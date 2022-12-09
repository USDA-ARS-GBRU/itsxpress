#!/usr/bin/env python
"""ITSxpress: A python module to rapidly trim ITS amplicon sequences from FASTQ files.

Authors: Adam Rivers, Kyle weber, USDA Agricultural Research Service

The internally transcribed spacer region is a region between the highly conserved small
subunit (SSU) of rRNA and the large subunit (LSU) of the rRNA. The eukaryotic ITS contains
the 5.8s gene and two variable length spacer regions. In amplicon sequencing studies it is
common practice to trim off the conserved (SSU, 5,8S or LSU) regions. Bengtsson-Palme
et al. (2013) published software the software package ITSx to do this.

ITSxpress is a high-speed implementation of the methods in ITSx than also allows FASTQ
files to be processed. Processing FASTQ files Which is essential for analyzing
sequences using the newer exact Sequence Variant methods in Qiime2, Dada2, Deblur
and Unoise that are replacing OTU clustering.

ITSxpress is also available as a QIIME Plugin. See
https://github.com/USDA-ARS-GBRU/q2_itsxpress for details.

Process:
    * Merges and error corrects reads using bbduk if reads are paired-end
    * Deduplicates reads using Vmatch to eliminate redundant hmm searches
    * Searches for conserved regions using the ITSx hmms, using HMMsearch:
    * Parses everything in python returning (optionally gzipped) fastq files.

Reference:
    Johan Bengtsson-Palme, Vilmar Veldre, Martin Ryberg, Martin Hartmann, Sara Branco,
    Zheng Wang, Anna Godhe, Yann Bertrand, Pierre De Wit, Marisol Sanchez,
    Ingo Ebersberger, Kemal Sanli, Filipe de Souza, Erik Kristiansson, Kessy Abarenkov,
    K. Martin Eriksson, R. Henrik Nilsson. (2013). ITSx: Improved software detection
    and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other
    eukaryotes for use in environmental sequencing. Methods in Ecology and Evolution,
    4: 914-919, 2013 (DOI: 10.1111/2041-210X.12073)
"""



import gzip
import tempfile
import argparse
import subprocess
import logging
import time
import os
import shutil
import math
from itertools import tee

from numpy import empty

from Bio import SeqIO

from itsxpress.definitions import ROOT_DIR, taxa_choices, taxa_dict, maxmismatches, maxratio,maxee

def restricted_float(x):
    x = float(x)
    if x < 0.99 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.99, 1.0]"%(x,))
    return x


def myparser():
    parser = argparse.ArgumentParser(description='ITSxpress: A python module to rapidly \
        trim ITS amplicon sequences from Fastq files.')
    parser.add_argument('--fastq', '-f', type=str, required=True,
                        help='A .fastq, .fq, .fastq.gz or .fq.gz file. Interleaved or not.')
    parser.add_argument('--single_end', '-s', action='store_true', default=False,
                        help='A flag to specify that the FASTQ file is single-ended (not paired). Default is false.')
    parser.add_argument('--fastq2', '-f2', type=str, default=None,
                        help='A .fastq, .fq, .fastq.gz or .fq.gz file. representing read 2 (optional)')
    parser.add_argument('--outfile', '-o', type=str, help="the trimmed Fastq file, if it \
                        ends in 'gz' it will be gzipped", required=True)
    parser.add_argument('--outfile2', '-o2', type=str, help="the trimmed read 2 Fastq file, if it \
                            ends in 'gz' it will be gzipped. If provided, reads will be returned unmerged.", default=None)
    parser.add_argument('--tempdir', help='The temp file directory', default=None)
    parser.add_argument('--keeptemp' ,help="Should intermediate files be kept?", action='store_true')
    parser.add_argument('--region', help='', choices=["ITS2", "ITS1", "ALL"], required=True)
    parser.add_argument('--taxa', help='The taxonomic group sequenced.', choices=taxa_choices, default="Fungi")
    parser.add_argument('--cluster_id', help='The percent identity for clustering reads range [0.99-1.0], set to 1 for exact dereplication.', type=restricted_float, default=1.0)
    parser.add_argument('--reversed_primers', '-rp',  help="Primers are in reverse orientation as in Taylor et al. 2016, DOI:10.1128/AEM.02576-16. If selected ITSxpress returns trimmed reads flipped to the forward orientation", action='store_true')
    parser.add_argument('--log' ,help="Log file", default="ITSxpress.log")
    parser.add_argument('--threads' ,help="Number of processor threads to use.", type=int, default=1)
    return parser




class ItsPosition:
    """Class for ITS positional information derived from hmmserach domtable files.

    Args:
        domtable (str):  the path locating the domtable file from HMMER 3 hmmsearch.
        region (str): The region of the ITS to extract choices: ["ITS1", "ITS2", "ALL"].

    Attributes:
        ddict (dict): A dictionary holding the scores and start and stop
            positions for the selected segment of each sequence.
            Example: {sample:{left:{score:31, pos:15}, right:{score:32, pos:354}, tlen:449}
        leftprefix (str): the left prefix to search for (set by type variable).
        rightprefix (str): the right prefix to search for (set by type variable).

    """
    def _score(self, sequence, stype, score, from_pos, to_pos, tlen):
        """Evaluates scores and positions from the new line of a domtable file and
            updates ddict if necessary.

        Args:
            sequence (str): The name of the sequence.
            stype (str): {'left', 'right'}
            score (int): The bit score from HMMSearch.
            to_pos (int): the ending position of the left sequence.
            from_pos (int): The beginning position of the right sequence.
            tlen (int): the length of the sequence

        """

        if stype in self.ddict[sequence]:
            if score > self.ddict[sequence][stype]["score"]:
                self.ddict[sequence][stype]["score"] = score
                self.ddict[sequence][stype]["to_pos"] = to_pos
                self.ddict[sequence][stype]["from_pos"] = from_pos
        else:
            self.ddict[sequence][stype] = {}
            self.ddict[sequence][stype]["score"] = score
            self.ddict[sequence][stype]["to_pos"] = to_pos
            self.ddict[sequence][stype]["from_pos"] = from_pos
            self.ddict[sequence]["tlen"] = tlen


    def parse(self):
        """Parses dom table from HMMsearch.

        The dom table is parsed and the start and stop position from the top scoring
        hmm math is saved. The start and stop positions of reach sequence are added to the ddict attribute.

        """
        try:
            with open(self.domtable, 'r') as f:
                for num, line in enumerate(f):
                    if not line.startswith("#"):
                        ll = line.split()
                        sequence = ll[0]
                        hmmprofile = ll[3]
                        score = float(ll[13])
                        from_pos = int(ll[19])
                        to_pos = int(ll[20])
                        tlen = int(ll[2])
                        if sequence not in self.ddict:
                            self.ddict[sequence] = {}
                        if hmmprofile.startswith(self.leftprefix):
                            self._score(sequence, 'left', score, from_pos, to_pos, tlen)
                        elif hmmprofile.startswith(self.rightprefix):
                            self._score(sequence, 'right', score, from_pos, to_pos, tlen)
        except Exception as e:
            logging.error("Exception occurred when parsing HMMSearh results")
            raise e

    def __init__(self, domtable, region):
        self.domtable = domtable
        self.ddict = {}
        if region == "ITS2":
            self.leftprefix = '3_'
            self.rightprefix = '4_'
        elif region == "ITS1":
            self.leftprefix = '1_'
            self.rightprefix = '2_'
        elif region == "ALL":
            self.leftprefix = '1_'
            self.rightprefix = '4_'
        self.parse()


    def get_position(self, sequence):
        """ Returns the start and stop positions for a given sequence.

        Args:
            sequence (str): The name of the sequence.

        Returns:
            (tuple): (start position, end position) zero indexed

        Raises:
            KeyError: If input sequence is not present in dictionary (no ITS start or stop sites were found)

        """

        try:
            if "left" in self.ddict[sequence]:
                start = int(self.ddict[sequence]["left"]["to_pos"])
            else:
                start = None
            if "right" in self.ddict[sequence]:
                stop = int(self.ddict[sequence]["right"]["from_pos"]) - 1
            else:
                stop = None
            if "tlen" in self.ddict[sequence]:
                tlen = int(self.ddict[sequence]["tlen"])
            else:
                tlen = None
            return(start, stop, tlen)
        except KeyError:
            logging.debug("No ITS stop or start sites were identified for sequence {}, skipping.".format(sequence))
            raise KeyError

class Dedup:
    """A class to handle deduplicated sequence data.

    To speed processing Vmatch is used to remove duplicate amplicons so that the
    start and stop sites are determined only once.

    Attributes:
        matchdict (dict): a dictionary of each sequence ID as a key and
            its representative sequence ID as a value {seq1:rep1, seq2:rep1, seq3:rep2}.
        uc_file (str): the location of the .uc file containing matching information.
        rep_file (str): The location of the representative sequences file.
        seq_file (str): The location of the complete sequences file.
        fastq (str): the location of the input fastq
        fastq2 (str) the location of optional Read2 input fastq if paired

    """


    def parse(self):
        """Parse the uc data file to populate the matchdict attribute.

        Raises:
            Exception: General exception if uc file is not parsed properly

        """
        try:
            with open(self.uc_file, 'r') as f:
                self.matchdict = {}
                for line in f:
                    ll = line.split()
                    datatype = ll[0]
                    ref = ll[9]
                    seq = ll[8]
                    if datatype == 'S':
                        self.matchdict[seq] = seq
                    elif datatype == 'H':
                        self.matchdict[seq] = ref
        except Exception as e:
            logging.exception("Could not parse the Vsearch '.uc' file.")
            raise e

    def __init__(self, uc_file, rep_file, seq_file, fastq=None, fastq2=None):
        self.matchdict = None
        self.uc_file = uc_file
        self.rep_file = rep_file
        self.seq_file = seq_file
        self.fastq = fastq
        self.fastq2 = fastq2
        self.parse()


    def _get_paired_seq_generator(self, zipseqgen, itspos):
        """This function takes a zipped object of two Biopython SeqIO sequence generators, and
        returns a two generators of Biopython SeqRecords for Dada2. Sequences where the ITS ends could
        not be determined are omitted.

        Args:
            zipseqgen (obj): A zipped object with two Biopython SeqIO generators
                          for the forward and reverse input sequences
            ispos (obj): An itsxpress ItsPosition object

        Returns:
            (obj): A two Python SeqRecord generators that yield filtered, trimmed sequence records.

        """
        def _filterfunc(ziprecord):
            """ Filters records down to those that contain a valid ITS start and stop position

            Args:
                record (obj): a Biopython SeqRecord object

            Returns:
                bool: True if an ITS start and stop positions are present false otherwise

            """
            try:
                record1, record2 = ziprecord
                if record1.id in self.matchdict:
                    repseq = self.matchdict[record1.id]
                    start, stop, tlen = itspos.get_position(repseq)
                    if start and stop:
                        if start < stop:
                            return True
                else:
                    return False
            except KeyError:
                return False



        def _map_func(ziprecord):
            """Trims the record down to the selected ITS region

            Args:
                record (obj): a Biopython SeqRecord object

            Returns:
                obj: two Biopython SeqRecord objects; forward and reverse reads trimmed to the ITS region
            """
            record1, record2 = ziprecord
            repseq = self.matchdict[record1.id]
            start, stop, tlen = itspos.get_position(repseq)
            r2start = tlen - stop
            return record1[start:], record2[r2start:]

        def _split_gen(gen):
            gen_a, gen_b = tee(gen, 2)
            return (a for a, b in gen_a), (b for a, b in gen_b)

        filt = filter(_filterfunc, zipseqgen)
        gen1 = map(_map_func, filt)
        gen1_split_a, gen1_split_b = _split_gen(gen1)
        #print(list(gen1_split_a))
        zeroseqctr = 0
        for i in list(gen1_split_a):
            if i.seq == "":
                zeroseqctr=zeroseqctr+1
                print(i.id)
        print("Total number of sequences that are empty Split A: ",zeroseqctr)
        print("  ")
        #print(list(gen1_split_b))
        zeroseqctr = 0
        for i in list(gen1_split_b):
            if i.seq == "":
                zeroseqctr=zeroseqctr+1
                print(i.id)
        print("Total number of sequences that are empty Split B: ",zeroseqctr)
        return gen1_split_a, gen1_split_b


    def create_paired_trimmed_seqs(self, outfile1, outfile2, gzipped, itspos):
        """Writes two FASTQ files, optionally gzipped, with the reads trimmed to the
            selected region.

        Args:
            outfile1 (str): The file to write the forward sequences to.
            outfile2 (str): The file to write the reverse sequences to.
            gzip (bool): Should the output files be gzipped?
            itspos (object): an ItsPosition object

        Returns:
            str: Name of the file written

        """

        def _write_seqs(seqs, outfile):
            """Helper function to optionally write sequences in compressed format

            Args:
                seqs (obj): A biopython SeqRecord generators
                outfile (str): A file name to writ the fastq data to.

            """
            if gzipped:
                with gzip.open(outfile, 'wt') as g:
                    SeqIO.write(seqs, g, "fastq")

            else:
                with open(outfile, 'w') as g:
                    SeqIO.write(seqs, g, "fastq")

        def _create_gen(f, g):
            """Create a sequence generator

            Args:
                f (str): a file name for read 1 fastq data
                g (str): a file name for read 2 fastq data

            """
            seqgen1 = SeqIO.parse(f, 'fastq')
            seqgen2 = SeqIO.parse(g, 'fastq')
            zipseqgen = zip(seqgen1, seqgen2)
            seqs1, seqs2 = self._get_paired_seq_generator(zipseqgen, itspos)
            _write_seqs(seqs1, outfile1)
            _write_seqs(seqs2, outfile2)

        try:
            if self.fastq.endswith(".gz") and self.fastq2.endswith(".gz"):
                with gzip.open(self.fastq, 'rt') as f:
                    with gzip.open(self.fastq2, 'rt') as g:
                        _create_gen(f, g)

            elif not (self.fastq2.endswith(".gz") or self.fastq2.endswith(".gz")):
                with open(self.fastq, 'r') as f:
                    with open(self.fastq2, 'r') as g:
                        _create_gen(f, g)
            else:
                raise ValueError("Fastq and Fastq2 files should both be gzipped or both be un-gzipped. Mixed input is not accepted.")

        except Exception as e:
            raise e


    def _get_trimmed_seq_generator(self, seqgen, itspos):
        """This function takes a Biopython SeqIO sequence generator, and
        returns a generator of trimmed sequences suitable for Deblur. Sequences where the ITS ends could
        not be determined are omitted.

        Args:
            seqgen (obj): A Biopython SeqIO generator of all input sequences
            ispos (obj): An itsxpress ItsPosition object

        Returns:
            (obj): A map object generator that yields filtered, trimmed sequence records.

        """
        def _filterfunc(record):
            """ Filters records down to those that contain a valid ITS start and stop position

            Args:
                record (obj): a Biopython SeqRecord object

            Returns:
                bool: True if an ITS start and stop positions are present false otherwise

            """
            try:
                if record.id in self.matchdict:
                    repseq = self.matchdict[record.id]
                    start, stop, tlen = itspos.get_position(repseq)
                    if start and stop:
                        if start < stop:
                            return True
                else:
                    return False
            except KeyError:
                return False



        def map_func(record):
            """Trims the record down to the selected ITS region

            Args:
                record (obj): a Biopython SeqRecord object

            Returns:
                obj: a Biopython SeqRecord object trimmed to the ITS region
            """
            repseq = self.matchdict[record.id]
            start, stop, tlen = itspos.get_position(repseq)
            return record[start:stop]

        filt = filter(_filterfunc, seqgen)
        r1 = map(map_func, filt)
        zeroseqctr = 0
        for i in list(r1):
            if i.seq == "":
                zeroseqctr=zeroseqctr+1
                print(i.id)
        print("Total number of sequences that are empty: ",zeroseqctr)
        return r1


    def create_trimmed_seqs(self, outfile, gzipped, itspos):
        """Creates a FASTQ file, optionally gzipped, with the reads trimmed to the
            selected region.

        Args:
            outfile (str): The file to write the sequences to.
            gzip (bool): Should the files be gzipped?
            itspos (object): an ItsPosition objectconda activate /project/gbru_fy21_tomato_ralstonia/ITSxpress/software/qiime2-2022.8_ITSxpressV2

        """
        def _write_seqs():
            if gzipped:
                tempf = os.path.join('./','temp.fa')
                with open(tempf, 'w') as g:
                    SeqIO.write(seqs, g, "fastq")
                    #print(seqs)
                with open(tempf,'rb') as f_in:
                    with gzip.open(outfile,'wb') as f_out:
                        f_out.writelines(f_in)

            else:
                with open(outfile, 'w') as g:
                    #print(seqs)
                    SeqIO.write(seqs, g, "fastq")


        if self.seq_file.endswith(".gz"):
            with gzip.open(self.seq_file, 'rt') as f:
                seqgen = SeqIO.parse(f, 'fastq')
                seqs = self._get_trimmed_seq_generator(seqgen, itspos)
                #print(seqs)
                _write_seqs()

        else:
            with open(self.seq_file, 'r') as f:
                seqgen = SeqIO.parse(f, 'fastq')
                seqs = self._get_trimmed_seq_generator(seqgen, itspos)
                #print(seqs)
                _write_seqs()


class SeqSample:
    """The class for processing sequence data into trimmed sequences.

    Attributes:
        tempdir (obj): A temporary directory object
        fastq (str): The path to the input FASTQ file
        uc_file (str): The path to the Vsearch uc mapping file
        rep_file: (str) the path to the representative sequences FASTA file created by Vsearch
        seq_file (str): the location of the fastq or fastq.gz sequence file used for analysis

    """


    def __init__(self, fastq, tempdir=None):
        if tempdir:
            if not os.path.exists(tempdir):
                logging.warning("Specified location for tempfile ({}) does not exist, using default location.".format(tempdir))
                self.tempdir = tempfile.mkdtemp(prefix='itsxpress_')
            else:
                self.tempdir = tempfile.mkdtemp(prefix='itsxpress_', dir=tempdir)
        else:
            self.tempdir = tempfile.mkdtemp(prefix='itsxpress_')
        self.fastq = fastq
        self.uc_file = None
        self.rep_file = None
        self.dom_file = None
        self.seq_file = None



    def deduplicate(self, threads=1):
        """Runs Vsearch dereplication to create a FASTA file of non-redundant sequences.

        Args:
            threads (int or str):the number of processor threads to use

        """
        try:
            self.uc_file=os.path.join(self.tempdir, 'uc.txt')
            self.rep_file=os.path.join(self.tempdir,'rep.fa')
            # parameters = ["vsearch",
            #               "--derep_fulllength",
            #               self.seq_file,
            #               "--output", self.rep_file,
            #               "--uc", self.uc_file,
            #               "--strand", "both",
            #               "--threads", str(threads)]
            parameters = ["vsearch",
                          "--fastx_uniques",
                          self.seq_file,
                          "--fastaout", self.rep_file,
                          "--uc", self.uc_file,
                          "--strand", "both"]
            p2 = subprocess.run(parameters, stderr=subprocess.PIPE)
            logging.info(p2.stderr.decode('utf-8'))
            p2.check_returncode()
        except subprocess.CalledProcessError as e:
            logging.exception("Could not perform dereplication with Vsearch. Error from Vsearch was:\n {}".format(p2.stderr.decode('utf-8')))
            raise e
        except FileNotFoundError as f:
            logging.error("Vsearch was not found, make sure Vsearch is installed and executable")
            raise f


    def cluster(self, threads=1, cluster_id=0.995):
        """Runs Vsearch clustering to create a FASTA file of non-redundant sequences.

        Args:
            threads (int or str):the number of processor threads to use

        """
        try:
            self.uc_file = os.path.join(self.tempdir, 'uc.txt')
            self.rep_file = os.path.join(self.tempdir,'rep.fa')
            parameters = ["vsearch",
                          "--cluster_size", self.seq_file,
                          "--centroids", self.rep_file,
                          "--uc", self.uc_file,
                          "--strand", "both",
                          "--id", str(cluster_id),
                          "--threads", str(threads)]
            p2 = subprocess.run(parameters, stderr=subprocess.PIPE)
            print(p2.stderr.decode('utf-8'))
            p2.check_returncode()
        except subprocess.CalledProcessError as e:
            logging.exception("Could not perform clustering with Vsearch. Error from Vsearch was:\n {}".format(p2.stderr.decode('utf-8')))
            raise e
        except FileNotFoundError as f:
            logging.error("Vsearch was not found, make sure Vsearch is installed and executable")
            raise f

    def _search(self, hmmfile, threads=1):
        try:
            self.dom_file = os.path.join(self.tempdir, 'domtbl.txt')
            #Run Hmmsearch
            parameters = ["hmmsearch",
                          "--domtblout",
                          self.dom_file,
                          "-T", "10",
                          "--cpu", str(threads),
                          "--tformat", "fasta",
                          "--F1", "1e-6",
                          "--F2", "1e-6",
                          "--F3", "1e-6",
                          hmmfile,
                          self.rep_file]
            p4 = subprocess.run(parameters, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
            p4.check_returncode()
        except subprocess.CalledProcessError as e :
            logging.exception("Could not perform ITS identification with hmmserach. The error was:\n {}".format(p4.stderr.decode('utf-8')))
            raise e
        except FileNotFoundError as f:
            logging.error("hmmsearch was not found, make sure HMMER3 is installed and executable")
            raise f


class SeqSamplePairedNotInterleaved(SeqSample):
    """SeqSample class extended to paired, two FASTQ file format.

    """
    def __init__(self, fastq, tempdir, fastq2, reversed_primers=False ):
        SeqSample.__init__(self, fastq, tempdir)
        if reversed_primers:
            self.r1 = fastq2
            self.fastq2 = fastq
        else:
            self.r1 = fastq
            self.fastq2 = fastq2

    def _merge_reads(self, threads):
        try:
            seq_file = os.path.join(self.tempdir, 'seq.fq')
            parameters = ['vsearch',
                          '--fastq_mergepairs' , self.r1,
                          '--reverse' , self.fastq2,
                          '--fastqout' ,seq_file,
                          '--fastq_maxdiffs' , str(maxmismatches),
                          '--fastq_maxee' , str(maxee),
                          '--threads'  ,str(threads)]
                          #'--fastq_maxee_rate' , str(maxratio),
            p1 = subprocess.run(parameters, stderr=subprocess.PIPE)
            self.seq_file = seq_file
            p1.check_returncode()
            logging.info(p1.stderr.decode('utf-8'))
        except subprocess.CalledProcessError as e:
            logging.exception("Could not perform read merging with vsearch. Error from vsearch was: \n  {}".format(p1.stderr.decode('utf-8')))
            raise e
        except FileNotFoundError as f:
            logging.error("vsearch was not found, make sure vsearch is installed on this environment")
            raise f

class SeqSampleNotPaired(SeqSample):
    """SeqSample class extended to unpaired format.

    """

    def __init__(self, fastq, tempdir):
        SeqSample.__init__(self, fastq, tempdir)
        self.seq_file = self.fastq
        self.r1 = self.fastq
        self.fastq2 = None

## Utility Functions

def _is_paired(fastq, fastq2, single_end):
    """Determines the workflow based on file inputs.

    Args:

    """
    if fastq and fastq2:
        paired_end = True
    elif single_end:
        paired_end = False
    else:
        paired_end = True
    return paired_end

def _logger_setup(logfile):
    """Set up logging to a logfile and the terminal standard out.

    Args:
        fastq (str): The path to a fastq or fastq.gz file
        fastq2 (str): The path to a fastq or fastq.gz file for the reverse sequences

    """
    try:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                            datefmt='%m-%d %H:%M',
                            filename=logfile,
                            filemode='w')
        # define a Handler which writes INFO messages or higher to the sys.stderr
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        # set a format which is simpler for console use
        formatter = logging.Formatter('%(asctime)s: %(levelname)-8s %(message)s')
        # tell the handler to use this format
        console.setFormatter(formatter)
        # add the handler to the root logger
        logging.getLogger('').addHandler(console)
    except Exception as e:
        print("An error occurred setting up logging")
        raise e

def _check_fastqs(fastq, fastq2=None):
    """Verifies the input files are valid fastq or fastq.gz files.
    Add interleaved explanation

    Args:
        fastq (str): The path to a fastq or fastq.gz file
        fastq2 (str): The path to a fastq or fastq.gz file for the reverse sequences

    Raises:
        ValueError: If Biopython detected invalid FASTQ files
    """
    def check_interleaved(file):
        try:
            if file.endswith('.gz'):
                f = gzip.open(file, 'rt')
            else:
                f = open(file, 'r')
            lines = f.readlines()
            L1 = lines[0:96:8]
            L2 = lines[4:96:8]
            #logging.info(L1)
            #logging.info(L2)
            L1_old = [i.strip().split('/', 1)[0] for i in L1]
            L2_old = [i.strip().split('/', 1)[0] for i in L2]
            L1_new = [i.strip().split(' ', 1)[0] for i in L1]
            L2_new = [i.strip().split(' ', 1)[0] for i in L2]
            assert L1_old != L2_old or L1_new != L2_new

        except AssertionError as a:
            logging.error("'File seems to be interleaved. ITSxpress will run with errors. Check BBmap reformat.sh to split interleaved files.")
        except IOError as f:
            logging.error("File may be wrong format for interleaved file check.")
            raise f


    def core(file):
        try:
            if file.endswith('.gz'):
                f = gzip.open(file, 'rt')
            else:
                f = open(file, 'r')
            n = 0
            for record in SeqIO.parse(f, 'fastq'):
                while n < 100:
                    n += 1
            check_interleaved(file)

            f.close()
            if fastq2:
                if fastq2.endswith('.gz'):
                    f = gzip.open(fastq2, 'rt')
                else:
                    f = open(fastq2, 'r')
                n = 0
                for record in SeqIO.parse(f, 'fastq'):
                    while n < 100:
                        n += 1
            f.close()
        except ValueError as e:
            logging.error("There appears to be an issue with the format of input file {}.".format(file))
            raise e
        except FileNotFoundError as f:
            logging.error("The input file {} could not be found.".format(file))
            raise f
        except Exception as g:
            logging.error("There appears to be an issue reading the input file {}.".format(file))
            raise g

    core(fastq)

def main(args=None):
    """Run Complete ITS trimming workflow.

    """
    # Set up logging
    t0 = time.time()
    parser = myparser()
    if not args:
        args = parser.parse_args()

    _logger_setup(args.log)
    try:
        logging.info("Verifying the input sequences.")
        _check_fastqs(args.fastq, args.fastq2)
        # Parse input types
        paired_end = _is_paired(args.fastq, args.fastq2, args.single_end)
        if paired_end:
            logging.info("Sequences are paired-end in two files. They will be merged using Vsearch.")
            sobj = SeqSamplePairedNotInterleaved(fastq=args.fastq, fastq2=args.fastq2, tempdir=args.tempdir, reversed_primers=args.reversed_primers)
            sobj._merge_reads(threads=str(args.threads))
        elif not paired_end:
            logging.info("Sequences are assumed to be single-end.")
            sobj = SeqSampleNotPaired(fastq=args.fastq, tempdir=args.tempdir)
        logging.info("Temporary directory is: {}".format(sobj.tempdir))
        # Deduplicate
        logging.info("Unique sequences are being written to a temporary FASTA file with Vsearch.")
        if math.isclose(args.cluster_id, 1, rel_tol=1e-05):
            sobj.deduplicate(threads=str(args.threads))
        else:
            sobj.cluster(threads=str(args.threads),cluster_id=args.cluster_id)
        # HMMSearch for ITS regions
        logging.info("Searching for ITS start and stop sites using HMMSearch. This step takes a while.")
        hmmfile = os.path.join(ROOT_DIR, "ITSx_db","HMMs", taxa_dict[args.taxa])
        sobj._search(hmmfile=hmmfile, threads=str(args.threads))
        # Parse Hmmsearch output
        logging.info("Parsing HMM results.")
        its_pos = ItsPosition(domtable=sobj.dom_file, region=args.region)
        # Create deduplication object
        dedup_obj = Dedup(uc_file=sobj.uc_file, rep_file=sobj.rep_file, seq_file=sobj.seq_file, fastq=sobj.r1, fastq2=sobj.fastq2)
        # Create trimmed sequences
        logging.info("Writing out sequences")
        if args.outfile2:
            if args.outfile.split('.')[-1] == 'gz' and args.outfile2.split('.')[-1] == 'gz':
                dedup_obj.create_paired_trimmed_seqs(args.outfile, args.outfile2, gzipped=True, itspos=its_pos)
            else:
                dedup_obj.create_paired_trimmed_seqs(args.outfile, args.outfile2, gzipped=False, itspos=its_pos)
        else:
            if args.outfile.split('.')[-1] =='gz':
                dedup_obj.create_trimmed_seqs(args.outfile, gzipped=True, itspos=its_pos)
            else:
                dedup_obj.create_trimmed_seqs(args.outfile, gzipped=False, itspos=its_pos)
        t1 = time.time()
        fmttime = time.strftime("%H:%M:%S", time.gmtime(t1-t0))
        logging.info("ITSxpress ran in {}".format(fmttime))
    except Exception as e:
        logging.error("ITSxpress terminated with errors. See the log file for details.")
        logging.error(e)
        #raise e
        raise SystemExit(1)
    finally:
        try:
            if not args.keeptemp:
                shutil.rmtree(sobj.tempdir)
        except UnboundLocalError:
            pass
        except AttributeError:
            pass


if __name__ == '__main__':
    main()
