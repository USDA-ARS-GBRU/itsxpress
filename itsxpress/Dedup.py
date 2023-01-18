import logging
import gzip
import pyzstd as zstd
import os
from itertools import tee

from Bio import SeqIO

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


    def _get_paired_seq_generator(self, zipseqgen, itspos,wri_file):
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
        if wri_file == False:
            zeroseqctr1 = 0
            seqlist1=[]
            seqlist2=[]
            for i in list(gen1_split_a):
                if i.seq == "":
                    zeroseqctr1=zeroseqctr1+1
                    seqlist1.append(i.id)
            zeroseqctr2 = 0
            for i in list(gen1_split_b):
                if i.seq == "":
                    zeroseqctr2=zeroseqctr2+1
                    seqlist2.append(i.id)
            if zeroseqctr1 or zeroseqctr2 !=0:
                print("Total number of sequences that are empty Split A: ",zeroseqctr1)
                print("Sequence IDs: ")
                print(seqlist1)
                print("Total number of sequences that are empty Split B: ",zeroseqctr2)
                print("Sequence IDs: ")
                print(seqlist2) 
        return gen1_split_a, gen1_split_b


    def create_paired_trimmed_seqs(self, outfile1, outfile2, gzipped,zstd_file, itspos,wri_file):
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
            elif zstd_file:
                with zstd.open(outfile, 'wt')as g:
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
            seqs1, seqs2 = self._get_paired_seq_generator(zipseqgen, itspos,wri_file)
            if wri_file:
                _write_seqs(seqs1, outfile1)
                _write_seqs(seqs2, outfile2)

        try:
            if self.fastq.endswith(".gz") and self.fastq2.endswith(".gz"):
                with gzip.open(self.fastq, 'rt') as f:
                    with gzip.open(self.fastq2, 'rt') as g:
                        _create_gen(f, g)
            elif self.fastq.endswith(".zst") and self.fastq2.endswith(".zst"):
                with zstd.open(self.fastq,'rt') as f:
                    with zstd.open(self.fastq2, 'rt') as g:
                        _create_gen(f,g)
            elif  (self.fastq2.endswith(".fastq") or self.fastq2.endswith(".fastq")):
                with open(self.fastq, 'r') as f:
                    with open(self.fastq2, 'r') as g:
                        _create_gen(f, g)
            else:
                raise ValueError("Fastq and Fastq2 files should both be gzipped (.gz), zstd compressed (.zst) or both be uncompressed. Mixed input is not accepted.")

        except Exception as e:
            raise e


    def _get_trimmed_seq_generator(self, seqgen, itspos,wri_file):
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
        if wri_file == False:
            zeroseqctr = 0
            seqlist=[]
            r1 = map(map_func, filt)
            for i in list(r1):
                if i.seq == "":
                    zeroseqctr=zeroseqctr+1
                    seqlist.append(i.id)
            
            if zeroseqctr !=0:
                print("Total number of sequences that are empty: ",zeroseqctr)
                print("Sequence IDs: ")
                print(seqlist)
        return map(map_func, filt)


    def create_trimmed_seqs(self, outfile, gzipped,zstd_file, itspos,wri_file):
        """Creates a FASTQ file, optionally gzipped, with the reads trimmed to the
            selected region.
        Args:
            outfile (str): The file to write the sequences to.
            gzip (bool): Should the files be gzipped?
            itspos (object): an ItsPosition object
            wri_file (bool): Should file be written or checked for empty sequences?
        """
        def _write_seqs():
            if gzipped:
                tempf = os.path.join('./','temp.fa')
                with open(tempf, 'w') as g:
                    SeqIO.write(seqs, g, "fastq")
                with open(tempf,'rb') as f_in:
                    with gzip.open(outfile,'wb') as f_out:
                        f_out.writelines(f_in)
            elif zstd_file:
                tempf = os.path.join('./','temp.fa')
                with open(tempf, 'w') as g:
                    SeqIO.write(seqs, g, "fastq")
                with open(tempf,'rb') as f_in:
                    with zstd.open(outfile,'wb') as f_out:
                        f_out.writelines(f_in)            
            else:
                with open(outfile, 'w') as g:
                    SeqIO.write(seqs, g, "fastq")
                    

        if self.seq_file.endswith(".gz"):
            with gzip.open(self.seq_file, 'rt') as f:
                seqgen = SeqIO.parse(f, 'fastq')
                seqs = self._get_trimmed_seq_generator(seqgen, itspos,wri_file)
                if wri_file:
                    _write_seqs()
        elif self.seq_file.endswith(".zst"):
            with zstd.open(self.seq_file, 'rt') as f:
                seqgen = SeqIO.parse(f, 'fastq')
                seqs = self._get_trimmed_seq_generator(seqgen, itspos,wri_file)
                if wri_file:
                    _write_seqs()
        else:
            with open(self.seq_file, 'r') as f:
                seqgen = SeqIO.parse(f, 'fastq')
                seqs = self._get_trimmed_seq_generator(seqgen, itspos,wri_file)
                if wri_file:
                    _write_seqs()