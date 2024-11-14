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
Process:
    * Merges and error corrects reads using Vsearch if reads are paired-end
    * Deduplicates reads using Vsearch to eliminate redundant hmm searches
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
import pyzstd as zstd
import argparse
import logging
import time
import os
import shutil
import math
import tempfile
from itertools import tee, islice
import contextlib

from numpy import empty

from Bio import SeqIO

from itsxpress.definitions import ROOT_DIR, taxa_choices, taxa_dict, maxmismatches, maxratio
from itsxpress.Dedup import Dedup
from itsxpress.ITSposition import ItsPosition
from itsxpress.SeqSamplePaired import SeqSamplePairedNotInterleaved
from itsxpress.SeqSampleNotPaired import SeqSampleNotPaired
from ._version import __version__

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
    parser.add_argument('--allow_staggered_reads', help='Allow merging of staggered reads with --fastq_allowmergestagger \
                        for Vsearch --fastq_mergepairs. See Vsearch documentation. (Optional) Default is true.', default=True)
    parser.add_argument('--keeptemp' ,help="Should intermediate files be kept?", action='store_true')
    parser.add_argument('--region', help='', choices=["ITS2", "ITS1", "ALL"], required=True)
    parser.add_argument('--taxa', help='The taxonomic group sequenced.', choices=taxa_choices, default="Fungi")
    parser.add_argument('--cluster_id', help='The percent identity for clustering reads range [0.99-1.0], set to 1 for exact dereplication.', type=restricted_float, default=1.0)
    parser.add_argument('--reversed_primers', '-rp',  help="Primers are in reverse orientation as in Taylor et al. 2016, DOI:10.1128/AEM.02576-16. If selected ITSxpress returns trimmed reads flipped to the forward orientation", action='store_true')
    parser.add_argument('--log' ,help="Log file", default="ITSxpress.log")
    parser.add_argument('--threads' ,help="Number of processor threads to use.", type=int, default=1)
    parser.add_argument('--version', '-v' , help="display version",  action="version", version="ITSxpress version: " + __version__)
    return parser


## Utility Functions

def _is_paired(fastq, fastq2, single_end):
    """Determines the workflow based on file inputs.
    Args:
    fastq (str): The path to a fastq or fastq.gz file
    fastq2 (str): The path to a fastq or fastq.gz file for the reverse sequences
    single_end (bool): A flag to specify that the FASTQ file is single-ended (not paired). Default is false.
    """
    paired_end = None
    if fastq and fastq2:
        paired_end = True
    elif single_end:
        paired_end = False
    elif fastq and not fastq2:
        paired_end = False
        logging.info("Only one fastq file provided. Assuming single-end.")
    else:
        try:
            assert paired_end != None
        except AssertionError as a:
            logging.error("ITSxpress requires either a single-end file or two paired-end files. If this is a single-end file, please use the --single_end flag.")
            raise a
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

@contextlib.contextmanager
def read_file(filename: str, mode: str ='r') -> contextlib.contextmanager:
    """ A context manager for opening gzipped, stdz compressed or uncompressed files.

        Args:
            filename: the filename to open
            mode: 'r' for read 'w' for write. 't' will be adeed automacially for compressed files.
        
        Returns a context manager so the file can be opened using 'with'
    """
    try:
        if filename.endswith('.gz'):
            file = gzip.open(filename, mode +'t')
            yield file
        elif filename.endswith('.zst'):
            file = zstd.open(filename, mode + 't')
            yield file
        else:
            file = open(filename, mode)
            yield file
    except FileNotFoundError as f:
        logging.error("The input file {} could not be found.".format(filename))
        raise f
    except Exception as g:
        logging.error("There appears to be an issue reading the input file {}.".format(filename))
        raise g
    finally:
        if 'file' in locals():
            file.close()
    
def _check_fastqs(fastq: str, fastq2: str=None) -> None:
    """Verifies the input files are valid fastq or fastq.gz files. Also checks for interleaved files which are no longer supported.
    Args:
        fastq (str): The path to a fastq or fastq.gz file
        fastq2 (str): The path to a fastq or fastq.gz file for the reverse sequences

    Raises:
        ValueError: If Biopython detected invalid FASTQ files
    """
    # Validation method from BBtools (Thanks Brian!)
    def test_pair_names_str(id1, id2):
        len1 = len(id1)
        len2 = len(id2)
        if len1 != len2:
            return False  # Can happen in PacBio names, but never Illumina
        idx_slash1 = id1.rfind('/')
        idx_slash2 = id2.rfind('/')
        idx_space1 = id1.find(' ')
        idx_space2 = id2.find(' ')
        
        if idx_space1 == idx_space2 and idx_space1 > 0 and len1 >= idx_space1 + 3 and len2 >= idx_space2 + 3:
            if id1[idx_space1 + 1] == '1' and id1[idx_space1 + 2] == ':' and id2[idx_space2 + 1] == '2' and id2[idx_space2 + 2] == ':':
                for i in range(idx_space1):
                    if id1[i] != id2[i]:
                        return False
                return True
        
        if idx_slash1 == idx_slash2 and idx_slash1 > 0 and len1 >= idx_slash1 + 2 and len2 >= idx_slash2 + 2:
            if id1[idx_slash1 + 1] == '1' and id2[idx_slash2 + 1] == '2':
                for i in range(idx_slash1):
                    if id1[i] != id2[i]:
                        return False
                for i in range(idx_slash1 + 2, len1):
                    if id1[i] != id2[i]:
                        return False
                return True
        
        return id1 == id2
    
    def core(filename):
        with read_file(filename) as handle:
            records = SeqIO.parse(handle, 'fastq')
            reclist = list(islice(records, 2))
            return test_pair_names_str(reclist[0].id, reclist[1].id)
        
    def warn_mess(filename):
        return logging.warning( "The file {} may be interleaved, which is not supported. Please verify your input file manually.")
    
    if core(fastq):
        warn_mess(fastq)
    if fastq2:
        if core(fastq2):
            warn_mess(fastq2)

def _check_total_reads(file, file2 = None):
    """Check the total number of reads in the input file(s).
    """
    #Count every fourth line in fastq file.
    def core(file):
        with read_file(file) as handle:
            n = 0
            for i, _ in enumerate(handle):
                if i % 4 == 0:
                    n += 1
            return n
    
    reads = core(file)
    logging.info("Total number of reads in file {} is {}.".format(file, reads))
    #Why is file2 False? 
    if file2:
        reads = core(file2)
        logging.info("Total number of reads in file {} is {}.".format(file2, reads))

def create_temp_directory(tempdir_arg=None):
    """
    Creates a temporary directory at a user-defined location or at the default location.
    The directory name is prefixed with 'itsxpress_'.
    
    Ensures no file with the same name exists before creating the directory.
    
    Parameters:
    - tempdir_arg (str): Path to a directory provided by the user. If None, 
      a new temporary directory is created at the default location.
    
    Returns:
    - str: Path to the temporary directory, or None if there was an error.
    """
    try:
        if tempdir_arg:
            if not os.path.exists(tempdir_arg):
                os.makedirs(tempdir_arg)
                logging.info(f"Directory '{tempdir_arg}' has been created.")
            else:
                if os.path.isfile(tempdir_arg):
                    logging.error(f"A file with the same name '{tempdir_arg}' already exists. Cannot create directory.")
                    return None
                logging.info(f"Directory '{tempdir_arg}' already exists.")
            
            # Try creating a unique temporary directory inside user-specified directory
            temp_dir =  tempfile.mkdtemp(prefix="itsxpress_", dir=tempdir_arg)
            logging.info(f"Temporary directory '{temp_dir}' has been created at the user-defined location.")
        else:
            # Create a temporary directory at default location
            temp_dir = tempfile.mkdtemp(prefix="itsxpress_")
            logging.info(f"Temporary directory '{temp_dir}' has been created at the default location.")
        
        return temp_dir
    
    except Exception as e:
        logging.error(f"Failed to create temporary directory: {e}")
        return None


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
        logging.info("Starting ITSxpress version  {}".format(__version__))
        logging.info("Verifying the input sequences.")
        _check_fastqs(args.fastq, args.fastq2)
        # Parse input types
        paired_end = _is_paired(args.fastq, args.fastq2, args.single_end)
        session_tempdir = create_temp_directory(tempdir_arg=args.tempdir)
        if paired_end:
            logging.info("Sequences are paired-end in two files. They will be merged using Vsearch.")
            sobj = SeqSamplePairedNotInterleaved(fastq=args.fastq, fastq2=args.fastq2, tempdir=session_tempdir, reversed_primers=args.reversed_primers)
            sobj._merge_reads(threads=str(args.threads), stagger=args.allow_staggered_reads)
        elif not paired_end:
            logging.info("Sequences are assumed to be single-end.")
            sobj = SeqSampleNotPaired(fastq=args.fastq, tempdir=session_tempdir)
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
        # Trim sequences
        if args.outfile2:
            if args.outfile.split('.')[-1] == 'gz' and args.outfile2.split('.')[-1] == 'gz':
                dedup_obj.create_paired_trimmed_seqs(args.outfile, args.outfile2, gzipped=True,zstd_file = False,  itspos=its_pos,wri_file=True)
                #dedup_obj.create_paired_trimmed_seqs(args.outfile, args.outfile2, gzipped=True,zstd_file = False,  itspos=its_pos,wri_file=False)
            elif args.outfile.split('.')[-1] == 'zst' and args.outfile2.split('.')[-1] == 'zst':
                dedup_obj.create_paired_trimmed_seqs(args.outfile, args.outfile2, gzipped=False, zstd_file = True, itspos=its_pos,wri_file=True)
                #dedup_obj.create_paired_trimmed_seqs(args.outfile, args.outfile2, gzipped=False, zstd_file = True, itspos=its_pos,wri_file=False)
            else:
                dedup_obj.create_paired_trimmed_seqs(args.outfile, args.outfile2, gzipped=False,zstd_file = False,  itspos=its_pos,wri_file=True)
                #dedup_obj.create_paired_trimmed_seqs(args.outfile, args.outfile2, gzipped=False,zstd_file = False,  itspos=its_pos,wri_file=False)

        else:
            if args.outfile.split('.')[-1] == 'gz':
                dedup_obj.create_trimmed_seqs(args.outfile, gzipped=True,zstd_file = False, itspos=its_pos,wri_file=True,tempdir=sobj.tempdir)
                #dedup_obj.create_trimmed_seqs(args.outfile, gzipped=True,zstd_file = False, itspos=its_pos,wri_file=False)
                #add function with above create_trimmed_seqs
                #use said function to check for 0 length seqs
            elif args.outfile.split('.')[-1] == 'zst':
                dedup_obj.create_trimmed_seqs(args.outfile, gzipped=False, zstd_file = True, itspos=its_pos,wri_file=True,tempdir=sobj.tempdir)
                #dedup_obj.create_trimmed_seqs(args.outfile, gzipped=False, zstd_file = True, itspos=its_pos,wri_file=False)
            else:
                dedup_obj.create_trimmed_seqs(args.outfile, gzipped=False,zstd_file = False, itspos=its_pos,wri_file=True,tempdir=sobj.tempdir)
                #dedup_obj.create_trimmed_seqs(args.outfile, gzipped=False,zstd_file = False, itspos=its_pos,wri_file=False)
        # Count reads after trimming
        logging.info("Counting reads after trimming.")
        if args.outfile2:
            if args.fastq2:
                _check_total_reads(args.fastq, args.fastq2)
            else:
                _check_total_reads(args.fastq)
            _check_total_reads(args.outfile, args.outfile2)
        else:
            if args.fastq2:
                _check_total_reads(args.fastq, args.fastq2)
            else:
                _check_total_reads(args.fastq)
            _check_total_reads(args.outfile)
        
        t1 = time.time()
        fmttime = time.strftime("%H:%M:%S", time.gmtime(t1-t0))
        logging.info("ITSxpress ran in {}".format(fmttime))
    except Exception as e:
        logging.error("ITSxpress terminated with errors. See the log file for details.")
        logging.error(e)
        raise SystemExit(1)
    finally:
        try:
            if not args.keeptemp:
                shutil.rmtree(session_tempdir)
        except UnboundLocalError:
            pass
        except AttributeError:
            pass

if __name__ == '__main__':
    main()
