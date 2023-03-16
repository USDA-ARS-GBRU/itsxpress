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

from itsxpress.definitions import ROOT_DIR, taxa_choices, taxa_dict, maxmismatches, maxratio
from itsxpress.Dedup import Dedup
from itsxpress.ITSposition import ItsPosition
from itsxpress.SeqSamplePaired import SeqSamplePairedNotInterleaved
from itsxpress.SeqSampleNotPaired import SeqSampleNotPaired

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
    return parser


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
            elif file.endswith('.zst'):
                f = zstd.open(file, 'rt')
            else:
                f = open(file, 'r')
            lines = f.readlines()
            L1 = lines[0:96:8]
            L2 = lines[3:96:8]
            L1_old = [i.strip().split('/', 1)[0] for i in L1]
            L2_old = [i.strip().split('/', 1)[0] for i in L2]
            L1_new = [i.strip().split(' ', 1)[0] for i in L1]
            L2_new = [i.strip().split(' ', 1)[0] for i in L2]
            assert L1_old != L2_old or L1_new != L2_new

        except AssertionError as a:
            logging.error("'File may be interleaved. ITSxpress will run with errors. Check BBmap reformat.sh to split interleaved files before using ITSxpress.")
            raise a 
        except IOError as f:
            logging.error("File may be wrong format for interleaved file check.")
            raise f
        
  
    def core(file):   
        try:
            if file.endswith('.gz'):
                f = gzip.open(file, 'rt')
            elif file.endswith('.zst'):
                f = zstd.open(file, 'rt')
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
                elif file.endswith('.zst'):
                    f = zstd.open(file, 'rt')
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
    if fastq2:
        core(fastq2)

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
            sobj._merge_reads(threads=str(args.threads), stagger=args.allow_staggered_reads)
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
                dedup_obj.create_paired_trimmed_seqs(args.outfile, args.outfile2, gzipped=True,zstd_file = False,  itspos=its_pos,wri_file=True)
                dedup_obj.create_paired_trimmed_seqs(args.outfile, args.outfile2, gzipped=True,zstd_file = False,  itspos=its_pos,wri_file=False)
            elif args.outfile.split('.')[-1] == 'zst' and args.outfile2.split('.')[-1] == 'zst':
                dedup_obj.create_paired_trimmed_seqs(args.outfile, args.outfile2, gzipped=False, zstd_file = True, itspos=its_pos,wri_file=True)
                dedup_obj.create_paired_trimmed_seqs(args.outfile, args.outfile2, gzipped=False, zstd_file = True, itspos=its_pos,wri_file=False)
            else:
                dedup_obj.create_paired_trimmed_seqs(args.outfile, args.outfile2, gzipped=False,zstd_file = False,  itspos=its_pos,wri_file=True)
                dedup_obj.create_paired_trimmed_seqs(args.outfile, args.outfile2, gzipped=False,zstd_file = False,  itspos=its_pos,wri_file=False)

        else:
            if args.outfile.split('.')[-1] =='gz':
                dedup_obj.create_trimmed_seqs(args.outfile, gzipped=True,zstd_file = False, itspos=its_pos,wri_file=True)
                dedup_obj.create_trimmed_seqs(args.outfile, gzipped=True,zstd_file = False, itspos=its_pos,wri_file=False)
                #add function with above create_trimmed_seqs
                #use said function to check for 0 length seqs
            if args.outfile.split('.')[-1] =='zst':
                dedup_obj.create_trimmed_seqs(args.outfile, gzipped=False, zstd_file = True, itspos=its_pos,wri_file=True)
                dedup_obj.create_trimmed_seqs(args.outfile, gzipped=False, zstd_file = True, itspos=its_pos,wri_file=False)
            else:
                dedup_obj.create_trimmed_seqs(args.outfile, gzipped=False,zstd_file = False, itspos=its_pos,wri_file=True)
                dedup_obj.create_trimmed_seqs(args.outfile, gzipped=False,zstd_file = False, itspos=its_pos,wri_file=False)
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
                shutil.rmtree(sobj.tempdir)
        except UnboundLocalError:
            pass
        except AttributeError:
            pass


if __name__ == '__main__':
    main()
