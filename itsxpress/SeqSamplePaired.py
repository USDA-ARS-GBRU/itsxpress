
import subprocess
import logging
import os

from itsxpress.SeqSample import SeqSample
from itsxpress.definitions import ROOT_DIR, taxa_choices, taxa_dict, maxmismatches, maxratio, vsearch_fastq_qmax

logger = logging.getLogger(__name__)

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

    def _merge_reads(self, threads,stagger):
        try:
            seq_file = os.path.join(self.tempdir, 'seq.fq')
            if self.r1.split('.')[-1] == 'zst' and self.fastq2.split('.')[-1] == 'zst':
                r1_temp = os.path.join(self.tempdir, 'r1_temp.fq')
                r2_temp = os.path.join(self.tempdir, 'r2_temp.fq')
                parameters = ['zstd','-d', self.r1,'-o', r1_temp]
                p1 = subprocess.run(parameters, stderr=subprocess.PIPE)
                p1.check_returncode()
                logging.info(p1.stderr.decode('utf-8'))
                parameters = ['zstd','-d', self.fastq2,'-o', r2_temp]
                p1 = subprocess.run(parameters, stderr=subprocess.PIPE)
                p1.check_returncode()
                logging.info(p1.stderr.decode('utf-8'))
                self.r1 = r1_temp
                self.fastq2 = r2_temp

            if stagger:
                parameters = ['vsearch',
                          '--fastq_mergepairs', self.r1,
                          '--reverse', self.fastq2,
                          '--fastqout', seq_file,
                          '--fastq_maxdiffs', str(maxmismatches),
                          '--fastq_maxee', str(2),
                          '--threads', str(threads),
                          '--fastq_allowmergestagger',
                          '--fastq_qmax', str(vsearch_fastq_qmax)]
            else:
                parameters = ['vsearch',
                          '--fastq_mergepairs', self.r1,
                          '--reverse', self.fastq2,
                          '--fastqout', seq_file,
                          '--fastq_maxdiffs', str(maxmismatches),
                          '--fastq_maxee', str(2),
                          '--threads', str(threads),
                          '--fastq_qmax', str(vsearch_fastq_qmax)]
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