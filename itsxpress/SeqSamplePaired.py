from itsxpress.SeqSample import SeqSample



import subprocess
import logging
import os
import pyzstd as zstd

from Bio import SeqIO

from itsxpress.definitions import ROOT_DIR, taxa_choices, taxa_dict, maxmismatches, maxratio

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
                self.r1 = zstd.decompress(self.r1)
                self.fastq2 = zstd.decompress(self.fastq2)
            if stagger:
                parameters = ['vsearch',
                          '--fastq_mergepairs' , self.r1,
                          '--reverse' , self.fastq2,
                          '--fastqout' ,seq_file,
                          '--fastq_maxdiffs' , str(maxmismatches),
                          '--fastq_maxee' , str(2),
                          '--threads'  ,str(threads),
                          '--fastq_allowmergestagger']
            else:
                parameters = ['vsearch',
                          '--fastq_mergepairs' , self.r1,
                          '--reverse' , self.fastq2,
                          '--fastqout' ,seq_file,
                          '--fastq_maxdiffs' , str(maxmismatches),
                          '--fastq_maxee' , str(2),
                          '--threads'  ,str(threads)]
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