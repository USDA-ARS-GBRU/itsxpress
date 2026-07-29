import subprocess
import logging
import os
from typing import Union, Optional

from itsxpress.SeqSample import SeqSample
from itsxpress.definitions import ROOT_DIR, taxa_choices, taxa_dict, maxmismatches, maxratio, vsearch_fastq_qmax

logger = logging.getLogger(__name__)

class SeqSamplePairedNotInterleaved(SeqSample):
    """SeqSample class extended to paired, two FASTQ file format."""

    def __init__(self, fastq: str, tempdir: str, fastq2: str, reversed_primers: bool = False) -> None:
        """Initializes a SeqSamplePairedNotInterleaved object.

        Args:
            fastq: Path to the input forward/R1 FASTQ file.
            tempdir: Path to a temporary directory.
            fastq2: Path to the input reverse/R2 FASTQ file.
            reversed_primers: Whether to treat R1 and R2 in reversed orientation.
        """
        SeqSample.__init__(self, fastq, tempdir)
        if reversed_primers:
            self.r1 = fastq2
            self.fastq2 = fastq
        else:
            self.r1 = fastq
            self.fastq2 = fastq2

    def _merge_reads(self, threads: Union[int, str], stagger: bool) -> None:
        """Merges paired reads using Vsearch.

        Args:
            threads: The number of processor threads to use.
            stagger: Whether to allow merging of staggered reads.

        Raises:
            subprocess.CalledProcessError: If Vsearch or zstd execution fails.
            FileNotFoundError: If Vsearch is not found.
        """
        try:
            seq_file = os.path.join(self.tempdir, 'seq.fq')
            if not os.path.exists(self.tempdir):
                logging.info(f"Expected {self.tempdir} to exist, but it does not. Creating it now.")
                os.makedirs(self.tempdir)

            # Decompress zstd if inputs are in .zst format
            if self.r1 and self.fastq2 and self.r1.endswith('.zst') and self.fastq2.endswith('.zst'):
                r1_temp = os.path.join(self.tempdir, 'r1_temp.fq')
                r2_temp = os.path.join(self.tempdir, 'r2_temp.fq')

                parameters = ['zstd', '-d', self.r1, '-o', r1_temp]
                p1 = subprocess.run(parameters, stderr=subprocess.PIPE)
                p1.check_returncode()
                logging.info(p1.stderr.decode('utf-8'))

                parameters = ['zstd', '-d', self.fastq2, '-o', r2_temp]
                p1 = subprocess.run(parameters, stderr=subprocess.PIPE)
                p1.check_returncode()
                logging.info(p1.stderr.decode('utf-8'))

                self.r1 = r1_temp
                self.fastq2 = r2_temp

            if self.r1 is None or self.fastq2 is None:
                raise ValueError("Both r1 and fastq2 paths must be defined to merge reads.")

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
