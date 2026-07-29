import os
import logging
import subprocess
from typing import Optional, Union, Any

logger = logging.getLogger(__name__)

class SeqSample:
    """The base class for processing sequence data into trimmed sequences.

    Attributes:
        tempdir (str): Path to a temporary directory.
        fastq (str): The path to the input FASTQ file.
        uc_file (Optional[str]): The path to the Vsearch uc mapping file.
        rep_file (Optional[str]): The path to the representative sequences FASTA file.
        dom_file (Optional[str]): The path to the HMMER domtbl file.
        seq_file (Optional[str]): The path to the FASTQ/gzipped/zst file used for analysis.
        r1 (Optional[str]): Path to the forward/R1 sequence file. Defined here to prevent LSP violations.
        fastq2 (Optional[str]): Path to the reverse/R2 sequence file. Defined here to prevent LSP violations.
    """

    def __init__(self, fastq: str, tempdir: str) -> None:
        """Initializes a SeqSample object.

        Args:
            fastq: The path to the input FASTQ file.
            tempdir: Path to a temporary directory.
        """
        self.tempdir: str = tempdir
        self.fastq: str = fastq
        self.uc_file: Optional[str] = None
        self.rep_file: Optional[str] = None
        self.dom_file: Optional[str] = None
        self.seq_file: Optional[str] = None
        self.r1: Optional[str] = None
        self.fastq2: Optional[str] = None

    def deduplicate(self, threads: Union[int, str] = 1) -> None:
        """Runs Vsearch dereplication to create a FASTA file of non-redundant sequences.

        Args:
            threads: The number of processor threads to use.

        Raises:
            subprocess.CalledProcessError: If the Vsearch process fails.
            FileNotFoundError: If the Vsearch binary cannot be located.
        """
        try:
            self.uc_file = os.path.join(self.tempdir, 'uc.txt')
            self.rep_file = os.path.join(self.tempdir, 'rep.fa')
            parameters = [
                "vsearch",
                "--fastx_uniques",
                self.seq_file,
                "--fastaout", self.rep_file,
                "--uc", self.uc_file,
                "--strand", "both"
            ]
            p2 = subprocess.run(parameters, stderr=subprocess.PIPE)
            logging.info(p2.stderr.decode('utf-8'))
            p2.check_returncode()
        except subprocess.CalledProcessError as e:
            logging.exception(
                "Could not perform dereplication with Vsearch. Error from Vsearch was:\n {}"
                .format(p2.stderr.decode('utf-8'))
            )
            raise e
        except FileNotFoundError as f:
            logging.error("Vsearch was not found, make sure Vsearch is installed and executable")
            raise f

    def cluster(self, threads: Union[int, str], cluster_id: float = 0.995) -> None:
        """Runs Vsearch clustering to create a FASTA file of non-redundant sequences.

        Args:
            threads: The number of processor threads to use.
            cluster_id: The percent identity threshold for clustering.

        Raises:
            subprocess.CalledProcessError: If the Vsearch process fails.
            FileNotFoundError: If the Vsearch binary cannot be located.
        """
        try:
            self.uc_file = os.path.join(self.tempdir, 'uc.txt')
            self.rep_file = os.path.join(self.tempdir, 'rep.fa')
            parameters = [
                "vsearch",
                "--cluster_size", self.seq_file,
                "--centroids", self.rep_file,
                "--uc", self.uc_file,
                "--strand", "both",
                "--id", str(cluster_id),
                "--threads", str(threads)
            ]
            p2 = subprocess.run(parameters, stderr=subprocess.PIPE)
            print(p2.stderr.decode('utf-8'))
            p2.check_returncode()
        except subprocess.CalledProcessError as e:
            logging.exception(
                "Could not perform clustering with Vsearch. Error from Vsearch was:\n {}"
                .format(p2.stderr.decode('utf-8'))
            )
            raise e
        except FileNotFoundError as f:
            logging.error("Vsearch was not found, make sure Vsearch is installed and executable")
            raise f

    def _search(self, hmmfile: str, threads: Union[int, str]) -> None:
        """Searches for ITS regions using hmmsearch.

        Args:
            hmmfile: Path to the target HMM profile file.
            threads: Number of processor threads to use.

        Raises:
            subprocess.CalledProcessError: If the hmmsearch process fails.
            FileNotFoundError: If the hmmsearch binary cannot be located.
        """
        try:
            self.dom_file = os.path.join(self.tempdir, 'domtbl.txt')
            parameters = [
                "hmmsearch",
                "--domtblout", self.dom_file,
                "-T", "10",
                "--cpu", str(threads),
                "--tformat", "fasta",
                "--F1", "1e-6",
                "--F2", "1e-6",
                "--F3", "1e-6",
                hmmfile,
                self.rep_file
            ]
            p4 = subprocess.run(parameters, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
            p4.check_returncode()
        except subprocess.CalledProcessError as e:
            logging.exception(
                "Could not perform ITS identification with hmmsearch. The error was:\n {}"
                .format(p4.stderr.decode('utf-8'))
            )
            raise e
        except FileNotFoundError as f:
            logging.error("hmmsearch was not found, make sure HMMER3 is installed and executable")
            raise f
