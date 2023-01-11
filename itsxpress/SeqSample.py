import os
import logging
import tempfile
import subprocess



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