import os
import subprocess
import logging
import gzip
from itertools import tee
from typing import Dict, Tuple, Optional, Any, Generator, Iterator, Union, List

import pyzstd as zstd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from itsxpress.definitions import ROOT_DIR, taxa_choices, taxa_dict, maxmismatches, maxratio, vsearch_fastq_qmax

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

    def orient_reads(self, threads: Union[int, str] = 1) -> None:
        """Orients reads using Vsearch --orient against the universal reference database.

        Args:
            threads: The number of processor threads to use.

        Raises:
            subprocess.CalledProcessError: If Vsearch execution fails.
            FileNotFoundError: If Vsearch binary cannot be located.
        """
        try:
            orient_ref = os.path.join(ROOT_DIR, "universal_orient_ref_clean.fasta.gz")
            oriented_fastq = os.path.join(self.tempdir, 'oriented.fq')
            parameters = [
                "vsearch",
                "--orient", self.fastq,
                "--db", orient_ref,
                "--fastqout", oriented_fastq,
                "--id", "0.35",
                "--threads", str(threads)
            ]
            p = subprocess.run(parameters, stderr=subprocess.PIPE)
            p.check_returncode()
            logging.info(p.stderr.decode('utf-8'))

            # Update the sequence files to point to the oriented fastq
            self.fastq = oriented_fastq
            self.seq_file = oriented_fastq
            self.r1 = oriented_fastq
        except subprocess.CalledProcessError as e:
            logging.exception("Could not orient reads with Vsearch. Error from Vsearch was:\n {}".format(p.stderr.decode('utf-8')))
            raise e
        except FileNotFoundError as f:
            logging.error("Vsearch was not found, make sure Vsearch is installed and executable")
            raise f

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


class SeqSampleNotPaired(SeqSample):
    """SeqSample class extended to unpaired format."""

    def __init__(self, fastq: str, tempdir: str) -> None:
        """Initializes a SeqSampleNotPaired object for single-end reads.

        Args:
            fastq: Path to the input single-end FASTQ file.
            tempdir: Path to a temporary directory.
        """
        SeqSample.__init__(self, fastq, tempdir)
        self.seq_file = self.fastq
        self.r1 = self.fastq
        self.fastq2 = None


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


class ItsPosition:
    """Class for ITS positional information derived from hmmsearch domtable files.

    Attributes:
        domtable (str): Path to the HMMER3 hmmsearch domtable output file.
        region (str): The region of the ITS to extract. Must be "ITS1", "ITS2", or "ALL".
        ddict (Dict[str, Any]): A dictionary holding the scores and start/stop positions.
            Example: {sample: {left: {score: 31, to_pos: 15, from_pos: 0}, right: {score: 32, to_pos: 354, from_pos: 200}, tlen: 449}}
        leftprefix (str): Prefix of the HMM profiles corresponding to the left boundary.
        rightprefix (str): Prefix of the HMM profiles corresponding to the right boundary.
    """

    def __init__(self, domtable: str, region: str) -> None:
        """Initializes ItsPosition, parses the domtable file, and determines boundaries.

        Args:
            domtable: Path to the HMMER3 hmmsearch domtable output file.
            region: The region of the ITS to extract (choices: "ITS1", "ITS2", "ALL").
        """
        self.domtable: str = domtable
        self.ddict: Dict[str, Any] = {}
        if region == "ITS2":
            self.leftprefix: str = '3_'
            self.rightprefix: str = '4_'
        elif region == "ITS1":
            self.leftprefix = '1_'
            self.rightprefix = '2_'
        elif region == "ALL":
            self.leftprefix = '1_'
            self.rightprefix = '4_'
        self.parse()

    def _score(self, sequence: str, stype: str, score: float, from_pos: int, to_pos: int, tlen: int) -> None:
        """Evaluates scores and positions from a domtable line and updates ddict.

        Args:
            sequence: The name/ID of the sequence.
            stype: Either 'left' or 'right' indicating the boundary type.
            score: The bit score of the match from hmmsearch.
            from_pos: The beginning coordinate of the HMM match.
            to_pos: The ending coordinate of the HMM match.
            tlen: The total length of the sequence.
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

    def parse(self) -> None:
        """Parses the domtable file generated by hmmsearch.

        Identifies matches for left and right HMM boundaries and updates `ddict`
        with the highest-scoring boundaries.

        Raises:
            Exception: If an error occurs when parsing HMMER hmmsearch results.
        """
        try:
            with open(self.domtable, 'r') as f:
                for line in f:
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
            logging.error("Exception occurred when parsing HMMSearch results")
            raise e

    def get_position(self, sequence: str) -> Tuple[Optional[int], Optional[int], Optional[int]]:
        """Returns the start and stop positions for a given sequence.

        Args:
            sequence: The name/ID of the sequence.

        Returns:
            A tuple of (start coordinate, stop coordinate, total sequence length).
            Coordinates are 0-indexed and can be None if not found.

        Raises:
            KeyError: If the input sequence was not parsed or has no identified boundaries.
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
            return (start, stop, tlen)
        except KeyError:
            logging.debug("No ITS stop or start sites were identified for sequence {}, skipping.".format(sequence))
            raise KeyError


class Dedup:
    """A class to handle deduplicated sequence data.

    To speed processing, vsearch is used to remove duplicate amplicons so that the
    start and stop sites are determined only once.

    Attributes:
        matchdict (Optional[Dict[str, str]]): A dictionary of each sequence ID as a key and
            its representative sequence ID as a value {seq1: rep1, seq2: rep1, seq3: rep2}.
        uc_file (str): The location of the .uc file containing matching information.
        rep_file (str): The location of the representative sequences file.
        seq_file (str): The location of the complete sequences file.
        fastq (Optional[str]): The location of the input FASTQ file.
        fastq2 (Optional[str]): The location of the optional Read 2 input FASTQ if paired.
    """

    def __init__(self, uc_file: str, rep_file: str, seq_file: str, fastq: Optional[str] = None, fastq2: Optional[str] = None) -> None:
        """Initializes the Dedup object and parses the uc file.

        Args:
            uc_file: The path of the .uc mapping file generated by Vsearch.
            rep_file: The path of the representative FASTA file generated by Vsearch.
            seq_file: The path of the complete sequences file.
            fastq: Optional path to the input FASTQ/R1 file.
            fastq2: Optional path to the input FASTQ/R2 file.
        """
        self.matchdict: Optional[Dict[str, str]] = None
        self.uc_file: str = uc_file
        self.rep_file: str = rep_file
        self.seq_file: str = seq_file
        self.fastq: Optional[str] = fastq
        self.fastq2: Optional[str] = fastq2
        self.parse()

    def parse(self) -> None:
        """Parse the uc data file to populate the matchdict attribute.

        Raises:
            Exception: General exception if the uc file is not parsed properly.
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

    def _get_paired_seq_generator(self, zipseqgen: Iterator[Tuple[SeqRecord, SeqRecord]], itspos: ItsPosition, wri_file: bool, trim_ccs: bool = False) -> Tuple[Generator[SeqRecord, None, None], Generator[SeqRecord, None, None]]:
        """This function takes a zipped object of two Biopython SeqIO sequence generators, and

        returns two generators of Biopython SeqRecords for Dada2. Sequences where the ITS ends could
        not be determined are omitted.

        Args:
            zipseqgen: A zipped iterator of two Biopython SeqIO generators
                for the forward and reverse input sequences.
            itspos: An itsxpress ItsPosition object.
            wri_file: Whether the sequences will be written to files.
            trim_ccs: Whether the trim_ccs mode is enabled (attaches synthetic primers).

        Returns:
            A tuple of two python generators yielding forward and reverse trimmed SeqRecords.
        """
        def _filterfunc(ziprecord: Tuple[SeqRecord, SeqRecord]) -> bool:
            """Filters records down to those that contain a valid ITS start and stop position."""
            try:
                record1, record2 = ziprecord
                if self.matchdict is not None and record1.id in self.matchdict:
                    repseq = self.matchdict[record1.id]
                    start, stop, tlen = itspos.get_position(repseq)
                    if start is not None and stop is not None:
                        if start < stop:
                            return True
                return False
            except KeyError:
                return False

        def _stitch_primers(record: SeqRecord) -> SeqRecord:
            fwd_seq = "GACAGGTACAAGAAGGA"
            rev_seq = "TTAACCCAGTCTCCAGT"

            new_seq = Seq(fwd_seq) + record.seq + Seq(rev_seq)
            fwd_qual = [93] * len(fwd_seq)
            rev_qual = [93] * len(rev_seq)

            if "phred_quality" in record.letter_annotations:
                new_qual = fwd_qual + list(record.letter_annotations["phred_quality"]) + rev_qual
            else:
                new_qual = fwd_qual + [93] * len(record) + rev_qual

            new_record = SeqRecord(
                new_seq,
                id=record.id,
                name=record.name,
                description=record.description
            )
            new_record.letter_annotations["phred_quality"] = new_qual
            return new_record

        def _map_func(ziprecord: Tuple[SeqRecord, SeqRecord]) -> Tuple[SeqRecord, SeqRecord]:
            """Trims the record down to the selected ITS region."""
            record1, record2 = ziprecord
            if self.matchdict is None:
                raise ValueError("matchdict must be parsed before calling sequence generators.")
            repseq = self.matchdict[record1.id]
            start, stop, tlen = itspos.get_position(repseq)
            if start is None or stop is None or tlen is None:
                raise ValueError(f"Could not retrieve valid positions for sequence {record1.id}")
            r2start = tlen - stop
            r2end = tlen - start #calculate end of R2
            try:
                if stop > tlen:
                    record1_return = record1[start:]
                elif stop <= tlen:
                    record1_return = record1[start:stop]
                else:
                    raise ValueError("An error occurred when trimming the forward read of {}".format(record1.id))
                if r2end > tlen:
                    record2_return = record2[r2start:]
                elif r2end <= tlen:
                    record2_return = record2[r2start:r2end]
                else:
                    raise ValueError("An error occurred when trimming the reverse read of {}".format(record2.id))
            except ValueError as e:
                logging.exception(e)
                raise e

            if trim_ccs:
                record1_return = _stitch_primers(record1_return)
                record2_return = _stitch_primers(record2_return)

            return record1_return, record2_return

        def _split_gen(gen: Iterator[Tuple[SeqRecord, SeqRecord]]) -> Tuple[Generator[SeqRecord, None, None], Generator[SeqRecord, None, None]]:
            gen_a, gen_b = tee(gen, 2)
            return (a for a, b in gen_a), (b for a, b in gen_b)

        filt = filter(_filterfunc, zipseqgen)
        filt_a, filt_b = tee(filt, 2)

        gen1_a = map(_map_func, filt_a)
        gen1_b = map(_map_func, filt_b)

        gen1_split_a = (a for a, b in gen1_a)
        gen1_split_b = (b for a, b in gen1_b)

        if not wri_file:
            check_gen_a, gen1_split_a = tee(gen1_split_a, 2)
            check_gen_b, gen1_split_b = tee(gen1_split_b, 2)

            zeroseqctr1 = 0
            seqlist1: List[str] = []
            for i in check_gen_a:
                if i.seq == "":
                    zeroseqctr1 += 1
                    seqlist1.append(i.id)
            zeroseqctr2 = 0
            seqlist2: List[str] = []
            for i in check_gen_b:
                if i.seq == "":
                    zeroseqctr2 += 1
                    seqlist2.append(i.id)
            if zeroseqctr1 != 0 or zeroseqctr2 != 0:
                print("Total number of sequences that are empty Split A: ", zeroseqctr1)
                print("Sequence IDs: ")
                print(seqlist1)
                print("Total number of sequences that are empty Split B: ", zeroseqctr2)
                print("Sequence IDs: ")
                print(seqlist2)

        return gen1_split_a, gen1_split_b

    def create_paired_trimmed_seqs(self, outfile1: str, outfile2: str, gzipped: bool, zstd_file: bool, itspos: ItsPosition, wri_file: bool, trim_ccs: bool = False) -> None:
        """Writes two FASTQ files, optionally compressed, with the reads trimmed to the selected region.

        Args:
            outfile1: The file to write the forward sequences to.
            outfile2: The file to write the reverse sequences to.
            gzipped: Should the output files be gzipped?
            zstd_file: Should the output files be zstd compressed? (.zst files)
            itspos: an ItsPosition object.
            wri_file: Should file be written or checked for empty sequences?
            trim_ccs: Whether the trim_ccs mode is enabled.

        Raises:
            ValueError: If file suffix combinations are mixed or unsupported.
        """
        if self.fastq is None or self.fastq2 is None:
            raise ValueError("Both fastq and fastq2 paths must be defined to create paired trimmed sequences.")

        def _write_seqs(seqs: Generator[SeqRecord, None, None], outfile: str) -> None:
            """Helper function to write sequences in compressed or plain format."""
            if gzipped:
                with gzip.open(outfile, 'wt') as g:
                    SeqIO.write(seqs, g, "fastq")
            elif zstd_file:
                with zstd.open(outfile, 'wt') as g:
                    SeqIO.write(seqs, g, "fastq")
            else:
                with open(outfile, 'w') as g:
                    SeqIO.write(seqs, g, "fastq")

        def _create_gen(f: Any, g: Any) -> None:
            """Create sequence generators and write them."""
            seqgen1 = SeqIO.parse(f, 'fastq')
            seqgen2 = SeqIO.parse(g, 'fastq')
            zipseqgen = zip(seqgen1, seqgen2)
            seqs1, seqs2 = self._get_paired_seq_generator(zipseqgen, itspos, wri_file, trim_ccs=trim_ccs)
            if wri_file:
                _write_seqs(seqs1, outfile1)
                _write_seqs(seqs2, outfile2)

        try:
            if self.fastq.endswith(".gz") and self.fastq2.endswith(".gz"):
                with gzip.open(self.fastq, 'rt') as f:
                    with gzip.open(self.fastq2, 'rt') as g:
                        _create_gen(f, g)
            elif self.fastq.endswith(".zst") and self.fastq2.endswith(".zst"):
                with zstd.open(self.fastq, 'rt') as f:
                    with zstd.open(self.fastq2, 'rt') as g:
                        _create_gen(f, g)
            elif self.fastq.endswith(".fq") and self.fastq2.endswith(".fq"):
                with open(self.fastq, 'r') as f:
                    with open(self.fastq2, 'r') as g:
                        _create_gen(f, g)
            elif (self.fastq.endswith(".fastq") or self.fastq.endswith(".fq")) and (self.fastq2.endswith(".fastq") or self.fastq2.endswith(".fq")):
                with open(self.fastq, 'r') as f:
                    with open(self.fastq2, 'r') as g:
                        _create_gen(f, g)
            else:
                raise ValueError("Fastq and Fastq2 files should both be gzipped (.gz), zstd compressed (.zst) or both be uncompressed. Mixed input is not accepted.")
        except Exception as e:
            raise e

    def _get_trimmed_seq_generator(self, seqgen: Iterator[SeqRecord], itspos: ItsPosition, wri_file: bool, trim_ccs: bool = False) -> Generator[SeqRecord, None, None]:
        """This function takes a Biopython SeqIO sequence generator, and

        returns a generator of trimmed sequences suitable for Deblur. Sequences where the ITS ends could
        not be determined are omitted.

        Args:
            seqgen: A Biopython SeqIO generator of all input sequences.
            itspos: An itsxpress ItsPosition object.
            wri_file: Whether the output will be written.
            trim_ccs: Whether the trim_ccs mode is enabled.

        Returns:
            A python generator yielding trimmed SeqRecords.
        """
        def _filterfunc(record: SeqRecord) -> bool:
            """Filters records down to those that contain a valid ITS start and stop position."""
            try:
                if self.matchdict is not None and record.id in self.matchdict:
                    repseq = self.matchdict[record.id]
                    start, stop, tlen = itspos.get_position(repseq)
                    if start is not None and stop is not None:
                        if start < stop:
                            return True
                return False
            except KeyError:
                return False

        def _stitch_primers(record: SeqRecord) -> SeqRecord:
            fwd_seq = "GACAGGTACAAGAAGGA"
            rev_seq = "TTAACCCAGTCTCCAGT"

            new_seq = Seq(fwd_seq) + record.seq + Seq(rev_seq)
            fwd_qual = [93] * len(fwd_seq)
            rev_qual = [93] * len(rev_seq)

            if "phred_quality" in record.letter_annotations:
                new_qual = fwd_qual + list(record.letter_annotations["phred_quality"]) + rev_qual
            else:
                new_qual = fwd_qual + [93] * len(record) + rev_qual

            new_record = SeqRecord(
                new_seq,
                id=record.id,
                name=record.name,
                description=record.description
            )
            new_record.letter_annotations["phred_quality"] = new_qual
            return new_record

        def map_func(record: SeqRecord) -> SeqRecord:
            """Trims the record down to the selected ITS region."""
            if self.matchdict is None:
                raise ValueError("matchdict must be parsed before calling sequence generators.")
            repseq = self.matchdict[record.id]
            start, stop, tlen = itspos.get_position(repseq)
            if start is None or stop is None:
                raise ValueError(f"Could not retrieve valid positions for sequence {record.id}")
            trimmed = record[start:stop]
            if trim_ccs:
                trimmed = _stitch_primers(trimmed)
            return trimmed

        filt = filter(_filterfunc, seqgen)
        filt_a, filt_b = tee(filt, 2)
        r1 = map(map_func, filt_a)

        if not wri_file:
            check_r1, r1 = tee(r1, 2)
            zeroseqctr = 0
            seqlist: List[str] = []
            for i in check_r1:
                if i.seq == "":
                    zeroseqctr += 1
                    seqlist.append(i.id)
            if zeroseqctr != 0:
                print("Total number of sequences that are empty: ", zeroseqctr)
                print("Sequence IDs: ")
                print(seqlist)

        return (map_func(rec) for rec in filt_b)

    def create_trimmed_seqs(self, outfile: str, gzipped: bool, zstd_file: bool, itspos: ItsPosition, wri_file: bool, tempdir: str, trim_ccs: bool = False) -> None:
        """Creates a FASTQ file, optionally compressed, with the reads trimmed to the selected region.

        Args:
            outfile: The file to write the sequences to.
            gzipped: Should the files be gzipped?
            zstd_file: Should the files be zstd compressed? (.zst files)
            itspos: An ItsPosition object.
            wri_file: Should file be written or checked for empty sequences?
            tempdir: Path to a temporary directory.
            trim_ccs: Whether the trim_ccs mode is enabled.
        """
        def _write_seqs(seqs: Generator[SeqRecord, None, None]) -> None:
            if gzipped:
                tempf = os.path.join(tempdir, 'temp.fa')
                with open(tempf, 'w') as g:
                    SeqIO.write(seqs, g, "fastq")
                with open(tempf, 'rb') as f_in, gzip.open(outfile, 'wb') as f_out:
                    f_out.writelines(f_in)
            elif zstd_file:
                tempf = os.path.join(tempdir, 'temp.fa')
                with open(tempf, 'w') as g:
                    SeqIO.write(seqs, g, "fastq")
                with open(tempf, 'rb') as f_in:
                    with zstd.open(outfile, 'wb') as f_out:
                        f_out.writelines(f_in)
            else:
                with open(outfile, 'w') as g:
                    SeqIO.write(seqs, g, "fastq")

        if self.seq_file.endswith(".gz"):
            with gzip.open(self.seq_file, 'rt') as f:
                seqgen = SeqIO.parse(f, 'fastq')
                seqs = self._get_trimmed_seq_generator(seqgen, itspos, wri_file, trim_ccs=trim_ccs)
                if wri_file:
                    _write_seqs(seqs)
        elif self.seq_file.endswith(".zst"):
            with zstd.open(self.seq_file, 'rt') as f:
                seqgen = SeqIO.parse(f, 'fastq')
                seqs = self._get_trimmed_seq_generator(seqgen, itspos, wri_file, trim_ccs=trim_ccs)
                if wri_file:
                    _write_seqs(seqs)
        else:
            with open(self.seq_file, 'r') as f:
                seqgen = SeqIO.parse(f, 'fastq')
                seqs = self._get_trimmed_seq_generator(seqgen, itspos, wri_file, trim_ccs=trim_ccs)
                if wri_file:
                    _write_seqs(seqs)
