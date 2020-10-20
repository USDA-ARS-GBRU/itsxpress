#!/usr/bin/env python
"""A python module to integrate ITSxpress into QIIME for the trimming of amplicon sequences.

Authors: Adam Rivers and Kyle Weber, USDA Agricultural Research Service

The internally transcribed spacer region is a region between the highly conserved small
subunit (SSU) of rRNA and the large subunit (LSU) of the rRNA. The eukaryotic ITS contains
the 5.8s gene and two variable length spacer regions. In amplicon sequencing studies it is
common practice to trim off the conserved (SSU, 5,8S or LSU) regions. Bengtsson-Palme
et al. (2013) published software the software package ITSx to do this.

ITSxpress is a high-speed implementation of the methods in ITSx than also allows FASTQ
files to be processed. Processing FASTQ files Which is essential for analyzing
sequences using the newer exact Sequence Variant methods in Qiime2, Dada2, Deblur
and Unoise that are replacing OTU clustering.

ITSxpress is also available as a stand-alone software package. See
https://github.com/USDA-ARS-GBRU/itsxpress for details.

"""
import os
import pathlib
import shutil
import math

import pandas as pd
from q2_types.per_sample_sequences import (SingleLanePerSamplePairedEndFastqDirFmt,
                                           SingleLanePerSampleSingleEndFastqDirFmt,
                                           CasavaOneEightSingleLanePerSampleDirFmt)
from itsxpress import main as itsxpress
from itsxpress.definitions import (taxa_dict,
                   ROOT_DIR)

default_cluster_id=0.995


def _set_fastqs_and_check(fastq: str,
                          fastq2: str,
                          sample_id: str,
                          single_end: bool,
                          reversed_primers: bool,
                          threads: int) -> object:
    """Checks and writes the fastqs as well as if they are paired end, interleaved and single end.

        Args:
            fastq (str): The path to the forward reads file.
            fastq2 (str): The path to the reverse reads file.
            sample_id (str): The Sample ID.
            single_end (bool): If the sequences are singled ended or not
            threads (int): The amount of threads to use

        Returns:
            (object): The sobj object

        Raises:
            ValueError1: for FASTQ format issue.

        """
    # checking fastqs
    try:
        itsxpress._check_fastqs(fastq=fastq, fastq2=fastq2)
        # Parse input types
        paired_end, interleaved = itsxpress._is_paired(fastq=fastq,
                                                       fastq2=fastq2,
                                                       single_end=single_end)
    except (NotADirectoryError,
            FileNotFoundError):

        raise ValueError("There is a problem with the fastq file(s) you selected")
        # Create SeqSample objects and merge if needed.
    if paired_end and interleaved:
        sobj = itsxpress.SeqSamplePairedInterleaved(fastq=fastq,
                                                    tempdir=None,
                                                    reversed_primers=reversed_primers)
        sobj._merge_reads(threads=threads)
        return sobj

    elif paired_end and not interleaved:
        sobj = itsxpress.SeqSamplePairedNotInterleaved(fastq=fastq,
                                                       fastq2=fastq2,
                                                       tempdir=None,
                                                       reversed_primers=reversed_primers)
        sobj._merge_reads(threads=threads)
        return sobj

    elif not paired_end and not interleaved:
        sobj = itsxpress.SeqSampleNotPaired(fastq=fastq,
                                            tempdir=None)
        return sobj


def _taxa_prefix_to_taxa(taxa_prefix: str) -> str:
    """Turns the taxa prefix letter into the taxa

        Args:
            taxa_prefix (str): The taxa prefix that will be converted to taxa.

        Returns:
            (str): The Taxa

    """
    taxa_dic = {"A": "Alveolata", "B": "Bryophyta", "C": "Bacillariophyta", "D": "Amoebozoa", "E": "Euglenozoa",
                "F": "Fungi", "G": "Chlorophyta", "H": "Rhodophyta", "I": "Phaeophyceae", "L": "Marchantiophyta",
                "M": "Metazoa", "O": "Oomycota", "P": "Haptophyceae", "Q": "Raphidophyceae", "R": "Rhizaria",
                "S": "Synurophyceae", "T": "Tracheophyta", "U": "Eustigmatophyceae", "ALL": "All"}
    taxa_choice = taxa_dic[taxa_prefix]
    return taxa_choice


# First command Trim for SingleLanePerSampleSingleEndFastqDirFmt
def trim_single(per_sample_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
                region: str,
                taxa: str = "F",
                threads: int = 1,
                cluster_id: float = default_cluster_id) -> CasavaOneEightSingleLanePerSampleDirFmt:
    results = main(per_sample_sequences=per_sample_sequences,
                   threads=threads,
                   taxa=taxa,
                   region=region,
                   paired_in=False,
                   paired_out=False,
                   reversed_primers=False,
                   cluster_id=cluster_id)
    return results


# Second command Trim for SingleLanePerSamplePairedEndFastqDirFmt
def trim_pair(per_sample_sequences: SingleLanePerSamplePairedEndFastqDirFmt,
              region: str,
              taxa: str = "F",
              threads: int = 1,
              reversed_primers: bool = False,
              cluster_id: float = default_cluster_id) -> CasavaOneEightSingleLanePerSampleDirFmt:
    results = main(per_sample_sequences=per_sample_sequences,
                   threads=threads,
                   taxa=taxa,
                   region=region,
                   paired_in=True,
                   paired_out=False,
                   reversed_primers=reversed_primers,
                   cluster_id=cluster_id)
    return results

# Second command Trim for SingleLanePerSamplePairedEndFastqDirFmt
def trim_pair_output_unmerged(per_sample_sequences: SingleLanePerSamplePairedEndFastqDirFmt,
              region: str,
              taxa: str = "F",
              threads: int = 1,
              reversed_primers: bool = False,
              cluster_id: float = default_cluster_id) -> CasavaOneEightSingleLanePerSampleDirFmt:
    results = main(per_sample_sequences=per_sample_sequences,
                   threads=threads,
                   taxa=taxa,
                   region=region,
                   paired_in=True,
                   paired_out=True,
                   reversed_primers=reversed_primers,
                   cluster_id=cluster_id)
    return results
# The ITSxpress handling
def main(per_sample_sequences,
         threads: int,
         taxa: str,
         region: str,
         paired_in: bool,
         paired_out: bool,
         reversed_primers: bool,
         cluster_id: float) -> CasavaOneEightSingleLanePerSampleDirFmt:
    """The main communication between the plugin and the ITSxpress program.

    Args:
        per_sample_sequences (SingleLanePerSampleSingleEndFastqDirFmt or SingleLanePerSamplePairedEndFastqDirFmt):
        the input sequences.
        threads (int) : The number of threads to use.
        taxa (str): The taxa to be used for the search.
        region (str) : The region to be used for the search.
        paired_in (bool): Declares if input files are paired.
        paired_out (bool): Declares if output files should be paired.
        cluster_id (float):The percent identity for clustering reads, set to 1 for exact dereplication.

    Returns:
        (CasavaOneEightSingleLanePerSampleDirFmt): A catch-all output type for
        both single and paired-end reads.

    Raises:
        ValueError1: hmmsearch error.

    """
    # Setting the taxa
    taxa = _taxa_prefix_to_taxa(taxa)
    samples = per_sample_sequences.manifest.view(pd.DataFrame)
    # Creating result dir
    results = CasavaOneEightSingleLanePerSampleDirFmt()
    # Running the for loop for each sample
    for sample in samples.itertuples():
        # writing fastqs and their attributes and checking the files
        sobj = _set_fastqs_and_check(
            fastq=sample.forward,
            fastq2=sample.reverse if paired_in else None,
            sample_id=sample.Index,
            single_end=paired_out,
            reversed_primers=reversed_primers,
            threads=threads)
        # Deduplicate
        if math.isclose(cluster_id, 1,rel_tol=1e-05):
            sobj.deduplicate(threads=threads)
        else:
            sobj.cluster(threads=threads, cluster_id=cluster_id)
        try:
            # HMMSearch for ITS regions
            hmmfile = os.path.join(ROOT_DIR, "ITSx_db", "HMMs", taxa_dict[taxa])
            sobj._search(hmmfile=hmmfile, threads=threads)
        except (ModuleNotFoundError,
                FileNotFoundError,
                NotADirectoryError):

            raise ValueError("hmmsearch was not found, make sure HMMER3 is installed and executable")

        # Parse HMMseach output.
        its_pos = itsxpress.ItsPosition(domtable=sobj.dom_file,
                                        region=region)
        # Create deduplication object.
        dedup_obj = itsxpress.Dedup(uc_file=sobj.uc_file,
                                    rep_file=sobj.rep_file,
                                    seq_file=sobj.seq_file,
                                    fastq=sobj.r1,
                                    fastq2=sobj.fastq2)

        # Copy the original filename, that way we preserve all filename fields.
        out_path_fwd = os.path.join(str(results),
                                    pathlib.Path(sample.forward).name)

        # Create trimmed sequences.
        if paired_out:
            # Copy the original filename, that way we preserve all filename fields.
            out_path_rev = os.path.join(str(results),
                                        pathlib.Path(sample.reverse).name)
            dedup_obj.create_paired_trimmed_seqs(out_path_fwd,
                                                 out_path_rev,
                                                 gzipped=True,
                                                 itspos=its_pos)
        else:
            dedup_obj.create_trimmed_seqs(out_path_fwd,
                                      gzipped=True,
                                      itspos=its_pos)
        # Deleting the temp files.
        shutil.rmtree(sobj.tempdir)
    # Writing out the results.
    return results
