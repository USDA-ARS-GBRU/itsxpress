#!/usr/bin/env python
"""A python module to integrate ITSxpress into QIIME for the trimming of amplicon sequences.

Authors: Adam Rivers, Kyle Weber, Sveinn Einarsson USDA Agricultural Research Service

The internally transcribed spacer region is a region between the highly conserved small
subunit (SSU) of rRNA and the large subunit (LSU) of the rRNA. The eukaryotic ITS contains
the 5.8s gene and two variable length spacer regions. In amplicon sequencing studies it is
common practice to trim off the conserved (SSU, 5,8S or LSU) regions. Bengtsson-Palme
et al. (2013) published software the software package ITSx to do this.

ITSxpress is a high-speed implementation of the methods in ITSx than also allows FASTQ
files to be processed. Processing FASTQ files Which is essential for analyzing
sequences using the newer exact Sequence Variant methods in Qiime2, Dada2, Deblur
and Unoise that are replacing OTU clustering.
"""

import os
import pathlib
import shutil
import math
import tempfile
from typing import Any

import pandas as pd
from q2_types.per_sample_sequences import (SingleLanePerSamplePairedEndFastqDirFmt,
                                           SingleLanePerSampleSingleEndFastqDirFmt,
                                           CasavaOneEightSingleLanePerSampleDirFmt)
from itsxpress import main as itsxpress
from itsxpress.definitions import taxa_dict, ROOT_DIR

default_cluster_id: float = 1.0


def _set_fastqs_and_check(fastq: str,
                          fastq2: str,
                          tempdir: str,
                          sample_id: str,
                          single_end: bool,
                          reversed_primers: bool,
                          allow_staggered_reads: bool,
                          threads: int) -> Any:
    """Checks and writes the fastqs as well as if they are paired end and single end.

    Args:
        fastq: The path to the forward reads file.
        fastq2: The path to the reverse reads file.
        tempdir: Path to a temporary directory.
        sample_id: The Sample ID.
        single_end: If the sequences are singled ended or not.
        reversed_primers: Primers are in reverse orientation.
        allow_staggered_reads: Allow merging staggered reads.
        threads: The amount of threads to use.

    Returns:
        The sequence sample object.

    Raises:
        ValueError: for FASTQ format issue.
    """
    try:
        itsxpress._check_fastqs(fastq=fastq, fastq2=fastq2)
        paired_end = itsxpress._is_paired(fastq=fastq,
                                          fastq2=fastq2,
                                          single_end=single_end)
    except (NotADirectoryError, FileNotFoundError):
        raise ValueError("There is a problem with the fastq file(s) you selected")

    if paired_end:
        sobj = itsxpress.SeqSamplePairedNotInterleaved(fastq=fastq,
                                                       fastq2=fastq2,
                                                       tempdir=tempdir,
                                                       reversed_primers=reversed_primers)
        sobj._merge_reads(threads=threads, stagger=allow_staggered_reads)
        return sobj
    else:
        sobj = itsxpress.SeqSampleNotPaired(fastq=fastq,
                                            tempdir=tempdir)
        return sobj


def _taxa_prefix_to_taxa(taxa_prefix: str) -> str:
    """Turns the taxa prefix letter into the taxa name.

    Args:
        taxa_prefix: The taxa prefix character.

    Returns:
        The taxonomical group name.
    """
    taxa_dic = {"A": "Alveolata", "B": "Bryophyta", "C": "Bacillariophyta", "D": "Amoebozoa", "E": "Euglenozoa",
                "F": "Fungi", "G": "Chlorophyta", "H": "Rhodophyta", "I": "Phaeophyceae", "L": "Marchantiophyta",
                "M": "Metazoa", "O": "Oomycota", "P": "Haptophyceae", "Q": "Raphidophyceae", "R": "Rhizaria",
                "S": "Synurophyceae", "T": "Tracheophyta", "U": "Eustigmatophyceae", "Y": "Parabasalia", "ALL": "All"}
    return taxa_dic[taxa_prefix]


def trim_single(per_sample_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
                region: str,
                taxa: str = "F",
                threads: int = 1,
                cluster_id: float = default_cluster_id,
                trim_ccs: bool = False) -> CasavaOneEightSingleLanePerSampleDirFmt:
    """Trims single-end sequence reads.

    Args:
        per_sample_sequences: Input single-end sequences.
        region: The target ITS region to trim.
        taxa: Taxonomic group of interest.
        threads: Number of threads to use.
        cluster_id: Dereplication or clustering threshold identity.
        trim_ccs: Whether trim-ccs mode is enabled for PacBio reads.

    Returns:
        The Casava 1.8 single-end directory output.
    """
    results = main(per_sample_sequences=per_sample_sequences,
                   threads=threads,
                   taxa=taxa,
                   region=region,
                   paired_in=False,
                   paired_out=False,
                   reversed_primers=False,
                   allow_staggered_reads=False,
                   cluster_id=cluster_id,
                   trim_ccs=trim_ccs)
    return results


def trim_pair(per_sample_sequences: SingleLanePerSamplePairedEndFastqDirFmt,
              region: str,
              taxa: str = "F",
              threads: int = 1,
              reversed_primers: bool = False,
              allow_staggered_reads: bool = True,
              cluster_id: float = default_cluster_id) -> CasavaOneEightSingleLanePerSampleDirFmt:
    """Trims paired-end sequence reads.

    Args:
        per_sample_sequences: Input paired-end sequences.
        region: The target ITS region to trim.
        taxa: Taxonomic group of interest.
        threads: Number of threads to use.
        reversed_primers: Primers are in reverse orientation.
        allow_staggered_reads: Allow merging staggered reads.
        cluster_id: Dereplication or clustering threshold identity.

    Returns:
        The Casava 1.8 directory output.
    """
    results = main(per_sample_sequences=per_sample_sequences,
                   threads=threads,
                   taxa=taxa,
                   region=region,
                   paired_in=True,
                   paired_out=False,
                   reversed_primers=reversed_primers,
                   allow_staggered_reads=allow_staggered_reads,
                   cluster_id=cluster_id,
                   trim_ccs=False)
    return results


def trim_pair_output_unmerged(per_sample_sequences: SingleLanePerSamplePairedEndFastqDirFmt,
                              region: str,
                              taxa: str = "F",
                              threads: int = 1,
                              reversed_primers: bool = False,
                              allow_staggered_reads: bool = True,
                              cluster_id: float = default_cluster_id) -> CasavaOneEightSingleLanePerSampleDirFmt:
    """Trims paired-end sequence reads, keeping them unmerged.

    Args:
        per_sample_sequences: Input paired-end sequences.
        region: The target ITS region to trim.
        taxa: Taxonomic group of interest.
        threads: Number of threads to use.
        reversed_primers: Primers are in reverse orientation.
        allow_staggered_reads: Allow merging staggered reads.
        cluster_id: Dereplication or clustering threshold identity.

    Returns:
        The Casava 1.8 directory output.
    """
    results = main(per_sample_sequences=per_sample_sequences,
                   threads=threads,
                   taxa=taxa,
                   region=region,
                   paired_in=True,
                   paired_out=True,
                   reversed_primers=reversed_primers,
                   allow_staggered_reads=allow_staggered_reads,
                   cluster_id=cluster_id,
                   trim_ccs=False)
    return results


def main(per_sample_sequences: Any,
         threads: int,
         taxa: str,
         region: str,
         paired_in: bool,
         paired_out: bool,
         reversed_primers: bool,
         allow_staggered_reads: bool,
         cluster_id: float,
         trim_ccs: bool = False) -> CasavaOneEightSingleLanePerSampleDirFmt:
    """The main communication between the plugin and the ITSxpress program.

    Args:
        per_sample_sequences: The input sequences.
        threads: The number of threads to use.
        taxa: The taxonomic group used for the search.
        region: The region used for the search.
        paired_in: Declares if input files are paired.
        paired_out: Declares if output files should be paired.
        reversed_primers: Primers are in reverse orientation.
        allow_staggered_reads: Allows merging of staggered reads.
        cluster_id: The percent identity for clustering reads.
        trim_ccs: Whether the trim-ccs mode is enabled.

    Returns:
        A directory output type for fastq files.

    Raises:
        ValueError: hmmsearch error or folder creation error.
    """
    taxa = _taxa_prefix_to_taxa(taxa)
    samples = per_sample_sequences.manifest.view(pd.DataFrame)
    try:
        tempdir = tempfile.mkdtemp(prefix='itsxpress_')
    except Exception:
        raise ValueError("Could not create temporary directory")

    results = CasavaOneEightSingleLanePerSampleDirFmt()

    for sample in samples.itertuples():
        sobj = _set_fastqs_and_check(
            fastq=sample.forward,
            fastq2=sample.reverse if paired_in else None,
            tempdir=tempdir,
            sample_id=sample.Index,
            single_end=False if paired_in else True,
            reversed_primers=reversed_primers,
            allow_staggered_reads=allow_staggered_reads,
            threads=threads
        )
        if trim_ccs:
            sobj.orient_reads(threads=threads)

        if math.isclose(cluster_id, 1, rel_tol=1e-05):
            sobj.deduplicate(threads=threads)
        else:
            sobj.cluster(threads=threads, cluster_id=cluster_id)

        try:
            from itsxpress.main import create_runtime_hmm
            hmmfile = create_runtime_hmm(taxa, region, tempdir)
            sobj._search(hmmfile=hmmfile, threads=threads)
        except (ModuleNotFoundError, FileNotFoundError, NotADirectoryError):
            raise ValueError("hmmsearch was not found, make sure HMMER3 is installed and executable")

        its_pos = itsxpress.ItsPosition(domtable=sobj.dom_file, region=region)
        dedup_obj = itsxpress.Dedup(uc_file=sobj.uc_file,
                                    rep_file=sobj.rep_file,
                                    seq_file=sobj.seq_file,
                                    fastq=sobj.r1,
                                    fastq2=sobj.fastq2)

        out_path_fwd = os.path.join(str(results), pathlib.Path(sample.forward).name)

        if paired_out:
            out_path_rev = os.path.join(str(results), pathlib.Path(sample.reverse).name)
            dedup_obj.create_paired_trimmed_seqs(out_path_fwd,
                                                 out_path_rev,
                                                 gzipped=True,
                                                 zstd_file=False,
                                                 itspos=its_pos,
                                                 wri_file=True,
                                                 trim_ccs=trim_ccs)
        else:
            dedup_obj.create_trimmed_seqs(out_path_fwd,
                                          gzipped=True,
                                          zstd_file=False,
                                          itspos=its_pos,
                                          wri_file=True,
                                          tempdir=sobj.tempdir,
                                          trim_ccs=trim_ccs)

    if trim_ccs:
        import logging
        parent_pid = os.getppid()
        marker_path = os.path.join(tempfile.gettempdir(), f"itsxpress_trim_ccs_{parent_pid}.tmp")
        if not os.path.exists(marker_path):
            try:
                with open(marker_path, 'w') as f:
                    f.write('printed')
            except Exception:
                pass
            msg = (
                "\n" + "=" * 80 + "\n"
                "PacBio CCS trimming complete. You can denoise these sequences in DADA2 using:\n\n"
                "qiime dada2 denoise-ccs \\\n"
                "  --i-demultiplexed-seqs pacbio-trimmed.qza \\\n"
                "  --p-front GACAGGTACAAGAAGGA \\\n"
                "  --p-adapter ACTGGAGACTGGGTTAA \\\n"
                "  --p-min-len 50 \\\n"
                "  --p-max-len 1600 \\\n"
                "  --p-n-threads 4 \\\n"
                "  --output-dir dada2-results\n"
                "*" * 80 + "\n"
            )
            print(msg)
            logging.info(msg)

    shutil.rmtree(tempdir)
    return results
