"""Pipelines.py: Parallel QIIME 2 pipelines for ITSxpress.

"""
from typing import Any

def parallel_trim_single(
    ctx: Any,
    per_sample_sequences: Any,
    region: str,
    taxa: str,
    threads: int = 1,
    cluster_id: float = 1.0,
    trim_ccs: bool = False,
    num_splits: int = None
) -> Any:
    """Runs ITSxpress trim_single across multiple samples in parallel using Parsl.

    Args:
        ctx: QIIME 2 pipeline execution context.
        per_sample_sequences: Input single-end sequence artifact.
        region: The target ITS region to search ("ITS2", "ITS1", or "ALL").
        taxa: Taxonomic group of interest.
        threads: Number of processor threads to use per parallel worker job.
        cluster_id: Dereplication or clustering threshold identity.
        trim_ccs: Whether to enable trim-ccs mode for PacBio CCS reads.

    Returns:
        The combined trimmed single-end sequence artifact.
    """
    split_demux = ctx.get_action('demux', 'split_samples')
    trim_method = ctx.get_action('itsxpress', 'trim_single')
    collate_demux = ctx.get_action('demux', 'collate_samples')

    # SPLIT STAGE
    split_results = split_demux(demux=per_sample_sequences)
    samples_dict = split_results.partitioned_demux

    # APPLY STAGE (Concurrently trims every sample)
    trimmed_futures = {}
    for sample_id, sample_artifact in samples_dict.items():
        try:
            future_result = trim_method(
                per_sample_sequences=sample_artifact,
                region=region,
                taxa=taxa,
                threads=threads,
                cluster_id=cluster_id,
                trim_ccs=trim_ccs
            )
            # Explicitly calling _result() on the output of trim_method within
            # the try-except block to block the main thread and catch errors properly.
            trimmed, = future_result._result()
        except Exception as e:
            raise e
        trimmed_futures[sample_id] = trimmed

    # COMBINE STAGE
    collated = collate_demux(results=trimmed_futures)
    return collated.collated_demux


def parallel_trim_pair(
    ctx: Any,
    per_sample_sequences: Any,
    region: str,
    taxa: str,
    threads: int = 1,
    reversed_primers: bool = False,
    allow_staggered_reads: bool = True,
    cluster_id: float = 1.0
) -> Any:
    """Runs ITSxpress trim_pair across multiple samples in parallel using Parsl.

    Args:
        ctx: QIIME 2 pipeline execution context.
        per_sample_sequences: Input paired-end sequence artifact.
        region: The target ITS region to search ("ITS2", "ITS1", or "ALL").
        taxa: Taxonomic group of interest.
        threads: Number of processor threads to use per parallel worker job.
        reversed_primers: Primers are in reverse orientation.
        allow_staggered_reads: Allow merging staggered reads.
        cluster_id: Dereplication or clustering threshold identity.

    Returns:
        The combined trimmed and merged paired-end sequence artifact.
    """
    split_demux = ctx.get_action('demux', 'split_samples')
    trim_method = ctx.get_action('itsxpress', 'trim_pair')
    collate_demux = ctx.get_action('demux', 'collate_samples')

    # SPLIT STAGE
    split_results = split_demux(demux=per_sample_sequences)
    samples_dict = split_results.partitioned_demux

    # APPLY STAGE
    trimmed_futures = {}
    for sample_id, sample_artifact in samples_dict.items():
        try:
            future_result = trim_method(
                per_sample_sequences=sample_artifact,
                region=region,
                taxa=taxa,
                threads=threads,
                reversed_primers=reversed_primers,
                allow_staggered_reads=allow_staggered_reads,
                cluster_id=cluster_id
            )
            trimmed, = future_result._result()
        except Exception as e:
            raise e
        trimmed_futures[sample_id] = trimmed

    # COMBINE STAGE
    collated = collate_demux(results=trimmed_futures)
    return collated.collated_demux


def parallel_trim_output_unmerged(
    ctx: Any,
    per_sample_sequences: Any,
    region: str,
    taxa: str,
    threads: int = 1,
    reversed_primers: bool = False,
    allow_staggered_reads: bool = True,
    cluster_id: float = 1.0
) -> Any:
    """Runs ITSxpress trim_pair_output_unmerged across multiple samples in parallel using Parsl.

    Args:
        ctx: QIIME 2 pipeline execution context.
        per_sample_sequences: Input paired-end sequence artifact.
        region: The target ITS region to search ("ITS2", "ITS1", or "ALL").
        taxa: Taxonomic group of interest.
        threads: Number of processor threads to use per parallel worker job.
        reversed_primers: Primers are in reverse orientation.
        allow_staggered_reads: Allow merging staggered reads.
        cluster_id: Dereplication or clustering threshold identity.

    Returns:
        The combined trimmed unmerged paired-end sequence artifact.
    """
    split_demux = ctx.get_action('demux', 'split_samples')
    trim_method = ctx.get_action('itsxpress', 'trim_pair_output_unmerged')
    collate_demux = ctx.get_action('demux', 'collate_samples')

    # SPLIT STAGE
    split_results = split_demux(demux=per_sample_sequences)
    samples_dict = split_results.partitioned_demux

    # APPLY STAGE
    trimmed_futures = {}
    for sample_id, sample_artifact in samples_dict.items():
        try:
            future_result = trim_method(
                per_sample_sequences=sample_artifact,
                region=region,
                taxa=taxa,
                threads=threads,
                reversed_primers=reversed_primers,
                allow_staggered_reads=allow_staggered_reads,
                cluster_id=cluster_id
            )
            trimmed, = future_result._result()
        except Exception as e:
            raise e
        trimmed_futures[sample_id] = trimmed

    # COMBINE STAGE
    collated = collate_demux(results=trimmed_futures)
    return collated.collated_demux
