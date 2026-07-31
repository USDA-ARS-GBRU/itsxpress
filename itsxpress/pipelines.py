"""Pipelines.py: Parallel QIIME 2 pipelines for ITSxpress.

"""
try:
    from qiime2.plugin import IContext
except ImportError:
    class IContext:
        pass

try:
    from qiime2.sdk import Artifact
except ImportError:
    class Artifact:
        pass


def trim_single(
    ctx: IContext,
    per_sample_sequences: Artifact,
    region: str,
    taxa: str,
    threads: int = 1,
    cluster_id: float = 1.0,
    trim_ccs: bool = False,
    num_splits: int = None
) -> Artifact:
    """Runs ITSxpress trim_single_sample across multiple samples in parallel using Parsl.

    Args:
        ctx: QIIME 2 pipeline execution context.
        per_sample_sequences: Input single-end sequence artifact.
        region: The target ITS region to search ("ITS2", "ITS1", or "ALL").
        taxa: Taxonomic group of interest.
        threads: Number of processor threads to use per parallel worker job.
        cluster_id: Dereplication or clustering threshold identity.
        trim_ccs: Whether to enable trim-ccs mode for PacBio CCS reads.
        num_splits: Unused.

    Returns:
        The combined trimmed single-end sequence artifact.
    """
    split_action = ctx.get_action('itsxpress', 'split_single_end')
    trim_method = ctx.get_action('itsxpress', 'trim_single_sample')
    combine_action = ctx.get_action('itsxpress', 'combine_single')

    # SPLIT STAGE
    splits, = split_action(per_sample_sequences=per_sample_sequences)

    # APPLY STAGE (Concurrently trims every sample)
    trimmed_futures = {}
    for sample_id, sample_artifact in splits.items():
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
    combined_trimmed, = combine_action(results=trimmed_futures)
    return combined_trimmed


def trim_pair(
    ctx: IContext,
    per_sample_sequences: Artifact,
    region: str,
    taxa: str,
    threads: int = 1,
    reversed_primers: bool = False,
    allow_staggered_reads: bool = True,
    cluster_id: float = 1.0,
    num_splits: int = None
) -> Artifact:
    """Runs ITSxpress trim_pair_sample across multiple samples in parallel using Parsl.

    Args:
        ctx: QIIME 2 pipeline execution context.
        per_sample_sequences: Input paired-end sequence artifact.
        region: The target ITS region to search ("ITS2", "ITS1", or "ALL").
        taxa: Taxonomic group of interest.
        threads: Number of processor threads to use per parallel worker job.
        reversed_primers: Primers are in reverse orientation.
        allow_staggered_reads: Allow merging staggered reads.
        cluster_id: Dereplication or clustering threshold identity.
        num_splits: Unused.

    Returns:
        The combined trimmed and merged paired-end sequence artifact.
    """
    split_action = ctx.get_action('itsxpress', 'split_paired_end')
    trim_method = ctx.get_action('itsxpress', 'trim_pair_sample')
    combine_action = ctx.get_action('itsxpress', 'combine_pair')

    # SPLIT STAGE
    splits, = split_action(per_sample_sequences=per_sample_sequences)

    # APPLY STAGE
    trimmed_futures = {}
    for sample_id, sample_artifact in splits.items():
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
    combined_trimmed, = combine_action(results=trimmed_futures)
    return combined_trimmed


def trim_pair_output_unmerged(
    ctx: IContext,
    per_sample_sequences: Artifact,
    region: str,
    taxa: str,
    threads: int = 1,
    reversed_primers: bool = False,
    allow_staggered_reads: bool = True,
    cluster_id: float = 1.0,
    num_splits: int = None
) -> Artifact:
    """Runs ITSxpress trim_pair_sample_unmerged across multiple samples in parallel using Parsl.

    Args:
        ctx: QIIME 2 pipeline execution context.
        per_sample_sequences: Input paired-end sequence artifact.
        region: The target ITS region to search ("ITS2", "ITS1", or "ALL").
        taxa: Taxonomic group of interest.
        threads: Number of processor threads to use per parallel worker job.
        reversed_primers: Primers are in reverse orientation.
        allow_staggered_reads: Allow merging staggered reads.
        cluster_id: Dereplication or clustering threshold identity.
        num_splits: Unused.

    Returns:
        The combined trimmed unmerged paired-end sequence artifact.
    """
    split_action = ctx.get_action('itsxpress', 'split_paired_end')
    trim_method = ctx.get_action('itsxpress', 'trim_pair_sample_unmerged')
    combine_action = ctx.get_action('itsxpress', 'combine_pair_unmerged')

    # SPLIT STAGE
    splits, = split_action(per_sample_sequences=per_sample_sequences)

    # APPLY STAGE
    trimmed_futures = {}
    for sample_id, sample_artifact in splits.items():
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
    combined_trimmed, = combine_action(results=trimmed_futures)
    return combined_trimmed
