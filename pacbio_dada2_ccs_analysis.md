# Supporting PacBio CCS Reads for DADA2 `denoise-ccs`

This document details the unique requirements of PacBio Circular Consensus Sequencing (CCS) reads, the needs of the DADA2 `--denoise-ccs` pipeline, and the concrete changes required in `ITSxpress` to support this workflow.

---

## 1. PacBio CCS & DADA2 Characteristics

### A. Sequence Orientation
* **The Problem:** Unlike Illumina sequencing where forward and reverse reads have fixed orientations (R1 and R2), PacBio CCS reads represent long single-pass consensus molecules. They are stored in a single FASTQ file but are in a **random mix of forward and reverse-complement orientations**.
* **DADA2 Requirement:** The DADA2 `denoise-ccs` method requires reads to be consistently oriented (5' to 3') before denoising. DADA2's primer-removal step handles this by checking each read against the user-supplied forward primer, and if the reverse complement is a better match, **re-orienting (flipping) the read to the forward strand**.

### B. Length and Variable Amplicon Sizes
* PacBio CCS reads for fungal ITS can cover the entire ITS region (ITS1 + 5.8S + ITS2), which typically ranges from 500bp to 1600bp.
* Trimming must be highly accurate to preserve the full internal spacer diversity while discarding conserved SSU and LSU primer binding sites.

### C. Quality Scores (Highly Accurate Long Reads)
* PacBio CCS reads have exceptionally high quality scores (typically $Q \ge 30$, accuracy $> 99.9\%$).
* Unlike Illumina reads, there is no steep quality drop-off towards the 3' end of the reads. Truncation parameters (like `--p-trunc-len`) should remain set to `0` to avoid discarding long sequences.
* `ITSxpress` must preserve these precise quality scores during trimming.

---

## 2. Technical Requirements for `itsxpress` Integration

To seamlessly process PacBio CCS reads for DADA2, `ITSxpress` needs to perform four operations:

```
[Raw PacBio FASTQ (Mixed Strands)]
               │
               ▼
   [Step 1: Orient all reads]  ───► Flipping Reverse Complement reads to Forward strand
               │
               ▼
   [Step 2: ITS identification] ───► Run HMM search on forward-oriented sequences
               │
               ▼
   [Step 3: Precise Trimming]  ───► Discard Conserved regions, slice ITS coordinates
               │
               ▼
[Output FASTQ (Consistently oriented & trimmed)] ───► Ready for DADA2 `--denoise-ccs`
```

---

## 3. Necessary Changes in `ITSxpress` Codebase

### A. Step 1: Sequence Orientation (Pre-processing)
Because HMM models are strand-specific, oriening reads first makes HMM searches more accurate. We need to implement a pre-processing re-orientation function.

1. **Orientation Logic:**
   - The user provides the forward primer sequence (e.g., `--front` or a CLI option).
   - For each input read, we check whether the forward primer aligns better to the sequence or its reverse complement.
   - If it aligns better to the reverse complement, we replace the sequence and its quality scores with their reverse complement.

```python
from Bio.Seq import Seq

def reorient_read(record, forward_primer: str) -> Tuple[bool, object]:
    """Checks if the record needs to be reverse-complemented to match the forward primer.

    Returns:
        (is_flipped, oriented_record)
    """
    seq_str = str(record.seq)
    rc_seq_str = str(record.seq.reverse_complement())

    # Check simple prefix/substring match or global sequence alignment score
    # (If reverse-complement matches the forward primer better, flip the read)
    if rc_seq_str.find(forward_primer) < seq_str.find(forward_primer):
        flipped = record.reverse_complement(id=True, name=True, description=True)
        return True, flipped
    return False, record
```

### B. Step 2: Preserving the Single-End Pipeline
PacBio CCS reads are single-end reads. `ITSxpress` must route them through `SeqSampleNotPaired` without attempts to perform read merging or paired mapping.
- Ensure that the `--single_end` flag is active.
- Verify that `SeqSampleNotPaired` utilizes the newly re-oriented files for its downstream vsearch clustering and HMMER search operations.

### C. Step 3: Exact Quality Score Preservation
When trimming, we must slice both the sequence and the quality scores cleanly using coordinate slices:
```python
# Biopython slicing automatically preserves the corresponding quality scores (letter_annotations)
trimmed_record = record[start:stop]
```
Ensure that no quality score scaling or modification occurs, as DADA2 error models are highly sensitive to raw quality score distributions.

### D. Step 4: CLI and QIIME2 Plugin Integration

#### 1. CLI Additions
Add options to the main argument parser in `itsxpress/main.py` to support PacBio CCS preparation:
- `--pacbio_ccs`: Flag indicating the input data consists of PacBio CCS reads.
- `--front`: User-supplied forward primer to guide the orientation step.

#### 2. QIIME2 Plugin Changes (`itsxpress/q2_itsxpress.py` & `itsxpress/plugin_setup.py`)
Add a dedicated method or parameter to handle PacBio inputs.

In `itsxpress/plugin_setup.py`, register a parameter to pass primer sequences, or set up a wrapper function:

```python
# In itsxpress/q2_itsxpress.py:
def trim_pacbio_ccs(
    per_sample_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
    region: str,
    front: str,
    taxa: str = "F",
    threads: int = 1
) -> CasavaOneEightSingleLanePerSampleDirFmt:
    # 1. Reorient files in per_sample_sequences using the forward primer (front)
    # 2. Run the single-end trimming pipeline on the oriented sequences
    # 3. Return results format expected by DADA2 denoise-ccs
```

Register this method as an action in `plugin_setup.py` under the QIIME2 registry:
```python
plugin.methods.register_function(
    function=q2_itsxpress.trim_pacbio_ccs,
    inputs={'per_sample_sequences': SampleData[SequencesWithQuality]},
    parameters={
        'region': Str % Choices(['ITS1', 'ITS2', 'ALL']),
        'front': Str,
        'taxa': Str % Choices(taxa_choices),
        'threads': Int
    },
    outputs=[('trimmed_sequences', SampleData[SequencesWithQuality])],
    input_descriptions={'per_sample_sequences': 'The single-end PacBio CCS sequences.'},
    parameter_descriptions={
        'front': 'The forward primer sequence to guide read orientation.',
        'region': 'The ITS region to trim.'
    },
    output_descriptions={'trimmed_sequences': 'The trimmed and oriented sequences.'},
    name='Trim PacBio CCS reads',
    description='Trims and re-orients single-end PacBio CCS reads for downstream DADA2 analysis.'
)
```

### E. Step 5: Consensus Flanking Anchors for DADA2 Compatibility
To satisfy DADA2 `denoise-ccs`'s mandatory requirement for primer sequences without re-introducing non-biological artifacts, `ITSxpress` will utilize a fixed-length trimming window rather than an absolute boundary slice.

1. **Junction Anchor Slicing:**
   Instead of cutting exactly at the coordinates identified by the HMMER search, the trimming function will expand the coordinate boundaries outward to include exactly **15 base pairs** of the conserved flanking structural genes.
   - **5' Edge (Forward Anchor):** Retain the last 15 bp of the upstream gene (e.g., 18S for ITS1, or 5.8S for ITS2).
   - **3' Edge (Reverse Anchor):** Retain the first 15 bp of the downstream gene (e.g., 5.8S for ITS1, or 28S for ITS2).

2. **Consensus Calculation:**
   Because `ITSxpress` pre-orients 100% of the reads to the forward strand in Step 1, these 15 bp flanking regions will be highly conserved and predictable based on the user's selected `--taxa` and `--region`. The plugin will internally compute these two 15 bp strings.

3. **User Output Logging:**
   At the successful completion of the pipeline, the CLI and QIIME 2 visualizer will print the exact 15 bp consensus sequences to the terminal log. This tells the user exactly what to type into their downstream DADA2 command, removing all user guesswork.

```text
Success: ITSxpress completed. Reads are 100% forward-oriented.
To run 'qiime dada2 denoise-ccs' next, use these exact parameters for your taxon:

   --p-front    [15-BP_CALCULATED_FORWARD_CONSENUS]
   --p-adapter  [15-BP_CALCULATED_REVERSE_CONSENUS]
```

### 💻 Required Code Changes for `itsxpress/q2_itsxpress.py`

Update your proposed wrapper registration to output these calculated consensus strings as metadata or via the standard logging channel so the user receives them:

```python
# In itsxpress/q2_itsxpress.py:
def trim_pacbio_ccs(
    per_sample_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
    region: str,
    front: str,
    taxa: str = "F",
    threads: int = 1
) -> CasavaOneEightSingleLanePerSampleDirFmt:
    # 1. Reorient files in per_sample_sequences using the forward primer (front)
    # 2. Run the single-end trimming pipeline, leaving 15 bp flanking anchor boundaries
    # 3. Compute the 15 bp consensus sequences based on the 'taxa' and 'region' matrix
    # 4. Log/Print the calculated --p-front and --p-adapter strings for DADA2 compatibility
    # 5. Return results format expected by DADA2 denoise-ccs
```

***

Would you like to double-check the **exact 15 bp consensus nucleotide sequences** for the fungal junctions (18S, 5.8S, and 28S) to ensure they are hardcoded accurately for the script?

---

## 4. Operational Pipeline Checklist

When executing this pipeline for DADA2 `denoise-ccs`:
1. **Consistently Oriented Output:** Outputs must be 100% forward-oriented (5' to 3' relative to the forward primer).
2. **Conserved Primers Removal:** Conserved primer regions must be completely sliced off, leaving only the variable spacer sequences.
3. **No Quality Distortions:** Phred quality scores must be accurately preserved and untouched.
4. **Valid Format:** The output FASTQ file must be a valid Single-End format (`SampleData[SequencesWithQuality]`) compatible with `qiime dada2 denoise-ccs`.
