# Codebase Analysis: Classes, Methods, Redundancies, and Inconsistencies

This document provides a detailed analysis of the current class structure, methods, redundancies, inconsistencies, and logical bugs within the `ITSxpress` software package.

---

## 1. Class-by-Class Breakdown

### A. `SeqSample` (Base Class)
* **File Location:** `itsxpress/SeqSample.py`
* **Purpose:** Acts as the parent class for processing sequence data into trimmed sequences by interfacing with external binaries like `vsearch` and `hmmsearch`.
* **Key Attributes:**
  - `tempdir`: Temporary directory path.
  - `fastq`: Path to the input FASTQ file.
  - `uc_file`, `rep_file`, `dom_file`, `seq_file`: Paths to intermediate files created during the pipeline.
* **Key Methods:**
  - `deduplicate(threads)`: Executes `vsearch --fastx_uniques` to dereplicate sequences.
  - `cluster(threads, cluster_id)`: Executes `vsearch --cluster_size` to cluster sequences.
  - `_search(hmmfile, threads)`: Executes `hmmsearch --domtblout` to search for ITS profiles.

### B. `SeqSampleNotPaired` (Subclass)
* **File Location:** `itsxpress/SeqSampleNotPaired.py`
* **Purpose:** Represents single-end/unpaired sequence data.
* **Key Attributes:**
  - `seq_file = self.fastq`
  - `r1 = self.fastq`
  - `fastq2 = None`
* **Key Methods:** Inherits all methods from `SeqSample` without modification.

### C. `SeqSamplePairedNotInterleaved` (Subclass)
* **File Location:** `itsxpress/SeqSamplePaired.py`
* **Purpose:** Represents paired-end sequences stored in two separate FASTQ files.
* **Key Attributes:**
  - `r1`, `fastq2`: Holds R1 and R2 files. Swaps them if `reversed_primers=True`.
* **Key Methods:**
  - `_merge_reads(threads, stagger)`: Calls `zstd` to decompress if needed, then executes `vsearch --fastq_mergepairs` to merge R1 and R2 into a temporary `seq.fq` file.

### D. `ItsPosition`
* **File Location:** `itsxpress/ItsPosition.py`
* **Purpose:** Parses HMMER `domtblout` (domain table) outputs to locate the starting and ending coordinates of the ITS region.
* **Key Attributes:**
  - `ddict`: Dictionary holding scores, start and stop positions for sequences.
  - `leftprefix`, `rightprefix`: String prefixes identifying target HMMs for the selected region.
* **Key Methods:**
  - `_score(sequence, stype, score, from_pos, to_pos, tlen)`: Evaluates and updates the highest-bit-score positions.
  - `parse()`: Iteratively reads the HMMER domtable file.
  - `get_position(sequence)`: Returns `(start, stop, tlen)` zero-indexed coordinates.

### E. `Dedup`
* **File Location:** `itsxpress/Dedup.py`
* **Purpose:** Manages the mapping of duplicate reads back to their unique/representative sequences and writes the trimmed sequences.
* **Key Attributes:**
  - `matchdict`: Maps each sequence ID to its representative unique sequence ID.
* **Key Methods:**
  - `parse()`: Parses the `.uc` mapping file from `vsearch`.
  - `_get_paired_seq_generator()`, `_get_trimmed_seq_generator()`: Python generator functions utilizing `Biopython` to yield trimmed FASTQ records.
  - `create_paired_trimmed_seqs()`, `create_trimmed_seqs()`: Reads input files, runs generators, and writes outputs in plain, `.gz`, or `.zst` formats.

---

## 2. Structural and Architectural Inconsistencies

### A. Liskov Substitution Principle (LSP) Violations
The class inheritance hierarchy is structurally unsound:
- The base class `SeqSample` constructor only defines `fastq` and does not define attributes for `r1` or `fastq2`.
- `SeqSampleNotPaired` sets dummy values `r1 = fastq` and `fastq2 = None`.
- `SeqSamplePairedNotInterleaved` defines real paths for `r1` and `fastq2`.
- In `Dedup.py` and `main.py`, the code accesses `sobj.r1` and `sobj.fastq2` polymorphically. If a caller uses the base class `SeqSample` directly, these property accesses will throw `AttributeError`. A subclass should be substitutable for its base class without causing such errors.

### B. Double Responsibility (Configuration vs. Execution State)
The `SeqSample` classes merge configuration details (e.g., input file paths, temp directories) with execution states (e.g., executing commands and mutating internal path variables like `self.seq_file`, `self.rep_file`, `self.uc_file`). This makes state tracking highly fragile and difficult to trace.

### C. Hardcoded File Management in Shared Directories
The filenames inside the temporary directories (`uc.txt`, `rep.fa`, `domtbl.txt`, `seq.fq`) are hardcoded as static string literals across `SeqSample` and its subclasses. If multiple threads, processes, or samples share a single temp directory, they will silently overwrite each other's intermediate files, causing severe race conditions and corrupted runs.

---

## 3. Redundancy and Code Duplication

### A. Vsearch Command Execution
The `deduplicate` and `cluster` methods in `SeqSample.py` duplicate almost all of their shell-out logic, error handling (`subprocess.CalledProcessError`, `FileNotFoundError`), logging, and output assignment block.

### B. Input File Compression and Parsing
`Dedup.py` duplicates file-open/decompression check loops twice—once in `create_trimmed_seqs` and once in `create_paired_trimmed_seqs`. It manually checks suffix types (`.gz`, `.zst`, `.fq`, `.fastq`) and opens them accordingly using separate library packages (`gzip`, `pyzstd`, built-in `open`).
This logic completely bypasses the reusable `read_file` context manager defined in `main.py`, creating redundant and brittle compression/decompression logic across the codebase.

### C. Duplicate Main Execution Path
The file `itsxpress/q2_itsxpress.py` duplicates massive portions of the workflow orchestration logic found in `itsxpress/main.py`. This includes:
- Checking input files.
- Configuring/running vsearch, clustering, and deduplication.
- Initializing `ItsPosition` and `Dedup`.
- Instantiating and managing temporary directories.
Any bug fix, parameter addition, or efficiency improvement made to `main.py` must be manually duplicated in `q2_itsxpress.py`.

---

## 4. Logical Bugs and Critical Issues

### A. The "Coordinate Zero" Falsy Bug (Critical)
In `Dedup.py`, both `_get_paired_seq_generator` and `_get_trimmed_seq_generator` contain the following filtering logic:
```python
if start and stop:
    if start < stop:
        return True
```
In Python, the integer `0` is falsy. If the identified ITS region begins exactly at the first nucleotide (index `0`), `start` is evaluated as `0`. Because `0` is falsy, the condition `if start and stop` evaluates to `False`. As a result, **any sequence where the ITS region starts at the first base is silently discarded.**
* **The Fix:** Use explicit identity/None checks:
```python
if start is not None and stop is not None:
```

### B. Generator Exhaustion Bug under `wri_file == False`
In `Dedup.py`, the parameters `wri_file=False` are designed to check for empty sequence outputs. However, if `wri_file == False`, the method executes:
```python
for i in list(gen1_split_a):
    ...
```
Because `gen1_split_a` is a Python generator, calling `list()` on it **fully consumes and exhausts the generator.** When the method later returns `gen1_split_a`, it is empty, and any downstream operations attempting to write or parse it will receive zero records.

### C. Index Error on Empty FASTQ Check
In `main.py`, the `_check_fastqs` function validates paired names by pulling the first two records using `islice`:
```python
records = SeqIO.parse(handle, 'fastq')
reclist = list(islice(records, 2))
return test_pair_names_str(reclist[0].id, reclist[1].id)
```
If the input FASTQ file is empty, `reclist` will be empty (`[]`). Accessing `reclist[0]` will immediately trigger an unhandled `IndexError`, crashing the entire application instead of failing with a meaningful validation error.

### D. Redundant and Extremely Slow IO Reads
The `_check_total_reads` utility function in `main.py` is called up to four times during a standard run to count sequence reads. It opens the entire FASTQ file and loops line-by-line:
```python
for i, _ in enumerate(handle):
    if i % 4 == 0:
        n += 1
```
For large input files (e.g., multi-gigabyte files), this line-by-line iteration is incredibly slow and duplicates disk I/O, negating the performance benefits gained from high-speed C-based tools like `vsearch`.

---

## 5. Testing and Best Practices Deficiencies

- **Executable Mocking:** Unit tests (`tests/test_main_pytest.py`) directly execute system-level binaries (`vsearch`, `hmmsearch`, `zstd`). If any of these are missing from the runner's path, the entire unit test suite fails, making the test suite highly environment-dependent.
- **Global Logging side effects:** `_logger_setup` in `main.py` modifies the root logger's basic configuration (`logging.basicConfig`), which causes side effects and output pollution when `itsxpress` is imported as a library/plugin inside larger frameworks like QIIME2.
- **Lack of Modern Type Hints:** The codebase is written in an older Python style without PEP 484 type hinting or static analysis support, making code paths difficult to trace.
