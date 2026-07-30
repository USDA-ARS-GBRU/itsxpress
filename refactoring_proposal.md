# Refactoring Proposal: Architectural Redesign for ITSxpress

This document proposes a comprehensive refactoring and architectural redesign for `ITSxpress`. It aims to resolve the class redundancies, logical flaws, validation inconsistencies, and platform coupling identified in the codebase analysis.

---

## 1. Architectural Vision

We propose moving away from an inheritance-based model to a **composition-based architecture** driven by dataclasses, isolated pipeline executors, and unified file/stream managers.

```
┌────────────────────────────────────────────────────────┐
│                   CLI / QIIME2 Plugin                  │
└───────────────────────────┬────────────────────────────┘
                            │
                            ▼
┌────────────────────────────────────────────────────────┐
│                   PipelineCoordinator                  │
└───────────────────────────┬────────────────────────────┘
                            │
                            ├────────────────────────────┐
                            ▼                            ▼
┌──────────────────────────────────────┐      ┌─────────────────────────────┐
│          IO & Format Manager         │      │     DependencyExecutor      │
│  (Unified gz, zst, raw file parser)   │      │  (Wrapper for CLI binaries) │
└──────────────────────────────────────┘      └─────────────────────────────┘
```

---

## 2. Redesigned Class & Module Hierarchy

### A. `itsxpress.models` (Domain Objects)
Instead of subclassing sequence containers, we introduce clean immutable state representations:

```python
from dataclasses import dataclass
from typing import Optional, Tuple

@dataclass(frozen=True)
class SequenceInput:
    r1_path: str
    r2_path: Optional[str] = None
    single_end: bool = False
    reversed_primers: bool = False

    @property
    def is_paired(self) -> bool:
        return self.r2_path is not None and not self.single_end
```

### B. `itsxpress.executor` (Execution Wrappers)
Abstract all command-line operations into an isolated executor class. This separates system commands from domain/sequence representations.

```python
import subprocess
import logging
import shutil
from typing import List

class DependencyExecutor:
    """Manages system paths, execution, and verification of external binary dependencies."""

    def __init__(self, vsearch_path: str = "vsearch", hmmer_path: str = "hmmsearch"):
        self.vsearch = vsearch_path
        self.hmmer = hmmer_path
        self._verify_dependencies()

    def _verify_dependencies(self):
        for cmd in [self.vsearch, self.hmmer]:
            if not shutil.which(cmd):
                raise FileNotFoundError(f"Required dependency '{cmd}' is not executable or in PATH.")

    def run_vsearch(self, args: List[str]) -> subprocess.CompletedProcess:
        try:
            return subprocess.run([self.vsearch] + args, capture_output=True, check=True, text=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"vsearch execution failed: {e.stderr}")
            raise e
```

### C. `itsxpress.trimmer` (Core Domain Logic)
A unified trimmer engine that relies on composition, receiving parsed positions and outputs.

```python
class ITSTrimmer:
    """Coordinates trimming coordinates on BioPython record streams."""

    def __init__(self, start: int, stop: int, seq_len: int):
        # Explicit identity checks prevent the falsy zero index bug
        self.start = start if start is not None else 0
        self.stop = stop if stop is not None else seq_len

    def trim_forward(self, record):
        return record[self.start:self.stop]

    def trim_reverse(self, record):
        # Calculate matching coordinates for R2
        r2_start = seq_len - self.stop
        r2_end = seq_len - self.start
        return record[r2_start:r2_end]
```

### D. `itsxpress.pipeline` (Orchestration Coordinator)
A single, high-level coordinator that manages temp files, sequence parsing, and calls the `DependencyExecutor` and `ITSTrimmer`. Both the command-line interface (`main.py`) and QIIME2 setup (`q2_itsxpress.py`) will import and call this coordinator, eliminating duplicate execution paths.

---

## 3. Advantages of the Refactoring Proposal

| Area | Refactored Architecture | Advantage |
| :--- | :--- | :--- |
| **Code Redundancy** | Replaces duplicated pipeline scripts in `main.py` and `q2_itsxpress.py` with a single reusable `PipelineCoordinator`. | Single source of truth. Maintenance, updates, or feature additions only need to be written once. |
| **Robustness** | Safe coordinates checking (`is not None`) instead of truthy evaluations (`if start and stop`). | Fixes the critical "Coordinate Zero" bug where valid ITS sequences starting at base 0 were silently dropped. |
| **Concurrency / Thread Safety** | Dynamic, unique filenames based on hashes or UUIDs inside the working directories. | Safe multi-threaded or multi-sample parallel execution without race conditions or files overwriting each other. |
| **Separation of Concerns** | Decoupling execution states, dependency calling, and output formatting. | Makes code much easier to mock and test. No environment-dependent system failures during testing. |
| **I/O Efficiency** | Utilizing shared file buffers and context managers (`read_file`). | Drastically reduces disk read-write amplification and prevents double-reading files. |

---

## 4. Disadvantages and Risks of the Proposal

1. **Backwards Compatibility (QIIME2 Interoperability):**
   - QIIME2 is historically highly sensitive to changes in method signatures and dependency namespaces. Refactoring internal class models could break external plugin bindings if not carefully bridged.
   - *Mitigation:* Retain the legacy function signatures in `q2_itsxpress.py` but internally map them to the new `PipelineCoordinator` instead of using the custom duplicated pipeline loops.

2. **Integration & Refactoring Effort:**
   - Redesigning the core class hierarchy requires rewriting some existing unit tests.
   - *Mitigation:* Implement the refactoring incrementally, keeping integration tests intact to verify that output FASTA/FASTQ sequence hashes match the old version.

---

## 5. Improved Testing Strategy

To eliminate brittle tests dependent on external binary presence, we propose:

1. **Unit Tests with Mocking:**
   - Use Python's `unittest.mock` to mock `subprocess.run` calls, returning simulated `.uc` mapping files and `domtbl.txt` logs. This allows running the test suite on any system without `vsearch` or `hmmer` pre-installed.
2. **Regression Integration Tests:**
   - Provide a set of locked sample input FASTQ files.
   - Execute the complete pipeline on these files and assert that the trimmed sequence hashes match expected checksums.
3. **Empty File / Malformed Header Assertions:**
   - Expand the test suite to verify graceful failures (with detailed log statements) instead of unhandled `IndexError` crashes on empty files or unmatched read headers.
