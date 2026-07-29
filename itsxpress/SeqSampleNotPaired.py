from itsxpress.SeqSample import SeqSample

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
