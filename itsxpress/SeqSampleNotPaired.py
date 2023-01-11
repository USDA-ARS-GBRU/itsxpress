from itsxpress.SeqSample import SeqSample

class SeqSampleNotPaired(SeqSample):
    """SeqSample class extended to unpaired format.
    """

    def __init__(self, fastq, tempdir):
        SeqSample.__init__(self, fastq, tempdir)
        self.seq_file = self.fastq
        self.r1 = self.fastq
        self.fastq2 = None