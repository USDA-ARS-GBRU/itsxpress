# -*- coding: utf-8 -*-
import os
import tempfile
import shutil
import gzip
import pyzstd as zstd
import subprocess
import filecmp
import sys

from Bio import SeqIO
import pytest

import itsxpress

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
from itsxpress.definitions import ROOT_DIR, taxa_dict
hmmfile = os.path.join(ROOT_DIR,"ITSx_db","HMMs", taxa_dict["Fungi"])


def test_check_fastqs():
	fastq = os.path.join(TEST_DIR, "test_data", "4774-1-MSITS3_R1.fastq")
	fastq2 = os.path.join(TEST_DIR, "test_data", "broken.fastq")
	pytest.raises(ValueError, itsxpress.main._check_fastqs, fastq, fastq2)

def test_check_fastq_gzs():
	fastq = os.path.join(TEST_DIR,"test_data", "4774-1-MSITS3_R1.fastq.gz")
	fastq2 = os.path.join(TEST_DIR,"test_data", "broken.fastq.gz")
	pytest.raises(ValueError, itsxpress.main._check_fastqs, fastq, fastq2)

def test_its_position_init():
	itspos = itsxpress.main.ItsPosition(os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "domtbl.txt"), "ITS2")
	exp1 = {'tlen': 341, 'right': {'score': 59.1, 'to_pos': 326, 'from_pos': 282}, 'left': {'score': 52.2, 'to_pos': 128, 'from_pos': 84}}
	print(itspos.ddict["M02696:28:000000000-ATWK5:1:1101:19331:3209"])
	assert (exp1 == itspos.ddict["M02696:28:000000000-ATWK5:1:1101:19331:3209"])
	exp2 = {'right': {'score': 34.0, 'to_pos': 370, 'from_pos': 327}, 'tlen': 385}
	print(itspos.ddict["M02696:28:000000000-ATWK5:1:1101:23011:4341"])
	assert (exp2 == itspos.ddict["M02696:28:000000000-ATWK5:1:1101:23011:4341"])
	assert (len(itspos.ddict) == 137)

def test_dedup():
	uc = os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "uc.txt")
	seq = os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "seq.fq.gz")
	rep = os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "rep.fa")
	dedup = itsxpress.main.Dedup( uc_file=uc, rep_file=rep, seq_file=seq)
	# Check length of records
	assert len(dedup.matchdict) == 227
	# Check that non-representative seqs are logged
	assert dedup.matchdict['M02696:28:000000000-ATWK5:1:1101:11740:1800'] == 'M02696:28:000000000-ATWK5:1:1101:10899:1561'
	# Check that representative seqs are logged
	assert dedup.matchdict["M02696:28:000000000-ATWK5:1:1101:23011:4341"] == 'M02696:28:000000000-ATWK5:1:1101:23011:4341'

def test_dedup_create_trimmed_seqs():
	tf = tempfile.mkdtemp()
	print(tf)
	uc = os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "uc.txt")
	seq = os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "seq.fq.gz")
	rep = os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "rep.fa")
	dedup = itsxpress.main.Dedup( uc_file=uc, rep_file=rep, seq_file=seq)
	itspos = itsxpress.main.ItsPosition(os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "domtbl.txt"), "ITS2")
	# Check non gzipped
	dedup.create_trimmed_seqs(os.path.join(tf,"testout.fastq"), gzipped=False, zstd_file=False, wri_file=True, itspos=itspos,tempdir=tf)
	with open(os.path.join(tf,"testout.fastq"), 'r') as f:
		recs = SeqIO.parse(f, "fastq")
		n = 0
		length = 0
		for rec in recs:
			n += 1
			length += len(rec)
	print("n: {}, length: {}".format(n, length))
	assert n == 226
	assert length == 42637
	shutil.rmtree(tf)

def test_dedup_create_trimmed_seqs_gzipped():
	tf = tempfile.mkdtemp()
	print(tf)
	uc = os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "uc.txt")
	seq = os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "seq.fq.gz")
	rep = os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "rep.fa")
	dedup = itsxpress.main.Dedup( uc_file=uc, rep_file=rep, seq_file=seq)
	itspos = itsxpress.main.ItsPosition(os.path.join(TEST_DIR, "test_data", "ex_tmpdir", "domtbl.txt"), "ITS2")
	# Check gzipped
	dedup.create_trimmed_seqs(os.path.join(tf,"testout.fastq.gz"), gzipped=True, zstd_file=False, wri_file=True, itspos=itspos,tempdir=tf)
	with gzip.open(os.path.join(tf,"testout.fastq.gz"), 'rt') as f:
		recs = SeqIO.parse(f, "fastq")
		n = 0
		length = 0
		for rec in recs:
			n += 1
			length += len(rec)
	print("n: {}, length: {}".format(n,length))
	assert n == 226
	assert length == 42637
	shutil.rmtree(tf)

def test_dedup_create_trimmed_seqs_zst():
	tf = tempfile.mkdtemp()
	print(tf)
	uc = os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "uc.txt")
	seq = os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "seq.fq.gz")
	rep = os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "rep.fa")
	dedup = itsxpress.main.Dedup( uc_file=uc, rep_file=rep, seq_file=seq)
	itspos = itsxpress.main.ItsPosition(os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "domtbl.txt"), "ITS2")
	# Check zstd compression
	dedup.create_trimmed_seqs(os.path.join(tf,"testout.fastq.zst"), gzipped=False, zstd_file=True, wri_file=True, itspos=itspos,tempdir=tf)
	with zstd.open(os.path.join(tf,"testout.fastq.zst"), 'rt') as f:
		recs = SeqIO.parse(f, "fastq")
		n = 0
		length = 0
		for rec in recs:
			n += 1
			length += len(rec)
	print("n: {}, length: {}".format(n,length))
	assert n == 226
	assert length == 42637
	shutil.rmtree(tf)

#Following test will fail if Vsearch version is less than 2.20 or hmmer is not installed
def test_seq_sample_not_paired():
	tf = tempfile.mkdtemp()
	fastq = os.path.join(TEST_DIR,"test_data", "4774-1-MSITS3_merged.fastq")
	sobj = itsxpress.main.SeqSampleNotPaired(fastq=fastq, tempdir=tf)
	sobj.deduplicate(threads=1)
	sobj._search(hmmfile=hmmfile, threads=1)
	shutil.rmtree(tf)

def test_seq_sample_not_paired_clustered():
	tf = tempfile.mkdtemp()
	fastq = os.path.join(TEST_DIR,"test_data", "4774-1-MSITS3_merged.fastq")
	sobj = itsxpress.main.SeqSampleNotPaired(fastq=fastq, tempdir=tf)
	sobj.cluster(threads=1, cluster_id=0.995)
	sobj._search(hmmfile=hmmfile, threads=1)
	shutil.rmtree(tf)

def test_seq_sample_paired_not_interleaved():
	tf = tempfile.mkdtemp()
	fastq = os.path.join(TEST_DIR, "test_data", "4774-1-MSITS3_R1.fastq")
	fastq2 = os.path.join(TEST_DIR, "test_data", "4774-1-MSITS3_R2.fastq")
	sobj = itsxpress.main.SeqSamplePairedNotInterleaved(fastq=fastq, tempdir=tf, fastq2=fastq2)
	sobj._merge_reads(stagger = True, threads=1)
	sobj.deduplicate(threads=1)
	sobj._search(hmmfile=hmmfile, threads=1)
	shutil.rmtree(tf)


def test_is_paired():
	paired_end= itsxpress.main._is_paired("fastq1.fq", "fastq2.fq", single_end=False)
	assert paired_end == True

	paired_end= itsxpress.main._is_paired("fastq1.fq", None, single_end=False)
	assert paired_end == False

	paired_end= itsxpress.main._is_paired("fastq1.fq", None, single_end=True)
	assert paired_end == False


def test_myparser():
	parser = itsxpress.main.myparser()
	args = parser.parse_args(['--fastq', 'test.fastq','--outfile', 'test.out','--tempdir', 'dirt','--region','ITS1','--taxa', 'Fungi'])
	assert (args.fastq == 'test.fastq')


def test_main_paired():
	parser = itsxpress.main.myparser()
	tf = tempfile.mkdtemp()
	fastq = os.path.join(TEST_DIR, "test_data", "4774-1-MSITS3_R1.fastq")
	fastq2 = os.path.join(TEST_DIR, "test_data", "4774-1-MSITS3_R2.fastq")
	outfile = os.path.join(tf, 'testout.fastq')
	validation = os.path.join(TEST_DIR, "test_data", "testout.fastq")
	args = parser.parse_args(['--fastq', fastq, '--fastq2', fastq2, '--outfile', outfile, '--region','ITS2', '--taxa',  'Fungi', '--threads', '1'])
	itsxpress.main.main(args=args)
	seqs = SeqIO.parse(outfile, 'fastq')
	n = sum(1 for _ in seqs)
	assert (n==235)
	shutil.rmtree(tf)

def test_main_paired_high_qual():
	""" Test Qual scores over 41
	"""
	parser = itsxpress.main.myparser()
	tf = tempfile.mkdtemp()
	fastq = os.path.join(TEST_DIR, "test_data", "high_qual_scores_R1.fastq.gz")
	fastq2 = os.path.join(TEST_DIR, "test_data", "high_qual_scores_R2.fastq.gz")
	outfile = os.path.join(tf, 'testout.fastq')
	validation = os.path.join(TEST_DIR,"test_data", "testout.fastq")
	args = parser.parse_args(['--fastq', fastq, '--fastq2', fastq2, '--outfile', outfile, '--region','ITS2', '--taxa',  'Fungi', '--threads', '1'])
	itsxpress.main.main(args=args)
	seqs = SeqIO.parse(outfile, 'fastq')
	n = sum(1 for _ in seqs)
	print(n)
	#assert (n==235)
	shutil.rmtree(tf)

def test_main_merged():
	parser = itsxpress.main.myparser()
	tf = tempfile.mkdtemp()
	fastq = os.path.join(TEST_DIR,"test_data", "4774-1-MSITS3_merged.fastq")
	outfile = os.path.join(tf,'testout.fastq')
	validation = os.path.join(TEST_DIR,"test_data", "testout.fastq")
	args = parser.parse_args(['--fastq', fastq, '--single_end', '--outfile', outfile, '--region','ITS2', '--taxa', 'Fungi', '--threads', '1'])
	itsxpress.main.main(args=args)
	seqs = SeqIO.parse(outfile, 'fastq')
	n = sum(1 for _ in seqs)
	print(n)
	assert (n == 226)
	shutil.rmtree(tf)

def test_main_paired_no_cluster():
	parser = itsxpress.main.myparser()
	tf = tempfile.mkdtemp()
	fastq = os.path.join(TEST_DIR, "test_data", "4774-1-MSITS3_R1.fastq")
	fastq2 = os.path.join(TEST_DIR,"test_data", "4774-1-MSITS3_R2.fastq")
	outfile = os.path.join(tf,'testout.fastq')
	validation = os.path.join(TEST_DIR, "test_data", "testout.fastq")
	args = parser.parse_args(['--fastq', fastq, '--fastq2', fastq2, '--outfile', outfile, '--region','ITS2', '--taxa',  'Fungi', '--cluster_id', '1', '--threads', '1'])
	itsxpress.main.main(args=args)
	seqs = SeqIO.parse(outfile, 'fastq')
	n = sum(1 for _ in seqs)
	#assert (n==227)
	assert (n==235)
	shutil.rmtree(tf)


def test_get_paired_seq_generator():
	uc = os.path.join(TEST_DIR,  "test_data", "ex_tmpdir", "uc.txt")
	seq = os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "seq.fq.gz")
	rep = os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "rep.fa")
	fastq = os.path.join(TEST_DIR,"test_data", "4774-1-MSITS3_R1.fastq")
	fastq2 = os.path.join(TEST_DIR, "test_data", "4774-1-MSITS3_R2.fastq")
	dedup = itsxpress.main.Dedup( uc_file=uc, rep_file=rep, seq_file=seq, fastq=fastq, fastq2=fastq2)
	itspos = itsxpress.main.ItsPosition(os.path.join(TEST_DIR, "test_data", "ex_tmpdir", "domtbl.txt"), "ITS2")
	f = open(fastq, 'r')
	g = open(fastq2, 'r')
	seqgen1 = SeqIO.parse(f, 'fastq')
	seqgen2 = SeqIO.parse(g, 'fastq')
	zipseqgen = zip(seqgen1, seqgen2)
	seqs1, seqs2  = dedup._get_paired_seq_generator(zipseqgen, itspos, wri_file = True)
	n1 = 0
	n2 = 0
	for rec in seqs1:
		n1 += 1
	for rec in seqs2:
		n2 += 1
	assert (n1==226)
	assert (n2==226)

def test_create_paired_trimmed_seqs():
	uc = os.path.join(TEST_DIR, "test_data", "ex_tmpdir", "uc.txt")
	seq = os.path.join(TEST_DIR, "test_data", "ex_tmpdir", "seq.fq.gz")
	rep = os.path.join(TEST_DIR, "test_data", "ex_tmpdir", "rep.fa")
	fastq = os.path.join(TEST_DIR, "test_data", "4774-1-MSITS3_R1.fastq")
	fastq2 = os.path.join(TEST_DIR, "test_data", "4774-1-MSITS3_R2.fastq")
	dedup = itsxpress.main.Dedup(uc_file=uc, rep_file=rep, seq_file=seq, fastq=fastq, fastq2=fastq2)
	itspos = itsxpress.main.ItsPosition(os.path.join(TEST_DIR,"test_data", "ex_tmpdir", "domtbl.txt"), "ITS2")
	tf = tempfile.mkdtemp()
	print(tf)
	t1 = os.path.join(tf,'t2_r1.fq')
	t2 = os.path.join(tf,'t2_r2.fq')
	dedup.create_paired_trimmed_seqs(t1, t2, False,False, itspos,True)
	assert (filecmp.cmp(t1, os.path.join(TEST_DIR, "test_data", "t2_r1.fq")))
	assert (filecmp.cmp(t2, os.path.join(TEST_DIR, "test_data", "t2_r2.fq")))
	shutil.rmtree(tf)


# --- Additional Unit Tests to improve coverage and verify fixed logical bugs ---

def test_check_fastqs_empty(tmp_path):
	"""Verifies that check_fastqs handles empty files gracefully or returns correct False check on empty streams."""
	empty_file = tmp_path / "empty.fastq"
	empty_file.write_text("")
	# Since it is empty, records list is empty and should return False inside core()
	# Let's ensure calling _check_fastqs does not raise IndexError
	itsxpress.main._check_fastqs(str(empty_file))

def test_coordinate_zero_not_falsy():
	"""Verifies that coordinate 0 (zero-index) is treated as a valid non-None value."""
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord

	class MockItsPosition:
		def get_position(self, seq_id):
			# Coordinate zero must be recognized!
			return (0, 10, 100)

	itspos = MockItsPosition()

	# Create a dummy Dedup object
	uc = os.path.join(TEST_DIR, "test_data", "ex_tmpdir", "uc.txt")
	seq = os.path.join(TEST_DIR, "test_data", "ex_tmpdir", "seq.fq.gz")
	rep = os.path.join(TEST_DIR, "test_data", "ex_tmpdir", "rep.fa")
	dedup = itsxpress.main.Dedup(uc_file=uc, rep_file=rep, seq_file=seq)
	dedup.matchdict = {"seq1": "seq1"}

	# Verify start=0 does not trigger falsy skip in generator filter logic
	seqgen = [SeqRecord(Seq("ATCG" * 50), id="seq1", description="")]
	res_gen = list(dedup._get_trimmed_seq_generator(iter(seqgen), itspos, wri_file=True))

	# Since start=0 is evaluated correctly as a valid position, sequence must be successfully trimmed and returned.
	assert len(res_gen) == 1
	assert str(res_gen[0].seq) == "ATCGATCGAT"

def test_generator_exhaustion_wri_file_false():
	"""Verifies that calling _get_trimmed_seq_generator with wri_file=False does not exhaust the returned generator."""
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord

	class MockItsPosition:
		def get_position(self, seq_id):
			return (5, 15, 100)

	itspos = MockItsPosition()

	uc = os.path.join(TEST_DIR, "test_data", "ex_tmpdir", "uc.txt")
	seq = os.path.join(TEST_DIR, "test_data", "ex_tmpdir", "seq.fq.gz")
	rep = os.path.join(TEST_DIR, "test_data", "ex_tmpdir", "rep.fa")
	dedup = itsxpress.main.Dedup(uc_file=uc, rep_file=rep, seq_file=seq)
	dedup.matchdict = {"seq1": "seq1"}

	seqgen = [SeqRecord(Seq("ATCG" * 50), id="seq1", description="")]
	# Get generator with wri_file=False
	res_gen = dedup._get_trimmed_seq_generator(iter(seqgen), itspos, wri_file=False)

	# Convert generator to list. Under the old buggy code, this would have been empty [].
	# Under the new robust code, it remains unexhausted and fully functional!
	res_list = list(res_gen)
	assert len(res_list) == 1
	assert len(res_list[0].seq) == 10

def test_trim_ccs_stitching():
	"""Verifies that enabling trim_ccs properly stitches fake primers and sets maximum quality scores."""
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord

	class MockItsPosition:
		def get_position(self, seq_id):
			return (5, 15, 100)

	itspos = MockItsPosition()

	uc = os.path.join(TEST_DIR, "test_data", "ex_tmpdir", "uc.txt")
	seq = os.path.join(TEST_DIR, "test_data", "ex_tmpdir", "seq.fq.gz")
	rep = os.path.join(TEST_DIR, "test_data", "ex_tmpdir", "rep.fa")
	dedup = itsxpress.main.Dedup(uc_file=uc, rep_file=rep, seq_file=seq)
	dedup.matchdict = {"seq1": "seq1"}

	original_seq = "N" * 50
	seqgen = [SeqRecord(Seq(original_seq), id="seq1", description="")]
	# Set simple mock quality scores
	seqgen[0].letter_annotations["phred_quality"] = [40] * len(original_seq)

	# Get generator with trim_ccs=True
	res_gen = dedup._get_trimmed_seq_generator(iter(seqgen), itspos, wri_file=True, trim_ccs=True)
	res_list = list(res_gen)

	assert len(res_list) == 1
	stitched = res_list[0]

	fwd_primer = "GACAGGTACAAGAAGGA"
	rev_primer_rc = "TTAACCCAGTCTCCAGT"
	expected_seq = fwd_primer + "N" * 10 + rev_primer_rc

	assert str(stitched.seq) == expected_seq

	# Quality scores should align perfectly: fwd_primer has 93s, trimmed sequence has 40s, rev_primer_rc has 93s
	expected_quals = [93] * 17 + [40] * 10 + [93] * 17
	assert stitched.letter_annotations["phred_quality"] == expected_quals
