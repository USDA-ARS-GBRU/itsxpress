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

TEST_DIR = os.path.dirname(os.path.abspath(__name__))
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

try:
	import qiime2
	from qiime2.util import redirected_stdio
	import os
	import unittest
	from sys import getsizeof

	#from nose.tools import eq_, raises
	from q2_types.per_sample_sequences import (SingleLanePerSampleSingleEndFastqDirFmt,
                                           SingleLanePerSamplePairedEndFastqDirFmt,
                                           FastqManifestFormat)

	import pandas as pd
	import itsxpress.q2_itsxpress as q2_itsxpress
	import itsxpress.plugin_setup

	tempdir = tempfile.mkdtemp()
	# The test data dir
	TEST_DIR = os.path.join(os.path.dirname(os.path.abspath(__name__)),"tests")
	# Test info 1
	TEST_FILE = os.path.join(TEST_DIR,
							"test_data",
							"paired",
							"445cf54a-bf06-4852-8010-13a60fa1598c",
							"data")

	TEST_DATA = SingleLanePerSamplePairedEndFastqDirFmt(TEST_FILE, "r")
	# Test info 2
	TEST_FILE_PBMD = os.path.join(TEST_DIR,
								"test_data",
								"pairedBrokenMissingData",
								"50d5f31a-a761-4c04-990c-e7668fe6bf00",
								"data")

	TEST_DATA_PBMD = SingleLanePerSamplePairedEndFastqDirFmt(TEST_FILE_PBMD, "r")
	# Test info 3
	TEST_FILE_PAF = os.path.join(TEST_DIR,
								"test_data",
								"pairedAllForward",
								"445cf54a-bf06-4852-8010-13a60fa1598c",
								"data")
	TEST_DATA_PAF = SingleLanePerSamplePairedEndFastqDirFmt(TEST_FILE_PAF, "r")
	# Test info 4
	TEST_FILE_OUT = os.path.join(TEST_DIR,
								"test_data",
								"out",
								"d9955749-00d5-44ae-a628-4b2da43000e1",
								"data")
	TEST_DATA_OUT = SingleLanePerSamplePairedEndFastqDirFmt(TEST_FILE_OUT, "r")
	# Test info 5
	TEST_FILE_SINGLEOUT = os.path.join(TEST_DIR,
									"test_data",
									"singleOut",
									"75aea4f5-f10e-421e-91d2-feda9fe7b2e1",
									"data")
	TEST_DATA_SINGLEOUT = SingleLanePerSamplePairedEndFastqDirFmt(TEST_FILE_SINGLEOUT, "r")
	# Test info 6
	TEST_FILE_SINGLEIN = os.path.join(TEST_DIR,
									"test_data",
									"singleIn",
									"cfd0e65b-05fb-4329-9618-15ecd0aec9b3",
									"data")
	TEST_DATA_SINGLEIN = SingleLanePerSampleSingleEndFastqDirFmt(TEST_FILE_SINGLEIN, "r")
	# Test artifact1
	ARTIFACT_TYPE_P = "SampleData[PairedEndSequencesWithQuality]"
	# Test artifact2
	ARTIFACT_TYPE_S = "SampleData[SequencesWithQuality]"


	def test_taxa_prefix_to_taxa():
		exp1 = q2_itsxpress._taxa_prefix_to_taxa(taxa_prefix="A")
		assert exp1, "Alveolata"


	class SetFastqsAndCheckTests(unittest.TestCase):
		def test_single(self):
			samples = TEST_DATA_SINGLEIN.manifest.view(pd.DataFrame)
			for sample in samples.itertuples():
				obs = q2_itsxpress._set_fastqs_and_check(fastq=sample.forward,
													fastq2=None,
													tempdir=tempdir,
													sample_id=sample.Index,
													single_end=True,
													reversed_primers=False,
													allow_staggered_reads=False,
													threads=1)
				self.assertTrue("4774-1-MSITS3" in obs.fastq)
		def test_paired(self):
			samples = TEST_DATA.manifest.view(pd.DataFrame)
			for sample in samples.itertuples():
				obs = q2_itsxpress._set_fastqs_and_check(fastq=sample.forward,
													fastq2=sample.reverse,
													tempdir=tempdir,
													sample_id=sample.Index,
													single_end=True,
													reversed_primers=False,
													allow_staggered_reads=False,
													threads=1)
				self.assertTrue("4774-1-MSITS3" in obs.fastq)

			def test_trim_single_no_cluster():
				threads = 1
				taxa = "F"
				sample = None  # replace with actual sample object
				# Fix for missing argument
				obs = q2_itsxpress._set_fastqs_and_check(fastq=sample.forward,
														 fastq2=None,
														 tempdir=tempdir,
														 sample_id=sample.Index,
														 single_end=True,
														 reversed_primers=False,
														 allow_staggered_reads=False,
														 threads=threads,
														 cluster=False)



	def test_trim_single_no_cluster():
		threads = 1
		taxa = "F"
		region = "ITS2"
		cluster_id = 1

		exp1 = q2_itsxpress.trim_single(per_sample_sequences=TEST_DATA_SINGLEIN,
									threads=threads,
									taxa=taxa,
									region=region,
									cluster_id=cluster_id)
		exp2 = getsizeof(exp1)
		exp3 = getsizeof(TEST_DATA_SINGLEOUT)
		assert exp2, exp3
# #Testing if there is no hmmer package?
# 	def test_trim_pair_no_hmmer():
# 		threads = 1
# 		taxa = "F"
# 		region = "ITS2"

# 		pytest.raises(ValueError, lambda: q2_itsxpress.trim_pair(per_sample_sequences=TEST_DATA,
# 														threads=threads,
# 														taxa=taxa,
# 														region=region))

	class TrimTests(unittest.TestCase):
		def setUp(self):
			self.plugin = qiime2.sdk.PluginManager().plugins['itsxpress']
			self.trim_single_fn = self.plugin.methods['trim_single']
			self.trim_paired_fn = self.plugin.methods['trim_pair']
			self.trim_paired_unmerged_fn = \
				self.plugin.methods['trim_pair_output_unmerged']

			self.se_seqs = qiime2.Artifact.import_data(
				'SampleData[SequencesWithQuality]',
				TEST_FILE_SINGLEIN,
				'SingleLanePerSampleSingleEndFastqDirFmt')
			self.pe_seqs = qiime2.Artifact.import_data(
				'SampleData[PairedEndSequencesWithQuality]',
				TEST_FILE,
				'SingleLanePerSamplePairedEndFastqDirFmt')

		def test_trim_single_success(self):
			with redirected_stdio(stderr=os.devnull):
				obs_artifact, = self.trim_single_fn(self.se_seqs, 'ITS2')
			self.assertEqual(str(obs_artifact.type),
							'SampleData[SequencesWithQuality]')

			obs_dir = obs_artifact.view(SingleLanePerSampleSingleEndFastqDirFmt)
			self.assertEqual(getsizeof(obs_dir), getsizeof(TEST_DATA_SINGLEOUT))

			obs = obs_artifact.view(SingleLanePerSampleSingleEndFastqDirFmt)
			obs_manifest = list(obs.manifest.view(FastqManifestFormat).open())
			exp_manifest = [
				'sample-id,filename,direction\n',
				'4774-1-MSITS3,4774-1-MSITS3_0_L001_R1_001.fastq.gz,forward\n']
			self.assertEqual(obs_manifest, exp_manifest)

		def test_trim_pair_success(self):
			with redirected_stdio(stderr=os.devnull):
				obs_artifact, = self.trim_paired_fn(self.pe_seqs, 'ITS2')
			self.assertEqual(str(obs_artifact.type),
							'SampleData[JoinedSequencesWithQuality]')

			obs_dir = obs_artifact.view(SingleLanePerSampleSingleEndFastqDirFmt)
			self.assertEqual(getsizeof(obs_dir), getsizeof(TEST_DATA_OUT))

			obs = obs_artifact.view(SingleLanePerSampleSingleEndFastqDirFmt)
			obs_manifest = list(obs.manifest.view(FastqManifestFormat).open())
			exp_manifest = [
				'sample-id,filename,direction\n',
				'4774-1-MSITS3,4774-1-MSITS3_0_L001_R1_001.fastq.gz,forward\n']
			self.assertEqual(obs_manifest, exp_manifest)

		def test_trim_pair_output_unmerged_success(self):
			with redirected_stdio(stderr=os.devnull):
				obs_artifact, = self.trim_paired_unmerged_fn(self.pe_seqs, 'ITS2')
			self.assertEqual(str(obs_artifact.type),
							'SampleData[PairedEndSequencesWithQuality]')
			obs = obs_artifact.view(SingleLanePerSamplePairedEndFastqDirFmt)
			obs_manifest = list(obs.manifest.view(FastqManifestFormat).open())
			exp_manifest = [
				'sample-id,filename,direction\n',
				'4774-1-MSITS3,4774-1-MSITS3_0_L001_R1_001.fastq.gz,forward\n',
				'4774-1-MSITS3,4774-1-MSITS3_1_L001_R2_001.fastq.gz,reverse\n']
			self.assertEqual(obs_manifest, exp_manifest)


except ModuleNotFoundError as e:
	#logging
	print("{}.Could not initialize the Qiime plugin portion of ITSxpress. Command line ITSxpress will still work normally. If you wish to use the Qiime2 ITSxpress plugin, you need to install Qiime2 first into your environment.\n".format(e))
	pass