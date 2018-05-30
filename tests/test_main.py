# -*- coding: utf-8 -*-
import os
import tempfile
import shutil
import gzip
import subprocess

from Bio import SeqIO
from nose.tools import assert_raises

import itsxpress

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
from itsxpress.definitions import ROOT_DIR, taxa_dict
hmmfile = os.path.join(ROOT_DIR,"ITSx_db","HMMs", taxa_dict["Fungi"])


def test_check_fastqs():
	fastq = os.path.join(TEST_DIR, "test_data", "4774-1-MSITS3_R1.fastq")
	fastq2 = os.path.join(TEST_DIR, "test_data", "broken.fastq")
	assert_raises(subprocess.CalledProcessError, itsxpress.main._check_fastqs, fastq, fastq2)

def test_its_position_init():
	itspos = itsxpress.main.ItsPosition(os.path.join(TEST_DIR, "test_data", "ex_temdir", "domtbl.txt"), "ITS2")
	exp1 = {'left': {'score': '53.7', 'pos': '128'}, 'right': {'score': '60.0', 'pos': '282'}}
	assert exp1 == itspos.ddict["M02696:28:000000000-ATWK5:1:1101:19331:3209"]
	exp2 = {'right': {'score': '35.1', 'pos': '327'}}
	assert exp2 == itspos.ddict["M02696:28:000000000-ATWK5:1:1101:23011:4341"]
	assert len(itspos.ddict) == 137

def test_dedup():
	uc = os.path.join(TEST_DIR, "test_data", "ex_temdir", "uc.txt")
	seq = os.path.join(TEST_DIR, "test_data", "ex_temdir", "seq.fq.gz")
	rep = os.path.join(TEST_DIR, "test_data", "ex_temdir", "rep.fa")
	dedup = itsxpress.main.Dedup( uc_file=uc, rep_file=rep, seq_file=seq)
	# Check length of records
	assert len(dedup.matchdict) == 227
	# Check that non-representative seqs are logged
	assert dedup.matchdict['M02696:28:000000000-ATWK5:1:1101:11740:1800'] == 'M02696:28:000000000-ATWK5:1:1101:10899:1561'
	# Check that representative seqs are logged
	assert dedup.matchdict["M02696:28:000000000-ATWK5:1:1101:23011:4341"] == 'M02696:28:000000000-ATWK5:1:1101:23011:4341'

def test_dedup_create_trimmed_seqs():
	tf = tempfile.mkdtemp()
	uc = os.path.join(TEST_DIR, "test_data", "ex_temdir", "uc.txt")
	seq = os.path.join(TEST_DIR, "test_data", "ex_temdir", "seq.fq.gz")
	rep = os.path.join(TEST_DIR, "test_data", "ex_temdir", "rep.fa")
	dedup = itsxpress.main.Dedup( uc_file=uc, rep_file=rep, seq_file=seq)
	itspos = itsxpress.main.ItsPosition(os.path.join(TEST_DIR, "test_data", "ex_temdir", "domtbl.txt"), "ITS2")
	# Check non gzipped
	dedup.create_trimmed_seqs(os.path.join(tf,"testout.fastq"), gzipped=False, itspos=itspos)
	with open(os.path.join(tf,"testout.fastq"), 'r') as f:
		recs = SeqIO.parse(f, "fastq")
		n = 0
		length = 0
		for rec in recs:
			n += 1
			length += len(rec)
	assert n == 226
	assert length == 42863
	shutil.rmtree(tf)

def test_dedup_create_trimmed_seqs_gzipped():
	tf = tempfile.mkdtemp()
	uc = os.path.join(TEST_DIR, "test_data", "ex_temdir", "uc.txt")
	seq = os.path.join(TEST_DIR, "test_data", "ex_temdir", "seq.fq.gz")
	rep = os.path.join(TEST_DIR, "test_data", "ex_temdir", "rep.fa")
	dedup = itsxpress.main.Dedup( uc_file=uc, rep_file=rep, seq_file=seq)
	itspos = itsxpress.main.ItsPosition(os.path.join(TEST_DIR, "test_data", "ex_temdir", "domtbl.txt"), "ITS2")
	# Check non gzipped
	dedup.create_trimmed_seqs(os.path.join(tf,"testout.fastq.gz"), gzipped=True, itspos=itspos)
	with gzip.open(os.path.join(tf,"testout.fastq.gz"), 'rt') as f:
		recs = SeqIO.parse(f, "fastq")
		n = 0
		length = 0
		for rec in recs:
			n += 1
			length += len(rec)
	assert n == 226
	assert length == 42863
	shutil.rmtree(tf)

def test_seq_sample_paired_interleaved():
	fastq = os.path.join(TEST_DIR, "test_data", "4774-1-MSITS3_interleaved.fastq")
	sobj = itsxpress.main.SeqSamplePairedInterleaved(fastq=fastq, tempdir=".")
	sobj._merge_reads(threads=1)
	sobj._deduplicate(threads=1)
	sobj._search(hmmfile=hmmfile, threads=1)
	shutil.rmtree(sobj.tempdir)

def test_seq_sample_not_paired():
	fastq = os.path.join(TEST_DIR, "test_data", "4774-1-MSITS3_merged.fastq")
	sobj = itsxpress.main.SeqSampleNotPaired(fastq=fastq, tempdir=".")
	sobj._deduplicate(threads=1)
	sobj._search(hmmfile=hmmfile, threads=1)
	shutil.rmtree(sobj.tempdir)
	
def test_seq_sample_paired_not_interleaved():
	fastq = os.path.join(TEST_DIR, "test_data", "4774-1-MSITS3_R1.fastq")
	fastq2 = os.path.join(TEST_DIR, "test_data", "4774-1-MSITS3_R2.fastq")
	sobj = itsxpress.main.SeqSamplePairedNotInterleaved(fastq=fastq, tempdir=".", fastq2=fastq2)
	sobj._merge_reads(threads=1)
	sobj._deduplicate(threads=1)
	sobj._search(hmmfile=hmmfile, threads=1)
	shutil.rmtree(sobj.tempdir)
	
	
def test_is_paired():
	paired_end, interleaved = itsxpress.main._is_paired("fastq1.fq", "fastq2.fq", single_end=False)
	assert paired_end == True and interleaved == False
	paired_end, interleaved = itsxpress.main._is_paired("fastq1.fq", None, single_end=False)
	assert paired_end == True and interleaved == True
	paired_end, interleaved = itsxpress.main._is_paired("fastq1.fq", None, single_end=True)
	assert paired_end == False and interleaved == False


