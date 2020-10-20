import os
import unittest
from sys import getsizeof

from nose.tools import eq_, raises
from q2_types.per_sample_sequences import (SingleLanePerSampleSingleEndFastqDirFmt,
                                           SingleLanePerSamplePairedEndFastqDirFmt,
                                           FastqManifestFormat)

import itsxpress._itsxpress as _itsxpress

import qiime2
from qiime2.util import redirected_stdio
import pandas as pd

# The test data dir
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
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
    exp1 = _itsxpress._taxa_prefix_to_taxa(taxa_prefix="A")
    eq_(exp1, "Alveolata")


class SetFastqsAndCheckTests(unittest.TestCase):
    def test_single(self):
        samples = TEST_DATA_SINGLEIN.manifest.view(pd.DataFrame)
        for sample in samples.itertuples():
            obs = _itsxpress._set_fastqs_and_check(fastq=sample.forward,
                                                   fastq2=None,
                                                   sample_id=sample.Index,
                                                   single_end=True,
                                                   reversed_primers=False,
                                                   threads=1)
            self.assertTrue("4774-1-MSITS3" in obs.fastq)

    def test_paired(self):
        samples = TEST_DATA.manifest.view(pd.DataFrame)
        for sample in samples.itertuples():
            obs = _itsxpress._set_fastqs_and_check(fastq=sample.forward,
                                                   fastq2=sample.reverse,
                                                   sample_id=sample.Index,
                                                   single_end=True,
                                                   reversed_primers=False,
                                                   threads=1)
            self.assertTrue("4774-1-MSITS3" in obs.fastq)

    def test_paired_unmerged(self):
        samples = TEST_DATA.manifest.view(pd.DataFrame)
        for sample in samples.itertuples():
            obs = _itsxpress._set_fastqs_and_check(fastq=sample.forward,
                                                   fastq2=sample.reverse,
                                                   sample_id=sample.Index,
                                                   single_end=False,
                                                   reversed_primers=False,
                                                   threads=1)
            self.assertTrue("4774-1-MSITS3" in obs.fastq)
            self.assertTrue("4774-1-MSITS3" in obs.fastq2)

    def test_trim_pair_no_bb(self):
        samples = TEST_DATA.manifest.view(pd.DataFrame)
        for sample in samples:
            raises(ValueError, lambda: _itsxpress._set_fastqs_and_check(fastq=sample.forward,
                                                                        fastq2=sample.reverse,
                                                                        sample_id=sample.Index,
                                                                        single_end=False,
                                                                        reversed_primers=False,
                                                                        threads=1))
            raises(ValueError, lambda: _itsxpress._set_fastqs_and_check(fastq=sample.forward,
                                                                        fastq2=sample.reverse,
                                                                        sample_id=sample.Index,
                                                                        single_end=True,
                                                                        reversed_primers=False,
                                                                        threads=1))


def test_trim_single_no_cluster():
    threads = 1
    taxa = "F"
    region = "ITS2"
    cluster_id = 1

    exp1 = _itsxpress.trim_single(per_sample_sequences=TEST_DATA_SINGLEIN,
                                threads=threads,
                                taxa=taxa,
                                region=region,
                                cluster_id=cluster_id)
    exp2 = getsizeof(exp1)
    exp3 = getsizeof(TEST_DATA_SINGLEOUT)
    eq_(exp2, exp3)

def test_trim_pair_no_hmmer():
    threads = 1
    taxa = "F"
    region = "ITS2"

    raises(ValueError, lambda: _itsxpress.trim_pair(per_sample_sequences=TEST_DATA,
                                                    threads=threads,
                                                    taxa=taxa,
                                                    region=region))

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
