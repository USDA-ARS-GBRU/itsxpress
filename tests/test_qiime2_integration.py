# -*- coding: utf-8 -*-
import os
import tempfile
import unittest
from sys import getsizeof

import pytest
import pandas as pd

# Check if qiime2 is available
try:
    import qiime2
    from qiime2.util import redirected_stdio
    from q2_types.per_sample_sequences import (
        SingleLanePerSampleSingleEndFastqDirFmt,
        SingleLanePerSamplePairedEndFastqDirFmt,
        FastqManifestFormat,
    )
    import itsxpress.q2_itsxpress as q2_itsxpress

    HAS_QIIME2 = True
except ModuleNotFoundError:
    HAS_QIIME2 = False

# Skip all tests in this module if QIIME2 is not installed
pytestmark = pytest.mark.skipif(
    not HAS_QIIME2,
    reason="QIIME2 is not installed in the current environment. Skipping QIIME2 integration tests.",
)

if HAS_QIIME2:
    TEST_DIR = os.path.dirname(os.path.abspath(__file__))
    tempdir = tempfile.mkdtemp()

    # Test file paths
    TEST_FILE = os.path.join(
        TEST_DIR, "test_data", "paired", "445cf54a-bf06-4852-8010-13a60fa1598c", "data"
    )
    TEST_DATA = SingleLanePerSamplePairedEndFastqDirFmt(TEST_FILE, "r")

    TEST_FILE_PBMD = os.path.join(
        TEST_DIR,
        "test_data",
        "pairedBrokenMissingData",
        "50d5f31a-a761-4c04-990c-e7668fe6bf00",
        "data",
    )
    TEST_DATA_PBMD = SingleLanePerSamplePairedEndFastqDirFmt(TEST_FILE_PBMD, "r")

    TEST_FILE_PAF = os.path.join(
        TEST_DIR,
        "test_data",
        "pairedAllForward",
        "445cf54a-bf06-4852-8010-13a60fa1598c",
        "data",
    )
    TEST_DATA_PAF = SingleLanePerSamplePairedEndFastqDirFmt(TEST_FILE_PAF, "r")

    TEST_FILE_OUT = os.path.join(
        TEST_DIR, "test_data", "out", "d9955749-00d5-44ae-a628-4b2da43000e1", "data"
    )
    TEST_DATA_OUT = SingleLanePerSamplePairedEndFastqDirFmt(TEST_FILE_OUT, "r")

    TEST_FILE_SINGLEOUT = os.path.join(
        TEST_DIR,
        "test_data",
        "singleOut",
        "75aea4f5-f10e-421e-91d2-feda9fe7b2e1",
        "data",
    )
    TEST_DATA_SINGLEOUT = SingleLanePerSamplePairedEndFastqDirFmt(
        TEST_FILE_SINGLEOUT, "r"
    )

    TEST_FILE_SINGLEIN = os.path.join(
        TEST_DIR,
        "test_data",
        "singleIn",
        "cfd0e65b-05fb-4329-9618-15ecd0aec9b3",
        "data",
    )
    TEST_DATA_SINGLEIN = SingleLanePerSampleSingleEndFastqDirFmt(
        TEST_FILE_SINGLEIN, "r"
    )

    def test_taxa_prefix_to_taxa():
        exp1 = q2_itsxpress._taxa_prefix_to_taxa(taxa_prefix="A")
        assert exp1 == "Alveolata"

    class SetFastqsAndCheckTests(unittest.TestCase):
        def test_single(self):
            samples = TEST_DATA_SINGLEIN.manifest.view(pd.DataFrame)
            for sample in samples.itertuples():
                obs = q2_itsxpress._set_fastqs_and_check(
                    fastq=sample.forward,
                    fastq2=None,
                    tempdir=tempdir,
                    sample_id=sample.Index,
                    single_end=True,
                    reversed_primers=False,
                    allow_staggered_reads=False,
                    threads=1,
                )
                self.assertTrue("4774-1-MSITS3" in obs.fastq)

        def test_paired(self):
            samples = TEST_DATA.manifest.view(pd.DataFrame)
            for sample in samples.itertuples():
                obs = q2_itsxpress._set_fastqs_and_check(
                    fastq=sample.forward,
                    fastq2=sample.reverse,
                    tempdir=tempdir,
                    sample_id=sample.Index,
                    single_end=True,
                    reversed_primers=False,
                    allow_staggered_reads=False,
                    threads=1,
                )
                self.assertTrue("4774-1-MSITS3" in obs.fastq)

    def test_trim_single_no_cluster():
        threads = 1
        taxa = "F"
        region = "ITS2"
        cluster_id = 1

        exp1 = q2_itsxpress.trim_single(
            per_sample_sequences=TEST_DATA_SINGLEIN,
            threads=threads,
            taxa=taxa,
            region=region,
            cluster_id=cluster_id,
        )
        exp2 = getsizeof(exp1)
        exp3 = getsizeof(TEST_DATA_SINGLEOUT)
        assert exp2 == exp3

    class TrimTests(unittest.TestCase):
        def setUp(self):
            self.plugin = qiime2.sdk.PluginManager().plugins["itsxpress"]
            self.trim_single_fn = self.plugin.methods["trim_single"]
            self.trim_paired_fn = self.plugin.methods["trim_pair"]
            self.trim_paired_unmerged_fn = self.plugin.methods[
                "trim_pair_output_unmerged"
            ]

            self.se_seqs = qiime2.Artifact.import_data(
                "SampleData[SequencesWithQuality]",
                TEST_FILE_SINGLEIN,
                "SingleLanePerSampleSingleEndFastqDirFmt",
            )
            self.pe_seqs = qiime2.Artifact.import_data(
                "SampleData[PairedEndSequencesWithQuality]",
                TEST_FILE,
                "SingleLanePerSamplePairedEndFastqDirFmt",
            )

        def test_trim_single_success(self):
            with redirected_stdio(stderr=os.devnull):
                (obs_artifact,) = self.trim_single_fn(self.se_seqs, "ITS2")
            self.assertEqual(str(obs_artifact.type), "SampleData[SequencesWithQuality]")

            obs_dir = obs_artifact.view(SingleLanePerSampleSingleEndFastqDirFmt)
            self.assertEqual(getsizeof(obs_dir), getsizeof(TEST_DATA_SINGLEOUT))

            obs = obs_artifact.view(SingleLanePerSampleSingleEndFastqDirFmt)
            obs_manifest = list(obs.manifest.view(FastqManifestFormat).open())
            exp_manifest = [
                "sample-id,filename,direction\n",
                "4774-1-MSITS3,4774-1-MSITS3_0_L001_R1_001.fastq.gz,forward\n",
            ]
            self.assertEqual(obs_manifest, exp_manifest)

        def test_trim_pair_success(self):
            with redirected_stdio(stderr=os.devnull):
                (obs_artifact,) = self.trim_paired_fn(self.pe_seqs, "ITS2")
            self.assertEqual(
                str(obs_artifact.type), "SampleData[JoinedSequencesWithQuality]"
            )

            obs_dir = obs_artifact.view(SingleLanePerSampleSingleEndFastqDirFmt)
            self.assertEqual(getsizeof(obs_dir), getsizeof(TEST_DATA_OUT))

            obs = obs_artifact.view(SingleLanePerSampleSingleEndFastqDirFmt)
            obs_manifest = list(obs.manifest.view(FastqManifestFormat).open())
            exp_manifest = [
                "sample-id,filename,direction\n",
                "4774-1-MSITS3,4774-1-MSITS3_0_L001_R1_001.fastq.gz,forward\n",
            ]
            self.assertEqual(obs_manifest, exp_manifest)

        def test_trim_pair_output_unmerged_success(self):
            with redirected_stdio(stderr=os.devnull):
                (obs_artifact,) = self.trim_paired_unmerged_fn(self.pe_seqs, "ITS2")
            self.assertEqual(
                str(obs_artifact.type), "SampleData[PairedEndSequencesWithQuality]"
            )
            obs = obs_artifact.view(SingleLanePerSamplePairedEndFastqDirFmt)
            obs_manifest = list(obs.manifest.view(FastqManifestFormat).open())
            exp_manifest = [
                "sample-id,filename,direction\n",
                "4774-1-MSITS3,4774-1-MSITS3_0_L001_R1_001.fastq.gz,forward\n",
                "4774-1-MSITS3,4774-1-MSITS3_1_L001_R2_001.fastq.gz,reverse\n",
            ]
            self.assertEqual(obs_manifest, exp_manifest)
