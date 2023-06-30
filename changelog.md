2.0.0 (2023-06-28)
------------------
Release Highlights
- Removed BBmap dependency
    - BBmap scripts are no longer used in the pipeline, including:
        - reformat.sh (interleaved files no longer supported)
        - bbmerge.sh (merging of paired-end reads now done with Vsearch --fastq_mergepairs)
             - merging of paired-end reads is different between Vsearch and BBmerge, so results may differ
- Updated dereplication step for newer versions of Vsearch
    - Dereplication step now done using Vsearch --fastx_uniques (derep_fulllength command no longer supports fastq files)
- Qiime2 plugin version of ITSxpress is now part of the standalone package of ITSxpress
    - No longer need to install q2-itsxpress separately and will be installed if Qiime2 is already installed, otherwise only the standalone version will be installed
    - Qiime2 plugin version of ITSxpress is no longer maintained after q2-itsxpress v1.8.1
- Pacbio sequences are now supported if fastq file scores are in Illumina format
Bug Fixes
- Fixed bug where the q2-itsxpress plugin was not handling single-end reads correctly, and was looking for a reverse read file

1.8.1 (2023-06-02)
------------------
Release Highlights

- This is the final version that uses BBmap scripts. The next version will remove the BBmap dependency and use Vsearch for all steps.
    - This version can still be found on the EOL-1.8.1 branch of the repository
     
- Replaced derep_fulllength with fastx_uniques command compatible with Vsearch>=2.21.1
- Fixed version of dependencies
- Updated pip install config file to pyproject.toml
- Updated readme usage section to reference compatible Qiime2 version
- Added read count output to log file


1.8.0 (2019-12-9)
-----------------
- Added support for primer sets in the reverse orientation
- Fixed a bug that could cause crashes when an intermediate file was empty

1.7.2 (2018-11-8)
-----------------
- This release fixes issue [#8](https://github.com/USDA-ARS-GBRU/itsxpress/issues/8)
    - This issue caused ITSxpress to incorrectly trim about 0.2% of read pairs. Sometimes this would result in it writing blank fastq records which would cause Qiime to detect an error and stop processing.

1.7.0 (2018-09-12)
------------------
New Features:

- Support for the output of unmerged paired end files. This allows users to use Dada2 for sequence variant calling.
- The API is now documented at ReadTheDocs

1.6.4 (2018-7-26)
-----------------
- Fix for issue validating fastq.gz files that was not solved by v1.6.3

1.6.3 (2018-7-25)
-----------------
- Fixed issue validating fastq.gz and added tests.

1.6.2 (2018-7-25)
-----------------
- This release fixes an error that occasionally occurred when validating FASTQ files. ITSxpress used BBtools reformat.sh which occasionally threw an exception when validating FASTQ files due to a race condition. FASTQ file validation is now done with Biopython instead.

1.6.1 (2018-7-19)
-----------------
- Changed the default clustering identity to 99.5%.
- Experiments with fungal soil samples showed that ITSxpress and ITSx trimmed 99.822% of reads in the ITS1 region within 2 bases of each other and 99.099% of reads in the ITS2 region within 2 bases of each other at 99.5% identity. For higher accuracy, dereplication can be run at at 100% identity.


1.6.0 (2018-7-13)
-----------------
- This release adds a new feature to cluster merged sequences at less than 100% identity. This speeds up typical dataset trimming by about 10x over previous versions depending on the sample, without major effects on trimming accuracy. This feature is controlled with the --cluster_id flag. Default behavior is now to cluster at 0.987 identity.

1.5.6 (2018-7-3)
-----------------
- Database is now included in the release

1.5.4 (2018-6-27)
-----------------
- Fixed bug in handling of temporary files files specified with the --tempfile flag
- updated readme
- fixed issues with exception handling


1.5.2 (2018-6-21)
-----------------
- Fixed an indexing error causing ITS trimming to be off by 1 base.
- Fixed error when raising file not found exception
- removed old readme