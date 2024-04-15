ITSxpress: Software to rapidly trim  the Internally transcribed spacer (ITS) region of FASTQ files
==================================================================================================

.. image:: https://readthedocs.org/projects/itsxpress-package/badge/?version=latest
    :target: https://itsxpress-package.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://github.com/USDA-ARS-GBRU/itsxpress/actions/workflows/python-package-conda.yml/badge.svg
   :target: https://github.com/USDA-ARS-GBRU/itsxpress/actions/workflows/python-package-conda.yml
   :alt: Build Status

.. image:: https://anaconda.org/bioconda/itsxpress/badges/downloads.svg
   :target: https://anaconda.org/bioconda/itsxpress
   :alt: Anaconda-Server Badge
   
.. image:: https://img.shields.io/github/v/release/USDA-ARS-GBRU/itsxpress?style=social
   :target: https://github.com/USDA-ARS-GBRU/itsxpress/releases/latest
   :alt: GitHub release (latest by date)

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1304349.svg
  :target: https://doi.org/10.5281/zenodo.1304349

Author
-------
* Adam R. Rivers, US Department of Agriculture, Agricultural Research Service
* Sveinn V. Einarsson, US Department of Agriculture, Agricultural Research Service

Citation
--------
Rivers AR, Weber KC, Gardner TG et al. ITSxpress: Software to rapidly trim
internally transcribed spacer sequences with quality scores for marker gene
analysis [version 1; referees: awaiting peer review]. F1000Research 2018, 7:1418
(doi: `10.12688/f1000research.15704.1`_)

.. _`10.12688/f1000research.15704.1`: https://doi.org/10.12688/f1000research.15704.1

Introduction
-------------

The internally transcribed spacer region is a region between highly conserved the small
subunit (SSU) of rRNA and the large subunit (LSU) of the rRNA. In Eukaryotes it contains
the 5.8s genes and two variable length spacer regions. In amplicon sequencing studies it is
common practice to trim off the conserved (SSU, 5,8S or LSU) regions. `Bengtsson-Palme
et al. (2013)`_ published software the software package ITSx_ to do this.

ITSxpress is designed to support the calling of exact sequence variants rather than OTUs_.
This newer method of sequence error-correction requires quality score data from each
sequence, so each input sequence must be trimmed. ITSXpress makes this possible by
taking FASTQ data, de-replicating the sequences then identifying the start and stop
sites using HMMSearch.  Results are parsed and the trimmed files are returned. The ITS 1,
ITS2 or the entire ITS region including the 5.8s rRNA gene can be selected. ITSxpress
uses the hmm model from ITSx so results are comparable.

ITSxpress is also a QIIME2 plugin. Starting from version 2.0.0 of ITSxpress, the QIIME2 plugin is included with
the command line version of ITSxpress. The installation method will be slightly different depending on whether 
QIIME2 is being used.

.. _`Bengtsson-Palme et al. (2013)`: https://doi.org/10.1111/2041-210X.12073
.. _ITSx: http://microbiology.se/software/itsx/
.. _OTUs: https://doi.org/10.1038/ismej.2017.119
.. _`QIIME2 Plugin`: https://github.com/USDA-ARS-GBRU/q2_itsxpress


Installing ITSxpress for use as a QIIME2 Plugin
----------------------------------------------------

To install ITSxpress as a plugin for QIIME 2 first install QIIME 2 as a separate Conda/Mamba environemnt using thier instructions 
https://docs.qiime2.org/2024.2/install/ then add ITSxress to the QIIME 2 Conda environment. The examples below are for QIIME2 2 
version 2024.2 an so please update the commands if you want a newer release.
 

For Linux:

.. code-block:: bash

    wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.2-py38-linux-conda.yml
    mamba create -n qiime2-amplicon-2024.2 --file qiime2-amplicon-2024.2-py38-linux-conda.yml
    mamba activate qiime2-amplicon-2024.2
    mamba install ITSxpress
    qiime dev refresh-cache

For maxOS (Intel) and OS X:

.. code-block:: bash

    wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.2-py38-osx-conda.yml
    mamba create -n qiime2-amplicon-2024.2 --file qiime2-amplicon-2024.2-py38-osx-conda.yml
    mamba activate qiime2-amplicon-2024.2
    mamba install ITSxpress
    qiime dev refresh-cache

For macOS (Arm / Apple Silicon):

.. code-block:: bash

    wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.2-py38-osx-conda.yml
    CONDA_SUBDIR=osx-64 mamba create -n qiime2-amplicon-2024.2 --file qiime2-amplicon-2024.2-py38-osx-conda.yml
    mamba activate qiime2-amplicon-2024.2
    mamba config --env --set subdir osx-64
    mamba install ITSxpress
    qiime dev refresh-cache


Installing ITSxpress for standalone use
-------------------------------------------

For Linux, maxOS (Intel), and OS X:

.. code-block:: bash

    mamba create -n itsxpressenv -c bioconda -c conda-forge itsxpress
    mamba activate itsxpressenv 


For macOS (Arm/Apple Silicon):

.. code-block:: bash

    CONDA_SUBDIR=osx-64 mamba create -n itsxpressenv -c bioconda -c conda-forge itsxpress
    mamba activate itsxpressenv
    conda  config --env --set subdir osx-64


Running ITSxpress as a Docker container
-------------------------------------------

.. code-block:: bash
    
    docker pull ghcr.io/usda-ars-gbru/itsxpress
    docker run [Options...] itsxpress


Dependencies
-------------
The software requires Vsearch, Hmmer and Biopython. Bioconda
takes care of this for you so it is the preferred installation method.


Usage
---------


+-------------------------+---------------------------------------------------------------+
| Option                  | Description                                                   |
+=========================+===============================================================+
| -h, --help              | Show this help message and exit.                              |
+-------------------------+---------------------------------------------------------------+
| --fastq                 | A ``.fastq``, ``.fq``, ``.fastq.gz`` or ``.fq.gz`` file.      |
|                         | Required.                                                     |
+-------------------------+---------------------------------------------------------------+
| --single_end            | A flag to specify that the fastq file is single-ended (not    |
|                         | paired). Default is false.                                    |
+-------------------------+---------------------------------------------------------------+
| --fastq2                | A ``.fastq``, ``.fq``, ``.fastq.gz`` or ``.fq.gz`` file       |
|                         | representing read 2 if present, optional.                     |
+-------------------------+---------------------------------------------------------------+
| --outfile               | The trimmed FASTQ file, if it ends in ``gz`` it will be       |
|                         | gzipped.                                                      |
+-------------------------+---------------------------------------------------------------+
| --outfile2              | The trimmed FASTQ read 2 file, if it ends in ``gz`` it will   |
|                         | be gzipped. If used, reads will be retuned as unmerged pairs  |
|                         | rather than than merged.                                      |
+-------------------------+---------------------------------------------------------------+
| --tempdir               | Specify the temp file directory. Default is None.             |
+-------------------------+---------------------------------------------------------------+
| --keeptemp              | Should intermediate files be kept? Default is false.          |
+-------------------------+---------------------------------------------------------------+
| --region                | Options : {ITS2, ITS1, ALL}                                   |
+-------------------------+---------------------------------------------------------------+
| --taxa                  | Select the taxonomic group sequenced: {Alveolata, Bryophyta,  |
|                         | Bacillariophyta, Amoebozoa, Euglenozoa, Fungi, Chlorophyta,   |
|                         | Rhodophyta, Phaeophyceae, Marchantiophyta, Metazoa, Oomycota, |
|                         | Haptophyceae, Raphidophyceae, Rhizaria, Synurophyceae,        |
|                         | Tracheophyta, Eustigmatophyceae, Parabasalia, All}.           |
|                         | Default Fungi.                                                |
+-------------------------+---------------------------------------------------------------+
| --cluster_id            | The percent identity for clustering reads range [0.99-1.0],   |
|                         | set to 1 for exact de-replication. Default 1.0.               |
+-------------------------+---------------------------------------------------------------+
| --log                   | Log file. Default is ITSxpress.log.                           |
+-------------------------+---------------------------------------------------------------+
| --threads               | Number of processor threads to use. Default is 1.             |
+-------------------------+---------------------------------------------------------------+
| --reversed_primers      | Primers are in reverse orientation as in Taylor et al. 2016,  |
|                         | DOI:10.1128/AEM.02576-16. If selected ITSxpress returns       |
|                         | trimmed reads flipped to the forward orientation              |
+-------------------------+---------------------------------------------------------------+
| --allow_staggered_reads | Allow merging staggered reads with --fastq_allowmergestagger  |
|                         | for Vsearch --fastq_mergepairs. See Vsearch documentation.    |
|                         | (Optional) Default is true.                                   |
+-------------------------+---------------------------------------------------------------+



Examples
---------

Use case 1: Trimming the ITS2 region from a fungal amplicon sequencing dataset with
forward and reverse gzipped FASTQ files using two cpu threads. Return a single merged file for use in Deblur.

.. code-block:: bash

    itsxpress --fastq r1.fastq.gz --fastq2 r2.fastq.gz --region ITS2 \
    --taxa Fungi --log logfile.txt --outfile trimmed_reads.fastq.gz --threads 2

ITSxpress can take uncompressed, gzipped or zstd compressed FASTQ files and it can write uncompressed, gzipped or
zstd compressed FASTQ files. It expects FASTQ files to end in: .fq, .fastq, .fq.gz, fastq.gz, .fq.zst or fastq.zst.

Use case 2: Trimming the ITS2 region from a fungal amplicon sequencing dataset with
forward and reverse gzipped FASTQ files using two cpu threads. Return a forward
and reverse read files  for use in Dada2.

.. code-block:: bash

    itsxpress --fastq r1.fastq.gz --fastq2 r2.fastq.gz --region ITS2 \
    --taxa Fungi --log logfile.txt --outfile trimmed_reads.fastq.gz --threads 2

ITSxpress can take uncompressed, gzipped or zstd compressed FASTQ files and it can write uncompressed, gzipped or
zstd compressed FASTQ files. It expects FASTQ files to end in: .fq, .fastq, .fq.gz, fastq.gz, .fq.zst or fastq.zst.


Use case 3: Trimming the ITS2 region from a fungal amplicon sequencing dataset with
an single-ended gzipped FASTQ files using two cpu threads.

.. code-block:: bash

    itsxpress --fastq single-end.fastq.gz --single_end --region ITS2 --taxa Fungi \
    --log logfile.txt --outfile trimmed_reads.fastq.gz --threads 2

Single ended data is less common and may come from a dataset where the reads have already
been merged.

License information
--------------------
This software is a work of the United States Department of Agriculture,
Agricultural Research Service and is released under a Creative Commons CC0
public domain attribution.

