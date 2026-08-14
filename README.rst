ITSxpress: Software to rapidly trim  the internal transcribed spacer (ITS) region of FASTQ files
==================================================================================================

.. image:: https://readthedocs.org/projects/itsxpress-package/badge/?version=latest
    :target: https://itsxpress-package.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://github.com/USDA-ARS-GBRU/itsxpress/actions/workflows/ci-cd.yml/badge.svg
   :target: https://github.com/USDA-ARS-GBRU/itsxpress/actions/workflows/ci-cd.yml
   :alt: CI/CD Pipeline (Test & Publish Docker)


.. image:: https://anaconda.org/bioconda/itsxpress/badges/downloads.svg
   :target: https://anaconda.org/bioconda/itsxpress
   :alt: Anaconda-Server Badge
   
.. image:: https://img.shields.io/github/v/release/USDA-ARS-GBRU/itsxpress?style=social
   :target: https://github.com/USDA-ARS-GBRU/itsxpress/releases/latest
   :alt: GitHub release (latest by date)

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1304348.svg
   :target: https://doi.org/10.5281/zenodo.1304348

Authors
-------

* Adam R. Rivers, US Department of Agriculture, Agricultural Research Service, Gainesville, FL, USA
* Sveinn V. Einarsson, Queens University, Charlotte, NC, USA


Citations
---------

Einarsson SV, Rivers AR. 2024. ITSxpress version 2: software to rapidly trim internal transcribed 
spacer sequences with quality scores for amplicon sequencing. Microbiol Spectr 0:e00601-24.
(doi: `10.1128/spectrum.00601-24`_)

.. _`10.1128/spectrum.00601-24`: https://doi.org/10.1128/spectrum.00601-24

Rivers AR, Weber KC, Gardner TG, Liu, S, Armstrong, SD. ITSxpress: Software to rapidly trim
internally transcribed spacer sequences with quality scores for marker gene
analysis [version 1]. F1000Research 2018, 7:1418
(doi: `10.12688/f1000research.15704.1`_)

.. _`10.12688/f1000research.15704.1`: https://doi.org/10.12688/f1000research.15704.1

Introduction
-------------

The internally transcribed spacer region is a region between highly conserved the small
subunit (SSU) of rRNA and the large subunit (LSU) of the rRNA. In Eukaryotes it contains
the 5.8s genes and two variable length spacer regions. In amplicon sequencing studies it is
common practice to trim off the conserved (SSU, 5,8S or LSU) regions. `Bengtsson-Palme
et al. (2013)`_ published software the software package ITSx_ to do this.

ITSxpress is designed to support the calling of exact amplicon sequence variants (ASVs) rather than OTUs_.
This newer method of sequence error-correction requires quality score data from each
sequence, so each input sequence must be trimmed. ITSXpress makes this possible by
taking FASTQ data, de-replicating the sequences with Vsearch then identifying the start and stop
sites using HMMER hmmsearch.  Results are parsed and the trimmed files are returned. The ITS1,
ITS2 or the entire ITS region including the 5.8s rRNA gene can be selected. ITSxpress
uses the hmm models from ITSx so results are comparable.

ITSxpress is also a QIIME2 plugin. Starting from version 2.0.0 of ITSxpress, the QIIME2 plugin is included with
the command line version of ITSxpress. The installation method will be slightly different depending on whether 
QIIME2 is being used.

.. _`Bengtsson-Palme et al. (2013)`: https://doi.org/10.1111/2041-210X.12073
.. _ITSx: https://microbiology.se/software/itsx/
.. _OTUs: https://doi.org/10.1038/ismej.2017.119



Installing ITSxpress for use as a QIIME2 plugin
----------------------------------------------------

To install ITSxpress as a plugin for QIIME 2, first install QIIME 2 Amplicon Distribution as a separate Conda/Mamba environemnt using their instructions 
https://library.qiime2.org/quickstart/qiime2 then add ITSxpress to the QIIME 2 Conda environment. The examples below are for QIIME2 
version 2026.7 an so please update the commands if you want a newer release.
 

For Linux:

.. code-block:: bash

    conda env create \
      --name rachis-qiime2-2026.7 \
      --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.7/qiime2/released/rachis-qiime2-linux-64-conda.yml
    conda activate rachis-qiime2-2026.7
    conda install -c bioconda -c conda-forge itsxpress
    qiime dev refresh-cache

For maxOS (Intel) and OS X:

.. code-block:: bash

    conda env create \
        --name rachis-qiime2-2026.7 \
        --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.7/qiime2/released/rachis-qiime2-osx-64-conda.yml
    conda activate rachis-qiime2-2026.7
    conda install -c bioconda -c conda-forge itsxpress
    qiime dev refresh-cache

For macOS (Arm / Apple Silicon):

.. code-block:: bash

    CONDA_SUBDIR=osx-64 conda env create \
     --name rachis-qiime2-2026.7 \
     --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.7/qiime2/released/rachis-qiime2-osx-64-conda.yml
    conda activate rachis-qiime2-2026.7
    conda config --env --set subdir osx-64
    CONDA_SUBDIR=osx-64 conda install -c bioconda -c conda-forge itsxpress
    qiime dev refresh-cache


Installing ITSxpress for standalone use
-------------------------------------------

For Linux, maxOS (Intel), and OS X:

.. code-block:: bash

    conda create -n itsxpressenv -c bioconda -c conda-forge itsxpress
    conda activate itsxpressenv 


For macOS (Arm/Apple Silicon):

.. code-block:: bash

    CONDA_SUBDIR=osx-64 conda create -n itsxpressenv -c bioconda -c conda-forge itsxpress
    conda activate itsxpressenv
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
------


+--------------------------+------------------------------------------------------------------------------------+
| Option                   | Description                                                                        |
+==========================+====================================================================================+
| -h, --help               | Show this help message and exit.                                                   |
+--------------------------+------------------------------------------------------------------------------------+
|| --fastq                 || A ``.fastq``, ``.fq``, ``.fastq.gz`` or ``.fq.gz`` file.                          |
||                         || Required.                                                                         |
+--------------------------+------------------------------------------------------------------------------------+
|| --single_end            || A flag to specify that the fastq file is single-ended (not                        |
||                         || paired). Default is false.                                                        |
+--------------------------+------------------------------------------------------------------------------------+
|| --fastq2                || A ``.fastq``, ``.fq``, ``.fastq.gz`` or ``.fq.gz`` file                           |
||                         || representing read 2 if present, optional.                                         |
+--------------------------+------------------------------------------------------------------------------------+
|| --outfile               || The trimmed FASTQ file, if it ends in ``gz`` it will be                           |
||                         || gzipped.                                                                          |
+--------------------------+------------------------------------------------------------------------------------+
|| --outfile2              || The trimmed FASTQ read 2 file, if it ends in ``gz`` it will                       |
||                         || be gzipped. If used, reads will be retuned as unmerged pairs                      |
||                         || rather than than merged.                                                          |
+--------------------------+------------------------------------------------------------------------------------+
| --tempdir                | Specify the temp file directory. Default is None.                                  |
+--------------------------+------------------------------------------------------------------------------------+
| --keeptemp               | Should intermediate files be kept? Default is false.                               |
+--------------------------+------------------------------------------------------------------------------------+
| --region                 | Options : {ITS2, ITS1, ALL}                                                        |
+--------------------------+------------------------------------------------------------------------------------+
|| --taxa                  || Select the taxonomic group sequenced: {Alveolata, Bryophyta,                      |
||                         || Bacillariophyta, Amoebozoa, Euglenozoa, Fungi, Chlorophyta,                       |
||                         || Rhodophyta, Phaeophyceae, Marchantiophyta, Metazoa, Oomycota,                     |
||                         || Haptophyceae, Raphidophyceae, Rhizaria, Synurophyceae,                            |
||                         || Tracheophyta, Eustigmatophyceae, Parabasalia, All}.                               |
||                         || Default Fungi.  For QIIME single letter codes are used, see help menu for list.   |
+--------------------------+------------------------------------------------------------------------------------+
|| --cluster_id            || The percent identity for clustering reads range [0.99-1.0],                       |
||                         || set to 1 for exact de-replication. Default 1.0.                                   |
+--------------------------+------------------------------------------------------------------------------------+
| --log                    | Log file. Default is ITSxpress.log.                                                |
+--------------------------+------------------------------------------------------------------------------------+
| --threads                | Number of processor threads to use. Default is 1. There is not much benefit over 4.|
+--------------------------+------------------------------------------------------------------------------------+
|| --reversed_primers      || Primers are in reverse orientation as in Taylor et al. 2016,                      |
||                         || DOI:10.1128/AEM.02576-16. If selected ITSxpress returns                           |
||                         || trimmed reads flipped to the forward orientation                                  |
+--------------------------+------------------------------------------------------------------------------------+
|| --allow_staggered_reads || Allow merging staggered reads with --fastq_allowmergestagger                      |
||                         || for Vsearch --fastq_mergepairs. See Vsearch documentation.                        |
||                         || (Optional) Default is true.                                                       |
+--------------------------+------------------------------------------------------------------------------------+
|| --trim-ccs              || PacBio CCS reads mode. Stitch fake primers (GACAGGTACAAGAAGGA                     |
||                         || and adapter ACTGGAGACTGGGTTAA) with uniform quality scores.                       |
+--------------------------+------------------------------------------------------------------------------------+

QIIME2 plugin Examples
-----------------------

Example 1: trimming fungal paired-end Illumina reads (outputting unmerged reads)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When processing Illumina  paired-end fungal data intended for downstream analysis with DADA2, you must maintain separate forward and reverse read structures. The ``trim-pair-output-unmerged`` 
command trims the ITS region from both reads but keeps them unmerged.

.. code-block:: bash

   qiime itsxpress trim-pair-output-unmerged \
     --i-per-sample-sequences paired-end-demux.qza \
     --p-region ITS1 \
     --p-taxa F \
     --p-threads 4 \
     --o-trimmed trimmed-paired-unmerged.qza \
     --verbose

**Key Parameters:**

* ``--p-region ITS1``: Targets the ITS1 sub-region specifically.
* ``--p-taxa F``: Restricts HMM searching to Fungi (selected by default).
* ``--o-trimmed``: Generates an output artifact of type ``SampleData[PairedEndSequencesWithQuality]``.

Example 2: trimming fungal single-end reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your dataset consists of single-end reads use the ``trim-single`` command. 
Using this for merged paired reads is not recommended by DADA2. It uses the information from both reads in its error modeling.

.. code-block:: bash

   qiime itsxpress trim-single \
     --i-per-sample-sequences single-end-demux.qza \
     --p-region ITS2 \
     --p-taxa F \
     --p-threads 4 \
     --o-trimmed trimmed-single.qza \
     --verbose

**Key Parameters:**

* ``--i-per-sample-sequences``: Accepts a ``SampleData[SequencesWithQuality]`` artifact containing single-end data.
* ``--p-region ITS2``: Targets the ITS2 sub-region.

Example 3: trimming PacBio CCS long reads (single-end)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


PacBio Circular Consensus Sequencing (CCS) generates long, high-fidelity reads that can span the 
entire ITS region (ITS1, 5.8S, and ITS2). To handle these long reads correctly, you must pass the ``--p-trim-ccs`` flag to the ``trim-single`` command.
This feature orients input sequence reads against a universal reference database with Vsearch. It outputs trimmed reads with synthetic 
dummy primers on the ends (forward: `GACAGGTACAAGAAGGA`, reverse complement of reverse: `ACTGGAGACTGGGTTAA`) 
with uniform Phred quality scores of 93 (`~`) to satisfy downstream QIMME2 DADA2  plugin requirements.
At completion, ITSxpress provides a log message displaying the example `qiime dada2 denoise-ccs` arguments for denoising.


.. code-block:: bash

   qiime itsxpress trim-single \
     --i-per-sample-sequences pacbio-ccs-demux.qza \
     --p-region ALL \
     --p-taxa F \
     --p-trim-ccs \
     --p-threads 4 \
     --o-trimmed trimmed-pacbio-ccs.qza \
     --verbose

**Key Parameters:**

* ``--p-region ALL``: Instructs the plugin to extract the full ITS composite sequence (ITS1 + 5.8S + ITS2).
* ``--p-trim-ccs``: Explicitly enables the specialized trimming mode optimized for PacBio CCS reads.


Stand Alone Command Line Examples
---------------------------------

Use case 1: Trimming the ITS2 region from a fungal amplicon sequencing dataset with
forward and reverse gzipped FASTQ files using two cpu threads. Return a single merged file for use in Deblur.

.. code-block:: bash

    itsxpress --fastq r1.fastq.gz --fastq2 r2.fastq.gz --region ITS2 \
    --taxa Fungi --log logfile.txt --outfile trimmed_reads.fastq.gz --threads 2

ITSxpress can take uncompressed, gzipped or zstd compressed FASTQ files and it can write uncompressed, gzipped or
zstd compressed FASTQ files. It expects FASTQ files to end in: .fq, .fastq, .fq.gz, fastq.gz, .fq.zst or fastq.zst.

Use case 2: Trimming the ITS2 region from a fungal amplicon sequencing dataset with
forward and reverse gzipped FASTQ files using two cpu threads. Return forward
and reverse read files for use in Dada2.

.. code-block:: bash

    itsxpress --fastq r1.fastq.gz --fastq2 r2.fastq.gz --region ITS2 \
    --taxa Fungi --log logfile.txt --outfile trimmed_reads_r1.fastq.gz --outfile2 trimmed_reads_r2.fastq.gz --threads 2

ITSxpress can take uncompressed, gzipped or zstd compressed FASTQ files and it can write uncompressed, gzipped or
zstd compressed FASTQ files. It expects FASTQ files to end in: .fq, .fastq, .fq.gz, fastq.gz, .fq.zst or fastq.zst.


Use case 3: Trimming the ITS2 region from a fungal amplicon sequencing dataset with
an single-ended gzipped FASTQ files using two cpu threads.

.. code-block:: bash

    itsxpress --fastq single-end.fastq.gz --single_end --region ITS2 --taxa Fungi \
    --log logfile.txt --outfile trimmed_reads.fastq.gz --threads 2

Single ended data is less common and may come from a platform like Oxford Nanopore or PacBio with long single-end reads. 
Pre-merging paired-end reads is not recommended because the next step for most people is DADA2 which prefers paired-end reads.


ITS HMM Models
--------------------

The HMM models for identifying ITS regions were built by Henrik Nilsson at the University of Gothenburg. 
Originally, they were built for the ITSx software package; Dr. Nilsson has continued to update the models, but not on the same schedule as ITSx. 
To provide clarity about the reference data used, I have created a separate, versioned repository for the 
ITS sequence models: `ITS_HMMs on GitHub <https://github.com/USDA-ARS-GBRU/ITS_HMMs>`_ and assigned it a `DOI <https://doi.org/10.5281/zenodo.13285214>`_. 

ITSXpress currently uses ITS_HMMS version 2.0.0.


API Documentation
-----------------

The API documentation is available at https://itsxpress.readthedocs.io

License information
--------------------

This software is a work of the United States Department of Agriculture, Agricultural Research Service.
Under 17 U.S. Code § 105 the work is not copyrightable.  It is released under a Creative Commons CC0
public domain attribution.

