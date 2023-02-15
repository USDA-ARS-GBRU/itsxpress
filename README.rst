ITSxpress: Software to rapidly trim  the Internally transcribed spacer (ITS) region of FASTQ files
==================================================================================================
.. image:: https://travis-ci.org/USDA-ARS-GBRU/itsxpress.svg?branch=master
    :target: https://travis-ci.org/USDA-ARS-GBRU/itsxpress

.. image:: https://readthedocs.org/projects/itsxpress/badge/?version=latest
    :target: https://itsxpress.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://codecov.io/gh/USDA-ARS-GBRU/itsxpress/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/USDA-ARS-GBRU/itsxpress

.. image:: https://api.codacy.com/project/badge/Grade/7e2a4c97cde74bccb3e84b706d7a2aa5
  :target: https://www.codacy.com/app/GBRU/itsxpress?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=USDA-ARS-GBRU/itsxpress&amp;utm_campaign=Badge_Grade

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

#####

**This is the end of life version 1 of q2_itsxpress and the command line version of ITSxpress.
The new version 2 of ITSxpress, has the Qiime2 plugin built in with command line version of ITSxpress. See 
master branch of ITSxpress.**

#####

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

ITSxpress is also available as a `QIIME2 Plugin`_

.. _`Bengtsson-Palme et al. (2013)`: https://doi.org/10.1111/2041-210X.12073
.. _ITSx: http://microbiology.se/software/itsx/
.. _OTUs: https://doi.org/10.1038/ismej.2017.119
.. _`QIIME2 Plugin`: https://github.com/USDA-ARS-GBRU/q2_itsxpress


Installation
-------------

This is the installation of the final iteration of ITSxpress version 1: (BBmap is no longer used in ITSxpress version 2):

	- This version should primarily be used for reproducability with other datasets, which used ITSxpress =<1.8.1
	- The new version 2 is compatible with the newer versions of Qiime2

ITSxpress can be installed from:

1. Bioconda: (preferred method because it handles dependencies):

.. code-block:: bash

    conda install -c bioconda itsxpress==1.8.1

2. Pip: https://pypi.org/project/itsxpress/:

.. code-block:: bash

    conda install -y -c conda-forge hmmer==3.1b2
    conda install -y -c bioconda bbmap==38.69
    conda install -y -c bioconda vsearch==2.21.1
    pip install itsxpress


3. The Github repository: https://github.com/USDA-ARS-GBRU/itsxpress

.. code-block:: bash

    git clone https://github.com/USDA-ARS-GBRU/itsxpress.git


Dependencies
-------------
The software requires Vsearch, BBtools, Hmmer >= 3.1b and Biopython. Bioconda
takes care of this for you so it is the preferred installation method.


Usage
---------

-h, --help            	Show this help message and exit.

--fastq 				A ``.fastq``, ``.fq``, ``.fastq.gz`` or ``.fq.gz`` file. Interleaved
                        	or not. Required.

--single_end 			A flag to specify that the fastq file is single-ended (not paired).
                        	single-ended (not paired). Default is false.

--fastq2 				A ``.fastq``, ``.fq``, ``.fastq.gz`` or ``.fq.gz`` file representing read 2 if present, optional.

--outfile				The trimmed FASTQ file, if it ends in ``gz`` it will be gzipped.

--outfile2			The trimmed FASTQ read 2 file, if it ends in ``gz`` it will be gzipped. If used, reads will be retuned as unmerged pairs rather than than merged.

--tempdir				Specify the temp file directory. Default is None.

--keeptemp				Should intermediate files be kept? Default is false.

--region 				Options : {ITS2, ITS1, ALL}

--taxa					Select the taxonomic group sequenced: {Alveolata, Bryophyta,
						Bacillariophyta, Amoebozoa, Euglenozoa, Fungi, Chlorophyta,
						Rhodophyta, Phaeophyceae, Marchantiophyta, Metazoa,
						Oomycota, Haptophyceae, Raphidophyceae, Rhizaria, Synurophyceae,
						Tracheophyta, Eustigmatophyceae, All}. Default Fungi.

--cluster_id            The percent identity for clustering reads range [0.99-1.0], set to 1
                        for exact de-replication. Default 1.0.

--log		          	Log file. Default is ITSxpress.log.

--threads		     	Number of processor threads to use. Default is 1.

--reversed_primers  Primers are in reverse orientation as in Taylor et al. 2016,
                    DOI:10.1128/AEM.02576-16. If selected ITSxpress returns
                    trimmed reads flipped to the forward orientation



Examples
---------

Use case 1: Trimming the ITS2 region from a fungal amplicon sequencing dataset with
forward and reverse gzipped FASTQ files using two cpu threads. Return a single merged file for use in Deblur.

.. code-block:: bash

    itsxpress --fastq r1.fastq.gz --fastq2 r2.fastq.gz --region ITS2 \
    --taxa Fungi --log logfile.txt --outfile trimmed_reads.fastq.gz --threads 2

ITSxpress can take gzipped or un-gzipped FASTQ files and it can write gzipped or
un-gzipped FASTQ files. It expects FASTQ files to end in: .fq, .fastq, .fq.gz or fastq.gz.

Use case 2: Trimming the ITS2 region from a fungal amplicon sequencing dataset with
forward and reverse gzipped FASTQ files using two cpu threads. Return a forward
and reverse read files  for use in Dada2.

.. code-block:: bash

    itsxpress --fastq r1.fastq.gz --fastq2 r2.fastq.gz --region ITS2 \
    --taxa Fungi --log logfile.txt --outfile trimmed_reads.fastq.gz --threads 2

ITSxpress can take gzipped or un-gzipped FASTQ files and it can write gzipped or
un-gzipped FASTQ files. It expects FASTQ files to end in: .fq, .fastq, .fq.gz or fastq.gz.


Use case 3: Trimming the ITS2 region from a fungal amplicon sequencing dataset with
an interleaved gzipped FASTQ files using two cpu threads. Return a single merged file for use in Deblur.

.. code-block:: bash

    itsxpress --fastq interleaved.fastq.gz  --region ITS2 --taxa Fungi \
    --log logfile.txt --outfile trimmed_reads.fastq.gz --threads 2


Use case 4: Trimming the ITS2 region from a fungal amplicon sequencing dataset with
an single-ended gzipped FASTQ files using two cpu threads.

.. code-block:: bash

    itsxpress --fastq single-end.fastq.gz --single_end --region ITS2 --taxa Fungi \
    --log logfile.txt --outfile trimmed_reads.fastq.gz --threads 2

Single ended data is less common and may come from a dataset where the reads have already
been merged.

Use case 5: Trimming the ITS1 region from a Alveolata amplicon sequencing dataset with
an interleaved gzipped FASTQ files using 8 cpu threads.

.. code-block:: bash

    itsxpress --fastq interleaved.fastq.gz --region ITS1 --taxa Alveolata \
    --log logfile.txt --outfile trimmed_reads.fastq.gz --threads 8


License information
--------------------
This software is a work of the United States Department of Agriculture,
Agricultural Research Service and is released under a Creative Commons CC0
public domain attribution.
