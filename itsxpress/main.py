#!/usr/bin/env python
"""ITSxpress: A python module to rapidly trim ITS amplicon sequences from Fastq files.
Author: Adam Rivers, USDA Agricultural Reseach Service

The internally transcribed spacer region is a region between highly conserved the small
subunit (SSU) of rRNA and the large subunit (LSU) of the rRNA. In Eukaryotes it contains
the 5.8s genes and two variable length spacer regions. In amplicon sequening studies it is
common practice to trim off the conserved (SSU, 5,8S or LSU) regions. Bengtsson-Palme
et al. (2013) published software the software package ITSx to do this.

ITSxpress is a high-speed implementation of the methods in	ITSx. It can process a typical
ITS amplicon sample with 100,000 read pairs in about 5 minutes, aproxamatly 100x faster.
It also trims fastq files rather than just fasta files.

Process:
	* Merges and error corrects reads using bbduk if reade are paired-end
	* Deduplicates reads using Vmatch to eliminate redundant hmm searches
	* Searches for conserved regions using the ITSx hmms, useing HMMsearch:
	  https://cryptogenomicon.org/2011/05/27/hmmscan-vs-hmmsearch-speed-the-numerology/
	* Parses everyting in python returning (optionally gzipped) fastq files.

Refernce:
	Johan Bengtsson-Palme, Vilmar Veldre, Martin Ryberg, Martin Hartmann, Sara Branco,
	Zheng Wang, Anna Godhe, Yann Bertrand, Pierre De Wit, Marisol Sanchez,
	Ingo Ebersberger, Kemal Sanli, Filipe de Souza, Erik Kristiansson, Kessy Abarenkov,
	K. Martin Eriksson, R. Henrik Nilsson. (2013). ITSx: Improved software detection
	and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other
	eukaryotes for use in environmental sequencing. Methods in Ecology and Evolution,
	4: 914-919, 2013 (DOI: 10.1111/2041-210X.12073)
"""

import gzip
import tempfile
import argparse
import subprocess
import logging
import time
import os
import shutil

from Bio import SeqIO

from itsxpress.definitions import ROOT_DIR, taxa_choices, taxa_dict

def _myparser():
	parser = argparse.ArgumentParser(description='ITSxpress: A python module to rapidly \
		trim ITS amplicon sequences from Fastq files.')
	parser.add_argument('--fastq', '-f', type=str, required=True,
						help='A .fastq, .fq, .fastq.gz or .fq.gz file. Interleaved or not.')
	parser.add_argument('--single_end', '-s', action='store_true', default=False,
						help='A flag to specify if the fastq file is inteleaved single-ended (not paired). Default is false.')
	parser.add_argument('--fastq2', '-f2', type=str, default=None,
						help='A .fastq, .fq, .fastq.gz or .fq.gz file. representing read 2 (optional)')
	parser.add_argument('--outfile', '-o', type=str, help="the trimmed Fastq file, if it \
						ends in 'gz' it will be gzipped")
	parser.add_argument('--tempdir', help='Specify the temp file directory', default=None)
	parser.add_argument('--keeptemp' ,help="Should intermediate files be kept?", action='store_true')
	parser.add_argument('--region', help='', choices=["ITS2", "ITS1", "ALL"], required=True)
	parser.add_argument('--taxa', help='Select the taxonomic group sequenced',
						choices=taxa_choices, default="Fungi")
	parser.add_argument('--log' ,help="Log file", default="ITSxpress.log")
	parser.add_argument('--threads' ,help="Number of processor threads to use", default="1")
	try:
		args = parser.parse_args()
		return args
	except:
		parser.print_help()
		raise SystemExit(0)




class ItsPosition:
	"""Class for ITS positional information derived from hmmserach domtable files.


	Args:
		domtable (str):	 the path locating the domtable file from HMMER 3 hmmsearch.
		region (str): The region of the ITS to extract choises: ["ITS1", "ITS2", "ALL"].

	Attributes:
		ddict (dict): A dictionary holding the scores and start and stop
			positions for the selected segment of each sequence.
			Example: {sample:{left:{score:31, pos:15}, right:{score:32, pos:354}}
		leftprefix (str): the left prefix to search for (set by type variable).
		rightprefix (str): the right prefix to search for (set by type variable).

	Todo:
		* Add additional ITS regions.

	"""
	
	def _left_score(self, sequence, score, to_pos):
		"""Evaluates left scores and positions from the new line of a domtable file and
			updates ddict if neccicary.
		
		Args:
			sequence (str): The name of the sequence.
			score (int): The bit score from HMMSearch.
			to_pos (int): the ending position of the left seqeunce.
		
		"""
		
		if "left" in self.ddict[sequence]:
			if score > self.ddict[sequence]["left"]["score"]:
				 self.ddict[sequence]["left"]["score"]=score
				 self.ddict[sequence]["left"]["pos"]=to_pos
		else:
			self.ddict[sequence]["left"]={}
			self.ddict[sequence]["left"]["score"]=score
			self.ddict[sequence]["left"]["pos"]=to_pos
		
		
	def _right_score(self, sequence, score, from_pos):
		"""Evaluates right scores and positions form the new line of a domtable file and
		updates ddict if neccicary.
	
		Args:
			sequence (str): The name of the sequence.
			score (int): The bit score from HMMSearch.
			from_pos (int): The beginning position of the right seqeunce.
	
		"""
		if "right" in self.ddict[sequence]:
			if score > self.ddict[sequence]["right"]["score"]:
				 self.ddict[sequence]["right"]["score"]=score
				 self.ddict[sequence]["right"]["pos"]=from_pos
		else:
			self.ddict[sequence]["right"]={}
			self.ddict[sequence]["right"]["score"]=score
			self.ddict[sequence]["right"]["pos"]=from_pos
		

	def parse(self):
		"""Parses dom table from HMMsearch.
		
		The dom table is parsed to record the start and stop position from the top scoring
		hmm mathces. This populates the ddict attribute containing the positions at
		which to trim each sequence.
		
		"""
		try:
			with open(self.domtable , 'r') as f:
				for line in f:
					if not line.startswith("#"):
						ll=line.split()
						sequence=ll[0]
						hmmprofile=ll[3]
						score=ll[7]
						from_pos=ll[19]
						to_pos=ll[20]
						if sequence not in self.ddict:
							self.ddict[sequence]={}
						if hmmprofile.startswith(self.leftprefix):
							self._left_score(sequence, score, to_pos)
						elif hmmprofile.startswith(self.rightprefix):
							self._right_score(sequence, score, from_pos)
		except Exception as e:
			logging.error("Exception occured when parsing HMMSearh results")
			raise e
			
	def __init__(self, domtable, region):
		self.domtable = domtable
		self.ddict = {}
		if region=="ITS2":
			self.leftprefix='3_'
			self.rightprefix='4_'
		elif region=="ITS1":
			self.leftprefix='1_'
			self.rightprefix='2_'
		elif region=="ALL":
			self.leftprefix='1_'
			self.rightprefix='4_'
		self.parse()
	
			
	def get_position(self, sequence):
		""" Retuns the start and stop positions for a given sequence
	
		Args:
			sequence (str): The name of the sequence.
		
		Returns:
			(tuple): (start position, end position) zero indexed
			
		Raises:
			KeyError: If input sequence is not present in dictionary (no ITS start or stop sites were found)
		
		"""
	
		try:
			if "left" in self.ddict[sequence]:
				start = int(self.ddict[sequence]["left"]["pos"])
			else:
				start = None
			if "right" in self.ddict[sequence]:
				stop =	int(self.ddict[sequence]["right"]["pos"]) - 1
			else:
				stop = None
			return(start, stop)
		except KeyError:
			logging.debug("No ITS stop or start sites were identified for sequence {}, skipping.".format(sequence))
			raise KeyError
						
class Dedup:
	"""A class to handle deduplicated sequence data.
	
	To speed processing Vmatch is used to remove duplicate amplicons so that the
	start ansd	stop sites are determened only once.
	
	Attributes:
		matchdict (dict): a dictionary of each sequence ID as a key and
			its representative sequence ID as a value {seq1:rep1, seq2:rep1, seq2:rep2}.
		uc_file (str): the location of the .uc file contianing matching information.
		rep_file (str): The location of the representative sequences file.
		seq_file (str): Teh location of the complete sequences file.
		
	
	"""
	
		
	def parse(self):
		"""
		Parse the uc data file to populate the matchcdict attribute.
		
		Rasies:
			SystemExit: Exits if uc file cannot be parsed
	
		"""
		try:
			with open(self.uc_file, 'r') as f:
				self.matchdict = {}
				for line in f:
					ll = line.split()
					datatype = ll[0]
					ref = ll[9]
					seq = ll[8]
					if datatype == 'S':
						self.matchdict[seq] = seq
					elif datatype == 'H':
						self.matchdict[seq] = ref
		except Exception as e:
			logging.exception("Could not parse the Vsearch '.uc' file.")
			raise e
	
	def __init__(self, uc_file, rep_file, seq_file):
		self.matchdict = None
		self.uc_file = uc_file
		self.rep_file = rep_file
		self.seq_file = seq_file
		self.parse()
	
		
	
	def _get_trimmed_seq_generator(self, seqgen, itspos):
		"""This function takes a Biopython SeqIO sequence generator of sequences, and
		returns a generator of trimmed sequecnes. Sequences where the ITS ends could
		not be determed are ommited.
		
		Args:
		
			seqgen (obj): A Biopython SeqIO generator of all input sequences
			ispos (obj): a itsxpress ItsPosition object
		
		Returns:
			(obj): A map object generator that yeilds filtered, trimmed sequence records.
			
		"""
		
	
		def _filterfunc(record):
			""" Filters records down to those that contain a valid ITS start and stop position
			
			Args:
				record (obj): a Biopython SeqRecord object
			
			Returns:
				bool: True if an ITS start and stop positions are present false otherwise
			  
			"""
			try:
				if record.id in self.matchdict:
					repseq = self.matchdict[record.id]
					start, stop = itspos.get_position(repseq)
					if start and stop:
						return True
				else:
					return False
			except KeyError:
				return False

		def map_func(record):
			"""Trims the record down to to correct region
			
			Args:
				record (obj): a Biopython SeqRecord object
			
			Returns:
				obj: a Biopython SeqRecord object trimmed to the ITS region
			"""
			repseq = self.matchdict[record.id]
			start, stop = itspos.get_position(repseq)
			return record[start:stop]
			
		filt = filter(_filterfunc, seqgen)
		return map(map_func, filt)

				
	def create_trimmed_seqs(self, outfile, gzipped, itspos):
		"""Creates a fastq file, optionally gzipped, with the reads trimmed to the
			selected region.
		
		Args:
			outfile (str): the file to write the sequences to.
			gzip (bool): Should the files be gzipped?
			itspos (object): an ItsPosition object
		
		Returns:
			str: name of the file written
		
		"""
		
		def _write_seqs():
			if gzipped:
				with gzip.open(outfile, 'wt') as g:
					SeqIO.write(seqs, g, "fastq")
			else:
				with open(outfile, 'w') as g:
					SeqIO.write(seqs, g, "fastq")
				
		if self.seq_file.endswith(".gz"):
			with gzip.open(self.seq_file, 'rt') as f:
				seqgen = SeqIO.parse(f, 'fastq')
				seqs = self._get_trimmed_seq_generator(seqgen, itspos)
				_write_seqs()

		else:
			with open(seld.seq_file, 'r') as f:
				seqgen = SeqIO.parse(f, 'fastq')
				seqs = self._get_trimmed_seq_generator(seqgen, itspos)
				_write_seqs()
				
				
class SeqSample:
	"""The class for processing sequence data into trimmed sequences.
	
	Attributes:
		tempdir (obj): A temporary directory object
		fastq (str): The path to the input fastq file
		uc_file (str): The path to the Vsearch uc mapping file
		rep_file: (str) the path to the representitive seequences fasta file created by Vsearch
		seq_file (str): the location of the fastq or fastq.gz sequence file used for analysis
		
	"""
	

	def __init__(self, fastq, tempdir=None):
		self.tempdir = tempfile.mkdtemp(prefix='itsxpress_', dir=tempdir)
		self.fastq = fastq
		self.uc_file = None
		self.rep_file =	None
		self.dom_file = None
		self.seq_file = None
		
	
	def _deduplicate(self, threads=1):
		"""Runs Vsearch dereplication to create a fasta file of nonredundant sequences.
		
		Args:
			threads (int or str):the number of prosessor threads to use
		
		"""
		try:
			self.uc_file=os.path.join(self.tempdir, 'uc.txt')
			self.rep_file=os.path.join(self.tempdir,'rep.fa')
			parameters = ["vsearch",
						  "--derep_fulllength",
						  self.seq_file,
						  "--output", self.rep_file,
						  "--uc", self.uc_file,
						  "--strand", "both",
						  "--threads", str(threads)]
			p2 = subprocess.run(parameters, stderr=subprocess.PIPE)
			logging.info(p2.stderr.decode('utf-8'))
			p2.check_returncode()
		except subprocess.CalledProcessError as e:
			logging.exception("Could not perform dereplication with Vsearch. Error from Vsearch was:\n {}".format(p2.stderr.decode('utf-8')))
			raise e
		except FileNotFoundError as f:
			logging.error("Vsearch was not found, make sure Vsearch is installed and executible")
			raise f

	def _search(self, hmmfile, threads=1):
		try:
			self.dom_file=os.path.join(self.tempdir, 'domtbl.txt')
			#Run hmmsearch
			parameters = ["hmmsearch",
			  			  "--domtblout",
						  self.dom_file,
						  "-T", "10",
						  "--cpu", str(threads),
						  "--tformat", "fasta",
						  "--F1", "1e-6",
						  "--F2", "1e-6",
						  "--F3", "1e-6",
						  hmmfile,
						  self.rep_file]
			p4 = subprocess.run(parameters,stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
			p4.check_returncode()
		except subprocess.CalledProcessError as e :
			logging.exception("Could not perform ITS identificaton with hmmserach. The error was:\n {}".format(p4.stderr.decode('utf-8')))
			raise e
		except FileNotFoundError as f:
			logging.error("hmmsearch was not found, make sure HMMER3 is installed and executible")
			raise f
			
class SeqSamplePairedInterleaved(SeqSample):
	"""SeqSample class extended to paired, interleaved format.
	
	"""
	def __init__(self, fastq, tempdir):
		SeqSample.__init__(self, fastq, tempdir)
		
	def _merge_reads(self, threads):
		try:
			seq_file = os.path.join(self.tempdir, 'seq.fq.gz')
			parameters = ['bbmerge.sh',
					  'in=' + self.fastq,
					  'out=' + seq_file,
					  't=' + str(threads)]
			p1 = subprocess.run(parameters, stderr=subprocess.PIPE)
			self.seq_file = seq_file
			p1.check_returncode()
			logging.info(p1.stderr.decode('utf-8'))
		except subprocess.CalledProcessError as e:
			logging.exception("could not perform read merging with BBmerge. Error from BBmerge was: \n  {}".format(p1.stderr.decode('utf-8')))
			raise e
		except FileNotFoundError as f:
			logging.error("BBmerge was not found, make sure BBmerge is executible")
			raise f

class SeqSamplePairedNotInterleaved(SeqSample):
	"""SeqSample class extended to paired, two fastq file format.
	
	"""
	def __init__(self, fastq, tempdir, fastq2):
		SeqSample.__init__(self, fastq, tempdir)
		self.fastq2 = fastq2
	
	def _merge_reads(self, threads):
		try:
			seq_file = os.path.join(self.tempdir, 'seq.fq.gz')
			parameters = ['bbmerge.sh',
					  'in=' + self.fastq,
					  'in2=' + self.fastq2,
					  'out=' + seq_file,
					  't=' + str(threads)]
			p1 = subprocess.run(parameters, stderr=subprocess.PIPE)
			self.seq_file = seq_file
			p1.check_returncode()
			logging.info(p1.stderr.decode('utf-8'))
		except subprocess.CalledProcessError as e:
			logging.exception("could not perform read merging with BBmerge. Error from BBmerge was: \n  {}".format(p1.stderr.decode('utf-8')))
			raise e
		except FileNotFoundError as f:
			logging.error("BBmerge was not found, make sure BBmerge is executible")
			raise f

class SeqSampleNotPaired(SeqSample):
	"""SeqSample class extended to unpaired format.
	
	"""

	def __init__(self, fastq, tempdir):
		SeqSample.__init__(self, fastq, tempdir)
		self.seq_file = self.fastq

	

def _is_paired(fastq, fastq2, single_end):
	"""Dertermines workflow based on file inputs.
	
	Args:
		
	"""
	if fastq and fastq2:
		paired_end = True
		interleaved = False
	elif single_end:
		paired_end = False
		interleaved = False
	else:
		paired_end=True
		interleaved=True
	return paired_end, interleaved

def _logger_setup(logfile):
	"""Set up logging to a logfile and the terminal standard out.
	
	Args:
		fastq (str): The path to a fastq or fastq.gz file
		fastq2 (str): The path to a fastq or fastq.gz file for the reverese sequences
	
	"""
	try:
		logging.basicConfig(level=logging.DEBUG,
						format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
						datefmt='%m-%d %H:%M',
						filename=logfile,
						filemode='w')
		# define a Handler which writes INFO messages or higher to the sys.stderr
		console = logging.StreamHandler()
		console.setLevel(logging.INFO)
		# set a format which is simpler for console use
		formatter = logging.Formatter('%(asctime)s: %(levelname)-8s %(message)s')
		# tell the handler to use this format
		console.setFormatter(formatter)
		# add the handler to the root logger
		logging.getLogger('').addHandler(console)
	except Exception as e:
		print("An error occured setting up logging")
		raise e

def _check_fastqs(fastq, fastq2=None):
	"""Verifies the input files are valid fastq or fastq.gz files.
	
	Args:
		fastq (str): The path to a fastq or fastq.gz file
		fastq2 (str): The path to a fastq or fastq.gz file for the reverese sequences
	
	Raises:
		SystemExit if invalid input sequences are found.
	"""
	try:
		parameters=['reformat.sh', 'in='+fastq]
		if fastq2:
			parameters.append('in2=' + fastq2)
		p1 = subprocess.run(parameters, stderr=subprocess.PIPE)
		p1.check_returncode()
	except subprocess.CalledProcessError as e:
		logging.error("There appears to be an issue with your input fastq or fastq.gz file(s).")
		raise e
	except FileNotFoundError as e:
		logging.error("BBtools was not found. check that the BBtools reformat.sh package is executible")
		raise e


def main():
	"""Run Complete ITS trimming workflow
	
	"""
	# Set up logging
	t0 = time.time()
	args = _myparser()
	_logger_setup(args.log)
	try:
		logging.info("Verifing the input sequecnes.")
		_check_fastqs(args.fastq, args.fastq2)
		# Parse input types
		paired_end, interleaved = _is_paired(args.fastq,args.fastq2, args.single_end)
		# create SeqSample objects and merge if needed
		if paired_end and interleaved:
			logging.info("Sequences are paired-end and interleaved. They will be merged using BBmerge.")
			sobj = SeqSamplePairedInterleaved(fastq=args.fastq, tempdir=args.tempdir)
			sobj._merge_reads(threads=args.threads)
		elif paired_end and not interleaved:
			logging.info("Sequences are paired-end in two files. They will be merged using BBmerge.")
			sobj = SeqSamplePairedNotInterleaved(fastq=args.fastq, fastq2=args.fastq2, tempdir=args.tempdir)
			sobj._merge_reads(threads=args.threads)
		elif not paired_end and not interleaved:
			logging.info("Sequences are assumed to be single-end.")
			sobj = SeqSampleNotPaired(fastq=args.fastq, tempdir=args.tempdir)
		logging.info("Temporary directory is: {}".format(sobj.tempdir))
		# Deduplicate
		logging.info("Unique sequences are being written to a temporary fasta file with Vsearch.")
		sobj._deduplicate(threads=args.threads)
		# HMMSearch for ITS regions
		logging.info("Searching for ITS start and stop sites using HMMSearch. This step takes a while.")
		hmmfile = os.path.join(ROOT_DIR,"ITSx_db","HMMs", taxa_dict[args.taxa])
		sobj._search(hmmfile=hmmfile, threads=args.threads)
		# Parse HMMseach output
		logging.info("Parsing HMM results.")
		its_pos = ItsPosition(domtable=sobj.dom_file, region=args.region)
		# Create deduplication object
		dedup_obj = Dedup(uc_file=sobj.uc_file, rep_file=sobj.rep_file, seq_file=sobj.seq_file)
		# Create trimmed sequences
		logging.info("Writing out sequences")
		if args.outfile.split('.')[-1] =='gz':
			dedup_obj.create_trimmed_seqs(args.outfile, gzipped=True, itspos=its_pos)
		else:
			dedup_obj.create_trimmed_seqs(args.outfile, gzipped=False, itspos=its_pos)
		t1 = time.time()
		fmttime = time.strftime("%H:%M:%S",time.gmtime(t1-t0))
		logging.info("ITSxpress ran in {}".format(fmttime))
	except Exception:
		logging.error("ITSXpress terminated with errors. see the log file fo details.")
		raise SystemExit(1)
	finally:
		try:
			if not args.keeptemp:
				shutil.rmtree(sobj.tempdir)
		except UnboundLocalError:
			pass
		except AttributeError:
			pass


if __name__ == '__main__':
	main()
