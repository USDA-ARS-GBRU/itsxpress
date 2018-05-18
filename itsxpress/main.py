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
	* Merges and error corrects reads using bbduk
	* Deduplicates reads using Vmatch to eliminate redundant hmm searches
	* Searches for conserved regions using the ITSx hmms, but it uses HMMsearch rather 
	  HMMscan which is much faster for large numbers of sequences:
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
	parser.add_argument('--region', help='', choices=["ITS2", "ITS1", "ALL"], required=True)
	parser.add_argument('--taxa', help='Select the taxonomic group sequenced', 
						choices=taxa_choices, default="Fungi")
	parser.add_argument('--log' ,help="Log file", default="ITSxpress.log")
	parser.add_argument('--threads' ,help="Number of processor threads to use", default="1")
	args = parser.parse_args()
	return args



	
class ItsPosition:
	"""Class for ITS positional information derived from hmmserach domtable files.
	
	
	Args:
		domtable (str):	 the path locating the domtable file from HMMER 3 hmmsearch.
		type (str): The region of the ITS to extract choises: ["ITS2"].
	
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
		"""Parses domtable from HMMsearch.
		
		The dom table is parsed to record the start and stop position from the top scoring
		hmm mathces. This results in the ddict attribute containging the positions at 
		which to trim each sequence.
			
			
		""" 
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

	def __init__(self, domtable, type):
		self.domtable = domtable
		self.ddict = {}
		if type=="ITS2":
			self.leftprefix='3_'
			self.rightprefix='4_'
		elif type=="ITS1":
			self.leftprefix='1_'
			self.rightprefix='2_'
		elif type=="ALL":
			self.leftprefix='1_'
			self.rightprefix='4_'
		self.parse()
	
			
	def get_position(self, sequence):
		""" Retuns the start and stop positions for a given sequence
	
		Args:
			sequence (str): The name of the sequence.
		
		Returns:
			(tuple): (start position, end position) zero indexed
		
		"""
	
		try:
			if "left" in self.ddict[sequence]:
				start = int(self.ddict[sequence]["left"]["pos"]) - 1
			else:
				start = None
			if right in self.ddict[sequence]:
				stop =	int(self.ddict[sequence]["right"]["pos"]) - 1
			else:
				stop = None
			return(start, stop)
		except KeyError:
			logging.error("Could not return position for sequence {}.".format(sequence))

						
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
		Parse the uc data file to populat the matchcdict attribute.
	
		"""
	
		with open(self.uc_file, 'r') as f:
			self.matchdict = {}
			for line in f:
				ll = line.split()
				type = ll[0]
				ref = ll[9]
				seq = ll[8]
				if type == 'S':
					self.matchdict[seq] = seq
				elif type == 'H':
					self.matchdict[seq] = ref
	
	def __init__(self, uc_file, rep_file, seq_file):
		self.matchdict = None
		self.uc_file = uc_file
		self.rep_file = rep_file
		self.seq_file = seq_file
		self.parse()
	
	
	def _get_trimmed_seq_list(self, handle, itspos):
		seqs = []
		for record in SeqIO.parse(handle , "fastq"):
			try:
				repseq = self.matchdict[record.id]
				start, stop = itspos.get_position(repseq)
				seqs.append(record[start:stop])
			except:
				logging.warn("Could not parse {}, continuing.".format(record.id))
				continue
		return seqs
		
	def _get_trimmed_seq_generator(self, handle, itspos):
		pass
		seqgen = SeqIO.parse(handle , "fastq")
		for record in seqgen:
			if record.id in self.matchdict:
				repseq = self.matchdict[record.id]
			else: 
				continue
			start, stop = itspos.get_position(repseq)
		
		
				
	def create_trimmed_seqs(self, outfile, gzipped, itspos):
		"""Creates a fastq file, optionally gzipped, with the reads trimmed to the 
			selected region.
		
		Args:
			outfile (str): the file to write the sequences to.
			gzip (bool): Should the files be gzipped?
			itspos (object): an ItsPosition object
		
		Returns:
			(str): name of the file written
		
		"""
		if self.seq_file.endswith(".gz"):
			with gzip.open(self.seq_file, 'rt') as f:
				seqs = self._get_trimmed_seq_list(f,itspos)
		else:
			with open(seld.seq_file, 'r') as f:
				seqs = self._get_trimmed_seq_list(f, itspos)
		if gzipped:
			with gzip.open(outfile, 'wt') as g:
				SeqIO.write(seqs, g, "fastq")
		else:
			with open(outfile, 'w') as g:
				SeqIO.write(seqs, g, "fastq")
				
				
class SeqSample:
	"""The class for processing sequence data into trimmed sequences.
	
	Attributes:
		seq_file (str): the Location of the fastq or fastq.gz sequence file for the object
	
	"""
	

	def __init__(self, fastq, tempdir=None):
		self.tempdir = tempfile.mkdtemp(prefix='itsxpress_', dir=tempdir)
		self.fastq = fastq
		self.uc_file = None
		self.rep_file =	None
		self.dom_file = None
		self.seq_file = None
		
	
	def _deduplicate(self, threads=1):
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
			print(parameters)
			p2 = subprocess.run(parameters, stderr=subprocess.PIPE)
			logging.info(p2.stderr.decode('utf-8'))
			p2.check_returncode()
		except subprocess.CalledProcessError:
			logging.info("Could not perform dereplication with Vsearch")
			logging.info(p2.stderr.decode('utf-8'))
	
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
		except subprocess.CalledProcessError:
			logging.error("Could not perform ITS identificaton with hmmserach")
			logging.error(p3.stderr.decode('utf-8'))
			
			
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
			logging.info(p1.stderr.decode('utf-8'))
		except:
			logging.error("could not perform read merging with bbmerge")
			logging.error(p1.stderr.decode('utf-8'))

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
			logging.info(p1.stderr.decode('utf-8'))
		except:
			logging.error("could not perform read merging with bbmerge")
			logging.error(p1.stderr.decode('utf-8'))

class SeqSampleNotPaired(SeqSample):
	"""SeqSample class extended to unpaired format.
	
	"""

	def __init__(self, *args):
		SeqSample.__init__(self, *args)
		self.seq_file = self.fastq

	

def _is_paired(fastq, fastq2, single_end):
	"""dertermines workflow based on fill inputs
	
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
	

def main():
	"""Run Complete ITS trimming workflow
	
	"""
	#Set up logging
	args = _myparser()
	logging.basicConfig(filename=args.log, level=logging.INFO)
	# Parse input types
	paired_end, interleaved = _is_paired(args.fastq,args.fastq2, args.single_end)
	# create SeqSample objects and merge if needed
	if paired_end and interleaved:
		sobj = SeqSamplePairedInterleaved(fastq=args.fastq, tempdir=args.tempdir)
		sobj._merge_reads(threads=args.threads)
	elif paired_end and not interleaved:
		sobj = SeqSamplePairedNotInterleaved(fastq=args.fastq, fastq2=args.fastq2, tempdir=args.tempdir)
		sobj._merge_reads(threads=args.threads)
	elif not paired_end and not interleaved:
		sobj = SeqSampleNotPaired(fastq=args.fastq, tempdir=args.tempdir)
	#Deduplicate
	sobj._deduplicate(threads=args.threads)
	#HMMSearch for ITS regions
	hmmfile = os.path.join(ROOT_DIR,"ITSx_db","HMMs", taxa_dict[args.taxa])
	print("hmmfilepath = {}".format(hmmfile))
	sobj._search(hmmfile=hmmfile, threads=args.threads)
	# Parse HMMseach output
	its_pos = ItsPosition(domtable=sobj.dom_file, type=args.region)
	# Create deduplication object
	print("sobj uc_file is {}".format(sobj.uc_file))
	dedup_obj = Dedup(uc_file=sobj.uc_file, rep_file=sobj.rep_file, seq_file=sobj.seq_file)
	# create trimmed sequences
	if args.outfile.split('.')[-1] =='gz':
		dedup_obj.create_trimmed_seqs(args.outfile, gzipped=True, itspos=its_pos)
	else:
		dedup_obj.create_trimmed_seqs(args.outfile, gzipped=False, itspos=its_pos)
	#shutil.rmtree(sobj.tempdir)


if __name__ == '__main__':
	main()
