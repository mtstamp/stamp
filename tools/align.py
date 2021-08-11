#!/usr/bin/python
#############################################################################
#Package: stamp
#Updated by Yiqin Wang
#Usage:
#
#
##############################################################################

#the path to the executables
samtools = "samtools-1.6" 
bwa ="bwa-0.7.17"
bamleftalign="bamleftalign-1.1.0"

import sys, os
import re

from argparse import ArgumentParser
from collections import defaultdict
from itertools import chain

#basic file readers and parsers
from FileIO import FileIO, execute, pipe_input, pipe_output, read_arguments

#sequencing file reads and parsers
from SeqIO import phred, prob, trim, flip, FastaIO, FastqIO, SAMFlag, SAMCigar

from math import exp, log, log10, factorial

import scipy.optimize as optimize

#debug
debug = False
#logging handler
LOG = None

#change the name of the "chr" function to "char" 
char = chr

#global parameters
#the maximum shift allowed for the mapped location of the captured fragments
align_shift_max = 20

#insert size for ndna; obsolete
insert_size_min = 300
insert_size_max = 500

#the maximum proportion of soft clip size in the read alignment
soft_clip_ratio = 0.4

#the minimum major allele count and ratio
major_allele_ratio = 0.6
major_allele_count = 10

#the total size of mtDNA (rCRS)
mtdna_len = 16569

#a ternary tree structure implemented for fast searching barcodes
class TernaryTree:
	"""
	a ternary tree structure implemented for fast searching barcodes
	
	Attributes
	----------
	nalphabet: the number of possible choices for each letter
	nletter: the number of possible letters of a string (use fixed-length string)
	tree: the root node of the ternary tree
	
	"""
	def __init__(self, nalphabet, nletter):
		""" 
		the __init__ method
		
		Arguments
		----------
		see class Attributes 
		"""
		self.nalphabet = nalphabet
		self.nletter = nletter
		self.tree = [None, ]*self.nalphabet
	
	@staticmethod
	def __search__(string, offset, nalphabet, node):
		""" 
		search a fixed-length string starting from a node of the tree
		
		Arguments
		----------
		offset: starting position in the string
		nalphabet: number of possible choices for each letter
		node: starting node
		
		Returns
		----------
		a list of matched nodes
		
		"""
		ret_all = []
		for i in range(offset, len(string)):
			letter = string[i]
			if (letter >= nalphabet):
				#treat as a wildcard, matching all possible alphabets
				for j in range(nalphabet):
					node = node[j]
					if (node):
						#recursively search the nodes at the next level 
						ret = TernaryTree.__search__(string, i+1, nalphabet, node)
						if (ret):
							ret_all.extend(ret)
				return ret_all
			else:
				node = node[letter]
			if (node is None):
				return None
		return [node,]
	
	@staticmethod
	def __insert__(string, offset, nalphabet, nletter, node):
		""" 
		insert a fixed-length string into a tree starting from a node of the tree
		
		Arguments
		----------
		offset: starting position in the string
		nalphabet: number of possible choices for each letter
		nletter: number of possible letters in the string
		node: starting node
		
		Returns
		----------
		a list of nodes inserted
		
		"""
		ret_all = []
		for i in range(offset, len(string)):
			letter = string[i]
			if (letter >= nalphabet):
				#treat as a wildcard, matching all possible alphabets
				if (i < nletter - 1):
					for j in range(nalphabet):
						if (node[j] is None):
							node[j] = [None,]*nalphabet
						ret_all.extend(TernaryTree.__insert__(string, i+1, nalphabet, nletter, node[j]))
				else:
					#the final letter
					for j in range(nalphabet):
						if (node[j] is None):
							node[j] = []
						ret_all.append(node[j])
				return ret_all
			else:
				if (node[letter] is None):
					if (i < nletter - 1):
						node[letter] = [None,]*nalphabet
					else:
						#the final letter
						node[letter] = []
				node = node[letter]
		ret_all.append(node)
		return ret_all
	
	def search(self, string):
		#search for a fixed-length string from the root
		assert len(string) == self.nletter
		node = self.tree
		return TernaryTree.__search__(string, 0, self.nalphabet, self.tree)
		
	def insert(self, string):
		#insert a  fixed-length string from the root
		assert len(string) == self.nletter
		node = self.tree
		return TernaryTree.__insert__(string, 0, self.nalphabet, self.nletter, self.tree)

#compute the hamming distance between two aligned strings with wildcards
def string_dist(s1, s2, wildcard = "N", letters = ["A","T","C","G"]):
	"""
	compute the hamming distance between two aligned strings with wildcards
	
	Arguments
	----------
	s1: the first string
	s2: the second string
	wildcard: the letter of the wildcard 
	
	s1 and s2 should be of the same length and their positions should be aligned.
	A wildcard does not count into the the calculation of the distance.
	
	Returns
	----------
	int distance
	
	"""
	dist = 0
	for i, j in zip(s1, s2):
		if (i != wildcard and j != wildcard and i in letters and j in letters and i!=j):
			dist += 1
	return dist

#mask low-quality letters within a string

def string_mask(s, q, qmin):
	"""
	mask low-quality letters within a string
	
	Arguments
	----------
	s: the original string
	q: the base quality 
	qmin: the minimum quality score allowed
	
	Returns
	----------
	the revised string with all low-quality letters (<qmin) being changed to "N"
	
	"""
	#using "N"
	return "".join([(i if j >= qmin else "N") for i,j in zip(s, q)])

#read probe file

def retrieve_probes(probe_file):
	"""
	read probe file
	
	Arguments
	----------
	probe_file: the path to the probe file
	
	Returns
	----------
	a list of probes
	each item is a tuple of
	name, chromosome, start position, end position, length of the segment, \
	strand of R1 read, strand of R2 read, probe sequence of R1, probe sequence of R2, length of barcode
		
	"""
	probe = []
	n = 0
	with open(probe_file) as fh:
		for line in fh:
			line = line.rstrip("\r\n")
			n += 1
			if (line):
				info = line.split("\t")
				if (len(info) < 8):
					raise ValueError("Unknown probe information %s at line %d" % (line, n))
				name, chr, start, end, r1_strand, r2_strand, r1_probe, r2_probe, barcode_len = info[:9]
				try:
					start = int(start)+1
					end = int(end)+1
				except:
					start = end = 0
				if (end > start):
					if (r1_strand == "+"):
						#ligation arm ==> <== extension arm
						start += len(r1_probe)
						end -= len(r2_probe)
						#remove the gap-filling and extension start sites due to high error rates
						#start += 1
						#end -= 1 
					else:
						#extension arm ==> <== ligation arm
						start += len(r2_probe)
						end -= len(r1_probe)
						#remove the gap-filling and extension start sites due to high error rates
						#start += 1 
						#end -= 1
					insert = (end-start)+1
				else:
					start = end = insert = 0
				r1_probe = flip(r1_probe)
				barcode_len = int(barcode_len)
				probe.append((name, chr, start, end, insert, r1_strand, r2_strand, r1_probe, r2_probe, barcode_len))
	return probe

#parse allele information

def parse_allele(string):
	"""
	parse allele information
	
	Arguments
	----------
	string: format (position+allele)
	i.e. 102C; 102d; 102CT;

	Returns
	----------
	a tuple of an integer position and a string allele
	or a tuple of None, None if the format is incorrect. 
	"""
	string = string.strip()
	ret = re.match("(^\d+)", string)
	try:
		pos = ret.group()
		allele = string[len(pos):]
		if (not allele):
			return int(pos), None
	except:
		#print >>sys.stderr, "Unknown allele pattern: " + string
		return None, None
	return int(pos), allele

#hamming distance between two sequences

def sequence_dist(s1, s2):
	"""
	the hamming distance between two sequences
	
	Arguments
	----------
	s1: the first sequence 
	s2: the second sequence
	Both s1 and s2 are lists of tuples with the integer position and the allele
	s1 and s2 can be the returns of parse_allele(...)
	
	
	Returns
	----------
	int distance
	Ns are counted as missing alleles not mismatches.
	 
	"""
	#only polymorphisms recorded in s1 and s2 {pos:allele, ..} 
	missing = 0
	match = 0
	for p, a1 in s1.items():
		a2 = s2.get(p)
		if (a2 is not None):
			if (a2 == a1):
				match += 1
		else:
			missing += 1
	#unique s2 positions + unmatched positions + s1 unique positions
	return len(s2) - match + missing

#read hsd file for alternate mtdna and numt sequences
def retrieve_alternate_sequences(hsd_file, offset = 0):
	"""
	read hsd file for alternate mtdna and numt sequences
	
	Arguments
	----------
	hsd_file: the path to the hsd file 
	offset: int position offset, used to adjust mtdna positions if using a revised version of rCRS
	
	Returns
	----------
	a dict with keys representing names of the sequences
	Values are lists with the start and end of the sequences and alleles in it.
	Alleles are dicts of int positions and nucleotides {position:nucleotide,...}
	"""
	#polymorphisms should be annotated according to rCRS
	#head information
	#SampleId(NameID)\tRange\tHaplogroup(MappingID)\tPolymorphisms (delimited with tabs)
	fh = FileIO(hsd_file, "r")
	sequences = {}
	for line in fh:
		line = line.rstrip("\r\n").split("\t")
		if (len(line) < 4):
			continue
		if (line[1] == "Range"):
			#skip head line
			continue
		id = line[0]
		#assuming range is start-end
		try:
			start, end = line[1].split("-")
			start = int(start)+offset
			end = int(end.replace(";",""))+offset
		except:
			raise ValueError("Unknown mtDNA range: %s" % line[1])
		if (start > end):
			#split the sequence into two if it spans the D-loop region: start ==> mtdna_length+offset and 1==>end+offset
			range = [(id + "_1", start, mtdna_len+offset),(id + "_2", 1, end)]
		else:
			range = [(id, max(start,1), min(end,mtdna_len+offset)),]
		for id, start, end in range:
			polymorphisms = {}
			assert id not in sequences, "Duplicated sequence id <%s> in file <%s>" % (id, hsd_file)
			for s in line[3:]:
				if (s.strip()):
					pos, allele = parse_allele(s)
					pos += offset
					if (pos >= 1 and pos <= mtdna_len + offset):
						if (pos >= start and pos <= end):
							polymorphisms[pos] = allele
						if (pos > mtdna_len):
							#shift the last <offset> bps to the head
							pos -= mtdna_len
							if (pos >= start and pos <= end):
								polymorphisms[pos] = allele
			sequences[id] = [start, end, polymorphisms]
			#print id, sequences[id]
	return sequences

#retrieve alternate sequences within a given region
def retrieve_alternate_sequences_by_region(sequences, start, end, overlap_min = 0.85):
	"""
	retrieve alternate sequences within a given region
	
	Arguments
	----------
	sequences: return of retrieve_alternate_sequences(...)
	start: int start of the region to retrieve
	end: int end of the region to retrieve
	overlap_min: the minimum overlapping percentage between the sequence and the region
	
	Returns
	----------
	a dict with keys representing names of the sequences
	Values are dicts of int positions and nucleotides {position:nucleotide,...} in the sequences
	
	"""
	region = {}
	for id in sequences:
		seq_start, seq_end, polymorphisms = sequences[id]
		#must overlap with >overlap_min of the given region 
		if (seq_start >= start + (end-start)*overlap_min):
			continue
		if (seq_end <= end - (end-start)*(1-overlap_min)):
			continue
		if ((min(seq_end, end) - max(seq_start, start))/float(end-start) < overlap_min):
			continue
		region[id] = dict([[pos, allele] for pos, allele in polymorphisms.items() if (pos >= start and pos <= end)])
	return region

#identify the sequence with the smallest hamming distance to the query sequence
def retrieve_alternate_sequences_by_similarity(sequences, query):
	"""
	identify the sequence with the smallest hamming distance to the query sequence
	
	Arguments
	----------
	sequences: return of retrieve_alternate_sequences(...)
	query: a sequence to compare with the ones in the sequences
	
	Hamming distances will be computed using sequence_dist(...) between the query and each of the sequences.
	
	Returns
	----------
	a tuple of the distance and the sequence id with the smallest hamming distance to the query sequence
	
	"""
	dist_min = sys.maxint
	match = None
	for id, polymorphisms in sequences.items():
		dist = sequence_dist(query, polymorphisms)
		if (dist < dist_min):
			dist_min = dist
			match = id
	return dist_min, match

#identify probe and barcode sequences within read1 and read2
def retrieve_tags(sample_name, fastq_r1, fastq_r2, probe, probe_reverse = False, outpath = None, output_reads = False, output_offtarget_reads = True, min_qual = 13, trim_qual = 0):
	"""
	identify probe and barcode sequences within read1 and read2
	
	Arguments
	----------
	sample_name: the name of the sample
	fastq_r1: the path to the R1 fastq file
	fastq_r2: the path to the R2 fastq file
	probe: probe information returned from retrieve_probes(...)
	probe_reverse: switch R1 and R2 reads
	outpath: the path to store the resulting files
	output_reads: just to check the files if turned off 
	output_offtarget_reads: write R1 and R2 pairs with incorrect pairs of arm sequences to the output files; discard if turned off
	min_qual: minimum quality of barcode bases
	trim_qual: minium quality used to trim ends of R1 and R2
		
	Returns
	----------
	the paths to the processed R1 and R2 files
	
	Outputs
	----------
	${outpath}/${sample_name}_R1.fastq.gz
	${outpath}/${sample_name}_R2.fastq.gz
	${outpath}/${sample_name}.summary
	${outpath}/${sample_name}.barcode
		
	"""
	
	tag_start = 0 #skip arms
	tag_len = 50 #50bp to include barcode and arm
	read_len = 80 #the miminum read length after trimming R1 and R2
	probe_mismatch = 3 #the maximum mismatches in the probe arm
	barcode_dlen = 13 #the length of the barcode used
	barcode_minq = min_qual
	barcode_minq_char = char(barcode_minq+33)
	barcode_minlen = 9 #the minimum length of high-quality bases in the barcode
	barcode_index = [] #
	barcode_cache = []
	npr = len(probe)
	#for every combination of the probe arms
	for i in range(npr*npr):
		#a ternary tree for each
		barcode_index.append(TernaryTree(4, barcode_dlen))
		#store the search result for each barcode
		barcode_cache.append({})
	
	#nucleotide decoder
	nt = {"a":0,"A":0,"t":1,"T":1,"c":2,"C":2,"g":3,"G":3}
	
	#variables to store read count
	id = 0
	lowqual_read = 0
	lowqual_barcode = 0
	offtarget_read = 0
	mp_r1 = 0
	mp_r2 = [0,]*npr
	mp_num = [0,]*(npr*npr)
	
	if (not outpath):
		outpath = "."
	outname = outpath + os.path.sep + sample_name
	if (sample_name and output_reads):
		#zip the output files
		r1_output = FileIO(outname + "_R1.fastq.gz", "w")
		r2_output = FileIO(outname + "_R2.fastq.gz", "w")
		ret = (r1_output.name, r2_output.name)
	else:
		output_reads = False
		ret = ("","")
	
	if (trim_qual > 0):
		trim_limit = prob(trim_qual)
	else:
		trim_limit = None
	
	#for each read pairs of R1 and R2
	for r1_input, r2_input in zip(fastq_r1, fastq_r2):
		if (not probe_reverse):
			input = FastqIO(r1_input, r2_input)
		else:
			input = FastqIO(r2_input, r1_input)
		for r1, r2 in input:
			id += 1
			name_r1, read_r1, quality_r1 = r1
			name_r2, read_r2, quality_r2 = r2
			if (trim_limit is not None):
				l1 = trim(phred(quality_r1),trim_limit)
				l2 = trim(phred(quality_r2),trim_limit)
				if (l1 < read_len or l2 < read_len):
					lowqual_read += 1
					continue
				read_r1 = read_r1[:l1]; quality_r1 = quality_r1[:l1]
				read_r2 = read_r2[:l2]; quality_r2 = quality_r2[:l2]
			#mask low-quality bases in R1 and R2 with N
			tag_r1 = string_mask(read_r1[tag_start:(tag_start+tag_len)], quality_r1[tag_start:(tag_start+tag_len)], barcode_minq_char)
			tag_r2 = string_mask(read_r2[tag_start:(tag_start+tag_len)], quality_r2[tag_start:(tag_start+tag_len)], barcode_minq_char)
			
			dmin = len(read_r1)
			pid1 = -1
			for i, v in enumerate(probe):
				r1_probe, r2_probe, barcode_len = v[-3:]
				d = string_dist(tag_r1, r1_probe)
				if (d < dmin):
					#find the closest r1 arm
					r1_probe_len = len(r1_probe)
					dmin = d
					pid1 = i
			if (dmin <= probe_mismatch):
				#assign pid1 to r1 and then check r2
				r1_probe, r2_probe, barcode_len = probe[pid1][-3:]
				barcode_offset = -3
				dmin = len(read_r2)
				#for R2 consider the posibility of shifted barcode sequence caused by demultiplexing individual barcode
				for offset in range(-3,4):
					if (barcode_len+offset < 0):
						continue
					d = string_dist(tag_r2[barcode_len+offset:], r2_probe)
					#compute distance between r2 arm and probe arm
					if (d < dmin):
						dmin = d
						barcode_offset = offset
				if (dmin >= probe_mismatch):
					barcode_offset = 0
					#find the closest r2 arm
					for pid in range(npr):
						r1_probe, r2_probe, barcode_len = probe[pid][-3:]
						d = string_dist(tag_r2[barcode_len:], r2_probe)
						if (d <= probe_mismatch):
							r2_probe_len = len(r2_probe)
							pid2 = pid
							break
					else:
						#cannot determine r2 probe
						mp_r2[pid1] += 1
						pid2 = -1
				else:
					pid2 = pid1
				r2_probe_len = len(r2_probe)
				if (pid2 >= 0):
					if (barcode_len):
						barcode_len += barcode_offset
						#concatenate N to let the barcode have a fixed length (default:13)
						if (barcode_len < barcode_dlen):
							barcode = "N"*(barcode_dlen-barcode_len) + tag_r2[:barcode_len]
						else:
							#use the last barcode-dlen bases 
							barcode = tag_r2[(barcode_len-barcode_dlen):barcode_len]
					else:
						#assign a dummy barcode
						barcode = "N%d"%id
					nleft = len(barcode) - barcode.count("N")
					if (barcode_len and nleft < barcode_minlen):
						#determined based in barcode less than barcode_minlen; discard
						lowqual_barcode += 1
					else:
						#count the number of probe pairs
						mp_num[pid1*npr+pid2]+=1
						if (barcode_len):
							#update the barcode tree for the pid1 and pid2 probe pairs
							barcode_groups = barcode_index[pid1*npr+pid2].insert([nt.get(i,4) for i in barcode])
							#return the same barcode group; N is treated as a wildcard in the string
							nfam = 0
							for i in barcode_groups:
								i.append((id, barcode))
						else:
							#dummy barcode group
							barcode_groups = [[(id,barcode,)],]
						barcode_cache[pid1*npr+pid2][barcode] = barcode_groups
						if (output_reads):
							if (output_offtarget_reads or pid1 == pid2):
								#output trimmed reads
								r1_output.write(name_r1.split(" ")[0] + " XM:Z:1:%s:%s:%s\n" % (barcode,probe[pid1][0],probe[pid2][0]) + read_r1[(tag_start+r1_probe_len):] + "\n+\n" + quality_r1[(tag_start+r1_probe_len):] + "\n")
								r2_output.write(name_r2.split(" ")[0] + " XM:Z:2:%s:%s:%s\n" % (barcode,probe[pid1][0],probe[pid2][0]) + read_r2[(tag_start+barcode_len+r2_probe_len):] + "\n+\n" + quality_r2[(tag_start+barcode_len+r2_probe_len):] + "\n")
								#remove the gap-filling site at read1
								#r1_output.write(name_r1.split(" ")[0] + " XM:Z:1:%s:%s:%s\n" % (barcode,probe[pid1][0],probe[pid2][0]) + read_r1[(tag_start+r1_probe_len+1):] + "\n+\n" + quality_r1[(tag_start+r1_probe_len+1):] + "\n")
								#remove the extension-start site at read2
								#r2_output.write(name_r2.split(" ")[0] + " XM:Z:2:%s:%s:%s\n" % (barcode,probe[pid1][0],probe[pid2][0]) + read_r2[(tag_start+barcode_len+r2_probe_len+1):] + "\n+\n" + quality_r2[(tag_start+barcode_len+r2_probe_len+1):] + "\n")
							if (pid1 != pid2):
								offtarget_read += 1
			else:
				#cannot determine r1 probe
				mp_r1 += 1
			#if (id % 10000 == 0):
			#	print id-lowqual_read-lowqual_barcode-mp_r1-sum(mp_r2)-offtarget_read, id
			#	for pid in range(npr):
			#		print probe[pid][0], mp_num[pid*npr+pid]
	if (output_reads):
		#close file handles
		r1_output.close()
		r2_output.close()
		barcode_output = open(outname + ".barcode", "w")
	
	id_retained = id-lowqual_read-lowqual_barcode-mp_r1-sum(mp_r2) #read pairs retained
	print >>LOG, "total_reads", sample_name, id
	print >>LOG, "proper_reads", sample_name, id_retained-offtarget_read, offtarget_read
	print >>LOG, "improper_reads", sample_name, lowqual_read, lowqual_barcode, mp_r1, sum(mp_r2)
	
	barcode_hist_all = defaultdict(lambda: 0) #count the occurance of barcodes
	nfam_all = 0
	for pid1 in range(npr):
		for pid2 in range(npr):
			if (not output_offtarget_reads and pid1 != pid2):
				continue
			name1 = probe[pid1][0]
			name2 = probe[pid2][0]
			#print name1, name2
			barcode_hist = defaultdict(lambda: 0)
			for barcode, barcode_groups in barcode_cache[pid1*npr+pid2].items():
				groups = defaultdict(lambda: {})
				for group in barcode_groups:
					for i, j in group:
						groups[j][i] = 1
				nbarcode = len(groups[barcode])
				ntotal = sum([len(i) for i in groups.values()])
				other_barcodes = list(set([j for i, j in group if j != barcode for group in barcode_groups ]))
				if (len(other_barcodes) >= 2):
					dmin = -1
					for i in range(len(other_barcodes)):
						for j in range(i+1,len(other_barcodes)):
							dmin = max(dmin, string_dist(other_barcodes[i], other_barcodes[j]))
					if (dmin > 0):
						barcode_hist[1] += nbarcode
						if (output_reads):
							barcode_output.write("%s\t%s\t%s\t%s\t%d\n" % (name1, name2, barcode, "", -1))
					else:
						barcode_hist[ntotal] += nbarcode
						if (output_reads):
							#write barcode and similar barcodes with Ns
							barcode_output.write("%s\t%s\t%s\t%s\t%d\n" % (name1, name2, barcode, ",".join(other_barcodes), nbarcode))
				else:
					barcode_hist[ntotal] += nbarcode
					if (output_reads):
						barcode_output.write("%s\t%s\t%s\t%s\t%d\n" % (name1, name2, barcode, ",".join(other_barcodes), nbarcode))
			rfam = [(i, j, j/float(i)) for i, j in sorted(barcode_hist.items())]
			nfam = sum([i[2] for i in rfam])
			nfam_all += nfam
			#output summary of barcode counts for each probe pair
			for fam_size, read_num, fam_num in rfam:
				print >>LOG, "tag_barcode", sample_name, name1, name2, mp_num[pid1*npr+pid2], mp_r2[pid1], nfam, fam_size, read_num, fam_num
			for i, j in barcode_hist.items():
				barcode_hist_all[i] += j
	#output summary of barcode counts for all probe pairs
	for i, j in sorted(barcode_hist_all.items()):
		print >>LOG, "tag_barcode", sample_name, "all", "all", sum(mp_num), sum(mp_r2), nfam_all,i, j, j/float(i)
	if (output_reads):
		barcode_output.close()
	#ret = (r1_output.name, r2_output.name)
	return ret

def retrieve_XM_info(alignment, paired = True):
	"""
	extract the XM tag information from the alignment file
	
	Arguments
	----------
	sequences: one line of the sam/bam file
	paired: true if paired-end
	
	Returns
	----------
	a tuple of
	(1/2, barcode, probe_r1, probe_r2)
	or
	(barcode, probe_r1, probe_r2)
	"""
	for i in range(len(alignment)-1, 10, -1):
		if (alignment[i].startswith("XM:Z:")):
			info = alignment[i][5:].split(":")
			if (paired):
				if (len(info) != 4):
					raise ValueError("Unknown information column \"%s\"" % info)
			else:
				if (len(info) != 3):
					raise ValueError("Unknown information column \"%s\"" % info)
			return info
	return ""

def retrieve_allele_info(val, remover_n=False):
	"""
	extract allele information from val
	
	Arguments
	----------
	val: alleles "," separated
	remover_n: remove N alleles
	
	Returns
	----------
	a dict of int positions and alleles
	"""
	#allele format: num_alleles:pos1_a1,pos2_a2
	info = {}
	if (not val):
		return None
	val = val.split(",")
	try:
		n = int(val[0])
		if (n != len(val)-1):
			raise
		else:
			for i in range(1, n+1):
				pos, allele = parse_allele(val[i])
				if (pos is None):
					raise
				if (not remover_n or allele != "N"):
					info[pos] = allele
	except:
		raise ValueError("Unknown allele information <%s>" % val)
	return info

def retrieve_tag_info(alignment, tag, type):
	"""
	extract a certain tag from the alignment line
	
	Arguments
	----------
	alignment: one line of the sam/bam file
	tag: a tag string
	type a tag type string
	
	Returns
	----------
	a string of the tag value
	
	"""
	for i in range(len(alignment)-1, 10, -1):
		key = "%s:%s:"%(tag, type)
		if (alignment[i].startswith(key)):
			return alignment[i][len(key):]
	return None

def retrieve_tags_info(alignment):
	"""
	extract all tag information from the alignment line
	
	Arguments
	----------
	alignment: one line of the sam/bam file
	
	Returns
	----------
	a dict of all tags and their tag values
	
	"""
	tags = {}
	for i in range(len(alignment)-1, 10, -1):
		val = alignment[i].split(":", 2)
		if (len(val) == 3):
			tags[val[0]] = val[1:]
	return tags

def align_reads_to_reference(fastq_r1, fastq_r2, reference, probe_info):
	"""
	a wrapper function to call bwa mem
	
	Arguments
	----------
	fastq_r1: the path to the R1 fastq file
	fastq_r2: the path to the R2 fastq file
	reference: the path to the reference genome
	probe_info: probe information obtained from retrieve_probes(...)
	
	Returns
	----------
	an iterator of paired-end reads in the alignment file
	
	"""
	if (not os.path.exists(fastq_r1)):
		raise IOError("Cannot open %s" % fastq_r1)
	if (not os.path.exists(fastq_r2)):
		raise IOError("Cannot open %s" % fastq_r2)
	if (not os.path.exists(reference)):
		raise IOError("Cannot open %s" % reference)
	
	#mapping
	#map filtered reads to the reference genome
	#-L 100,5 forces matching after the clipped probe sequences
	cmd = "%s mem -L 100,5 -M %s -C %s %s" % (bwa, reference, fastq_r1, fastq_r2)
	pf = pipe_output(cmd)
	head = []
	alignments = []
	read_name = False
	n = 0
	for line in pf.stdout:
		n += 1
		line = line.rstrip("\r\n")
		if (not line):
			continue
		if (head is not None):
			if (line[0] == "@"):
				if (line.startswith("@PG")):
					#replace the actual file links with the variable names 
					reference_file = reference[reference.rfind(os.path.sep)+1:] 
					line = line.replace(reference, "$reference_path/"+reference_file)
					line = line.replace(fastq_r1, "$fastq_r1_file")
					line = line.replace(fastq_r2, "$fastq_r2_file")
				head.append(line)
				continue
			else:
				yield head
				head = None
		line = line.split("\t")
		info = retrieve_XM_info(line)
		if (not info):
			raise ValueError("Unknown information at line %d" % n)
		read, barcode, probe_r1, probe_r2 = info 
		if (read not in ("1","2") or 
			probe_r1 not in probe_info or
			probe_r2 not in probe_info):
			raise ValueError("Unknown information column \"%s\" at line %d" % (info, n))
		if (probe_r1 != probe_r2):
			#exclude off-target alignments
			continue
		name = line[0]
		if (read_name != name):
			if (alignments):
				#output the alignment group
				yield alignments
			#initialize a new group
			read_name = name
			alignments = []
		flag = SAMFlag(line[1]) #flag
		alignments.append((info, flag, line))
	if (alignments):
		#output the last group
		yield alignments
	pf.stdout.close()
	if (head is not None):
		raise ValueError("Bam head not loaded from <%s>" % cmd)

def retrieve_alignments(sample_name, fastq_r1, fastq_r2, probe, outpath, genome_refseq, mtdna_refseq, mtdna_cord_offset):
	"""
	align the r1 and r2 reads to the reference genome and mitochondrial genome
	
	Arguments
	----------
	sample_name: the name of the sample
	fastq_r1: r1 fastq file
	fastq_r2: r2 fastq file
	probe: probe informatioin obtained from retrieve_probe
	outpath: the path to store the output alignment files
	genome_refseq: the path to the reference genome
	mtdna_refseq: the path to the reference mtdna
	mtdna_cord_offset: the position shift in parsing the mtdna
	
	Returns
	----------
	a tuple of the paths to the alignment files (mtdna and ndna)
	
	"""
	#variables for counting reads output
	#genome_read_count = {}
	read_count = {}
	
	#create a dict for fast retrieval of probe info 
	probe_info = {}
	probe_ndna = {}
	#initialize
	for i in probe:
		name, chr = i[:2]
		probe_info[name] = i
		if (chr != "chrM"):
			probe_ndna[name] = 1
		#reads mapped to reference genome
		#genome_read_count[name] = {}
		#reads after QC, umapped reads, reads mapped to wrong chr, strand error, position error or low mapq (<20), numts
		read_count[name] = [0, 0, 0, 0, 0, 0]
	
	#keep track of mtdna fragments can be aligned to ndna
	numts = defaultdict(lambda: 0)
	numts_alignments = {}
	#temporaily store ndba fragments with proper insert length
	ndna = []
	#parse complete genome mapping results
	ndna_alignments = align_reads_to_reference(fastq_r1, fastq_r2, genome_refseq, probe_info)
	head = ndna_alignments.next() #head information
	for alignments in ndna_alignments:
		alignment_pair = [None,None]
		for info, flag, alignment in alignments:
			r, barcode, r1_probe, r2_probe = info
			#count = genome_read_count[r1_probe]
			if (not flag.primary):
				#exclude secondary alignments
				continue
			if (flag.unmapped or flag.mate_unmapped):
				if (r1_probe in probe_ndna):
					read_count[r1_probe][1] += 1
				#unmapped
				continue
			r = int(info[0]) #1/2
			if (alignment_pair[r-1]):
				for i in alignments:
					print >>sys.stderr, "#", i[2]
				raise ValueError("Ambiguous alignment %s (%s)" % (alignment[0], ":".join(info)))
			alignment_pair[r-1] = alignment #temporarily store reads until both are found
		#found both reads in a pair
		if (alignment_pair[0] and alignment_pair[1]):
			alignment = alignment_pair[0]
			if (probe_info[r1_probe][1] != "chrM"):
				#temperarily store proper nDNA alignments 
				ndna.append(alignments)
			else:
				numts_coord = []
				for alignment in alignment_pair:
					chr = alignment[2]
					mapq = int(alignment[4])
					if (chr != "chrM" and mapq > 10):
						#treat it as numts if it can be mapped to nDNA with mapping quanlity >10 (<10% mapping errors)
						numts[alignment[0]] += 1
						#store numts alignment information
						numts_coord.append("%s:%s:%d"%(chr, alignment[3],mapq))
						#print "numts", alignment[0], numts[alignment[0]]
						if (numts[alignment[0]] > 1):
							#remove the entire read family if a pair of it was mapped to nDNA
							numts[r1_probe + "-" +  barcode] = 2
							#R1 and R2 mapping stats: chr:pos:mapq_chr:pos:mapq
							numts_alignments[r1_probe + "-" +  barcode] = "_".join(numts_coord)
	
	#output bam file (pipe to samtools)
	if (not outpath):
		outpath = "."
	outfile = outpath + os.path.sep + sample_name + ".%s.bam"
	cmd = "%s view -Sb -o %s - " % (samtools, outfile) 
	
	#output data on nDNA and mtDNA separately
	out_ndna = pipe_input(cmd % "ndna")
	out_mtdna = pipe_input(cmd % "mtdna")
	
	#return output file names
	ret = (outfile % "mtdna", outfile % "ndna")
	
	#retain head informtion
	for line in head:
		out_ndna.stdin.write(line+"\n")
	
	#parse mtDNA mapping results
	mtdna_alignments = align_reads_to_reference(fastq_r1, fastq_r2, mtdna_refseq, probe_info)
	head = mtdna_alignments.next()
	#retain head information
	for line in head:
		out_mtdna.stdin.write(line+"\n")
	#iterate all alignments
	for alignments in chain(ndna, mtdna_alignments):
		alignment_pair = [None,None]
		for info, flag, alignment in alignments:
			r, barcode, r1_probe, r2_probe = info
			count = read_count[r1_probe]
			if (not flag.primary):
				#exclude secondary alignments
				continue
			if (flag.unmapped or flag.mate_unmapped):
				#unmapped
				if (r1_probe not in probe_ndna):
					count[1] += 1
				continue
			name, chr, start, end, insert, r1_strand, r2_strand = probe_info[info[2]][:7]
			if (alignment[2] != chr or alignment[6] != "="):
				#mapped to a wrong chromosome
				count[1] += 1
				continue
			r = int(info[0]) #1/2
			if (r == 1):
				#check strand consistency
				if (flag.strand_reverse != (r1_strand == "-")):
					count[3] += 1
					continue
				if (flag.mate_strand_reverse != (r2_strand == "-")):
					count[3] += 1
					continue
			else:
				if (flag.strand_reverse == (r1_strand == "-")):
					count[3] += 1
					continue
				if (flag.mate_strand_reverse == (r2_strand == "-")):
					count[3] += 1
					continue
			if ((r1_strand == "+") == (r == 1)):
				#check reference position
				align_start = int(alignment[3])
				if (chr == "chrM"):
					align_start -= mtdna_cord_offset
				align_end = align_start + int(alignment[8])
				#allow a certain degree of shift in the aligned position
				if (align_start < start - align_shift_max or align_start > start + align_shift_max):
					#print r1_strand, r, r1_probe, align_start, start, align_end, end
					count[4] += 2 #exclude this pair
					continue
				if (align_end < end - align_shift_max or align_end > end + align_shift_max):
					#print r1_strand, r, r1_probe, align_start, start, align_end, end
					count[4] += 2 #exclude this pair
					continue
			if (alignment_pair[r-1]):
				#for i in alignments:
				#	print i[2]
				raise ValueError("Ambiguous alignment %s (%s)" % (alignment[0], ":".join(info)))
			alignment_pair[r-1] = alignment
		if (alignment_pair[0] and alignment_pair[1]):
			#the minimum alignment score is 20 for both reads 
			if (int(alignment_pair[0][4]) >= 20 and int(alignment_pair[1][4]) >= 20):
				if (alignment_pair[0][2] == "chrM"):
					if (numts.get(alignment_pair[0][0], 0) > 1 or numts.get(r1_probe + "-" +  barcode, 0) > 1):
						#read pairs or other pairs with the same barcode can be mapped to nDNA
						#print r1_probe, numts.get(alignment_pair[0][0])
						count[5] += 2
						#get coordinates
						coord = numts_alignments.get(r1_probe + "-" +  barcode)
						if (coord):
							coord = "BWA:%s"%coord
						else:
							coord = "BWA"
						#mark potential numts in reads annotation
						alignment_pair[0].append("XQ:Z:NUMTS(%s)"%coord)
						alignment_pair[1].append("XQ:Z:NUMTS(%s)"%coord)
					else:
						count[0] += 2
					#output if both reads pass all of the filters
					out_mtdna.stdin.write("\t".join(alignment_pair[0]) + "\n")
					out_mtdna.stdin.write("\t".join(alignment_pair[1]) + "\n")
				else:
					count[0] += 2
					#output if both reads pass all of the filters
					out_ndna.stdin.write("\t".join(alignment_pair[0]) + "\n")
					out_ndna.stdin.write("\t".join(alignment_pair[1]) + "\n")
					
			else:
				#only one read found
				count[4] += 2
	
	out_ndna.stdin.close()
	out_mtdna.stdin.close()
	for i in sorted(read_count.keys()):
		print >>LOG, "alignment", sample_name, i, " ".join(map(str, read_count[i]))
	
	#ret = (outfile % "mtdna", outfile % "ndna")
	return ret

def process_alignments(alignment_file, refseq, outfile = None, recalibrate_base_qual = True):
	"""
	a wrapper function to call local realignment and base quality recalibration
	
	Arguments
	----------
	alignment_file: the path to the original alignment file
	refseq: the path to the reference genome
	outfile:  the name of the output file
	recalibrate_base_qual: true if recalibratin gbase qualily
		
	Returns
	----------
	the path to the processed alignment file
	
	"""
	
	global debug
	#debug = True
	
	if (alignment_file.endswith(".bam")):
		alignment_file = alignment_file[:-4]
	if (not outfile):
		outfile = alignment_file
	elif (outfile.endswith(".bam")):
		outfile = outfile[:-4]
	infile = alignment_file + ".bam"
	if (not os.path.exists(infile)):
		raise IOError("Cannot open %s" % infile)
	execute("%s sort -o %s.sorted.bam %s" % (samtools, outfile, infile), run = not debug)
	#execute("%s index %s.bam" % (samtools, outfile))
	infile = outfile + ".sorted"
	if (recalibrate_base_qual):
		outfile = outfile + ".sorted.realign.recal"
		#delete the existing file
		execute("rm -f %s.bam" % outfile, run = not debug)
		#multiple alignment to align indels and recalculate md
		execute("%s -c -f %s < %s.bam | %s calmd -EArb - %s >%s.bam" % (bamleftalign, refseq, infile, samtools, refseq, outfile), run = not debug)
	else:
		outfile = outfile + ".sorted.realign"
		#delete the existing file
		execute("rm -f %s.bam" % outfile, run = not debug)
		#multiple alignment to align indels and recalculate md
		execute("%s -c -f %s < %s.bam >%s.bam" % (bamleftalign, refseq, infile, outfile), run = not debug)
	#create bam index
	execute("%s index %s.bam" % (samtools,outfile), run = not debug)
	#delete the temperary sorted bam file
	execute("rm -f %s.bam" % infile, run = not debug)
	return outfile

def retrieve_read_family(bam_file, prb_name, chr, start, end, barcode_index):
	"""
	group reads according to their attached barcode
	
	Arguments
	----------
	bam_file: the path to the alignment file
	prb_name: probe to be processed
	chr: chromosome
	start: start position; used to extract reads
	end: end position; used to extract reads
	barcode_index: a dict of barcode (grouped based on possible Ns) 
		
	Returns
	----------
	a dict of read families (barcode as the key) and the associated paired-end reads
	
	"""
	#use samtools to extract reads
	cmd = "%s view %s.bam %s:%d-%d" % (samtools, bam_file, chr, start, end)
	pf = pipe_output(cmd)
	family = defaultdict(lambda: defaultdict(lambda: [None, None]))
	for line in pf.stdout:
		line = line.rstrip("\r\n").split("\t")
		read, barcode, r1_probe, r2_probe = retrieve_XM_info(line)
		if (r1_probe == prb_name):
			#mask barcode if there are similar barcodes with Ns (low-quality base treated as a wildcard)
			barcode = barcode_index.get(barcode, barcode)
			qname = line[0] #qname as the key
			alignment = family[barcode][qname]
			r = int(read)-1 #0:r1; 1:r2
			if (alignment[r]):
				raise ValueError("Ambiguous alignment %s (%s)" % (alignment[0], ":".join([read, barcode, r1_probe, r2_probe])))
			alignment[r] = line
	pf.stdout.close()
	return family

def estimate_family_num(prb_name, sample_name, family_count):
	"""
	estimate the average number of captured fragments being sequenced
	
	not the average number of reads in each read family
	
	Arguments
	----------
	prb_name: probe to be processed
	sample_name: name of the sample
	family_count: a dict of the number of read families with a specific family size (the number of reads in a read family) 
		
	Returns
	----------
	a tuple of the number of read family and the estimated number of fragments
	
	"""
	def poission_count(l):
		#compute a log-likelihood of the poisson distribution with the parameter lambda equal to l
		lnl = log(l)
		lnc = log(1.0-exp(-l))
		ret = 0.0 #-1000*exp(-l)*l
		for i, j in family_distr:
			ret += j*(i*lnl-l-log(factorial(i))-lnc)
		return -ret
	family_distr = [] 
	nfamily = float(sum(family_count.values()))
	nsample = 0
	for i, j in sorted(family_count.items()):
		print >>LOG, "family_size",sample_name, prb_name, i, j
		nsample += i*j
		family_distr.append([float(i),int(j/nfamily*1000)])
	ret = optimize.minimize_scalar(poission_count, bounds=(0.0, 20.0), method='bounded')
	print >>LOG, "family_num", sample_name, prb_name, nfamily, nsample, ret.x, nsample/ret.x, 3*nsample/ret.x
	return(nfamily, nsample/ret.x)

def call_consensus_reads(alignments):
	"""
	compute the base of a sequence in a read family using the Bayesian rule
	
	Arguments
	----------
	alignments: a list of base and base quality of reads in a read family at a certain aligned position
		
	Returns
	----------
	a tuple of the most likely base and its estimated base quality
	
	"""
	#combining qulaity score using the bayesian rule
	#used 1/3 the error rate
	merge_corrected = True
	#for insertions and deletions
	deletion  = 0
	insertion = defaultdict(lambda: [])
	#the minimum rate of reads to call a indel
	indel_frac = 0.7
	
	#the minimum rate of reads to call a single base
	base_frac = 0.3
	
	#likelihood of A,T,C,G 
	r = [1.0,1.0,1.0,1.0]
	n = 0
	for s, q in alignments:
		if (q is not None and q < 4):
			continue
		if (s == "A"):
			#phred to probability
			n += 1
			q = prob(q); r[0] *= 1-q; #true
			if (merge_corrected): q = q/3.0
			r[1] *= q; r[2] *= q; r[3] *= q #errors
		elif (s == "T"):
			n += 1
			q = prob(q); r[1] *= 1-q; #true
			if (merge_corrected): q = q/3.0
			r[0] *= q; r[2] *= q; r[3] *= q #errors
		elif (s == "C"):
			n += 1
			q = prob(q); r[2] *= 1-q; #true
			if (merge_corrected): q = q/3.0
			r[0] *= q; r[1] *= q; r[3] *= q #errors
		elif (s == "G"):
			n += 1
			q = prob(q); r[3] *= 1-q; #true
			if (merge_corrected): q = q/3.0
			r[0] *= q; r[1] *= q; r[2] *= q #errors
		elif (s == "."):
			deletion += 1
		elif (len(s) > 1):
			insertion[s].append(float(sum(q))/len(q)) #use average quality score
		elif (s == "N"):
			pass
		else:
			raise ValueError("Unknown nucleotide \"%s\"" % s)
	#no effective reads
	if (n == 0 and deletion == 0 and not insertion):
		return "N", 0
	#deletion
	if (deletion > len(alignments)*indel_frac):
		return ".", 0
	#insertion
	if (insertion):
		major_insert = sorted([(len(j),i, j) for i, j in insertion.items()], reverse = True)[0]
		insertion = sum([len(i) for i in insertion.values()])
		if (major_insert[0] > len(alignments)*indel_frac):
			#quality average over all supports
			return major_insert[1], [int(round(sum(major_insert[2])/major_insert[0])),]*len(major_insert[1])
	rtotal = sum(r)
	#N if more than base_frac of bases were removed due to low-quality or indels
	if (rtotal == 0 or n < len(alignments)*base_frac):
		return "N", 0
	rml = sorted([(j,i) for i, j in enumerate(r)], reverse = True)[0]
	#take the first one if two nts have the same quality. Its quality will be in the range of 0.25~0.5, and will be filtered out in following steps.
	base = ["A","T","C","G"][rml[1]]
	if (merge_corrected):
		seq_error = (1-rml[0]/rtotal)
	else:
		seq_error = (1-rml[0]/rtotal)/3.0
	if (seq_error <= 0):
		#too small; use the largest value
		return base, 223
	else:
		#rounded to the nearest integer
		return base, int(round(-10*log10(seq_error)))

def build_consensus_reads(family, strand, clip_end, ref_seq, ref_seq_offset, min_base_qual):
	"""
	estimate the consensus read sequence of the paired-end reads in a read family using the Bayesian rule
	
	Arguments
	----------
	family: a dict of all paired-end reads in a read family
	strand: stand of the mtdna reads used to determine start and end positions
	clip_end: clip the overlapping region if turned on
	ref_seq: the reference mtdna sequence
	ref_seq_offset: the position offset used in parsing mtdna variants
	min_base_qual: minimum base quality score for determining high-quality bases
		
	Returns
	----------
	a list containing
	1. alignment line (sam format) of the consensus read
	2. number of mismatch bases
	3. number of paired-end reads in the read family
	4. number of low-quality bases
	5. sum of the number of low-quality bases and the length of gap (Ns)
	6. existing annotation of the paired-end reads
	
	"""
	fam_size = 0
	#fetch coordinates of consensus sequence
	#reads on the left-heand side
	rl_cord = [sys.maxint, -1]
	rl_alignments = []
	#reads on the right-hand side
	rr_cord = [sys.maxint, -1]
	rr_alignments = []
	
	#directly copy other information from the alignments provided
	#set the name of the consensus read as that of the first paired-end read
	chr = family[0][0][2]
	name = family[0][0][0]
	
	#use the maximum mapq of the family for consensus read
	mapq = 0
	
	#treat the consensus read as a single-end read
	flag = 0x42
	
	#existing xq information
	annot = {}
	
	#rl ===> <=== rr
	for alignment in family:  
		r1, r2 = alignment
		
		#retain existing quality annotation
		r1_annot = retrieve_tag_info(r1, "XQ", "Z")
		if (r1_annot):
			for i in r1_annot.split(","):
				if (i):
					annot[i] = 1
		r2_annot = retrieve_tag_info(r2, "XQ", "Z")
		if (r2_annot):
			for i in r2_annot.split(","):
				if (i):
					annot[i] = 1
		
		#treat as single-ended reads
		if (strand == "+"):
			#r1 ===> <=== r2
			rl_alignment = SAMCigar(r1[5], int(r1[3]), r1[9], r1[10])
			rr_alignment = SAMCigar(r2[5], int(r2[3]), r2[9], r2[10])
		else:
			#r2 ===> <=== r1
			rl_alignment = SAMCigar(r2[5], int(r2[3]), r2[9], r2[10])
			rr_alignment = SAMCigar(r1[5], int(r1[3]), r1[9], r1[10])
		#check soft clip length
		soft_clip = 0
		for i, j in chain(rl_alignment.cigar, rr_alignment.cigar):
			if (j == "S"):
				soft_clip += i
		if (soft_clip >= soft_clip_ratio*(len(r1[9])+len(r2[9]))):
			continue
		
		fam_size += 1
		
		gap_size = rr_alignment.ref_start - rl_alignment.ref_end
		#left and right reads overlap 
		if (gap_size < 0 and clip_end):
			#clip the read end with the lower quality
			qual_left = 0
			qual_left_num = 0 
			for i in range(rl_alignment.ref_len+gap_size, rl_alignment.ref_len):
				q = rl_alignment.align_qual[i]
				if (isinstance(q, list)):
					qual_left += sum(q)
					qual_left_num += len(q)
				elif (q >= 0):
					qual_left += q
					qual_left_num += 1
			qual_right = 0
			qual_right_num = 0 
			for i in range(-gap_size):
				q = rr_alignment.align_qual[i]
				if (isinstance(q, list)):
					qual_right += sum(q)
					qual_right_num += len(q)
				elif (q > 0):
					qual_right += q
					qual_right_num += 1
			avg_qual_left = float(qual_left)/qual_left_num
			avg_qual_right = float(qual_right)/qual_right_num
			if (avg_qual_left > avg_qual_right):
				#clip head of read right
				rr_alignment.clip_head(-gap_size)
			elif (avg_qual_right > avg_qual_left):
				#clip tail of read left
				rl_alignment.clip_tail(-gap_size)
			else:
				#if equal, clip read 2
				if (strand == "+"):
					rr_alignment.clip_head(-gap_size)
				else:
					rl_alignment.clip_tail(-gap_size)
			#raise NotImplementedError("Overlapping reads are not supported now;\n Use other tools to clip before building consensus.")
		rl_cord[0] = min(rl_cord[0], rl_alignment.ref_start)
		rl_cord[1] = max(rl_cord[1], rl_alignment.ref_end)
		rr_cord[0] = min(rr_cord[0], rr_alignment.ref_start)
		rr_cord[1] = max(rr_cord[1], rr_alignment.ref_end)
		rl_alignments.append(rl_alignment)
		rr_alignments.append(rr_alignment)
		#take the highest mapping score
		mapq = max(mapq, int(r1[4]))
	
	#no reads left
	if (not fam_size):
		return [None, {}, 0, 0, 0, {}]
	
	if (rl_cord[1] < rr_cord[0]):
		#no overlapping region: rl/gap/rr 
		segments = ((rl_cord[0], rl_cord[1], rl_alignments), (rl_cord[1], rr_cord[0], []), (rr_cord[0], rr_cord[1], rr_alignments))
	else:
		#overlapping region: rl/rl+rr/rr
		segments = ((rl_cord[0], rr_cord[0], rl_alignments), (rr_cord[0], rl_cord[1], rl_alignments + rr_alignments), (rl_cord[1], rr_cord[1], rr_alignments))

	consensus_cigar=[]
	consensus_seq = []
	consensus_qual = []
	cur_stat = ""
	mismatch = {} #keep track of mismatches
	lowqual = 0 #number of low-quality bases
	ref_len = 0 #length of the aligned reference sequence
	gap_size = 0 #length of the gap (Ns)
	for start, end, alignments in segments:
		if (alignments):
			ref_len += (end - start)
			for pos in xrange(start, end):
				s, q = call_consensus_reads([i.fetch(pos) for i in alignments])
				ref_pos = pos-ref_seq_offset
				if (ref_pos >= len(ref_seq)):
					ref = "N"
				else:
					ref = ref_seq[ref_pos]
				if (len(s) > 1):
					#insertion
					if (pos != end-1):
						consensus_cigar.append([len(s)-1, "I"])
						consensus_cigar.append([1, "M"])
						cur_stat = "M"
					else:
						consensus_cigar.append([len(s), "I"])
						cur_stat = "I"
					consensus_seq.extend(s)
					consensus_qual.extend(q)
					mismatch[pos] = s
					if (float(sum(q))/len(q) < min_base_qual):
						lowqual += len(s)
				elif (s == "."):
					#deletion
					if (cur_stat != "D"):
						consensus_cigar.append([1,"D"])
					else:
						consensus_cigar[-1][0] += 1
					#mismatch[pos] = ""
					cur_stat = "D"
				else:
					#match or mismatch
					if (cur_stat != "M"):
						consensus_cigar.append([1,"M"])
					else:
						consensus_cigar[-1][0] += 1
					consensus_seq.append(s)
					consensus_qual.append(q)
					if (q < min_base_qual):
						lowqual += 1
					if (s != "N" and s != ref):
						mismatch[pos] = s
					cur_stat = "M"
		else:
			#fill gap with Ns
			gap_size = end - start
			if (cur_stat == "M"):
				consensus_cigar[-1][0] += gap_size
			else:
				consensus_cigar.append([gap_size, "M"])
				cur_stat = "M"
			consensus_seq.extend(["N",]*gap_size)
			consensus_qual.extend([0,]*gap_size)
	
	consensus_cigar = "".join([str(i[0])+i[1] for i in consensus_cigar])
	
	consensus_seq = "".join(consensus_seq)
	#printable characters have the largest ASCII code of 126 (equivalent to a phred-like error rate < 1e-9)
	consensus_qual = "".join([char(33+min(93,i)) for i in consensus_qual])
	
	return [[name, str(flag), chr, str(rl_cord[0]), str(mapq), consensus_cigar, "*", "0", "0", consensus_seq, consensus_qual], mismatch, fam_size, lowqual, lowqual+gap_size, annot]

def retrieve_re_alignments(sample_name, alignments_mtdna, alignments_ndna, outpath, genome_refseq, mtdna_refseq, recalibrate_base_qual):
	"""
	realign and recalibrate reads
	
	Arguments
	----------
	sample_name: name of the sample
	alignments_mtdna: path to alignment file of mtdna reads returned from retrieve_alignment(...)
	alignments_ndna: path to alignment file of ndna reads returned from retrieve_alignment(...)
	outpath: path to store the resulting consensus alignment files
	genome_refseq: path to the reference genome
	mtdna_refseq: path to the reference mtdna
	recalibrate_base_qual: recalibrate base quality if turned on (default: yes)
		
	Returns
	----------
	a tuple of paths to the re-alignment files
	
	"""
	
	#mtdna
	alignments_mtdna = process_alignments(alignments_mtdna, mtdna_refseq, outpath + os.path.sep + sample_name + ".mtdna", recalibrate_base_qual)
	#ndna
	alignments_ndna = process_alignments(alignments_ndna, genome_refseq, outpath + os.path.sep + sample_name + ".ndna", recalibrate_base_qual)
	
	return (alignments_mtdna + ".bam", alignments_ndna + ".bam") 
	
def retrieve_consensus_reads(sample_name, alignments_mtdna, alignments_ndna, probe, barcode_file, outpath, genome_refseq, mtdna_refseq, mtdna_cord_offset, mtdna_alt_seq, recalibrate_base_qual, clip_end, consensus_nm, consensus_xf, consensus_gap, consensus_qual, consensus_lqbase):
	"""
	call consensus reads from paired-end reads
	
	Arguments
	----------
	sample_name: name of the sample
	alignments_mtdna: path to alignment file of mtdna reads returned from retrieve_re_alignment(...)
	alignments_ndna: path to alignment file of ndna reads returned from retrieve_re_alignment(...)
	probe: probe information returned from retrieve_probe(...)
	barcode_file: path to the barcode file output in retrieve_tags(...)
	outpath: path to store the resulting consensus alignment files
	genome_refseq: path to the reference genome
	mtdna_refseq: path to the reference mtdna
	mtdna_cord_offset: position offset used in parsing mtdna variants
	mtdna_alt_seq: alternate mtDNA sequences. e.g., sequences of numts, alignment artifacts
	recalibrate_base_qual: recalibrate base quality if turned on (default: yes)
	clip_end: clip the overlapping region of r1 and r2 (default: no)
	consensus_nm: the maximum rate of mismatches in consensus reads (default: 0.025 >=8 in most target regions)
	consensus_xf: the minimum family size of consensus reads
	consensus_gap: the maximum length of gap allowed in each consensus read (output as N; default: 50)
	consensus_qual: the minimum base quality score used to count low quality bases (default: 0; not used)
	consensus_lqbase: the maximum number of low quality bases allowed in each consensus read (default:200; not used)
	
	Returns
	----------
	a tuple of paths to the alignment files of consensus reads
	
	"""
	if (not outpath):
		outpath = "."
	#retrieve barcode information processed in retrieve_tags
	barcode_cache = defaultdict(lambda:list())
	#default barcode file: *.barcode
	if (os.path.exists(barcode_file)):
		#read barcode group information if it exists
		#barcode groups were defined during tag retrieval where barcodes diff nt at a low-quality site were considered from the same read family
		#retrieve barcode information processed in retrieve_tags
		barcode_dmin = 3
		with open(barcode_file) as fh:
			for line in fh:
				line = line.rstrip("\r\n").split("\t")
				if (len(line) != 5):
					continue
				r1_probe, r2_probe, barcode, other_barcode, num = line
				if (r1_probe == r2_probe):
					#only consider on-target fragments
					N = barcode.count("N")
					other_barcode = [i for i in other_barcode.strip().split(",") if i and i.count("N") < barcode_dmin]
					if (other_barcode and N < barcode_dmin):
						barcode_cache[r1_probe].append([barcode.count("N"), -int(num), barcode, other_barcode])
						
			#sort barcode based on base uncertainty
			for r1_probe in barcode_cache:
				barcode_cache[r1_probe].sort()
	
	#build universal barcode families
	barcode_index = defaultdict(lambda:dict())
	for r1_probe in barcode_cache:
		for N, n, barcode, other_barcode in barcode_cache[r1_probe]:
			if (barcode not in barcode_index[r1_probe]):
				#use barcode to represent all other barcodes
				for i in other_barcode:
					if (i not in barcode_index[r1_probe]):
						barcode_index[r1_probe][i] = barcode
					#else:
					#	print "Ambigious barcode", r1_probe, N, n, barcode, other_barcode, barcode_index[r1_probe][i]
		#print r1_probe, len(barcode_index[r1_probe])
	del barcode_cache
	#print len(barcode_index)
	
	#alignments_mtdna = process_alignments(alignments_mtdna, mtdna_refseq, outpath + os.path.sep + sample_name + ".mtdna", recalibrate_base_qual)
	#alignments_ndna = process_alignments(alignments_ndna, genome_refseq, outpath + os.path.sep + sample_name + ".ndna", recalibrate_base_qual)
	
	#create index file for fast retrieval
	if (alignments_mtdna.endswith(".bam")):
		alignments_mtdna = alignments_mtdna[:-4]
	if (not os.path.exists(alignments_mtdna + ".bam.bai")):
		execute("%s index %s.bam" % (samtools,alignments_mtdna), run = not debug)
	if (alignments_ndna.endswith(".bam")):
		alignments_ndna = alignments_ndna[:-4]
	if (not os.path.exists(alignments_ndna + ".bam.bai")):
		execute("%s index %s.bam" % (samtools,alignments_ndna), run = not debug)
	
	#retain head information
	cmd = samtools + " view -H %s.bam"
	head_mtdna = pipe_output(cmd % alignments_mtdna)
	head_ndna = pipe_output(cmd % alignments_ndna)
	
	#output consensus reads in bam format
	if (not outpath):
		outpath = "."
	outfile = outpath + os.path.sep + sample_name + ".%s.consensus.bam"
	cmd = "%s view -hbS -o %s - " % (samtools, outfile)
	
	#output data on nDNA and mtDNA separately
	out_mtdna_consensus = pipe_input(cmd % "mtdna")
	out_mtdna_consensus.stdin.write(head_mtdna.stdout.read())
	head_mtdna.stdout.close()
	
	out_ndna_consensus = pipe_input(cmd % "ndna")
	out_ndna_consensus.stdin.write(head_ndna.stdout.read())
	head_ndna.stdout.close()
	
	#return consensus alignment file names
	ret = (outfile % "mtdna", outfile % "ndna")
	
	#keep track of family size distribution
	mtdna_family_count = defaultdict(lambda: 0)
	ndna_family_count = defaultdict(lambda: 0)
	mtdna_nc = 0
	ndna_nc = 0
	
	#reference sequences are used for counting mismatches
	seq_ndna = FastaIO(genome_refseq)
	seq_mtdna = FastaIO(mtdna_refseq)
	
	#process one probe at a time
	for i in probe:
		prb_name, chr, start, end, amp_size, r1_strand = i[:6]
		#if (prb_name not in ["B10"]):
		#	continue
		if (chr == "chrM"):
			mtdna_nc += 1
			read_family = retrieve_read_family(alignments_mtdna, prb_name, chr, start+mtdna_cord_offset, end+mtdna_cord_offset, barcode_index[prb_name])
			if (not read_family):
				estimate_family_num(prb_name, sample_name, {})
				continue
			#relax the mismatch filter for sites in the D-loop region due to a higher mutation rate
			if (start >= 16000 or end < 560):
				nm_factor = 1.2 #default: 2.5% ==> 3%
			else:
				nm_factor = 1.0
			start += mtdna_cord_offset
			end += mtdna_cord_offset
			#retrieve reference sequence
			chr, offset, ref_seq = seq_mtdna.fetch(chr, start-align_shift_max, end+align_shift_max)
			#retrieve alternate sequences, which are considered as numts or alignment artifacts
			alt_seq = retrieve_alternate_sequences_by_region(mtdna_alt_seq, start, end)
		else:
			ndna_nc += 1
			read_family = retrieve_read_family(alignments_ndna, prb_name, chr, start, end, barcode_index[prb_name])
			if (not read_family):
				estimate_family_num(prb_name, sample_name, {})
				continue
			nm_factor = 1.0
			#retrieve reference sequence
			chr, offset, ref_seq = seq_ndna.fetch(chr, start-align_shift_max, end+align_shift_max)
			#no alternate sequences
			alt_seq = {}
		#print ref_seq, start, end, offset, 
		#find the most common cigar string for R1 and R2
		r1_cigar = defaultdict(lambda: 0)
		r2_cigar = defaultdict(lambda: 0)
		#record the family size distribution
		family_count = defaultdict(lambda: 0)
		#keep track of alternative major alleles
		consensus_allele = defaultdict(lambda:defaultdict(lambda:0))
		#total reads before quality control
		total_reads = 0
		#number of reads after quality control
		usable_reads = 0
		#temparily store consensus sequences for building major allelles
		consensus_info = []
		for barcode, family in read_family.items():
			for r1, r2 in family.values():
				r1_cigar[r1[5]] += 1
				r2_cigar[r2[5]] += 1
			#retrieve consensus read information
			consensus = build_consensus_reads(family.values(), r1_strand, clip_end, ref_seq, offset+1, consensus_qual)
			#retain barcode information
			consensus.append(barcode)
			if (consensus[0] is not None):
				consensus_info.append(consensus)
				#record the position and allele of each mismatch 
				for pos, allele in consensus[1].items():
					consensus_allele[pos][allele] += 1
		#print out the most abundant cigar string
		#R1
		for j,i in sorted([[j,i] for i,j in r1_cigar.items()], reverse = True):
			print prb_name, "r1", i, j, float(j)/sum(r1_cigar.values())
			break
		#R2
		for j,i in sorted([[j,i] for i,j in r2_cigar.items()], reverse = True):
			print prb_name, "r2", i, j, float(j)/sum(r2_cigar.values())
			break
		#determine the major allele for each mismatch found
		major_allele = {}
		for pos in consensus_allele:
			#print (major/minor) allele at each consensus read site to stdout (used for debugging, to be removed)
			print prb_name, pos, dict(consensus_allele[pos]) 
			for allele, count in consensus_allele[pos].items():
				#must account for over half of the consensus sequences found and occur at least major_allele_count times 
				if (count > major_allele_ratio*len(consensus_info) and count > major_allele_count):
					major_allele[pos] = allele
					print prb_name, "major", pos, allele
				else:
					print prb_name, "minor", pos, allele, count, major_allele_ratio, major_allele_count, len(consensus_info) 
		#keep track of #mismatches observed
		nm_max = consensus_nm*nm_factor*(end-start)
		nm_window_max = 2*consensus_nm*nm_factor*100
		for consensus, mismatch, fam_size, lowqual, ns, annot, barcode in consensus_info:
			#quality filters on family size, reads quality and gap size
			total_reads += fam_size
			if (fam_size >= consensus_xf and lowqual <= consensus_lqbase and ns <= consensus_gap):
				#count total #mismatches and #mismatches in 100 bp sliding windows
				nm = []
				nm_window = [0,]*((end-start)/50+1)
				for pos, allele in mismatch.items():
					if (major_allele.get(pos) != allele):
						##check for remaining NuMTs @ some known sites
						nm.append([pos, allele])
						window_id = (pos-start)/50
						if (window_id < len(nm_window)):
							nm_window[window_id] += 1
							if (window_id > 0):
								nm_window[window_id-1] += 1
				#same as the reference sequence but different from the consensus (i.e. at homoplasmic sites)
				#thus, nm may be larger than len(mismatch)
				for pos, allele in major_allele.items():
					if (mismatch.get(pos) is None):
						nm.append([pos,allele])
						window_id = (pos-start)/50
						if (window_id < len(nm_window)):
							nm_window[window_id] += 1
							if (window_id > 0):
								nm_window[window_id-1] += 1
				#if (nm != len(mismatch)):
				#	print len(mismatch), nm, nm_window, mismatch
				#remove reads with clutered mismatches in sliding windows of a breath of 100bp and overlaping regions of 50bp
				if (len(nm) > nm_max or sorted(nm_window, reverse=True)[0] > nm_window_max):
					annot["EXMISMATCH"] = 1
				#calculate the distances between read sequences and known nuclear mitochondrial DNA sequences
				#compare the distances between read  
				if (alt_seq):
					nm_alt, id = retrieve_alternate_sequences_by_similarity(alt_seq, mismatch)
					if (nm_alt < len(nm)):
						#print nm_alt, id, alt_seq, nm 
						annot["NUMTS(%s)"%id] = 1
				if (not annot):
					#no detectable quality issues
					annot["PASS"] = 1
					family_count[fam_size] += 1
					usable_reads += fam_size
				#output alternate alleles and corresponding positions
				alt_allele = [str(p)+a for p, a in sorted(mismatch.items())]
				#output minor alleles and corresponding positions
				minor_allele = [str(p)+a for p, a in sorted(nm)]
				
				#info = ["NM:i:%d"%len(mismatch), "XA:Z:%s"%",".join(map(str,[nm,]+sorted(nm_pos))), "XN:Z:%s"%",".join(map(str,[len(numts),]+numts)), "XQ:Z:%s"%",".join(annot.keys()), XF:i:%d" % fam_size, "XM:Z:%s:%s:%s"%(barcode,prb_name,prb_name)]
				info = ["NM:i:%d"%len(mismatch), "XP:Z:%s"%",".join(map(str,[len(alt_allele),]+alt_allele)), "XA:Z:%s"%",".join(map(str,[len(minor_allele),]+minor_allele)), "XQ:Z:%s"%",".join(annot.keys()), "XF:i:%d" % fam_size, "XM:Z:%s:%s:%s"%(barcode,prb_name,prb_name)]
				#output the consensus sequences for each read family as a single-end read
				#XP: number of mismatches compared to the reference sequence: mismatched position and allele
				#XA: number of mismatches compared to the major mtdna allele: mismatched position and allele
				if (chr == "chrM"):
					out_mtdna_consensus.stdin.write("\t".join(consensus + info) + "\n")
				else:
					out_ndna_consensus.stdin.write("\t".join(consensus + info) + "\n")
		
		#estimate average family size
		estimate_family_num(prb_name, sample_name, family_count)
		#pooling faimly size counts for mtDNA and nDNA probes, respectively
		if (chr == "chrM"):
			#sum up all read families mapped to mtdna
			for i, j in family_count.items():
				mtdna_family_count[i] += j	
		else:
			#sum up all read families mapped to ndna
			for i, j in family_count.items():
				ndna_family_count[i] += j
		
		#print out reads before and after quality filtering
		print prb_name, "S", total_reads, usable_reads, float(usable_reads)/total_reads if (total_reads > 0) else "-"
	#estimate average family size for all possible mtDNA fragments
	#fitting a poisson distribution 
	mtdna_cn, mtdna_pool = estimate_family_num("mtDNA", sample_name, mtdna_family_count)
	#estimate average family size for all possible nDNA fragments
	ndna_cn, ndna_pool = estimate_family_num("nDNA", sample_name, ndna_family_count)
	#mtDNA copy number as the ratio between mtDNA coverage/family number and nDNA coverage/family number
	#this mtDNA copy number method does not take capture efficiency differences into consideration. (obsolete)
	if (ndna_cn > 0 and ndna_pool > 0):
		#this is a beta function for mtDNA copy number calculation, not used for real analyses at this point
		print >>LOG, "mito_copy_number", sample_name, mtdna_cn*ndna_nc/mtdna_nc/ndna_cn*2, mtdna_pool*ndna_nc/mtdna_nc/ndna_pool*2
	#close file handles
	out_ndna_consensus.stdin.close()
	out_mtdna_consensus.stdin.close()
	
	return ret

#example code: python align.py  -r1 fastq/mt_R1.fastq.gz -p probe_set.txt -o fastq --genome ~/data/ref_seq/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa --mtdna ~/data/ref_seq/chrMT/rCRS_chrM_16449-1-16569.fa --mtdna-offset 120  mt

def run(prog, args):
	"""
	run the main program
	
	Arguments
	----------
	prog: prog name
	args: arguments list
		
	Returns
	----------
	None
	
	"""
	#parse arguments
	parser = ArgumentParser(prog = prog, usage="%(prog)s [-options] sample\n", version="%(prog)s v0.1.1")
	#parser.add_argument("-r1", "--read1", type = str, required = True, nargs='+', dest = "r1", help="fastq file(s) for read 1")
	parser.add_argument("-r1", "--read1", type = str, nargs='+', dest = "r1", help="fastq file(s) for read 1")
	parser.add_argument("-r2", "--read2", type = str, nargs='+', dest = "r2", help="fastq file(s) for read 2")
	parser.add_argument("-p", "--probe", type = str, dest = "probe", help="the file with the probe information")
	parser.add_argument("-o", "--outpath", type = str, default = ".", dest = "outpath", help="path where to store the temporary alignment files")
	parser.add_argument("--consensus-outpath", type = str, dest = "consensus_outpath", help="path where to store the consensus alignment files (default: outpath)")
	parser.add_argument("--override", action = "store_true" , default = False, dest = "override", help="override output files (default: skip exiting output files)")
	parser.add_argument("--genome", type=str, dest="genome", help="the complete genome reference")
	parser.add_argument("--mtdna", type=str, dest="mtdna", help="the mtDNA reference sequence")
	parser.add_argument("--mtdna-offset", type=int, dest="mtdna_offset", help="the position offset used in parsing mtDNA read alignments")
	parser.add_argument("--numts", type=str, dest="numts", help="known nuclear mitochondrial DNA segments (HSD format)")
	parser.add_argument("--qual", type=int, default = 13, dest="min_qual", help="the minimum base quality for barcode sequence")
	parser.add_argument("--trim", type=int, default = 0, dest="trim", help="trim the low-quality ends in R1 and R2 with the specified base quality")
	parser.add_argument("-pr", "--probe-reverse", action = "store_true", default = False, dest = "probe_reverse", help = "reverse R1 and R2")
	parser.add_argument("--output-off-target", action = "store_true", default = False, dest = "output_off_target", help="output off-target R1 and R2 pairs")
	parser.add_argument("--no-recalibration", action = "store_true", default = False, dest = "no_recalibration", help="turn off base quality re-calibration")
	parser.add_argument("--clip-end", action = "store_true", default = False, dest = "clip_end", help="clip the overlapping region of R1 and R2 instead of merging it") #by default, two reads will be treated as a single-ended read
	parser.add_argument("--consensus-nm", type=float, default = 0.025, dest = "consensus_nm", help="the maximum rate of mismatches in consensus reads")
	parser.add_argument("--consensus-xf", type=int, default = 1, dest = "consensus_xf", help="the minimum family size of consensus reads")
	parser.add_argument("--consensus-gap", type = int, default = 50, dest = "consensus_gap", help = "the maximium length of gap allowed in each consensus read (output as N)")
	parser.add_argument("--consensus-qual", type = int, default = 0, dest = "consensus_qual", help = "the minimum base quality score used to count low quality bases")
	parser.add_argument("--consensus-lqbase", type = int, default = 200, dest = "consensus_lqbase", help = "the maximum number of low quality bases allowed in each consensus read")
	parser.add_argument("sample", type = str, help="sample name")
	options = parser.parse_args(args)
	
	#set default values for global variables
	logging_file = ""
	global LOG
	
	#get defaults
	#package path
	pack_path = os.path.dirname(os.path.abspath(__file__))
	arg_file = pack_path + os.path.sep + "stamp.arg"
	if (os.path.exists(arg_file)):
		default_arguments = read_arguments(arg_file)
	else:
		default_arguments = {} 
	
	global samtools
	samtools = default_arguments.get("samtools", samtools)
	
	global bwa 
	bwa = default_arguments.get("bwa", bwa)
	
	global bamleftalign
	bamleftalign = default_arguments.get("bamleftalign", bamleftalign)
	
	#set default values for options if not provided
	if (options.genome is None):
		options.genome = default_arguments.get("genome")
	
	if (options.mtdna is None):
		options.mtdna = default_arguments.get("mtdna")
		
	if (options.mtdna_offset is None):
		if (default_arguments.get("mtdna_offset")):
			options.mtdna_offset = int(default_arguments.get("mtdna_offset"))
	
	if (options.probe is None):
		options.probe = default_arguments.get("probe")
	
	#read probe information
	probe = retrieve_probes(options.probe)
	
	fastq_r1 = options.r1
	fastq_r2 = None
	if (fastq_r1 is not None):
		if (options.r2 is None):
			#use the same prefix of the R1 file for the R2 file
			fastq_r2 = [r1.replace("_R1","_R2") for r1 in fastq_r1]
		else:
			fastq_r2 = options.r2
		assert len(fastq_r1) == len(fastq_r2), "R1 and R2 do not match in numbers"
	
	if (options.outpath):
		outpath = options.outpath
	else:
		outpath = "."
	
	#output files
	r1 = outpath + os.path.sep + options.sample + "_R1.fastq.gz"
	r2 = outpath + os.path.sep + options.sample + "_R2.fastq.gz"
	
	#step 1
	#check the existence of the output files
	#skip existing files if override is turned off
	"""
	if (options.override or not os.path.exists(r1) or not os.path.exists(r2) or not os.path.exists(alignments_mtdna) or not os.path.exists(alignments_ndna)):
		#logging
		logging_file = outpath + os.path.sep + options.sample + ".summary"
		LOG = open(logging_file, "w")
		#filter reads according to known probe sequences, output file names containing filtered reads
		r1, r2 = retrieve_tags(options.sample, fastq_r1, fastq_r2, probe, options.probe_reverse, outpath, True, options.output_off_target, options.min_qual, options.trim)
		#map filtered reads to the reference genomes
		alignments_mtdna, alignments_ndna = retrieve_alignments(options.sample, r1, r2, probe, outpath, options.genome, options.mtdna, options.mtdna_offset)
	"""
	#skip the read de-multiplexing step if the output file already exists
	if (fastq_r1 is not None and str(fastq_r1) != r1 and (options.override or not os.path.exists(r1) or not os.path.exists(r2))):
		#logging
		logging_file = outpath + os.path.sep + options.sample + ".summary"
		LOG = open(logging_file, "w")
		#filter reads according to known probe sequences, output file names containing filtered reads
		r1, r2 = retrieve_tags(options.sample, fastq_r1, fastq_r2, probe, options.probe_reverse, outpath, True, options.output_off_target, options.min_qual, options.trim)
	
	alignments_mtdna = outpath + os.path.sep + options.sample + ".mtdna.bam"
	alignments_ndna = outpath + os.path.sep + options.sample + ".ndna.bam"
	#step 2
	#skip the alignment step if the output file already exists
	if (os.path.exists(r1) and (options.override or not os.path.exists(alignments_mtdna) or not os.path.exists(alignments_ndna))):
		#logging
		if (LOG is None):
			logging_file = outpath + os.path.sep + options.sample + ".summary"
			LOG = open(logging_file, "a")
		#map filtered reads to the reference genomes
		alignments_mtdna, alignments_ndna = retrieve_alignments(options.sample, r1, r2, probe, outpath, options.genome, options.mtdna, options.mtdna_offset)
	
	barcode_file = outpath + os.path.sep + options.sample + ".barcode"
	#build consensus sequence
	#output files for the mtDNA and nDNA consensus reads
	if (options.consensus_outpath):
		outpath = options.consensus_outpath
	
	#step 3
	#check the existence of the output files
	#skip existing files if override is turned off
	re_alignments_mtdna = outpath + os.path.sep + options.sample + ".mtdna.sorted.realign"+(".recal.bam" if not options.no_recalibration else ".bam")
	re_alignments_ndna = outpath + os.path.sep + options.sample + ".ndna.sorted.realign"+(".recal.bam" if not options.no_recalibration else ".bam")
	if (os.path.exists(alignments_mtdna) and (options.override or not os.path.exists(re_alignments_mtdna) or not os.path.exists(re_alignments_ndna))):
			re_alignments_mtdna, re_alignments_ndna = retrieve_re_alignments(options.sample, alignments_mtdna, alignments_ndna, outpath, options.genome, options.mtdna, not options.no_recalibration)
	
	#read alternate mtDNA sequences. e.g., sequences of numts, alignment artifacts in hsd format
	if (options.numts):
		mtdna_alt_seq = retrieve_alternate_sequences(options.numts, options.mtdna_offset)
	else:
		mtdna_alt_seq = {}
	consensus_mtdna = outpath + os.path.sep + options.sample + ".mtdna.consensus.bam"
	consensus_ndna = outpath + os.path.sep + options.sample + ".ndna.consensus.bam"
	
	#step 4
	#check the existence of the output files
	#skip existing files if override is turned off
	if (options.override or not os.path.exists(consensus_mtdna) or not os.path.exists(consensus_ndna)):
		#logging
		if (logging_file != outpath + os.path.sep + options.sample + ".summary"):
			logging_file = outpath + os.path.sep + options.sample + ".summary"
			if (LOG):
				#close the current fh for logging
				LOG.close()
			LOG = open(logging_file, "w")
		
		consensus_mtdna, consensus_ndna = retrieve_consensus_reads(options.sample, re_alignments_mtdna, re_alignments_ndna, probe, barcode_file, outpath, options.genome, options.mtdna, options.mtdna_offset, mtdna_alt_seq, not options.no_recalibration, options.clip_end, options.consensus_nm, options.consensus_xf, options.consensus_gap, options.consensus_qual, options.consensus_lqbase)
	
	#close the fh for logging
	if (LOG):
		LOG.close()

def main():
	run(sys.argv[0], sys.argv[1:])
	
if __name__ == "__main__":
	main()