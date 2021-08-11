#############################################################################
#
#Updated by: Yiqin Wang
#
##############################################################################

import os
import re
from FileIO import FileIO
#from collections import defaultdict

def flip(x):
	#return a reverve complementary string
	return "".join(reversed([{"A":"T","T":"A","G":"C","C":"G","[":"]","]":"["}.get(i,i) for i in x]))

def phred(score):
	#return the integer phred score
	return map(lambda i: ord(i)-33, score)

#set up the table for converting phred score to the corresponding error probability
phred2prob = []
for i in range(0, 223):
	phred2prob.append(10**(-i/10.0))

def prob(phred_score):
	#convert a phred score to the error probability 
	global phred2prob
	return phred2prob[phred_score]
	
def trim(qual, limit = 0.05):
	#return the left n base with the highest quaily
	#modified Richard Mott's trimming algorithm and Li He's seqtk
	s = s_max = 0.0
	end = 0
	for n, q in enumerate(qual):
		p = prob(q)
		s += limit-p
		if (s < 0):
			s = 0
		if (s >= s_max):
			s_max = s
			end = n + 1
	return end

class FastaIO:
	""" 
	This is a class for reading fasta files
	
	Attributes
	----------
	name: the names of the sequences in the fasta file
	__sequence__: a dict of sequence read from the fasta file
	__seq_index__: a dict of sequence index from the fasta .fai file
	__seq_file__: file handle of the fasta file
	
	"""
	
	def __init__(self, fasta_file, index_file = None):
		""" 
		the __init__ method
		
		read sequence information from the fasta_file or sequence indices from the index file
		
		Arguments
		----------
		fasta_file: the path to the fasta file
		index_file: the path to the fasta index file; if not specified, will look for file names ${fasta_file}.fai
		
		"""
		#read index if exists
		if (index_file is None):
			index_file = fasta_file + ".fai"
		if (os.path.exists(index_file)):
			#read index into memory
			with open(index_file, "rb") as fai_fh:
				self.__load_index__(fai_fh)
			self.__seq_file__ = open(fasta_file, "rb")
		else:
			#read all sequence data into memory
			with open(fasta_file, "rb") as fasta_fh:
				self.__load_sequence__(fasta_fh)
			self.__seq_file__ = None
	
	def __load_index__(self, fai_fh):
		""" 
		read sequence indices from the index file
		
		update __seq_index__ with a dict of sequence index ({name:sequence_index})
		
		Arguments
		----------
		fai_fh: file handle of the fasta index file
		
		"""
		n = 0
		seq_index = {}
		seq_name = []
		for line in fai_fh:
			n += 1
			i = line.strip().split()
			#name, length, offset, number of nt per line, number of character per line
			assert len(i) == 5, "Invalid index file line '%s' at line %d" % (line,n)
			name, seq_len, offset, line_nt, line_length = i[0], int(i[1]), int(i[2]), int(i[3]), int(i[4])
			assert name and seq_index.get(name) is None, "Invalid sequence name '%s' at line %d" % (name,n)
			seq_index[name] = [seq_len, offset, line_nt, line_length]
			seq_name.append(name)
		self.name = seq_name
		self.__sequence__ = None
		self.__seq_index__ = seq_index
	
	def __load_sequence__(self, fasta_fh, skip_line = 0):
		""" 
		read sequence information from the fasta file
		
		update __sequence__ with a dict of sequence information read from the file ({name:sequence})
		
		Arguments
		----------
		fasta_fh: file handle of the fasta file
		skip_line: no. of lines to skip in the fasta file
		
		"""
		sequence = {}
		seq_name = []
		seq_line = []
		#align_state = ""
		name=""
		n = 0
		for line in fasta_fh:
			n += 1
			if (skip_line > 0):
				skip_line -= 1
				continue
			line = line.strip()
			if (not line):
				continue
			if (line[0] == ">"):
				#head line
				if (seq_line):
					sequence[name] = "".join(seq_line)
					seq_line = []
				name = line.split()[0][1:]
				assert name and sequence.get(name) is None, "Invalid sequence name '%s' at line %d" % (name,n)
				seq_name.append(name)
			else:
				seq_line.append(line)
		if (seq_line):
			sequence[name] = "".join(seq_line)
		self.name = seq_name
		self.__sequence__ = sequence
		self.__seq_index__ = None
	
	def extract(self, cord):
		""" 
		extract a segment of a certain sequence 
		
		Arguments
		----------
		cord: sequence name and positions, e.g., name, name:start-end (start and end are 1-based, represent the region of [start,end))
		
		Returns
		----------
		see fetch(...)
		
		"""
		
		seq = cord.split(":")
		if (len(seq) == 1):
			seq = seq[0]
			start = end = None
		else:
			seq, pos = seq[:2]
			pos = pos.split("-")
			if (len(pos) == 1):
				start = int(pos[0])
				end = None
			else:
				start = int(pos[0])
				end = int(pos[1])
		return self.fetch(seq, start, end)
	
	def fetch(self, seq, start = None, end = None):
		""" 
		extract a segment of a certain sequence 
		
		Arguments
		----------
		seq: name of the sequence
		start: start position
		end: end position
		start and end are 1-based, represent the region of [start,end).
		
		Returns
		----------
		a tuple of
		1. sequence name (the same as the argument seq)
		2. start (the offset position in the original sequence, 0-based)
		3. the sequence of the region (None if the region is not valid or not included)
		
		"""
		
		if (self.__sequence__ is not None):
			sequence = self.__sequence__.get(seq)
			if (start is None and end is None):
				start = 0
				end = len(sequence) 
			elif (end is None):
				start = max(0, start-1)
				end = len(sequence) 
			else:
				start = max(0, start-1)
				end = min(end,len(sequence))
			sequence = sequence[start:end]
		elif (self.__seq_index__ is not None):
			index = self.__seq_index__.get(seq)
			seq_len, offset, line_nt, line_length = index
			if (start is None or start < 1):
				start = 1
			if (end is None or end > seq_len):
				end = seq_len + 1
			if (end <= start):
				end = start + 1
			start -= 1
			end -= 1
			pos = start/line_nt*line_length + start%line_nt
			length = end/line_nt*line_length+end%line_nt - pos
			self.__seq_file__.seek(offset + pos)
			sequence = self.__seq_file__.read(length)
			sequence = "".join(sequence.split())
			assert len(sequence) == (end-start), "Invalid sequence '%s' read from fasta file <%s>" % (seq,self.__seq_file__.name)
		else:
			sequence = None
		#start is zero-based
		return (seq, start, sequence)
	
	def list(self):
		#list all sequence names
		return self.__sequence__.keys()

class FastqIO:
	""" 
	This is a class for reading fastq files
	
	Attributes
	----------
	r1: file handle of R1 fastq file
	r2: file handle of R2 fastq file
	r1_phred: the cumulative phred score at each position in R1; list of length 500; record while reading R1
	r2_phred: the cumulative phred score at each position in R2; list of length 500; record while reading R2
	r1_length: the length distribution of R1; dict of read length and count; record while reading R1
	r2_length: the length distribution of R2; dict of read length and count; record while reading R2
	qc: true/false if recording read information of R1 and R2
		
	"""
	
	def __init__(self, fastq_r1, fastq_r2, quality_check = False):
		""" 
		the __init__ method
		
		Arguments
		----------
		fastq_r1: the path to the R1 fastq file
		fastq_r2: the path to the R2 fastq file
		quality_check: true/false if recording read information of R1 and R2
		
		"""
		self.r1 = FileIO(fastq_r1, "r")
		self.r2 = FileIO(fastq_r2, "r")
		self.r1_phred = [0,]*500 #500 bp, long enough for most NGSs
		self.r2_phred = [0,]*500
		self.r1_length = {}
		self.r2_length = {}
		#self.n = 0
		self.qc = quality_check
	
	def load(self, fastq_r1, fastq_r2):
		""" 
		open new R1 and R2 fastq files
		
		Arguments
		----------
		fastq_r1: the path to the R1 fastq file
		fastq_r2: the path to the R2 fastq file
		
		"""
		if (self.r1):
			self.r1.close()
		self.r1 = FileIO(fastq_r1, "r")
		if (self.r2):
			self.r2.close()
		self.r2 = FileIO(fastq_r2, "r2")
	
	def __iter__(self):
		return self
	
	def __quality__(self, phred, length):
		""" 
		compute the average base quality for each position in a list
		
		Arguments
		----------
		phred: a list of cumulative phred scores
		length: a dict of read length and count
		
		Returns
		----------
		a list of the average base quality for each position in the phred list
		
		"""
		length = [(i, length[i]) for i in sorted(length.keys())]
		phred_avg = [0,]*500
		for i, v in enumerate(length):
			for j in range(i):
				phred_avg[j] += v[1]
		for i in range(len(phred_avg)):
			if (phred_avg):
				phred_avg[i] = phred[i]/float(phred_avg[i])
		return phred_avg
	
	def quality(self):
		""" 
		compute the average base quality for each position in R1 and R2
		
		use __quality__(...) and the read quality and length information were recorded. 
				
		Returns
		----------
		a tuple of
		1. a list of the average base quality for each position in R1
		2. a list of the average base quality for each position in R2
		
		"""
		r1_quality = self.__quality__(self.r1_phred, self.r1_length)
		r2_quality = self.__quality__(self.r2_phred, self.r2_length)
		return (r1_quality, r2_quality)
	
	def next(self):
		""" 
		the next function for iterating each read pair in the fastq files provided
		
		Returns
		----------
		a tuple of
		1. R1: a tuple of read name, read sequence, read quality
		2. R2: a tuple of read name, read sequence, read quality
		
		"""
		state = 0
		n_r1 = 0
		n_r2 = 0
		while (True):
			#Reads in R1 and R2 should be in pair and in the same order.
			while (True):
				line_r1 = self.r1.readline()
				if (not line_r1):
					raise StopIteration
				line_r1 = line_r1.rstrip("\r\n")
				n_r1 += 1
				if (line_r1):
					break
			while (True):
				line_r2 = self.r2.readline()
				if (not line_r2):
					raise StopIteration
				line_r2 = line_r2.rstrip("\r\n")
				n_r2 += 1
				if (line_r2):
					break
			if (state == 0):
				#read name line
				if (line_r1[0] != "@"):
					raise ValueError("Error: R1 read name %s @ line %d" % (line_r1, n_r1))
				if (line_r2[0] != "@"):
					raise ValueError("Error: R2 read name %s @ line %d" % (line_r2, n_r2))
				name_r1 = line_r1
				name_r2 = line_r2
				state = 1
			elif (state == 1):
				#read line
				read_r1 = line_r1
				read_r2 = line_r2
				state = 2
			elif (state == 2):
				#info line
				if (line_r1[0] != "+"):
					raise ValueError("Error: R1 info %s @ line %d" % (line_r1, n_r1))
				if (line_r2[0] != "+"):
					raise ValueError("Error: R2 info %s @ line %d" % (line_r2, n_r2))
				state = 3
			elif (state == 3):
				state = 0
				#quality line
				if (len(read_r1) != len(line_r1)):
					raise ValueError("Error: R1 read length (%d) and quality length (%d) do not match" % (len(read_r1), len(line_r1)))
				if (len(read_r2) != len(line_r2)):
					raise ValueError("Error: R2 read length (%d) and quality length (%d) do not match" % (len(read_r2), len(line_r2)))
				quality_r1 = line_r1
				quality_r2 = line_r2
				if (self.qc):
					#obsolete
					#compute cumulative quality in R1 and R2
					phred_q1 = phred(quality_r1)
					for i,q in enumerate(phred_q1):
						self.r1_phred[i] += q
					phred_q2 = phred(quality_r2)
					for i,q in enumerate(phred_q2):
						self.r2_phred[i] += q
					#R1 and R2 length distribution
					self.r1_length[len(quality_r1)] += 1
					self.r2_length[len(quality_r2)] += 1
					return ((name_r1, read_r1, quality_r1, self.r1_phred), (name_r2, read_r2, quality_r2, self.r2_phred))
				else:
					return ((name_r1, read_r1, quality_r1), (name_r2, read_r2, quality_r2))

#class to parse sam files

class SAMFlag:
	""" 
	a class to parse the mapping flag in the sam/bam file
	
	Attributes
	----------
	stat: int stat
	true/false status for
	proper_pair
	unmapped
	mate_unmapped
	strand_reverse
	mate_strand_reverse
	first_segment
	last_segment
	primary
	duplicate
		
	"""
	def __init__(self, stat):
		""" 
		the __init__ method
		
		Arguments
		----------
		stat: string or int stat of the mapping flag in the sam/bam file 
		"""
		stat = int(stat)
		self.stat = stat
		self.proper_pair = (stat & 0x2 > 0)
		self.unmapped = (stat & 0x4 > 0)
		self.mate_unmapped = (stat & 0x8 > 0) 
		self.strand_reverse = (stat & 0x10 > 0)
		self.mate_strand_reverse = (stat & 0x20 > 0)
		self.first_segment = (stat & 0x40 > 0)
		self.last_segment = (stat & 0x80 > 0)
		self.primary = (stat & 0x100 == 0 and stat & 0x800 == 0)
		self.duplicate = (stat & 0x400 > 0)

CIGAR_PATTERN = re.compile("(\d*)(.)")
class SAMCigar:
	""" 
	a class to parse the cigar string
	
	Attributes
	----------
	cigar: a tuple representing the cigar string
	align_seq: base of aligned reads at each position of the reference
	align_qual: base quality of aligned reads at each position of the reference
	ref_start: the start position in the reference
	ref_len: the mapped length in the reference
	ref_end: the end position in the reference
		
	"""
	def __init__(self, cigar, ref_start = 0, seq = None, qual = None):
		""" 
		a class to parse the cigar string
		
		Attributes
		----------
		cigar: a tuple representing the cigar string
		ref_start: the start position in the reference
		seq: sequence string
		qual: quality string
			
		"""
		self.ref_start = ref_start
		self.cigar = [(int(i), j) for i, j in re.findall(CIGAR_PATTERN, cigar)]
		cur = 0
		if (seq and qual):
			if (len(seq) != len(qual)):
				raise ValueError("seq string and qual string do not have the same length <%s> <%s>" % (seq, qual))
			#rebuild alignment according to the cigar string
			align_seq = []
			align_qual = []
			last_stat = ""
			for i, j in self.cigar:
				if (j in ("M","X","=")):
					#match or mismatch
					if (last_stat == "I"):
						#modify the last aligned position
						align_seq[-1] += seq[cur]
						align_qual[-1].extend(phred(qual[cur]))
						cur += 1
						i -= 1
					align_seq.extend(list(seq[cur:(cur+i)]))
					align_qual.extend(phred(qual[cur:(cur+i)]))
					cur += i
				elif (j == "H"):
					#hard clip
					pass
				elif (j == "S"):
					#soft clip
					cur += i
				elif (j == "I"):
					#insertion
					align_seq.append(seq[cur:(cur+i)])
					align_qual.append(phred(qual[cur:(cur+i)]))
					cur += i
				elif (j == "D"):
					#deletion
					align_seq.extend([".",]*i)
					align_qual.extend([None,]*i)
				else:
					#other cigar stats, such as N in mRNA alignment
					raise NotImplementedError("The current version does not support <%s> in cigar string" % j)
				last_stat = j
			self.align_seq = align_seq
			self.align_qual = align_qual
			self.ref_len = len(align_seq)
			if (last_stat == "I" and len(align_qual[-1])==1):
				#convert 1-length insertion at the end to a mismatch 
				align_qual[-1] = align_qual[-1][0]
		else:
			self.align_seq = None
			self.align_qual = None
			self.ref_len = 0
		self.ref_end = self.ref_start + self.ref_len
		
	def fetch(self, ref_pos):
		""" 
		get the base mapped to the position in the reference
		
		Attributes
		----------
		ref_pos: the position in the reference
		
		Returns
		----------
		a tuple of the base and the base quality
		"""
		if (ref_pos >= self.ref_start and ref_pos < self.ref_end):
			p = ref_pos - self.ref_start
			return (self.align_seq[p], self.align_qual[p])
		else:
			return ("N", 0)
	
	def clip_head(self, size):
		""" 
		trim the aligned reads for a certain length from the start
		
		Attributes
		----------
		size: number of aligned positions to trim
		
		Returns
		----------
		self
		"""
		if (size >= self.ref_len):
			raise ValueError("Cannot clip segments longer than alignment length")
		#update the aligned reads and qualty
		self.align_seq = self.align_seq[size:]
		self.align_qual = self.align_qual[size:]
		#update the reference start and length
		self.ref_start += size
		self.ref_len -= size
		return self
	
	def clip_tail(self, size):
		""" 
		trim the aligned reads for a certain length from the end
		
		Attributes
		----------
		size: number of aligned positions to trim
		
		Returns
		----------
		self
		"""
		if (size >= self.ref_len):
			raise ValueError("Cannot clip segments longer than alignment length")
		#update the aligned reads and qualty
		self.align_seq = self.align_seq[:-size]
		self.align_qual = self.align_qual[:-size]
		#update the reference start and length
		self.ref_end -= size
		self.ref_len -= size
		return self
