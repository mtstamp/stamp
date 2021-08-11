#!/usr/bin/python
#############################################################################
#Package: stamp
#Updated By Yiqin Wang 
#Usage:
#main
#call mtDNA variation from mapped reads
#
##############################################################################

import sys
import re
import copy

from collections import Counter, defaultdict
from argparse import ArgumentParser

#math libraries
import random
from math import sqrt, log10
import scipy.stats as stats
import scipy.optimize as optimize

#basic file readers and parsers
from FileIO import FileIO

def formatAllele(allele):
	#format allele string; obsolete
	return ";".join(["%s=%d" % (j, i) for i, j in allele.items()])

def addCounts(c1, c2):
	#sum up values in two dicts; obsolete
	for i, n in c2.items():
		if (c1.get(i) is not None):
			c1[i] += n
		else:
			c1[i] = n
	return c1

def stringList(l):
	#convert all items in a list to strings
	return [str(i) for i in l]

def readRQ(readqualities):
	"""
	convert read quality string to a list of phred scores
	
	Arguments
	----------
	a string of read quality stored in the alignment file
	
	Returns
	----------
	a list of phred scores
	
	"""
	return [ord(i)-33 for i in readqualities]

def readRB(readbases):
	"""
	a function to read base information in the pileup file
	
	Arguments
	----------
	a string of read bases stored in the pileup file
	
	Returns
	----------
	a tuple of a string of single bases, a dict of deletions
	
	"""
	readbases_old = readbases 
	readbases=re.sub('\^.?','',readbases)
	readbases=re.sub('\$','',readbases)
	readbases=readbases.replace("$","")
	#indel string
	indel_pat = re.compile(r"([-\+])(\d+)\D?")
	#indels
	deletions = {}
	insertions = {}
	singlent = ""
	#current position
	offset = 0
	while (True):
		ret = indel_pat.search(readbases, offset)
		if (ret is not None):
			type = ret.groups()[0]
			length = int(ret.groups()[1])
			start = ret.start() + 1 + len(ret.groups()[1])
			#store the indel string
			seq = readbases[start:(start+length)]
			if (type == "-"):
				#deletion
				if (deletions.get(seq) is None):
					deletions[seq] = 1
				else:
					deletions[seq] += 1
			else:
				#insertion
				if (insertions.get(seq) is None):
					insertions[seq] = 1
				else:
					insertions[seq] += 1
			#single bases left
			singlent += readbases[offset:ret.start()]
			offset = start+length
		else:
			singlent += readbases[offset:]
			break
	return singlent, insertions, deletions

def readMpileup(fh, sample_idx = 0, sample_name = "s", quality_min = 20, rm_rate = 0.5, retain_quality = True):
	"""
	parse mpileup file for base information
	
	Arguments
	----------
	fh: file handle of the pileup file
	sample_idx: index of the samples to extract (start from the 4th column)
	sample_name: the names of the samples to extract (list)
	quality_min: the minimum of the BAQ used to count bases
	rm_rate: the minimum rate of bases retained after BAQ filtering
	retain_quality: return the quality list
	
	Returns
	----------
	an iterator of tuples containing:
	1. sample index
	2. sample name
	3. chromosome
	4. position
	5. letter of base on the forward strand
	6. total depth of reads
	7. base count (a dict of keys representing bases)
	8. insertion count
	9. deletion count
	10. associated BAQ list (optional)
	
	"""
	if (isinstance(sample_idx,int)):
		sample_idx = [sample_idx,]
	if (isinstance(sample_name,str)):
		sample_name = [sample_name,]
	assert len(sample_idx) == len(sample_name), "Sample ids and names don't match."
	line = fh.readline()
	line = line.rstrip("\r\n").split("\t")
	#number of samples in the mpileup file
	nsample = (len(line)-3)/3 
	assert nsample > max(sample_idx) and len(line) % 3 == 0, "There is not enough columns in the mpileup file <%s>" % fh.name
	fh.seek(0)
	for line in fh:
		line = line.rstrip("\r\n")
		if (not line):
			continue
		line = line.split("\t")
		assert len(line) == nsample*3 + 3, "The number of columns in the mpileup file <%s> doesn't match" % fh.name
		chr = line[0]
		pos = int(line[1])
		ref_fwd = line[2].upper() #reference allele on the forward strand
		ref_rev = line[2].lower() #reference allele on the reverse strand
		for i in range(len(sample_idx)):
			id = sample_idx[i]
			name = sample_name[i]
			depth = int(line[3*id+3])
			readqualities = readRQ(line[3*id+5])
			readbases, insertions, deletions = readRB(line[3*id+4])
			#print len(readbases) ,len(insertions), len(deletions), len(readqualities)
			#if (len(readbases) == depth and (insertions or deletions)):
			#	print depth, len(readbases) ,len(insertions), len(deletions)
			if (len(readbases) != len(readqualities)):
				print >>sys.stderr, "Mismatched read base and quality", line
				print >>sys.stderr, "Mismatched read base and quality", len(readbases) ,len(insertions), len(deletions), len(readqualities)
			#assert len(readbases) == len(readqualities)
			#assert len(readbases)+len(insertions)+len(deletions) == len(readqualities)
			#print len(readbases)+len(insertions)+len(deletions), len(readqualities)
			nt_count = {"A":0,"T":0,"G":0,"C":0,"a":0,"t":0,"c":0,"g":0}
			rq = []
			for n, q in zip(readbases, readqualities):
				if (q >= quality_min):
					if (n == "."):
						nt_count[ref_fwd] += 1
					elif (n == ","):
						nt_count[ref_rev] += 1
					elif (nt_count.get(n) is not None):
						nt_count[n] += 1
					else:
						continue
					rq.append((n,q))
			ins_count = dict(Counter(insertions)) #to update
			del_count = dict(Counter(deletions))
			if (depth and sum(nt_count.values())/float(depth) >= rm_rate):
				#do not return if no reads are found or the rate of bases left is lower than rm_rate
				yield (i, name, chr, pos, ref_fwd, depth, nt_count, ins_count, del_count, rq if (retain_quality) else None)

class MTSite:
	""" 
	This is a class for processing a single variant site
	
	Attributes
	----------
	chr: chromosome (string)
	position: int position
	depth: int depth of reads
	ref: reference allele (string)
	allele_count: a dict of allele count (readMpileup)
	ins_count: a dict of insertion count (readMpileup)
	del_count: a dict of deletion count (readMpileup)
	readsquality: a list of base quality (readMpileup)
	is_heteroplasmy: 0/1 status
	is_substitution: 0/1 status
	heteroplasmy: a list of heteroplasmy status (see __str__)
	"""
	
	def __init__(self, chr, pos, depth, ref, allele_count, ins_count, del_count, readsquality = None):
		""" 
		the __init__ method
		
		Arguments
		----------
		see class Attributes 
		"""
	
		self.chr = chr
		self.pos = pos
		self.depth = depth
		self.ref = ref
		self.allele_count = allele_count
		self.ins_count = ins_count
		self.del_count = del_count
		self.readsquality = readsquality
		self.is_heteroplasmy = 0
		self.is_substitution = 0
		self.heteroplasmy = None
	
	def __str__(self):
		"""
		returns the string representation of the object
		
		Returns
		----------
		a tab-separated string containing
		1. chromosome
		2. position
		3. total depth of reads
		4. depth of reads on the forward strand (after QC)
		5. depth of reads on the reverse strand (after QC)
		6. the major allele
		7. the major allele on the forward strand
		8. the major allele on the reverse strand
		9. is a heteroplasmy
		10. is a substitution
		11-18. heteroplasmy information
		11. minor allele
		12. the fraction of the minor allele
		13. the MLE of the fraction of the minor allele
		14. the Log-likelihood of the MLE fraction
		15. 2.5% lower bound of the MLE fraction
		16. 97.5% upper bound of the MLE fraction
		17. P value from Fisher's test for difference from the minimum fraction
		18. P value for strand bias of variant fractions
		"""
		heteroplasmy = ["",]*8 if self.heteroplasmy is None else self.heteroplasmy
		return "\t".join(stringList([self.chr, self.pos, self.ref, self.depth, self.depth_fwd, self.depth_rev, self.allele] + \
						 self.allele_fwd + self.allele_rev + [self.is_heteroplasmy, self.is_substitution] + heteroplasmy))
	
	def string(self):
		#same as __str__
		return self.__str__() 
	
	def callAllele(self, min_depth, min_depth_fwd, min_depth_rev, min_minor_depth, min_minor_depth_fwd, min_minor_depth_rev, min_het_freq):
		""" 
		identify the major and the minor alleles at a position
		
		Arguments
		----------
		the requirements to call major and minor alleles
		min_depth: the minimum read depth
		min_depth_fwd: the minimum read depth on the forward strand
		min_depth_rev: the minimum read depth on the reverse strand
		min_minor_depth: the minimum read depth of the minor allele
		min_minor_depth_fwd: the minimum read depth of the minor allele on the forward strand
		min_minor_depth_rev: the minimum read depth of the minor allele on the reverse strand
		min_het_freq: the minimum fraction of the minor allele
		
		Returns
		----------
		None
		
		"""
		#get allele count
		allele_fwd = self.allele_fwd = [self.allele_count.get(i, 0) for i in ["A", "T", "C", "G"]]
		allele_rev = self.allele_rev = [self.allele_count.get(i, 0) for i in ["a", "t", "c", "g"]]
		self.allele_all = [allele_fwd[i] + allele_rev[i] for i in range(4)]
		#read depth
		depth_fwd = sum(allele_fwd)
		depth_rev = sum(allele_rev)
		self.depth_fwd = depth_fwd
		self.depth_rev = depth_rev
		self.heteroplasmy = None
		if (depth_fwd + depth_rev < min_depth):
			self.allele = "N"
		else:
			allele_all =  [[self.allele_all[i],i] for i in range(4)]
			allele_all.sort(reverse = True)
			major_allele = ["A","T","C","G"][allele_all[0][1]]
			minor_allele = ["A","T","C","G"][allele_all[1][1]]
			if (allele_all[1][0] > 0):
				p_fisher = None
				p_sb = None
				a1 = allele_all[0][0]
				a2 = allele_all[1][0]
				if (depth_fwd >= min_depth_fwd and depth_rev >= min_depth_rev):
					if (a2 >= min_minor_depth):
						a12 = allele_fwd[allele_all[1][1]]
						a22 = allele_rev[allele_all[1][1]]
						if (a12 >= min_minor_depth_fwd and a22 >= min_minor_depth_rev):
							#estimate heteroplasmy frequency
							freq_ratio, freq_mle, freq_llr, freq_low, freq_high = self.estimateHeteroplamy(major_allele, a1, minor_allele, a2)
							
							#test whether freq is significantly hight than the minimum frequency.
							odds, p_fisher = stats.fisher_exact([[(a1+a2)*(1-min_het_freq), (a1+a2)*min_het_freq],[a1,a2]], "greater")
							
							#test strand bias
							a11 = allele_fwd[allele_all[0][1]]
							a21 = allele_rev[allele_all[0][1]]
							odds, p_sb = stats.fisher_exact([[a11,a12],[a21,a22]])
							#output heteroplasmy
							self.heteroplasmy = [minor_allele, freq_ratio, freq_mle, freq_llr, freq_low, freq_high, p_fisher, p_sb]
							if (freq_ratio > min_het_freq):
								self.is_heteroplasmy = 1
					else:
						pass
			self.allele = major_allele
			if (self.allele != self.ref):
				self.is_substitution = 1
	
	def downSample(self, depth):
		""" 
		Down sample reads to the depth specified.
		Reads were stored in the reads quality attribute.
		
		Arguments
		----------
		depth: int depth of read; should be smaller than the current read depth
		
		Returns
		----------
		None
		
		"""
		assert self.readsquality is not None, "No reads and qualities stored for down sampling" 
		if (len(self.readsquality) <= depth):
			#return a copy of self if the current depth is already smaller than the given depth
			return copy.deepcopy(self)
		#sampling depth reads
		rq = random.sample(self.readsquality, depth)
		nt_count = {"A":0,"T":0,"G":0,"C":0,"a":0,"t":0,"c":0,"g":0}
		ref_fwd = self.ref.upper()
		ref_rev = self.ref.lower()
		for i, j in rq:
			if (i == "."):
				nt_count[ref_fwd] += 1
			elif (i == ","):
				nt_count[ref_rev] += 1
			elif (nt_count.get(i) is not None):
				nt_count[i] += 1
			else:
				continue
		#copy the current attributes
		ret = copy.deepcopy(self)
		#update the read depth
		ret.allele_count = nt_count
		#update the reads quality list
		ret.readsquality = rq
		return ret
	
	def estimateHeteroplamy(self, major, a1, minor, a2):
		""" 
		estimate the minor allele fraction of a heteroplasmy
		if the original bases and BAQ are stored, a log-likelihood estimation of the fraction is computed
		
		Arguments
		----------
		major: major allele
		a1: int major allele count
		minor: minor allele
		a2: int minor allele count
		
		Returns
		----------
		a tuple of 
		1. fraction: a2/(a1+a2)
		2. mle estimated fraction
		3. log-likelihood of the mle estimation
		4 and 5. 95% CI of the fraction
		
		"""
		freq_ratio = a2 / float(a1+a2)
		freq_mle = ""
		llr = ""
		if (self.readsquality):
			quality_major = defaultdict(lambda:0)
			quality_minor = defaultdict(lambda:0)
			for a, q in self.readsquality:
				if (a == "," or a == "."):
					a = self.ref
				else:
					a = a.upper()
				#convert phred score to a probality
				q = 10**(-q/10.0)
				if (a == major):
					quality_major[q] += 1
				elif (a == minor):
					quality_minor[q] += 1
			def likelihood(f):
				est = 1.0
				for q, n in quality_major.items():
					est *= (((1-f)*q + f*(1-q)))**n
				for q, n in quality_minor.items():
					est *= (((1-f)*(1-q) + f*q))**n
				#-log(ext)
				return(-est)
			def loglikelihood(f):
				est = 0.0
				for q, n in quality_major.items():
					est += n*log10(((1-f)*q + f*(1-q)))
				for q, n in quality_minor.items():
					est += n*log10(((1-f)*(1-q) + f*q))
				#-log(ext)
				return(-est)
			ret = optimize.minimize_scalar(loglikelihood, bounds=(0.0, 1.0), method='bounded')
			if (ret.success):
				freq = freq_mle = 1-ret.x
				l0 = max(-loglikelihood(1.0),-loglikelihood(0.0))
				l1 = -loglikelihood(ret.x)
				if (l0 == 0):
					llr = float('Inf')
				else:
					llr = l1-l0
					#llr = log10(l1)-log10(l0) #already log10 transformed
			else:
				freq = freq_ratio
		else:
			freq = freq_ratio
		num = a1+a2
		sd = sqrt((freq*(1-freq))/num)
		#95% CI of the fraction
		freq_low = max(freq - 1.96*sd,0.0)
		freq_high = min(freq + 1.96*sd,1.0)
		return (freq_ratio, freq_mle, llr, freq_low, freq_high) 
	
	def compareAllele(self, sample2, ref_allele = None):
		""" 
		compare alleles with another sample 
		
		Arguments
		----------
		sample2: (MTSite)
		ref_allele: give a reference allele (optional)
		
		Returns
		----------
		a status string of 
		Allele difference/Same Allele/Comparable Heteroplasmy/Decreased Heteroplasmy/Increased Heteroplasmy
		
		"""
		if (self.allele == "N" or sample2.allele == "N"):
			return ""
		if (self.allele != sample2.allele and self.is_heteroplasmy == 0 and sample2.is_heteroplasmy == 0):
			return "Allele difference"
		p_fisher = 0.05
		if (ref_allele is not None and ref_allele in ["A","T","C","G"]):
			ref_allele = ["A","T","C","G"].index(ref_allele)
		else:
			ref_allele = None
		if (self.is_heteroplasmy):
			a1 = ["A","T","C","G"].index(self.allele)
			a2 = ["A","T","C","G"].index(self.heteroplasmy[0])
			if (ref_allele is not None and a2 == ref_allele):
				a2 = a1
				a1 = ref_allele
			if (sample2.is_heteroplasmy):
				ret = self.compareHeteroplasmy(sample2, a1, a2, p_fisher)
			else:
				#a22 = ["A","T","C","G"].index(sample2.heteroplasmy[0])
				#if (a2 != a22):
				if (sample2.allele == ["A","T","C","G"][a1]):
					ret = -1
				elif (sample2.allele == ["A","T","C","G"][a2]):
					ret = 1
				else:
					return "Allele difference"
				#ret = self.compareHeteroplasmy(sample2, a1, a2, p_fisher)
				#else:
				#ret = self.compareHeteroplasmy(sample2, a1, a2, p_fisher)
		elif (sample2.is_heteroplasmy):
			a1 = ["A","T","C","G"].index(sample2.allele)
			a2 = ["A","T","C","G"].index(sample2.heteroplasmy[0])
			if (ref_allele is not None and a2 == ref_allele):
				a2 = a1
				a1 = ref_allele
			if (self.allele == ["A","T","C","G"][a1]):
				ret = 1
			elif (self.allele == ["A","T","C","G"][a2]):
				ret = -1
			else:
				return "Allele difference"
			#ret = self.compareHeteroplasmy(sample2, a1, a2, p_fisher)
		else:
			return "Same Allele"
		if (ret == 0):
			return "Comparable Heteroplasmy"
		elif (ret == -1):
			return "Decreased Heteroplasmy"
		else:
			return "Increased Heteroplasmy"
	
	def compareHeteroplasmy(self, sample2, a1, a2, p = 0.05):
		""" 
		compare heteroplasmic fraction with another sample 
		
		Arguments
		----------
		sample2: (MTSite)
		a1: major allele
		a2: minor allele'
		p: significance level
		
		Returns
		----------
		-1: sample2 has a higher fraction
		1: sample1 has a higher fraction
		0: no difference
		"""
		odds, p_fisher = stats.fisher_exact([[self.allele_all[a1], self.allele_all[a2]],\
									   [sample2.allele_all[a1], sample2.allele_all[a2]]], "less")
		if (p_fisher < p):
			return -1
		odds, p_fisher = stats.fisher_exact([[self.allele_all[a1], self.allele_all[a2]],\
									   [sample2.allele_all[a1], sample2.allele_all[a2]]], "greater")
		if (p_fisher < p):
			return 1
		return 0

class MTScan:
	""" 
	This is a class for processing an mpileup file
	
	Attributes
	----------
	fh: file handle
	fh_to_close: true/false
	sample: sample index in the mpileup file to extract (list)
	name: names of the samples (list)
	
	QC filters to call heteroplasmies
	arguments used in MTSite.callAllele(...)
	base_quality: the minimum base quality
	min_reads_rate: the minimum rate of bases with BAQ >= base_quality
	min_depth: the minimum read depth
	min_depth_fwd: the minimum read depth on the forward strand
	min_depth_rev: the minimum read depth on the reverse strand
	min_minor_depth: the minimum read depth of the minor allele
	min_minor_depth_fwd: the minimum read depth of the minor allele on the forward strand
	min_minor_depth_rev: the minimum read depth of the minor allele on the reverse strand
	min_het_freq: the minimum fraction of the minor allele
	
	"""
	def __init__(self, pileup_file, sample = 0, name = "s", base_quality = 20, min_reads_rate = 0.5, min_depth = 10, min_depth_fwd = 1, min_depth_rev = 1, min_minor_depth = 1, min_minor_depth_fwd = 1, min_minor_depth_rev = 1, min_het_freq = 0.01):
		""" 
		the __init__ method
		
		Arguments
		----------
		see class Attributes 
		"""
		if (isinstance(pileup_file, str)):
			self.fh = FileIO(pileup_file, "r")
			self.fh_to_close = True
		else:
			self.fh = pileup_file
			self.fh_to_close = False
		if (isinstance(name, str)):
			name = [name,]
		if (sample is None):
			sample = range(len(name))
		elif (isinstance(sample, int)):
			sample = [sample,]
		assert len(sample) == len(name), "Sample and Name should be of same length."
		self.sample = sample
		self.name = name
		self.base_quality = base_quality
		self.min_reads_rate = min_reads_rate
		self.min_depth = min_depth
		self.min_depth_fwd = min_depth_fwd
		self.min_depth_rev = min_depth_rev
		self.min_minor_depth = min_minor_depth
		self.min_minor_depth_fwd = min_minor_depth_fwd
		self.min_minor_depth_rev = min_minor_depth_rev
		self.min_het_freq = min_het_freq
	
	def __del__(self):
		#close the file handle
		if (self.fh_to_close):
			self.fh.close()
	
	def allSites(self, end = 16570):
		""" 
		call mtDNA variants at all sites
		
		Arguments
		----------
		end: the end position
		
		Returns
		----------
		an iterator of tuples containing
		1. lists of MTSite for samples indicated
		2. int position
		3. is a variant (true/false) 
		"""
		cur = 1
		is_var = False
		var = [None,]*len(self.sample)
		#iterate all pileup lines
		for line in  readMpileup(self.fh, self.sample, self.name, self.base_quality, self.min_reads_rate):
			idx, name, chr, pos, ref, depth, allele_count, ins_count, del_count, rq = line
			if (pos > cur):
				#new line
				yield var, cur, is_var
				#reset var
				var = [None,]*len(self.sample)
				is_var = False
				cur += 1
				#output empty lines
				for p in range(cur, pos):
					yield var, p, is_var
				#move cur to pos
				cur = pos
			#new site information
			site = MTSite(chr, pos, depth, ref, allele_count, ins_count, del_count, rq)
			#determine variant alleles
			site.callAllele(self.min_depth, self.min_depth_fwd, self.min_depth_rev, self.min_minor_depth, self.min_minor_depth_fwd, self.min_minor_depth_rev, self.min_het_freq)
			#determine variant status
			if (site.is_heteroplasmy or site.is_substitution):
				is_var = True
			#temporarily store this variant
			var[idx] = site
		yield var, cur, is_var
		var = [None,]*len(self.sample)
		is_var = False
		#output empty lines for the remaining sites
		for p in range(cur+1, end):
			yield var, p, is_var
	
	def varSites(self):
		""" 
		call mtDNA variants at only variant sites
		same as allSites(...) but only return sites where variants are found

		Returns
		----------
		an iterator of tuples containing
		1. lists of MTSite for samples indicated
		2. int position
		"""
		cur = -1
		is_var = False
		var = [None,]*len(self.sample)
		for line in  readMpileup(self.fh, self.sample, self.name, self.base_quality, self.min_reads_rate):
			idx, name, chr, pos, ref, depth, allele_count, ins_count, del_count, rq = line
			if (cur == -1):
				cur = pos
			if (pos != cur):
				if (is_var):
					#return when it is a variant
					yield var, cur
				var = [None,]*len(self.sample)
				is_var = False
				cur = pos
			#new site information 
			site = MTSite(chr, pos, depth, ref, allele_count, ins_count, del_count, rq)
			#determine variant alleles
			site.callAllele(self.min_depth, self.min_depth_fwd, self.min_depth_rev, self.min_minor_depth, self.min_minor_depth_fwd, self.min_minor_depth_rev, self.min_het_freq)
			#determine variant status
			if (site.is_heteroplasmy or site.is_substitution):
				is_var = True
			#temporarily store this variant
			var[idx] = site
		if (is_var):
			yield var, cur
	
	def reset(self):
		#set the current position in the file handle to the beginning of the file
		self.fh.seek(0)

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
	parser = ArgumentParser(prog=prog, usage="%(prog)s [-options] input output\n", version="%(prog)s v0.1.1")
	parser.add_argument("--mind", type = int, default = 10, dest = "min_depth", help="the minimum read depth to call variants")
	parser.add_argument("--mind-fwd", type = int, default = 0, dest = "min_depth_fwd", help="the minimum read depth on the forward strand to call variants")
	parser.add_argument("--mind-rev", type = int, default = 0, dest = "min_depth_rev", help="the minimum read depth on the reverse strand to call variants")
	parser.add_argument("--minh", type = int, default = 2, dest = "min_minor_depth", help="the minimum number of minor alleles to call heteroplasmies")
	parser.add_argument("--minh-fwd", type = int, default = 0, dest = "min_minor_depth_fwd", help="the minimum number of minor alleles on the forward strand to call heteroplasmies")
	parser.add_argument("--minh-rev", type = int, default = 0, dest = "min_minor_depth_rev", help="the minimum number of minor alleles on the reverse strand to call heteroplasmies")
	parser.add_argument("--min-het", type = float, default = 0.005, dest = "min_het_freq", help="the minimum minor allele fraction of heteroplasmies")
	parser.add_argument("--min-qual", type = int, default = 20, dest = "min_qual", help="the minimum base quality")
	parser.add_argument("--max-qual", type = int, default = 10000, dest = "max_qual", help="the maximum base quality")
	parser.add_argument("--min-qual-rate", type = float, default = 0.0, dest = "min_qual_rate", help="the minimum proportion of bases passing the quality filter(s)")
	parser.add_argument("--mle", action = "store_true", default = False, dest = "mle", help="use maximum likelihood estimation to compute variant quality for heteroplasmies")
	parser.add_argument("--family", type = str, default = "a", dest = "family", help="family name")
	parser.add_argument("--sample", type = int, nargs="+", default = None, dest = "sample", help="sample columns to process in the mpileup file")
	parser.add_argument("--name", type = str, nargs="+", default = ["a",], dest = "name", help="the corresponding names of the samples to process")
	parser.add_argument("--batch", action="store_true", default = False, dest = "batch", help="proceed in the batch mode (read arguments from the batch file)")
	parser.add_argument("--all-sites", action="store_true", default = False, dest = "all_sites", help="output information for all sites instead of only variant sites")
	parser.add_argument("--summarize-coverage", action = "store_true", default = False, dest = "summarize_coverage", help="summarize the depth of read coverage at each mtDNA site")
	parser.add_argument("--harmonize-coverage", action = "store_true", default = False, dest = "harmonize_coverage", help="down sampling reads to achieve equal sequencing depth between samples")
	parser.add_argument("--pair-wise", action = "store_true", default = False, dest = "pair", help = "do pair-wise comparison of the heteroplasmic fractions between samples")
	parser.add_argument("--ref-sample", type = int, dest = "ref_sample", help="the reference sample used in the comparison")
	parser.add_argument("input", type=str, help = "a mpileup file or a batch file")
	parser.add_argument("output", type=str, help = "the prefix of output files")
	options = parser.parse_args(args)
	
	batch = []
	name = []
	if (options.batch):
		#batch mode
		batch_parser = ArgumentParser(usage="[--sample id] [--name id] batchfile")
		#batch_parser.add_argument("--family", type = str, default = None, dest = "family", help="specify family name")
		batch_parser.add_argument("--sample", type = int, nargs="+", default = None, dest = "sample", help="samples to process")
		batch_parser.add_argument("--name", type = str, nargs="+", default = ["a",], dest = "name", help="specify sample names")
		batch_parser.add_argument("input", type=str, help = "mpileup file")
		with open(options.input) as fh:
			for line in fh:
				line = line.rstrip("\r\n")
				if (not line):
					continue
				batch_options = batch_parser.parse_args(line.split())
				scanner = MTScan(batch_options.input, batch_options.sample, batch_options.name, options.min_qual, options.min_qual_rate, options.min_depth, options.min_depth_fwd, options.min_depth_rev, options.min_minor_depth, options.min_minor_depth_fwd, options.min_minor_depth_rev, options.min_het_freq)
				batch.append(scanner.allSites())
				name.extend(batch_options.name)
	else:
		scanner = MTScan(options.input, options.sample, options.name, options.min_qual, options.min_qual_rate, options.min_depth, options.min_depth_fwd, options.min_depth_rev, options.min_minor_depth, options.min_minor_depth_fwd, options.min_minor_depth_rev, options.min_het_freq)
		batch.append(scanner.allSites())
		name = options.name
	assert len(set(name)) == len(name), "Duplicated sample name in %s ." % ";".join(name)
	if (options.output == "-"):
		out = sys.stdout
	else:
		out = open(options.output + ".var", "w")
	if (options.summarize_coverage and options.output != "-"):
		out_coverage =  open(options.output + ".coverage", "w")
		out_coverage.write("family\tchr\tpos")
		for i in name:
			out_coverage.write("\t")
			out_coverage.write("\t".join([i+".depth",i+".depth_fwd", i+".depth_rev", i+".depth_hq"]))
		out_coverage.write("\n")
	else:
		out_coverage = None
	n = 1
	head = "\t".join(["family","sample","chr","pos","ref","depth","depth_fwd","depth_rev","allele","A","T","C","G","a","t","c","g",\
					  "heteroplasmy","substitution","het_allele","het_freq","het_freq_mle","het_freq_llr","het_low","het_high","het_p_fisher","het_p_sbias"])
	out.write(head + "\n")
	
	cur = 1
	end = 16570
	while (cur < end):
		sites_total = []
		is_var_total = False
		for scanner in batch:
			sites, pos, is_var = scanner.next()
			assert pos == cur, "Mismatched position %d, %d" % (pos, cur)
			sites_total.extend(sites)
			if (is_var):
				is_var_total = True
		if (out_coverage):
			out_coverage.write("%s\tchrM\t%s" % (options.family,cur))
			for i, site in enumerate(sites_total):
				out_coverage.write("\t")
				if (site is not None):
					out_coverage.write("\t".join(map(str, [site.depth, site.depth_fwd, site.depth_rev, (site.depth_fwd+site.depth_rev)/float(site.depth) if site.depth else 0.0])))
				else:
					out_coverage.write("\t".join(map(str, [0, 0, 0, 0])))
			out_coverage.write("\n")
		#skip invariant sites if all_sites has not been turned on
		if (options.all_sites or is_var_total):
			if (options.harmonize_coverage):
				min_depth = min([len(site.readsquality) if site is not None else 0 for site in sites_total])
				if (min_depth > 0 and min_depth >= options.min_depth):
					for i, site in enumerate(sites_total):
						#output site
						if (len(site.readsquality) > min_depth):
							site = site.downSample(min_depth)
							site.callAllele(options.min_depth, options.min_depth_fwd, options.min_depth_rev, options.min_minor_depth, options.min_minor_depth_fwd, options.min_minor_depth_rev, options.min_het_freq)
						print >>out, options.family + "\t" + name[i] + "\t" + site.string()
			else:
				for i, site in enumerate(sites_total):
					#output site
					if (site is not None):
						print >>out, options.family + "\t" + name[i] + "\t" + site.string()
					else:
						print >>out, options.family + "\t" + name[i] + "\t" + "\t".join(["chrM",str(cur)]+["",]*23)
		cur += 1
	
	if (out != sys.stdout):
		out.close()
		if (out_coverage):
			out_coverage.close()

def main():
	run(sys.argv[0], sys.argv[1:])
	
if __name__ == "__main__":
	main()