#!/usr/bin/python
#############################################################################
#Package: stamp
#Updated By Yiqin Wang 
#Usage:
#main
#retrieve variant information
#require .var file output from scan.py
#
##############################################################################

import sys

from argparse import ArgumentParser

def median(v):
	#median value in the v
	v = sorted(v)
	if (len(v) % 2):
		return v[len(v)/2]
	else:
		return float(v[len(v)/2-1] + v[len(v)/2]) / 2.0

def mean(v):
	#mean value in the v
	return sum(v)/float(len(v))

def countif(v, f):
	#count the number of items in v filling the requirement of f
	n = 0
	for i in v:
		if (f(i)):
			n += 1
	return n

def toInt(s, d = None):
	#int convertion
	try:
		s = int(s)
		return s
	except:
		return d

def toFloat(s, d = None):
	#float convertion
	try:
		s = float(s)
		return s
	except:
		return d

#globe variables
llr_min = 5
sbias_min = 0.01
reference_seq = {}

class MtVariant:
	""" 
	This is a class for storing information of one mtDNA variant
	
	Attributes
	----------
	family: family id
	sample: sample id
	pos: position
	depth: depth of reads
	ref: the reference allele
	depth_fwd: depth of reads on the forward strand
	depth_rev: depth of reads on the reverse strand
	depth_qc: depth of reads passing QC
	depth_ratio: the ratio of reads passing QC
	alt_freq: alternate allele fraction (minor allele fraction if not specified) 0 if this variant does not pass QC (low_qual=true)
	allele: the major allele
	alt_allele: the minor allele 
	alt_freq_llr: the log-likelihood estimation of the minor allele fraction
	alt_freq_sbias: the strand bias P value
	low_qual: true/false indicator of whether the variant pass all the QC filters
	alt_freq_raw: the original minor allele fraction regardless of the QC status
	
	"""
	
	def __init__(self, family, sample, pos, ref, depth, depth_fwd, depth_rev, allele, alt_allele, alt_freq, alt_freq_llr, alt_freq_sbias):
		""" 
		the __init__ method
		
		Arguments
		----------
		see class Attributes 
		"""
		
		self.family = family
		self.sample = sample
		self.pos = int(pos)
		self.depth = int(depth)
		if (self.pos not in reference_seq):
			reference_seq[self.pos] = ref
		self.ref = ref
		self.depth_fwd = int(depth_fwd)
		self.depth_rev = int(depth_rev)
		self.depth_qc = self.depth_fwd + self.depth_rev
		self.depth_ratio = float(self.depth_qc)/self.depth if (self.depth > 0) else 0.0
		self.alt_freq = toFloat(alt_freq, 0.0)
		self.allele = allele
		self.alt_allele = alt_allele 
		if (self.alt_freq > 0.5):
			#self.alt_freq = 1.0 - self.alt_freq
			#self.allele = alt_allele
			#self.alt_allele = allele
			pass
		self.alt_freq_llr = toFloat(alt_freq_llr, 10000.0)
		self.alt_freq_sbias = toFloat(alt_freq_sbias, 1.0)
		self.alt_freq_raw = self.alt_freq
		if (self.alt_freq_llr <= llr_min or self.alt_freq_sbias <= sbias_min):
			#set frequency to zero
			self.alt_freq = 0.0
			self.low_qual = True
			#self.alt_allele = ""
			#if (self.alt_freq_sbias <= sbias_min):
			#	 self.alt_freq_raw = 0 
			#	self.depth_qc = 0 
			#self.alt_freq_llr = 10000.0
			#self.alt_freq_sbias = 1.0
		else:
			self.low_qual = False
		
		self.idrev = self.id = self.ref + str(self.pos) + self.ref
		self.dev_freq = self.alt_freq
		#determine the derived allele
		if (self.allele and self.allele != self.ref):
			self.id = self.ref + str(self.pos) + self.allele
			self.idrev = self.allele + str(self.pos) + self.ref
			self.dev_freq = 1.0-self.alt_freq
		elif (self.alt_allele and self.alt_allele != self.ref):
			self.id = self.ref + str(self.pos) + self.alt_allele
			self.idrev = self.alt_allele + str(self.pos) + self.ref

def getMtVariant(data, depth_min = 10, depth_ratio_min = 0.0):
	""" 
	extract variants in the data that can pass QC filters specified 
	
	Attributes
	----------
	data: a dict of variants grouped by family and sample
	depth_min: the minimum depth of a variant
	depth_ratio_min: the minimum ratio of reads that pass QC
	
	Returns
	----------
	a dict of extracted variants grouped by family and sample
	
	"""
	var = {}
	n = 0
	for family in data:
		var[family] = {}
		for sample in data[family]:
			var[family][sample] = {}
			for variant in data[family][sample].values():
				if (variant.ref == "N"):
					continue
				if (variant.depth_qc >= depth_min and variant.depth_ratio >= depth_ratio_min):
					var[family][sample][variant.pos] = variant
					n += 1
				else:
					var[family][sample][variant.pos] = None
	print "Read %d mitochondrial DNA variants" % n
	return var

def getMtMajorAllele(data, depth_min = 10, depth_ratio_min = 0.0):
	""" 
	extract variants that have the major allele different from the reference allele 
	
	Attributes
	----------
	data: a dict of variants grouped by family and sample
	depth_min: the minimum depth of a variant
	depth_ratio_min: the minimum ratio of reads that pass QC
	
	Returns
	----------
	a dict of extracted variants grouped by family and sample
	
	"""
	var = {}
	n = 0
	for family in data:
		var[family] = {}
		for sample in data[family]:
			var[family][sample] = {}
			for variant in data[family][sample].values():
				if (variant.ref == "N"):
					continue
				if (variant.depth_qc >= depth_min and variant.depth_ratio >= depth_ratio_min):
					if (variant.allele != variant.ref):
						var[family][sample][variant.pos] = variant
						n += 1
	print "Read %d mitochondrial DNA variants with major allele different than the reference allele" % n
	return var
	
def getMtHomoplasmy(data, depth_min = 10, depth_ratio_min = 0.0, freq_min = 0.01):
	""" 
	extract homoplasmies
	
	Attributes
	----------
	data: a dict of variants grouped by family and sample
	depth_min: the minimum depth of a variant
	depth_ratio_min: the minimum ratio of reads that pass QC
	freq_min: the maximum frequency of a possible minor allele
	 
	Returns
	----------
	a dict of extracted variants grouped by family and sample
	
	"""
	homo = {}
	n = m = 0
	for family in data:
		homo[family] = {}
		for sample in data[family]:
			homo[family][sample] = {}
			for variant in data[family][sample].values():
				if (variant.ref == "N"):
					continue
				n += 1
				if (variant.depth_qc >= depth_min and variant.depth_ratio >= depth_ratio_min):
					if (variant.allele != variant.ref and variant.alt_freq < freq_min):
						homo[family][sample][variant.pos] = variant
						m += 1
				else:
					#include missing
					homo[family][sample][variant.pos] = None
	print "Read %d mitochondrial homoplasmies" % m
	return homo

def isHeteroplasmy(variant, depth_min = 40, depth_strand = 0, depth_ratio_min = 0.0, freq_min = 0.01):
	""" 
	determine whether a variant is a heteroplasmy according to the filers specified
	
	Attributes
	----------
	variant: MTVariant
	depth_min: the minimum depth
	depth_strand: the minimum depth on either the forward and the reverse strand
	depth_ratio_min: the minimum ratio of reads passing QC
	freq_min: the minimum minor allele frequency
	
	Returns
	----------
	True: is a heteroplasmy
	False: not a heteroplasmy
	None: does not pass QC
	
	"""
	if (variant.depth_qc >= depth_min and variant.depth_fwd >= depth_strand and variant.depth_rev >= depth_strand and variant.depth_ratio >= depth_ratio_min):
		if (variant.alt_freq >= freq_min):
			return True
		else:
			return False
	else:
		return None

def getMtHeteroplasmy(data, depth_min = 40, depth_strand = 0, depth_ratio_min = 0.0, freq_min = 0.01, freq_min_raw = 0):
	""" 
	extract heteroplasmies
	
	Attributes
	----------
	data: a dict of variants grouped by family and sample
	depth_min: the minimum depth
	depth_strand: the minimum depth on either the forward and the reverse strand
	depth_ratio_min: the minimum ratio of reads passing QC
	freq_min: the minimum minor allele frequency (after QC)
	freq_min_raw: the minimum minor allele frequency (before QC; see MTVariant for details)
	
	Returns:
	----------
	a dict of extracted variants grouped by family and sample
	
	"""
	het = {}
	n = 0
	for family in data:
		het[family] = {}
		for sample in data[family]:
			het[family][sample] = {}
			for variant in data[family][sample].values():
				if (variant.ref == "N"):
					continue
				if (variant.depth_qc >= depth_min and variant.depth_fwd >= depth_strand and variant.depth_rev >= depth_strand and variant.depth_ratio >= depth_ratio_min):
					if (variant.alt_freq >= freq_min and variant.alt_freq_raw >= freq_min_raw):
						het[family][sample][variant.pos] = variant
						n += 1
	print "Read %d mitochondrial heteroplasmies" % n
	return het   

def readMtVariant(variant_file, fam_excl = {}, pos_excl = {}):
	""" 
	read variants stored in the variant file
	
	Attributes
	----------
	variant_file: the path to the variant file
	fam_excl: family ids to be excluded (dict, list or tuple)
	pos_excl: positions to be excluded (dict, list or tuple)
	
	Returns
	----------
	a tuple of
	1. file head information
	2. a dict of variants grouped by family and sample
	
	"""
	data = {}
	n = 0
	if (variant_file == "-"):
		#use standard input instead
		fh = sys.stdin
	else:
		fh = open(variant_file)
	head = fh.readline()
	head = head.rstrip("\r\n").split("\t")
	assert len(head) >= 27, "Truncated head line of the variant file"
	for line in fh:
		line = line.rstrip("\r\n").split("\t")
		family,sample,chr,pos,ref,depth,depth_fwd,depth_rev,allele,A,T,C,G,a,t,c,g,\
		heteroplasmy,substitution,het_allele,het_freq,het_freq_mle,het_freq_llr,het_low,het_high,het_p_fisher,het_p_sbias = line[:27]
		if (family in fam_excl):
			#exclude this family
			continue
		if (family == "family"):
			#skip the head line
			continue
		#a new variant
		variant = MtVariant(family,sample,pos,ref,depth,depth_fwd,depth_rev,allele,het_allele,het_freq_mle,het_freq_llr,het_p_sbias)
		#temporarily store the original line
		variant.line_cache = line[:]
		if (family not in data):
			data[family] = {}
		if (sample not in data[family]):
			data[family][sample] = {}
		pos = variant.pos
		assert pos not in data[family][sample], "Duplicated vairant at position %d in sample %s." % (variant.pos, variant.sample)
		if (pos in pos_excl):
			#exclude this position
			continue
		data[family][sample][pos] = variant
		n += 1
	print "Read %d mitochondrial DNA variants" % n
	return head, data

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
	parser = ArgumentParser(prog = prog, usage="%(prog)s [-options] [input] [output]\n", version="%(prog)s v0.1.1")
	parser.add_argument("--exclude", type = str, dest="exclude", help="exclude mtDNA sites from analysis")
	parser.add_argument("--remove", type = str, dest="remove", help="remove families from analysis")
	parser.add_argument("--keep", type = str, dest="keep", help="keep only the families for analysis")
	parser.add_argument("--depth", type = float, default = 10, dest="depth", help="the minimum read depth of all variants")
	parser.add_argument("--depth-min", type = float, default = 40, dest="depth_min", help="the minimum read depth of heteroplasmies")
	parser.add_argument("--hq-min", type = float, default = 0.7, dest="hq_min", help="the minimum ratio of high-quality reads of heteroplasmies")
	parser.add_argument("--llr-min", type = float, default = 5, dest="llr_min", help="the minimum quality score of heteroplasmies")
	parser.add_argument("--sbias-min", type = float, default = 0.001, dest="sbias_min", help="the minimum P value for strand bias analysis of heteroplasmies")
	parser.add_argument("--frac-min", type = float, default = 0.01, dest="frac_min", help="the minimum minor allele fraction of heteroplasmies")
	parser.add_argument("--dev-frac-min", type = float, default = 0.90, dest="dev_frac_min", help="the minimum variant allele fraction of homoplasmies")
	parser.add_argument("--annotate", type = str, dest="annotate", help="annotate variants according to the file specified")
	parser.add_argument("--output-ped", default=False, action="store_true", dest="output_ped", help="output the variants detected to a ped file")
	parser.add_argument("--output-hsd", default=False, action="store_true", dest="output_hsd", help="output major allele to the hsd file")
	parser.add_argument("--output-minor-hsd", default=False, action="store_true", dest="output_minor_hsd", help="output minor allele to the hsd file")
	parser.add_argument("input", help="the variant file output from scan")
	parser.add_argument("output", help="the prefix of output files")
	options = parser.parse_args(args)
	
	#initialize globle variables
	global llr_mim, sbias_min
	llr_min = options.llr_min
	sbias_min = options.sbias_min
	
	pos_excl = {}
	if (options.exclude):
		with open(options.exclude) as fh:
			for line in fh:
				line = line.strip()
				try:
					pos_excl[int(line)] = 1
				except:
					continue
	fam_excl = {}
	if (options.remove):
		with open(options.remove) as fh:
			for line in fh:
				line = line.strip()
				try:
					fam_excl[int(line)] = 1
				except:
					continue
	head, data = readMtVariant(options.input, fam_excl, pos_excl)
	#head = "family\tsample\tchr\tpos\tref\tdepth\tdepth_fwd\tdepth_rev\tallele\tA\tT\tC\tG\ta\tt\tc\tg\theteroplasmy\tsubstitution\thet_allele\thet_freq\thet_freq_mle\thet_freq_llr\thet_low\thet_high\thet_p_fisher\thet_p_sbias".split("\t")
	head.append("stat")
	
	annot = {}
	annot_len = 0
	if (options.annotate):
		#read annotation file
		#build annotation table
		if (options.annotate.endswith(".csv")):
			delim = ","
		else:
			delim = "\t"
		with open(options.annotate, "r") as fh:
			line = fh.readline()
			line = line.rstrip("\r\n").split(delim)
			n = line.index("id")
			head.extend(line[n+1:])
			annot_len = len(line) - n - 1
			for line in fh:
				line = line.rstrip("\r\n").split(delim)
				if (not line):
					continue
				id = line[n]
				annot[id] = line[n+1:]
	annot_null = ["",]*annot_len
	
	if (options.output_hsd):
		out_hsd = open(options.output + ".qc.hsd", "wb")
		out_hsd.write("SampleId\tRange\tHaplogroup\tPolymorphisms (delimited with tabs)\n")
		if (options.output_minor_hsd):
			out_minor_hsd = open(options.output + "minor.qc.hsd", "wb")
			out_minor_hsd.write("SampleId\tRange\tHaplogroup\tPolymorphisms (delimited with tabs)\n")
		else:
			out_minor_hsd = None
	else:
		out_hsd = None
		out_minor_hsd = None
	
	var_all = getMtVariant(data, depth_min = 0, depth_ratio_min = 0)
	sample_all = []
	if (options.keep):
		with open(options.keep, "rb") as fh:
			for line in fh:
				line = line.rstrip("\r\n")
				if (not line):
					continue
				family, sample = line.split("\t")
				if (family not in fam_excl):
					sample_all.append([family, sample])
	else:
		for family in sorted(var_all.keys()):
			for sample in sorted(var_all[family].keys()):
				sample_all.append([family, sample])
	
	#output sample names
	#order corresponds to that of the samples in the ped file and the hsd file
	with open(options.output + ".qc.tfam", "wb") as out_fam:
		for family, sample in sample_all:
			#use the default phenotype value -9
			out_fam.write("\t".join([family, sample, "0", "0", "-9", "-9"])+"\n")
	
	sites_all = {}
	with open(options.output + ".qc.annot", "wb") as out:
		#output the head line
		out.write("\t".join(head) + "\n")
		idx = 0 #sample idx
		for family, sample in sample_all:
			if (family in var_all and sample in var_all[family]):
				var = var_all[family][sample]
			else:
				var = {}
			homoplasmy = []
			heteroplasmy = []
			for pos in sorted(var.keys()):
				v = var[pos]
				if (v):
					if (isHeteroplasmy(v, depth_min = options.depth_min, depth_strand = 0, depth_ratio_min = options.hq_min, freq_min = options.frac_min)):
						stat = "heteroplasmy"
						heteroplasmy.append(v)
						add_var = True
						a1 = v.allele
						a2 = v.alt_allele
					elif (v.depth >= options.depth and v.allele !=  v.ref and v.dev_freq >= options.dev_frac_min):
						stat = "homoplasmy"
						homoplasmy.append(v)
						add_var = True
						a1 = a2 = v.allele
					elif (v.alt_freq_raw > options.frac_min):
						#variant does not pass the filters of variant quality and strand bias (see MTVariant)
						stat = "heteroplasmy possible"
						add_var = False
						a1 = a2 = None
					else:
						stat = "unkown"
						add_var = False
						a1 = a2 = None
					out.write("\t".join(v.line_cache + [stat,] + annot.get(v.id,annot_null))+"\n")
				else:
					stat = "unkown"
					add_var = True
					a1 = a2 = "N"
				if (add_var):
					if (pos not in sites_all):
						sites_all[pos] = [0, 0, 0, {}] ##homoplamy, #heteroplasmy, #missing, #{sample: allele}
					site = sites_all[pos]
					site[3][idx]= a1+"\t"+a2
					if (a1 != a2):
						site[1] += 1
					elif (a2 == "N"):
						site[2] += 1
					else:
						site[0] += 1
			idx += 1
			if (out_hsd):
				#use sample index (one-based) instead of the real sample name
				#output major alleles
				major_allele = [[v.pos,str(v.allele)] for v in homoplasmy] + [[v.pos,str(v.allele)] for v in heteroplasmy if v.dev_freq >= options.dev_frac_min]
				if (not major_allele):
					major_allele = ["1G",]
				else:
					major_allele.sort()
					major_allele = [str(p)+str(a) for p, a in major_allele]
				out_hsd.write("\t".join([str(idx),"1-16569;","?"] + major_allele)+"\n")
				if (out_minor_hsd):
					minor_allele = [[v.pos,str(v.alt_allele)] for v in heteroplasmy]
					if (minor_allele):
						minor_allele.sort()
						minor_allele = [str(p)+str(a) for p, a in minor_allele]
					out_minor_hsd.write("\t".join([str(idx),"1-16569;","?"] + minor_allele)+"\n")

	if (out_hsd):
		out_hsd.close()
		if (out_minor_hsd):
			out_minor_hsd.close()
	
	if (options.output_ped):
		sites = sorted(sites_all.keys())
		out_ped = open(options.output + ".qc.tped", "wb")
		out_map = open(options.output + ".qc.map", "wb")
		for i in sites:
			site = sites_all[i]
			out_ped.write("\t".join(["26",str(i)+reference_seq[i],"0",str(i)]))
			out_ped.write("\t")
			site_sample = site[3]
			ref = reference_seq[i]+"\t"+reference_seq[i]
			out_ped.write("\t".join([site_sample.get(j,ref) for j in range(len(sample_all))]))
			out_ped.write("\n")
			out_map.write("\t".join(["26",str(i)+reference_seq[i],"0",str(i)]+list(map(str,site[:3])))+"\n")
	"""
	as plink does not handle multi-allelic variants
	to remove these variants in R
	nallele <- apply(tped[,5:ncol(tped)],1,function(x){length(unique(x))})
	write.table(tped[nallele==2,], "biallele.tped", sep="\t", quote = F,col.names = F, row.names = F)
	"""

def main():
	run(sys.argv[0], sys.argv[1:])
	
if __name__ == "__main__":
	main()
