#!/usr/bin/python
#############################################################################
#Package: stamp
#Updated by Yiqin Wang
#Usage:
#summarize read coverage for whole-genome sequencing data
#
##############################################################################

import os, sys
import time
from argparse import ArgumentParser

from FileIO import pipe_output, read_arguments
#from MTFactory import MTFactory, MTThreading

samtools = "samtools-1.3"

#length of each chromosome
chr_hg38 = {'chr1':248956422 ,
'chr2':242193529 ,
'chr3':198295559 ,
'chr4':190214555 ,
'chr5':181538259 ,
'chr6':170805979 ,
'chr7':159345973 ,
'chr8':145138636 ,
'chr9':138394717 ,
'chr10':133797422 ,
'chr11':135086622 ,
'chr12':133275309 ,
'chr13':114364328 ,
'chr14':107043718 ,
'chr15':101991189 ,
'chr16':90338345 ,
'chr17':83257441 ,
'chr18':80373285 ,
'chr19':58617616 ,
'chr20':64444167 ,
'chr21':46709983 ,
'chr22':50818468 ,
'chrX':156040895 ,
'chrY':57227415 ,
'chrM':16569
}

def median(v):
	#median value in the v
	if (not v):
		return None 
	v = sorted(v)
	if (len(v) % 2):
		return v[len(v)/2]
	else:
		return float(v[len(v)/2-1] + v[len(v)/2]) / 2.0

def mean(v):
	#mean value in the v
	if (not v):
		return None 
	return sum(v)/len(v)

def timeDiff(start, end):
	#compute running time
	time_elapsed = end - start
	hours = int(time_elapsed/60/60)
	minutes = int((time_elapsed - hours*60*60)/60)
	seconds =  time_elapsed - hours*60*60 - minutes*60
	return "%d hours, %d minutes, and %f seconds" % (hours, minutes, seconds) 
	
def calculateCoverage(alignment_file, ref_seq, chr, qual_filter, proper_pair, win_size = 100000, sliding_size = 50000):
	"""
	Compute depth of read coverage

	Arguments
	----------
	alignment file: input alignment file
	ref_seq: the reference genome sequence used for read alignment
	chr: the chromosome to be processed
	qual_filter: quality filter string
	proper_pair: retain only unique, properly paired reads
	win_size: compute read coverage across the reference genome in sliding windows
	sliding_size: sliding size

	Returns
	----------
	a tuple of
	1. a list of read coverage in each window
	2. a lit of sites covered with reads in each window
	"""

	assert chr in chr_hg38, "Cannot find chromosome %s" % chr
	#depth file handle
	if (proper_pair):
		ph = pipe_output("%s view -h -T %s -F 0xF08 -f 0x2 %s %s|%s depth %s -a -" % (samtools, ref_seq, alignment_file, chr, samtools, qual_filter))
	else:
		ph = pipe_output("%s depth %s --reference %s -r %s -a %s" % (samtools, qual_filter, ref_seq, chr, alignment_file))

	#assume the window size is an integer times the sliding size
	sliding_ratio = win_size/sliding_size
	win_size = sliding_ratio*sliding_size

	if (chr != "chrM"):
		#chr 1-Y
		chr_length = chr_hg38[chr]
		blocks = (chr_length/sliding_size)+1 
		win_coverage = [0,]*blocks
		win_sites = [0,]*blocks
		#chr_coverage = [0,]*blocks
		#chr_sites = [0,]*blocks
		pos = 0
		win_end = sliding_size
		win_cur = 0
		total_cov = 0
		total_sites = 0
		#count read numbers in windows of the sliding size so that depth at each site is counted once.
		for line in ph.stdout:
			#line = line.strip().split("\t")
			#name = line[0]
			#pos = int(line[1])-1
			#cov = int(line[2])
			cov = int(line[line.rfind("\t"):])
			total_cov += cov
			if (cov > 0):
				total_sites += 1 
			#win_end = pos/sliding_size+1
			#win_start = max(0,(pos-win_size)/sliding_size+1)
			#for i in xrange(max(0,win_start), win_end):
			#	chr_coverage[i] += cov
			#	chr_sites[i] += 1
			pos += 1
			if (pos == win_end):
				win_coverage[win_cur] = total_cov
				win_sites[win_cur] = total_sites
				win_end += sliding_size
				win_cur += 1
				total_cov = 0
				total_sites = 0
		if (win_cur < len(win_coverage)):
			win_coverage[win_cur] = total_cov
			win_sites[win_cur] = total_sites
		chr_coverage = []
		chr_sites = []
		#compute read depth in windows of the window size based on the "sliding_ratio" windows
		for win_cur in xrange(len(win_coverage)-sliding_ratio+1):
			total_coverage = 0
			total_sites = 0
			for sliding_cur in  range(sliding_ratio):
				total_coverage += win_coverage[win_cur+sliding_cur]
				total_sites += win_sites[win_cur+sliding_cur]
			chr_coverage.append(total_coverage)
			chr_sites.append(total_sites)
	else:
		#chrM
		chr_length = chr_hg38[chr]
		#all sites, coding sites outside the D-loop region
		#record coverage for all sites and sites in the coding region
		#for mtDNA content calculation, use mtDNA coverage from the .coverage file output by wgs-align
		chr_coverage = [0,0]
		chr_sites = [0,0]
		pos = 0
		for line in ph.stdout:
			#line = line.strip().split("\t")
			#name = line[0]
			#pos = int(line[1])-1
			#cov = int(line[2])
			cov = int(line[line.rfind("\t"):])
			chr_coverage[0] += cov
			chr_sites[0] += 1
			if (pos >= 560 and pos <= 16000):
				chr_coverage[1] += cov
				chr_sites[1] += 1
			pos += 1
	return chr_coverage, chr_sites
	
def calculateMitoCN(coverage, output, win_size = 500000, sliding_size = 250000, min_sites = 0.8):
	"""
	output read depth for each chromosome

	Arguments
	----------
	coverage: a dict of coverage information with keys of chromosomes and values of results from calculateCoverage
	output: output path and file prefix
	win_size: compute read coverage across the reference genome in sliding windows
	sliding_size: sliding size
	min_sites: the minimum ratio of sites in a window covered with reads

	Returns
	----------
	None
	"""

	coverage_avg = {}
	chrM_cov_all = 0
	chrM_cov_coding = 0

	#output coverage information of all windows to the file
	with open(output + ".bin.coverage", "wb") as out:
		for chr in sorted(chr_hg38.keys()):
			coverage_bin = []
			if (chr not in coverage):
				continue
			chr_cov, chr_sites = coverage[chr]
			for n in xrange(len(chr_cov)):
				total_cov = chr_cov[n]
				total_sites = chr_sites[n]
				if (total_sites):
					avg_cov = total_cov/float(total_sites)
				else:
					avg_cov = 0.0
				out.write("%s\t%d\t%d\t%d\t%f\n" % (chr, n, total_cov, total_sites, avg_cov))
				if (chr != "chrM"):
					if (n == len(chr_cov)-1):
						#the last window
						if (total_sites >= int((chr_hg38[chr] % sliding_size)*min_sites)):
							coverage_bin.append(avg_cov)
					else:
						#retain only windows with >win_size*min_sites sites covered with reads
						if (total_sites >= int(win_size*min_sites)):
							coverage_bin.append(avg_cov)
				else:
					coverage_bin.append(avg_cov)
			if (chr != "chrM"):
				coverage_avg[chr] = (len(coverage_bin), mean(coverage_bin), median(coverage_bin))
			else:
				chrM_cov_all = coverage_bin[0]
				chrM_cov_coding = coverage_bin[1]

	#output coverage information for each chromosome to the file
	with open(output + ".bin.coverage.info", "wb") as out:
		for chr in sorted(chr_hg38.keys()):
			if (chr not in coverage_avg):
				continue
			chr_cov_num, chr_cov_avg, chr_cov_median = coverage_avg[chr]
			#each line:chromosome, window_number, average window coverage, median window coverage
			out.write("%s\t%d\t%f\t%f" % (chr, chr_cov_num, chr_cov_avg, chr_cov_median))
			#obsolete: compute mtDNA copy number as 2 times the read coverage in mtDNA and that in each chromosome 1-Y
			#for mtDNA content calculation, use mtDNA coverage from the .coverage file output by wgs-align instead
			if (chr_cov_avg):
				out.write("\t%f\t%f" % (chrM_cov_all/chr_cov_avg*2, chrM_cov_coding/chr_cov_avg*2))
			else:
				out.write("\tNA\tNA")
			if (chr_cov_median):
				out.write("\t%f\t%f\n" % (chrM_cov_all/chr_cov_median*2, chrM_cov_coding/chr_cov_median*2))
			else:
				out.write("\tNA\tNA\n")

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

	parser = ArgumentParser(prog=prog, usage="%(prog)s [-options] [outprefix]\n", version="%(prog)s v0.1")
	#parser.add_argument("--thread", "-t", type=int, default=1, dest="thread", help="specify number of parallele processes")
	parser.add_argument("--refseq", "-r", dest="refseq", help="specify the reference genome sequences")
	parser.add_argument("--win-size", type=int, default = 100000, dest="win_size", help="specify the size of the sliding window (default:100kb)")
	parser.add_argument("--sliding-size", type=int, default = 50000, dest="sliding_size", help="specify the step size of the sliding window (default:50kb)")
	parser.add_argument("--min-sites", type=float, default = 0.8, dest="min_sites", help="specify the minimum proportion of sites with read coverage allowed in a sliding window")
	parser.add_argument("--qual-filter", type=str, default = "-q 23 -Q 20", dest="qual_filter", help="specify quality filters")
	parser.add_argument("--proper-pair", action="store_true", dest="proper_pair", help="only include unique reads that are properly paired")
	parser.add_argument("--chr-incl", type=str, dest="chr_incl", help="specify chromosomes (default is all)")
	parser.add_argument("input", type=str, help = "input bam/cram file")
	parser.add_argument("output", type=str, help = "output file")
	options = parser.parse_args(args)

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

	#initialize parameters
	#enable multi-threading
	#njobs = options.thread
	#if (njobs < 1):
	#	njobs = 1
	input = options.input
	output = options.output

	#multi-thread handle
	#engine = MTThreading()
	#factory = MTFactory(engine, thread = options.thread, time=1)
	jobs = []
	
	if (options.chr_incl):
		chr_incl = {}
		for i in options.chr_incl.split(","):
			i = i.strip()
			if (not i):
				continue
			chr_incl[i] = 1
	else:
		chr_incl = None
	
	start = time.time()
	coverage = {}
	for chr in chr_hg38.keys():
		if (chr_incl is None or chr in chr_incl):
			#if (njobs > 1):
			#	#allocate one thread for each chromosome
			#	engine.setOption(func = calculateCoverage, args = (input, options.refseq, chr, options.qual_filter, options.proper_pair, options.win_size, options.sliding_size))
			#	jobs.append(factory.analyse(chr))
			#else:
				coverage[chr] = calculateCoverage(input, options.refseq, chr, options.qual_filter, options.proper_pair, options.win_size, options.sliding_size)

	#if (jobs):
	#	for chr in jobs:
	#		#wait until all calculation completes
	#		factory.waitResult(chr)
	#	for chr in jobs:
	#		#get return from calculateCoverage
	#		coverage[chr] = factory.getResult(chr)

	calculateMitoCN(coverage, output, options.win_size, options.sliding_size, options.min_sites)
	
	print >>sys.stderr, "Completed in %s!" % timeDiff(start, time.time())

def main():
	run(sys.argv[0], sys.argv[1:])

if __name__ == "__main__":
	main()

