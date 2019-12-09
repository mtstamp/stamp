#############################################################################
#
#Updated by: Yiqin Wang
#
##############################################################################

import sys, os
import numpy

sites = {}
coverage = {}
mtdna_len = 16569
probe_set = ["A%d"% i for i in xrange(1,13)]
probe_set = probe_set + ["B%d"% i for i in xrange(1,13)]
probe_set = probe_set + ["C%d"% i for i in xrange(1,13)]
probe_set = probe_set + ["D%d"% i for i in xrange(1,11)]

def countif(n, func):
	"""
	return number of items in n having a non-zero return from func
	"""
	n = [1 for i in n if (func(i))]
	return sum(n)

def count_coverage(depth):
	"""
	return values of mean, median, 25% percentile, 75% percentile, count of over 100, count of 500, count of 1000, \
	count of over 0.5xmean, count of over 0.2xmean, count of over 0.1xmean
	"""
	cov_mean = sum(depth)/float(mtdna_len)
	cov_median = numpy.percentile(depth, 50)
	cov_25p = numpy.percentile(depth, 25)
	cov_75p = numpy.percentile(depth, 75)
	cov_gt1000 = countif(depth, lambda x: x >= 1000)
	cov_gt500 = countif(depth, lambda x: x >= 500)
	cov_gt100 = countif(depth, lambda x: x >= 100)
	cov_gt05m = countif(depth, lambda x: x >= 0.5*cov_mean)
	cov_gt02m = countif(depth, lambda x: x >= 0.2*cov_mean)
	cov_gt01m = countif(depth, lambda x: x >= 0.1*cov_mean)
	return [cov_mean, cov_median, cov_25p, cov_75p, cov_gt100, cov_gt500, cov_gt1000, cov_gt05m, cov_gt02m, cov_gt01m]

#out = open(out_file, "w")
if __name__ == "__main__":
	if (len(sys.argv) < 2):
		print >>sys.stderr, "Usage: %s file1.coverage [file2.coverage] ..." %sys.argv[0]
		exit(-1)
	#write summary to stdout
	out = sys.stdout
	head = ["name", "file"]
	#file header
	for mtdna_cov in ["cov_mean","cov_median","cov_25p","cov_75p","cov_gt100","cov_gt500","cov_gt1000","cov_gt05xmean","cov_gt02xmean","cov_gt01xmean"]:
		head.append(mtdna_cov)
		head.append(mtdna_cov + "_hq")
	for probe in probe_set:
		head.append("cov_probe_%s" % probe)
		head.append("cov_probe_%s_hq" % probe)
	out.write("\t".join(head) + "\n")
	#for fname in glob.glob(var_path + os.path.sep + "*.coverage"):
	#read coverage files from pileup one by one
	for fname in sys.argv[1:]:
		if (fname.endswith(".coverage")):
			name = fname[(fname.rfind(os.path.sep)+1):-len(".coverage")]
			print >>sys.stderr, "Processing %s: %s" % (name, fname)
			#intialize
			depth = [0,]*mtdna_len
			depth_hq = [0,]*mtdna_len
			depth_probe = {}
			depth_probe_hq = {}
			#open file
			with open(fname, "r") as fh:
				fh.readline() #skip head line
				for line in fh:
					line = line.rstrip("\r\n").split("\t")
					pos = int(line[2])-1
					"""
					the 4th to 10th columns in coverage file record
					depth, depth with BAQ between 0 and 10, between 10 and 20, between 20 and 30, between 30 and 40, over 40
					"""
					d, q0, q1, q2, q3, q4 = map(int, line[4:10])
					depth[pos] = d
					#high-quality read (BAQ > 30)
					depth_hq[pos] = q3+q4
					#probe information
					probe = line[-1]
					#skip sites covered by two probes
					if (probe.find("|") == -1):
						probe = probe.split("(")[0]
						if (probe not in depth_probe):
							depth_probe[probe] = []
							depth_probe_hq[probe] = []
						#store depths of all sites from this probe
						depth_probe[probe].append(d)
						depth_probe_hq[probe].append(q3+q4)
			#data to output
			sample_cov = [name,fname]
			#stats of all reads
			depth = count_coverage(depth)
			#stats of high-quality reads 
			depth_hq = count_coverage(depth_hq)
			#rearrange based on the head
			for i in range(len(depth)):
				sample_cov.append(depth[i])
				sample_cov.append(depth_hq[i])
			for probe in probe_set:
				cov = depth_probe.get(probe, [0,])
				#use the mean depth for each probe (may only diff a little at the ends of the reads)
				sample_cov.append(float(sum(cov))/len(cov))
				cov_hq = depth_probe_hq.get(probe, [0,])
				sample_cov.append(float(sum(cov_hq))/len(cov_hq))
			#output
			out.write("\t".join(map(str, sample_cov))+"\n")
	out.close()
