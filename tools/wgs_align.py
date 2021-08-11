#!/usr/bin/python
#############################################################################
#Package: stamp
#Updated by Yiqin Wang
#Usage:
#generate mtDNA alignments for whole-genome sequencing data
#
##############################################################################

from argparse import ArgumentParser
import os, sys
import pysam

from SeqIO import phred
from FileIO import FileIO, execute, pipe_output, read_arguments

#executables
samtools = "samtools-1.6" 
bwa ="bwa-0.7.17"
bamleftalign="bamleftalign-1.1.0"
picard = "java -jar  ~/bin/picard-2.16.jar"

override = False

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
	return sum(v)/float(len(v))

def quantile(v, cutoff):
	#percentile in the v is larger than or equal to the cutoff
	return sum([i >= cutoff for i in v])/float(len(v))

def retrieve_alignment(name, infile, outfile, fastq = None, remap_genome = True, genome_refseq = "", remap_mtdna = True, mtdna_refseq ="", mtdna_offset = 0, recal = True):
	"""
	perform genome alignment and quality filtering

	Arguments
	----------
	name: sample name
	infile: input path and file name prefix
	outfile: output path and file name prefix
	fastq: files of read 1 and read 2 (tuple (r1, r2))
	remap_genome: perform whole-genome alignment (true or false; default is true)
	genome_refseq: the reference genome sequences
	remap_mtdna: perform mitochondrial genome alignment (true or false; default is true)
	mtdna_refseq: the reference mitochondrial genome sequence
	mtdna_offset: the position offset used in parsing mtDNA read alignments
	recal: perform base quality recalibration (true or false; default is true)

	Returns
	----------
	the name of the alignment file after quality filtering (string)
	"""

	global override
	if (not override):
		if (recal):
			if (remap_mtdna):
				fname = "%s.mtdna.realign.recal.dedup.sorted.qc.bam" % outfile
			else:
				fname = "%s.realign.recal.dedup.sorted.qc.bam" % outfile
		else:
			fname = "%s.qc.bam" % outfile
		if (os.path.exists(fname)):
			#return file name if the output file has already been created.
			return fname

	#defaults of quality filters
	#maximum rate of nucleotide mismatches allowed
	error_rate = 0.05
	#minimum length of aligned read
	min_query_len = 40
	#minimum mapping quality score
	mapq_min = 20

	if (recal):
		if (remap_genome):
			#perform whole-genome alignment
			if (fastq):
				#mapping using fastq files
				execute("%s mem -M %s %s %s | %s view -hb - >%s.realign.bam" % (bwa, genome_refseq, fastq[0], fastq[1], samtools, outfile))
			elif (remap_genome):
				#convert bam to fastq
				execute("%s SamToFastq VALIDATION_STRINGENCY=LENIENT I=%s.bam F=%s.1.fastq F2=%s.2.fastq" % (picard, infile, outfile, outfile))
				#remap to the reference whole genome sequences
				execute("%s mem -M %s %s.1.fastq %s.2.fastq | %s view -hb - >%s.realign.bam" % (bwa, genome_refseq, outfile, outfile, samtools, outfile))
				#remove temporary fastq files
				execute("rm -rf %s.1.fastq %s.2.fastq" % (outfile, outfile))
			infile = outfile + ".realign"
			#sort and index
			execute("%s sort -o %s.realign.sorted.bam %s.bam" % (samtools, outfile, infile))
			execute("%s index %s.realign.sorted.bam" % (samtools, outfile))
			#mark duplicates
			execute("%s MarkDuplicates I=%s.realign.sorted.bam  O=%s.realign.dedup.sorted.bam METRICS_FILE=%s.realign.dedup.sorted.metrics" % (picard, outfile, outfile, outfile))
			#delete temporary files
			execute("rm -rf %s.realign.bam" % outfile)
			execute("rm -rf %s.realign.sorted.bam" % outfile)
			execute("rm -rf %s.realign.sorted.bam.bai" % outfile)
			infile = outfile + ".realign.dedup.sorted"
		
		if (remap_mtdna):
			#perform mitochondrial genome alignment
			#input alignment file
			pf = pysam.AlignmentFile(infile + ".bam", "rb")
			#output alignment file for reads after quality filtering
			out = pysam.AlignmentFile(outfile +".mtdna.bam", "wb", template = pf)
			#retain only read pairs mapped to chrM
			#record name of read 1 and read 2 and mapping quality
			r1 = {}
			r2 = {}
			mapq = {}
			for i in pf:
				if (i.is_unmapped):
					i.mapping_quality = 0
				if (i.reference_name == "chrM" and i.next_reference_name == "chrM" and i.is_paired): 
					#mate mapped to chrMT
					if (i.is_read1):
						r1[i.qname] = 1
					elif (i.is_read2):
						r2[i.qname] = 1
					if (i.qname not in mapq):
						mapq[i.qname] = []
					mapq[i.qname].append(i.mapq)
			pf.reset()
			for i in pf:
				if (i.qname in r1 and i.qname in r2 and max(mapq[i.qname]) >= mapq_min):
					#output read pairs with MAPQ>=mapq_min
					out.write(i)
			out.close()
			#remove whole genome alignment
			#if (not remap_genome):
			#	execute("rm -rf %s.realign.dedup.sorted.bam" % outfile)
			#convert bam to fastq
			execute("%s SamToFastq VALIDATION_STRINGENCY=LENIENT I=%s.mtdna.bam F=%s.mtdna.1.fastq F2=%s.mtdna.2.fastq" % (picard, outfile, outfile, outfile))
			#remap to the reference mitochondrial genome sequence
			execute("%s mem -M %s %s.mtdna.1.fastq %s.mtdna.2.fastq | %s view -hb - >%s.mtdna.realign.bam" % (bwa, mtdna_refseq, outfile, outfile, samtools, outfile))
			#remove temporary fastq files
			execute("rm -rf %s.mtdna.bam %s.mtdna.1.fastq %s.mtdna.2.fastq" % (outfile, outfile, outfile))
			#sort and index
			execute("%s sort -o %s.mtdna.realign.sorted.bam %s.mtdna.realign.bam" % (samtools, outfile, outfile))
			execute("%s index %s.mtdna.realign.sorted.bam" % (samtools,outfile))
			#mark duplicates
			execute("%s MarkDuplicates I=%s.mtdna.realign.sorted.bam  O=%s.mtdna.realign.dedup.sorted.bam METRICS_FILE=%s.mtdna.realign.dedup.sorted.metrics" % (picard, outfile, outfile, outfile))
			#remove temporary files
			execute("rm -rf %s.mtdna.realign.bam" % outfile)
			execute("rm -rf %s.mtdna.realign.sorted.bam" % outfile)
			execute("rm -rf %s.mtdna.realign.sorted.bam.bai" % outfile)

			#perform local realignment and base quality recalibration
			execute("rm -rf %s.mtdna.realign.recal.dedup.sorted.bam" % outfile)
			execute("%s -f %s < %s.mtdna.realign.dedup.sorted.bam | %s calmd -EArb - %s >>%s.mtdna.realign.recal.dedup.sorted.bam" % (bamleftalign, mtdna_refseq, outfile, samtools, mtdna_refseq, outfile))
			execute("rm -rf %s.mtdna.realign.dedup.sorted.bam" % outfile)

			#input and output alignment files for quality filtering
			pf_rm = "%s.mtdna.realign.recal.dedup.sorted.bam"%outfile
			pf = pysam.AlignmentFile(pf_rm, "rb")
			out_qc = pysam.AlignmentFile("%s.mtdna.realign.recal.dedup.sorted.qc.bam"%outfile, "wb", template = pf)
			out_rm = pysam.AlignmentFile("%s.mtdna.realign.recal.dedup.sorted.rm.bam"%outfile, "wb", template = pf) 
		else:
			#perform local realignment and base quality recalibration
			execute("rm -rf %s.realign.recal.dedup.sorted.bam" % outfile)
			execute("%s -f %s < %s.bam | %s calmd -EArb - %s >>%s.realign.recal.dedup.sorted.bam" % (bamleftalign, genome_refseq, infile, samtools, genome_refseq, outfile))

			#input and output alignment files for quality filtering
			pf_rm = "%s.realign.recal.dedup.sorted.bam"%outfile
			pf = pysam.AlignmentFile(pf_rm, "rb")
			out_qc = pysam.AlignmentFile("%s.realign.recal.dedup.sorted.qc.bam"%outfile, "wb", template = pf)
			out_rm = pysam.AlignmentFile("%s.realign.recal.dedup.sorted.rm.bam"%outfile, "wb", template = pf)
	else:
		#input and output alignment files for quality filtering
		pf_rm = None
		pf = pysam.AlignmentFile(infile + ".bam", "rb")
		out_qc = pysam.AlignmentFile(outfile +".qc.bam", "wb", template = pf)
		out_rm = pysam.AlignmentFile(outfile +".rm.bam", "wb", template = pf)
	
	#only keep paired reads, proper pairs
	#record the number of reads that pass and failed in quality filering.
	#p: pass; q: short read; n: excessive mismatches; m: strand error; r: wrong location; d: duplicate.
	p = q = n = m = r = d = 0
	#name of reads to be removed
	rm_qseq = {}
	for i in pf:
		if (i.is_secondary or i.reference_id < 0 or i.next_reference_id < 0):
			continue
		if (i.reference_name != "chrM" or i.next_reference_name != "chrM"):
			r += 1 
			rm_qseq[i.qname] = 1
			continue
		if (i.is_paired and (i.is_reverse != i.mate_is_reverse) and (not i.mate_is_unmapped)):
			#if (i.is_paired and (i.is_reverse != i.mate_is_reverse) and (not i.is_duplicate) and (not i.mate_is_unmapped)):
			if (i.is_duplicate):
				rm_qseq[i.qname] = 1
				d += 1
				continue
			if (i.query_length < min_query_len):
				rm_qseq[i.qname] = 1
				q+=1
				continue
			if (i.has_tag("NM")):
				nm = i.get_tag("NM")
				#only count match_len as the number of matched positions, excluding insertions and soft clips
				match_len = sum([j[1] for j in i.cigar if j[0] == 0])
				#count indel as one event
				indel_len = sum([j[1]-1 for j in i.cigar if j[0] in (1,2)])
				nm -= indel_len
				#if (nm > error_rate*i.query_length):
				#print >>sys.stderr, nm, match_len, error_rate*match_len
				if (nm > error_rate*match_len):
					rm_qseq[i.qname] = 1
					n += 1
					print >>sys.stderr, "Excessive mismatch", nm, match_len, i.qname, i.cigar
					continue
			#coding region reads must be in proper pair
			if ((min(i.next_reference_start,i.pos) > 560 + mtdna_offset or max(i.next_reference_start,i.pos)) < 16000 + mtdna_offset and not i.is_proper_pair):
				rm_qseq[i.qname] = 1
				m+=1
				continue
			p += 1
		else:
			rm_qseq[i.qname] = 1
			r += 1
			#count number of mismatches NM:i:%d
	print >>sys.stderr, "Removing low quality reads (%s):\t%d\t%d\t%d\t%d\t%d\t%d" % (name, p, d, q, n, m, r)
	pf.reset()
	for i in pf:
		if (not i.is_secondary and i.reference_id >= 0 and i.next_reference_id >= 0 and not rm_qseq.get(i.qname)):
			#output reads that pass the quality check.
			out_qc.write(i)
		else:
			#output reads removed to another file.
			out_rm.write(i)

	#close file handles and remove the temporary alignment file
	pf.close()
	out_qc.close()
	out_rm.close()
	if (pf_rm):
		execute("rm -rf %s" % pf_rm)

	return out_qc.filename

def generate_pileup_file(family, sample_name, outfile, gzip, alignment_file, mtdna_refseq, mtdna_offset, mapq_min = 20, qual_min = 0):
	"""
	perform base summarization

	Arguments
	----------
	family: family id
	sample_name: sample ids; samples will be output in the same order in the pileup file (list)
	outfile: output path and file name prefix
	gzip: compress file (true or false)
	alignment_file: names of the alignment files from retrieve_alignment
	mtdna_refseq: the reference mitochondrial genome sequence
	mtdna_offset: the position offset used in parsing mtDNA read alignments
	mapq_min: minimum mapping quality
	qual_min: minimum base alignment quality

	Returns
	----------
	None
	"""

	global override
	if (not override):
		if (os.path.exists(outfile + ".coverage") and os.path.exists(outfile + ".adj.pileup" + (".gz" if (gzip) else ""))):
			#return if the output file has already been created
			return

	#length of mitochondrial genome
	mtdna_len = 16569

	#pileup file handle
	pf = pipe_output("%s mpileup -q %d -Q %d -B -d 500000 -f %s %s" % (samtools, mapq_min, qual_min, mtdna_refseq, " ".join(alignment_file)))
	
	#summarize site coverage
	out_coverage = open(outfile + ".coverage", "w")
	#coverage file head line
	head = ["family","name","chr", "pos", "pos.adj", "ref", "depth", "Q0", "Q1", "Q2", "Q3", "Q4"]
	out_coverage.write("\t".join(head) + "\n")
	
	#trim and move the shifted reads to the correct positions of the original mitochondrial genome
	out_name =  outfile + ".adj.pileup"
	if (gzip):
		#compress output file
		out_name += ".gz"
	out_pileup = FileIO(out_name, "w", compresslevel=3)
	
	#temporarily store reads mapped to the end of the shifted mtDNA (the last mtdna_offset bps)
	tmp_line = [{} for i in range(len(sample_name))]
	
	#iterate reads in the pileup file generated
	nline = 0
	for line in pf.stdout:
		nline += 1
		line = line.rstrip("\r\n")
		if (not line):
			continue
		line = line.split("\t")
		assert len(line) == 3*len(sample_name)+3, "Invalid pileup file at line %d" % nline
		chr, pos, ref = line[:3]
		if (chr != "chrM"):
			continue
		pos = int(pos)
		#original position
		pos_adj = pos - mtdna_offset
		if (pos_adj > 0):
			#directly output if the site is not located in the shifted region
			out_pileup.write("\t".join([chr, str(pos_adj), ref]))
			for s in range(len(sample_name)):
				depth, r, q = line[(3*s+3):(3*s+6)]
				depth = int(depth)
				qual = [0,0,0,0,0]
				#group quals into <10, 10-20, 20-30, 30-40, >=40
				for i in phred(q):
					i = int(i)/10
					if (i >= 4):
						i = 4
					qual[i] += 1
				l = tmp_line[s].get(pos_adj)
				if (l):
					#output reads at the corresponding site of the shifted region
					chr1, ref1, depth1, qual1, r1, q1 = l
					assert ref1 == ref, "the reference allele does not match at position %d" % pos_adj
					#pileup reads if they aligned to the same positions
					depth = int(depth) + int(depth1)
					#concatenate reads and read qualities 
					r += r1
					q += q1
					#sum up quality stats
					qual = [i+j for i,j in zip(qual, qual1)]
					#delete temp records for the position
					del tmp_line[s][pos_adj]
				#output coverage and quality stats
				out_coverage.write("\t".join(map(str, [family, sample_name[s], chr, pos, pos_adj, ref, depth,] + qual)) + "\n")
				#output reads
				out_pileup.write("\t"+"\t".join([str(depth), r, q]))
			out_pileup.write("\n")
		else:
			#the site is not located in the shifted region
			for s in range(len(sample_name)):
				depth, r, q = line[(3*s+3):(3*s+6)]
				#temporarily store reads aligned to the last mtdna_offset bps
				depth = int(depth)
				qual = [0,0,0,0,0]
				#group quals into <10, 10-20, 20-30, 30-40, >40
				for i in phred(q):
					i = int(i)/10
					if (i >= 4):
						i = 4
					qual[i] += 1
				pos = mtdna_len + pos_adj
				tmp_line[s][pos] = [chr, ref, depth, qual, r, q]
	if (tmp_line):
		#output reads aligned to the last mtdna_offset bps
		for pos_adj in sorted(tmp_line[0].keys()):
			chr, ref, depth, qual, r, q = tmp_line[0][pos_adj]
			out_pileup.write("\t".join([chr, str(pos_adj), ref]))
			#for s in range(len(samples)):
			for s in range(len(sample_name)):
				chr, ref, depth, qual, r, q = tmp_line[s][pos_adj]
				#output coverage and quality stats
				out_coverage.write("\t".join(map(str, [family, sample_name[s], chr, pos_adj-mtdna_len, pos_adj, ref, depth,] + qual))+"\n")
				#output reads
				out_pileup.write("\t"+"\t".join([str(depth), r, q]))
		out_pileup.write("\n")
	
	#close file handles
	pf.stdout.close()
	out_pileup.close()
	out_coverage.close()

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

	parser = ArgumentParser(prog=prog, usage="%(prog)s [-options] [outpath]\n", version="%(prog)s v0.1")
	parser.add_argument("--fastq", type=str, nargs = "+", dest = "fastq_file", help="fastq files (in the format: name:R1:R2 ...)")
	parser.add_argument("--bam", type=str, nargs = "+", dest = "bam_file", help="bam files (in the format: name:bam ...)")
	parser.add_argument("--genome", type=str, dest="genome", help="specify the complete reference genome sequences")
	parser.add_argument("--mtdna", type=str, dest="mtdna", help="specify the mtDNA reference sequence")
	parser.add_argument("--mtdna-offset", type=int, default=0, dest="mtdna_offset", help="specify the coordinate offset when parsing mtDNA alignments")
	parser.add_argument("--remap-genome", action="store_true", default=False, dest="remap_genome", help="re-map reads from bam files to the reference genome sequences")
	parser.add_argument("--remap-mtdna", action="store_true", default=False, dest="remap_mtdna", help="re-map reads from bam files to the mtDNA sequence")
	parser.add_argument("--recal", action="store_true", default=False, dest="recal", help="enable BAQ re-calibration and indel re-alignment")
	parser.add_argument("--family", type=str, default="", dest="family", help="specify the family name of output sample(s)")
	parser.add_argument("--pileup", action = "store_true", default=False, dest="pileup", help="convert bam file to pileup file")
	parser.add_argument("-q","--mapq-min", type=int, dest="mapq_min", help="minimum mapping quality")
	parser.add_argument("-Q","--qual-min", type=int, dest="qual_min", help="minimum base alignment quality")
	parser.add_argument("-r", "--realign-recal", action = "store_true", default=False, dest="realign_recal", help="enable remap-genome, remap-mtdna recal and pileup")
	parser.add_argument("-z", "--gzip", action="store_true", default = False, dest = "gzip", help="compress the resulting pileup file with gzip")
	parser.add_argument("-s","--suffix", type=str, default="mtdna", dest="suffix", help="name suffix of the resulting pileup file (family.suffix.type)")
	parser.add_argument("--override", action="store_true", default = False, dest = "override", help="regenerate all files without checking file existence")
	parser.add_argument("outpath", help="specify the prefix of the output files")
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

	global bwa
	bwa = default_arguments.get("bwa", bwa)

	global bamleftalign
	bamleftalign = default_arguments.get("bamleftalign", bamleftalign)

	global picard
	picard = default_arguments.get("picard", picard)

	#retrieve file names
	alignment_file_all = []
	sample_name_all = []
	if (options.outpath.strip()):
		outpath = options.outpath
	else:
		outpath = "."
	outpath = outpath + os.path.sep
	outsuffix = ""
	if (options.suffix.strip()):
		outsuffix = options.suffix.strip()
		if (outsuffix[0] != "."):
			outsuffix = "." + outsuffix

	#initialize global variables
	global override
	override = options.override

	#initialize mapping parameters
	if (options.realign_recal):
		if (options.genome):
			options.remap_genome = True
			options.pileup = True
		if (options.mtdna):
			options.remap_mtdna = True
			options.pileup = True
		options.recal = True

	if (options.fastq_file):
		#obtain reads from fastq files
		for i in options.fastq_file:
			try:
				name, read1, read2 = i.split(":")
			except:
				print >>sys.stderr, "Skipping sample %s" % i
				continue
			if (name.strip()):
				alignment_file = retrieve_alignment(name, outpath + name, outpath + name, [read1, read2], options.remap_genome, options.genome, options.remap_mtdna, options.mtdna, options.mtdna_offset, options.recal)
				alignment_file_all.append(alignment_file)
				sample_name_all.append(name)
	elif (options.bam_file):
		#obtain reads from bam files
		for i in options.bam_file:
			try:
				name, bam = i.split(":")
			except:
				print >>sys.stderr, "Skipping sample %s" % i
				continue
			if (name.strip()):
				if (bam.endswith(".bam")):
					bam = bam[:-4]
				alignment_file = retrieve_alignment(name, bam, outpath + name, None, options.remap_genome, options.genome, options.remap_mtdna, options.mtdna, options.mtdna_offset, options.recal)
				alignment_file_all.append(alignment_file)
				sample_name_all.append(name)

	#initialize parameters for base summarization
	family = options.family
	if (not family):
		family = sample_name_all[0]
	#default minimum MAPQ is 20
	if (options.mapq_min is None or options.mapq_min < 0):
		if (options.recal and options.remap_mtdna):
			mapq_min = 0
		else:
			mapq_min = 20
	else:
		mapq_min = options.mapq_min
	#default minimum BAQ is zero
	if (not options.qual_min or options.qual_min < 0):
		qual_min = 0
	else:
		qual_min = options.qual_min
	#base summarization
	if (options.pileup and alignment_file_all):
		if (options.remap_mtdna):
			generate_pileup_file(family, sample_name_all, outpath + family + outsuffix, options.gzip, alignment_file_all, options.mtdna, options.mtdna_offset, mapq_min, qual_min)
		else:
			generate_pileup_file(family, sample_name_all, outpath + family + outsuffix, options.gzip, alignment_file_all, options.genome, options.mtdna_offset, mapq_min, qual_min)

def main():
	run(sys.argv[0], sys.argv[1:])

if __name__ == "__main__":
	main()