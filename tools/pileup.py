#!/usr/bin/python
#############################################################################
#Package: stamp
#Updated by Yiqin Wang
#Usage:
#
#
##############################################################################


import sys, os
from argparse import ArgumentParser

#basic file readers and parsers
from SeqIO import phred, SAMCigar
from FileIO import FileIO, execute, pipe_input, pipe_output, read_arguments

#executables, updated Jan, 2018
samtools = "samtools-1.6" 
#bwa ="bwa-0.7.17"
#bamleftalign="bamleftalign-1.1.0"

#alignment parsers
from align import retrieve_XM_info, retrieve_tags_info, retrieve_allele_info

#sequence distance calculation
from align import retrieve_alternate_sequences, retrieve_alternate_sequences_by_region, retrieve_alternate_sequences_by_similarity

#output alignment information
def output_alignment_info(out, info, tags_incl = ("XF","XQ","XM","XA","XP")):
	"""
	output the information of the alignment line
	
	Arguments
	----------
	out: the output file handle
	info: an alignment line (string) or a list of alignment line items 
	tags_incl: alignment annotation tags to be included for output
		
	Returns
	----------
	None
	
	"""

	if (isinstance(info, str)):
		info = info.rstrip("\r\n").split("\t")
	seq = info[:11] #QNAME to QUAL
	#exclude temporary alignment annotations
	annot = []
	for i in range(11,len(info)):
		s = info[i]
		if (s.startswith("X")):
			for j in tags_incl:
				if (s.startswith(j)):
					annot.append(s)
					break
		else:
			annot.append(s)
	out.write("\t".join(seq+annot)+"\n")

#filter consensus alignments according to the provided criteria

def summarize_alignment_file(alignment_file, tag_excl, tag_incl, mtdna_alt_seq, alt_fmin, alt_fmax):
	"""
	summarize the consensus read information without generating the pileup file.
	
	Arguments
	----------
	alignment_file: the consensus read alignment file(s) output from align
	tag_excl: exclude consensus reads with the tags specified
	tag_incl: include consensus reads with the tags specified
	mtdna_alt_seq: alternate mtdna sequences. i.e. retrieve_alternate_sequences(...)
	alt_fmin: the minimum fraction to exclude reads from NUMTS (default:0%)
	alt_fmax: the maximum fraction to exclude reads from NUMTS (default:100%)

	Returns
	----------
	a dict of
	alleles detected by each probe pair
	
	"""
	
	#count numbers of aligned consensus reads after QC
	alignment_summary = {}
	#alternate mtDNA sequence alleles marked by regions ("start-end")
	mtdna_alt_seq_by_region = {}
	#alternate mDNA sequences marked by amplicons ("probe")
	mtdna_alt_seq_by_probe = {}
	#store potential numts
	numts_seq = {}
	#store probe family size distribution
	probe_count = {}
	#store minor and alternate allele count
	allele_count = {}
	
	#iterate through all alignment files
	for i in alignment_file:
		mtdna_alignment = pipe_output("%s view %s" % (samtools, i))
		#iterate through lines of the consensus bam file
		for line in mtdna_alignment.stdout:
			info = line.rstrip("\r\n").split("\t")
			#retrieve barcode and amplicon information
			probe_info = retrieve_XM_info(info, paired = False)
			if (not probe_info):
				continue
			barcode, r1_probe, r2_probe = probe_info
			if (r1_probe != r2_probe):
				continue
			#retrieve tag information from mapping
			tags = retrieve_tags_info(info)
			#retrieve family size
			family_size = tags.get("XF")
			if (family_size is not None):
				family_size = int(family_size[1])
			else:
				continue
			#parse tag information
			xq_info = tags.get("XQ")
			if (xq_info is None or len(xq_info) < 2):
				qual_info = []
			else:
				qual_info = [i.upper() for i in xq_info[1].split(",")]
			#check appropriate quality or tag information for including and excluding reads
			if (tag_incl):
				n = 0
				for i in qual_info:
					for j in tag_incl:
						if (i.startswith(j)):
							n += 1
				if (n == 0):
					continue
			if (tag_excl):
				n = 0
				for i in qual_info:
					for j in tag_excl:
						if (i.startswith(j)):
							n += 1
				if (n > 0):
					#print info[0],qual_info, tag_excl
					continue
			
			#count the number of each probe retained
			if (r1_probe not in probe_count):
				probe_count[r1_probe] = 0
			probe_count[r1_probe] += 1
			
			#parse minor and alternate allele information
			minor_allele = retrieve_allele_info(tags.get("XA", ["",""])[1], remover_n = True)
			alt_allele = retrieve_allele_info(tags.get("XP",["",""])[1], remover_n = True)
			if (minor_allele):
				#mismatches compared to the mtDNA haplogroup
				nm = len(minor_allele)
				#calculate the length of mtDNA reference seq aligned with the read based on the Cigar string
				ref_length = 0
				for s, t in SAMCigar(info[5]).cigar:
					if (t in ["M","X","=","D"]):
						ref_length += s
				start = int(info[3])
				end = start+ref_length
				region = "%d-%d" % (start, end)
				alt_seq = mtdna_alt_seq_by_region.get(region)
				if (alt_seq is None):
					#temporarily store regional sequence information and amplicon information
					alt_seq = mtdna_alt_seq_by_region[region] = retrieve_alternate_sequences_by_region(mtdna_alt_seq, start, end)
					for id in alt_seq:
						if (id not in mtdna_alt_seq_by_probe):
							mtdna_alt_seq_by_probe[id] = {}
						if (r1_probe not in mtdna_alt_seq_by_probe[id]):
							mtdna_alt_seq_by_probe[id][r1_probe] = []
				if (alt_seq):
					nm_alt, id = retrieve_alternate_sequences_by_similarity(alt_seq, alt_allele)
					if (nm_alt < len(minor_allele)):
						#append NUMTS tag to quality information
						#qual_info.append("NUMTS(%s)"%id)
						#count numbers of alternate sequences found for "probe"
						#temporarily store the alignment information
						mtdna_alt_seq_by_probe[id][r1_probe].append([family_size, len(minor_allele), len(alt_allele)])
						continue
				else:
					is_numts = False
					for i in qual_info:
						if (i.startswith("NUMTS")):
							#trim off the NUMTS tag 
							id = i[6:-1]
							if (id not in mtdna_alt_seq_by_probe):
								mtdna_alt_seq_by_probe[id] = {}
							if (r1_probe not in mtdna_alt_seq_by_probe[id]):
								mtdna_alt_seq_by_probe[id][r1_probe] = []
							#temporarily store the alignment information
							mtdna_alt_seq_by_probe[id][r1_probe].append(info)
							is_numts = True
					if (is_numts):
						continue
			if (r1_probe not in allele_count):
				allele_count[r1_probe] = {}
			if (family_size not in allele_count[r1_probe]):
				allele_count[r1_probe][family_size] = [[],[]]
			allele_count[r1_probe][family_size][0].append(len(minor_allele))
			allele_count[r1_probe][family_size][1].append(len(alt_allele))
		mtdna_alignment.stdout.close()
	
	for id in mtdna_alt_seq_by_probe:
		#print id, mtdna_alt_seq_by_probe[id]
		for r1_probe in mtdna_alt_seq_by_probe[id]:
			n = len(mtdna_alt_seq_by_probe[id][r1_probe])
			freq = n/float(probe_count[r1_probe])
			if (freq >= alt_fmin and freq <= alt_fmax):
				if (n > 0):
					print id, r1_probe, n
			else: 
				#print id, r1_probe, n, "output"
				for family_size, minor_allele_nm, alt_allele_nm in mtdna_alt_seq_by_probe[id][r1_probe]:
					if (r1_probe not in allele_count):
						allele_count[r1_probe] = {}
					if (family_size not in allele_count[r1_probe]):
						allele_count[r1_probe][family_size] = [[],[]]
					allele_count[r1_probe][family_size][0].append(minor_allele_nm)
					allele_count[r1_probe][family_size][1].append(alt_allele_nm)
	
	return allele_count

def update_qual_info(info, add_val, remove_val):
	"""
	update the XQ tag in the alignment
	
	Arguments
	----------
	info: one line in sam/bam files (list)
	add_val: new annotations to add (list or string)
	remove_val: annotations to remove (list of string)
	
	Returns
	----------
	updated alignment (list)
	
	"""
	
	if (not add_val):
		add_val = []
	elif (not isinstance(add_val, list) or not isinstance(add_val, set)):
		add_val = [add_val,]
	
	if (not remove_val):
		remove_val = []
	elif (not isinstance(remove_val, list) or not isinstance(remove_val, set)):
		remove_val = [remove_val,]
	
	for i in range(11,len(info)):
		if (info[i].startswith("XQ")):
			xq_info = info[i].split(":", 2)
			qual_info = xq_info[2].split(",")
			for j in remove_val:
				if (j in qual_info):
					qual_info.remove(j)
			for j in add_val:
				if (j not in add_val):
					qual_info.append(j)
			info[i] = xq_info[0]+":"+xq_info[1]+":"+",".join(qual_info)
			break
	else:
		if (add_val):
			info.append("XQ:Z:"+",".join(add_val))
	return info

def annotate_alignment_file(sample_name, outpath, alignment_file, family_size_min, family_size_max, nm_max, nm_max_dloop, mtdna_offset, tag_excl, tag_incl, mtdna_alt_seq, alt_fmin, alt_fmax):
	"""
	perform quality control annotation on the alignments of consensus reads.
	all consensus reads are retained in the resulting alignment file.
	Reads that pass quality control will be marked as "PASS" in the XQ tag, which can be extracted using "--tag-incl PASS"
	
	Arguments
	----------
	sample_name: the name of the sample
	outpath: the path to store the output files
	alignment_file: the consensus read alignment file(s) output from align
	family_size_min: the minimum read family size of consensus reads (read family size refers the number of paired-end reads used to construct the consensus read)
	family_size_max: the maximum read family size of consensus reads
	nm_max: the maximum number of mismatches of consensus reads to the major mtDNA sequence in the coding region
	nm_max_dloop: the maximum number of mismatches of consensus reads to the major mtDNA sequence in the dloop region
	mtdna_offset: the position offset used in parsing mtDNA read alignments
	tag_excl: exclude consensus reads with the tags specified
	tag_incl: include consensus reads with the tags specified
	mtdna_alt_seq: alternate mtdna sequences. i.e. retrieve_alternate_sequences(...)
	alt_fmin: the minimum fraction to exclude reads from NUMTS (default:0%)
	alt_fmax: the maximum fraction to exclude reads from NUMTS (default:100%)
	
	Returns
	----------
	a tuple of
	1. the path to the processed alignment file
	2. a dict of reads passing QC 
	"""
	#dump QC+ reads to a new bam file 
	cmd = "%s view -hbS -o %s - " % (samtools, outpath + os.path.sep + sample_name + ".%s.consensus.bam")
	
	#retain head information
	head_mtdna = pipe_output("%s view -H %s" % (samtools, alignment_file[0]))
	
	#output consensus alignments in bam format
	if (not outpath):
		outpath = "."
	
	#pipe to samtools for bam format conversion
	#out_alignment_file = outpath + os.path.sep + sample_name + ".mtdna.consensus.bam"
	out_alignment_file = outpath + os.path.sep + sample_name + ".unsorted.mtdna.consensus.bam"
	out_mtdna_alignment = pipe_input("%s view -hbS -o %s - " % (samtools, out_alignment_file))
	
	#retain head information
	out_mtdna_alignment.stdin.write(head_mtdna.stdout.read())
	head_mtdna.stdout.close()
	
	#set an accetable range of family size
	fs_check = True
	if (family_size_min is None):
		if (family_size_max is None):
			fs_check = False
		else:
			family_size_min = 1
	else:
		if (family_size_max is None):
			family_size_max = sys.maxint
	
	#set an accetable range of nucleotide mismatches
	if (nm_max is None or nm_max < 0):
		nm_max = sys.maxint
		
	if (nm_max_dloop is None or nm_max_dloop < 0):
		nm_max_dloop = sys.maxint
	
	#count numbers of aligned consensus reads after QC
	alignment_summary = {}
	#alternate mtDNA sequence alleles marked by regions ("start-end")
	mtdna_alt_seq_by_region = {}
	#aleternate mDNA sequences marked by amplicons ("probe")
	mtdna_alt_seq_by_probe = {}
	#store potential numts
	numts_seq = {}
	#store probe family size distribution
	probe_count = {}
	#iterate through all alignment files
	ntotal = 0
	noutput = 0
	for i in alignment_file:
		mtdna_alignment = pipe_output("%s view %s" % (samtools, i))
		#iterate through lines of the consensus bam file
		for line in mtdna_alignment.stdout:
			ntotal += 1
			info = line.rstrip("\r\n").split("\t")
			#retrieve barcode and amplicon information
			probe_info = retrieve_XM_info(info, paired = False)
			if (not probe_info):
				update_qual_info(info, "NOBARCODE", "PASS")
				output_alignment_info(out_mtdna_alignment.stdin, info)
				continue
			barcode, r1_probe, r2_probe = probe_info
			if (r1_probe != r2_probe):
				update_qual_info(info, "OFFTARGET", "PASS")
				output_alignment_info(out_mtdna_alignment.stdin, info)
				continue
			#retrieve tag information from mapping
			tags = retrieve_tags_info(info)
			#check for appropriate family size
			if (fs_check):
				family_size = tags.get("XF")
				if (family_size is not None):
					family_size = int(family_size[1])
					if (family_size > family_size_max or family_size < family_size_min):
						update_qual_info(info, "IMPROPERFS", "PASS")
						output_alignment_info(out_mtdna_alignment.stdin, info)
						continue
				else:
					update_qual_info(info, "IMPROPERFS", "PASS")
					output_alignment_info(out_mtdna_alignment.stdin, info)
					continue
			#parse tag information
			xq_info = tags.get("XQ")
			if (xq_info is None or len(xq_info) < 2):
				qual_info = []
			else:
				qual_info = [i.upper() for i in xq_info[1].split(",")]
			
			#check appropriate quality or tag information for including and excluding reads
			if (tag_incl):
				#only consider reads with tags provided
				n = 0
				for i in qual_info:
					for j in tag_incl:
						if (i.startswith(j)):
							n += 1
				if (n == 0):
					update_qual_info(info, "", "PASS")
					output_alignment_info(out_mtdna_alignment.stdin, info)
					continue
			else:
				if (tag_excl):
					n = 0
					for i in qual_info:
						for j in tag_excl:
							if (i.startswith(j)):
								n += 1
					if (n > 0):
						update_qual_info(info, "EXCLTAG", "PASS")
						output_alignment_info(out_mtdna_alignment.stdin, info)
						#print info[0],qual_info, tag_excl
						continue
				
				#count the number of each probe retained
				if (r1_probe not in probe_count):
					probe_count[r1_probe] = 0
				probe_count[r1_probe] += 1
				
				#parse minor and alternate allele information
				minor_allele = retrieve_allele_info(tags.get("XA", ["",""])[1], remover_n = True)
				alt_allele = retrieve_allele_info(tags.get("XP",["",""])[1], remover_n = True)
				if (minor_allele):
					#mismatches compared to the major mtDNA sequence
					nm = len(minor_allele)
					#calculate the length of mtDNA reference seq aligned with the read based on the Cigar string
					ref_length = 0
					for s, t in SAMCigar(info[5]).cigar:
						if (t in ["M","X","=","D"]):
							ref_length += s
					start = int(info[3])
					end = start+ref_length
					if (start >= 16000 + mtdna_offset or end < 560 + mtdna_offset):
						#within D-loop region
						if (nm_max_dloop is not None):
							if (nm > nm_max_dloop):
								update_qual_info(info, "EXMISMATCH", "PASS")
								output_alignment_info(out_mtdna_alignment.stdin, info)
								continue
					else:
						#within coding region
						if (nm_max is not None):
							if (nm > nm_max):
								update_qual_info(info, "EXMISMATCH", "PASS")
								output_alignment_info(out_mtdna_alignment.stdin, info)
								continue
					region = "%d-%d" % (start, end)
					alt_seq = mtdna_alt_seq_by_region.get(region)
					if (alt_seq is None):
						#temporarily store regional sequence information and amplicon information
						alt_seq = mtdna_alt_seq_by_region[region] = retrieve_alternate_sequences_by_region(mtdna_alt_seq, start, end, overlap_min = 0.9)
						for id in alt_seq:
							if (id not in mtdna_alt_seq_by_probe):
								mtdna_alt_seq_by_probe[id] = {}
							if (r1_probe not in mtdna_alt_seq_by_probe[id]):
								mtdna_alt_seq_by_probe[id][r1_probe] = []
					if (alt_seq):
						nm_alt, id = retrieve_alternate_sequences_by_similarity(alt_seq, alt_allele)
						if (nm_alt < len(minor_allele)):
							#append NUMTS tag to quality information
							#qual_info.append("NUMTS(%s)"%id)
							#count numbers of alternate sequences found for "probe"
							#temporarily store the alignment information
							mtdna_alt_seq_by_probe[id][r1_probe].append(info)
							update_qual_info(info, "NUMTS(%s)"%id, "PASS")
							#for i in range(11,len(info)):
							#	if (info[i].startswith("XQ")):
							#		xq_info = xq_info[1].split(",")
							#		if ("PASS" in xq_info):
							#			xq_info.remove("PASS")
							#		xq_info.append("NUMTS(%s)"%id)
							#		info[i] = "XQ:Z:%s"%",".join(xq_info)
							#		break
							#else:
							#	info.append("XQ:Z:NUMTS(%s)"%id)
							continue
					else:
						is_numts = False
						for i in qual_info:
							if (i.startswith("NUMTS")):
								#trim off the NUMTS tag 
								id = i[6:-1]
								if (id not in mtdna_alt_seq_by_probe):
									mtdna_alt_seq_by_probe[id] = {}
								if (r1_probe not in mtdna_alt_seq_by_probe[id]):
									mtdna_alt_seq_by_probe[id][r1_probe] = []
								#temporarily store the alignment information
								mtdna_alt_seq_by_probe[id][r1_probe].append(info)
								is_numts = True
						if (is_numts):
							continue
				#elif alt_allele:
				#	for pos in alt_allele:
				#		if (pos not in major_allele):
				#			major_allele[pos] = alt_allele[pos]
				
			#remaining reads were ouput to a bam file for generating the pileup file
			#out_mtdna_alignment.stdin.write(line)
			output_alignment_info(out_mtdna_alignment.stdin, info)
			noutput += 1
			if (r1_probe not in alignment_summary):
				alignment_summary[r1_probe] = 0
			alignment_summary[r1_probe] += 1
		mtdna_alignment.stdout.close()
	
	info_out = {}
	#alignment not passing QC of NUMTS
	info_numts = {}
	#iterate all potential NUMTS reads
	for id in mtdna_alt_seq_by_probe:
		#print id, mtdna_alt_seq_by_probe[id]
		for r1_probe in mtdna_alt_seq_by_probe[id]:
			n = len(mtdna_alt_seq_by_probe[id][r1_probe])
			freq = min(1.0, n/float(probe_count[r1_probe]))
			if (freq >= alt_fmin and freq <= alt_fmax):
				if (n > 0):
					print id, r1_probe, n
				for info in mtdna_alt_seq_by_probe[id][r1_probe]:
					#temporarity store numts reads
					info_numts[info[0]]=info
				#discard this consensus read
			else: 
				#the frequency of the minor allele does not fall into the fraction range of NUMTS
				#retain this consensus read
				#print id, r1_probe, n, "output"
				for info in mtdna_alt_seq_by_probe[id][r1_probe]:
					#out_mtdna_alignment.stdin.write(line)
					if (info[0] not in info_out):
						#set QC status as PASS
						update_qual_info(info, "PASS", "")
						#output each consensus read only once
						output_alignment_info(out_mtdna_alignment.stdin, info)
						info_out[info[0]] = 1
						noutput += 1
	
	#output consensus reads with NUMTS annotation
	for info in info_numts.values():
		if (info[0] not in info_out):
			#output each consensus read only once
			update_qual_info(info, "", "PASS")
			output_alignment_info(out_mtdna_alignment.stdin, info)
	
	out_mtdna_alignment.stdin.close()
	
	#sort reads according to the aligned mtDNA positions
	out_alignment_sorted_file = outpath + os.path.sep + sample_name + ".mtdna.consensus.bam"
	execute("%s sort -o %s %s " % (samtools, out_alignment_sorted_file, out_alignment_file))
	
	#remove temporary unsorted alignment files to save space
	execute("rm -f %s" % out_alignment_file)	

	print ntotal, noutput
	
	#return out_alignment_file, alignment_summary
	return out_alignment_sorted_file, alignment_summary

def filter_alignment_file(sample_name, outpath, alignment_file, family_size_min, family_size_max, nm_max, nm_max_dloop, mtdna_offset, tag_excl, tag_incl, mtdna_alt_seq, alt_fmin, alt_fmax):
	"""
	perform quality control filtering on the alignments of consensus reads.
	only consensus reads passing quality control are retained in the resulting alignment file.
	
	Arguments
	----------
	sample_name: the name of the sample
	outpath: the path to store the output files
	alignment_file: the consensus read alignment file(s) output from align
	family_size_min: the minimum read family size of consensus reads (read family size refers the number of paired-end reads used to construct the consensus read)
	family_size_max: the maximum read family size of consensus reads
	nm_max: the maximum number of mismatches of consensus reads to the major mtDNA sequence in the coding region
	nm_max_dloop: the maximum number of mismatches of consensus reads to the major mtDNA sequence in the dloop region
	mtdna_offset: the position offset used in parsing mtDNA read alignments
	tag_excl: exclude consensus reads with the tags specified
	tag_incl: include consensus reads with the tags specified
	mtdna_alt_seq: alternate mtdna sequences. i.e. retrieve_alternate_sequences(...)
	alt_fmin: the minimum fraction to exclude reads from NUMTS (default:0%)
	alt_fmax: the maximum fraction to exclude reads from NUMTS (default:100%)
	
	Returns
	----------
	a tuple of
	1. the path to the processed alignment file
	2. a dict of reads output
	
	"""
	#dump QC+ reads to a new bam file 
	cmd = "%s view -hbS -o %s - " % (samtools, outpath + os.path.sep + sample_name + ".%s.consensus.bam")
	
	#retain head information
	head_mtdna = pipe_output("%s view -H %s" % (samtools, alignment_file[0]))
	
	#output consensus alignments in bam format
	if (not outpath):
		outpath = "."
	
	#pipe to samtools for bam format conversion
	#out_alignment_file = outpath + os.path.sep + sample_name + ".mtdna.consensus.bam"
	out_alignment_file = outpath + os.path.sep + sample_name + ".unsorted.mtdna.consensus.bam"
	out_mtdna_alignment = pipe_input("%s view -hbS -o %s - " % (samtools, out_alignment_file))
	#retain head imformation
	out_mtdna_alignment.stdin.write(head_mtdna.stdout.read())
	head_mtdna.stdout.close()
	
	#set an accetable range of family size
	fs_check = True
	if (family_size_min is None):
		if (family_size_max is None):
			fs_check = False
		else:
			family_size_min = 1
	else:
		if (family_size_max is None):
			family_size_max = sys.maxint
	
	#set an accetable range of nucleotide mismatches
	if (nm_max is None or nm_max < 0):
		nm_max = sys.maxint
		
	if (nm_max_dloop is None or nm_max_dloop < 0):
		nm_max_dloop = sys.maxint

	#count numbers of aligned consensus reads after QC
	alignment_summary = {}
	#alternate mtDNA sequence alleles marked by regions ("start-end")
	mtdna_alt_seq_by_region = {}
	#aleternate mDNA sequences marked by amplicons ("probe")
	mtdna_alt_seq_by_probe = {}
	#store potential numts
	numts_seq = {}
	#store probe family size distribution
	probe_count = {}
	#iterate through all alignment files
	ntotal = 0
	noutput = 0
	for i in alignment_file:
		mtdna_alignment = pipe_output("%s view %s" % (samtools, i))
		#iterate through lines of the consensus bam file
		for line in mtdna_alignment.stdout:
			ntotal += 1
			info = line.rstrip("\r\n").split("\t")
			#retrieve barcode and amplicon information
			probe_info = retrieve_XM_info(info, paired = False)
			if (not probe_info):
				continue
			barcode, r1_probe, r2_probe = probe_info
			if (r1_probe != r2_probe):
				continue
			#retrieve tag information from mapping
			tags = retrieve_tags_info(info)
			#check for appropriate family size
			if (fs_check):
				family_size = tags.get("XF")
				if (family_size is not None):
					family_size = int(family_size[1])
					if (family_size > family_size_max or family_size < family_size_min):
						continue
				else:
					continue
			#parse tag information
			xq_info = tags.get("XQ")
			if (xq_info is None or len(xq_info) < 2):
				qual_info = []
			else:
				qual_info = [i.upper() for i in xq_info[1].split(",")]
			
			#check appropriate quality or tag information for including and excluding reads 	
			if (tag_incl):
				#only consider reads with tags provided
				n = 0
				for i in qual_info:
					for j in tag_incl:
						if (i.startswith(j)):
							n += 1
				if (n == 0):
					continue
			else:
				if (tag_excl):
					n = 0
					for i in qual_info:
						for j in tag_excl:
							if (i.startswith(j)):
								n += 1
					if (n > 0):
						#print info[0],qual_info, tag_excl
						continue
				
				#count the number of each probe retained
				if (r1_probe not in probe_count):
					probe_count[r1_probe] = 0
				probe_count[r1_probe] += 1
				
				#parse minor and alternate allele information
				minor_allele = retrieve_allele_info(tags.get("XA", ["",""])[1], remover_n = True)
				alt_allele = retrieve_allele_info(tags.get("XP",["",""])[1], remover_n = True)
				if (minor_allele):
					#mismatches compared to the major mtDNA sequence
					nm = len(minor_allele)
					#calculate the length of mtDNA reference seq aligned with the read based on the Cigar string
					ref_length = 0
					for s, t in SAMCigar(info[5]).cigar:
						if (t in ["M","X","=","D"]):
							ref_length += s
					start = int(info[3])
					end = start+ref_length
					if (start >= 16000 + mtdna_offset or end < 560 + mtdna_offset):
						#within D-loop region
						if (nm_max_dloop is not None):
							if (nm > nm_max_dloop):
								continue
					else:
						#within coding region
						if (nm_max is not None):
							if (nm > nm_max):
								continue
					region = "%d-%d" % (start, end)
					alt_seq = mtdna_alt_seq_by_region.get(region)
					if (alt_seq is None):
						#temporarily store regional sequence information and amplicon information
						alt_seq = mtdna_alt_seq_by_region[region] = retrieve_alternate_sequences_by_region(mtdna_alt_seq, start, end, overlap_min = 0.9)
						for id in alt_seq:
							if (id not in mtdna_alt_seq_by_probe):
								mtdna_alt_seq_by_probe[id] = {}
							if (r1_probe not in mtdna_alt_seq_by_probe[id]):
								mtdna_alt_seq_by_probe[id][r1_probe] = []
					if (alt_seq):
						nm_alt, id = retrieve_alternate_sequences_by_similarity(alt_seq, alt_allele)
						#print nm_alt, r1_probe, region,id, minor_allele, alt_allele, alt_seq
						if (nm_alt < len(minor_allele)):
							#append NUMTS tag to quality information
							#qual_info.append("NUMTS(%s)"%id)
							#count numbers of alternate sequences found for "probe"
							#temporarily store the alignment information
							mtdna_alt_seq_by_probe[id][r1_probe].append(info)
							for i in range(11,len(info)):
								if (info[i].startswith("XQ")):
									xq_info = xq_info[1].split(",")
									if ("PASS" in xq_info):
										xq_info.remove("PASS")
									xq_info.append("NUMTS(%s)"%id)
									info[i] = "XQ:Z:%s"%",".join(xq_info)
									break
							else:
								info.append("XQ:Z:NUMTS(%s)"%id)
							continue
					else:
						is_numts = False
						for i in qual_info:
							if (i.startswith("NUMTS")):
								#trim off the NUMTS tag 
								id = i[6:-1]
								if (id not in mtdna_alt_seq_by_probe):
									mtdna_alt_seq_by_probe[id] = {}
								if (r1_probe not in mtdna_alt_seq_by_probe[id]):
									mtdna_alt_seq_by_probe[id][r1_probe] = []
								#temporarily store the alignment information
								mtdna_alt_seq_by_probe[id][r1_probe].append(info)
								is_numts = True
						if (is_numts):
							continue
				#elif alt_allele:
				#	for pos in alt_allele:
				#		if (pos not in major_allele):
				#			major_allele[pos] = alt_allele[pos]
			
			#remaining reads were ouput to a bam file for generating the pileup file
			#out_mtdna_alignment.stdin.write(line)
			output_alignment_info(out_mtdna_alignment.stdin, info)
			noutput += 1
			if (r1_probe not in alignment_summary):
				alignment_summary[r1_probe] = 0
			alignment_summary[r1_probe] += 1
		mtdna_alignment.stdout.close()
	
	#print probe_count
	info_out = {}
	#iterate all potential NUMTS reads
	for id in mtdna_alt_seq_by_probe:
		#print id, mtdna_alt_seq_by_probe[id]
		for r1_probe in mtdna_alt_seq_by_probe[id]:
			n = len(mtdna_alt_seq_by_probe[id][r1_probe])
			freq = min(1.0, n/float(probe_count[r1_probe]))
			if (freq >= alt_fmin and freq <= alt_fmax):
				if (n > 0):
					#discard this consensus read
					print id, r1_probe, n
			else: 
				#the frequency of the minor allele does not fall into the fraction range of NUMTS
				#retain this consensus read
				#print id, r1_probe, n, "@@"
				for info in mtdna_alt_seq_by_probe[id][r1_probe]:
					#print line 
					#out_mtdna_alignment.stdin.write(line)
					if (info[0] not in info_out):
						#output each consensus read only once
						output_alignment_info(out_mtdna_alignment.stdin, info)
						info_out[info[0]] = 1
						noutput += 1
	
	out_mtdna_alignment.stdin.close()
	
	#sort reads according to the aligned mtDNA positions
	out_alignment_sorted_file = outpath + os.path.sep + sample_name + ".mtdna.consensus.bam"
	execute("%s sort -o %s %s " % (samtools, out_alignment_sorted_file, out_alignment_file))
	
	#remove temporary unsorted alignment files to save space
	execute("rm -f %s" % out_alignment_file)	

	print ntotal, noutput
	
	#return out_alignment_file, alignment_summary
	return out_alignment_sorted_file, alignment_summary


def generate_pileup_file(sample_name, outpath, gzip_pileup, alignment_file, alignment_summary, probe_file, mtdna_refseq, qual_min, mtdna_offset):
	"""
	a function wrapper to generate the pileup file from the alignment file using samtools
	mtdna positions are corrected to those of rCRS
	
	Arguments
	----------
	sample_name: the name of the sample
	outpath: the path to store the output files
	gzip_pileup: compress the output file
	alignment_file: the processed alignment file returned from filter_alignment_file(...)
	alignment_summary: a dict of read summary returned from filter_alignment_file(...)
	probe_file: the path to the probe file
	mtdna_refseq: the mtdna reference sequence
	qual_min: the minimum quality score to output (used in samtools mpileup -Q )
	mtdna_offset: the position offset used to parse mtdna sites.
	
	Returns
	----------
	None
	
	Outputs
	----------
	${output}/${sample_name}.mtdna.consensus.adj.pileup(.gz): the resulting pileup file
	${output}/${sample_name}.coverage: the read coverage information for all mtdna sites (tsv file)  
	
	"""
	
	#parse amplicon information
	mtdna_len = 16569
	read_len = 250
	
	#amplicon information for each position in mtdna
	amp_info = [[] for i in xrange(mtdna_len+mtdna_offset+1)] #amplicon covered at each position
	amp_cov = [0,]*(mtdna_len+mtdna_offset+1) #amplicons sequencing depth at each position
	amp_r1_pos = [mtdna_len,]*(mtdna_len+mtdna_offset+1) #relative position at read 1
	amp_r2_pos = [mtdna_len,]*(mtdna_len+mtdna_offset+1) #relative position at read 2
	amp_r1_probe = [mtdna_len,]*(mtdna_len+mtdna_offset+1) #relative position to the r1 probe
	amp_r2_probe = [mtdna_len,]*(mtdna_len+mtdna_offset+1) #relative position to the r2 probe
	
	#parse probe file for amplicon imformation
	with open(probe_file, "r") as fh:
		for line in fh:
			line = line.rstrip("\r\n")
			if (not line):
				continue
			name, chr, start, end, s1, s2, r1_probe, r2_probe, blen = line.split("\t")
			start = int(start)
			end = int(end)
			#length of r1 and r2 probes
			if (s1 == "+"):
				p1 = len(r1_probe.strip())
				p2 = len(r2_probe.strip())
			else:
				p1 = len(r2_probe.strip())
				p2 = len(r1_probe.strip())
			#barcode length
			blen = int(blen)
			if (chr == "chrM"):
				#number of amplicons in the QC+ bam file
				cov = alignment_summary.get(name, 0)
				amp = []
				if (start < 0):
					#split the amplicon into halves in the D-loop region
					amp.append([mtdna_len+start, mtdna_len, s1, p1, 0])
					amp.append([1, end, s1, 0, p2])
				else:
					amp.append([start, end, s1, p1, p2])
				#positions in r1 and r2 reads
				#positions in probe
				for start, end, s1, p1, p2 in amp:
					for i in range(start+p1, end-p2+1):
						amp_info[i].append("%s(%s)"%(name,s1))
						amp_cov[i] += cov
					if (s1 == "+"):
						for i in xrange(p1):
							amp_r1_probe[i+start] = p1-i #position in R1 probe
						for i in xrange(start+p1, end-p2+1):
							amp_r1_pos[i] = min(amp_r1_pos[i], i-start) #position in R1
							amp_r2_pos[i] = min(amp_r2_pos[i], end+1-i+blen) #position in R2
						for i in xrange(p2): #position in R2 probe
							amp_r2_probe[end-i] = p2-i
					else:
						for i in xrange(p1):
							amp_r2_probe[i+start] = p1-i
						for i in xrange(start+p1, end-p2+1):
							amp_r2_pos[i] = min(amp_r2_pos[i], i-start+blen)
							amp_r1_pos[i] = min(amp_r1_pos[i], end+1-i)
						for i in xrange(p2):
							amp_r1_probe[end-i] = p2-i
	
	if (alignment_file.endswith(".bam")):
		alignment_file = alignment_file[:-4]
	
	#sort reads according to the aligned mtDNA positions
	#execute("%s sort -o %s.sorted.bam %s.bam " % (samtools, alignment_file, alignment_file))
	
	#pileup reads using samtools
	#mapq >= 20 & baseq >= qual_min
	#pf = pipe_output("%s mpileup -q 20 -Q %d -B -d 500000 -f %s %s.sorted.bam" % (samtools, qual_min, mtdna_refseq, alignment_file))
	pf = pipe_output("%s mpileup -q 20 -Q %d -B -d 500000 -f %s %s.bam" % (samtools, qual_min, mtdna_refseq, alignment_file))
	
	#summarize site coverage
	out_coverage = open(outpath + os.path.sep + sample_name + ".coverage", "w")
	head = ["chr", "pos", "pos.adj", "ref", "depth", "Q0", "Q1", "Q2", "Q3", "Q4", "amps", "amp.r1.pos","amp.r2.pos","amp.r1.probe","amp.r2.probe","amp.info"]
	out_coverage.write("\t".join(head) + "\n")
	
	#trim and move the shifted reads to the correct rCRS positions
	out_name =  outpath + os.path.sep + sample_name + ".mtdna.consensus.adj.pileup"
	if (gzip_pileup):
		out_name += ".gz"
	out_pileup = FileIO(out_name, "w", compresslevel=3)
	
	#temporarily store amplicons mapped to the end of the shifted mtDNA (the last mtdna_offset bps)
	tmp_line = {}
	
	#iterate reads in the pileup file generated
	for line in pf.stdout:
		line = line.rstrip("\r\n")
		if (not line):
			continue
		chr, pos, ref, depth, r, q = line.split("\t")
		qual = [0,0,0,0,0]
		#group quals into <10, 10-20, 20-30, 30-40, >40
		for i in phred(q):
			i = int(i)/10
			if (i >= 4):
				i = 4
			qual[i] += 1
		pos = int(pos)
		depth = int(depth)
		pos_adj = pos - mtdna_offset
		if (pos_adj > 0):
			l = tmp_line.get(pos_adj)
			if (l):
				chr1, ref1, depth1, qual1, r1, q1 = l
				assert ref1 == ref, "the reference allele does not match at position %d" % pos_adj
				#pileup reads if they aligned to the same positions
				depth = int(depth) + int(depth1)
				#concatenate reads and read qualities, respectively
				r += r1
				q += q1
				#sum up quality stats
				qual = [i+j for i,j in zip(qual, qual1)]
				#delete temp records for the position
				del tmp_line[pos_adj]
			#output coverage and quality stats
			out_coverage.write("\t".join(map(str, [chr, pos, pos_adj, ref, depth] + qual + [amp_cov[pos_adj], amp_r1_pos[pos_adj], amp_r2_pos[pos_adj], amp_r1_probe[pos_adj], amp_r2_probe[pos_adj],"|".join(amp_info[pos_adj])]))+"\n")
			#output reads
			out_pileup.write("\t".join([chr, str(pos_adj), ref, str(depth), r, q])+"\n")
		else:
			#temporarily store reads aligned to the last mtdna_offset bps
			pos = mtdna_len + pos_adj
			tmp_line[pos] = [chr, ref, depth, qual, r, q]
	if (tmp_line):
		#output reads aligned to the last mtdna_offset bps
		for pos_adj in sorted(tmp_line.keys()):
			chr, ref, depth, qual, r, q = tmp_line[pos_adj]
			out_coverage.write("\t".join(map(str, [chr, pos_adj-mtdna_len, pos_adj, ref, depth] + qual + [amp_cov[pos_adj], amp_r1_pos[pos_adj], amp_r2_pos[pos_adj], amp_r1_probe[pos_adj], amp_r2_probe[pos_adj],"|".join(amp_info[pos_adj])]))+"\n")
			out_pileup.write("\t".join([chr, str(pos_adj), ref, str(depth), r, q])+"\n")
	
	#close file handles
	pf.stdout.close()
	out_pileup.close()
	out_coverage.close()
	
	#remove temporary sorted alignment files to save space
	#execute("rm -f %s.sorted.bam" % alignment_file)	

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
	parser = ArgumentParser(prog=prog, usage="%(prog)s [-options] sample\n", version="%(prog)s v0.1.1")
	parser.add_argument("-a", "--alignment", type = str, nargs='+', dest = "alignment_file", help="the consensus read alignment file(s) output from align")
	parser.add_argument("-p", "--probe", type = str, dest = "probe", help="the file with the probe information")
	parser.add_argument("-o", "--outpath", type = str, default = ".", dest = "outpath", help="path where to store the pileup file")
	parser.add_argument("-z", "--gzip", action="store_true", default = False, dest = "gzip", help="compress the pileup file with gzip")
	parser.add_argument("--mtdna", type=str, dest="mtdna", help="the mtDNA reference sequence")
	parser.add_argument("--mtdna-offset", type=int, dest="mtdna_offset", help="the position offset used in parsing mtDNA read alignments")
	parser.add_argument("--fs-min", type=int, default = None, dest = "fs_min", help="the minimum read family size of consensus reads")
	parser.add_argument("--fs-max", type=int, default = None, dest = "fs_max", help="the maximum read family size of consensus reads")
	parser.add_argument("--nm-max", type=float, default = None, dest = "nm_max", help="the maximum number of mismatches of consensus reads to the major mtDNA sequence in the coding region")
	parser.add_argument("--nm-max-dloop", type=float, default = None, dest = "nm_max_dloop", help="the maximum number of mismatches of consensus reads to the major mtDNA sequence in the Dloop region")
	parser.add_argument("--tag-excl", type=str, default="NUMTS(BWA,EXMISMATCH", dest = "tag_excl", help = "exclude consensus reads with the tags specified")
	parser.add_argument("--tag-incl", type=str, default="", dest="tag_incl", help="include consensus reads with the tags specified")
	parser.add_argument("--numts-excl", type=str, dest="numts", help="exclude reads from nuclear mitochondrial DNA segments specified by the HSD file")
	parser.add_argument("--numts-fmin", type=float, default = 0.0, dest="numts_fmin", help="the minimum fraction to exclude reads from NUMTS")
	parser.add_argument("--numts-fmax", type=float, default = 1.0, dest="numts_fmax", help="the maximum fraction to exclude reads from NUMTS")
	parser.add_argument("--qual", type=int, default = 0, dest="qual_min", help="the minimum base quality to output")
	parser.add_argument("--count-nm", action="store_true", default = False, dest = "count_nm", help="summarize NM information for reads in the consensus read alignment file(s)")
	parser.add_argument("--output-alignment", action="store_true", default = False, dest = "output_alignment", help="output all alignment reads along with QC annotation")
	parser.add_argument("--retain-alignment", action="store_true", default = False, dest = "retain_alignment", help="retain the resulting alignment file after QC filtering")
	parser.add_argument("sample", type = str, help="sample name")
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
	
	#set default values for options if not provided
	if (options.mtdna is None):
		options.mtdna = default_arguments.get("mtdna")
	
	if (options.mtdna_offset is None):
		if (default_arguments.get("mtdna_offset")):
			options.mtdna_offset = int(default_arguments.get("mtdna_offset"))
	
	if (options.probe is None):
		options.probe = default_arguments.get("probe")
	
	if options.alignment_file:
		alignment_file = options.alignment_file
	else:
		alignment_file = [options.sample + ".mtdna.consensus.bam",]
	for i in alignment_file:
		assert os.path.exists(i), "Cannot find the consensus bam file <%s>." % i
	
	tag_excl = {}
	if (options.tag_excl is not None):
		for i in options.tag_excl.split(","):
			i = i.strip().upper()
			if (i):
				tag_excl[i] = 1
	
	tag_incl = {}
	if (options.tag_incl is not None):
		for i in options.tag_incl.split(","):
			i = i.strip().upper()
			if (i):
				tag_incl[i] = 1	
	
	#read alternate mtDNA sequences. e.g., seqeuences of numts, alignment artifacts in hsd format
	if (options.numts):
		mtdna_alt_seq = retrieve_alternate_sequences(options.numts, options.mtdna_offset)
	else:
		mtdna_alt_seq = {}
	
	if (options.count_nm):
		#summarize the consensus read information without generating the pileup file
		allele_count = summarize_alignment_file(alignment_file, tag_excl, tag_incl, mtdna_alt_seq, options.numts_fmin, options.numts_fmax)
		with open(options.outpath+os.path.sep+options.sample+".mtdna.consensus.summary", "wb") as out:
			out.write("ID\tPROBE\tFAMSIZE\tMINOR_NM\tALT_NM\n")
			for r1_probe in sorted(allele_count.keys()):
				for family_size in sorted(allele_count[r1_probe].keys()):
					minor_nm, alt_nm = allele_count[r1_probe][family_size]
					for i, j in zip(minor_nm, alt_nm):
						out.write("\t".join(map(str,[options.sample, r1_probe, family_size, i, j]))+"\n")
	elif (options.output_alignment):
		#annotate mtDNA consensus read alignments with the filters specified
		alignment_file_after_qc, alignment_reads_after_qc = annotate_alignment_file(options.sample, options.outpath, alignment_file, options.fs_min, options.fs_max, options.nm_max, options.nm_max_dloop, options.mtdna_offset, tag_excl, tag_incl, mtdna_alt_seq, options.numts_fmin, options.numts_fmax)
	else:
		#process mtDNA consensus read alignments with the filters specified
		alignment_file_after_qc, alignment_reads_after_qc = filter_alignment_file(options.sample, options.outpath, alignment_file, options.fs_min, options.fs_max, options.nm_max, options.nm_max_dloop, options.mtdna_offset, tag_excl, tag_incl, mtdna_alt_seq, options.numts_fmin, options.numts_fmax)
		
		#generate the consensus pileup file for calling variants
		generate_pileup_file(options.sample, options.outpath, options.gzip, alignment_file_after_qc, alignment_reads_after_qc, options.probe, options.mtdna, options.qual_min, options.mtdna_offset)
		
		if (not options.retain_alignment):
			#delete temporary alignment files
			execute("rm -f %s" % alignment_file_after_qc)

def main():
	run(sys.argv[0], sys.argv[1:])
	
if __name__ == "__main__":
	main()
