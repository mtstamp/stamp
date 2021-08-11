#############################################################################
#
#Updated by: Yiqin Wang
#
##############################################################################

import os, sys

if (__name__ == "__main__"):
	if (len(sys.argv) != 5):
		print >>sys.stderr, "Usage: %s combined_variant_file individual_variant_file_path individual_variant_file_suffix output_file" % sys.argv[0]
		print >>sys.stderr, "Example: %s all.q30.var.combined consensus.realign.recal/var.f2 .f2.q40 consensus.realign.recal/all.q30.f2_q40.var.combined" % sys.argv[0]
		exit(-1)
	
	het_file = sys.argv[1]
	var_path = sys.argv[2]
	var_suffix = sys.argv[3]
	out_file = sys.argv[4]
	
	sites = {}
	
	head = []
	#read all variant sites in het_file
	with open(het_file, "r") as fh:
		head = fh.readline()
		head = head.rstrip("\r\n").split("\t")
		for line in fh:
			line = line.rstrip("\r\n").split("\t")
			family, name, chr, pos = line[:4]
			pos = int(pos)
			if (name not in sites):
				sites[name] = {}
			#use sample name and variant position as ids
			sites[name][pos] = [family, name, chr, pos]
	
	with open(out_file, "w") as out:
		#retain head information
		out.write("\t".join(head)+"\n")
		for name in sorted(sites.keys()):
			out_pos = sites[name]
			var = {}
			#individual variant file must be present in the folder
			with open(var_path + os.path.sep + name + var_suffix + ".var", "r") as fh:
				fh.readline()
				for line in fh:
					line = line.rstrip("\r\n").split("\t")
					family, name, chr, pos = line[:4]
					pos = int(pos)
					#record line
					if (pos in out_pos):
						var[pos] = line
			#print out number of variants identified in this sample to stdout
			print name, len(out_pos), len(var)
			#order variants according to their positions
			for i in sorted(out_pos.keys()):
				line = var.get(i, out_pos[i])
				#fill blanks to match the number of items specified by the file header
				if (len(line) < len(head)):
					line = line + ["",]*(len(head)-len(line))
				out.write("\t".join(map(str,line))+"\n")
	
