#############################################################################
#
#Updated by: Yiqin Wang
#
##############################################################################

import sys, os

#initialize probe names
probe_all = []
for i,n in [("A",12),("B",12),("C",12),("D",10)]:
	for j in range(1,n+1):
		probe_all.append(i+str(j))
probe_all.extend(["b2m","SERPINA1","emc1","wrn","axl"])
probe_all.extend(["mtDNA","nDNA"])
	
if __name__ == "__main__":
	if (len(sys.argv) != 2):
		print >>sys.stderr, "Usage: %s input" % sys.argv[0]
		print >>sys.stderr, "Example: %s consensus.realign.recal/all" % sys.argv[0]
		exit(-1)
	n = 0
	
	head = ["name",]
	#usable reads information after alignment and consensus read calling for each probe
	for i in probe_all:
		head.extend([i+".amps",i+".reads",i+".fam_size"])
	#total read information
	head.extend(["total.reads","proper.reads.1","proper.reads.2","improper.reads.1","improper.reads.2","improper.reads.3","improper.reads.4"])
	#family size distribution
	head.extend(["mtDNA.amps.s%d" % (i+1) for i in range(20)])
	head.extend(["nDNA.amps.s%d" % (i+1) for i in range(20)])
	print >>sys.stdout, "\t".join(head)
	family_num = {}
	with open("%.family_num.summary" % sys.argv[1], "r") as fh:
		for line in fh:
			line = line.rstrip("\r\n").split(" ")
			#sample id
			id = line[1]
			
			#probe id
			probe = line[2]
			
			#number of consensus reads
			amps = line[3]
			
			#number of paired-end reads
			reads = line[4]
			
			#average number of paired-end reads among total reads captured (fitting a poisson ditribution)
			fam_size = line[5]
			if (id not in family_num):
				family_num[id] = {}
			#no duplicate probe
			assert probe not in family_num[id]
			family_num[id][probe] = [amps, reads, fam_size]
	reads = {}
	with open("%s.all.reads.summary"  % sys.argv[1], "r") as fh:
		for line in fh:
			line = line.strip().split()
			id = line[1]
			if (id not in reads):
				reads[id] = [0,0,0,0,0,0,0]
			if (line[0].find("total_reads") > 0):
				reads[id][0] = int(line[2])
			elif (line[0].find("proper_reads") > 0):
				reads[id][1] = int(line[2])
				reads[id][2] = int(line[3])
			elif (line[0].find("improper_reads") > 0):
				reads[id][3] = int(line[2])
				reads[id][4] = int(line[3])
				reads[id][5] = int(line[4])
				reads[id][6] = int(line[5])
	family_size_ndna = {}
	family_size_mtdna = {}
	if (os.path.exists("%s.all.family_size.summary" % sys.argv[1])):
		with open("%s.all.family_size.summary" % sys.argv[1], "r") as fh:
			for line in fh:
				line = line.rstrip("\r\n").split(" ")
				#probe id
				probe = line[2]
				
				#number of paired-end reads
				id = line[1]
				#round to a maximum of 20
				n = min(20,int(line[3]))
				
				#number of consensus reads
				amps = int(line[4])
				
				if (probe == "mtDNA"):
					if (id not in family_size_mtdna):
						#record up to family size of 20 
						family_size_mtdna[id] = [0,]*20
					#n-1: 1 to 0 based list conversion
					family_size_mtdna[id][n-1] += amps
				elif (probe == "nDNA"):
					if (id not in family_size_ndna):
						#record up to family size of 20 
						family_size_ndna[id] = [0,]*20
					family_size_ndna[id][n-1] += amps
	for id in sorted(family_num.keys()):
		line = [id, ]
		#set 0 as default
		for probe in probe_all:
			line.extend(family_num[id].get(probe, ["0",]*3))
		line.extend(reads.get(id, [0,]*7))
		line.extend(family_size_mtdna.get(id, [0,]*20))
		line.extend(family_size_ndna.get(id, [0,]*20))
		#output
		print >>sys.stdout, "\t".join(map(str,line))

