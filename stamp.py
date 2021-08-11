##!/usr/bin/python
#############################################################################
#Package: stamp
#Coded by:  Yiqin wang
#Usage:
#
#
##############################################################################

import sys, os

prog = "stamp"
version = "0.1.2"

def print_help():
	help_msg = """
	Program: %(prog)s
	Version: %(version)s
	
	Usage: %(prog)s <command> [options]
	Command: align      generate consensus read alignments
	         pileup     summarize consensus read bases
	         wgs-align  generate mtDNA alignments (using whole-genome sequencing data)
	         wgs-depth  summarize read coverage (using whole-genome sequencing data)
	         scan       variant identification
	         annot      variant annotation
	""" % {"prog":prog, "version":version}
	
	print >>sys.stderr, help_msg

def check_packages(packages):
	#check availability of packages
	for package in packages:
		#sys.stderr.write("%s: " % package)
		try:
			test = __import__(package)
		except ImportError, e:
			print >>sys.stderr, "%s: missing!" % package
			exit(-1)

def init_run(prog, args):
	from argparse import ArgumentParser
	from tools.FileIO import read_arguments, write_arguments

	parser = ArgumentParser(prog = prog, usage="%(prog)s [-options] sample\n", version="%(prog)s v0.1.1")
	parser.add_argument("--samtools", type = str, dest = "samtools", help="full path of the samtools executable")
	parser.add_argument("--bwa", type = str, dest = "bwa", help="full path of the bwa executable")
	parser.add_argument("--bamleftalign", type = str, dest = "bamleftalign", help="full path of the bamleftalign/freebayes executable")
	parser.add_argument("--picard", type=str, dest="picard", help="full path of the picard executable")
	parser.add_argument("-p", "--probe", type = str, dest = "probe", help="full path of the file containing the STAMP probe information")
	parser.add_argument("--genome", type=str, dest="genome", help="full path of the complete reference genome sequences")
	parser.add_argument("--mtdna", type=str, dest="mtdna", help="full path of the mtDNA reference sequence")
	parser.add_argument("--mtdna-offset", type=int, dest="mtdna_offset", help="the position offset used in parsing mtDNA read alignments")
	options = parser.parse_args(args)
	
	print >>sys.stderr, "Initializing stamp ..."
	check_packages(["scipy",])
	
	arguments = []
	arguments.append(["samtools", options.samtools, parser._option_string_actions["--samtools"].help])
	arguments.append(["bwa", options.bwa, parser._option_string_actions["--bwa"].help])
	arguments.append(["bamleftalign", options.bamleftalign, parser._option_string_actions["--bamleftalign"].help])
	arguments.append(["picard", options.picard, parser._option_string_actions["--picard"].help])
	arguments.append(["probe", options.probe, parser._option_string_actions["--probe"].help])
	arguments.append(["genome", options.genome, parser._option_string_actions["--genome"].help])
	arguments.append(["mtdna", options.mtdna, parser._option_string_actions["--mtdna"].help])
	arguments.append(["mtdna_offset", options.mtdna_offset, parser._option_string_actions["--mtdna-offset"].help])
	
	#package path
	pack_path = os.path.dirname(os.path.abspath(__file__))
	arg_file = pack_path + os.path.sep + "tools" + os.path.sep + "stamp.arg"
	
	if (os.path.exists(arg_file)):
		#overwrite existing values
		default_arguments = read_arguments(arg_file)
		for i in xrange(len(arguments)):
			if (arguments[i][1] is None and arguments[i][0] in default_arguments):
				arguments[i][1] = default_arguments.get(arguments[i][0])
	
	write_arguments(arguments, arg_file)
	print >>sys.stderr, "Done!"

def main():
	if (len(sys.argv) < 2):
		print_help()
		exit(-1)
	command = sys.argv[1]
	if (command == "align"):
		import tools.align as align
		#run stamp align
		align.run("%s %s" % (prog,command), sys.argv[2:])
	elif (command == "pileup"):
		import tools.pileup as pileup
		#run stamp pileup
		pileup.run("%s %s" % (prog,command), sys.argv[2:])
	elif (command == "scan"):
		import tools.scan as scan
		#run stamp scan
		scan.run("%s %s" % (prog,command), sys.argv[2:])
	elif (command == "annot"):
		import tools.annot as annot
		#run stamp annot
		annot.run("%s %s" % (prog,command), sys.argv[2:])
	elif (command == "wgs-align"):
		import tools.wgs_align as wgs_align
		#run stamp wgs-align
		#pysam is currently required to parse alignment file in wgs-align
		check_packages(["pysam", ])
		wgs_align.run("%s %s" % (prog, command), sys.argv[2:])
	elif (command == "wgs-depth"):
		import tools.wgs_depth as wgs_depth
		#run stamp wgs-depth
		wgs_depth.run("%s %s" % (prog, command), sys.argv[2:])
	elif (command == "init"):
		#run stamp init
		init_run("%s %s" % (prog,command), sys.argv[2:])
	else:
		#unknown command; print the help page and exit.
		print_help()

if __name__ == "__main__":
	main()
