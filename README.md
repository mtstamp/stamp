# STAMP tool kit

**`stamp`** is a set of computational tools developed for processing sequencing data from **mt-STAMP** (**s**equencing by **t**argeted **a**mplification of **m**ultiplex **p**robes). **mt-STAMP** is a labor- and cost-effective sequencing method designed for assessing human mitochondrial DNA (mtDNA) variations, including mtDNA homoplasmies, mtDNA heterpolasmies and mtDNA content, in large-scale population studies. 

The **mt-STAMP** manuscript is in preparation. If you have questions about **mt-STAMP**, please contact [Zhenglong Gu](mailto:zg27@cornell.edu) (at Cornell University) for details.
If you have questions about **`stamp`** and related tools in this repository, please contact [Yiqin Wang](mailto:yw729@cornell.edu).

## Table of Contents

1. [Install and set up prerequisites](#install-and-set-up-pre-requisites)
2. [Lists of stamp arguments](#lists-of-stamp-arguments)
3. [Examples of stamp](#examples-of-stamp)
4. [Other tools](#other-tools)
5. [References](#references)
6. [Licensing](#licensing)

## Install and set up prerequisites

At the root of your project, use **`git clone`** to download the current **`stamp`** repository.
```shell
git clone https://github.com/mtstamp/stamp.git
```
**`stamp`** relies on [bwa-mem](https://github.com/lh3/bwa) [1] and bamleftalign from [freebayes](https://github.com/ekg/freebayes) [2] for aligning paired-end reads, and [samtools](https://github.com/samtools/) [3] for processing alignment file and summarizing base information. So you need to have access to the executables of these tools. 
**`stamp`** uses the reference human genome containing both nuclear DNA (genome assembly GRCh38) and mitochondrial DNA (mtDNA; Revised Cambridge Reference Sequence, rCRS) sequences for read alignment, you can download the reference human genome (GRCh38 plus rCRS and other sequences) and the pre-built index files from the ftp site of the 1000 Genomes project.
```bash
cd stamp
cd refseq
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.*

```
**`stamp`** also aligns paired-end reads to a revised rCRS with the final 120bp copied to the start to accommodate alignment of reads in the D-loop region of the circular mtDNA sequence. This revised rCRS and its pre-built index files are provided in the folder refseq.

Once the download is complete, use **`stamp init`** to check the availability of necessary python packages and set the default values for a number of commonly used arguments in **`stamp`**.

```bash
python stamp.py init \
--bwa $HOME/bin/bwa \
--samtools $HOME/bin/samtools \
--bamleftalign $HOME/bin/bamleftalign \
--probe $PWD/refseq/stamp_el_probe.txt \
--genome $PWD/refseq/GRCh38_full_analysis_set_plus_decoy_hla.fa \
--mtdna $PWD/refseq/rCRS_chrM_16449-1-16569.fa \
--mtdna-offset 120

```

The default values of these arguments can also be manually edited in the argument file under the "tools" folder of your project.
```bash
#stamp.arg
#full path of the samtools executable
samtools=/home/user/bin/samtools
#full path of the bwa executable
bwa=/home/user/bin/bwa
#full path of the bamleftalign/freebyes executable
bamleftalign=/home/user/bin/bamleftalign
#full path of the file with the probe information
probe=/home/user/stamp/refseq/stamp_el_probe.txt
#full path of the complete genome reference
genome=/home/user/stamp/refseq/GRCh38_full_analysis_set_plus_decoy_hla.fa
#full path of the mtDNA reference sequence
mtdna=/home/user/stamp/refseq/rCRS_chrM_16449-1-16569.fa
#the position offset used in parsing mtDNA read alignments
mtdna_offset=120
```

## Lists of stamp arguments

There are four functional modules in stamp, including `algin`, `pileup`, `scan`, and `annot`.
Here is the command to list all stamp modules:

```console
stamp@mito:~$ python stamp.py -h

        Usage: stamp.py <command> [options]
        Command: align                  generate the consensus read alignments
                 pileup                 summarize the consensus read bases
                 scan                   variant identification
                 annot                  variant annotation

```

You can print out the list of arguments in each module by specifying the name of the module along with the "**--help**" or "**-h**" argument.

**1. print out all arguments of **`stamp align`****
```console
stamp@mito:~$ python stamp.py align -h

usage: stamp align [-options] sample

positional arguments:
  sample                                sample name

optional arguments:
  -h, --help                            show this help message and exit
  -v, --version                         show program's version number and exit
  -r1 R1 [R1 ...], --read1 R1 [R1 ...]  fastq file(s) for read 1
  -r2 R2 [R2 ...], --read2 R2 [R2 ...]  fastq file(s) for read 2
  -p PROBE, --probe PROBE               the file with the probe information
  -o OUTPATH, --outpath OUTPATH         path where to store the temporary alignment files
 --consensus-outpath CONSENSUS_OUTPATH  path where to store the consensus alignment files (default: outpath)
  --override                            override output files (default: skip exiting output files)
  --genome GENOME                       the complete genome reference
  --mtdna MTDNA                         the mtDNA reference sequence
  --mtdna-offset MTDNA_OFFSET           the position offset used in parsing mtDNA read alignments
  --numts NUMTS                         known nuclear mitochondrial DNA segments (HSD file)
  ...
```


**2. print out all arguments of **`stamp pileup`****
```console
stamp@mito:~$ python stamp.py pileup --help

usage: stamp pileup [-options] sample

positional arguments:
  sample                                sample name

optional arguments:
  -h, --help                            show this help message and exit
  -v, --version                         show program's version number and exit
  -a ALIGNMENT_FILE [ALIGNMENT_FILE ...], --alignment ALIGNMENT_FILE [ALIGNMENT_FILE ...]
                                        the consensus read alignment file(s) output from align
  -p PROBE, --probe PROBE
                                        the file with the probe information
  -o OUTPATH, --outpath OUTPATH         path where to store the pileup file
  -z, --gzip                            compress the pileup file with gzip
  --mtdna MTDNA                         the mtDNA reference sequence
  --mtdna-offset MTDNA_OFFSET           the position offset used in parsing mtDNA read alignments
  --fs-min FS_MIN                       the minimum read family size of consensus reads
  --fs-max FS_MAX                       the maximum read family size of consensus reads
  --nm-max NM_MAX                       the maximum number of mismatches of consensus reads compared to the major mtDNA sequence in the coding region
  --nm-max-dloop NM_MAX_DLOOP           the maximum number of mismatches of consensus reads compared to the major  mtDNA sequence in the Dloop region
  --tag-excl TAG_EXCL                   exclude consensus reads with the tags specified
  --tag-incl TAG_INCL                   include consensus reads with the tags specified
  --numts-excl NUMTS                    exclude reads from nuclear mitochondrial DNA segments specified by the HSD file
   ...
```


**3. print out all arguments of **`stamp scan`****

```console
stamp@mito:~$ python stamp.py scan --help

usage: stamp scan [-options] input output

positional arguments:
  input                                 a mpileup file or a batch file
  output                                the prefix of output files

optional arguments:
  -h, --help                            show this help message and exit
  -v, --version                         show program's version number and exit
  --mind MIN_DEPTH                      the minimum read depth to call variants
  --mind-fwd MIN_DEPTH_FWD              the minimum read depth on the forward strand to call variants
  --mind-rev MIN_DEPTH_REV              the minimum read depth on the reverse strand to call variants
  --minh MIN_MINOR_DEPTH                the minimum number of minor alleles to call heteroplasmies
  --minh-fwd MIN_MINOR_DEPTH_FWD        the minimum number of minor alleles on the forward strand to call heteroplasmies
  --minh-rev MIN_MINOR_DEPTH_REV        the minimum number of minor alleles on the reverse strand to call heteroplasmies
  --min-het MIN_HET_FREQ                the minimum minor allele fraction of heteroplasmies
  --min-qual MIN_QUAL                   the minimum base quality
  --max-qual MAX_QUAL                   the maximum base quality
  --min-qual-rate MIN_QUAL_RATE         the minimum proportion of bases passing the quality filter(s)
  --mle                                 use maximum likelihood estimation to compute variant quality for heteroplasmies
  --family FAMILY                       family name
  --sample SAMPLE [SAMPLE ...]          sample columns to process in the mpileup file
  --name NAME [NAME ...]                the corresponding names of the samples to process
  --batch                               proceed in the batch mode 
  --all-sites                           output information for all sites instead of only variant sites
  ...
```


**4. print out all arguments of **`stamp annot`****

```console
stamp@mito:~$ python stamp.py annot --help

usage: stamp annot [-options] input output

positional arguments:
  input                                 the variant file output from scan
  output                                the prefix of output files

optional arguments:
  -h, --help                            show this help message and exit
  -v, --version                         show program's version number and exit
  --exclude EXCLUDE                     exclude mtDNA sites from analysis
  --remove REMOVE                       remove families from analysis
  --keep KEEP                           keep only the families for analysis
  --depth DEPTH                         the minimum read depth of all variants
  --depth-min DEPTH_MIN                 the minimum read depth of heteroplasmies
  --hq-min HQ_MIN                       the minimum rate of high-quality reads of heteroplasmies
  --llr-min LLR_MIN                     the minimum quality score of heteroplasmies
  --sbias-min SBIAS_MIN                 the minimum P value for strand bias analysis of heteroplasmies
  --frac-min FRAC_MIN                   the minimum minor allele fraction of heteroplasmies
  --dev-frac-min DEV_FRAC_MIN           the minimum variant allele fraction of homoplasmies
  --annotate ANNOTATE                   annotate variants according to the file specified
  --output-ped                          output the variants detected to a ped file
  --output-hsd                          output the major alleles to a hsd file
  --output-minor-hsd                    output the minor alleles to a hsd file
  ...
```

## Examples of stamp

**1. alignment of paired-end reads using `stamp align`**

**Example:** perform paired-end read alignment and call consensus reads for sample1
```bash
python stamp.py align -r1 sample1.R1.fastq -r2 sample1.R2.fastq \
--genome GRCh38_full_analysis_set_plus_decoy_hla.fa \
--mtdna rCRS_chrM_16449-1-16569.fa --mtdna-offset 120 \
--probe stamp_el_probe.txt --numts numts_sites_hg38.rCRS.hsd \
--outpath bam --consensus-outpath consensus.realign.recal/bam sample1 

```

**Important arguments:**

- **-r1**, **-r2**: the R1 and R2 files containing paired-end reads; **`stamp align`** can read and combine multiple R1 and R2 files of a sample.
- **--genome**: the complete human genome sequence in FASTA format with bwa index.
- **--mtdna**: the mtDNA sequence in FASTA format with bwa index.
- **--mtdna-offset**: the default offset in relation to the rCRS is 120bp.
- **--numts**: the file with a collection of known NUMTS sequences in the reference genome in HSD format (http://haplogrep.uibk.ac.at/blog/tag/hsd/).
- **--probe**: the file containing the EL probe information.
- **--output**: the output directory for the alignment files of the paired-end reads.

**Output files:**
- **bam/sample1_R(1/2).fastq.gz**: the fastq files with barcode and probe information trimmed and stored as an annotation in read description (i.e. XM:Z:1:barcode:probe1:probe2).
- **consensus.realign.recal/bam/sample1.mtdna.sorted.realign.recal.bam**: the alignment file of paired-end reads generated after base recalibration.  
- **consensus.realign.recal/bam/sample1.mtdna.consensus.bam**: the alignment file of consensus reads. Consensus reads are output as single-end reads constructed from paired-end reads with the same barcode (read family). Read names are assigned using that of the first paired-end read in a read family. The number of paired-end reads and the barcode is recorded in the SAM TAGs "XF" and "XM" of the consensus read.

***
**2. Post-alignment processing and base summarization using **`stamp pileup`****

**Example:** exclude consensus reads with excessive mismatches or NUMTS annotations
```bash
python stamp.py pileup -a sample1.mtdna.consensus.bam \
--mtdna rCRS_chrM_16449-1-16569.fa --mtdna-offset 120 --probe stamp_el_probe.txt \
--numts-excl numts_sites_hg38.rCRS.hsd --nm-max 5 --nm-max-dloop 8 \
--tag-excl NUMTS,EXMISMATCH -o consensus.realign.recal/var sample1 
```

**Example:** only retain consensus reads constructed from duplicate paired-end reads
```bash
python stamp.py pileup -a sample1.mtdna.consensus.bam \
--mtdna rCRS_chrM_16449-1-16569.fa --mtdna-offset 120 --probe stamp_el_probe.txt \
--numts-excl numts_sites_hg38.rCRS.hsd  --nm-max 5 --nm-max-dloop 8 --fs-min 2 \
--tag-excl NUMTS,EXMISMATCH -o consensus.realign.recal/var.f2 sample1
```
**Important arguments:**
- **-a**: the consensus read alignment file(s) that are output from align; **`stamp pileup`** can read and combine multiple alignment files of a sample.
- **-o**: the file path where to store the pileup file.
- **--mtdna**: the mtDNA sequence in FASTA format with bwa index.
- **--mtdna-offset**: the default offset in relation to the rCRS is 120bp.
- **--numts-excl**: the file with a collection of NUMTS sequences in HSD format.
- **--probe**: the file containing the EL probe information.
- **--tag-excl**/**--tag-incl**: filter out or retain consensus reads using quality information in SAM TAG "XQ"
**Output file:**
- **consensus.realign.recal/var/sample1.mtdna.consensus.adj.pileup**: the resulting pileup file, which is compressed if "**-z**" is turned on.
- **consensus.realign.recal/var/sample1.coverage**: a table recording consensus read coverage at each site. Data include: (1) chromosome; (2) position; (3) original position in rCRS; (4) the reference allele; (5) total depth of consensus reads; (6) the number of reads with BAQ between 0 and 10; (7) the number of reads with BAQ between 10 and 20; (8) the number of reads with BAQ between 20 and 30; (9) the number of reads with BAQ between 30 and 40; (10) the number of reads with BAQ over 40; (11) total read number; (12)  the relative distance to the nearest ligation probe; (13) the relative distance to the nearest extension probe; (14) location on the ligation probe starting from the 3' end; (15) location on the extension probe starting from the 3' end; (16) probe information.

***
**3. mtDNA variant detection using **`stamp scan`****

**Example:** call mtDNA homoplasmies and heteroplasmies which had a minor allele fraction >=1% with at least 5 minor allele and a total of 500 reads with BAQ>=30
```bash
python stamp.py scan --name sample1 \
--mind 100 --mindh 500 --min-qual 30 --minh 5 --min-het 0.01 --mle \
consensus.realign.recal/var/sample1.mtdna.consensus.adj.pileup.gz consensus.realign.recal/var/sample1.q30

```

**Important arguments:**
- **input:** the input mpileup/pileup file; when "**--batch**" is turned on, **`stamp scan`** will read the arguments and the path of the pileup file from each line in the batch file provided.
- **--sample:** indicates the order of the column(s) to process in the mpileup file.
- **--name:** the name(s) of the sample(s) in the mpileup file.
- **--family:** a family identifier of the sample(s).

**Output file:**
- **consensus.realign.recal/var/sample1.q30.var:** a text file (tab-separated) which contains information on the variants identified; it contains information for all mtDNA sites if "**--all-sites**" is turned on. Data include: (1) family id; (2) sample id; (3) chromosome; (4) position in rCRS; (5) reference allele (rCRS); (6) total depth of unique (consensus) reads; (7) depth of reads on the forward strand; (8) depth of reads on the reverse strand; (9) major allele; (10-13) numbers of A, T, C and G on the forward strand passing quality filtering; (14-17) numbers of A, T, C and G on the reverse strand passing quality filtering; (18) binary status of heteroplasmy; (19) binary status of homoplasmy; (20) minor allele; (21) minor allele fraction; (22) minor allele fraction computed using maximum likelihood estimation; (23) minor allele fraction log likelihood ratio; (24-25) lower and upper bounds of the 95% confidence interval of the minor allele fraction; (26) Fisher's exact test on the minor allele count; (27) strand difference in allele fractions (Fisher's exact test P value).


***
**4. mtDNA variant annotation using **`stamp annot`****

**Example:** mtDNA variants for all samples
```bash
python stamp.py annot --exclude low_quality_sites.txt --remove low_quality_samples.txt \
--frac-min 0.01 --output-ped --output-hsd \
--annotate refseq/mtDNA_mutation_all.annovar.tsv all.q30.var.combined all.q30.var
```

**Important arguments:**
- **input:** the var file from stamp scan, or the file obtained by merging multiple var files.
- **--exclude**; **--remove**; **--keep**: the text file containing a list of mtDNA sites or sample IDs to be excluded or included.
- **--output-ped**; **--output-hsd**; **--output-hsd-minor**: convert information in the var file to a ped or hsd file.
- **--annotate**: a text file with annotation information for mtDNA. An example is provided under the folder refseq.

**Output files:**
- **all.q30.var.qc.annot:** a text file (tab-separated) contains information on the variants passing all quality filters and the related mtDNA annotation. The status column in this file indicates the type of the variants: "homoplasmy", "heteroplasmy", and "possible heteroplasmy". "Possible heteroplasmy" indicates a variant which has VAF great than the provided cutoff but does not pass other quality control filters. 
- **all.q30.var.qc.hsd:** the hsd file recording the major allele information different from the reference mtDNA sequence.
- **all.q30.var.qc.{tped/map/fam}:** the plink files recording the major allele information different from the reference mtDNA sequence.


## Other tools
Several python or R scripts (only with limited comments now) used for processing results generated from **`stamp`** are provied in the folder tools.
Some functionality of those scripts will be incorperated into **`stamp`** in the future.

**1. extract and combine alignment information in the ".summary" files output from `stamp align`**

```bash
#bam/*.summary is a summary of the results of paired-end read demultiplexing,
#including the number of read pairs included (proper reads) and excluded (improper reads).
grep _reads bam/*.summary >consensus.realign.recal/all.reads.summary

#${align}/bam/*.summary contains information on read alignement and consensus read calling
#extract the lines with a "family_num" tag in the summary file \
#which recorde the number of paired-end reads and consensus reads from each probe of STAMP
grep family_num consensus.realign.recal/bam/*.summary \
>consensus.realign.recal/all.family_num.summary

#extract the lines with a "family_size" tag in the summary file \
#which record the number of paried-end reads used to construct consensus reads for each probe of STAMP
grep family_size consensus.realign.recal/bam/*.summary|grep DNA \
>consensus.realign.recal/all.family_size.summary

#combine information extracted above into a table with each row representing one sample and \
#each column representing one type of information extracted
python tools/summarize_reads.py consensus.realign.recal/all >consensus.realign.recal/all.probe.combined
```

**2. extract site covarage information in the ".coverage" file output from `stamp pileup`**

```bash
#summarize coverage of consensus reads for each probe of STAMP after quality filtering
python tools/summarize_coverage_probe.py consensus.realign.recal/var/*.coverage \
>consensus.realign.recal/all.coverage.summary

#summarize coverage of only consensus reads created with duplicate paried-end reads\
#for each probe of STAMP after quality filtering
python tools/summarize_coverage_probe.py ${batch}/${align}/var.f2/*.coverage \
>consensus.realign.recal/all.f2.coverage.summary
```

**3. compare read depths and allele fractions at corresonding sites of two samples or in two ".var" files** 

```bash
#extract all possible mtDNA variants from the ".var" file output from `stamp scan`
#retain head information
head -n 1 consensus.realign.recal/var/*.q30.var | tail -n 1 \
>consensus.realign.recal/var/all.q30.var.combined
#the 18th and 19th columns in the ".var" file record a binary variable indicating \
#whether a variant is detected (homoplasmy or heteroplasmy) before quality filtering
#combine information from all samples (all ".var" files under this folder)
gawk '$18 == 1 || $19 == 1 {print $0}' consensus.realign.recal/var/*.q30.var \
>consensus.realign.recal/all.q30.var.combined

#extract information in corresponding samples and sites in ".var" files under the "var.f1" folder
#according to the variants recorded in all.q30.var.combined 
python tools/summarize_variation.py consensus.realign.recal/all.q30.var.combined \
consensus.realign.recal/var.f1 .f1.q23 consensus.realign.recal/all.q30.f1_q23.var.combined
#extract information in corresponding samples and sites in ".var" files under the "var.f2" folder
#according to the variants recorded in all.q30.var.combined   
python tools/summarize_variation.py consensus.realign.recal/all.q30.var.combined \
consensus.realign.recal/var.f2 .f2.q40 consensus.realign.recal/all.q30.f2_q40.var.combined

#compare allele fractions of corresponding sites in two ".var" files
#i.e. compare allele fractions between consensus reads with and without duplications
#stored in the files under the "var.f2" folder and the "var.f1" folder respectively
#consensus.realign.recal/all.q30.f1andf2.var.combined can also be read by **`stamp annot`** for annotation
Rscript tools/summarize_variation_freq.R consensus.realign.recal/all.q30.var.combined \
consensus.realign.recal/all.q30.f1_q23.var.combined \
consensus.realign.recal/all.q30.f2_q40.var.combined \
consensus.realign.recal/all.q30.f1andf2.var.combined
```

## References
1) Li, H. & Durbin, R. Fast and accurate long-read alignment with Burrows-Wheeler transform. Bioinformatics **26**, 589-95 (2010).
2) Garrison, E. & Marth, G. Haplotype-based variant detection from short-read sequencing. arXiv preprint arXiv:1207.3907 (2012).
3) Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics **25**, 2078-9 (2009).

## Licensing

This project is licensed under the terms of the MIT license.