
**WORKFLOW**:

MAPPED DATA > FILTERING

In this section, we will show how to perform a basic filtering of sites, after the reads have been mapped or aligned.

For most of the examples, we will use the program [ANGSD](http://popgen.dk/wiki/index.php/ANGSD) (Analysis of Next Generation Sequencing Data) developed by Thorfinn Korneliussen and Anders Albrechtsen at the University of Copenhagen. 
More information about its rationale and implemented methods can be found [here](http://www.ncbi.nlm.nih.gov/pubmed/25420514).

We will use 60 BAM files of human samples (of African, European, and Native American descent), a reference genome, and putative ancestral sequence.
The human data represents a small genomic region (1MB on chromosome 11) extracted from the 1000 Genomes Project data set.

## Preparation

Please set the path for all programs and data we will be using.
As an example these are my paths.
```
ANGSD=/data/data/Software/angsd
SAMTOOLS=/data/data/Software/samtools-1.3
NGSDIST=/data/Software/ngsDist
NGSTOOLS=/data/Software/ngsTools
NGSADMIX=/data/data/Software/NGSadmix/NGSadmix
FASTME=/data/data/Software/fastme-2.1.4/src/fastme
```

If you downloaded the data using the provided script, this is what you should specify.
```
REF=Data/hs37d5.fa.gz
ANC=Data/hg19ancNoChr.fa.gz
```

#### ANGSD

Here we will use ANGSD to analyse our data. To see a full list of options type
```
$ANGSD/angsd
```
and you should see something like
```
...
Overview of methods:
	-GL		Estimate genotype likelihoods
	-doCounts	Calculate various counts statistics
	-doAsso		Perform association study
	-doMaf		Estimate allele frequencies
	-doError	Estimate the type specific error rates
	-doAncError	Estimate the errorrate based on perfect fastas
	-doHWE		Est inbreedning per site
	-doGeno		Call genotypes
	-doFasta	Generate a fasta for a BAM file
	-doAbbababa	Perform an ABBA-BABA test
	-sites		Analyse specific sites (can force major/minor)
	-doSaf		Estimate the SFS and/or neutrality tests genotype calling
	-doHetPlas	Estimate hetplasmy by calculating a pooled haploid frequency

	Below are options that can be usefull
	-bam		Options relating to bam reading
	-doMajorMinor	Infer the major/minor using different approaches
	-ref/-anc	Read reference or ancestral genome
	-doSNPstat	Calculate various SNPstat
	many others

For information of specific options type: 
	./angsd METHODNAME eg 
		./angsd -GL
		./angsd -doMaf
		./angsd -doAsso etc
		./angsd sites for information about indexing -sites files
Examples:
	Estimate MAF for bam files in 'list'
		'./angsd -bam list -GL 2 -doMaf 2 -out RES -doMajorMinor 1'
```

We will see later of to perform SNP and genotype calling (and many other things) with ANGSD.

ANGSD can accept several input files, as described [here](http://popgen.dk/angsd/index.php/Input):

* BAM/CRAM
* Pileup
* Genotype likelihood/probability files
* VCF

#### Basic filtering post-mapping

Here we show how ANGSD can also perform some basic filtering of the data.
These filters are based on:

* quality and depth, see [here](http://www.popgen.dk/angsd/index.php/Filters)
* SNP quality, see [here](http://popgen.dk/angsd/index.php/SnpFilters)
* sites, see [here](http://popgen.dk/angsd/index.php/Sites)

If the input file is in BAM format, the possible options are:
```
$ANGSD/angsd -bam
...
---------------
parseArgs_bambi.cpp: bam reader:
	-bam/-b		(null)	(list of BAM/CRAM files)
	-i		(null)	(Single BAM/CRAM file)
	-r		(null)	Supply a single region in commandline (see examples below)
	-rf		(null)	Supply multiple regions in a file (see examples below)
	-remove_bads	1	Discard 'bad' reads, (flag >=256) 
	-uniqueOnly	0	Discards reads that doesnt map uniquely
	-show		0	Mimic 'samtools mpileup' also supply -ref fasta for printing reference column
	-minMapQ	0	Discard reads with mapping quality below
	-minQ		13	Discard bases with base quality below
	-trim		0	Number of based to discard at both ends of the reads
	-only_proper_pairs 1	Only use reads where the mate could be mapped
	-C		0	adjust mapQ for excessive mismatches (as SAMtools), supply -ref
	-baq		0	adjust qscores around indels (as SAMtools), supply -ref
	-checkBamHeaders 1	Exit if difference in BAM headers
	-doCheck	1	Keep going even if datafile is not suffixed with .bam/.cram
	-downSample	0.000000	Downsample to the fraction of original data
	-nReads		50	Number of reads to pop from each BAM/CRAMs
	-minChunkSize	250	Minimum size of chunk sent to analyses

Examples for region specification:
		chr:		Use entire chromosome: chr
		chr:start-	Use region from start to end of chr
		chr:-stop	Use region from beginning of chromosome: chr to stop
		chr:start-stop	Use region from start to stop from chromosome: chr
		chr:site	Use single site on chromosome: chr
```


Some basic filtering consists in removing, for instance, read with low quality and/or with multiple hits, and this can be achieved using the parameters ```-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1```.

Also, you may want to remove reads with low mapping quality and sites with low quality or covered by few reads (low depth).
Under these circumnstances, the assignment of individual genotypes and SNPs is problematics, and can lead to errors. 

However, it is necessary to know the overall distribution of per-site depth, in order to avoid filtering too many sites.
We first derive the distribution of quality scores and depth on our data set using ```-doQsDist 1 -doDepth 1```.

```
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL.qc \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 1200
```

As an illustration here, ```-maxDepth 1200``` corresponds to a per-sample average depth of 20.
This option means that all sites with depth equal or greater than 1200 will be binned together.

Have a look at the files generated:
```
...
-> Output filenames:
		->"Results/ALL.qc.arg"
		->"Results/ALL.qc.qs"
		->"Results/ALL.qc.depthSample"
		->"Results/ALL.qc.depthGlobal"
...
```
```
# counts of quality scores
less -S Results/ALL.qc.qs
# counts of per-sample depth
less -S Results/ALL.qc.depthSample 
wc -l Results/ALL.qc.depthSample # 60 Results/ALL.qc.depthSample
# counts of global depth
less -S Results/ALL.qc.depthGlobal 
```

It is convenient to compute the percentiles of these distributions (and visualize them) in order to make an informative decision on the threshold values we will use for our filtering.

```
Rscript Scripts/plotQC.R Results/ALL.qc 2> /dev/null
```
Have a look at the output files:
```
less -S Results/ALL.qc.info
evince Results/ALL.qc.pdf
``` 

Which values would you choose as sensible thresholds on quality score and global depth (min and max)?
We may also want to remove sites where half of the individual have no data. This is achieved by the -minInd option.
A possible command line would contain the following filtering:
```
...
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
...
```
which corresponds to the following scenario:

Parameter | Meaning |
--- | --- |
-minInd 30 | use only sites with data from at least N individuals |
-setMinDepth 210 | minimum total depth |
-setMaxDepth 700 | minimum total depth |

-------------------

Finally, for the rest of the analyses on SNP and genotype calling we are going to use only a small fraction of the entire data set, namely only the first 100k sites.
Typically you can achieve these by setting `-rf` option in ANGSD but since our BAM files do not have a proper header, we have to specify each site we want to analyse, and create a BED file.
This file can be generated with (knowning that our region is on chromosome 11 from 61M to 62M):
```
Rscript -e 'write.table(cbind(rep(11,100000), seq(61000001,61100000,1)), file="sites.txt", sep="\t", quote=F, row.names=F, col.names=F)'
```

---------------------

OPTIONAL

ANGSD can also compute more sophisticated metrics to filter out SNPs, as described [here](http://popgen.dk/angsd/index.php/SnpFilters), mostly based on:

* strand bias
* deviation from HWE
* quality score bias

The strand bias models are described [here](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3532123/). Some examples of strand biases, taken from the previously cited study, can be found [here](https://github.com/mfumagalli/WoodsHole/blob/master/Files/strand_bias.png).
Different models are implemented, as seen [here](https://github.com/mfumagalli/WoodsHole/blob/master/Files/strand_bias_eq.png).

These options are still under development and thus they will not be used during this session.
Furthermore, since our global sample is a mix of populations, the HWE filtering is not appropriate.
As a general guideline, a typical command line to report SNP statistics is:

```
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
	-doMaf 1 -doMajorMinor 1 -GL 1 -SNP_pval 1e-2 -hwe_pval 1 -doSnpStat 1
```

The output files are:
```
less -S Results/ALL.hwe.gz
less -S Results/ALL.snpStat.gz
```

----------------------------

OPTIONAL - DO NOT RUN

After we made our filtering choices we can extract the valid sites we will use in the forthcoming analyses.
ANGSD can analyse a set of predefined sites, which we can derived from the .mafs.gz file.
We will discuss these options later.

```
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 \
	-minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
	-doMaf 1 -doMajorMinor 1 -GL 1 &> /dev/null
```	

How many sites?
```
zcat Results/ALL.mafs.gz | wc -l
# 955533
```
This number includes the header.

Write the coordinates of these sites in a file. The coordinates are the first 2 columns.
```
zcat Results/ALL.mafs.gz | tail -n+2 | cut -f 1,2 > sites.txt
```

We can also analyse only the first 100k sites.
Also, ANGSD requires this file to be indexed.
```
zcat Results/ALL.mafs.gz | tail -n+2 | cut -f 1,2 | head -n 100000 > sites.txt
$ANGSD/angsd sites index sites.txt
```
Now we can simply tell ANGSD to analyse only these sites using the option `-sites`.

----------------------------
