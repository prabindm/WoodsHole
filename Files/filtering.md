
**WORKFLOW**:

MAPPED DATA > SNP CALLING

In this section, we will go through some examples on how to assign variable sites from BAM files, once the data has been filtered.
We will show how to call SNPs with different methods, and we will compare their results.

We will mainly use the program [ANGSD](http://popgen.dk/wiki/index.php/ANGSD) (Analysis of Next Generation Sequencing Data) developed by Thorfinn Korneliussen and Anders Albrechtsen at the University of Copenhagen. 
More information about its rationale and implemented methods can be found [here](http://www.ncbi.nlm.nih.gov/pubmed/25420514).

We will also utilise the widely-used SAMtools at some point.

As an illustration, we will use 60 BAM files of human samples (of African, European, and Native American descent), a reference genome, and putative ancestral sequence.
The human data represents a small genomic region (1MB on chromosome 11) extracted from the 1000 Genomes Project data set.
More information on this project can be found [here](http://www.1000genomes.org/), including their last publication available [here](http://www.nature.com/nature/journal/v526/n7571/full/nature15393.html).

## Preparation

Please set the path for all programs and data we will be using.

```
ANGSD=/data/data/Software/angsd
SAMTOOLS=/data/data/Software/samtools-1.3
NGSDIST=/data/Software/ngsDist
NGSTOOLS=/data/Software/ngsTools
NGSADMIX=/data/data/Software/NGSadmix/NGSadmix
FASTME=/data/data/Software/fastme-2.1.4/src/fastme

REF=/data/data/hg19/chr11.fa
ANC=Data/hg19ancNoChr.fa.gz
```

Create a folder for results.
```
mkdir Results
```

### Estimating allele frequencies and calling SNPs

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

ANGSD can also perform some basic filtering of the data, as described [here](http://www.popgen.dk/angsd/index.php/Filters). 

-------------

## Basic filtering post-mapping

```
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $ANC -out Results/ALL.qc \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 1200
```

As an illustration here, -maxDepth 1200 corresponds to a per-sample average depth of 20. 


----------------------------


As a first step we want to estimate **allele frequencies**:

```
ngsTools/angsd/angsd -doMaf
...
-doMaf	0 (Calculate persite frequencies '.mafs.gz')
	1: Frequency (fixed major and minor)
	2: Frequency (fixed major unknown minor)
	4: Frequency from genotype probabilities
	8: AlleleCounts based method (known major minor)
	NB. Filedumping is supressed if value is negative
...
```

Therefore our command line could be:
```
ngsTools/angsd/angsd sites index input/lyca/sites.angsd.bed # to index file with sites to keep
ngsTools/angsd/angsd -b input/lyca/bams.list -GL 1 -doMajorMinor 1 -doMaf 2 -sites input/lyca/sites.angsd.bed -out output/lyca
```
Results are stored in the `output/` folder. 

Which are the output files?
```
output/lyca.arg
output/lyca.mafs.gz
```
`.args` file is a summary of all options used, while `.mafs.gz` file shows the allele frequencies computed at each site.

We can have a look at the estimated **Minor Allele Frequency** (MAF) data:
```
gunzip -c output/lyca.mafs.gz | head
```
and you should see something like
```
chromo	position	major	minor	unknownEM	nInd
ref_contig	126	T	A	0.000000	20
ref_contig	129	A	C	0.000000	20
ref_contig	131	C	A	0.000000	20
ref_contig	133	T	C	0.000000	20
ref_contig	135	A	C	0.000000	20
ref_contig	136	T	A	0.000000	20
ref_contig	138	C	A	0.000000	20
ref_contig	139	A	C	0.000000	20
ref_contig	141	T	A	0.000000	20
```
for a total of
```
gunzip -c output/lyca.mafs.gz | tail -n+2 | wc -l
   99196
```
sites.

The first and second column indicate the position of each site, then we have major and minor alleles (based on the reference sequence), the estimated allele frequency, and the number of samples with data.

To generate this file we used some options in ANGSD, `-GL 1 -doMajorMinor 1 -doMaf 2`. 
What do they mean? 

Let us check the online help:
```
ngsTools/angsd/angsd -GL
	...	
	-GL=0: 
	1: SAMtools
	2: GATK
	3: SOAPsnp
	...
```
This means that with `-GL` you can control how to compute genotype likelihoods, which method to use.

```
ngsTools/angsd/angsd -doMajorMinor
	...
	-doMajorMinor	0
	1: Infer major and minor from GL
	2: Infer major and minor from allele counts
	3: use major and minor from a file (requires -sites file.txt)
	4: Use reference allele as major (requires -ref)
	5: Use ancestral allele as major (requires -anc)
	-skipTriallelic	0
	...
```
With -doMajorMinor you can set how to define the 2 allelic states.

```
ngsTools/angsd/angsd -doMaf
	...
	-doMaf	0 (Calculate persite frequencies '.mafs.gz')
	1: Frequency (fixed major and minor)
	2: Frequency (fixed major unknown minor)
	4: Frequency from genotype probabilities
	8: AlleleCounts based method (known major minor)
	NB. Filedumping is supressed if value is negative
	...
```
With `-doMaf` you can choose the method to estimate allele frequencies.
Please note that you can even combine different methods in the same output by simply adding the `-doMaf` parameters.
Explanations of these methods can be found [here](http://popgen.dk/angsd/index.php/Allele_Frequency_estimation).

--------

We may be interested in looking at allele frequencies only for sites that are actually variable in our sample. 
Therefore we want to perform a **SNP calling**. 
There are several ways to call SNPs using ANGSD, for instance by using these options:
```
	-minMaf  	0.000000	(Remove sites with MAF below)
	-SNP_pval	1.000000	(Remove sites with a pvalue larger)
```
Therefore we can consider assigning as SNPs sites whose estimated allele frequency is above a certain threhsold (e.g. the frequency of a singleton) or whose probability of being variable is above a specified value.

As an illustration, let us call SNPs by computing:
 - genotype likelihoods using GATK method;
 - major and minor alleles from allele counts (you need to specify -doCounts 1);
 - frequency from known major allele;
 - SNPs as those having MAF>0.05.

The command line is:
```
ngsTools/angsd/angsd -b input/lyca/bams.list -sites input/lyca/sites.angsd.bed -GL 2 -doMajorMinor 2 -doMaf 1 -minMaf 0.05 -doCounts 1 -out output/lyca
```
Please note that not all the combinations of parameters are possible, since they might be in conflict or require additional steps or flags.

You can have a look at the results:
```
gunzip -c output/lyca.mafs.gz | less

 chromo  position        major   minor   knownEM nInd
 ref_contig      153     A       T       0.180146        20
 ref_contig      199     G       A       0.156008        20
 ref_contig      213     A       G       0.127374        20
 ref_contig      275     A       T       0.178517        20
 ref_contig      321     G       A       0.157150        20
 ...
```

As a general guidance, `-GL 1`, `-doMaf 1/2` and `-doMajorMinor 1` should be the preferred choice when data uncertainty is high.
For accurate and recalibrate data users can choose `-GL 2` as well.
If interested in analyzing very low frequency SNPs, then `-doMaf 2` should be selected.
When accurate information on reference sequence or outgroup are available, one can use `-doMajorMinor` to 4 or 5.
Also, detecting variable sites based on their probability of being SNPs is generally a better choice than defining a threshold on the allele frequency. 
However, various cutoffs and a dedicated filtering should be perform to assess robustenss of your called SNPs.

**EXERCISE**
Estimate allele frequencies and call SNPs using the example dataset (or your own) using ANGSD.
Try varying the cutoff for SNP calling and record how many sites are predicted to be variable for each scenario.
Identify which sites are not predicted to be variable anymore with a more stringent cutoff (e.g. between a pair of scenario), and plot their allele frequencies.
A possible **SOLUTION** for this exercise is given [here](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/solutions.txt).

**ADDITIONAL MATERIAL**
Command lines for SNP calling using SAMtools, and comparison with ANGSD, are given [here](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/snpcall_samtools.txt).





