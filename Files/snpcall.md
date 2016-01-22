
**WORKFLOW**:

MAPPED DATA > FILTERING > SNP CALLING

In this section, we will go through some examples on how to assign variable sites from BAM files, once the data has been filtered.
We will show how to call SNPs with different methods, and we will compare their results.

We will mainly use the program [ANGSD](http://popgen.dk/wiki/index.php/ANGSD) and the widely-used [SAMtools](http://samtools.sourceforge.net/) at some point.


### Estimating allele frequencies and calling SNPs

One of the first quantities of interest is the identification of which sites are variable, or polymorphic, in our sample.
This search can be firstly translated in the estimation of the allele frequency for each site.
In other words, at each site we want to to estimate (or count) how many copies of different alleles (two in case of biallelic SNPs) we observe in our sample (across all sequenced individuals).

ANGSD has an option to estimate **allele frequencies**:

```
$ANGSD/angsd -doMaf
...
-doMaf	0 (Calculate persite frequencies '.mafs.gz')
	1: Frequency (fixed major and minor)
	2: Frequency (fixed major unknown minor)
	4: Frequency from genotype probabilities
	8: AlleleCounts based method (known major minor)
	NB. Filedumping is supressed if value is negative
-doPost	0	(Calculate posterior prob 3xgprob)
	1: Using frequency as prior
	2: Using uniform prior
Filters:
	-minMaf  	-1.000000	(Remove sites with MAF below)
	-SNP_pval	1.000000	(Remove sites with a pvalue larger)
	-rmTriallelic	0.000000	(Remove sites with a pvalue lower)
Extras:
	-ref	(null)	(Filename for fasta reference)
	-anc	(null)	(Filename for fasta ancestral)
	-eps	0.001000 [Only used for -doMaf &8]
	-beagleProb	0 (Dump beagle style postprobs)
	-indFname	(null) (file containing individual inbreedcoeficients)
NB These frequency estimators requires major/minor -doMajorMinor
```

Therefore, the estimation of allele frequencies requires the specification of how to assign the major and minor alleles (if biallelic).
```
$ANGSD/angsd -doMajorMinor
...
	-doMajorMinor	0
	1: Infer major and minor from GL
	2: Infer major and minor from allele counts
	3: use major and minor from a file (requires -sites file.txt)
	4: Use reference allele as major (requires -ref)
	5: Use ancestral allele as major (requires -anc)
	-rmTrans: remove transitions 0
	-skipTriallelic	0
```

Finally, one needs to specify with genotype likelihood model to use.
```
$ANGSD/angsd -GL
...
	-GL=0: 
	1: SAMtools
	2: GATK
	3: SOAPsnp
	4: SYK
	5: phys
	-trim		0		(zero means no trimming)
	-tmpdir		angsd_tmpdir/	(used by SOAPsnp)
	-errors		(null)		(used by SYK)
	-minInd		0		(0 indicates no filtering)

Filedumping:
	-doGlf	0
	1: binary glf (10 log likes)	.glf.gz
	2: beagle likelihood file	.beagle.gz
	3: binary 3 times likelihood	.glf.gz
	4: text version (10 log likes)	.glf.gz
```

From these observations, our command line could be:
```
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 -sites sites.txt\
	-GL 1 -doMajorMinor 4 -doMaf 1 -skipTriallelic 1 &> /dev/null
```
where we specify:
* -GL 1: genotype likelihood model is in SAMtools
* -doMajorMinor 4: force the major allele to be the reference (the minor is inferred)
* -doMaf 1: major and minor are fixed

Which are the output files?
```
->"Results/ALL.arg"
->"Results/ALL.mafs.gz"
```
`.args` file is a summary of all options used, while `.mafs.gz` file shows the allele frequencies computed at each site.

Have a look at this file which contains estimates of **Minor Allele Frequency** (MAF) values.
```
zcat Results/ALL.mafs.gz | head
```
and you should see
```
chromo	position	major	minor	ref	knownEM	nInd
11	61005556	C	A	C	0.000001	49
11	61005584	C	A	C	0.000001	49
11	61005597	C	T	C	0.001915	49
11	61005994	A	C	A	0.000001	46
11	61005995	G	A	G	0.000001	46
11	61005996	C	A	C	0.000001	46
11	61005997	C	A	C	0.000001	46
11	61005998	T	A	T	0.000001	46
11	61005999	G	A	G	0.000000	46
```

The first and second column indicate the position of each site, then we have major and minor alleles (based on the reference sequence), the estimated allele frequency, and the number of samples with data.

To see all sites you can type:
```
less -S Results/ALL.mafs.gz
```

To generate this file we used some options in ANGSD, `-GL 1 -doMajorMinor 4 -doMaf 1`.
To summarise:
* with `-GL` you can control how to compute genotype likelihoods, which method to use, more info [here](http://popgen.dk/angsd/index.php/Genotype_likelihoods)
* with `-doMajorMinor` you can set how to define the two allelic states, more info [here](http://popgen.dk/angsd/index.php/Inferring_Major_and_Minor_alleles)
* with `-doMaf` you can choose the method to estimate allele frequencies, more info [here](http://popgen.dk/angsd/index.php/Allele_Frequency_estimation).

------------

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
 - major and minor alleles from allele counts (you need to specify -doCounts 1; not recommended in many cases);
 - frequency from known major allele but unknown minor;
 - SNPs as those having MAF=>0.01.

```
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $ANC -out Results/ALL \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 -sites sites.txt \
	-GL 2 -doMajorMinor 2 -doMaf 2 -skipTriallelic 1  \
	-minMaf 0.01 &> /dev/null
```

You can have a look at the results:
```
zcat Results/ALL.mafs.gz | head

chromo	position	major	minor	ref	unknownEM	nInd
11	61006040	G	A	G	0.010250	53
11	61007218	T	C	C	0.444894	58
11	61007710	C	G	C	0.303258	59
11	61007781	A	G	G	0.010488	59
11	61007786	G	T	T	0.025118	59
11	61007804	T	C	C	0.329408	55
11	61014040	C	A	C	0.051804	60
11	61014194	C	T	C	0.014751	59
11	61014749	A	T	A	0.010908	60
```

How many SNPs?
```
zcat Results/ALL.mafs.gz | tail -n+2 | wc -l
```

As a general guidance, `-GL 1`, `-doMaf 1/2` and `-doMajorMinor 1` should be the preferred choice when data uncertainty is high.
If interested in analyzing very low frequency SNPs, then `-doMaf 2` should be selected.
When accurate information on reference sequence or outgroup are available, one can use `-doMajorMinor` to 4 or 5.
Also, detecting variable sites based on their probability of being SNPs is generally a better choice than defining a threshold on the allele frequency.
However, various cutoffs and a dedicated filtering should be perform to assess robustenss of your called SNPs.

-------------------------

**EXERCISE**

Try varying the cutoff for SNP calling and record how many sites are predicted to be variable for each scenario.
Identify which sites are not predicted to be variable anymore with a more stringent cutoff (e.g. between a pair of scenario), and plot their allele frequencies.

```
# iterate over some cutoffs
for PV in 0.05 1e-2 1e-4 1e-6
do
        if [ $PV == 0.05 ]; then echo SNP_pval NR_SNPs; fi
        $ANGSD/angsd -P 4 -b ALL.bamlist -ref $ANC -out Results/ALL.$PV \
		-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
		-minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 -sites sites.txt \
		-GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
		-SNP_pval $PV &> /dev/null
	echo $PV `zcat Results/ALL.$PV.mafs.gz | tail -n+2 | wc -l`
done
```

A possible output is (your numbers may be different):
```
SNP_pval NR_SNPs
0.05 384
1e-2 344
1e-4 277
1e-6 241
```

Which sites differ from 0.05 and 0.01? What is their frequency?
This script will also print out the first 20 discordant sites (pK.EM is the p-value for the SNP calling test).
```
Rscript -e 'mafs1=read.table(gzfile("Results/ALL.1e-2.mafs.gz"), he=T, strings=F); mafs5=read.table(gzfile("Results/ALL.0.05.mafs.gz"), header=T, stringsAsFact=F); mafs5[!(mafs5[,2] %in% mafs1[,2]),][1:20,]; pdf(file="Results/diff_snpcall.pdf"); par(mfrow=c(1,2)); hist(as.numeric(mafs5[!(mafs5[,2] %in% mafs1[,2]),][,6]), main="Discordant SNPs", xlab="MAF (DAF)"); hist(as.numeric(mafs5[(mafs5[,2] %in% mafs1[,2]),][,6]), main="Concordant SNPs", xlab="MAF"); dev.off();'

evince Results/diff_snpcall.pdf
```

Can you draw some conclusions from these results?
Which frequencies are more difficult to estimate and therefore affect SNP calling?


### ADDITIONAL MATERIAL

The following material is provided as a pure indication, and not all command lines have been tested for compatibility with the most recent version of used programs.

#### SAMtools

We also provide command lines [here](https://github.com/mfumagalli/EvoGen_course/blob/master/Files/snpcall_samtools.md) to call SNPs using SAMtools, and to compare results with ANGSD.

----------------------------

[HOME](https://github.com/mfumagalli/WoodsHole)



