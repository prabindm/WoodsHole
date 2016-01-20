
**WORKFLOW**:

FILTERED DATA > SNP CALLING > GENOTYPE CALLING

Here we will explore several ways to call genotypes from sequencing data.
We will use ANGSD (and SAMtools as additional material), and compare results using different options.

We now see how to use ANGSD to call genotypes.
The specific option is `-doGeno`:
```
$ANGSD/angsd -doGeno
...
-doGeno	0
	1: write major and minor
	2: write the called genotype encoded as -1,0,1,2, -1=not called
	4: write the called genotype directly: eg AA,AC etc 
	8: write the posterior probability of all possible genotypes
	16: write the posterior probability of called gentype
	32: write the posterior probability of called gentype as binary
	-> A combination of the above can be choosen by summing the values, EG write 0,1,2 types with majorminor as -doGeno 3
	-postCutoff=0.333333 (Only genotype to missing if below this threshold)
	-geno_minDepth=-1	(-1 indicates no cutof)
	-geno_maxDepth=-1	(-1 indicates no cutof)
	-geno_minMM=-1.000000	(minimum fraction af major-minor bases)
	-minInd=0	(only keep sites if you call genotypes from this number of individuals)

	NB When writing the posterior the -postCutoff is not used
	NB geno_minDepth requires -doCounts
	NB geno_maxDepth requires -doCounts
```

Therefore, if we set `-doGeno 2`, genotypes are coded as 0,1,2, as the number of alternate alleles.
If we want to print the major and minor alleles as well then we set `-doGeno 3`.

To calculate the posterior probability of genotypes we need to define a model.
```
$ANGSD/angsd -doGeno

...
-doPost	0	(Calculate posterior prob 3xgprob)
	1: Using frequency as prior
	2: Using uniform prior
...
```
`-doPost 1` uses the estimate per-site allele frequency as a prior for genotype proportions, assuming Hardy Weinberg Equilibrium.
We will see later what to do when the assumption of HWE is not valid.

A typical command for genotype calling assuming HWE is:

```
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $ANC -out Results/ALL \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 -sites sites.txt\
	-GL 1 -doMajorMinor 1 -doMaf 2 -skipTriallelic 1 \
	-SNP_pval 1e-3 \
	-doGeno 3 -doPost 1 -postCutoff 0 &> /dev/null
```

Have a look at the output file:
```
less -S Results/ALL.geno.gz
```

How many sites have at least one missing genotype?
```
zcat Results/ALL.geno.gz | grep -1 - | wc -l
# 0
```

Why is that?

You can control how to set missing genotype when their confidence is low with `-postCutoff`.
For instance, we can set as missing genotypes when their (highest) genotype posterior probability is below 0.95:

```
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $ANC -out Results/ALL \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 -sites sites.txt\
	-GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
	-SNP_pval 1e-2 \
	-doGeno 3 -doPost 1 -postCutoff 0.95 &> /dev/null
```

How many sites do we have in total?
How many sites have at least one missing genotype now?
```
zcat Results/ALL.geno.gz | wc -l
# 574
zcat Results/ALL.geno.gz | grep -1 - | wc -l
# 510
```

Why are there some many sites with missing genotypes?

The mean depth per sample is around 7/8, therefore genotypes cannot be assigned with very high confidence.

Setting this threshold depends on the mean sequencing depth of your data, as well as your application.
For some analyses you need to work only with high quality genotypes (e.g. measure of proportion of shared SNPs for gene flow estimate), while for others you can be more relaxed (e.g. estimate of overall nucleotide diversity).
We will show later how to accurately estimate summary statistics with low-depth data.

If we use a uniform prior, then the command line requires `-doPost 2`:

```
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $ANC -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 -sites sites.txt\
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-2 \
        -doGeno 3 -doPost 2 -postCutoff 0.95 &> /dev/null
```

How many sites have at least one missing genotype in this case?
```
zcat Results/ALL.geno.gz | grep -1 - | wc -l
# 10904
```

Did you expect such difference compared to the case of HWE-based prior?

Optionally, if you want to investigate some of these differences more in detail, you can look at posterior probabilities and eventually to the raw data (BAM and mpileup files) using samtools.

###ADDITIONAL MATERIAL

The following material is provided as a pure indication, and not all command lines have been tested for compatibility with the most recent version of used programs.

#### Inbred species

For some studies (domesticated or self-pollinated species), it is important to consider any deviation from HWE when performing genotype or SNP calling.
We provide some command lines ([here](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/inbreeding.md)) to achieve this using [ngsF](https://github.com/fgvieira/ngsF).

#### SAMtools

We also provide command lines ([here](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/genocall_samtools.md)) to call genotypes using SAMtools, and to compare results with ANGSD.

#### BEAGLE

You can also use BEAGLE to increase accuracy of your genotype calling.
Several examples using BEAGLE to impute data are given [here](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/imputation.md).

#### FreeBayes

Freebayes is another tool for SNP and Genotype calling, available [here](https://github.com/ekg/freebayes).

#### GATK

Alternatively, one can use GATK, which runs slower and requires more steps. Here is an example to generate a VCF file:
```
# GATK does not support .gz compressed references
gunzip ref/hg19.fa.gz
# GATK requires a dictionary for the reference
# remember to change the path to the program
java -jar CreateSequenceDictionary.jar R= ref/hg19.fa O=ref/hg19.dict
# run gatk
./gatk  -R ref/hg19.fa -T UnifiedGenotyper -I bams.list -L 1 -nt 4 -o gatk.vcf
# zip the file again for safety
gzip ref/hg19.fa
```





