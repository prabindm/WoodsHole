
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

# hwe no filt
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $ANC -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-2 \
        -doGeno 3 -doPost 1 -postCutoff 0 &> /dev/null

zcat Results/ALL.geno.gz | grep -1 - | wc -l
# 0

enotypes are coded as 0,1,2, as the number of alternative alleles.

# hwe filt
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $ANC -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-2 \
        -doGeno 3 -doPost 1 -postCutoff 0.95 &> /dev/null

zcat Results/ALL.geno.gz | wc -l
# 12219
zcat Results/ALL.geno.gz | grep -1 - | wc -l
# 7171

# unif filt
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $ANC -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-2 \
        -doGeno 3 -doPost 2 -postCutoff 0.95 &> /dev/null

zcat Results/ALL.geno.gz | grep -1 - | wc -l
# 10904

# investigate some differences? diff file1 file2? then look at postprobs and then eventually to raw data (bam to mpileup) with samtools?



Genotypes are coded as 0,1,2, as the number of alternative alleles. 

If we want to print the major and minor allele then we set `-doGeno 3`:

Values coded as -1 represent missing data. Indeed, genotypes with a posterior probability less than a specified threshold are set as missing. 
You can vary this cutoff by setting the option -postCutoff.
For instance, if we set `-postCutoff 0`:
we get
sites with missing genotypes, while if we are more stringent and use 0.95 as cutoff
we get
sites with at least one individual with missing genotype (all!).
Why is that?

Indeed, the mean depth per sample is around 4, therefore genotypes cannot be assigned with very high confidence.

Setting this threshold depends on the mean sequencing depth of your data, as well as your application. 
For some analyses you need to work only with high quality genotypes (e.g. measure of proportion of shared SNPs for gene flow estimate), while for others you can be more relaxed (e.g. estimate of overall nucleotide diversity). 
We will show later how to accurately estimate summary statistics with low-depth data.


### Inbred species

**ADDITIONAL MATERIAL**
For some studies (domesticated or self-pollinated species), it is important to consider any deviation from HWE when performing genotype or SNP calling.
We provide some command lines ([here](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/inbreeding.md)) to achieve this using ngsTools.


### SAMtools

**ADDITIONAL MATERIAL**
We also provide command lines ([here](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/genocall_samtools.md)) to call genotypes using SAMtools, and to compare results with ANGSD.


### BEAGLE

**ADDITIONAL MATERIAL**
You can also use BEAGLE to increase accuracy of your genotype calling.
Several examples using BEAGLE to impute data are given [here](https://github.com/mfumagalli/EvoGen_course/tree/master/Files/imputation.md).

### FreeBayes

Freebayes is another tool for SNP and Genotype calling, available [here](https://github.com/ekg/freebayes).
It is especially suitable for indels and CNVs detection.

### GATK

Alternatively, one can use GATK, which runs slower and requires more steps. Here are commands to generate a VCF file from the previous examples:
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

Look at the generated VCF file (gatk.vcf). How many SNPs does the program predict?





