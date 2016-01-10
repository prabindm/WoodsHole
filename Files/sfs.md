

**WORKFLOW**:

FILTERED DATA > GENOTYPE AND SNP CALLING > POPULATION GENETICS (SFS)

Another important aspect of data analysis for population genetics is the estimate of the Site Frequency Spectrum (SFS). SFS records the proportions of sites at different allele frequencies. It can be folded or unfolded, and the latter case implies the use of an outgroup species to define the ancestral state. SFS is informative on the demography of the population or on selective events (when estimated at a local scale).

We use ANGSD to estimate SFS using on example dataset, using the methods described [here](http://www.ncbi.nlm.nih.gov/pubmed/22911679).
Details on the implementation can be found [here](http://popgen.dk/angsd/index.php/SFS_Estimation).
Briefly, from sequencing data one computes genotype likelihoods (as previously described). 
From these quantities ANGSD computes posterior probabilities of Sample Allele Frequency (SAF), for each site. 
Finally, an estimate of the SFS is computed.

Sequence data -> Genotype likelihoods -> Posterior probabilities of allele frequencies -> SFS

These steps can be accomplished in ANGSD using `-doSaf 1/2` options and the program `realSFS`.

```
$ANGSD/angsd -doSaf
...
-doSaf		0
	1: perform multisample GL estimation
	2: use an inbreeding version
	3: calculate genotype probabilities
	4: Assume genotype posteriors as input (still beta) 
	-doThetas		0 (calculate thetas)
	-underFlowProtect	0
	-fold			0 (deprecated)
	-anc			(null) (ancestral fasta)
	-noTrans		0 (remove transitions)
	-pest			(null) (prior SFS)
	-isHap			0 (is haploid beta!)
NB:
	  If -pest is supplied in addition to -doSaf then the output will then be posterior probability of the sample allelefrequency for each site
```

The SFS is typically computed for each population separately.
We need to slightly modify the filtering options as now each population has 20 samples. So now we set `-minInd 10 -setMinDepth 70 -setMaxDepth 235`.
Moreover, we want to estimate the unfolded SFS and we use a putative ancestral sequence to polarise our alleles (in ancestral and derived states).

```
for POP in LWK TSI PEL
do
        echo $POP
        $ANGSD/angsd -P 4 -b $POP.bamlist -ref $REF -anc $ANC -out Results/$POP \
        	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        	-minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 70 -setMaxDepth 235 -doCounts 1 \
        	-GL 1 -doSaf 1 &> /dev/null
done
```

Have a look at the output file.
```
$ANGSD/misc/realSFS print Results/LWK.saf.idx | less -S
```
These values represent the sample allele frequency likelihoods at each site.

The next step would be to use these likelihoods and estimate the overall SFS.
This is achieved by the program `realSFS`.

```
$ANGSD/misc/realSFS
	-> ---./realSFS------
	-> EXAMPLES FOR ESTIMATING THE (MULTI) SFS:

	-> Estimate the SFS for entire genome??
	-> ./realSFS afile.saf.idx 

	-> 1) Estimate the SFS for entire chromosome 22 ??
	-> ./realSFS afile.saf.idx -r chr22 

	-> 2) Estimate the 2d-SFS for entire chromosome 22 ??
	-> ./realSFS afile1.saf.idx  afile2.saf.idx -r chr22 

	-> 3) Estimate the SFS for the first 500megabases (this will span multiple chromosomes) ??
	-> ./realSFS afile.saf.idx -nSites 500000000 

	-> 4) Estimate the SFS around a gene ??
	-> ./realSFS afile.saf.idx -r chr2:135000000-140000000 

	-> Other options [-P nthreads -tole tolerence_for_breaking_EM -maxIter max_nr_iterations -bootstrap number_of_replications]

	-> See realSFS print for possible print options
	-> Use realSFS print_header for printing the header

	->------------------
	-> NB: Output is now counts of sites instead of log probs!!
	-> NB: You can print data with ./realSFS print afile.saf.idx !!
	-> NB: Higher order SFS can be estimated by simply supplying multiple .saf.idx files!!
	-> NB: Program uses accelerated EM, to use standard EM supply -m 0
```

This will estimate the SFS for each population

```
for POP in LWK TSI PEL
do
        echo $POP
        $ANGSD/misc/realSFS Results/$POP.saf.idx 2> /dev/null > Results/$POP.sfs
done
```

How many values do you expect?
```
cat Results/LWK.sfs
```

Let us plot the SFS for each pop using this simple R script.
```
Rscript Scripts/plotSFS.R Results/LWK.sfs
evince Results/LWK.sfs.pdf
```

Now we compare the 3 SFS.
```
Rscript Scripts/plotSFS.R Results/LWK.sfs Results/TSI.sfs Results/PEL.sfs
evince Results/LWK_TSI_PEL.pdf
```
Do they behave like expected? 
Which population has more SNPs?
Which population has a higher proportion of common (not rare) variants?

It is sometimes convenient to generate bootstrapped replicates of the SFS, by sampling with replacements genomic segments.
This could be used for instance to get confidence intervals when using the SFS for demographic inferences.
This can be achieved in ANGSD using:
```
$ANGSD/misc/realSFS Results/LWK.saf.idx -bootstrap 10  2> /dev/null > Results/LWK.boots.sfs
```
This command may take some time.
The output file has one line for each boostrapped replicate.

--------------------------

We discussed how this method does not rely on genotype calling.
How does it compare against methods that assign individual genotypes?

Let us make this test on the TSI data set.
We now estimate the SFS from called genotypes using either a HWE-based or uniform prior.
```

$ANGSD/angsd -P 4 -b $POP.bamlist -ref $REF -anc $ANC -out Results/TSI_hwe \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 70 -setMaxDepth 235 -doCounts 1 \
		-doMaf 2 -SNP_pval 0.01 \
                -GL 1 -doGeno 2 -doPost 1 -doMajorMinor 5 &> /dev/null

$ANGSD/angsd -P 4 -b $POP.bamlist -ref $REF -anc $ANC -out Results/TSI_unif \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 70 -setMaxDepth 235 -doCounts 1 \
		-doMaf 2 -SNP_pval 0.01 \
                -GL 1 -doGeno 2 -doPost 2 -doMajorMinor 5 &> /dev/null
```

From called genotypes, we can calculate the SFS by summing all observations of derived alleles across sites.
```
Rscript Scripts/getSFS.R Results/TSI_hwe.geno.gz > Results/TSI_hwe.sfs
Rscript Scripts/getSFS.R Results/TSI_unif.geno.gz > Results/TSI_unif.sfs
```
Plot the comparison without normalisation but looking at actual counts.
```
Rscript Scripts/plotSFS_abs.R Results/TSI.sfs Results/TSI_hwe.sfs Results/TSI_unif.sfs
evince Results/TSI_TSI_hwe_TSI_unif.pdf
```
Comment on what observed.

Based on these considerations, let us compute the SFS without SNP calling.

```
$ANGSD/angsd -P 4 -b $POP.bamlist -ref $REF -anc $ANC -out Results/TSI_hwe \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 70 -setMaxDepth 235 -doCounts 1 \
                -GL 1 -doGeno 2 -doPost 1 -doMajorMinor 5 -doMaf 2 &> /dev/null

$ANGSD/angsd -P 4 -b $POP.bamlist -ref $REF -anc $ANC -out Results/TSI_unif \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 70 -setMaxDepth 235 -doCounts 1 \
                -GL 1 -doGeno 2 -doPost 2 -doMajorMinor 5 -doMaf 2 &> /dev/null

Rscript Scripts/getSFS.R Results/TSI_hwe.geno.gz > Results/TSI_hwe.sfs
Rscript Scripts/getSFS.R Results/TSI_unif.geno.gz > Results/TSI_unif.sfs

Rscript Scripts/plotSFS_abs.R Results/TSI.sfs Results/TSI_hwe.sfs Results/TSI_unif.sfs
evince Results/TSI_TSI_hwe_TSI_unif.pdf
```
Comment on the result.

-----------------------

It is very useful to estimate a multi-dimensional SFS, for instance the joint SFS between 2 populations (2D).
This can be used for making inferences on their divergence process (time, migration rate and so on).

An important issue when doing this is to be sure that we are comparing the exactly same corresponding sites between populations.
ANGSD does that automatically and considers only a set of overlapping sites.
The 2D-SFS between LWK and TSI is computed with:
```
$ANGSD/misc/realSFS -P 4 Results/LWK.saf.idx Results/TSI.saf.idx 2> /dev/null > Results/LWK.TSI.sfs
```

The output file is a flatten matrix, where each value is the count of sites with the corresponding joint frequency ordered as [0.0] [0.1] and so on.
```
less -S Results/LWK.TSI.sfs
```
You can plot it, but you need to define how many samples you have per population.
```
Rscript Scripts/plot2DSFS.R Results/LWK.TSI.sfs 20 20
evince Results/LWK.TSI.sfs.pdf
```

Finally, you can even estimate SFS with higher order of magnitude:
```
$ANGSD/misc/realSFS -P 4 Results/LWK.saf.idx Results/TSI.saf.idx Results/PEL.saf.idx 2> /dev/null > Results/LWK.TSI.PEL.sfs
```

--------------------------------------------------

Under the same rationale, summary statistics and indexes of nucleotide diversity can be calculated without relying on called genotypes in ANGSD.
Briefly, expectations of such statistics are estimated from the sample allele frequency probabilities.

Here we show a typical pipeline, assuming we have already estimate the SFS for each population (see above).
The pipeline works as follow: 
-doSaf (likelihoods) -> misc/realSFS (SFS) -> -doSaf (posterior probabilities) -> -doThetas (summary statistics) -> misc/thetaStat (sliding windows)

```
for POP in LWK TSI PEL
do
        echo $POP

        # compute saf posterior probabilities
        $ANGSD/angsd -P 4 -b $POP.bamlist -ref $ANC -anc $ANC -out Results/$POP \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 70 -setMaxDepth 235 -doCounts 1 \
        -GL 1 -doSaf 1 \
        -doThetas 1 -pest Results/$POP.sfs &> /dev/null

        # summary statistics per site
        $ANGSD/misc/thetaStat make_bed Results/$POP.thetas.gz &> /dev/null

        # sliding windows
        $ANGSD/misc/thetaStat do_stat Results/$POP.thetas.gz -nChr 11 -win 50000 -step 10000  -outnames Results/$POP.thetas &> /dev/null
        mv Results/$POP.theta.thetasWindow.gz.pestPG Results/$POP.thetas.win.txt

done
```

Notes from ANGSD website:
* Output in the thetas.gz are the log scaled per site estimates of the thetas
* Output in the pestPG file are the sum of the per site estimates for a region

Have a look at output file:
```
less -S Results/PEL.thetas.win.txt
```
Columns are:
 - (indexStart,indexStop)(posStart,posStop)(regStart,regStop) chrom window_center; <br>
 - 5 estimators of theta: Watterson, pairwise, Fu & Li, Fay H, L; <br>
 - 5 neutrality test statistics: Tajima D, Fu&Li F, Fu&Li D, Fay H, Zeng E. <br>
 - The final column is the effetive number of sites with data in the window. <br>

Plot the results:
```
Rscript Scripts/plotSS.R
evince Results/all.ss.pdf
```





