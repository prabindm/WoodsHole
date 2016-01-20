
## CASE STUDY

*MOTIVATION*

Detecting signatures of natural selection in the genome has the twofold meaning of (i) understanding which adaptive processes shaped genetic variation and (ii) identifying putative functional variants.
In case of humans, biological pathways enriched with selection signatures include pigmentation, immune-system regulation and metabolic processes. 
The latter may be related to human adaptation to different diet regimes, depending on local food availability (e.g. the case of lactase persistence in dairy-practicing populations).

In this regard, fatty acid desatureses gene family (FADS) include three different genes (1,2,3) located on chromosome 11.
They are important players in regulation of metabolite levels, and are responsible for processing of omega-3 (e.g. from fish oil) and omega-6 (e.g. from land mammals) polyunsaturated fatty acids.
FADS gene have been found to be under positive selection in several human populations, including [Europeans](http://www.ncbi.nlm.nih.gov/pubmed/26595274) and Greenlandic [Inuit](http://www.ncbi.nlm.nih.gov/pubmed/26383953).
Interestingly, selection seems to act on different genes and variants within FADS in these populations.

Selection acting on FADS genes in Inuit was estimated to start 20k years ago, possibly predating the split between the ancestors of Inuit and Native Americans.

*HYPOTHESIS*

Is positive selection on FADS genes targeting Native American populations as well?

If so, is selection acting on the same variants as in Inuit?

*PLAN OF ACTION*

- Retrieve genomic data for American population (low-depth)
- Check population structure of American samples relating to Europeans and Africans
- Select individuals with high Native American ancestry
- Compute allele frequencies for SNPs of interest in Native Americans, Europeans and Africans
- Compare their allele frequencies as a hint of selection

---------------

1) Investigate population structure (or clustering) of PEL samples (Peruvians), EUR (Europeans) and LWK (Africans).

One solution would be to perform a PCA, MDS or some clustering based on genetic distances among samples.
Then we can check whether some PEL fall within EUR/LWK or all PEL form a separate clade. 

First, compute genotype posterior probabilities for all samples.
 
```
# Assuming HWE, without filtering based on probabilities, with SNP calling
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGeno 8 -doPost 1 &> /dev/null
```

Record how many sites we retrieve.
```
N_SITES=`zcat Results/ALL.mafs.gz | tail -n+2 | wc -l`
echo $N_SITES
```

Create a file with labels indicating the population of interest for each sample.
```
Rscript -e 'cat(paste(rep(c("LWK","TSI","PEL"),each=20), rep(1:20, 3), sep="_"), sep="\n", file="pops.label")'
cat pops.label
```

With [ngsDist](https://github.com/fgvieira/ngsDist) we can computer pairwise genetic distances without relying on individual genotype calls.
```
$NGSDIST/ngsDist -verbose 1 -geno Results/ALL.geno.gz -probs -n_ind 60 -n_sites $N_SITES -labels pops.label -o Results/ALL.dist -n_threads 4
less -S Results/ALL.dist
```

We can visualise the pairwise genetic distances in form of a tree.
```
$FASTME -D 1 -i Results/ALL.dist -o Results/ALL.tree -m b -n b &> /dev/null
cat Results/ALL.tree
```

Plot the tree.
```
Rscript -e 'library(ape); library(phangorn); pdf(file="Results/ALL.tree.pdf"); plot(read.tree("Results/ALL.tree"), cex=0.5); dev.off();' &> /dev/null
evince Results/ALL.tree.pdf
```

Therefore, PEL samples appear very close to EUR.
The next step would be to compute admixture proportions in order to select a subset of putative Native American (unadmixed) samples.

------------------------

2) Compute admixture proportions across samples.

We use ngsAdmix, which again works on genotype probabilities and not on individual calls.
This is suitable for low-depth data.

ngsAdmix requires genotype likelihoods in BEAGLE format as input.
We can compute these quantities with ANGSD with `-doGlf 2`.
```
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGlf 2 &> /dev/null
```

We assume 3 ancestral populations (Europeans, Africans and Native Americans) making up the genetic diversity of our samples.
Therefore we compute admixture proportions with 3 ancestral components.
```
K=3
$NGSADMIX -likes Results/ALL.beagle.gz -K $K -outfiles Results/ALL.admix.K$K -P 4 -minMaf 0.02
```

Combine samples IDs with admixture proportions and inspect the results.
```
paste ALL.bamlist Results/ALL.admix.K$K.qopt > Results/ALL.admix.K$K.txt
less -S Results/ALL.admix.K$K.txt
```

From these quantities we can extract how many samples (and which ones) have a high proportion of Native American ancestry (e.g. >0.90).
We can also plot the individual ancestral proportions for PEL samples.
We need to specify the order of ancestral populations in the admixture file `Results/ALL.admix.K$K.txt`. 
In my case these are PEL TSI LWK
``` 
Rscript Scripts/getUnadmixed.R 0.90 PEL TSI LWK
```
Inspect the results.
```
less -S Results/PEL_unadm.BAMs.txt
evince Results/ALL.admix.PEL.pdf
```

Now we have a subset of putative Native American samples.
We can calculate allele frequencies for only these samples.

------------------------

3) Estimate allele frequencies for SNPs in FADS genes of interest

In ANGSD we can restrict our analyses on a subset of positions of interest using the `-sites` option.
The positions we are looking at are the one found under selection in Inuit, shown [here]((https://github.com/mfumagalli/WoodsHole/blob/master/Files/snps_inuit.png)): 
- 11 61627960 <br>
- 11 61631510 <br>
- 11 61632310 <br>
- 11 61641717 <br>
- 11 61624414 <br>
- 11 61597212 <br>

The file with these positions need to be formatted as (chromosome positions).
```
> snps.txt
echo 11 61627960 >> snps.txt
echo 11 61631510 >> snps.txt
echo 11 61632310 >> snps.txt
echo 11 61641717 >> snps.txt
echo 11 61624414 >> snps.txt
echo 11 61597212 >> snps.txt
```
Inspect the file.
```
cat snps.txt
```

We need to index this file in order for ANGSD to process it.
```
$ANGSD/angsd sites index snps.txt
```

We are interested in calculating the derived allele frequencies, so are using the ancestral sequence to polarise the alleles.
Create new lists of BAM files.
```
head -n 20 ALL.bamlist > LWK.sub.bamlist
tail -n 20 ALL.bamlist > TSI.sub.bamlist
cp Results/PEL_unadm.BAMs.txt PEL.sub.bamlist
```

Run ANGSD to compute allele frequencies.
Here we change the filtering (more relaxed) since we are interested in outputting all sites.
```
for POP in LWK TSI PEL
do
        echo $POP
        $ANGSD/angsd -P 4 -b $POP.sub.bamlist -ref $REF -anc $ANC -out Results/$POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 1 -setMinDepth 10 -setMaxDepth 500 -doCounts 1 \
                -GL 1 -doMajorMinor 5 -doMaf 1 -skipTriallelic 1 \
                -sites snps.txt &> /dev/null
done
```

Inspect the results.
```
zcat Results/LWK.mafs.gz Results/TSI.mafs.gz Results/PEL.mafs.gz
```

Do you see any allele frequency differentiation?

----------------------------

[HOME](https://github.com/mfumagalli/WoodsHole)









