
# IGNORE THIS FILE
This is only for testing and debugging.

## path to programs

ANGSD=/data/data/Software/angsd
# -> angsd version: 0.910-48-g3ecc812 (htslib: 1.2.1-69-gb79f40a)
SAMTOOLS=/data/data/Software/samtools-1.3
NGSDIST=/data/Software/ngsDist
NGSTOOLS=/data/Software/ngsTools

## Data

bash getBams.sh
# in Data/PEL.BAMs/* and TSI and LWK
rm *.bai

wget http://dna.ku.dk/~thorfinn/hg19ancNoChr.fa.gz
mv *.fa.gz Data/.
#zcat hg19ancNoChr.fa.gz > ancHg19.fa

$SAMTOOLS/samtools faidx Data/hg19ancNoChr.fa.gz

ANC=Data/hg19ancNoChr.fa.gz

# USE REF! ANC ONLY FOR SFS and after

zcat hs37d5.fa.gz > hs37d5.fa
bgzip hs37d5.fa
$SAMTOOLS/samtools faidx hs37d5.fa.gz

ls Data/LWK.BAMs/*.bam Data/TSI.BAMs/*.bam Data/PEL.BAMs/*.bam > ALL.bamlist
ls Data/LWK.BAMs/*.bam > LWK.bamlist
ls Data/TSI.BAMs/*.bam > TSI.bamlist
ls Data/PEL.BAMs/*.bam > PEL.bamlist

## Basic filtering post-mapping

mkdir Results

# -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -r 11:61000000-62000000

$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL.qc \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 1200
# 1200 is per-sample 20X

# look at files... less -S

Rscript $NGSTOOLS/scripts/plotQC.R Results/ALL.qc 2> /dev/null

less -S Results/ALL.qc.info

evince Results/ALL.qc.pdf

# mention also hwe and strand bias filters? maybe print out stats!
# see http://popgen.dk/angsd/index.php/SnpFilters

## Allele frequencies

#`-minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 250 -setMaxDepth 750 -doCounts 1`

$ANGSD/angsd -P 4 -b ALL.bamlist -ref $ANC -out Results/ALL \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
	-GL 1 -doMajorMinor 4 -doMaf 1 -skipTriallelic 1

zcat Results/ALL.mafs.gz | head

zcat Results/ALL.mafs.gz | less -S

# how many sites?
zcat Results/ALL.mafs.gz | tail -n+2 | wc -l


## SNP calling

$ANGSD/angsd -P 4 -b ALL.bamlist -ref $ANC -out Results/ALL \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
	-GL 1 -doMajorMinor 4 -doMaf 1 -skipTriallelic 1 \
	-minMaf 0.01

# how many sites?
zcat Results/ALL.mafs.gz | tail -n+2 | wc -l

# iterate over some cutoffs
for PV in 0.05 1e-2 1e-4 1e-6
do
	if [ $PV == 0.05 ]; then echo SNP_pval NR_SNPs; fi
	$ANGSD/angsd -P 4 -b ALL.bamlist -ref $ANC -out Results/ALL.$PV \
	        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
	        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
		-SNP_pval $PV &> /dev/null
	echo $PV `zcat Results/ALL.$PV.mafs.gz | tail -n+2 | wc -l`
done

# which sites differ from 0.05 and 0.01? and what is their frequency?
Rscript -e 'mafs1=read.table(gzfile("Results/ALL.1e-2.mafs.gz"), he=T, strings=F); mafs5=read.table(gzfile("Results/ALL.0.05.mafs.gz"), header=T, stringsAsFact=F); mafs5[!(mafs5[,2] %in% mafs1[,2]),][1:20,]; pdf(file="Results/diff_snpcall.pdf"); par(mfrow=c(1,2)); hist(as.numeric(mafs5[!(mafs5[,2] %in% mafs1[,2]),][,6]), main="Discordant SNPs", xlab="MAF (DAF)"); hist(as.numeric(mafs5[(mafs5[,2] %in% mafs1[,2]),][,6]), main="Concordant SNPs", xlab="MAF"); dev.off();'
evince Results/diff_snpcall.pdf

# why this bump at freq=1? this is daf! and freq=1 is considered as SNP while freq=0 not; use -doMajorMinor 1 and you won't see this

## Genotype calling

$ANGSD/angsd -doGeno

$ANGSD/angsd -doPost

# hwe no filt
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $ANC -out Results/ALL \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
	-GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
	-SNP_pval 1e-2 \
	-doGeno 3 -doPost 1 -postCutoff 0 &> /dev/null

zcat Results/ALL.geno.gz | grep -1 - | wc -l
# 0

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


## SFS

$ANGSD/angsd -doSaf

# reset the depth filtering to /3 since previously it was the global depth

# use anc instead of ref

# 1D for each pop
for POP in LWK TSI PEL
do
	echo $POP
	$ANGSD/angsd -P 4 -b $POP.bamlist -ref $ANC -anc $ANC -out Results/$POP \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 70 -setMaxDepth 235 -doCounts 1 \
	-GL 1 -doSaf 1 &> /dev/null
done

# discuss output files

$ANGSD/misc/realSFS print Results/LWK.saf.idx | less -S

for POP in LWK TSI PEL
do
	echo $POP
	$ANGSD/misc/realSFS Results/$POP.saf.idx 2> /dev/null > Results/$POP.sfs
done

# how many bins?
cat Results/LWK.sfs

# plot the 3 1d-plots
Rscript plotSFS.R Results/LWK.sfs
evince Results/LWK.sfs.pdf

Rscript plotSFS.R Results/LWK.sfs Results/TSI.sfs Results/PEL.sfs
evince Results/LWK_TSI_PEL.pdf

# Bootstrapping
$ANGSD/misc/realSFS Results/LWK.saf.idx -bootstrap 10  2> /dev/null > Results/LWK.boots.sfs

## 2D, 3D so show using single pops
$ANGSD/misc/realSFS -P 4 Results/LWK.saf.idx Results/TSI.saf.idx 2> /dev/null > Results/LWK.TSI.sfs
# mention issue of comparing the same sites!

# output file is flatten matrix (like in dadi?)

# plot 2d and 3d?
$ANGSD/misc/realSFS -P 4 Results/LWK.saf.idx Results/TSI.saf.idx Results/PEL.saf.idx 2> /dev/null > Results/LWK.TSI.PEL.sfs

## SUMMARY STATS (pi, tajima)

# in Europeans, as high DAF?

# doSaf (likes) -> SFS -> doSaf (post probs)

for POP in LWK TSI PEL
do
	echo $POP

	# compute saf posterior probabilities
	# $ANGSD/angsd -P 4 -b $POP.bamlist -ref $ANC -anc $ANC -out Results/$POP \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 70 -setMaxDepth 235 -doCounts 1 \
	-GL 1 -doSaf 1 \
	-doThetas 1 -pest Results/$POP.sfs &> /dev/null

	# summary stats per site
	$ANGSD/misc/thetaStat make_bed Results/$POP.thetas.gz &> /dev/null

	# sliding windows
	$ANGSD/misc/thetaStat do_stat Results/$POP.thetas.gz -nChr 11 -win 50000 -step 10000  -outnames Results/$POP.thetas &> /dev/null
	mv Results/$POP.theta.thetasWindow.gz.pestPG Results/$POP.thetas.win.txt

done

# note:
#- Output in the thetas.gz are the log scaled per site estimates of the thetas
#- Output in the pestPG file are the sum of the per site estimates for a region

less -S Results/PEL.thetas.win.txt

Rscript plotSS.R
evince Results/all.ss.pdf

## PCA, ngsDIST, ngsAdmixture

## dist

# hwe no filt, geno post probs
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $ANC -out Results/ALL \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
	-GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
	-SNP_pval 1e-2 \
	-doGeno 8 -doPost 1 &> /dev/null

N_SITES=`zcat Results/ALL.mafs.gz | tail -n+2 | wc -l`
echo $N_SITES

Rscript -e 'cat(paste(rep(c("LWK","TSI","PEL"),each=20), rep(1:20, 3), sep="_"), sep="\n", file="pops.label")'

$NGSDIST/ngsDist -verbose 1 -geno Results/ALL.geno.gz -probs -n_ind 60 -n_sites $N_SITES -labels pops.label -o Results/ALL.dist -n_threads 4

/data/data/Software/fastme-2.1.4/src/fastme -D 1 -i Results/ALL.dist -o Results/ALL.tree -m b -n b &> /dev/null

Rscript -e 'library(ape); library(phangorn); pdf(file="Results/ALL.tree.pdf"); plot(read.tree("Results/ALL.tree"), cex=0.5); dev.off();' &> /dev/null

evince Results/ALL.tree.pdf

# so PEL are admixed, so use ngsadmix to get unadmixed or select most unadmixed

## ngsAdmix

NGSADMIX=/data/data/Software/NGSadmix/NGSadmix

# hwe no filt, geno post probs
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $ANC -out Results/ALL \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
	-GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
	-SNP_pval 1e-2 \
	-doGlf 2 &> /dev/null

K=3
$NGSADMIX -likes Results/ALL.beagle.gz -K $K -outfiles Results/ALL.admix.K$K -P 4 -minMaf 0.02

paste ALL.bamlist Results/ALL.admix.K$K.qopt > Results/ALL.admix.K$K.txt

Rscript getUnadmixed.R 0.90

less -S Results/PEL_unadm.txt

evince Results/ALL.admix.PEL.pdf


## Estimate allele frequencies for FADS snps

# use -sites (File containing sites to keep (chr pos))

snps.txt

11 61627960
11 61631510
11 61632310
11 61641717
11 61624414
11 61597212

$ANGSD/angsd sites index snps.txt

# use ANC, and so doMajorMinor

#cat LWK.BAMs.txt TSI.BAMs.txt Results/PEL_unadm.txt > ALL_unadm.filelist

head -n 20 ALL.bamlist > LWK.bamlist
tail -n 20 ALL.bamlist > TSI.bamlist
cp Results/PEL_unadm.BAMs.txt PEL.bamlist

# change filtering?

for POP in LWK TSI PEL
do
	echo $POP
	$ANGSD/angsd -P 4 -b $POP.bamlist -ref $ANC -anc $ANC -out Results/$POP \
		-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
		-minMapQ 20 -minQ 20 -minInd 1 -setMinDepth 10 -setMaxDepth 500 -doCounts 1 \
		-GL 1 -doMajorMinor 5 -doMaf 1 -skipTriallelic 1 \
		-sites snps.txt &> /dev/null
done

zcat Results/LWK.mafs.gz Results/TSI.mafs.gz Results/PEL.mafs.gz








































