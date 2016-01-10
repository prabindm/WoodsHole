

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








