
Pipeline to download and process the data to be used for this workshop.

```
$SAMTOOLS=samtools
mkdir Data
mkdir Results
```

```
bash Scripts/getBams.sh
# in Data/PEL.BAMs/* and TSI and LWK
rm *.bai
```

```
ls Data/LWK.BAMs/*.bam Data/TSI.BAMs/*.bam Data/PEL.BAMs/*.bam > ALL.bamlist
ls Data/LWK.BAMs/*.bam > LWK.bamlist
ls Data/TSI.BAMs/*.bam > TSI.bamlist
ls Data/PEL.BAMs/*.bam > PEL.bamlist
```

```
wget http://dna.ku.dk/~thorfinn/hg19ancNoChr.fa.gz
mv *.fa.gz Data/.
zcat Data/hg19ancNoChr.fa.gz > Data/ancHg19.fa
$SAMTOOLS/samtools faidx Data/hg19ancNoChr.fa.gz
```

```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
zcat hs37d5.fa.gz > hs37d5.fa
bgzip hs37d5.fa
$SAMTOOLS/samtools faidx hs37d5.fa.gz
mv hs37d5* Data/.
```







