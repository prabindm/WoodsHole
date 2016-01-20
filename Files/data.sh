
## Pipeline to download and process the data to be used for this workshop.

# set path
SAMTOOLS=/data/data/Software/samtools-1.3
echo Is this your path to samtools? $SAMTOOLS

# get unrelated samples iDs
Rscript Scripts/getUnrelated.R

# download BAM files
bash Scripts/getBams.sh
# this creates files and folders in Data/PEL.BAMs/* and TSI and LWK

# create file with list of BAMs
ls Data/LWK.BAMs/*.bam Data/TSI.BAMs/*.bam Data/PEL.BAMs/*.bam > ALL.bamlist
ls Data/LWK.BAMs/*.bam > LWK.bamlist
ls Data/TSI.BAMs/*.bam > TSI.bamlist
ls Data/PEL.BAMs/*.bam > PEL.bamlist

# download ancestral sequence
echo Downloading and processing ancestral sequence...
wget http://dna.ku.dk/~thorfinn/hg19ancNoChr.fa.gz &>/dev/null
mv *.fa.gz Data/.
# zcat Data/hg19ancNoChr.fa.gz > Data/ancHg19.fa
$SAMTOOLS/samtools faidx Data/hg19ancNoChr.fa.gz

# download reference sequence
echo Downloading and processing reference sequence...
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz &>/dev/null
gunzip hs37d5.fa.gz &>/dev/null
echo Bgzipping... this may take some time...
bgzip hs37d5.fa
$SAMTOOLS/samtools faidx hs37d5.fa.gz
mv hs37d5* Data/.

echo Done!
ls -lh Data/* > Data/download.log
echo Open Data/download.log to see what files have been generated.



