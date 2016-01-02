# WoodsHole

Material for Workshop on Population and Speciation Genomics

27th January 2016

## Material

To download all the material in this web page use [git](http://git-scm.com/):

	git clone https://github.com/mfumagalli/EvoGen_course

Enter the folder and create a subfolder "Work" where you will be running the exercises:

	cd WoodsHole
	mkdir Work
	cd Work
	mkdir Results # you results will be saved here

## Data

We will use some example data sets:

 - Raw reads data of **chipmunks** (optional) are from [this](http://www.ncbi.nlm.nih.gov/pubmed/24118668) paper, and are kindly provided by Tyler Linderoth and Ke Bi.
In this paper, Authors sequenced around 4Mbp from early 20th century samples. We will use a very small subset of the original dataset.

 - BAM files of **butterflies** are a small sample subset of the original study and we will analyse around 1 Mbp.
More information about the original aligned files can be obtained at the original paper [here](http://www.ncbi.nlm.nih.gov/pubmed/22759293).
In this study, authors analysed a butterfly Next-Generation Sequencing dataset sequenced using Illumina GAII technology.
The dataset consists of a total of 381 samples of Lycaeides idas, Lycaeides melissa and an hybrid collection (Jackson Hole).
DNA resequencing was conducted on custom reduced genomic complexity libraries (RAD-seq).
We will analysed only a subset of 20 samples.

 - Other BAM files that will be used are a subset of **human** sequencing data from the 1000 Genomes Project.
More information on this project can be found [here](http://www.1000genomes.org/).
This dataset comprises 33 individuals of European descent.

 - Genotype likelihoods from inbred samples (optional) will be generated on-the-run.

## Agenda

### Lecture:

* SNP and genotype calling
* Estimation of summary statistics from NGS data
* Paper discussion

### Practical

* SNP and genotype calling using ANGSD
* Advanced methods to estimate SFS and summary statistics
* Additional material:
	+ Basic data filtering
	+ Dealing with inbred sample
	+ Population structure and admixture from low-depth data

## Credits

Some materials have been borrowed (and then adapted) from [Thorfinn Korneliussen](http://scholar.google.co.uk/citations?user=-YNWF4AAAAAJ&hl=en), [Anders Albrechtsen](http://popgen.dk/albrecht/web/WelcomePage.html), [Tyler Linderoth](http://scholar.google.com/citations?user=dTuxmzkAAAAJ&hl=en), [Filipe G. Vieira](http://scholar.google.com/citations?user=gvZmPNQAAAAJ&hl=en).




 







