# WoodsHole

Material for Workshop on Population and Speciation Genomics

27th January 2016

## Material

To download all the material in this web page use [git](http://git-scm.com/):

	git clone https://github.com/mfumagalli/EvoGen_course

`Data` and `Results` folder will be created automatically.

## Data

As an illustration, we will use 60 BAM files of human samples (of African, European, and Native American descent), a reference genome, and putative ancestral sequence.
The human data represents a small genomic region (1MB on chromosome 11) extracted from the 1000 Genomes Project data set.
More information on this project can be found [here](http://www.1000genomes.org/), including their last publication available [here](http://www.nature.com/nature/journal/v526/n7571/full/nature15393.html).

All data is publicly available.
A pipeline to retrieve such data is provided [here](https://github.com/mfumagalli/WoodsHole/blob/master/Files/data.md).
You need to have samtools installed in your /usr/bin to run this.
Data will be saved (but not pushed to git main repository) in `Data` folder.

Additional scripts are be provided in the `Scripts/` folder.


## Agenda

### Lecture:

* Basics of data handling and filtering
* SNPs and genotypes calling
* Estimation of summary statistics from low-depth data
* Paper discussion

### Practical

* Basic [filtering](https://github.com/mfumagalli/WoodsHole/blob/master/Files/filtering.md)
* Estimation of allele frequencies and [SNP calling](https://github.com/mfumagalli/WoodsHole/blob/master/Files/snpcall.md)
* [Genotype calling](https://github.com/mfumagalli/WoodsHole/blob/master/Files/genocall.md)
* Advanced methods to estimate [SFS](https://github.com/mfumagalli/WoodsHole/blob/master/Files/sfs.md)
* [Exercise](https://github.com/mfumagalli/WoodsHole/blob/master/Files/exercise.md): identification of allele frequency differentation between, with admixture assessment and quantification, from low-depth data: the case of FADS genetic variation in Native Americans

## Credits

Some materials have been borrowed (and then adapted) from [Thorfinn Korneliussen](http://scholar.google.co.uk/citations?user=-YNWF4AAAAAJ&hl=en), [Anders Albrechtsen](http://popgen.dk/albrecht/web/WelcomePage.html), [Tyler Linderoth](http://scholar.google.com/citations?user=dTuxmzkAAAAJ&hl=en), [Filipe G. Vieira](http://scholar.google.com/citations?user=gvZmPNQAAAAJ&hl=en), Dean Ousby.


