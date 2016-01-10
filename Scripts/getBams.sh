
echo Retrieving file names...
NS=20
# get IDs 
for POP in LWK TSI PEL;
do

	INPUT=$POP.txt
	INDS=`cat $INPUT`
	OUTPUT=$POP.BAMs.txt
	> $OUTPUT
	> tmp

	for i in $INDS;
	do
		grep $i Files/phase3_bamlist.txt >> tmp
	done
	# wc -l tmp
	head -n $NS tmp > $OUTPUT 
	rm tmp
done

echo Downloading BAM files...
# download and index bams
#mkdir Data
for POP in LWK TSI PEL;
do
	mkdir Data/$POP.BAMs
	echo $POP
	INDLIST=`cat $POP.BAMs.txt`
	for i in $INDLIST;
	do
		NAME=`echo -n $i | tail -c 58`
		echo $NAME
		samtools view -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/$i 11:61000000-62000000 > Data/$POP.BAMs$NAME 2> /dev/null
		#samtools index $POP.BAMs/$NAME
	done
done

echo Removing index files...
rm *.bai

exit



