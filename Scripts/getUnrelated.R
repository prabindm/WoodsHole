
cat("Extracting unrelated IDs...\n")

fin="File/20130606_g1k.ped"

sam=read.table(fin, stringsAsFac=F, head=T, sep="\t")

#table(sam$Relationship)

#NSAMPLES=20

#pops to extract are
# LWK TSI PEL
pops=c("LWK", "TSI", "PEL")

cat("Written these files: ")
for (pop in pops) {
	pid=sam$Individual.ID[which(sam$Pop==pop & sam$Rela %in% c("mother","father","unrel","unrels"))]
	cat(pid, sep="\n", file=paste(pop, ".txt", sep="",collapse=""))
	#cat(pop, ":", length(pid), "\n")
	cat(paste(pop, ".txt", sep="",collapse=""), "")
}
cat("\n")



