
args=commandArgs(T)

fin="Results/ALL.admix.K3.txt"

th=as.numeric(args[1])

pops=c(args[2], args[3], args[4])
inds=c(which(pops=="LWK"), which(pops=="TSI"), which(pops=="PEL"))
rm(pops)

pops=c("AFR", "EUR", "NAM")
cols=c("blue","red","green")

pops=pops[inds]
cols=cols[inds]

pdf(file="Results/ALL.admix.PEL.pdf")

admix<-t(as.matrix(read.table(fin)[21:40,-1]))
barplot(admix,col=cols,space=0,border=NA,xlab="Individuals",ylab="admixture",legend=pops)

dev.off()

res=read.table(fin, stringsAsFactors=F, head=F)[21:40,]

# V1 is NAM, V2 is AFR, V3 is EUR

ii=which(res$V2>=th)
cat("Samples retained:", length(ii), "\n")

cat(res$V1[ii], sep="\n", file="Results/PEL_unadm.BAMs.txt")

cat("Output files:", "Results/ALL.admix.PEL.pdf", "Results/PEL_unadm.BAMs.txt", "\n")




