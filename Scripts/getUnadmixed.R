
th=as.numeric(commandArgs(T))

fin="Results/ALL.admix.K3.txt"


pops=c("AFR", "EUR", "NAM")
cols=c("blue","red","green")

pops=pops[c(3,1,2)]
cols=cols[c(3,1,2)]

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




