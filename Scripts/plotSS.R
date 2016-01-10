
#

win=5e4

fout="Results/all.ss.pdf"

fins=c("Results/LWK.thetas.win.txt", "Results/TSI.thetas.win.txt", "Results/PEL.thetas.win.txt")

theta1 <- read.table(fins[1], fill=TRUE)
theta2 <- read.table(fins[2], fill=TRUE)
theta3 <- read.table(fins[3], fill=TRUE)

xl=c(61e6, 62e6)

pdf(file=fout)

par(mfrow=c(2,2))

plot(x=theta1$V3, y=theta1$V4, pch=16, ty="l", col="blue", ylim=c(0,max(c(theta1$V4, theta2$V4, theta3$V4), na.rm=T)), xlab="Chromosome 11", ylab="Theta (Watt.)", xlim=xl)
points(x=theta2$V3, y=theta2$V4, pch=16, col="black", ty="l")
points(x=theta3$V3, y=theta3$V4, pch=16, col="green", ty="l")
legend("topleft", legend=c("LWK", "TSI", "PEL"), col=c("blue","black","green"), pch=16)

plot(x=theta1$V3, y=theta1$V5, pch=16, col="blue", ylim=c(0,max(c(theta1$V5, theta2$V5, theta3$V5), na.rm=T)), xlab="Chromosome 11", ylab="Theta (Pi)", xlim=xl, ty="l")
points(x=theta2$V3, y=theta2$V5, pch=16, col="black", ty="l")
points(x=theta3$V3, y=theta3$V5, pch=16, col="green", ty="l")
legend("topleft", legend=c("LWK", "TSI", "PEL"), col=c("blue","black","green"), pch=16)

plot(x=theta1$V3, y=theta1$V9, pch=16, col="blue", ylim=c(min(c(theta1$V9, theta2$V9, theta3$V9), na.rm=T), max(c(theta1$V9, theta2$V9, theta3$V9), na.rm=T)), xlab="Chromosome 11", ylab="Tajima's D (~Pi-Watt)", xlim=xl, ty="l")
points(x=theta2$V3, y=theta2$V9, pch=16, col="black", ty="l")
points(x=theta3$V3, y=theta3$V9, pch=16, col="green", ty="l")
legend("topleft", legend=c("LWK", "TSI", "PEL"), col=c("blue","black","green"), pch=16)

plot(x=theta1$V3, y=theta1$V14, pch=16, col="blue", ylim=c(min(c(theta1$V14, theta2$V14, theta3$V14), na.rm=T),win), xlab="Chromosome 11", ylab="Nr. of sites", xlim=xl, ty="l")
points(x=theta2$V3, y=theta2$V14, pch=16, col="black", ty="l")
points(x=theta3$V3, y=theta3$V14, pch=16, col="green", ty="l")
legend("bottomright", legend=c("LWK", "TSI", "PEL"), col=c("blue","black","green"), pch=16)

dev.off()


cat("Output file:", fout, "\n")




