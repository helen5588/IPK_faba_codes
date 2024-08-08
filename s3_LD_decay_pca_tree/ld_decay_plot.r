read.table("Spring.ld.decay.bin")->data
##Dist   Mean_r^2        Mean_D' Sum_r^2 Sum_D'  NumberPairs
pdf("Spring.ld.decay.pdf")
plot(data[,1]/1000,data[,2],type="l",col="blue",main="LD decay",xlab="Distance(Kb)",xlim=c(0,300),ylim=c(0,0.4),ylab=expression(r^{2}),bty="n",lwd=2)
dev.off()
png("Spring.ld.decay.png")
plot(data[,1]/1000,data[,2],type="l",col="blue",main="LD decay",xlab="Distance(Kb)",xlim=c(0,300),ylim=c(0,0.4),ylab=expression(r^{2}),bty="n",lwd=2)
dev.off()
