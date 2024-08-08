library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
mydata<-read.table(args[1],header = T) 

#dose.labs <- c("chr1_partA", "chr1_partB", "chr2","chr3","chr4","chr5","chr6")
#names(dose.labs) <- c("1", "2", "3","4","5","6","7") 
mydata$start <- ifelse(mydata$chr == "chr1_partB", mydata$start+2000000000,mydata$start)
#mydata$end <- ifelse(mydata$chr == "chr1_partB", mydata$end+2000000000,mydata$end)
mydata$chr[mydata$chr == "chr1_partA"]<-"chr1"
mydata$chr[mydata$chr == "chr1_partB"]<-"chr1"
mydata1<-na.omit(mydata)
#pdf(args[2], width = 8, height = 8)
png(args[2], width = 1200, height = 1300)
ggplot(mydata1, aes(x=start/1000000,y=value,group=factor(chr)))+geom_point()+facet_grid(chr ~ .)+#,labeller = labeller(chr = dose.labs))
xlab("Physical distance(Mb)")+ylab("Fst value")+theme(legend.position = "none")+
 theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(color = 'black'),text=element_text(family="Times"), #title=element_text(size=16,hjust=0.2,lineheight=0.2),
       axis.title.x=element_text(size=30,face="bold",color="black",hjust=0.5),
       axis.title.y=element_text(size=30,face="bold",color="black",hjust=0.5),
       axis.text.x=element_text(family="Times",size=30,color="black",angle=45,hjust=1),
       axis.text.y=element_text(family="Times",size=30,color="black"),legend.text=element_text(size=30),strip.text = element_text(
           size = 25)) 
library(qqman)
#pdf(args[3], width = 15, height = 6)
#mydata1$snp<-paste(mydata1$chr,mydata1$ps,sep='_')
mydata1$chr[mydata1$chr == "chr1"]<-"1"
mydata1$chr[mydata1$chr == "chr2"]<-"2"
mydata1$chr[mydata1$chr == "chr3"]<-"3"
mydata1$chr[mydata1$chr == "chr4"]<-"4"
mydata1$chr[mydata1$chr == "chr5"]<-"5"
mydata1$chr[mydata1$chr == "chr6"]<-"6"
mydata1$chr <- as.numeric(mydata1$chr)
mydata1$start <- as.numeric(mydata1$start)
mydata1$value <- as.numeric(mydata1$value)
mydata1$snp<-paste(mydata1$chr,mydata1$start,sep='_')
#manhattan(mydata1, chr="chr",bp="start", snp="Gene", p="value",col = c("blue4", "orange3"),suggestiveline = F, genomewideline = 0.5,logp=F,annotateTop=T)
#dev.off()

png(args[3], width = 1000, height = 600)
manhattan(mydata1, cex=1,chr="chr",bp="start", snp="Gene", p="value",col = c("blue4", "orange3"),suggestiveline = F, genomewideline = 0.5,logp=F,annotateTop=T,annotatePval=0.5,ylab="Fst",las=1, bty='l',chrlabs=NULL)
dev.off()
