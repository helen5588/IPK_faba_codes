if(!require(gplots)) install.packages("gplots")
if(!require(ggplot2)) install.packages("ggplot2")
library(gplots)
library(ggplot2)
read.table("all.id.eigenvector.xls",h=T)->PCA

grc = factor(PCA[["Group"]])
shape_level <- nlevels(grc)
if (shape_level < 15){ shapes = (0:shape_level) %% 15}else{ shapes = c(0:14,c((15:shape_level) %% 110 + 18)) }


pdf("new.PC3_PC4.r.pdf",width=7,height=6)
if (shape_level < 6 ){
df=PCA[order(PCA$Cluster),]
ggplot(data = df, mapping = aes(x = PC3, y = PC4,colour=Cluster)) + geom_point(size = 3, alpha = 0.6 ) + #scale_color_brewer(palette="Set1")+#,direction=-1)+#+
#scale_x_continuous(limits=c(-0.1,0.1))+
#scale_y_continuous(limits=c(-0.2,0.2))+
scale_color_manual(values = c("#00008B","#bd0817"))+
#scale_color_manual(values = c("#c9c9c7","#bd0817"))+    #rev(RColorBrewer::brewer.pal(9,Set1)))+ 
theme(plot.title = element_text(hjust = 0.5),legend.key=element_blank(),legend.text=element_text(size=16))+ ##,legend.position = c(0.9,0),legend.justification=c(0.9,0))+
labs(x="PC3 (2.00%)",y="PC4 (1.62%)")+theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(color = 'black'),text=element_text(family="Times"), #title=element_text(size=16,hjust=0.2,lineheight=0.2),
 axis.title.x=element_text(size=16,face="bold",color="black",hjust=0.5),
 axis.title.y=element_text(size=16,face="bold",color="black",hjust=0.5),
 axis.text.x=element_text(family="Times",size=14,color="black",hjust=0.5), #angle=45,hjust=1),
 axis.text.y=element_text(family="Times",size=14,color="black"),legend.title=element_blank())+
 #scale_color_brewer(breaks=c("ICARDA", "IPK_Genebank", "NPZ","PROFABA","Winter_panel"),   labels=c("ICARDA genebank", "IPK genebank", "Elite lines (NPZ)","ProFaba","Winter panel"),palette="Set1")+
guides(fill=guide_legend(title=NULL))
 #+guides(color=guide_legend(reverse=TRUE))
} else {
ggplot(data = PCA, mapping = aes(x = PC3, y = PC4,colour=Cluster)) + geom_point(size = 3 ,alpha = 0.6) +ggtitle("PCA") +scale_fill_brewer(palette="Set1",direction=-1) +theme(plot.title = element_text(hjust = 0.5))+  scale_shape_manual(values=seq(0,15))+labs(x="PC3 (2.00%)",y="PC4 (1.62%)") 	}

dev.off();

