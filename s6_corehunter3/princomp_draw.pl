#!/usr/bin/perl -w 
use strict;
use Cwd 'abs_path';
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw/max min sum/;
=head1 Version
  Author: Hailin Zhang,zhangh@ipk-gatersleben.de
  data : 2023.10.24

=head1 Description
  generate R script for Principal Components Analysis to draw several figures displaying the analysis result.
  input: vcf of snp among samples. 

=head1 Usage
  perl princomp_draw.pl -vcf <vcf file> [--core <core in file> --o <output file>]
        --vcf <str>   vcf format
        --core <str>      core numbers information,optional
        --o <str>          Out directory(default ".")
	--exclude <str>    Excluding ID list 
        --include <str>   including ID list
	--help             get help information and exit
	--phe <str>   phenotype format like the following:

  ID,trait1,trait2,trait3,trait4
  TYPE,RF,RF,RF,RF
  GPID_00002,-0.265441456402348,0.175503202460647,-1.75236527888222,2.99999999999992
  GPID_00003,-0.265441456402348,0.450335468307614,-2.74228860055942,3.99999999999992
  GPID_00004,-0.265441456402348,0.175503202460647,-0.983815758403017,4.99999999999992
  
  #####The meaining of the second row#####
  Variable type Code Default data type(first charactor)
  Nominal	N String
  Ordinal	O Integer
  Interval I Integer
  Ratio R	Double
  
  Data type Code(second charactor)
  Boolean B
  Short	T
  Integer I
  Long L
  Float	F
  Double D
  Big Integer R
  Big Decimal M
  Date A
  String S
  None X
=head1 Example
  eg 1: perl princomp_draw.pl --vcf sample.vcf.gz
  eg 2: perl princomp_draw.pl --vcf sample.vcf.gz --core 12 --o workdir
  eg 3: perl princomp_draw.pl --vcf sample.vcf.gz --core 12 --phe PROFABA_GBS.nan.csv  --o workdir

=cut

my ($vcf, $Outdir, $core_info, $Help, $Rscript, $phe,$exclude,$include);
GetOptions(
	"vcf:s"=>\$vcf,
	"o:s"=>\$Outdir,
	"help"=>\$Help,
	"R:s"=>\$Rscript,
	"core:s" => \$core_info,
	"phe:s" => \$phe,
	"exclude:s" =>\$exclude,
	"include:s" =>\$include
	
);
die `pod2text $0` if (@ARGV != 0 || !$vcf || $Help);
$Outdir ||= "." if (!defined $Outdir);
$core_info ||= 0.2 if (!defined $core_info); #Corehunter default parameters for getting 20% samples
$core_info=int($core_info);
$Rscript = "/opt/Bio/R/3.5.1-conda/bin/Rscript" unless ($Rscript and -f $Rscript);
$phe ||= "" if (!defined $phe);
$exclude ||= "" if (!defined $exclude);
$include ||= "" if (!defined $include);
#######################Set up the environment#####################
#export LD_LIBRARY_PATH=/usr/lib64:/usr/lib:/filer-dg/agruppen/dg7/fengj/genome/hbrevisubulatum/hic_1M/lib:/opt/Bio/jdk/1.8.0_66/jre/lib/amd64/server/:/opt/Bio/jdk/1.8.0_66/jre/lib/amd64;

#######################PCA_plink#######################
my $ttm="\$0";
system "/opt/Bio/bcftools/1.15.1/bin/bcftools view -H $vcf | cut -f 1 | uniq | awk '{print $ttm,$ttm}'|sed 's/ /\t/g' > filename.chrom-map.txt;"; #20240102 add
system "/opt/Bio/vcftools/0.1.16/bin/vcftools --gzvcf $vcf --plink --chrom-map filename.chrom-map.txt  --out $vcf.pca ";
system "/opt/Bio/plink/1.90b6.9/bin/plink --noweb --file $vcf.pca --make-bed --out $vcf.pca_bfile --allow-extra-chr ;";
system "/opt/Bio/plink/1.90b6.9/bin/plink --threads 1 --bfile  $vcf.pca_bfile --allow-extra-chr --pca 100 --out $vcf.pca_bfile;" ;

#######################distance file for the input of Corehunter3#######################
system "/filer-dg/agruppen/dg6/zhangh/software/VCF2Dis/bin/VCF2Dis -InPut $vcf -OutPut $vcf.distance.tmp";
system "sed -i '1d' $vcf.distance.tmp;cut -f 1 $vcf.distance.tmp|sed '1i\\ID' >$vcf.distance.tmp2;python3 /filer-dg/agruppen/dg6/zhangh/software/BIN/transfer.py $vcf.distance.tmp2 $vcf.distance.tmp3;cat $vcf.distance.tmp3 $vcf.distance.tmp |sed 's/ //g'|sed 's/\t\$//g'|sed 's/\t/,/g' >$vcf.distance.csv;rm -rf $vcf*.tmp*";
#######################IBS distance file for the input of Corehunter3#######################
system "/opt/Bio/plink/1.90b6.9/bin/plink --bfile $vcf.pca_bfile --allow-extra-chr --distance-matrix --out $vcf.pca_bfile;paste -d '\t' $vcf.pca_bfile.mdist.id $vcf.pca_bfile.mdist |cut -f 2- >$vcf.pca_bfile.mdist.tmp;head -1 $vcf.distance.csv >$vcf.pca_bfile.mdist.tmp2;cat $vcf.pca_bfile.mdist.tmp2 $vcf.pca_bfile.mdist.tmp |sed 's/\t/,/g'|sed 's/ /,/g'|sed 's/,\$//g' >$vcf.pca_bfile.mdist.csv;";
#######################generate R script#####################
my $keyname = basename($vcf);
open RSC,"> $Outdir/$keyname.R" or die $!;
print RSC <<EOF;
.libPaths(c("/filer-dg/agruppen/DG/zhangh/path/to/Rlib/folder","/opt/Bio/R/3.5.1-conda/lib/R/library","/opt/Bio/R_LIBS/3.5.1","/opt/Bio/jdk/1.8.0_66/jre/lib/amd64/server"))
options(java.parameters = "-Xmx10G")
library(rJava)
library(corehunter)

#######################corehunter from genotype#######################
if (nchar("$exclude") >0 ){
excl <-read.table(file="$exclude",header = F,stringsAsFactors = F,sep="\t")
excl <- c(excl[,1])
}
if (nchar("$include") >0 ){
inc <-read.table(file="$include",header = F,stringsAsFactors = F,sep="\t")
inc <- c(inc[,1])
}
my.geno <- read.csv(file="$vcf.pca.ped",header=F,stringsAsFactors = F,sep="\t")
my.geno.dd <- my.geno[,-(2:6)]

tt <-  paste("Marker", 0:(ncol(my.geno.dd)-1),sep="")
tt[1] <-"ID"
names(my.geno.dd) <- tt
write.table(my.geno.dd, file = "$vcf.pca.ped.csv", sep = ",", row.names = F, col.names = T, quote = F)
my.geno.data <- genotypes( file="$vcf.pca.ped.csv",format = "default")
if (nchar("$exclude") ==0 && nchar("$include")==0){
geno.core.number<-sampleCore(my.geno.data,size=$core_info)
} else if (nchar("$exclude") > 0 && nchar("$include")==0){
geno.core.number<-sampleCore(my.geno.data,size=$core_info,never.selected=c(excl))
} else if (nchar("$exclude") == 0 && nchar("$include")>0){
geno.core.number<-sampleCore(my.geno.data,size=$core_info,always.selected=c(inc))
} else {
geno.core.number<-sampleCore(my.geno.data,size=$core_info,always.selected=c(inc),never.selected=c(excl))
}

eigvec <- read.table("$vcf.pca_bfile.eigenvec", header = F, stringsAsFactors = F)
eigvec <-data.frame(eigvec[,1],eigvec)
tt <-  paste("PC", -2:(ncol(eigvec)-3),sep="")
tt[1] <- "SampleName"
tt[2] <- "Group"
tt[3] <- "Cluster"
names(eigvec) <- tt

eigvec\$Group[eigvec\$Group %in% geno.core.number\$sel] <-"Core"
eigvec\$Cluster[eigvec\$Cluster %in% geno.core.number\$sel] <-"Core"
`%!in%` <- Negate(`%in%`)
eigvec\$Group[eigvec\$Group %!in% geno.core.number\$sel] <-"Others"
geno.core.number\$sel<-c("Core",geno.core.number\$sel)
eigvec\$Cluster[eigvec\$Cluster %!in% geno.core.number\$sel] <-"Others"
write.table(eigvec[1:ncol(eigvec)], file = "$vcf.pca_bfile.geno.eigenvector.xls", sep = "\t", row.names = F, col.names = T, quote = F)

eigval <- read.table("$vcf.pca_bfile.eigenval", header = F)
pcs <- paste("PC", 1:nrow(eigval),sep="")
eigval[nrow(eigval),1] <- 0
percentage <- eigval\$V1/sum(eigval\$V1)*100
eigval_df <- as.data.frame(cbind(pcs, eigval[,1], percentage), stringsAsFactors = F)
names(eigval_df) <- c("PCs", "variance", "proportion")
eigval_df\$variance <- as.numeric(eigval_df\$variance)
eigval_df\$proportion <- as.numeric(eigval_df\$proportion)
names(eigval_df) <- c("PCs", "#eval", "eval%")
write.table(eigval_df[,2:3], file = "$vcf.pca_bfile.geno.eigenvaltor.xls", sep = "\t", quote = F, row.names = F, col.names = T)


#######################corehunter from distance#######################


my.data <- distances (file = "$vcf.distance.csv")
if (nchar("$exclude") ==0 && nchar("$include")==0){
 core.number<-sampleCore(my.data,size=$core_info)
 } else if (nchar("$exclude") > 0 && nchar("$include")==0){
 core.number<-sampleCore(my.data,size=$core_info,never.selected=c(excl))
 } else if (nchar("$exclude") == 0 && nchar("$include")>0){
 core.number<-sampleCore(my.data,size=$core_info,always.selected=c(inc))
 } else {
 core.number<-sampleCore(my.data,size=$core_info,always.selected=c(inc),never.selected=c(excl))
 }

eigvec <- read.table("$vcf.pca_bfile.eigenvec", header = F, stringsAsFactors = F)
eigvec <-data.frame(eigvec[,1],eigvec)
tt <-  paste("PC", -2:(ncol(eigvec)-3),sep="")
tt[1] <- "SampleName" 
tt[2] <- "Group"
tt[3] <- "Cluster"
names(eigvec) <- tt

eigvec\$Group[eigvec\$Group %in% core.number\$sel] <-"Core"
eigvec\$Cluster[eigvec\$Cluster %in% core.number\$sel] <-"Core"
`%!in%` <- Negate(`%in%`)
eigvec\$Group[eigvec\$Group %!in% core.number\$sel] <-"Others"
core.number\$sel<-c("Core",core.number\$sel)
eigvec\$Cluster[eigvec\$Cluster %!in% core.number\$sel] <-"Others"
write.table(eigvec[1:ncol(eigvec)], file = "$vcf.pca_bfile.eigenvector.xls", sep = "\t", row.names = F, col.names = T, quote = F)

eigval <- read.table("$vcf.pca_bfile.eigenval", header = F)
pcs <- paste("PC", 1:nrow(eigval),sep="")
eigval[nrow(eigval),1] <- 0
percentage <- eigval\$V1/sum(eigval\$V1)*100
eigval_df <- as.data.frame(cbind(pcs, eigval[,1], percentage), stringsAsFactors = F)
names(eigval_df) <- c("PCs", "variance", "proportion")
eigval_df\$variance <- as.numeric(eigval_df\$variance)
eigval_df\$proportion <- as.numeric(eigval_df\$proportion)
names(eigval_df) <- c("PCs", "#eval", "eval%")
write.table(eigval_df[,2:3], file = "$vcf.pca_bfile.eigenvaltor.xls", sep = "\t", quote = F, row.names = F, col.names = T)



#20231103 add the IBS_distance
#
my.data <- distances (file = "$vcf.pca_bfile.mdist.csv")
if (nchar("$exclude") ==0 && nchar("$include")==0){
  core.number<-sampleCore(my.data,size=$core_info)
   } else if (nchar("$exclude") > 0 && nchar("$include")==0){
     core.number<-sampleCore(my.data,size=$core_info,never.selected=c(excl))
      } else if (nchar("$exclude") == 0 && nchar("$include")>0){
        core.number<-sampleCore(my.data,size=$core_info,always.selected=c(inc))
	 } else {
	  core.number<-sampleCore(my.data,size=$core_info,always.selected=c(inc),never.selected=c(excl))
	   }

	   eigvec <- read.table("$vcf.pca_bfile.eigenvec", header = F, stringsAsFactors = F)
	   eigvec <-data.frame(eigvec[,1],eigvec)
	   tt <-  paste("PC", -2:(ncol(eigvec)-3),sep="")
	   tt[1] <- "SampleName"
	   tt[2] <- "Group"
	   tt[3] <- "Cluster"
	   names(eigvec) <- tt

	   eigvec\$Group[eigvec\$Group %in% core.number\$sel] <-"Core"
	   eigvec\$Cluster[eigvec\$Cluster %in% core.number\$sel] <-"Core"
	   `%!in%` <- Negate(`%in%`)
	   eigvec\$Group[eigvec\$Group %!in% core.number\$sel] <-"Others"
	   core.number\$sel<-c("Core",core.number\$sel)
	   eigvec\$Cluster[eigvec\$Cluster %!in% core.number\$sel] <-"Others"
	   write.table(eigvec[1:ncol(eigvec)], file = "$vcf.pca_bfile.ibs.eigenvector.xls", sep = "\t", row.names = F, col.names = T, quote = F)

	   eigval <- read.table("$vcf.pca_bfile.eigenval", header = F)
	   pcs <- paste("PC", 1:nrow(eigval),sep="")
	   eigval[nrow(eigval),1] <- 0
	   percentage <- eigval\$V1/sum(eigval\$V1)*100
	   eigval_df <- as.data.frame(cbind(pcs, eigval[,1], percentage), stringsAsFactors = F)
	   names(eigval_df) <- c("PCs", "variance", "proportion")
	   eigval_df\$variance <- as.numeric(eigval_df\$variance)
	   eigval_df\$proportion <- as.numeric(eigval_df\$proportion)
	   names(eigval_df) <- c("PCs", "#eval", "eval%")
	   write.table(eigval_df[,2:3], file = "$vcf.pca_bfile.ibs.eigenvaltor.xls", sep = "\t", quote = F, row.names = F, col.names = T)
EOF
close RSC;
#######################corehunter from genotype and phenotypes#######################
if ($phe ne ""){
  open RSC2,"> $Outdir/$keyname.phe.R" or die $!;
  print RSC2 <<EOF;
.libPaths(c("/filer-dg/agruppen/DG/zhangh/path/to/Rlib/folder","/opt/Bio/R_LIBS/3.5.1","/opt/Bio/jdk/1.8.0_66/jre/lib/amd64/server"))
options(java.parameters = "-Xmx10G")
library(rJava)
library(corehunter)
#######################corehunter from genotype and phenotype#######################
if (nchar("$exclude") >0 ){
 excl <-read.table(file="$exclude",header = F,stringsAsFactors = F,sep="\t")
 excl <- c(excl[,1])
 }
 if (nchar("$include") >0 ){
 inc <-read.table(file="$include",header = F,stringsAsFactors = F,sep="\t")
 inc <- c(inc[,1])
 }
my.geno <- read.csv(file="$vcf.pca.ped",header=F,stringsAsFactors = F,sep="\t")
my.phe <- read.csv(file="$phe",header=F,stringsAsFactors = F,sep=",")
my.phe <- na.omit(my.phe)
start_row <- 3
end_row <- nrow(my.phe)
sort_index <- order(my.phe\$V1[start_row:end_row])
sorted_phe <- my.phe
sorted_phe[start_row:end_row, ] <- my.phe[start_row:end_row, ][sort_index, ]
write.table(sorted_phe, file = "$phe.csv", sep = ",", row.names = F, col.names = F, quote = F)
my.geno.dd <- my.geno[my.geno\$V1 %in% my.phe\$V1[start_row:end_row],-(2:6)]
#tt <-  paste("Marker", 0:(ncol(my.geno.dd)-1),sep="")
#tt[1] <-"ID"
#names(my.geno.dd) <- tt
sort_index <- order(my.geno.dd\$V1)
sorted_geno <- my.geno.dd
sorted_geno<-my.geno.dd[sort_index, ]
tt <-  paste("Marker", 0:(ncol(my.geno.dd)-1),sep="")
tt[1] <-"ID"
names(sorted_geno) <- tt
write.table(sorted_geno, file = "$vcf.pca.ped.phe.csv", sep = ",", row.names = F, col.names = T, quote = F)
my.geno.data <- genotypes( file="$vcf.pca.ped.phe.csv",format = "default")
my.pheo.data <-phenotypes(file = "$phe.csv")
my.data<-coreHunterData(my.geno.data,my.pheo.data)
if (nchar("$exclude") ==0 && nchar("$include")==0){
 phe.core.number<-sampleCore(my.data,size=$core_info)
 } else if (nchar("$exclude") > 0 && nchar("$include")==0){
 excla<-excl[excl %in% sorted_geno\$ID]
 phe.core.number<-sampleCore(my.data,size=$core_info,never.selected=c(excla))
 } else if (nchar("$exclude") == 0 && nchar("$include")>0){
 inca<-inc[inc %in% sorted_geno\$ID]
 `%!in%` <- Negate(`%in%`)
 ddt<-inc[inc %!in% sorted_geno\$ID]
 phe.core.number<-sampleCore(my.data,size=$core_info-length(ddt),always.selected=c(inca))
 } else {
 excla<-excl[excl %in% sorted_geno\$ID]
 inca<-inc[inc %in% sorted_geno\$ID]
 `%!in%` <- Negate(`%in%`)
 ddt<-inc[inc %!in% sorted_geno\$ID]
 phe.core.number<-sampleCore(my.data,size=$core_info-length(ddt),always.selected=c(inca),never.selected=c(excla))
 
 }

eigvec <- read.table("$vcf.pca_bfile.eigenvec", header = F, stringsAsFactors = F)
eigvec <-data.frame(eigvec[,1],eigvec)
tt <-  paste("PC", -2:(ncol(eigvec)-3),sep="")
tt[1] <- "SampleName"
tt[2] <- "Group"
tt[3] <- "Cluster"
names(eigvec) <- tt
if (!exists('ddt')){
eigvec\$Group[eigvec\$Group %in% phe.core.number\$sel] <-"Core"
eigvec\$Cluster[eigvec\$Cluster %in% phe.core.number\$sel] <-"Core"
`%!in%` <- Negate(`%in%`)
eigvec\$Group[eigvec\$Group %!in% phe.core.number\$sel] <-"Others"
phe.core.number\$sel<-c("Core",phe.core.number\$sel)
eigvec\$Cluster[eigvec\$Cluster %!in% phe.core.number\$sel] <-"Others"
} else{
phe.core.number\$sel<-c(phe.core.number\$sel,ddt)
eigvec\$Group[eigvec\$Group %in% phe.core.number\$sel] <-"Core"
eigvec\$Cluster[eigvec\$Cluster %in% phe.core.number\$sel] <-"Core"
`%!in%` <- Negate(`%in%`)
eigvec\$Group[eigvec\$Group %!in% phe.core.number\$sel] <-"Others"
phe.core.number\$sel<-c("Core",phe.core.number\$sel)
eigvec\$Cluster[eigvec\$Cluster %!in% phe.core.number\$sel] <-"Others"
}

write.table(eigvec[1:ncol(eigvec)], file = "$vcf.pca_bfile.phe.eigenvector.xls", sep = "\t", row.names = F, col.names = T, quote = F)

eigval <- read.table("$vcf.pca_bfile.eigenval", header = F)
pcs <- paste("PC", 1:nrow(eigval),sep="")
eigval[nrow(eigval),1] <- 0
percentage <- eigval\$V1/sum(eigval\$V1)*100
eigval_df <- as.data.frame(cbind(pcs, eigval[,1], percentage), stringsAsFactors = F)
names(eigval_df) <- c("PCs", "variance", "proportion")
eigval_df\$variance <- as.numeric(eigval_df\$variance)
eigval_df\$proportion <- as.numeric(eigval_df\$proportion)
names(eigval_df) <- c("PCs", "#eval", "eval%")
write.table(eigval_df[,2:3], file = "$vcf.pca_bfile.phe.eigenvaltor.xls", sep = "\t", quote = F, row.names = F, col.names = T)

EOF
close RSC2;
}
print STDERR "Rscript generated successfully!\nThen wait a moment please!\n";
die "Error in Rscript execution!\n Please check your R's PATH!\n" if system("$Rscript $Outdir/$keyname.R >$Outdir/$keyname.Rout");
if ($phe ne ""){
die "Error in Rscript execution!\n Please check your R's PATH!\n" if system("$Rscript $Outdir/$keyname.phe.R >$Outdir/$keyname.phe.Rout");
}
system("rm -rf  $Outdir/$keyname.R $Outdir/$keyname.Rout");
system("/filer-dg/agruppen/dg6/zhangh/software/BIN/Plot2Deig3 -InFile $vcf.pca_bfile.geno.eigenvector.xls -OutPut $Outdir/PCA.geno -Columns 4:5 -ShowEval -ColorBrewer Set1");
system("/filer-dg/agruppen/dg6/zhangh/software/BIN/Plot2Deig3 -InFile $vcf.pca_bfile.geno.eigenvector.xls -OutPut $Outdir/PCA.geno -Columns 5:6 -ShowEval -ColorBrewer Set1");
system("/filer-dg/agruppen/dg6/zhangh/software/BIN/Plot2Deig3 -InFile $vcf.pca_bfile.eigenvector.xls -OutPut $Outdir/PCA.distance -Columns 4:5 -ShowEval -ColorBrewer Set1");
system("/filer-dg/agruppen/dg6/zhangh/software/BIN/Plot2Deig3 -InFile $vcf.pca_bfile.eigenvector.xls -OutPut $Outdir/PCA.distance -Columns 5:6 -ShowEval -ColorBrewer Set1");
system("/filer-dg/agruppen/dg6/zhangh/software/BIN/Plot2Deig3 -InFile $vcf.pca_bfile.ibs.eigenvector.xls -OutPut $Outdir/PCA.ibs -Columns 4:5 -ShowEval -ColorBrewer Set1");
system("/filer-dg/agruppen/dg6/zhangh/software/BIN/Plot2Deig3 -InFile $vcf.pca_bfile.ibs.eigenvector.xls -OutPut $Outdir/PCA.ibs -Columns 5:6 - ShowEval -ColorBrewer Set1");
if ($phe  ne ""){
system("/filer-dg/agruppen/dg6/zhangh/software/BIN/Plot2Deig3 -InFile $vcf.pca_bfile.phe.eigenvector.xls -OutPut $Outdir/PCA.phe -Columns 4:5 -ShowEval -ColorBrewer Set1");
system("/filer-dg/agruppen/dg6/zhangh/software/BIN/Plot2Deig3 -InFile $vcf.pca_bfile.phe.eigenvector.xls -OutPut $Outdir/PCA.phe -Columns 5:6 -ShowEval -ColorBrewer Set1");
system("grep Core $vcf.pca_bfile.phe.eigenvector.xls |cut -f 1 >Core.lst");
}

system("rm -rf *.svg *.png");
#system("convert -density 120 $Outdir/PCA.Com3-Com4.pdf $Outdir/PCA.Com3-Com4.png");
print STDERR "Finished!\nPlease enjoy the beautiful graphics using a viewer!\nYou can modify the R script and re-run it to get more suitable figure(s)~~\nDefault is keeping core list from phenotype and genotype.You can also choose p_distance, ibs or just genotype core samples.\n";
