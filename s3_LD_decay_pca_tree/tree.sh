#https://github.com/BGI-shenzhen/VCF2Dis.git

VCF2Dis -InPut all.snp.vcf.gz -OutPut All.mat -KeepMF
_fneighbor -datafile All.mat -outfile All.mat.fneighbor -matrixtype s -treetype n -outtreefile ./All.mat.tree
