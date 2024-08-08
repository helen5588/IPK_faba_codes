library(dplyr)
library(tidyr)
library(stringr)
###Put the eggNOG result and iprscan result 
egg<-read.csv("eggnog.anno",header=T,sep="\t")
egg[egg==""] <- NA
egg2<-read.csv("iprscan.go.keep",header=T,sep="\t")
#GID     GENENAME
#Vfaba.Hedin2.R2.1g000654.1      Vfaba.Hedin2.R2.1g000654.1
#Vfaba.Hedin2.R2.1g002437.1      Vfaba.Hedin2.R2.1g002437.1
gene_info <- egg %>%dplyr::select(GID = query, GENENAME = query) %>% na.omit() 
gene_info<-rbind(gene_info,egg2)
goterms <- egg %>%dplyr::select(query, GOs) %>% na.omit() %>% filter(str_detect(GOs,"GO")) 

all_go_list=str_split(goterms$GOs,",") 
gene2go <- data.frame(GID = rep(goterms$query, times = sapply(all_go_list, length)), GO = unlist(all_go_list), EVIDENCE = "IEA") %>% filter(str_detect(GO,"GO")) 
###save the ouput to gene2go.txt file
write.table(gene2go,file="gene2go.txt",sep="\t",row.names=F,quote=F) 
read.table("gene2go.txt",header=T,sep="\t")

koterms <- egg %>%dplyr::select(GID = query, KO=KEGG_ko)%>%na.omit()%>% filter(str_detect(KO,"ko")) 


if(!file.exists('kegg_info.RData')){
  
  library(jsonlite)
  library(purrr)
  library(RCurl)
  
  update_kegg <- function(json = "ko00001.json",file=NULL) {
    pathway2name <- tibble(Pathway = character(), Name = character())
    ko2pathway <- tibble(Ko = character(), Pathway = character())
    
    kegg <- fromJSON(json)
    
    for (a in seq_along(kegg[["children"]][["children"]])) {
      A <- kegg[["children"]][["name"]][[a]]
      
      for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
        B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
        
        for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
          pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
          
          pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
          pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
          pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
          
          kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
          
          kos <- str_match(kos_info, "K[0-9]*")[,1]
          
          ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
        }
      }
    }
    
    save(pathway2name, ko2pathway, file = file)
  }
  
  update_kegg(json = "ko00001.json",file="kegg_info.RData")
  
}
load("kegg_info.RData")
head(ko2pathway) 
head(pathway2name)
write.table(pathway2name,file="pathway2name.txt",sep="\t",row.names=F,quote=F) 


library(stringr)
colnames(ko2pathway)=c("KO",'Pathway') 
koterms$KO=str_replace_all(koterms$KO,"ko:","") 
gene2pathway <- koterms %>% left_join(ko2pathway, by = "KO") %>%dplyr::select(GID, Pathway) %>%na.omit() 

##optional
gene2pathway_name<-left_join(gene2pathway,pathway2name,by="Pathway")
write.table(gene2pathway_name,file="gene2pathway_name.txt",sep="\t",row.names=F,quote=F) 

library(clusterProfiler)
library(AnnotationForge)
makeOrgPackage(gene_info=gene_info, go=gene2go, ko=koterms,  pathway=gene2pathway, version="0.0.1", maintainer='Hailin Zhang <zhanghailin5588@gmail.com>', author='Hailin Zhang <zhanghailin5588@gmail.com>',outputDir=".", tax_id="3906", genus="Vicia", species="faba",goTable="go")
install.packages('org.Vfaba.eg.db',repos = NULL, type="source")=scientific_name

####Annotate gene.list using the faba bean's gene database
library(org.Vfaba.eg.db)
library(clusterProfiler)
library(enrichplot)
data <- read.table("gene.list",header=F) 
genes <- as.character(data$V1) 
ego <- enrichGO(gene          = genes, 
                OrgDb         = org.Vfaba.eg.db, 
                keyType       = 'GID', 
                ont           = "ALL", #  "BP", "MF", "CC", "ALL"
                pAdjustMethod = "BH", #  "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                pvalueCutoff  = 0.05, # pvalue，default pvalueCutoff = 0.05
                qvalueCutoff  = 0.3, # qvalue，default qvalueCutoff = 0.2
                readable      = FALSE ) 

write.table(as.data.frame(ego),"go_enrich.csv",sep="\t",row.names =F,quote=F) 

data<-read.table("go_enrich.csv",sep="\t",header=T,quote="")
geneID_all <- unlist(apply(as.matrix(data$geneID),1,function(x) unlist(strsplit(x,'/'))))

ego2<-new("enrichResult", result=data, gene=geneID_all, pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,ontology="ALL",keytype="GID",universe='Unknown',geneSets=list(),organism="Unknown",readable=FALSE) 
###barplot and dotplot for the GO
barplot(ego2, showCategory=30, title="EnrichmentGO")
dotplot(ego2, showCategory=30,color = "pvalue")

kk<-read.table("pathway2gene.txt",sep="\t",header=T,quote="")
KEGG <- enricher(genes, TERM2GENE = kk, TERM2NAME = pathway2name,gson = NULL, pvalueCutoff = 1, qvalueCutoff = 0.3,  pAdjustMethod = "BH", minGSSize = 1)
write.table(as.data.frame(KEGG),"kegg_enrich.csv",sep="\t",row.names =F,quote=F) 

data<-read.table("kegg_enrich.csv",sep="\t",header=T,quote="")
geneID_all <- unlist(apply(as.matrix(data$geneID),1,function(x) unlist(strsplit(x,'/'))))
ego2<-new("enrichResult", result=data, gene=geneID_all, pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,ontology="ALL",keytype="GID",universe='Unknown',geneSets=list(),organism="Unknown",readable=FALSE) 
###barplot and dotplot for the KEGG
barplot(ego2, showCategory=30, title="",color = "pvalue")
dotplot(ego2, showCategory=30,color = "pvalue")
