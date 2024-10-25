
rm(list = ls())

library(clusterProfiler)
library(GSVA)
library(tidyverse)
library(GSEABase)

load("~/scABE_filter_anno_day31421_fullgene_public.RData")


sce <- sce@assays$originalexp@counts
go_list <- getGmt("~/KEGG_P53_SIGNALING_PATHWAY.v2023.1.Hs.gmt")
gsva_mat <- gsva(expr=as.matrix(sce) , 
                 gset.idx.list=go_list, 
                 method = 'ssgsea',
                 kcdf="Gaussian" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                 min.sz=1,
                 abs.ranking=F,
                 verbose=T, 
                 ssgsea.norm=TRUE,
                 parallel.sz = parallel::detectCores())#调用所有核
write.csv(gsva_mat,"~/gsva_go_matrix.csv")



