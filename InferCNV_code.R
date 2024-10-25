
rm(list = ls())
library(tidyverse)
library(Seurat)
library(infercnv)
library(magrittr)
library(RColorBrewer)
options("Seurat.object.assay.version" = "v3")


load('~/scABE_filter_anno_day31421_fullgene_public.RData')

NTC_path = '~/NTC_subset_day31421_PDCD1/'
# blank_path = '/data01/meisongshi/scABE/result/blank/'

if (!dir.exists(NTC_path)) {dir.create(NTC_path)}
# if (!dir.exists(blank_path)) {dir.create(blank_path)}


for (cell in c('B')){#,'CD8',,'B','Gamma_delta_T'
  sce <- subset(sce, subset=(anno_level1 == cell))
  print(table(sce@meta.data$anno_level1))
  
  raw_meta = sce@meta.data
  sam_info = raw_meta %>% 
  distinct(sample_group,condition,target,time)
  print(sam_info)
  
  load('/data01/huangpinzheng/scABE/data/001_gtf_raw.RData')
  
  counts <- GetAssayData(sce, slot = 'counts')
  #counts[1:10,1:100]
  
  gene.chr = gtf_raw %>%subset(.,type == "gene"& gene_name %in% rownames(counts)) %>% #&gene_type=='protein_coding' 
    as.data.frame() %>%
    dplyr::select(gene_name, seqnames, start, end)%>%
    dplyr::distinct(gene_name, .keep_all=T)%>%
    set_rownames(.$gene_name)%>%
    dplyr::select(-gene_name)
  print('------------------gene.chr----------------------')
  head(gene.chr,10)
  
  options(scipen = 100)
  
  
  ###################################### NTC as ref #############################################################
  for (cond in c('ABE','Cas9')) {#
  # cond = 'ABE'
    for (tim in c('day3','day7','day21')) {
      # tim='day3'
      se_info = sam_info %>% 
        filter(.,condition==cond&time==tim)
      
      se_sample = as.character(se_info$sample_group)
      se_sample = se_sample[str_split(se_sample,'-',simplify=T)[,2]!='B2M']
      print(se_sample)
      ref = se_sample[str_split(se_sample,'-',simplify = T)[,2]=='NTC']
      print('----------------------ref-------------------------')
      print(ref)
      
      
      subsce_obj_1 <- subset(sce, sample_group%in%se_sample[1])
      subsce_obj_1 <- subset(subsce_obj_1, downsample=600)
      subsce_obj_2 <- subset(sce, sample_group%in%se_sample[2])
      subsce_obj_2 <- subset(subsce_obj_2, downsample=600)
      subsce_obj <- merge(subsce_obj_1, subsce_obj_2)
      
      counts <- GetAssayData(subsce_obj, slot = 'counts')
      
      raw_count = counts[rownames(gene.chr),]
      
      anno_df = raw_meta %>% 
        dplyr::select(sample_group)
      
      anno_df_se = anno_df %>% 
        filter(.,rownames(.)%in%colnames(counts)) %>% 
        mutate(sample_group=as.character(sample_group))
      
      
      print('------------------anno_df_se--------------------------------')
      print(head(anno_df_se,10))
      print(table(anno_df_se$sample_group))
      
      count_se = counts[,rownames(anno_df_se)]
      print('------------------count_se--------------------------------')
      #print(count_se)
      
      
      
      infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=count_se,
                                           annotations_file=anno_df_se,
                                           delim="\t",
                                           gene_order_file=gene.chr,
                                           ref_group_names=ref)
      
      infercnv_obj = infercnv::run(infercnv_obj,
                                   min_cells_per_gene = 10,
                                   cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                   out_dir=paste0(NTC_path,cell,'_',cond,'_',tim),  # result_blank/
                                   #analysis_mode="samples" ,
                                   cluster_by_groups=T,  # 根据细胞类型对肿瘤细胞执行单独的聚类，如 cell annotations 文件中定义的那样
                                   # hclust_method="ward.D2",
                                   denoise=T,  # 去噪处理
                                   no_prelim_plot = F,
                                   write_expr_matrix = T,  #重要！！
                                   num_threads=60, # 设置线程数, 多线程运行，加快计算速度
                                   HMM=T)
      
      
  #    infercnv::plot_cnv(infercnv_obj, #上两步得到的infercnv对象
  #                       plot_chr_scale = T, #画染色体全长，默认只画出（分析用到的）基因
  #                       output_filename = paste0(NTC_path,cond,'_',tim,'/better_plot'),output_format = "pdf", #保存为pdf文件
  #                       custom_color_pal =  color.palette(c("steelblue","white","firebrick"), c(2, 2))) #改颜色
      
    }
    
  }

}











