### object make ################
```
gene.order.file <-"/home/starjjbang/2022/Single_cell/1_data_analysis/2_figure_data/hg19.inferCNV.gtf"
rds_dir = '/home/starjjbang/2022/Single_cell/1_data_analysis/3_rda/'
cnv_re_dir = '/home/starjjbang/2022/Single_cell/1_data_analysis/5_infer_CNV/input/'
infercnv.obj <- CreateInfercnvObject(
  raw_counts_matrix = paste0(cnv_re_dir,'ct_exp.csv'),
  annotations_file = paste0(cnv_re_dir,'ct_annotation.csv'),
  gene_order_file = gene.order.file,
  ref_group_names = c("T cells","B cells"),delim = "\t")
saveRDS(infercnv.obj, paste0(rds_dir,"8_epi_inferCNV_obj.rds"))

infercnv.obj <- CreateInfercnvObject(
  raw_counts_matrix = paste0(cnv_re_dir,'all_ct_exp.csv'),
  annotations_file = paste0(cnv_re_dir,'all_ct_annotation.csv'),
  gene_order_file = gene.order.file,
  ref_group_names = c("T cells","B cells"),delim = "\t")
saveRDS(infercnv.obj, paste0(rds_dir,"9_all_inferCNV_obj.rds"))
```
### running ################
```
rm(list = ls())
library(infercnv)
all_dir = paste0(rds_dir,"9_all_inferCNV_obj.rds")
epi_dir = paste0(rds_dir,"8_epi_inferCNV_obj.rds")
infer_cnv_run(epi_dir,"epi") #### inferCNV run
infer_cnv_run(all_dir,"all")
```

infer_cnv_run <-function(in_dat_dir,fn){
  cnv_re_dir = paste0('/home/starjjbang/2022/Single_cell/1_data_analysis/5_infer_CNV/',fn)
  dir.create(file.path('/home/starjjbang/2022/Single_cell/1_data_analysis/5_infer_CNV/', fn,'/'))
  infercnv.obj = readRDS(in_dat_dir)
  infercnv_obj_default <- infercnv::run(infercnv.obj,cutoff=0.1, 
                                        out_dir = cnv_re_dir,
                                        window_length=101,
                                        max_centered_threshold=3,
                                        cluster_by_groups=F,
                                        plot_steps=FALSE,
                                        denoise=TRUE,
                                        sd_amplifier=1.3,
                                        analysis_mode = "samples",
                                        HMM=FALSE,
                                        no_prelim_plot=TRUE,
                                        png_res=300)
}


### ture epithelial cell filtering ################
```
include_group_annotation <- FALSE
adjust_normal_thresholds <- TRUE
cancer_x_threshold_sd_multiplier <- 2
cancer_y_threshold_sd_multiplier <- 1.5
normal_x_threshold_sd_multiplier <- 1
normal_y_threshold_sd_multiplier <- 1.25

library(Seurat)
library(infercnv)
library(RColorBrewer)
library(ComplexHeatmap)#BiocManager::install("ComplexHeatmap")
library(circlize)
library(reshape2)
library(ggplot2)
library(scales)
library(fpc)
library(dplyr)
library(naturalsort)
library(cowplot)

func_dir = '/home/starjjbang/2022/code/inferCNV/functions/'
rds_dir = '/home/starjjbang/2022/Single_cell/1_data_analysis/3_rda/'
```
