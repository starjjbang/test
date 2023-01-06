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

### 1. Load InferCNV output and create heatmap and metadata dfs ###
```{r}
library(Seurat)
paste0('/home/starjjbang/2022/Single_cell/1_data_analysis/5_infer_CNV/epi/')
infercnv_output <- data.frame(t(read.table('/home/starjjbang/2022/Single_cell/1_data_analysis/5_infer_CNV/epi/infercnv.observations.txt',sep=' ')))
infercnv_output[1:3,1:3]
harmony.obj = readRDS(paste0(rds_dir,"7_1_cell_ano_garnett_with_hamony.rds"))
metadata <- prepare_infercnv_metadata(harmony.obj, subset_data = F, as.data.frame(t(infercnv_output)), for_infercnv=F)
epithelial_metadata <- metadata$metadata
epithelial_ids <- epithelial_metadata$cell_ids[grep("pithelial", epithelial_metadata$cell_type)]
epithelial_heatmap <- infercnv_output[rownames(infercnv_output) %in% epithelial_ids,]
epithelial_heatmap[1:3,1:3]
```

### 2. Add QC metadata ###
```{r}
# create epithelial_metadata df and only include epithelial cells in epithelial_heatmap:
epithelial_metadata <- epithelial_metadata[rownames(epithelial_heatmap),]
# order epithelial metadata cell type cluster levels:
epithelial_metadata$cell_type <- factor(
  epithelial_metadata$cell_type,
  levels = naturalsort(unique(epithelial_metadata$cell_type))
)
QC <- data.frame(
  row.names = names(Idents(harmony.obj)),
  nUMI = harmony.obj@meta.data$nCount_RNA,
  nGene = harmony.obj@meta.data$nFeature_RNA
)
QC <- QC[rownames(epithelial_metadata),]
epithelial_metadata <- cbind(epithelial_metadata, QC)

rds_dir = '/home/starjjbang/2022/Single_cell/1_data_analysis/3_rda/'
saveRDS(epithelial_heatmap, paste0(rds_dir,"10_epithelial_heatmap_with_cell_type_and_QC.rds"))
saveRDS(epithelial_metadata, paste0(rds_dir, "11_epithelial_metadata_with_cell_type_and_QC.Rdata"))
```

### 3. Add CNA and correlation value metadata ###
```{r}
source('/home/starjjbang/2022/code/paper_code/single_function.R')
rds_dir = '/home/starjjbang/2022/Single_cell/1_data_analysis/3_rda/'

# determine CNA values and add to epithelial_metadata:
#print("Determining CNA values and adding to epithelial metadata df...")
# scale infercnv values to -1:1, square values and take the mean:
scaled_df <- as.data.frame(rescale(as.matrix(epithelial_heatmap), c(-1,1)))
dim(scaled_df)

#the mean of the squares of these values was used to define a genomic instability score for each cell.
CNA_values <- apply(scaled_df, 1, function(y) {
  #y[is.na(y)] <- 0
  #scaled_y <- rescale(y, c(-1, 1))
  return(mean(y^2))
})
CNA_value_df <- data.frame(row.names = names(CNA_values),CNA_value = CNA_values)
dim(CNA_value_df)
epithelial_metadata <- cbind(epithelial_metadata, CNA_value_df)
saveRDS(epithelial_metadata, paste0(rds_dir,"12_epithelial_metadata_with_cell_type_and_QC_CNV_value.rds"))

# determine top 5% cancer cells:
CNA_order <- order(CNA_value_df$CNA_value, decreasing=T)
ordered_CNA_values  <- data.frame(row.names = rownames(CNA_value_df)[CNA_order],CNA_value = CNA_value_df[CNA_order,])
top_cancer <- head(ordered_CNA_values, nrow(ordered_CNA_values)*0.05)

# find average genome-wide CNV predictions across genome:
top_cancer_CNV_average <- apply(epithelial_heatmap[rownames(top_cancer),], 2, mean)
epithelial_heatmap[1:3,1:3]
dim(epithelial_heatmap)
saveRDS(epithelial_heatmap, paste0(rds_dir,"13_epithelial_heatmap_with_cell_type_and_QC_CNV_avg.rds"))

# find correlations of each cell's CNVs with top_GIN_CNV_average:
for(i in 1:nrow(epithelial_heatmap)){
  x = epithelial_heatmap[i,]
  if (length(unique(as.numeric(x))) == 1) {
    cor_result <- data.frame(cor.estimate="no_CNVs_recorded", 
                             cor.p.value="no_CNVs_recorded")
  } else {
    cor <- cor.test(as.numeric(x), top_cancer_CNV_average, method = "kendall")
    cor_result <- data.frame(cor$estimate, cor$p.value)
  }
  rownames(cor_result)=rownames(epithelial_heatmap)[i]
  write.table(cor_result[1,],'/home/starjjbang/2022/Single_cell/1_data_analysis/5_infer_CNV/epi/normal_filter/cor.txt',
              col.names = FALSE,quote=FALSE,sep='\t', append=TRUE)
  print(i)
}
```

### 4a. Call normals and add to epithelial_metadata ###
```{r}
rds_dir = '/home/starjjbang/2022/Single_cell/1_data_analysis/3_rda/'
epithelial_metadata = readRDS(paste0(rds_dir,"12_epithelial_metadata_with_cell_type_and_QC_CNV_value.rds"))
correlation_df <-unique(as.data.frame(read.table('/home/starjjbang/2022/Single_cell/1_data_analysis/5_infer_CNV/epi/normal_filter/cor.txt',sep='\t')))
colnames(correlation_df) = c('cell_ids','cor.estimate','cor.p.value')
epithelial_metadata <- merge(epithelial_metadata,correlation_df,by="cell_ids")
epithelial_metadata = epithelial_metadata[,-5]
saveRDS(epithelial_metadata, paste0(rds_dir,"14_epithelial_metadata_with_cell_type_and_QC_CNV_value_cor.rds"))


quad_df <- data.frame(
  row.names = rownames(epithelial_metadata),
  CNA_value = epithelial_metadata$CNA_value, 
  cor.estimate = epithelial_metadata$cor.estimate
)
library(dplyr)
library(fpc)
library(cluster )
scaled_quad_df <- scale(quad_df) %>% as.data.frame()
# run silhouette cluster analysis to determine clusters and thresholds:
pamk_result <- pamk(scaled_quad_df, krange=1:4)
pamk_result$nc
silhouette_result <- pam(scaled_quad_df, pamk_result$nc)
saveRDS(silhouette_result, paste0(rds_dir,"15_silhouette_result.rds"))

sil_values <- as.data.frame(silhouette_result$silinfo$widths)
sil_result <- data.frame(row.names=names(silhouette_result$clustering),
                         cluster=silhouette_result$clustering,
                         sil_width=sil_values$sil_width)

# add sil_result to epithelial_metadata:
epithelial_metadata = readRDS(paste0(rds_dir,"14_epithelial_metadata_with_cell_type_and_QC_CNV_value_cor.rds"))
epithelial_metadata <- cbind(epithelial_metadata, sil_result)
saveRDS(epithelial_metadata,paste0(rds_dir,"16_epithelial_metadata_with_cell_type_and_QC_CNV_value_cor_sil_result.rds"))


