###################################
#     Running DESeq2 analysis     #
#      on head and body data      #
#       at the tissue- and        #
#          cluster-level          #
###################################

#### Load libraries ####
library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)
library(scuttle)
library(DESeq2)
library(tidyverse)

#### Set working directory ####
setwd("./from_s_loom_filtered_sctransformed_h5seurat/")

#### Loop over all tissue samples ####
tissues <- c('Head', 'Body')

#### Run DE analysis ####
### At the tissue level ###
for (tissue in tissues) {
  ## Read in Seurat file ##
  tissue_seurat <- LoadH5Seurat(paste0("r_", tolower(tissue),
                                       "_cell_type_filt_decontx_stringent_filt_sctransformed.h5Seurat"))
  
  ## Remove further sex-specific cells from carcass ##
  if (tissue == 'Body') {
    tissue_seurat <- subset(x = tissue_seurat,
                            subset = s_annotation %in% 
                                         c('follicle cell',
                                           'polar follicle cell'),
                            invert = TRUE)
  }
  gc()
  
  ## Convert Seurat object to Single Cell Experiment object ##
  tissue_sce <- as.SingleCellExperiment(tissue_seurat)
  sum_by <- c("sex", "sample_id")
  
  ## Aggregate counts ##
  summed_sce <- scuttle::aggregateAcrossCells(tissue_sce, id = SummarizedExperiment::colData(tissue_sce)[, sum_by])
  colnames(summed_sce) <- apply(SummarizedExperiment::colData(summed_sce)[, sum_by], 1,
                                function(x) paste(x, collapse = "_"))
  data_for_deseq2 <- SummarizedExperiment::assay(summed_sce, "counts")
  metadata_for_deseq2 <- data.frame(id = apply(SummarizedExperiment::colData(summed_sce)[, sum_by], 1,
                                               function(x) paste(x, collapse = "_")),
                                    sex = SummarizedExperiment::colData(summed_sce)[, 'sex'],
                                    sample_id = SummarizedExperiment::colData(summed_sce)[, 'sample_id'])
  metadata_for_deseq2$sex <- as.factor(metadata_for_deseq2$sex)
  metadata_for_deseq2$sample_id <- as.factor(metadata_for_deseq2$sample_id)
  
  ## Run DESeq2 ##
  deseq_obj <- DESeq2::DESeqDataSetFromMatrix(countData = data_for_deseq2,
                                     colData = metadata_for_deseq2,
                                     design = ~ sex)
  deseq_analysis <- DESeq2::DESeq(deseq_obj)
  comparisons <- c("sex",
                   levels(metadata_for_deseq2$sex)[2],
                   levels(metadata_for_deseq2$sex)[1])
  deseq_results <- DESeq2::results(deseq_analysis,
                           contrast = comparisons,
                           alpha = 0.05)
  deseq_results_df <- deseq_results %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  
  ## Write results to file ##
  write.table(deseq_results_df,
            paste0("deseq_results_from_raw_relaxed_", tolower(tissue), "_per_tissue.csv"),
            sep = ",",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
  rm(tissue_seurat)
  rm(tissue_sce)
  rm(sum_by)
  rm(summed_sce)
  rm(data_for_deseq2)
  rm(metadata_for_deseq2)
  rm(deseq_obj)
  rm(deseq_analysis)
  rm(comparisons)
  rm(deseq_results)
  rm(deseq_results_df)
  gc()
}


### At the cluster level ###
## Function to run DESeq2 from counts aggregated for each cell type cluster ##
get_deseq2_cluster_results <- function(cluster, metadata_df, data_df) {
  metadata_for_deseq2 <- metadata_df[which(metadata_df$s_broad_annotation ==
                                             cluster),]
  data_for_deseq2 <- data_df[, unique(metadata_for_deseq2$id)]
  deseq_obj<- DESeq2::DESeqDataSetFromMatrix(countData = data_for_deseq2,
                                     colData = metadata_for_deseq2,
                                     design = ~ sex)
  deseq_analysis <- DESeq2::DESeq(deseq_obj)
  comparisons <- c("sex",
                   levels(metadata_for_deseq2$sex)[2],
                   levels(metadata_for_deseq2$sex)[1])
  deseq_results <- DESeq2::results(deseq_analysis,
                           contrast = comparisons,
                           alpha = 0.05)
  deseq_results_df <- deseq_results %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  
  deseq_results_df$cluster <- cluster
  return(deseq_results_df)
}

## Use broadly annotated cell types ##
# Read in annotation file
annotation <- read.csv("./fly_cell_atlas_proj/data/cell_types_and_broad_annotation2.csv",
                       header = TRUE)

colnames(annotation) <- c("s_annotation", "broad_annotation")


## Run analysis on all clusters within a tissue ##
for (tissue in tissues) {
  ## Read in Seurat file ##
  tissue_seurat <- LoadH5Seurat(paste0("r_", tolower(tissue),
                                       "_cell_type_filt_decontx_stringent_filt_sctransformed.h5Seurat"))
  gc()
  
  ## Remove further sex-specific cells from carcass ##
  if (tissue == 'Body') {
    tissue_seurat <- subset(x = tissue_seurat,
                            subset = s_annotation %in% 
                              c('follicle cell',
                                'polar follicle cell'),
                            invert = TRUE)
  }
  gc()
  
  ## Merge with stringent annotation ##
  tissue_annotation_df <- left_join(tissue_seurat@meta.data[, c('cellID',
                                                                's_annotation')],
                                    annotation)
  
  rownames(tissue_annotation_df) <- tissue_annotation_df$cellID
  
  tissue_seurat <- AddMetaData(object = tissue_seurat,
                               metadata = tissue_annotation_df$broad_annotation,
                               col.name = 's_broad_annotation')
  
  ## Convert Seurat object to Single Cell Experiment object ##
  tissue_sce <- as.SingleCellExperiment(tissue_seurat)
  rm(tissue_seurat)
  rm(tissue_annotation_df)
  gc()
  sum_by <- c("sex", "sample_id", "s_broad_annotation")
  
  ## Aggregate counts for each cell type cluster ##
  summed_sce <- scuttle::aggregateAcrossCells(tissue_sce, id = SummarizedExperiment::colData(tissue_sce)[, sum_by])
  colnames(summed_sce) <- apply(SummarizedExperiment::colData(summed_sce)[, sum_by], 1,
                                function(x) paste(x, collapse = "_"))
  data_for_deseq2 <- SummarizedExperiment::assay(summed_sce, "counts")
  metadata_for_deseq2 <- data.frame(id = apply(SummarizedExperiment::colData(summed_sce)[, sum_by], 1,
                                               function(x) paste(x, collapse = "_")),
                                    sex = SummarizedExperiment::colData(summed_sce)[, 'sex'],
                                    sample_id = SummarizedExperiment::colData(summed_sce)[, 'sample_id'],
                                    s_broad_annotation = SummarizedExperiment::colData(summed_sce)[, 's_broad_annotation'])
  rm(tissue_sce)
  rm(sum_by)
  rm(summed_sce)
  gc()
  metadata_for_deseq2$sex <- as.factor(metadata_for_deseq2$sex)
  metadata_for_deseq2$sample_id <- as.factor(metadata_for_deseq2$sample_id)
  metadata_for_deseq2$s_broad_annotation <- as.factor(metadata_for_deseq2$s_broad_annotation)
  clusters <- levels(metadata_for_deseq2$s_broad_annotation)
  
  ## Run DESeq2 and write output to file ##
  for (cluster in clusters) {
    deseq_results_df <- get_deseq2_cluster_results(cluster = cluster,
                                                   metadata_df = metadata_for_deseq2,
                                                   data_df = data_for_deseq2)
    if (match(cluster, clusters) == 1) {
      write.table(deseq_results_df,
                paste0("deseq_results_from_raw_relaxed_", tolower(tissue), "_per_cluster.csv"),
                sep = ",",
                quote = FALSE,
                col.names = TRUE,
                row.names = FALSE)
    } else {
      write.table(deseq_results_df,
                paste0("deseq_results_from_raw_relaxed_", tolower(tissue), "_per_cluster.csv"),
                sep = ",",
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE,
                append = TRUE)
    }
    rm(deseq_results_df)
    gc()
  }
  
  rm(metadata_for_deseq2)
  rm(data_for_deseq2)
  rm(cluster)
  rm(clusters)
  gc()  
}

