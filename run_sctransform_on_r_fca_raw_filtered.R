#################################
#	Run ScTransform on	#
#	raw s_fca filtered	#
#	     counts 	        #
#################################

#### Load libraries ####
library(tidyverse)
library(ape)
library(Seurat)
library(SeuratDisk)
library(celda)
library(glmGamPoi)
library(limma)
library(sctransform)

#### Set working directory ####
setwd("./fca_relaxed_loom/")

#### Loop over all tissue samples ####
tissues <- c('Antenna', 'Malpighian_tubule', 'Body_wall',
             'Oenocyte', 'Fat_body', 'Gut', 'Ovary',
             'Proboscis_and_maxillary_palps', 'Haltere',
             'Heart', 'Trachea', 'Leg', 'Wing', 'Head', 'Body',
             'Testis', 'Male_reproductive_glands')

#### Read in annotation file ####
gff_data <- read.gff(file = "./fly_cell_atlas_proj/data/dmel-genes-r6.31.gff")
gff_data$gene <- str_match(gff_data$attributes, 'Name=(.*?);Alias=')[1,2]
gff_data <- gff_data %>%
  mutate(gene = str_match(attributes, 'Name=(.*?);')[,2])
mit_gff_data <- gff_data[gff_data$seqid == "mitochondrion_genome",]
rm(gff_data)
gc()

#### Run normalisation step ####
for (tissue in tissues) {
  #### Read in Seurat file ####
  tissue_seurat <- LoadH5Seurat(paste0("./r_raw_", tolower(tissue), 
                                       "_filtered.h5Seurat"))
  if ("mix" %in% unique(tissue_seurat@meta.data$sex)) {
    tissue_seurat <- subset(x = tissue_seurat,
                            subset = sex != "mix",
                            invert = FALSE)
  }
  tissue_sce <- as.SingleCellExperiment(tissue_seurat)
  tissue_sce <- decontX(x = tissue_sce,
                        z = tissue_sce$annotation,
                        batch = tissue_sce$batch)
  tissue_seurat[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(tissue_sce)))
  rm(tissue_sce)
  gc()
  #### Add percent.mito to metadata ####
  #### Get mit-linked data from Seurat file ####
  mit_linked_genes <- as.numeric(which(rownames(tissue_seurat@assays$decontXcounts@counts) %in% mit_gff_data[, "gene"]))
  mit_linked_data <-  GetAssayData(object = tissue_seurat,
                                   assay = 'decontXcounts',
                                   slot = 'counts')[mit_linked_genes,]
  mit_linked_data <- as.data.frame(mit_linked_data)
  mit_linked_data <- mit_linked_data %>%
    tibble::rownames_to_column(var = "gene") %>%
    gather("cellID", "counts", 2:(ncol(tissue_seurat@assays$decontXcounts@counts) + 1))
  mit_linked_df <- left_join(mit_linked_data,
                             tissue_seurat@meta.data[, c('cellID', 'sex')],
                             by = "cellID")
  rm(mit_linked_genes)
  rm(mit_linked_data)
  gc()
  #### Run stringent filtering ####
  mit_linked_df <- mit_linked_df %>%
    group_by(cellID, sex) %>%
    summarise(total_mit = sum(counts))
  mit_linked_df <- left_join(mit_linked_df,
                             data.frame(cellID = names(colSums(tissue_seurat@assays$decontXcounts@counts)),
                                        total_counts = colSums(tissue_seurat@assays$decontXcounts@counts),
                                        total_genes = colSums(tissue_seurat@assays$decontXcounts@counts != 0),
                                        total_ones = colSums(tissue_seurat@assays$decontXcounts@counts == 1)))
  mit_linked_df$percent.mito <- (mit_linked_df$total_mit/mit_linked_df$total_counts)*100
  mit_linked_df <- as.data.frame(mit_linked_df)
  rownames(mit_linked_df) <- mit_linked_df$cellID
  #### Add outlier status ####
  median_counts <- median(mit_linked_df$total_counts)
  median_genes <- median(mit_linked_df$total_genes)
  mit_linked_df <- mit_linked_df %>%
    mutate(outlier_status = case_when(abs(total_counts - median_counts) > 3*median_counts ~ 'outlier',
                                      abs(total_genes - median_genes) > 3*median_genes ~ 'outlier',
                                      percent.mito > 5 ~ 'outlier',
                                      total_counts < 500 ~ 'outlier',
                                      total_genes < 200 ~ 'outlier',
                                      .default = 'keeper'))
  mit_linked_df <- mit_linked_df %>%
    filter(outlier_status == 'keeper')
  tissue_seurat <- subset(x = tissue_seurat,
                          subset = cellID %in% mit_linked_df$cellID,
                          invert = FALSE)
  tissue_seurat <- AddMetaData(tissue_seurat,
                               mit_linked_df[, c('percent.mito',
                                                 'total_genes',
                                                 'total_ones')])
  rm(median_counts)
  rm(median_genes)
  rm(mit_linked_df)
  gc()
  #### Remove genes expressed in less than 3 cells ####
  cells_per_gene <- rowSums(tissue_seurat@assays$decontXcounts@counts != 0)
  cells_per_gene <- data.frame(gene = names(cells_per_gene),
                               ncells = cells_per_gene)
  cells_per_gene <- cells_per_gene %>%
    filter(ncells >= 3)
  tissue_seurat <- subset(x = tissue_seurat,
                          features = cells_per_gene$gene)
  rm(cells_per_gene)
  gc()
  
  #### Filter to keep cells in Fly Cell Atlas stringent dataset ####
  stringent_seurat <- LoadH5Seurat(paste0("./fly_cell_atlas_proj/data/h5seurat_filtered_files/s_raw_",
                                          tolower(tissue), "_filtered.h5Seurat"))
  
  tissue_seurat <- subset(x = tissue_seurat,
                          subset = cellID %in% stringent_seurat$cellID,
                          invert = FALSE)
  
  ## And add stringent annotation to metadata ##
  stringent_annotation <- stringent_seurat@meta.data[,c(4,6,7)]
  colnames(stringent_annotation) <- c('cellID', 's_annotation', 's_broad_annotation')
  stringent_annotation <- left_join(x = tissue_seurat@meta.data[,c(4,6,7)],
                                    y = stringent_annotation)
  rownames(stringent_annotation) <- stringent_annotation$cellID
  tissue_seurat <- AddMetaData(tissue_seurat,
                               stringent_annotation[,c(4,5)])
  
  #### Remove spermatids from male accessory glands Seurat object ####
  if (tissue == 'Male_reproductive_glands') {
    tissue_seurat <- subset(x = tissue_seurat,
                            subset = annotation != 'spermatid',
                            invert = FALSE)
  }
  ### Set decontX counts to default assay ###
  DefaultAssay(tissue_seurat) <- 'decontXcounts'
  
  #### Run SCTransform v2 on counts ####
  tissue_seurat <- SCTransform(object = tissue_seurat,
                               assay = 'decontXcounts',
                               vst.flavor = 'v2',
                               vars.to.regress = c('nCount_decontXcounts',
                                                   'percent.mito'),
                               verbose = TRUE)
  SaveH5Seurat(object = tissue_seurat,
               filename = paste0('r_', tolower(tissue),
                                 "_cell_type_filt_decontx_stringent_filt_sctransformed.h5Seurat"))
  rm(tissue_seurat)
  gc()
}
