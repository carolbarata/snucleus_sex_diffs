###########################################
#     Analysing Majane et al's (2022)     #
#       dataset for X:A differences       #
#                 in Dmel                 #
##########################################

#### Load libraries ####
library(Matrix)
library(tidyverse)
library(tidyr)
library(ape)
library(SingleCellExperiment)
library(data.table)
library(magrittr)
library(celda)
library(glmGamPoi)
library(Seurat)
library(SeuratDisk)
library(limma)
library(sctransform)
library(ggpubr)
library(ggridges)

#### Set working directory ####
setwd("./majane_2022_data")

#### Read in Dmel raw data ####
## Matrix with counts ##
raw_dmel <- readMM(file = "../majane_2022_data/matrix_raw.mtx")
clean_dmel <- readMM(file = "matrix.mtx")
clean_yak <- readMM(file = "matrix_yak.mtx")
clean_sim <- readMM(file = "matrix_sim.mtx")

## Cell IDs ##
raw_cellIDs <- read.table("../majane_2022_data/barcodes_raw.txt")
clean_cellIDs <- read.table("barcodes.txt")
clean_cellIDs_yak <- read.table("barcodes_yak.txt")
clean_cellIDs_sim <- read.table("barcodes_sim.txt")

## Gene IDs ##
raw_geneIDs <- read.table("../majane_2022_data/geneIDs_raw.txt")
clean_geneIDs <- read.table("geneIDs.txt")
clean_geneIDs_yak <- read.table("geneIDs_yak.txt")
clean_geneIDs_sim <- read.table("geneIDs_sim.txt")

## Attach cell and gene IDs to counts matrix ##
colnames(raw_dmel) <- raw_cellIDs$V1
rownames(raw_dmel) <- raw_geneIDs$V1
colnames(clean_dmel) <- clean_cellIDs$V1
rownames(clean_dmel) <- clean_geneIDs$V1
colnames(clean_yak) <- clean_cellIDs_yak$V1
rownames(clean_yak) <- clean_geneIDs_yak$V1
colnames(clean_sim) <- clean_cellIDs_sim$V1
rownames(clean_sim) <- clean_geneIDs_sim$V1

rm(raw_cellIDs)
rm(raw_geneIDs)
rm(clean_cellIDs)
rm(clean_geneIDs)
rm(clean_cellIDs_yak)
rm(clean_geneIDs_yak)
rm(clean_cellIDs_sim)
rm(clean_geneIDs_sim)
gc()

#### Read in annotation file ####
gff_data <- read.gff(file = "./fly_cell_atlas_proj/data/dmel-genes-r6.33.gff")
#gff_data <- read.gff(file = "E:/fly_cell_atlas_proj/data/dmel-genes-r6.33.gff")
gff_data$gene <- str_match(gff_data$attributes, 'Name=(.*?);Alias=')[,2]
gff_data <- gff_data %>%
  mutate(gene = str_match(attributes, 'Name=(.*?);')[,2])
gff_data$fb_id <- str_match(gff_data$attributes, 'ID=(.*?);Name=')[,2]
mit_gff_data <- gff_data[gff_data$seqid == "mitochondrion_genome",]

# D.yak
gff_yak <- read.gff(file = "./fly_cell_atlas_proj/data/dyak-genes-r1.05.gff")
gff_yak$gene <- str_match(gff_yak$attributes, 'Name=(.*?);Alias=')[,2]
gff_yak <- gff_yak %>%
  mutate(gene = str_match(attributes, 'Name=(.*?);')[,2])
gff_yak$fb_id <- str_match(gff_yak$attributes, 'ID=(.*?);Name=')[,2]

# D.sim
gff_sim <- read.gff(file = "./fly_cell_atlas_proj/data/dsim-genes-r2.02.gff")
gff_sim$gene <- str_match(gff_sim$attributes, 'Name=(.*?);Alias=')[,2]
gff_sim <- gff_sim %>%
  mutate(gene = str_match(attributes, 'Name=(.*?);')[,2])
gff_sim$fb_id <- str_match(gff_sim$attributes, 'ID=(.*?);Name=')[,2]

#### Change counts matrix rownames to actual gene names ####
# Raw
fb_genes <- left_join(data.frame(fb_id = rownames(raw_dmel)),
                      gff_data[, c('gene', 'fb_id', 'seqid')])

rownames(raw_dmel) <- fb_genes$gene

rm(fb_genes)

# Clean D.mel
fb_genes <- left_join(data.frame(fb_id = rownames(clean_dmel)),
                      gff_data[, c('gene', 'fb_id', 'seqid')])

rownames(clean_dmel) <- fb_genes$gene

rm(fb_genes)

# Clean D.yak
fb_genes <- left_join(data.frame(fb_id = rownames(clean_yak)),
                      gff_yak[, c('gene', 'fb_id', 'seqid')])

rownames(clean_yak) <- fb_genes$gene

rm(fb_genes)

# Clean D.sim
fb_genes <- left_join(data.frame(fb_id = rownames(clean_sim)),
                      gff_sim[, c('gene', 'fb_id', 'seqid')])

rownames(clean_sim) <- fb_genes$gene

rm(fb_genes)

#### Run DecontX on counts data ####
mag_sce <- SingleCellExperiment(assays = list(counts = raw_dmel))
mag_sce <- decontX(x = mag_sce,
                   z = mag_sce$annotation)

mag_sce <- SingleCellExperiment(assays = list(counts = clean_dmel))
#mag_sce <- decontX(x = mag_sce,
#                   z = mag_sce$annotation)

mag_sce_yak <- SingleCellExperiment(assays = list(counts = clean_yak))

mag_sce_sim <- SingleCellExperiment(assays = list(counts = clean_sim))

#### Convert table to Seurat object ####
# Raw
raw_seurat <- CreateSeuratObject(raw_dmel)

metadata_df <- data.frame(sex = 'male',
                          cellID = rownames(raw_seurat@meta.data))

rownames(metadata_df) <- metadata_df$cellID

raw_seurat <- AddMetaData(raw_seurat,
                          metadata_df)

rm(metadata_df)

# Clean Dmel
clean_seurat <- CreateSeuratObject(clean_dmel)

metadata_df <- data.frame(sex = 'male',
                          cellID = rownames(clean_seurat@meta.data))

rownames(metadata_df) <- metadata_df$cellID

clean_seurat <- AddMetaData(clean_seurat,
                            metadata_df)

rm(metadata_df)

# Clean D.yak
clean_yak_seurat <- CreateSeuratObject(clean_yak)

metadata_df <- data.frame(sex = 'male',
                          cellID = rownames(clean_yak_seurat@meta.data))

rownames(metadata_df) <- metadata_df$cellID

clean_yak_seurat <- AddMetaData(clean_yak_seurat,
                            metadata_df)

rm(metadata_df)

# Clean D.sim
clean_sim_seurat <- CreateSeuratObject(clean_sim)

metadata_df <- data.frame(sex = 'male',
                          cellID = rownames(clean_sim_seurat@meta.data))

rownames(metadata_df) <- metadata_df$cellID

clean_sim_seurat <- AddMetaData(clean_sim_seurat,
                                metadata_df)

rm(metadata_df)

#### Add decontX counts to Seurat object ####
raw_seurat[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(mag_sce)))

#### Run stringent filtering on cells ####
## Calculate mitochondrial gene content ##
mit_linked_genes <- as.numeric(which(rownames(raw_seurat@assays$decontXcounts@counts) %in% mit_gff_data[, "gene"]))
mit_linked_data <-  GetAssayData(object = raw_seurat,
                                 assay = 'decontXcounts',
                                 slot = 'counts')[mit_linked_genes,]
mit_linked_data <- as.data.frame(mit_linked_data)
mit_linked_data <- mit_linked_data %>%
  tibble::rownames_to_column(var = "gene") %>%
  gather("cellID", "counts", 2:(ncol(raw_seurat@assays$decontXcounts@counts) + 1))
mit_linked_df <- left_join(mit_linked_data,
                           raw_seurat@meta.data[, c('cellID', 'sex')],
                           by = "cellID")
rm(mit_linked_genes)
rm(mit_linked_data)
gc()

mit_linked_genes <- as.numeric(which(rownames(clean_seurat@assays$RNA@counts) %in% mit_gff_data[, "gene"]))
mit_linked_data <-  GetAssayData(object = clean_seurat,
                                 assay = 'RNA',
                                 slot = 'counts')[mit_linked_genes,]
mit_linked_data <- as.data.frame(mit_linked_data)
mit_linked_data <- mit_linked_data %>%
  tibble::rownames_to_column(var = "gene") %>%
  gather("cellID", "counts", 2:(ncol(clean_seurat@assays$RNA@counts) + 1))
mit_linked_df <- left_join(mit_linked_data,
                           clean_seurat@meta.data[, c('cellID', 'sex')],
                           by = "cellID")
rm(mit_linked_genes)
rm(mit_linked_data)
gc()

#### Run stringent filtering ####
mit_linked_df <- mit_linked_df %>%
  group_by(cellID, sex) %>%
  summarise(total_mit = sum(counts))
mit_linked_df <- left_join(mit_linked_df,
                           data.frame(cellID = names(colSums(raw_seurat@assays$decontXcounts@counts)),
                                      total_counts = colSums(raw_seurat@assays$decontXcounts@counts),
                                      total_genes = colSums(raw_seurat@assays$decontXcounts@counts != 0),
                                      total_ones = colSums(raw_seurat@assays$decontXcounts@counts == 1)))
mit_linked_df$percent.mito <- (mit_linked_df$total_mit/mit_linked_df$total_counts)*100
mit_linked_df <- as.data.frame(mit_linked_df)
rownames(mit_linked_df) <- mit_linked_df$cellID

mit_linked_df <- mit_linked_df %>%
  group_by(cellID, sex) %>%
  summarise(total_mit = sum(counts))
mit_linked_df <- left_join(mit_linked_df,
                           data.frame(cellID = names(colSums(clean_seurat@assays$RNA@counts)),
                                      total_counts = colSums(clean_seurat@assays$RNA@counts),
                                      total_genes = colSums(clean_seurat@assays$RNA@counts != 0),
                                      total_ones = colSums(clean_seurat@assays$RNA@counts == 1)))
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

raw_seurat <- subset(x = raw_seurat,
                     subset = cellID %in% mit_linked_df$cellID,
                     invert = FALSE)
raw_seurat <- AddMetaData(raw_seurat,
                          mit_linked_df[, c('percent.mito',
                                            'total_genes',
                                            'total_ones')])

#clean_seurat <- subset(x = clean_seurat,
#                     subset = cellID %in% mit_linked_df$cellID,
#                     invert = FALSE)
#clean_seurat <- AddMetaData(clean_seurat,
#                          mit_linked_df[, c('percent.mito',
#                                            'total_genes',
#                                            'total_ones')])

rm(median_counts)
rm(median_genes)
rm(mit_linked_df)
gc()

#### Remove genes expressed in less than 3 cells ####
cells_per_gene <- rowSums(raw_seurat@assays$decontXcounts@counts != 0)
cells_per_gene <- data.frame(gene = names(cells_per_gene),
                             ncells = cells_per_gene)
cells_per_gene <- cells_per_gene %>%
  filter(ncells >= 3)

raw_seurat <- subset(x = raw_seurat,
                        features = cells_per_gene$gene)
rm(cells_per_gene)
gc()

#cells_per_gene <- rowSums(clean_seurat@assays$RNA@counts != 0)
#cells_per_gene <- data.frame(gene = names(cells_per_gene),
#                             ncells = cells_per_gene)
#cells_per_gene <- cells_per_gene %>%
#  filter(ncells >= 3)

#clean_seurat <- subset(x = clean_seurat,
#                     features = cells_per_gene$gene)
#rm(cells_per_gene)
#gc()

#### Set decontX as default assay ####
DefaultAssay(object = raw_seurat) <- 'decontXcounts'

DefaultAssay(object = clean_seurat) <- 'RNA'

DefaultAssay(object = clean_yak_seurat) <- 'RNA'

DefaultAssay(object = clean_sim_seurat) <- 'RNA'


#### Rerun Majane et al (2022) code for annotating cell types ####
### Data processing for cell clustering ###
## Normalisation ##
raw_seurat %<>% NormalizeData(normalization.method = "LogNormalize",
                              scale.factor = 10000)

clean_seurat %<>% NormalizeData(normalization.method = "LogNormalize",
                                scale.factor = 10000)

clean_yak_seurat %<>% NormalizeData(normalization.method = "LogNormalize",
                                scale.factor = 10000)

clean_sim_seurat %<>% NormalizeData(normalization.method = "LogNormalize",
                                scale.factor = 10000)

## Find variable genes ##
raw_seurat %<>% FindVariableFeatures(selection.method = "dispersion")

clean_seurat %<>% FindVariableFeatures(selection.method = "dispersion")

clean_yak_seurat %<>% FindVariableFeatures(selection.method = "dispersion")

clean_sim_seurat %<>% FindVariableFeatures(selection.method = "dispersion")

## Scale data to have equal variance ##
raw_seurat %<>% ScaleData(features = rownames(raw_seurat@assays$RNA@counts))

clean_seurat %<>% ScaleData(features = rownames(clean_seurat@assays$RNA@counts))

clean_yak_seurat %<>% ScaleData(features = rownames(clean_yak_seurat@assays$RNA@counts))

clean_sim_seurat %<>% ScaleData(features = rownames(clean_sim_seurat@assays$RNA@counts))

## Run PCA
raw_seurat %<>% RunPCA()

clean_seurat %<>% RunPCA()

clean_yak_seurat %<>% RunPCA()

clean_sim_seurat %<>% RunPCA()

## JackStraw permutation analysis of PC significance
raw_seurat <- JackStraw(object = raw_seurat,
                        dims = 40)
raw_seurat <- ScoreJackStraw(raw_seurat,
                             dims = 1:40)

clean_seurat <- JackStraw(object = clean_seurat,
                        dims = 40)
clean_seurat <- ScoreJackStraw(clean_seurat,
                             dims = 1:40)

clean_yak_seurat <- JackStraw(object = clean_yak_seurat,
                          dims = 40)
clean_yak_seurat <- ScoreJackStraw(clean_yak_seurat,
                               dims = 1:40)

clean_sim_seurat <- JackStraw(object = clean_sim_seurat,
                          dims = 40)
clean_sim_seurat <- ScoreJackStraw(clean_sim_seurat,
                               dims = 1:40)

## Choose PCs to use ##
# Determine percent of variation associated with each PC
pct <- raw_seurat@reductions$pca@stdev / sum(raw_seurat@reductions$pca@stdev) * 100

pct <- clean_seurat@reductions$pca@stdev / sum(clean_seurat@reductions$pca@stdev) * 100

pct_yak <- clean_yak_seurat@reductions$pca@stdev / sum(clean_yak_seurat@reductions$pca@stdev) * 100

pct_sim <- clean_sim_seurat@reductions$pca@stdev / sum(clean_sim_seurat@reductions$pca@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

cumu_yak <- cumsum(pct_yak)

cumu_sim <- cumsum(pct_sim)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1_yak <- which(cumu_yak > 90 & pct < 5)[1]

co1_sim <- which(cumu_sim > 90 & pct < 5)[1]

co1

co1_yak

co1_sim

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

co2_yak <- sort(which((pct_yak[1:length(pct_yak) - 1] - pct_yak[2:length(pct_yak)]) > 0.1), decreasing = T)[1] + 1

co2_sim <- sort(which((pct_sim[1:length(pct_sim) - 1] - pct_sim[2:length(pct_sim)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

co2_yak

co2_sim

pcs <- min(co1, co2)

pcs_yak <- min(co1_yak, co2_yak)

pcs_sim <- min(co1_sim, co2_sim)

# Minimum of the two calculation
pcs

pcs_yak

pcs_sim

JackStrawPlot(object = raw_seurat, 
              dims = 1:40)

JackStrawPlot(object = clean_seurat, 
              dims = 1:40)

JackStrawPlot(object = clean_yak_seurat, 
              dims = 1:40)

JackStrawPlot(object = clean_sim_seurat, 
              dims = 1:40)

use.pcs = c(1:4,6,8,10:11, 13)

use_clean.pcs = c(1:7)

use_yak.pcs = c(1:6)

use_sim.pcs = c(1:13)

## Run UMAP analysis ##
# Raw
raw_seurat <- RunUMAP(raw_seurat,
                      reduction = "pca",
                      dims = use.pcs)

raw_seurat <- FindNeighbors(raw_seurat,
                            reduction = "pca",
                            dims = use.pcs)

# Clean D.mel
clean_seurat <- RunUMAP(clean_seurat,
                      reduction = "pca",
                      dims = use_clean.pcs)

clean_seurat <- FindNeighbors(clean_seurat,
                            reduction = "pca",
                            dims = use_clean.pcs)

# Clean D.yak
clean_yak_seurat <- RunUMAP(clean_yak_seurat,
                        reduction = "pca",
                        dims = use_yak.pcs)

clean_yak_seurat <- FindNeighbors(clean_yak_seurat,
                              reduction = "pca",
                              dims = use_yak.pcs)

# Clean D.sim
clean_sim_seurat <- RunUMAP(clean_sim_seurat,
                            reduction = "pca",
                            dims = use_sim.pcs)

clean_sim_seurat <- FindNeighbors(clean_sim_seurat,
                                  reduction = "pca",
                                  dims = use_sim.pcs)

## Find clusters at multiple resolutions ##
raw_seurat <- FindClusters(object = raw_seurat,
                           resolution = seq(0.2, 1.4, 0.1),
                           verbose = FALSE)

clean_seurat <- FindClusters(object = clean_seurat,
                           resolution = seq(0.2, 1.4, 0.1),
                           verbose = FALSE)

clean_yak_seurat <- FindClusters(object = clean_yak_seurat,
                             resolution = seq(0.2, 1.4, 0.1),
                             verbose = FALSE)

clean_sim_seurat <- FindClusters(object = clean_sim_seurat,
                             resolution = seq(0.2, 1.4, 0.1),
                             verbose = FALSE)

## Pick appropriate resolution ##
# Raw
sapply(grep("res",colnames(raw_seurat@meta.data),value = TRUE),
       function(x) length(unique(raw_seurat@meta.data[,x])))

Idents(raw_seurat) <- "decontXcounts_snn_res.0.2"

DimPlot(raw_seurat, reduction = "umap", label = T, pt.size = 2) +
  theme_minimal(base_size = 15)

# Clean D.mel
sapply(grep("res",colnames(clean_seurat@meta.data),value = TRUE),
       function(x) length(unique(clean_seurat@meta.data[,x])))

Idents(clean_seurat) <- "RNA_snn_res.0.2"

DimPlot(clean_seurat, reduction = "umap", label = T, pt.size = 2) +
  theme_minimal(base_size = 15)

# Clean D.yak
sapply(grep("res",colnames(clean_yak_seurat@meta.data),value = TRUE),
       function(x) length(unique(clean_yak_seurat@meta.data[,x])))

Idents(clean_yak_seurat) <- "RNA_snn_res.0.2"

DimPlot(clean_yak_seurat, reduction = "umap", label = T, pt.size = 2) +
  theme_minimal(base_size = 15)

# Clean D.sim
sapply(grep("res",colnames(clean_sim_seurat@meta.data),value = TRUE),
       function(x) length(unique(clean_sim_seurat@meta.data[,x])))

Idents(clean_sim_seurat) <- "RNA_snn_res.0.2"

DimPlot(clean_sim_seurat, reduction = "umap", label = T, pt.size = 2) +
  theme_minimal(base_size = 15)


FeaturePlot(raw_seurat,
            features = c("SP", "Acp36DE", "Acp95EF", "lectin-46Ca",
                         "lectin-46Cb", "abd-A", "lncRNA:iab8", "vvl",
                         "Dup99B"))

FeaturePlot(clean_seurat,
            features = c("SP", "Acp36DE", "Acp95EF", "lectin-46Ca",
                         "lectin-46Cb", "abd-A", "lncRNA:iab8", "vvl",
                         "Dup99B"))

## Assign cell types ##
raw_seurat %<>% RenameIdents('0' = 'MCsp1',
                             '1' = 'MCsp2',
                             '2' = 'EDC',
                             '3' = 'SC')

clean_seurat %<>% RenameIdents('0' = 'MCsp1',
                             '1' = 'MCsp2',
                             '2' = 'EDC',
                             '3' = 'SC')

clean_yak_seurat %<>% RenameIdents('0' = 'MCsp1',
                               '1' = 'MCsp2',
                               '2' = 'EDC',
                               '3' = 'SC')

clean_sim_seurat %<>% RenameIdents('0' = 'MCsp1',
                               '1' = 'MCsp2',
                               '2' = 'EDC',
                               '3' = 'SC')




## Add annotation to metadata ##
# Raw
annotation <- data.frame(cell_type_id = Idents(object = raw_seurat))
rownames(annotation) <- colnames(x = raw_seurat)

annotation <- annotation %>%
  mutate(annotation = case_when((cell_type_id == "MCsp1") ~ "main cell 1",
                                (cell_type_id == "MCsp2") ~ "main cell 2",
                                (cell_type_id == "SC") ~ "secondary cell",
                                (cell_type_id == "EDC") ~ "ejaculatory duct cell"))

raw_seurat <- AddMetaData(object = raw_seurat,
                          metadata = annotation$annotation,
                          col.name = 'annotation')

rm(annotation)

# Clean D.mel
annotation <- data.frame(cell_type_id = Idents(object = clean_seurat))
rownames(annotation) <- colnames(x = clean_seurat)

annotation <- annotation %>%
  mutate(annotation = case_when((cell_type_id == "MCsp1") ~ "main cell 1",
                                (cell_type_id == "MCsp2") ~ "main cell 2",
                                (cell_type_id == "SC") ~ "secondary cell",
                                (cell_type_id == "EDC") ~ "ejaculatory duct cell"))

clean_seurat <- AddMetaData(object = clean_seurat,
                          metadata = annotation$annotation,
                          col.name = 'annotation')

rm(annotation)

# Clean D.yak
annotation <- data.frame(cell_type_id = Idents(object = clean_yak_seurat))
rownames(annotation) <- colnames(x = clean_yak_seurat)

annotation <- annotation %>%
  mutate(annotation = case_when((cell_type_id == "MCsp1") ~ "main cell 1",
                                (cell_type_id == "MCsp2") ~ "main cell 2",
                                (cell_type_id == "SC") ~ "secondary cell",
                                (cell_type_id == "EDC") ~ "ejaculatory duct cell"))

clean_yak_seurat <- AddMetaData(object = clean_yak_seurat,
                            metadata = annotation$annotation,
                            col.name = 'annotation')

rm(annotation)

# Clean D.sim
annotation <- data.frame(cell_type_id = Idents(object = clean_sim_seurat))
rownames(annotation) <- colnames(x = clean_sim_seurat)

annotation <- annotation %>%
  mutate(annotation = case_when((cell_type_id == "MCsp1") ~ "main cell 1",
                                (cell_type_id == "MCsp2") ~ "main cell 2",
                                (cell_type_id == "SC") ~ "secondary cell",
                                (cell_type_id == "EDC") ~ "ejaculatory duct cell"))

clean_sim_seurat <- AddMetaData(object = clean_sim_seurat,
                                metadata = annotation$annotation,
                                col.name = 'annotation')

rm(annotation)



## Run scTransform on decontX assay data ##
raw_seurat <- SCTransform(object = raw_seurat,
                          assay = 'decontXcounts',
                          vst.flavor = 'v2',
                          vars.to.regress = c('nCount_decontXcounts',
                                              'percent.mito'),
                          verbose = TRUE)

#clean_seurat <- SCTransform(object = clean_seurat,
#                          assay = 'RNA',
#                          vst.flavor = 'v2',
#                          vars.to.regress = c('nCount_RNA',
#                                              'percent.mito'),
#                          verbose = TRUE)


## Save Seurat object to file ##
SaveH5Seurat(object = raw_seurat,
             filename = paste0("majane_etal_dmel_male_reproductive_glands_decontx_stringent_filt_sctransformed.h5Seurat"))

#SaveH5Seurat(object = clean_seurat,
#             filename = paste0("majane_etal_dmel_male_reproductive_glands_clean_decontx_stringent_filt_sctransformed.h5Seurat"))

#### Run X:A analysis ####
### Calculate X/A expression ratio function ###
x_a_ratio_n_max <- function(cell, seurat_file, dmel_genes) {
  #### Extract genes that have expression > 0 ####
  expressed_data <- as.data.frame(seurat_file@assays$RNA@counts[which(as.data.frame(seurat_file@assays$RNA@counts[,cell]) > 0), cell])
  expressed_genes <- data.frame(expression = expressed_data[,1],
                                gene = as.character(rownames(expressed_data)),
                                cell_id = as.character(colnames(seurat_file@assays$RNA[, cell])),
                                sex = 'male',
                                cell_type = as.character(seurat_file$annotation[cell]))
  rm(expressed_data)
  
  expressed_genes <- merge(expressed_genes,
                           dmel_genes,
                           by = "gene")
  
  #### Calculate X/A expression ratio for each cell ####
  data_for_ratio <- expressed_genes %>%
    group_by(status) %>%
    summarise(total_exp = sum(expression),
              median_exp = median(expression),
              max_exp = quantile(expression, probs = 0.95),
              expressed_counts = n())
  if ("msl-1" %in% expressed_genes$gene){
    msl1_expression <- expressed_genes[expressed_genes$gene == "msl-1", "expression"]
  } else {
    msl1_expression <- 0
  }
  if ("msl-2" %in% expressed_genes$gene){
    msl2_expression <- expressed_genes[expressed_genes$gene == "msl-2", "expression"]
  } else {
    msl2_expression <- 0
  }
  if ("msl-3" %in% expressed_genes$gene){
    msl3_expression <- expressed_genes[expressed_genes$gene == "msl-3", "expression"]
  } else {
    msl3_expression <- 0
  }
  if ("mof" %in% expressed_genes$gene){
    mof_expression <- expressed_genes[expressed_genes$gene == "mof", "expression"]
  } else {
    mof_expression <- 0
  }
  if ("mle" %in% expressed_genes$gene){
    mle_expression <- expressed_genes[expressed_genes$gene == "mle", "expression"]
  } else {
    mle_expression <- 0
  }
  if ("lncRNA:roX1" %in% expressed_genes$gene){
    rox1_expression <- expressed_genes[expressed_genes$gene == "lncRNA:roX1", "expression"]
  } else {
    rox1_expression <- 0
  }
  if ("lncRNA:roX2" %in% expressed_genes$gene){
    rox2_expression <- expressed_genes[expressed_genes$gene == "lncRNA:roX2", "expression"]
  } else {
    rox2_expression <- 0
  }
  if ("dsx" %in% expressed_genes$gene){
    dsx_expression <- expressed_genes[expressed_genes$gene == "dsx", "expression"]
  } else {
    dsx_expression <- 0
  }
  if ("X" %in% data_for_ratio$status && "autosome" %in% data_for_ratio$status) {
    max_x <- as.numeric(data_for_ratio[data_for_ratio$status == "X", "max_exp"])
    counts_x <- as.numeric(data_for_ratio[data_for_ratio$status == "X", "expressed_counts"])
    max_autosomes <- as.numeric(data_for_ratio[data_for_ratio$status == "autosome", "max_exp"])
    counts_autosomes <- as.numeric(data_for_ratio[data_for_ratio$status == "autosome", "expressed_counts"])
    data_for_ratio <- data_for_ratio %>%
      summarise(total_x = data_for_ratio[status == "X", "total_exp"],
                total_aut = data_for_ratio[status == "autosome", "total_exp"],
                total_ratio = data_for_ratio[status == "X", "total_exp"]/data_for_ratio[status == "autosome", "total_exp"],
                median_ratio = data_for_ratio[status == "X", "median_exp"]/data_for_ratio[status == "autosome", "median_exp"],
                normalised_ratio = ((data_for_ratio[status == "X",
                                                    "total_exp"]/data_for_ratio[status == "X",
                                                                                "expressed_counts"])/(data_for_ratio[status == "autosome",
                                                                                                                     "total_exp"]/data_for_ratio[status == "autosome",
                                                                                                                                                 "expressed_counts"])))
    cell_ratio_df <- data.frame(cell_id = unique(expressed_genes$cell_id),
                                cell_type = unique(expressed_genes$cell_type),
                                sex = unique(expressed_genes$sex),
                                tissue = 'male_reproductive_glands',
                                total_x = as.numeric(data_for_ratio$total_x),
                                total_aut = as.numeric(data_for_ratio$total_aut),
                                total_ratio = as.numeric(data_for_ratio$total_ratio),
                                median_ratio = as.numeric(data_for_ratio$median_ratio),
                                normalised_ratio = as.numeric(data_for_ratio$normalised_ratio),
                                max_x = max_x,
                                max_autosomes = max_autosomes,
                                msl1_expression = msl1_expression,
                                msl2_expression = msl2_expression,
                                msl3_expression = msl3_expression,
                                mof_expression = mof_expression,
                                mle_expression = mle_expression,
                                rox1_expression = rox1_expression,
                                rox2_expression = rox2_expression,
                                dsx_expression = dsx_expression,
                                counts_x = counts_x,
                                counts_autosomes = counts_autosomes)
    
  } else if (!("X" %in% data_for_ratio$status)) {
    if (!("autosome" %in% data_for_ratio$status)) {
      max_autosomes <- NA
      counts_autosomes <- NA
      total_aut <- 0
    } else {
      max_autosomes <- as.numeric(data_for_ratio[data_for_ratio$status == "autosome", "max_exp"])
      counts_autosomes <- as.numeric(data_for_ratio[data_for_ratio$status == "autosome", "counts"])
      total_aut <- as.numeric(data_for_ratio[data_for_ratio$status == "autosome", "total_exp"])
    }
    cell_ratio_df <- data.frame(cell_id = unique(expressed_genes$cell_id),
                                cell_type = unique(expressed_genes$cell_type),
                                sex = unique(expressed_genes$sex),
                                tissue = 'male_reproductive_glands',
                                total_x = 0,
                                total_aut = total_aut,
                                total_ratio = 0,
                                median_ratio = 0,
                                normalised_ratio = 0,
                                max_x = NA,
                                max_autosomes = max_autosomes,
                                msl1_expression = msl1_expression,
                                msl2_expression = msl2_expression,
                                msl3_expression = msl3_expression,
                                mof_expression = mof_expression,
                                mle_expression = mle_expression,
                                rox1_expression = rox1_expression,
                                rox2_expression = rox2_expression,
                                dsx_expression = dsx_expression,
                                counts_x = 0,
                                counts_autosomes = counts_autosomes)
  }
  return(cell_ratio_df)
}


#### Find genes in GFF data ####
# Raw
dmel_genes <- as.data.frame(subset(gff_data, gff_data$gene 
                                     %in%
                                       rownames(raw_seurat@assays$SCT@counts))[,c('seqid', 'gene')])
colnames(dmel_genes) <- c('chromosome', 'gene')

# Clean D.mel
dmel_genes <- as.data.frame(subset(gff_data, gff_data$gene 
                                   %in%
                                     rownames(clean_seurat@assays$RNA@counts))[,c('seqid', 'gene')])
colnames(dmel_genes) <- c('chromosome', 'gene')

# Clean D.yak
dmel_genes <- as.data.frame(subset(gff_yak, gff_yak$gene 
                                   %in%
                                     rownames(clean_yak_seurat@assays$RNA@counts))[,c('seqid', 'gene')])
colnames(dmel_genes) <- c('chromosome', 'gene')

# Clean D.sim
dmel_genes <- as.data.frame(subset(gff_sim, gff_sim$gene 
                                   %in%
                                     rownames(clean_sim_seurat@assays$RNA@counts))[,c('seqid', 'gene')])
colnames(dmel_genes) <- c('chromosome', 'gene')

dmel_genes$chromosome <- gsub('Scf_', '',
                              dmel_genes$chromosome)

#### Set chromosome status ####
dmel_genes <- dmel_genes %>%
  mutate(status = case_when((chromosome == "X") ~ "X",
                            (chromosome == "Y") ~ "Y",
                            (chromosome == "mitochondrion_genome") ~ "mt",
                            (chromosome == "2L") ~ "autosome",
                            (chromosome == "2R") ~ "autosome",
                            (chromosome == "3L") ~ "autosome",
                            (chromosome == "3R") ~ "autosome",
                            (chromosome == "4") ~ "autosome",
                            TRUE ~ "other"))

for (cell in 1:dim(raw_seurat@assays$SCT@counts)[2]) {
  cell_data <- c()
  cell_data <- x_a_ratio_n_max(cell = cell,
                               seurat_file = raw_seurat)
  if (class(cell_data) == "NULL") {
    next
  }
  
  if (cell == 1){
    write.table(cell_data,
                paste0("x_a_normalised_expression_ratio_n_max_expression_data_majane_etal.csv"),
                row.names = FALSE,
                col.names = TRUE)
  } else {
    write.table(cell_data,
                paste0("x_a_normalised_expression_ratio_n_max_expression_data_majane_etal.csv"),
                row.names = FALSE,
                col.names = FALSE,
                append = TRUE)
  }
}


# NB: Replace with clean_seurat/clean_yak_seurat/clean_sim_seurat below
for (cell in 1:dim(clean_sim_seurat@assays$RNA@counts)[2]) {
  cell_data <- c()
  cell_data <- x_a_ratio_n_max(cell = cell,
                               seurat_file = clean_sim_seurat,
                               dmel_genes = dmel_genes) 
  if (class(cell_data) == "NULL") {
    next
  }
  
  if (cell == 1){
    write.table(cell_data,
                paste0("x_a_normalised_expression_ratio_n_max_og_clean_sim_expression_data_majane_etal.csv"),
                row.names = FALSE,
                col.names = TRUE)
  } else {
    write.table(cell_data,
                paste0("x_a_normalised_expression_ratio_n_max_og_clean_sim_expression_data_majane_etal.csv"),
                row.names = FALSE,
                col.names = FALSE,
                append = TRUE)
  }
}


#### Plot X:A data ####
## Read in csv file ##
dmel_xtoa <- read.table("x_a_normalised_expression_ratio_n_max_expression_data_majane_etal.csv",
                        header = TRUE)

rownames(dmel_xtoa) <- dmel_xtoa$cell_id

# Clean D.mel
clean_xtoa <- read.table("x_a_normalised_expression_ratio_n_max_og_clean_expression_data_majane_etal.csv",
                         header = TRUE)

rownames(clean_xtoa) <- clean_xtoa$cell_id

# Clean D.yak
yak_xtoa <- read.table("x_a_normalised_expression_ratio_n_max_og_clean_yak_expression_data_majane_etal.csv",
                        header = TRUE)

rownames(yak_xtoa) <- yak_xtoa$cell_id

# Clean D.sim
sim_xtoa <- read.table("x_a_normalised_expression_ratio_n_max_og_clean_sim_expression_data_majane_etal.csv",
                       header = TRUE)

rownames(sim_xtoa) <- sim_xtoa$cell_id

## Read in h5Seurat file ##
dmel_seurat <- LoadH5Seurat("majane_etal_dmel_male_reproductive_glands_decontx_stringent_filt_sctransformed.h5Seurat")

clean_seurat <- LoadH5Seurat("majane_etal_dmel_male_reproductive_glands_clean_decontx_stringent_filt_sctransformed.h5Seurat")

## Add ratio info as metadata ##
# Add total_ratio
dmel_seurat <- AddMetaData(object = dmel_seurat,
                           metadata = dmel_xtoa$total_ratio,
                           col.name = 'total_ratio')

clean_seurat <- AddMetaData(object = clean_seurat,
                           metadata = clean_xtoa$total_ratio,
                           col.name = 'total_ratio')

# Add normalised_ratio
dmel_seurat <- AddMetaData(object = dmel_seurat,
                           metadata = dmel_xtoa$normalised_ratio,
                           col.name = 'normalised_ratio')

clean_seurat <- AddMetaData(object = clean_seurat,
                           metadata = clean_xtoa$normalised_ratio,
                           col.name = 'normalised_ratio')

clean_yak_seurat <- AddMetaData(object = clean_yak_seurat,
                            metadata = yak_xtoa$normalised_ratio,
                            col.name = 'normalised_ratio')

clean_sim_seurat <- AddMetaData(object = clean_sim_seurat,
                            metadata = sim_xtoa$normalised_ratio,
                            col.name = 'normalised_ratio')

# Add rox1 expression
dmel_seurat <- AddMetaData(object = dmel_seurat,
                           metadata = dmel_xtoa$rox1_expression,
                           col.name = 'rox1_expression')

clean_seurat <- AddMetaData(object = clean_seurat,
                           metadata = clean_xtoa$rox1_expression,
                           col.name = 'rox1_expression')

# Add rox2 expression
dmel_seurat <- AddMetaData(object = dmel_seurat,
                           metadata = dmel_xtoa$rox2_expression,
                           col.name = 'rox2_expression')

clean_seurat <- AddMetaData(object = clean_seurat,
                           metadata = clean_xtoa$rox2_expression,
                           col.name = 'rox2_expression')

## Change annotation for clarity ##
Idents(dmel_seurat) <- 'annotation'
Idents(clean_seurat) <- 'annotation'
Idents(clean_yak_seurat) <- 'annotation'
Idents(clean_sim_seurat) <- 'annotation'

dmel_seurat@meta.data$annotation <- gsub(pattern = 'ejaculatory duct cell',
                                           replacement = 'ejaculatory duct',
                                           dmel_seurat@meta.data$annotation)

clean_seurat@meta.data$annotation <- gsub(pattern = 'ejaculatory duct cell',
                                         replacement = 'ejaculatory duct',
                                         clean_seurat@meta.data$annotation)

clean_yak_seurat@meta.data$annotation <- gsub(pattern = 'ejaculatory duct cell',
                                          replacement = 'ejaculatory duct',
                                          clean_yak_seurat@meta.data$annotation)

clean_sim_seurat@meta.data$annotation <- gsub(pattern = 'ejaculatory duct cell',
                                          replacement = 'ejaculatory duct',
                                          clean_sim_seurat@meta.data$annotation)

dmel_xtoa <- FeaturePlot(dmel_seurat,
                        features = 'normalised_ratio',
                        reduction = "umap",
                        label = TRUE,
                        repel = TRUE,
                        label.size = 2.5,
                        pt.size = 1.2,
                        order = TRUE,
                        min.cutoff = 0,
                        max.cutoff = 1.1,
                        cols = c('#6b008a', '#E69F00')) +
  ggtitle("mean X:A") +
  labs(x = "UMAP 1",
       y = "UMAP 2") +
  theme_minimal() +
  theme(text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5,
                                  face = 'bold'))

clean_xtoa_sim <- FeaturePlot(clean_sim_seurat,
                         features = 'normalised_ratio',
                         reduction = "umap",
                         label = TRUE,
                         repel = TRUE,
                         label.size = 2.5,
                         pt.size = 1.2,
                         order = TRUE,
                         min.cutoff = 0,
                         max.cutoff = 1.1,
                         cols = c('#6b008a', '#E69F00')) +
  ggtitle("mean X:A") +
  labs(x = "UMAP 1",
       y = "UMAP 2") +
  theme_minimal() +
  theme(text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5,
                                  face = 'bold'))

# Clean D.mel
x_a_data <- data.frame(cellID = clean_seurat@meta.data$cellID,
                       normalised_ratio = clean_seurat@meta.data$normalised_ratio,
                       annotation = clean_seurat@meta.data$annotation)

# Clean D.yak
x_a_data <- data.frame(cellID = clean_yak_seurat@meta.data$cellID,
                       normalised_ratio = clean_yak_seurat@meta.data$normalised_ratio,
                       annotation = clean_yak_seurat@meta.data$annotation)

# Clean D.sim
x_a_data <- data.frame(cellID = clean_sim_seurat@meta.data$cellID,
                       normalised_ratio = clean_sim_seurat@meta.data$normalised_ratio,
                       annotation = clean_sim_seurat@meta.data$annotation)

x_a_data$annotation <- factor(x_a_data$annotation,
                                levels = c('main cell 1',
                                           'main cell 2',
                                           'secondary cell',
                                           'ejaculatory duct'))

x_a_normalised_plot_sim <- ggplot(data = x_a_data,
                              aes(x = normalised_ratio,
                                  y = annotation,
                                  colour = annotation,
                                  fill = annotation)) + 
  ggridges::geom_density_ridges(jittered_points = TRUE,
                                position = position_points_jitter(width = 0.05, height = 0),
                                point_shape = '|', 
                                point_size = 3,
                                point_alpha = 1,
                                alpha = 0.3,
                                scale = 1) +
  labs(x = "X:A mean expression",
       y = "") +
  scale_colour_manual(values = c('#e69f00', '#f0e442',  '#9786e6', '#288e60')) +
  scale_fill_manual(values = c('#e69f00', '#f0e442', '#9786e6', '#288e60')) +
  lims(x = c(0.3, 1.1)) +
  theme(legend.position = "none",
        text = element_text(size = 11),
        panel.grid.major = element_line(color = "#999999",
                                        linetype = 'dashed'),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Clean D.mel
clean_seurat <- ScaleData(clean_seurat,
                          features = rownames(clean_seurat))

clean_heatmeap <- DoHeatmap(clean_seurat,
                            features = c("SP", "Acp36DE", "Acp95EF", "lectin-46Ca",
                                         "lectin-46Cb", "abd-A", "lncRNA:iab8",
                                         "vvl","Dup99B"),
                            slot = "scale.data",
                            group.colors = c('#e69f00', '#9786e6', '#288e60',
                                             '#abd20c', '#905395'),
                            size = 3)

# Clean D.yak
clean_yak_seurat <- ScaleData(clean_yak_seurat,
                          features = rownames(clean_yak_seurat))

clean_heatmeap <- DoHeatmap(clean_yak_seurat,
                            features = c("Dyak\\GE21959", "Dyak\\GE18098",
                                         "Dyak\\GE19302",
                                         "Dyak\\GE19303", "Dyak\\GE24883",
                                         "Dyak\\GE16054",
                                         "Dyak\\GE10438"),
                            slot = "scale.data",
                            group.colors = c('#e69f00', '#9786e6', '#288e60',
                                             '#abd20c', '#905395'),
                            size = 3)

# Clean D.sim
clean_sim_seurat <- ScaleData(clean_sim_seurat,
                              features = rownames(clean_sim_seurat))

clean_heatmeap <- DoHeatmap(clean_sim_seurat,
                            features = c("Dsim\\Acp70A", "Dsim\\Acp36DE",
                                         "Dsim\\lectin-46Ca",
                                         "Dsim\\GD10684", "Dsim\\GD20268",
                                         "Dsim\\GD22109",
                                         "Dsim\\GD17687"),
                            slot = "scale.data",
                            group.colors = c('#e69f00', '#9786e6', '#288e60',
                                             '#abd20c', '#905395'),
                            size = 3)

x_to_a_combi <- ggarrange(ggarrange(x_a_normalised_plot,
                                   clean_xtoa,
                                   labels = c('(A)', '(B)'),
                                   nrow = 1),
                          clean_heatmeap,
                          labels = c('', '(C)'),
                          nrow = 2)


ggsave(file = paste0("majane_og_clean_xtoa_ridges_umap_and_heatmap.pdf"),
       plot = x_to_a_combi,
       width = 210, height = 210, units = c("mm"), dpi = 300)

title_plot_yak <- as_ggplot(text_grob("Drosophila yakuba",
                                   color = '#999999',
                                   face = "bold.italic",
                                   size = 12,
                                   just = c(0.5,1))) + 
  theme(plot.margin = margin(0, 0, 0.5, 0, "cm"))

title_plot_sim <- as_ggplot(text_grob("Drosophila simulans",
                                      color = '#999999',
                                      face = "bold.italic",
                                      size = 12,
                                      just = c(0.5,1))) + 
  theme(plot.margin = margin(0, 0, 0.5, 0, "cm"))


yak_sim_combi <- ggarrange(title_plot_sim,
                           ggarrange(x_a_normalised_plot_sim,
                                     clean_xtoa_sim,
                                     labels = c('(A)', '(B)'),
                                     nrow = 1, ncol = 2),
                           title_plot_yak,
                           ggarrange(x_a_normalised_plot_yak,
                                     clean_xtoa_yak,
                                     labels =  c('(C)', '(D)'),
                                     nrow = 1, ncol = 2),
                           nrow = 4,
                           heights = c(0.15, 1, 0.15, 1))

ggsave(file = paste0("majane_og_clean_sim_yak_xtoa_ridges_umap.pdf"),
       plot = yak_sim_combi,
       width = 180, height = 160, units = c("mm"), dpi = 300)



dmel_dcc <- FeaturePlot(clean_seurat,
                        features = c('lncRNA:roX1',
                                     'lncRNA:roX2',
                                     'mle',
                                     'mof',
                                     'msl-1',
                                     'msl-2',
                                     'msl-3',
                                     'dsx'),
                        reduction = "umap",
                        slot = 'scale.data',
                        label = TRUE,
                        repel = TRUE,
                        label.size = 2.5,
                        pt.size = 1.2,
                        order = TRUE,
                        keep.scale = 'all',
                        cols = c('#6b008a', '#E69F00'))

dmel_other <- FeaturePlot(dmel_seurat,
                        features = c('Mtor',
                                     'JIL-1',
                                     'Top2',
                                     'Clamp'),
                        reduction = "umap",
                        slot = 'counts',
                        label = TRUE,
                        repel = TRUE,
                        label.size = 2.5,
                        pt.size = 1.2,
                        order = TRUE,
                        keep.scale = 'all',
                        cols = c('#6b008a', '#E69F00'))


#### Find highly expressed genes ####
## Subset each cell type ##
dmel_sc <- subset(dmel_seurat,
                  subset = annotation == 'secondary cell')

dmel_ed <- subset(dmel_seurat,
                  subset = annotation == 'ejaculatory duct')

dmel_mc1 <- subset(dmel_seurat,
                   subset = annotation == 'main cell 1')

dmel_mc2 <- subset(dmel_seurat,
                   subset = annotation == 'main cell 2')

## Calculate average gene expression across cell types ##
dmel_sc_mean_counts <- rowMeans(dmel_sc@assays$SCT@counts)

dmel_ed_mean_counts <- rowMeans(dmel_ed@assays$SCT@counts)

dmel_mc1_mean_counts <- rowMeans(dmel_mc1@assays$SCT@counts)

dmel_mc2_mean_counts <- rowMeans(dmel_mc2@assays$SCT@counts)

## Combine data ##
dmel_mean_counts <- data.frame(gene = names(dmel_sc_mean_counts),
                               main_cell_1 = dmel_mc1_mean_counts,
                               main_cell_2 = dmel_mc2_mean_counts,
                               secondary_cell = dmel_sc_mean_counts,
                               ejaculatory_duct = dmel_ed_mean_counts)
library(data.table)
dmel_mean_counts <- melt(dmel_mean_counts,
                         id.vars = 'gene')

rm(list = c('dmel_sc_mean_counts', 'dmel_ed_mean_counts',
            'dmel_mc1_mean_counts', 'dmel_mc2_mean_counts'))

## Find 20 most expressed genes per cell type ##
dmel_mean_counts %>%
  group_by(variable) %>%
  top_n(20, value)



### Read in Majane et al main cell 1 and 2 marker list ###
markers <- read.csv("G:/majane_2022_data/majane_etal_main_cell_markers.csv",
                    header = TRUE)

markers <- markers %>%
  mutate(marker_status = case_when(avg_logFC > 0 ~ "main cell 2",
                                   avg_logFC < 0 ~ "main cell 1"),
         significance = case_when(abs(avg_logFC) > 1 ~ 'significant',
                                  .default = 'n.significant'))

DoHeatmap(subset(dmel_seurat, annotation %in% c('main cell 1', 'main cell 2')),
          features = subset(markers, marker_status == 'main cell 2')$gene)

DoHeatmap(subset(dmel_seurat, annotation %in% c('main cell 1', 'main cell 2')),
          features = subset(markers, significance == 'significant')$gene)


#### Produce dotplot with expression of DCC genes ####
genes <- c("lncRNA:roX1", "lncRNA:roX2", "msl-1", "msl-2", "msl-3", "mle", "mof",
           "dsx")

Idents(dmel_seurat) <- 'annotation'

dcc_dotplot <- DotPlot(object = dmel_seurat,
                       features = genes,
                       cols = c('#6b008a', '#E69F00'),
                       scale = FALSE) &
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        text = element_text(size = 11))

ggsave(file = 'majane_mag_clusters_dcc_genes_dotplot_unscaled_data.pdf',
       plot = dcc_dotplot,
       width = 240, height = 90, units = c("mm"), dpi = 300)
