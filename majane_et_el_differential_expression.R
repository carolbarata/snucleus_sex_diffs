###########################################
#     Analysing Majane et al's (2022)     #
#       dataset for X:A differences       #
#                 in Dmel                 #
#                continued                #
###########################################

#### Load libraries ####
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(ggpubr)
library(ggridges)

#### Read in h5Seurat file ####
dmel_seurat <- LoadH5Seurat("majane_etal_dmel_male_reproductive_glands_decontx_stringent_filt_sctransformed.h5Seurat")

#clean_seurat <- LoadH5Seurat("majane_etal_dmel_male_reproductive_glands_clean_decontx_stringent_filt_sctransformed.h5Seurat")


#### Perform DE between main cell clusters ####
#mag_subsubset <- subset(x = clean_seurat,
#                        subset = (annotation %in% c('main cell 1', 'main cell 2')))
mag_subsubset <- subset(x = dmel_seurat,
                        subset = (annotation %in% c('main cell 1', 'main cell 2')))


### Get counts for DCC genes ###
genes <- c("lncRNA:roX1", "lncRNA:roX2", "msl-1", "msl-2", "msl-3", "mle",
           "dsx", "Mtor", "JIL-1", "Top2", "Clamp")

gene_counts <- GetAssayData(object = mag_subsubset,
                            slot = 'counts',
                            layer = 'SCT')[genes,]

gene_counts <- as.matrix(gene_counts)

gene_counts <- t(gene_counts)

gene_counts <- as.data.frame(gene_counts)

gene_counts$cellID <- rownames(gene_counts)

gene_counts <- left_join(gene_counts, mag_subsubset@meta.data[,c('cellID', 'annotation')],
                         by = c("cellID" = "cellID"))

gene_counts %>%
  group_by(annotation) %>%
  summarise(median = median(`lncRNA:roX1`))


### MAST analysis ###
library(SingleCellExperiment)
library(MAST)
mag_sce <- SingleCellExperiment(list(logcounts = GetAssayData(mag_subsubset,
                                                              layer = "data",
                                                              assay = "SCT"),
                                     counts = GetAssayData(mag_subsubset,
                                                           layer = "counts",
                                                           assay = "SCT")),
                                colData = mag_subsubset@meta.data)

mag_sce_filt <- mag_sce[rowSums(SingleCellExperiment::counts(mag_sce) > 1) >= 10,]
mag_sca <- SceToSingleCellAssay(mag_sce_filt)
print("Dimensions before subsetting:")
print(dim(mag_sca))
print("")
# keep genes that are expressed in more than 10% of all cells
mag_sca_filt <- mag_sca[freq(mag_sca) > 0.1,]
print("Dimensions after subsetting:")
print(dim(mag_sca_filt))
print("")
# add a column to the data which contains scaled number of genes that are expressed in each cell
cdr2 <- colSums(assay(mag_sca) > 0)
colData(mag_sca)$ngeneson <- scale(cdr2)
# store the columns that we are interested in as factors
annotation <- factor(colData(mag_sca)$annotation)
# set the reference level
annotation <- relevel(annotation,"main cell 1")
# define and fit the model
zlmCond <- zlm(formula = ~ngeneson + annotation, 
               sca = mag_sca)
# perform likelihood-ratio test for the condition that we are interested in    
summaryCond <- summary(zlmCond, doLRT = 'annotationmain cell 2')
# get the table with log-fold changes and p-values
summaryDt <- summaryCond$datatable
result_majane <- merge(summaryDt[contrast=='annotationmain cell 2' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                summaryDt[contrast=='annotationmain cell 2' & component=='logFC', .(primerid, coef)],
                by='primerid') # logFC coefficients
# MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
result_majane[,coef:=result_majane[,coef]/log(2)]
# do multiple testing correction
result_majane[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
result_majane = result_majane[result_majane$FDR<0.01,, drop=F]

result_majane <- stats::na.omit(as.data.frame(result_majane))
write.table(result_majane, "majane_dmel_mast_result_mag_main_cell_clusters.txt",
            row.names = FALSE)

result_majane <- read.table("majane_dmel_mast_result_mag_main_cell_clusters.txt",
                     header = TRUE)

result_majane$log10FDR <- -log10(result_majane$FDR)

result_majane <- result_majane %>%
  mutate(status = case_when(abs(coef) >= 1 ~ 'DE',
                            abs(coef) < 1 ~ 'unbiased'))
