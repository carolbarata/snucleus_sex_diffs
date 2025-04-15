##############################################
#             Investigating male             #
#             reproductive glands            #
#             inc accessory gland            #
#            and ejaculatory duct            #
##############################################

#### Load libraries ####
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(ggpubr)
library(ggridges)


#### Read in csv file ####
mag_xtoa <- read.table("x_a_normalised_expression_ratio_n_max_expression_data_male_reproductive_glands_relaxed_after_filtering.csv",
                       header = TRUE)

rownames(mag_xtoa) <- mag_xtoa$cell_id

### Read in h5Seurat file ###
mag_relaxed_seurat <- LoadH5Seurat("r_male_reproductive_glands_cell_type_filt_decontx_stringent_filt_sctransformed.h5Seurat")

### Read in anndata stringent data ###
mag_anndata <- read_h5ad('./fly_cell_atlas_proj/data/h5ad_initial_files/s_fca_biohub_male_reproductive_glands_10x.h5ad',
                         backed = 'r')

## Extract UMAP coordinates ##
mag_umap <- as.data.frame(mag_anndata$obsm[['X_umap']])

rownames(mag_umap) <- rownames(mag_anndata$obs)

mag_umap$cellID <- rownames(mag_umap)

# Only keep cells that are present in Seurat file
mag_umap <- subset(mag_umap,
                   cellID %in% mag_relaxed_seurat@meta.data$cellID)

## Extract metadata ##
mag_metadata <- as.data.frame(mag_anndata$obs)

rm(mag_anndata)
gc()

mag_metadata$cellID <- rownames(mag_metadata)

mag_metadata <- subset(mag_metadata,
                   cellID %in% mag_relaxed_seurat@meta.data$cellID)


### Check which cells are in restricted dataset ###
common_cells <- which(mag_relaxed_seurat@meta.data$cellID %in% mag_xtoa$cell_id)

mag_relaxed_seurat <- subset(x = mag_relaxed_seurat,
                     subset = cellID %in% mag_xtoa$cell_id,
                     invert = FALSE)

### Add ratio info as metadata ###
## Add total_ratio ##
mag_relaxed_seurat <- AddMetaData(object = mag_relaxed_seurat,
                           metadata = mag_xtoa$total_ratio,
                           col.name = 'total_ratio')

## Add normalised_ratio ##
mag_relaxed_seurat <- AddMetaData(object = mag_relaxed_seurat,
                           metadata = mag_xtoa$normalised_ratio,
                           col.name = 'normalised_ratio')


# Add rox1 expression
mag_relaxed_seurat <- AddMetaData(object = mag_relaxed_seurat,
                           metadata = mag_xtoa$rox1_expression,
                           col.name = 'rox1_expression')

# Add rox2 expression
mag_relaxed_seurat <- AddMetaData(object = mag_relaxed_seurat,
                           metadata = mag_xtoa$rox2_expression,
                           col.name = 'rox2_expression')

mag_relaxed_seurat[['umap']] <- CreateDimReducObject(embeddings = as.matrix(mag_umap[,c(1,2)]),
                                                     key = 'umap_',
                                                     global = TRUE,
                                                     assay = 'SCT')
#### Remove spermatid cells in s_annotation ####
mag_relaxed_seurat <- subset(x = mag_relaxed_seurat,
                             subset = s_annotation != 'spermatid')

#### Save Seurat object with UMAP projections to file #### 
SaveH5Seurat(object = mag_relaxed_seurat,
             filename = paste0("r_male_reproductive_glands_cell_type_filt_decontx_stringent_filt_sctransformed_with_umap_and_xtoa.h5Seurat"))

#### Read in Seurat file ####
mag_relaxed_seurat <- LoadH5Seurat("r_male_reproductive_glands_cell_type_filt_decontx_stringent_filt_sctransformed_with_umap_and_xtoa.h5Seurat")

Idents(mag_relaxed_seurat) <- 's_annotation'

DimPlot(mag_relaxed_seurat,
        label = TRUE)

mag_subset <- subset(mag_relaxed_seurat,
                     subset = s_annotation %in% c('male accessory gland main cell',
                                                  'ejaculatory bulb',
                                                  'ejaculatory bulb epithelium',
                                                  'anterior ejaculatory duct',
                                                  'male accessory gland secondary cell'))

mag_subset@meta.data$s_annotation <- gsub(pattern = 'male accessory gland main cell',
                                          replacement = 'main cell',
                                          mag_subset@meta.data$s_annotation)


mag_subset@meta.data$s_annotation <- gsub(pattern = 'male accessory gland secondary cell',
                                          replacement = 'secondary cell',
                                          mag_subset@meta.data$s_annotation)

Idents(mag_subset) <- 's_annotation'

mag_subset <- ScaleData(mag_subset,
                        features = rownames(mag_subset))

mag_xtoa <- FeaturePlot(mag_subset,
            features = 'normalised_ratio',
            reduction = "umap",
            label = TRUE,
            repel = TRUE,
            label.size = 2.5,
            pt.size = 1.2,
            order = TRUE,
            min.cutoff = 0,
            max.cutoff = 1.1) +
  scale_colour_gradient(low = '#6b008a',
                        high = '#E69F00',
                        name = 'high DC') +
  ggtitle("mean X:A") +
  labs(x = "UMAP 1",
       y = "UMAP 2",
       tag = 'low DC') +
  annotate("text", x = 11, y = -8.5,
           label = "main cell",
           size = 2.6) +
  theme_minimal() +
  theme(text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5,
                                  face = 'bold'),
        plot.tag.position = c(0.91, 0.29),
        plot.tag = element_text(size = 9,
                                face = 'bold'),
        legend.title = element_text(size = 9,
                                    face = 'bold',
                                    vjust = 2))


x_a_data <- data.frame(cellID = mag_subset@meta.data$cellID,
                       normalised_ratio = mag_subset@meta.data$normalised_ratio,
                       s_annotation = mag_subset@meta.data$s_annotation)

x_a_data$s_annotation <- factor(x_a_data$s_annotation,
                                levels = c('main cell',
                                           'secondary cell',
                                           'anterior ejaculatory duct',
                                           'ejaculatory bulb',
                                           'ejaculatory bulb epithelium'))

x_a_normalised_plot <- ggplot(data = x_a_data,
                                aes(x = normalised_ratio,
                                    y = s_annotation,
                                    colour = s_annotation,
                                    fill = s_annotation)) + 
  ggridges::geom_density_ridges(jittered_points = TRUE,
                                position = position_points_jitter(width = 0.05, height = 0),
                                point_shape = '|', 
                                point_size = 3,
                                point_alpha = 1,
                                alpha = 0.3,
                                scale = 1) +
  labs(x = "X:A mean expression",
       y = "") +
  scale_colour_manual(values = c('#e69f00', '#9786e6', '#288e60',
                                 '#abd20c', '#905395')) +
  scale_fill_manual(values = c('#e69f00', '#9786e6', '#288e60',
                                 '#abd20c', '#905395')) +
  lims(x = c(0.3, 1.1)) +
  theme(legend.position = "none",
        text = element_text(size = 11),
        panel.grid.major = element_line(color = "#999999",
                                        linetype = 'dashed'),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())


x_to_a_combi <- ggarrange(x_a_normalised_plot, mag_xtoa,
                          labels = c('(A)', '(B)'),
                          nrow = 1)

ggsave(file = paste0("mag_xtoa_ridges_and_umap.pdf"),
       plot = x_to_a_combi,
       width = 210, height = 90, units = c("mm"), dpi = 300)


#### Find highly expressed genes ####
### Get mean expression per gene per cluster ###
average_expression <- as.data.frame(AverageExpression(object = mag_subset,
                                                      assays = 'SCT',
                                                      #slot = 'counts',
                                                      slot = 'data',
                                                      group.by = 's_annotation'))

average_expression$gene <- rownames(average_expression)


#### Produce dotplot with expression of genes ####
## Run if plotting DCC genes ##
genes <- c("lncRNA:roX1", "lncRNA:roX2", "msl-1", "msl-2", "msl-3", "mle", "mof",
           "dsx")

## Run if plotting other sex determination pathway genes ##
genes <- c("tra", "tra2", "Sxl", "fru")

## Run if plotting other DCC genes ##
genes <- c("Mtor", "Clamp", "JIL-1", "Top2")

mag_subset@meta.data <- mag_subset@meta.data %>%
  mutate(extra_annotation = case_when(s_annotation == 'main cell' &
                                        normalised_ratio >= 0.6857191 ~ 'high DC main cell',
                                      s_annotation == 'main cell' &
                                        normalised_ratio < 0.6857191 ~ 'low DC main cell',
                                      s_annotation == 'anterior ejaculatory duct' ~
                                        'anterior ejaculatory duct',
                                      s_annotation == 'ejaculatory bulb' ~
                                        'ejaculatory bulb',
                                      s_annotation == 'ejaculatory bulb epithelium' ~
                                        'ejaculatory bulb epithelium',
                                      s_annotation == 'secondary cell' ~
                                        'secondary cell'))

Idents(mag_subset) <- 'extra_annotation'

dcc_dotplot <- DotPlot(object = subset(mag_subset,
                                       subset = extra_annotation %in%
                                         c('high DC main cell',
                                           'low DC main cell',
                                           'anterior ejaculatory duct',
                                           'secondary cell')),
                       features = genes,
                       cols = c('#6b008a', '#E69F00'),
                       scale = FALSE) &
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        text = element_text(size = 11))

ggsave(file = 'mag_clusters_other_dcc_genes_dotplot_unscaled_data.pdf',
       plot = dcc_dotplot,
       bg = "white",
       width = 160, height = 80, units = c("mm"), dpi = 300)


#### Produce UMAP plot with expression of DCC genes ####
genes <- c("lncRNA:roX1", "lncRNA:roX2", "msl-1", "msl-2", "msl-3", "mle", "mof",
           "dsx")
for (gene in genes) {
  umap_plot <- FeaturePlot(mag_subset,
                           features = gene,
                           reduction = "umap",
                           slot = 'counts',
                           pt.size = 0.25,
                           raster = F) +
    labs(colour = paste0(gene, ' expression')) + 
    scale_colour_gradientn(colours = alpha(c("#999999", "#E69F00", "#C48112"), 0.8)) +
    ggtitle('Male reproductive glands') +
    theme(plot.title = element_text(hjust = 0.5, size = 15),
          legend.position = "right",
          legend.justification = "center",
          legend.title = element_text(face = "bold", size = 8))  
  umap_plot <- LabelClusters(umap_plot, id = 'ident',
                             size = 2, box = TRUE,
                             label.padding = 0.12,
                             label.size = 0,
                             repel = TRUE)
  #### Save plot to file
  gene <- gsub("lncRNA:", "", gene)
  outfilename <- paste0(gene, '_male_rep_glands_umap_plot.png')
  ggsave(file = outfilename,
         plot = umap_plot,
         width = 210, height = 120, units = c("mm"), dpi = 300)
  rm(umap_plot)
  gc()
}


#### Get counts for DCC genes ####
genes <- c("lncRNA:roX1", "lncRNA:roX2", "msl-1", "msl-2", "msl-3", "mle", "mof",
           "dsx", "Mtor", "JIL-1", "Top2", "Clamp", "Sxl")

gene_counts <- GetAssayData(object = mag_subset,
                            layer = 'counts',
                            assay = 'SCT')[genes,]

gene_counts <- as.matrix(gene_counts)

gene_counts <- t(gene_counts)

gene_counts <- as.data.frame(gene_counts)

gene_counts$cellID <- rownames(gene_counts) 

gene_counts <- left_join(gene_counts, x_a_data,
                         by = c("cellID" = "cellID"))

normratio_density <- density(gene_counts$normalised_ratio)
optimize(approxfun(normratio_density$x,
                   normratio_density$y),
         interval = c(0.6,0.8))$minimum

gene_counts <- gene_counts %>%
  mutate(status = case_when(normalised_ratio >= 0.6857191 ~ 'high DC',
                            normalised_ratio < 0.6857191 ~ 'low DC'))

gene_counts <- pivot_longer(gene_counts,
                            cols = 1:13,
                            names_to = 'gene',
                            values_to = 'gene_counts')

gene_counts <- gene_counts %>% filter(s_annotation == 'main cell')

#### Set plot theme ####
my_theme <- theme(legend.position = "bottom",
                  text = element_text(size = 11),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.spacing = unit(0, "mm"),
                  panel.background = element_blank(),
                  axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                  axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))


### Plot Sxl, roX1 and roX2 for high/low DC peaks ###
wilcox_results1 <- subset(gene_counts, gene %in% c("lncRNA:roX1",
                                                   "lncRNA:roX2",
                                                   "Sxl")) %>%
  group_by(gene) %>%
  summarise(pval = wilcox.test(gene_counts ~ status)$p.value)

wilcox_results1$pval <- as.numeric(wilcox_results1$pval)

wilcox_results1$pval_bonf <- p.adjust(wilcox_results1$pval,
                                       method = "bonferroni")

wilcox_results1 <- wilcox_results1 %>%
  mutate(label = case_when(pval_bonf < 0.001 ~ '***',
                           pval_bonf < 0.01 ~ '**',
                           pval_bonf < 0.05 ~ '*',
                           is.na(pval_bonf) ~ '',
                           TRUE ~ 'N.S.'))

wilcox_results1$status <- 'high DC'

dcc_boxplots <- ggplot(data = subset(gene_counts, gene %in% c("lncRNA:roX1",
                                                              "lncRNA:roX2",
                                                              "Sxl")),
                       aes(y = log(gene_counts + 1),
                           x = gene,
                           fill = status)) +
  geom_boxplot(notch = TRUE,
               outlier.alpha = 0.6,
               outlier.shape = 95,
               outlier.size = 5) +
  geom_text(data = wilcox_results1,
            mapping = aes(y = 4.8,
                          label = label),
            colour = 'black',
            size = 5) +
  labs(x = "Gene",
       y = "log(counts + 1)") +
  scale_fill_manual(name = 'Dosage compensation status',
                    values = c('#E69F00', '#6b008a')) +
  my_theme +
  theme(plot.margin = unit(c(1, 3.5, 1, 3.5), 'cm'),
        legend.title = element_text(size = 8,
                                    hjust = 1.5))

ggsave(file = "mag_boxplot_dsx_rox1_rox2_for_high_and_low_dc_peaks.pdf",
       plot = dcc_boxplots,
       width = 120, height = 90, units = c("mm"), dpi = 300)


## Combine with X:A in MAG clusters and X:A UMAP
blank <- ggplot() + theme_void()
x_to_a_combi <- ggarrange(ggarrange(x_a_normalised_plot,
                                    mag_xtoa,
                                    labels = c('(A)', '(B)'),
                                    nrow = 1),
                          ggarrange(blank,
                                    dcc_boxplots,
                                    blank,
                                    nrow = 1,
                                    labels = c('', '   (C)', ''),
                                    widths = c(0.12, 1, 0.12)),
                          heights = c(1, 0.8),
                          nrow = 2)

ggsave(file = "mag_xtoa_ridges_and_umap_plus_boxplot_dsx_rox1_rox2_for_high_and_low_dc_peaks.pdf",
       plot = x_to_a_combi,
       width = 210, height = 180, units = c("mm"), dpi = 300)


#### Perform DE between main cell clusters ####
mag_subsubset <- subset(x = mag_subset,
                        subset = s_annotation == 'main cell')

mag_subsubset@meta.data <- mag_subsubset@meta.data %>%
  mutate(status = ifelse((mag_subsubset@meta.data$normalised_ratio >= 0.6857191),
                         "high DC",  "low DC"))

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
status <- factor(colData(mag_sca)$status)
# set the reference level
status <- relevel(status,"high DC")
# define and fit the model
zlmCond <- zlm(formula = ~ngeneson + status, 
               sca = mag_sca)
# perform likelihood-ratio test for the condition that we are interested in    
summaryCond <- summary(zlmCond, doLRT = 'statuslow DC')
# get the table with log-fold changes and p-values
summaryDt <- summaryCond$datatable
result <- merge(summaryDt[contrast=='statuslow DC' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                summaryDt[contrast=='statuslow DC' & component=='logFC', .(primerid, coef)],
                by='primerid') # logFC coefficients
# MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
result[,coef:=result[,coef]/log(2)]
# do multiple testing correction
result[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
result = result[result$FDR < 0.01,, drop = F]

result <- stats::na.omit(as.data.frame(result))
write.table(result, "mast_result_mag_main_cell_clusters.txt",
            row.names = FALSE)

result <- read.table("mast_result_mag_main_cell_clusters.txt",
                     header = TRUE)

result$log10FDR <- -log10(result$FDR)

result <- result %>%
  mutate(status = case_when(abs(coef) >= 1 ~ 'DE',
                            abs(coef) < 1 ~ 'unbiased'))


### Compare with Majane et al results ###
result_majane <- read.table("../majane_2022_data/majane_dmel_mast_result_mag_main_cell_clusters.txt",
                            header = TRUE)

result_majane$log10FDR <- -log10(result_majane$FDR)

result_majane <- result_majane %>%
  mutate(status = case_when(abs(coef) >= 1 ~ 'DE',
                            abs(coef) < 1 ~ 'unbiased'))


result_de <- result %>%
  filter(status == 'DE') %>%
  left_join(result_majane, by = 'primerid')

result_de <- result_de %>% 
  arrange(desc(coef.x))

write.csv(x = result_de,
          file = "mast_results_our_de_genes_in_majane_dmel.csv",
          row.names = FALSE)

fdr_vs_fdr <- ggplot(data = result_de,
                  aes(x = log10FDR.x,
                      y = log10FDR.y,
                      colour = status.y)) +
  geom_point() +
  geom_smooth(method = 'lm', se = TRUE, level = 0.95,
              colour = '#4c4c4c') +
  scale_colour_manual(name = 'Status in Majane et al (2022) dataset',
                      values = c('#E69F00', '#353535'),
                      na.translate = FALSE) +
  labs(x = "-log10(adjusted p-value our study)",
       y = expression("-log10(adjusted p-value Majane et al (2022))")) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 12),
        legend.text = element_text(size = 10))

fc_vs_fc <- ggplot(data = result_de,
                     aes(x = abs(coef.x),
                         y = abs(coef.y),
                         colour = status.y)) +
  geom_point() +
  geom_smooth(method = 'lm', se = TRUE, level = 0.95,
              colour = '#4c4c4c') +
  scale_colour_manual(name = 'Status in Majane et al (2022) dataset',
                      values = c('#E69F00', '#353535'),
                      na.translate = FALSE) +
  labs(x = "log2(Fold Change our study)",
       y = expression("log2(Fold Change Majane et al (2022))"),
       colour = "") +     
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 12),
        legend.text = element_text(size = 10))


fdr_and_fc <- ggarrange(fdr_vs_fdr, fc_vs_fc,
                        labels = c('(A)', '(B)'),
                        nrow = 1,
                        common.legend = TRUE,
                        legend = 'bottom',
                        hjust = 0.02)

ggsave(file = paste0("adjust_pval_and_log2fc_ours_vs_majane_clean.pdf"),
       plot = fdr_and_fc,
       width = 210, height = 120, units = c("mm"), dpi = 300)


#### Plot number of X and autosomal expressed genes ####
### Load libraries ###
library(ape)
library(Seurat)
library(SeuratDisk)
library(tidyverse)

### Read in annotation file ###
gff_data <- read.gff(file = "./fly_cell_atlas_proj/data/dmel-genes-r6.31.gff")
gff_data$gene <- str_match(gff_data$attributes, 'Name=(.*?);Alias=')[1,2]
gff_data <- gff_data %>%
  mutate(gene = str_match(attributes, 'Name=(.*?);')[,2])

### Extract number of expressed genes per cluster###
## Read in Seurat file ##
mag_seurat <- LoadH5Seurat(paste0('./r_male_reproductive_glands_cell_type_filt_decontx_stringent_filt_sctransformed.h5Seurat'))

mag_subset <- subset(mag_seurat,
                     subset = s_annotation %in% c('male accessory gland main cell',
                                                  'ejaculatory bulb',
                                                  'ejaculatory bulb epithelium',
                                                  'anterior ejaculatory duct',
                                                  'male accessory gland secondary cell'))

rm(mag_seurat)
gc()

## Get counts data ##
counts_data <- GetAssayData(object = mag_subset,
                            layer = 'counts',
                            #slot = 'counts',
                            assay = 'SCT')

metadata_df <- mag_subset@meta.data[, c('cellID', 'sex', 's_annotation')]

rm(mag_subset)
gc()

counts_data <- as.matrix(counts_data)
gc()

counts_data <- as.data.frame(counts_data)
gc()

counts_data <- rownames_to_column(counts_data,
                                  var = 'gene')
gc()

counts_data <- pivot_longer(counts_data,
                            cols = 2:length(colnames(counts_data)),
                            names_to = 'cellID',
                            values_to = 'counts')
gc()

counts_data <- counts_data %>%
  filter(counts >= 1)
gc()

#### Match cell ID with metadata ####
counts_data <- left_join(counts_data, metadata_df)
rm(metadata_df)
gc()

#### Find genes in GFF data ####
tissue_genes <- as.data.frame(subset(gff_data, gff_data$gene 
                                     %in%
                                       counts_data$gene)[,c('seqid', 'gene')])
colnames(tissue_genes) <- c('chromosome', 'gene')

#### Set chromosome status ####
tissue_genes <- tissue_genes %>%
  mutate(chr_status = case_when((chromosome == "X") ~ "X",
                                (chromosome == "Y") ~ "Y",
                                (chromosome == "mitochondrion_genome") ~ "mt",
                                TRUE ~ "autosome"))

counts_data <- left_join(counts_data, tissue_genes)
rm(tissue_genes)
gc()

## Read in X:A ratio data ##
x_a_ratio_df_plus <- read.csv("x_a_normalised_expression_ratio_n_max_expression_data_relaxed_after_filtering_all_tissues.csv",
                              header = TRUE)

x_a_ratio_df_plus <- subset(x_a_ratio_df_plus,
                            tissue == "male reproductive glands")

names(counts_data)[names(counts_data) == 'cellID'] <- 'cell_id'

counts_data <- left_join(counts_data,
                         x_a_ratio_df_plus[, c('cell_id',
                                               'normalised_ratio')])

counts_data <- counts_data %>%
  mutate(status = case_when(normalised_ratio >= 0.6857191 ~ 'high DC',
                            normalised_ratio < 0.6857191 ~ 'low DC'))

counts_data <- counts_data %>%
  mutate(extra_status = case_when(s_annotation == 'anterior ejaculatory duct' ~
                                    'anterior ejaculatory duct',
                                  s_annotation == 'ejaculatory bulb' ~
                                    'ejaculatory bulb',
                                  s_annotation == 'ejaculatory bulb epithelium' ~
                                    'ejaculatory bulb epithelium',
                                  s_annotation == 'male accessory gland secondary cell' ~
                                    'secondary cell',
                                  s_annotation == 'male accessory gland main cell' &
                                    status == 'high DC' ~ 'main cell high DC',
                                  s_annotation == 'male accessory gland main cell' &
                                    status == 'low DC' ~ 'main cell low DC'))
## Run if ##
# plotting the percentage of genes expressed per cluster #
counts_data_summary <- counts_data %>%
  filter(counts >= 1) %>%
  group_by(cell_id, chr_status, extra_status, sex) %>%
  summarise(expressed_genes = n())
#

rm(counts_data)
gc()


## Set plot theme ##
my_theme <- theme(legend.position = "bottom",
                  text = element_text(size = 11),
                  strip.text = element_text(face = 'bold'),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.spacing = unit(0, "mm"),
                  panel.background = element_blank(),
                  axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                  axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))


number_of_genes <- gff_data %>%
  filter(seqid %in% c('2L', '2R', '3L', '3R', '4', 'X') & type == 'gene') %>%
  mutate(chr_status = case_when((seqid == "X") ~ "X",
                                (seqid == "Y") ~ "Y",
                                (seqid == "mitochondrion_genome") ~ "mt",
                                TRUE ~ "autosome")) %>%
  group_by(chr_status) %>%
  summarise(number_of_genes = n())


counts_data_summary <- left_join(counts_data_summary,
                                 number_of_genes)

counts_data_summary$percentage_of_genes <- counts_data_summary$expressed_genes/counts_data_summary$number_of_genes


### Create boxplots ###
## Load libraries ##
library(ggpubr)

x_and_autosomes <- ggplot(data = subset(counts_data_summary,
                                        chr_status %in% c('X', 'autosome')),
                          aes(x = chr_status,
                              y = percentage_of_genes*100)) +
  geom_violin(aes(group = chr_status,
                  colour = extra_status,
                  fill = extra_status),
              alpha = 0.6,
              linewidth = 0.6) +
  geom_boxplot(aes(group = chr_status,
                   colour = extra_status),
               linewidth = 0.6,
               width = 0.2,
               notch = TRUE,
               outlier.alpha = 0.6,
               outlier.shape = 95,
               outlier.size = 5,
               fill = 'white') +
  stat_compare_means(comparisons = list(c('autosome', 'X')),
                     method = 'wilcox.test',
                     label = 'p.signif',
                     label.y = rep(12.5, times = 6),
                     tip.length = 0.005) +
  scale_colour_manual(guide = 'none',
                      values = c('#288e60', '#abd20c', '#905395',
                                 '#e69f00', '#f0e442', '#9786e6')) +
  scale_fill_manual(guide = 'none',
                    values = c('#288e60', '#abd20c', '#905395',
                               '#e69f00', '#f0e442', '#9786e6')) +
  labs(x = '',
       y = 'Percentage of expressed genes (%)') +
  lims(y = c(NA, 13.5)) +
  facet_wrap(vars(extra_status)) +
  my_theme

ggsave(file = "percent_x_and_autosomal_expressed_genes_per_cluster_mag.pdf",
       plot = x_and_autosomes,
       width = 140, height = 160, units = c("mm"), dpi = 300)


#### Produce dotplot with expression of male germline markers ####
### Read in marker data ###
germline_markers <- read.csv("fca_male_germline_cell_markers.csv")

germline_markers$marker <- gsub(" ", "",
                                germline_markers$marker)

germline_markers_summary <- germline_markers %>%
  group_by(marker) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts))

top_genes <- germline_markers_summary$marker[1:16]
top_genes <- sort(c("tHMG1", "tHMG2", "D1", "Phs", "Zif",
                    "HmgD", "CG12104", "CG34367", "kmg", "CG11970",
                    "CG3995", "CG12942", "CG13617", "p53", "stwl",
                    "Doc1", "vis", "dmrt93B", "CG11294", "dac"))

mag_subset@meta.data <- mag_subset@meta.data %>%
  mutate(extra_annotation = case_when(s_annotation == 'main cell' &
                                        normalised_ratio >= 0.6857191 ~ 'high DC main cell',
                                      s_annotation == 'main cell' &
                                        normalised_ratio < 0.6857191 ~ 'low DC main cell',
                                      s_annotation == 'anterior ejaculatory duct' ~
                                        'anterior ejaculatory duct',
                                      s_annotation == 'ejaculatory bulb' ~
                                        'ejaculatory bulb',
                                      s_annotation == 'ejaculatory bulb epithelium' ~
                                        'ejaculatory bulb epithelium',
                                      s_annotation == 'secondary cell' ~
                                        'secondary cell'))

Idents(mag_subset) <- 'extra_annotation'

## Read in testis data ##
testis_relaxed_seurat <- LoadH5Seurat("r_testis_cell_type_filt_decontx_stringent_filt_sctransformed.h5Seurat")

testis_subset <- subset(testis_relaxed_seurat,
                        subset = s_broad_annotation %in% c('male germline cell'))

testis_subset@meta.data <- testis_subset@meta.data %>%
  mutate(extra_annotation = case_when(s_annotation %in% c('spermatocyte 0',
                                                          'spermatocyte 1',
                                                          'spermatocyte 2',
                                                          'spermatocyte 3',
                                                          'spermatocyte 4',
                                                          'spermatocyte 5',
                                                          'spermatocyte 6',
                                                          'spermatocyte 7a',
                                                          'spermatocyte',
                                                          'late primary spermatocyte') ~
                                        'spermatocyte',
                                      s_annotation %in% c('early-mid elongation-stage spermatid',
                                                          'early elongation stage spermatid',
                                                          'mid-late elongation-stage spermatid',
                                                          'spermatid') ~
                                        'spermatid',
                                      s_annotation %in% c('spermatogonium',
                                                          'spermatogonium-spermatocyte transition',
                                                          'mid-late proliferating spermatogonia') ~
                                        'spermatogonia'))

Idents(testis_subset) <- 'extra_annotation'

merged_mag_testis <- merge(x = mag_subset, y = testis_subset)

merged_mag_testis <- subset(merged_mag_testis,
                            subset = extra_annotation %in% c('ejaculatory bulb',
                                                             'ejaculatory bulb epithelium'),
                            invert = TRUE)

germline_marker_dotplot <- DotPlot(object = merged_mag_testis,
                                   assay = "SCT",
                                   features = top_genes,
                                   cols = c('#6b008a', '#E69F00'),
                                   scale = FALSE) &
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        text = element_text(size = 11),
        axis.text.x = element_text(size = 7.5,
                                   angle = 10,
                                   vjust = 0.7))

ggsave(file = 'mag_clusters_male_germline_markers_dotplot_unscaled_data.png',
       plot = germline_marker_dotplot,
       width = 290, height = 120, units = c("mm"), dpi = 300)


#### Plot dsx expression across tissues ####
#### Read in compiled data from csv ####
x_a_ratio_df_plus <- read.csv("x_a_normalised_expression_ratio_n_max_expression_data_relaxed_after_filtering_all_tissues.csv",
                              header = TRUE)

x_a_ratio_df_plus$tissue <- gsub("^body$", 'carcass',
                                 x_a_ratio_df_plus$tissue)

x_a_ratio_df_plus$tissue <- gsub("proboscis and maxillary palps", 'proboscis',
                                 x_a_ratio_df_plus$tissue)

x_a_ratio_df_plus$broad_annotation <- gsub("female germline cell",
                                           "germline cell",
                                           x_a_ratio_df_plus$broad_annotation)


### Further cell filtering ###
x_a_ratio_df_plus <- x_a_ratio_df_plus %>%
  filter(!(tissue == 'carcass' &
             broad_annotation == 'female reproductive system')) %>%
  filter(!(tissue %in% c('gonad',
                         'male reproductive glands') &
             broad_annotation == 'muscle cell')) %>%
  filter(!(tissue %in% c('gonad',
                         'male reproductive glands') &
             broad_annotation == 'hemocyte')) %>%
  filter(!(tissue %in% c('body wall', 'heart', 'leg', 'malpighian tubule',
                         'oenocyte', 'proboscis',
                         'gonad') &
             broad_annotation == 'tracheolar cell')) %>%
  filter(!(tissue %in% c('body wall',
                         'heart') &
             broad_annotation == 'salivary gland')) %>%
  filter(!(tissue == 'gonad' &
             broad_annotation == 'somatic precursor cell')) %>%
  filter(!(tissue %in% c('gut',
                         'heart') &
             broad_annotation == 'excretory system')) %>%
  filter(!(tissue == 'male reproductive glands' &
             broad_annotation == 'reproductive system')) %>%
  filter(!(tissue == 'male reproductive glands' &
             broad_annotation == 'germline'))

### Filter out male reproductive glands ###
x_a_ratio_df_plus <- x_a_ratio_df_plus %>%
  filter(tissue != 'male reproductive glands')

### Calculate dsx median for each tissue ###
median_per_tissue <- x_a_ratio_df_plus %>%
  group_by(tissue, sex) %>%
  summarise(median_dsx = median(dsx_expression))


### Set colour status for head & carcass tissues ###
x_a_ratio_df_plus <- x_a_ratio_df_plus %>%
  mutate(colour_status = case_when(tissue %in% c('head', 'antenna',
                                                 'proboscis') ~
                                     'head colour',
                                   .default = 'carcass colour'))

median_per_tissue <- median_per_tissue %>%
  mutate(colour_status = case_when(tissue %in% c('head', 'antenna',
                                                 'proboscis') ~
                                     'head colour',
                                   .default = 'carcass colour'))

### Set colour name vector
cols <- c("head colour" = '#288e60',
          "carcass colour" = '#905395')

### Set plot theme ###
my_theme <- theme(legend.position = "bottom",
                  text = element_text(size = 11),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.spacing = unit(0, "mm"),
                  panel.background = element_blank(),
                  axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                  axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))


x_a_ratio_df_plus$tissue <- factor(x_a_ratio_df_plus$tissue,
                                   levels = c('head', 'antenna', 'proboscis',
                                              'carcass', 'trachea', 'heart',
                                              'malpighian tubule', 'gut',
                                              'oenocyte', 'fat body', 'wing',
                                              'haltere', 'leg', 'body wall',
                                              'gonad'))

median_per_tissue$tissue <- factor(median_per_tissue$tissue,
                                   levels = c('head', 'antenna', 'proboscis',
                                              'carcass', 'trachea', 'heart',
                                              'malpighian tubule', 'gut',
                                              'oenocyte', 'fat body', 'wing',
                                              'haltere', 'leg', 'body wall',
                                              'gonad'))

dsx_expression_plot <- ggplot(data = x_a_ratio_df_plus,
                              aes(x = tissue)) +
  geom_boxplot(aes(y = dsx_expression,
                   fill = colour_status),
               notch = TRUE,
               outliers = FALSE,
               width = 0.65) +
  geom_vline(xintercept = 3.45,
             colour = '#288e60',
             linetype = "dashed",
             linewidth = 0.25) +
  geom_vline(xintercept = 3.55,
             colour = '#905395',
             linetype = "dashed",
             linewidth = 0.25) +
  annotate(geom = "text",
           x = 3.25, y = Inf,
           hjust = 1,
           label = "In the head",
           color = '#288e60',
           angle = 90) +
  annotate(geom = "text",
           x = 3.75, y = Inf,
           hjust = -0.001,
           label = "In the carcass",
           color = '#905395',
           angle = 270) +
  geom_label(data = median_per_tissue,
             aes(y = median_dsx,
                 label = median_dsx,
                 colour = colour_status,
                 fontface = 'bold'),
             alpha = 0.75,
             size = 2.5,
             vjust = -0.4,
             label.size = NA) +
  facet_wrap(vars(sex),
             nrow = 2,
             labeller = as_labeller(c(`female` = '(A) Females',
                                      `male` = '(B) Males'))) +
  labs(x = '',
       y = 'dsx expression counts') +
  scale_fill_manual(values = cols,
                    guide = 'none') +
  scale_colour_manual(values = cols,
                      guide = 'none') +
  my_theme +
  theme(axis.text.x = element_text(size = 7.5,
                                   angle = 10,
                                   vjust = 0.7),
        strip.text.x = element_text(hjust = 0,
                                    size = 11,
                                    face = 'bold'))

ggsave(file = 'dsx_expression_across_tissues.png',
       plot = dsx_expression_plot,
       width = 297, height = 160, units = c("mm"), dpi = 300)
