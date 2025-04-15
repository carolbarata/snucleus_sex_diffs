# Plotting pseudo-bulk vs bulk RNAseq
# expression profiles for head
# and carcass

## Load libraries ##
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(ggpmisc)

## Set working directory ##
setwd("./fca_relaxed_loom/")

## Read in bulk RNAseq ##
### Set sample information ###
samples_info <- data.frame(sample = rep('dmel', times = 6),
                           replicate = c(1:3, 1:3),
                           tissue = c(rep('head', times = 3),
                                      rep('carc', times = 3)))

### Set function to read in data ###
read_n_add_sample <- function(sample, tissue, replicate) {
  filename1 <- paste0(sample, "_", tissue, "_", "f", replicate,
                      "_rsem_counts.out.genes.results")
  sample_data <- read.table(filename1, header = TRUE)
  filename2 <- paste0(sample, "_", tissue, "_", "m", replicate,
                      "_rsem_counts.out.genes.results")
  sample_data_m <- read.table(filename2, header = TRUE)
  sample_data <- sample_data %>%
    left_join(sample_data_m, join_by(gene_id))
  sample_data$sample <- sample
  sample_data$replicate <- replicate
  sample_data$tissue <- tissue
  return(sample_data)
}


### Read in results ###
sample_results <- pmap_df(.l = samples_info, .f = read_n_add_sample)

sample_results$tissue <- gsub('carc', 'carcass',
                              sample_results$tissue)

## Replace FB id with gene name ##
### Read in gff file ###
library(ape)
gff_data <- read.gff(file = "/media/cdecastr/Samsung USB1/fly_cell_atlas_proj/data/dmel-genes-r6.31.gff")
gff_data$gene <- str_match(gff_data$attributes, 'Name=(.*?);Alias=')[1,2]
gff_data$fb_id <- str_match(gff_data$attributes, 'ID=(.*?);Name=')[1,2]
gff_data <- gff_data %>%
  mutate(gene = str_match(attributes, 'Name=(.*?);')[,2],
         fb_id = str_match(attributes, 'ID=(.*?);')[,2])
colnames(gff_data)[1] <- 'chromosome'
gff_data <- gff_data %>%
  mutate(gene_length = end - start)

colnames(sample_results)[1] <- 'fb_id'

sample_results <- left_join(sample_results,
                            gff_data[, c('chromosome', 'gene', 'fb_id')])


## Read in single nucleus RNAseq ##
tissues <- c('Head', 'Body')

for (tissue in tissues) {
  ### Read in Seurat file ###
  tissue_seurat <- LoadH5Seurat(paste0('./r_', tolower(tissue),
                                       "_cell_type_filt_decontx_stringent_filt_sctransformed.h5Seurat"))
  
  
  ### Further cell type filtering ###
  to_remove <- c('artefact',
                 #'unannotated',
                 'female reproductive system',
                 'male reproductive system')
  
  tissue_seurat_final <- subset(tissue_seurat,
                                subset = s_broad_annotation %in% to_remove,
                                invert = TRUE)

  rm(tissue_seurat)
  rm(to_remove)
  gc()
  
  ### Get metadata ###
  metadata_df <- tissue_seurat_final@meta.data[, c('cellID',
                                                   'sex',
                                                   'sample_id')]
  
  metadata_df$cellID <- colnames(tissue_seurat_final@assays$decontXcounts$counts)
  
  replicate_number_df <- metadata_df %>%
    select(sex, sample_id) %>%
    distinct(.keep_all = TRUE) %>%
    group_by(sex) %>%
    mutate(replicate = 1:length(unique(sample_id)))
  
  metadata_df <- left_join(metadata_df,
                           replicate_number_df)
  
  rm(replicate_number_df)
  gc()

  ### Split gene list into sublists ###
  gene_lists <- split(rownames(tissue_seurat_final@assays$decontXcounts),
                               ceiling(seq_along(rownames(tissue_seurat_final@assays$decontXcounts))/50))

  #### Loop over intervals and output to file ####
  for (gene_list in gene_lists) {
    #### Unlist gene_list ###
    gene_list <- unlist(gene_list)
    
    ### Get summary data per gene per sex ###
    #### Get counts data ####
    counts_data <- GetAssayData(object = tissue_seurat_final,
                                layer = 'counts',
                                assay = 'decontXcounts')[gene_list,]

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
    
    #### Match cell ID with metadata ####
    counts_data <- left_join(x = counts_data,
                             y = metadata_df,
                             by = 'cellID',
                             relationship = "many-to-one")
    
    gc()
    
    ### Calculate summary stats for females and males ###
    #### Sum counts across cells within a replicate for each gene ####
    counts_data_summary <- counts_data %>%
      group_by(gene, sex, sample_id, replicate) %>%
      summarise(sum = sum(counts, na.rm = TRUE))
    
    rm(counts_data)
    gc()
    
    ### Split female and male data across cols ###
    counts_data_summary_mod <- counts_data_summary %>%
      filter(sex == 'female')
    
    counts_data_summary_m <- counts_data_summary %>%
      filter(sex == 'male')
    
    counts_data_summary_mod <- left_join(counts_data_summary_mod,
                                         counts_data_summary_m,
                                         by = c('gene', 'replicate'),
                                         relationship = "one-to-one")
    
    rm(counts_data_summary)
    rm(counts_data_summary_m)
    gc()
    
    if (tissue == 'Body') {
      tissue <- 'Carcass'
    }
    counts_data_summary_mod$tissue <- tolower(tissue)
    
    if (file.exists("decontX_expression_data_per_gene_sex_rep_head_n_carc.csv")){
      write.table(counts_data_summary_mod,
                  "decontX_expression_data_per_gene_sex_rep_head_n_carc.csv",
                  append = TRUE,
                  sep = ",",
                  col.names = FALSE,
                  row.names = FALSE)
    } else {
      write.csv(counts_data_summary_mod,
                "decontX_expression_data_per_gene_sex_rep_head_n_carc.csv",
                row.names = FALSE)
    }
    
    rm(counts_data_summary_mod)
    gc()
  }
  
  rm(tissue_seurat_final)
  rm(gene_lists)
  rm(metadata_df)
  gc()
}


### Read in single nucleus data from file ###
pseudo_bulk_df <- read.csv(file = "expression_data_per_gene_sex_rep_head_n_carc.csv")

## Calculate average expression per replicate ##
### For bulk RNAseq
sample_results <- sample_results %>%
  group_by(gene, chromosome, tissue) %>%
  summarise(bulk_mean_sum_f = mean(TPM.x, na.rm = TRUE),
            bulk_mean_sum_m = mean(TPM.y, na.rm = TRUE))

### For pseudo-bulk RNAseq
#### Calculate TPM ####
pseudo_bulk_df <- pseudo_bulk_df %>%
  left_join(y = gff_data %>% select(gene, type, gene_length)) %>%
  mutate(gene_length_kb = gene_length/1000,
         rpk_f = sum.x/gene_length_kb,
         rpk_m = sum.y/gene_length_kb)

pseudo_bulk_df_total <- pseudo_bulk_df %>%
  group_by(tissue, replicate) %>%
  summarise(total_rpk_f = sum(rpk_f, na.rm = TRUE)/1e6,
            total_rpk_m = sum(rpk_m, na.rm = TRUE)/1e6)

pseudo_bulk_df <- pseudo_bulk_df %>%
  left_join(y = pseudo_bulk_df_total,
            relationship = 'many-to-one')

pseudo_bulk_df <- pseudo_bulk_df %>%
  mutate(tpm_f = rpk_f/total_rpk_f,
         tpm_m = rpk_m/total_rpk_m)

#### Calculate average ####
pseudo_bulk_df <- pseudo_bulk_df %>%
  group_by(gene, tissue) %>%
  summarise(pseudo_mean_sum_f = mean(tpm_f, na.rm = TRUE),
            pseudo_mean_sum_m = mean(tpm_m, na.rm = TRUE))


## Join with bulk data ##
merged_df <- left_join(x = sample_results,
                       y = pseudo_bulk_df,
                       by = c('gene', 'tissue'))

merged_df <- merged_df %>%
  filter(chromosome %in% c('3R', '3L', '2R', '2L', '4', 'X', 'Y'))

## Plot bulk vs pseudo-bulk expression ##
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


### Set colour name vector
cols <- c("female-biased" = "#287C8E",
          "male-biased" = "#E69F00",
          "head" = '#288e60',
          "carcass" = '#905395')

## Generate plot ##
merged_df_final <- merged_df %>%
  filter(chromosome %in% c('2L', '2R', '3L', '3R', '4', 'X'))

bulk_vs_pseudo_m <- ggplot(data = merged_df_final,
                           aes(x = log(bulk_mean_sum_m),
                               y = log(pseudo_mean_sum_m),
                               colour = tissue)) +
  geom_point(alpha = 0.4) +
  stat_poly_line(colour = "#999999",
                 na.rm = TRUE) +
  stat_poly_eq(na.rm = TRUE,
               use_label(c("eq", "R2")),
               colour = "#4C4C4C") +
  facet_wrap(vars(tissue)) +
  labs(x = 'log(TPM) for bulk RNA-seq',
       y = 'log(TPM) for pseudo-bulk RNA-seq') +
  scale_colour_manual(guide = 'none',
                      values = cols) +
  my_theme +
  theme(panel.spacing.x = grid::unit(0.26, "lines"),
        strip.text = element_text(face = 'bold',
                                  size = 11))

ggsave(file = paste0("bulk_tpm_vs_sctransformed_tpm_pseudobulk_head_n_carcass_males.pdf"),
       plot = bulk_vs_pseudo_m,
       bg = 'white',
       width = 160, height = 100, units = c("mm"), dpi = 300)

bulk_vs_pseudo_f <- ggplot(data = merged_df_final,
                           aes(x = log(bulk_mean_sum_f),
                               y = log(pseudo_mean_sum_f),
                               colour = tissue)) +
  geom_point(alpha = 0.4) +
  stat_poly_line(colour = "#999999",
                 na.rm = TRUE) +
  stat_poly_eq(na.rm = TRUE,
               use_label(c("eq", "R2")),
               colour = "#4C4C4C") +
  facet_wrap(vars(tissue)) +
  labs(x = 'log(TPM) for bulk RNA-seq',
       y = 'log(TPM) for pseudo-bulk RNA-seq') +
  scale_colour_manual(guide = 'none',
                      values = cols) +
  my_theme +
  theme(panel.spacing.x = grid::unit(0.26, "lines"),
        strip.text = element_text(face = 'bold',
                                  size = 11))

ggsave(file = paste0("bulk_tpm_vs_sctransformed_tpm_pseudobulk_head_n_carcass_females.pdf"),
       plot = bulk_vs_pseudo_f,
       bg = 'white',
       width = 160, height = 100, units = c("mm"), dpi = 300)
