#######################################
#         DE analysis results:        #
#             head & body             #
#       DESEQ2 across the genome      #
#######################################

# Load libraries
library(tidyverse)

# Loop over all tissue samples
tissues <- c('Head', 'Body')

# Set function to read in data
read_n_add_tissue <- function(tissue, level) {
  filename <- paste0("deseq_results_from_raw_relaxed_", 
                     tissue,
                     "_per_",
                     level,
                     ".csv")
  tissue_data <- read.csv(filename, header = TRUE)
  tissue_data$tissue <- tolower(tissue)
  return(tissue_data)
}

# Read in tissue-level results
tissue_results <- tolower(tissues) %>%
  map_df(~ read_n_add_tissue(., level = 'tissue'))

tissue_results$log2FoldChange <- tissue_results$log2FoldChange*(-1) 

# Set significance status
tissue_results <- tissue_results %>%
  mutate(status = case_when(tissue_results$padj <= 0.05 &
                              tissue_results$log2FoldChange >= 1 ~ "female-biased",
                            tissue_results$padj <= 0.05 &
                              tissue_results$log2FoldChange <= -1 ~ "male-biased",
                            tissue_results$padj <= 0.05 &
                              tissue_results$log2FoldChange < 1 |
                              tissue_results$padj <= 0.05 &
                              tissue_results$log2FoldChange > -1 ~ "n.s. log2FC",
                            .default = "unbiased"))


# Read in cluster-level results
cluster_results <- tolower(tissues) %>%
  map_df(~ read_n_add_tissue(., level = 'cluster'))


cluster_results$log2FoldChange <- cluster_results$log2FoldChange*(-1) 

cluster_results <- cluster_results %>%
  mutate(status = case_when(cluster_results$padj <= 0.05 &
                              cluster_results$log2FoldChange >= 1 ~ "female-biased",
                            cluster_results$padj <= 0.05 &
                              cluster_results$log2FoldChange <= -1 ~ "male-biased",
                            cluster_results$padj <= 0.05 &
                              cluster_results$log2FoldChange < 1 |
                              cluster_results$padj <= 0.05 &
                              cluster_results$log2FoldChange > -1 ~ "n.s. log2FC",
                            .default = "unbiased"))


# Match each gene to its respective chromosome
## Read in gff file
library(ape)
gff_data <- read.gff(file = "./fly_cell_atlas_proj/data/dmel-genes-r6.31.gff")
gff_data$gene <- str_match(gff_data$attributes, 'Name=(.*?);Alias=')[1,2]
gff_data <- gff_data %>%
  mutate(gene = str_match(attributes, 'Name=(.*?);')[,2])
colnames(gff_data)[1] <- 'chromosome'

tissue_results <- left_join(tissue_results,
                            gff_data[, c('chromosome', 'gene')])
cluster_results <- left_join(cluster_results,
                             gff_data[, c('chromosome', 'gene')])


# Calculate observed/expected sex-biased genes per chromosome
## For tissue level results
tissue_chr <- tissue_results %>%
  filter(chromosome %in% c('2L', '2R', '3L', '3R', '4', 'X')) %>%
  group_by(tissue, status, chromosome) %>%
  summarise(sum = n())

tissue_chr_total <- tissue_results %>%
  filter(chromosome %in% c('2L', '2R', '3L', '3R', '4', 'X')) %>%
  group_by(tissue, chromosome) %>%
  summarise(total_chr = n()) %>%
  group_by(tissue) %>%
  mutate(total_genes = sum(total_chr)) %>%
  slice(rep(1:n(), each = length(unique(tissue_chr$status)))) %>%
  mutate(status = rep(unique(tissue_chr$status),
                      times = length(unique(tissue_chr$chromosome))))

tissue_chr <- full_join(tissue_chr, tissue_chr_total)

tissue_chr[is.na(tissue_chr)] <- 0

tissue_chr_sb <- tissue_results %>%
  filter(chromosome %in% c('2L', '2R', '3L', '3R', '4', 'X')) %>%
  group_by(tissue, status) %>%
  summarise(total_sb = n())

tissue_chr <- left_join(tissue_chr, tissue_chr_sb)

tissue_chr <- tissue_chr %>%
  mutate(total_percent = total_chr/total_genes,
         expected = total_sb*total_percent,
         obs_exp = sum/expected)

### Average obs/exp sex-biased genes for all autosomes
tissue_auto <- tissue_chr %>%
  filter(!(chromosome == 'X')) %>%
  group_by(tissue, status) %>%
  summarise(mean_ratio = mean(obs_exp))

### Obs/exp sex-biased genes for each chromosome
tissue_chr <- tissue_chr %>%
  filter(status %in% c('female-biased', 'male-biased'))

#### Save table as csv
write.csv(tissue_chr[,c('tissue', 'status', 'chromosome', 'obs_exp')],
          "observed_over_expected_sex_biased_genes_per_chr_per_tissue.csv",
          row.names = FALSE)

### Combine autosomal data
tissue_x_auto <- tissue_chr %>%
  mutate(chr_status = case_when(chromosome == 'X' ~ 'X',
                                chromosome != 'X' ~ 'autosome')) %>%
  group_by(tissue, chr_status, status) %>%
  summarise(sum = sum(sum),
            total_chr = sum(total_chr),
            total_sb = sum(total_sb),
            total_percent = sum(total_percent)) %>%
  mutate(expected = total_sb*total_percent,
         obs_exp = sum/expected,
         non_sb_sum = total_chr - sum)

#### Save table as csv
write.csv(tissue_x_auto[,c('tissue', 'status', 'chr_status', 'obs_exp')],
          "observed_over_expected_sex_biased_genes_autosomes_and_xchr_per_tissue.csv",
          row.names = FALSE)

### Fisher's Exact test
tissue_fet <- tissue_x_auto %>%
  group_by(tissue, status) %>%
  summarise(pvalue = fisher.test(x = matrix(c(sum, non_sb_sum),
                                            nrow = 2))$p.value,
            odds_ratio = fisher.test(x = matrix(c(sum, non_sb_sum),
                                                nrow = 2))$estimate)

tissue_fet <- tissue_fet %>%
  group_by(tissue) %>%
  mutate(padj = p.adjust(pvalue, method = 'BH'))

#### Save table as csv
write.csv(tissue_fet,
          "fet_sex_biased_genes_autosomes_and_xchr_per_tissue_plus_padj.csv",
          row.names = FALSE)


## For cluster-level results
cluster_chr <- cluster_results %>%
  filter(chromosome %in% c('2L', '2R', '3L', '3R', '4', 'X')) %>%
  group_by(tissue, cluster, status, chromosome) %>%
  summarise(sum = n())

cluster_chr_total <- cluster_results %>%
  filter(chromosome %in% c('2L', '2R', '3L', '3R', '4', 'X')) %>%
  group_by(tissue, cluster, chromosome) %>%
  summarise(total_chr = n()) %>%
  group_by(tissue, cluster) %>%
  mutate(total_genes = sum(total_chr))

cluster_chr <- left_join(cluster_chr, cluster_chr_total)

cluster_chr_sb <- cluster_results %>%
  filter(chromosome %in% c('2L', '2R', '3L', '3R', '4', 'X')) %>%
  group_by(tissue, cluster, status) %>%
  summarise(total_sb = n())

cluster_chr <- left_join(cluster_chr, cluster_chr_sb)

cluster_chr <- cluster_chr %>%
  mutate(total_percent = total_chr/total_genes,
         expected = total_sb*total_percent,
         obs_exp = sum/expected)

### Average obs/exp sex-biased genes for all autosomes
cluster_auto <- cluster_chr %>%
  filter(!(chromosome == 'X')) %>%
  group_by(tissue, cluster, status) %>%
  summarise(mean_ratio = mean(obs_exp))

### Obs/exp sex-biased genes for each chromosome
cluster_chr <- cluster_chr %>%
  filter(status %in% c('female-biased', 'male-biased'))

write.csv(cluster_chr[,c('tissue', 'cluster', 'status', 'chromosome', 'obs_exp')],
          "observed_over_expected_sex_biased_genes_per_chr_per_cluster.csv",
          row.names = FALSE)

### Combine autosomal data
cluster_x_auto <- cluster_chr %>%
  mutate(chr_status = case_when(chromosome == 'X' ~ 'X',
                                chromosome != 'X' ~ 'autosome')) %>%
  group_by(tissue, cluster, chr_status, status) %>%
  summarise(sum = sum(sum),
            total_chr = sum(total_chr),
            total_sb = sum(total_sb),
            total_percent = sum(total_percent)) %>%
  mutate(expected = total_sb*total_percent,
         obs_exp = sum/expected,
         non_sb_sum = total_chr - sum)

#### Save table as csv
write.csv(cluster_x_auto[,c('tissue', 'cluster', 'status', 'chr_status', 'obs_exp')],
          "observed_over_expected_sex_biased_genes_autosomes_and_xchr_per_cluster.csv",
          row.names = FALSE)

### Fisher's Exact test
cluster_fet <- cluster_x_auto %>%
  group_by(tissue, cluster, status) %>%
  summarise(pvalue = fisher.test(x = matrix(c(sum, non_sb_sum),
                                            nrow = 2))$p.value,
            odds_ratio = fisher.test(x = matrix(c(sum, non_sb_sum),
                                                nrow = 2))$estimate)

cluster_fet <- cluster_fet %>%
  group_by(tissue, cluster) %>%
  mutate(padj =  p.adjust(pvalue, method = 'BH'))

#### Save table as csv
write.csv(cluster_fet,
          "fet_sex_biased_genes_autosomes_and_xchr_per_cluster_plus_padj.csv",
          row.names = FALSE)


