#######################################
#         DE analysis results:        #
#             head & body             #
#             DESEQ2 plots            #
#######################################

#### Load libraries ####
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(lemon)

#### Loop over all tissue samples ####
tissues <- c('Head', 'Body')

#### Set function to read in data ####
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

#### Read in tissue-level results ####
tissue_results <- tolower(tissues) %>%
  map_df(~ read_n_add_tissue(., level = 'tissue'))

tissue_results$log2FoldChange <- tissue_results$log2FoldChange*(-1)

#### Set significance status ####
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

#### Replace 'body' with 'carcass' ####
tissue_results$tissue <- gsub("^body$", 'carcass', tissue_results$tissue)


#### Read in cluster-level results ####
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

#### Replace 'body' with 'carcass' ####
cluster_results$tissue <- gsub("^body$", 'carcass', cluster_results$tissue)


#### For stacked barplots ####
for_barplot_tissue <- tissue_results %>%
  group_by(tissue, status) %>%
  summarise(counts = n())

for_barplot_tissue$tissue <- factor(for_barplot_tissue$tissue,
                                    levels = c('head', 'carcass'))

tissue_barplot <- ggplot(for_barplot_tissue,
                         aes(fill = status,
                             y = counts,
                             x = tissue)) +
  geom_bar(position = "fill",
           stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(name = 'Status',
                    values = c("#287C8E", "#E69F00", "#999999", "#4C4C4C")) +
  labs(x = "Tissue",
       y = "Percentage of genes") +
  guides(fill = guide_legend(nrow = 2)) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 12),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank())

for_barplot_cluster <- cluster_results %>%
  group_by(tissue, cluster, status) %>%
  summarise(counts = n())

for_barplot_cluster$tissue <- factor(for_barplot_cluster$tissue,
                                     levels = c('head', 'body'))

for_barplot_cluster$cluster <- factor(for_barplot_cluster$cluster,
                                      levels = c("neuron",
                                                 "sensory neuron",
                                                 "glial cell",
                                                 "tracheolar cell",
                                                 "hemocyte",
                                                 "fat cell",
                                                 "epithelial cell",
                                                 "oenocyte",
                                                 "muscle cell",
                                                 "unannotated"))

cluster_barplot <- ggplot(for_barplot_cluster,
                         aes(fill = status,
                             y = counts,
                             x = cluster)) +
  geom_bar(position = "fill",
           stat = "identity") +
  scale_y_continuous(labels = scales::percent,
                     position = 'right') +
  scale_fill_manual(name = 'Status',
                    values = c("#287C8E", "#E69F00", "#999999", "#4C4C4C")) +
  labs(x = "Cell type",
       y = "Percentage of genes") +
  facet_grid(vars(tissue),
             switch = "y",
             scales = 'free',
             drop = TRUE) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 10),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        strip.text.y.left = element_text(angle = 0,
                                         hjust = 1),
        strip.background = element_blank(),
        plot.margin = unit(c(0,1.7,0,1.7), "cm"))



#### Match each gene to its respective chromosome ####
### Read in gff file ###
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


tissue_results <- tissue_results %>%
  filter(chromosome != 'mitochondrion_genome' &
           chromosome != 'Y')

cluster_results <- cluster_results %>%
  filter(chromosome != 'mitochondrion_genome' &
           chromosome != 'Y')

#### For % of expressed genes that are female- and male-biased ####
### Combine tissue and cluster results ###
tissue_results$cluster <- tissue_results$tissue
merged_df <- rbind(cluster_results,
                   tissue_results)

### Count total genes in head & carcass ###
total_genes <- merged_df %>%
  group_by(tissue) %>%
  distinct(gene) %>%
  summarise(total_expressed = n())

### Get genes that are both FB and MB ###
merged_biased <- merged_df %>%
  select(gene, tissue, status) %>%
  filter(status %in% c('female-biased', 'male-biased')) %>%
  group_by(tissue, gene) %>%
  summarise(number_status = n_distinct(status)) %>%
  filter(number_status > 1)

merged_biased %>%
  group_by(tissue) %>%
  summarise(number_of_genes = n())
  
### Split genes by tissue-only, cluster-only and shared ###
merged_df <- merged_df %>%
  filter(status %in% c('female-biased', 'male-biased')) %>%
  group_by(tissue, status, gene) %>%
  summarise(shared_status = case_when(n() >= 2 & 
                                        'head' %in% cluster ~ 'shared',
                                      n() >= 2 & 
                                        'carcass' %in% cluster ~ 'shared',
                                      n() == 1 &
                                        'head' %in% cluster ~ 'tissue-level only',
                                      n() == 1 &
                                        'carcass' %in% cluster ~ 'tissue-level only',
                                      .default = 'cluster-level only'))


for_barplot <-  left_join(x = merged_df,
                          y = total_genes,
                          by = 'tissue',
                          relationship = 'many-to-one') %>%
  group_by(tissue, status, shared_status, total_expressed) %>%
  summarise(counts = n()) %>%
  group_by(tissue) %>%
  mutate(percent_expressed = counts/total_expressed)

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

carc_data <- for_barplot %>%
  filter(tissue == 'carcass' & 
           status %in% c('female-biased', 'male-biased'))

carc_data$shared_status <- factor(x = carc_data$shared_status,
                                  levels = c('tissue-level only',
                                             'shared',
                                             'cluster-level only'))

carc_barplot <- ggplot(carc_data,
                         aes(fill = status,
                             y = percent_expressed,
                             x = shared_status)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = counts),
            position = position_stack(vjust = 0.5),
            colour = 'white',
            size = 2.5,
            fontface = "bold") +
  scale_y_continuous(labels = scales::label_percent(scale = 100,
                                                    suffix = "%")) +
  scale_fill_manual(name = 'Status',
                    values = c("female-biased" = "#287C8E",
                               "male-biased" = "#E69F00")) +
  labs(x = "Tissue",
       y = "Percentage of expressed genes") +
  guides(fill = guide_legend(nrow = 1)) +
  my_theme +
  theme(plot.margin = unit(c(1,0.5,0,0.5), 'cm'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 0.5, vjust = 0.65, 
                                   size = 8),
        legend.margin = margin(0,0,0,0))

head_data <- for_barplot %>%
  filter(tissue == 'head' & 
           status %in% c('female-biased', 'male-biased'))

head_data$shared_status <- factor(x = head_data$shared_status,
                                  levels = c('tissue-level only',
                                             'shared',
                                             'cluster-level only'))

head_barplot <- ggplot(head_data,
                       aes(fill = status,
                           y = percent_expressed,
                           x = shared_status)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = counts),
            position = position_stack(vjust = 0.5),
            colour = 'white',
            size = 2.5,
            fontface = "bold") +
  scale_y_continuous(labels = scales::label_percent(scale = 100,
                                                    suffix = "%")) +
  scale_fill_manual(name = 'Status',
                    values = c("female-biased" = "#287C8E",
                               "male-biased" = "#E69F00")) +
  labs(x = "Tissue",
       y = "Percentage of expressed genes") +
  guides(fill = guide_legend(nrow = 1)) +
  my_theme +
  theme(plot.margin = unit(c(1,0.5,0,0.5), 'cm'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 0.5, vjust = 0.65, 
                                   size = 8),
        legend.margin = margin(0,0,0,0))


#### For UpSet plots ####
### Load libraries ###
library(ComplexUpset)

### Set ggplot theme ###
my_theme <- theme(legend.position = "bottom",
                  text = element_text(size = 10),
                  strip.text = element_text(face = 'bold'),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.spacing = unit(0, "mm"),
                  panel.background = element_blank(),
                  axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                  axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))

### Combine tissue and cluster results ###
tissue_results$cluster <- tissue_results$tissue
merged_df <- rbind(cluster_results,
                   tissue_results)

### UpSet plot for head ###
## Subset sex biased genes in the head ##
merged_head <- merged_df %>%
   filter(status %in% c('female-biased',
                        'male-biased') &
            tissue == 'head' &
            cluster != 'unannotated') %>%
   mutate(val = 1)
 
cluster_order <- rev(c("head",
                       "neuron",
                       "sensory neuron",
                       "glial cell",
                       "hemocyte",
                       "fat cell",
                       "epithelial cell",
                       "muscle cell"))

merged_head <- merged_head %>%
  arrange(match(cluster, cluster_order))

merged_head$cluster <- factor(merged_head$cluster,
                              levels = c("head",
                                             "neuron",
                                             "sensory neuron",
                                             "glial cell",
                                             "hemocyte",
                                             "fat cell",
                                             "epithelial cell",
                                             "muscle cell"))


merged_head <- merged_head %>%
  pivot_wider(id_cols = gene,
              names_from = cluster,
              values_from = val,
              values_fill = list(n = 0)) %>%
  mutate(across(-gene, ~replace_na(.x, 0))) %>%
  mutate(across(-gene, ~ifelse(.x == 1, TRUE,FALSE)))

clusters <- colnames(merged_head)[2:length(colnames(merged_head))]


## UpSet plot ##
cluster_head_upset <- ComplexUpset::upset(merged_head, clusters,
                                          name = '',
                                          base_annotations = list('Number of sex-biased genes' = 
                                                                    intersection_size(text = list(size = 2.5,
                                                                                                  fontface = 'bold'))),
                                          set_sizes = FALSE,
                                          sort_sets = FALSE,
                                          min_degree = 2,
                                          #queries = list(upset_query(intersect = 'head',
                                          #                           color = '#8a377b',
                                          #                           fill = '#8a377b',
                                          #                           only_components = c('intersections_matrix',
                                          #                                               'Intersection size'))),
                                          themes = upset_default_themes(panel.grid.major = element_blank(),
                                                                        panel.grid.minor = element_blank()))



ggsave(file = paste0("upset_plot_head_deseq2_results_per_cluster.pdf"),
       plot = cluster_head_upset,
       width = 140, height = 100, units = c("mm"), dpi = 300,
       bg = 'white')


### UpSet plot for body ###
## Subset sex biased genes in the body ##
merged_body <- merged_df %>%
  filter(status %in% c('female-biased',
                       'male-biased') &
           tissue == 'carcass' &
           cluster != 'unannotated') %>%
  mutate(val = 1)

cluster_order <- rev(c("carcass",
                       "neuron",
                       "sensory neuron",
                       "glial cell",
                       "tracheolar cell",
                       "hemocyte",
                       "fat cell",
                       "epithelial cell",
                       "oenocyte",
                       "muscle cell"))

merged_body <- merged_body %>%
  arrange(match(cluster, cluster_order))

merged_body$cluster <- factor(merged_body$cluster,
                              levels = rev(c("carcass",
                                             "neuron",
                                             "sensory neuron",
                                             "glial cell",
                                             "tracheolar cell",
                                             "hemocyte",
                                             "fat cell",
                                             "epithelial cell",
                                             "oenocyte",
                                             "muscle cell")))


merged_body <- merged_body %>%
  pivot_wider(id_cols = gene,
              names_from = cluster,
              values_from = val,
              values_fill = list(n = 0)) %>%
  mutate(across(-gene, ~replace_na(.x, 0))) %>%
  mutate(across(-gene, ~ifelse(.x == 1, TRUE,FALSE)))

clusters <- colnames(merged_body)[2:length(colnames(merged_body))]


## UpSet plot ##
cluster_body_upset <- ComplexUpset::upset(merged_body, clusters,
                                          name = '',
                                          min_size = 8,
                                          base_annotations = list('Number of sex-biased genes' = 
                                                                    intersection_size(text = aes(size = 2,
                                                                                                 fontface = 'bold'))),
                                          set_sizes = FALSE,
                                          sort_sets = FALSE,
                                          min_degree = 2,
                                          # queries = list(upset_query(intersect = 'body',
                                          #                            color = '#8a377b',
                                          #                            fill = '#8a377b',
                                          #                            only_components = c('intersections_matrix',
                                          #                                                'Intersection size'))),
                                          themes = upset_default_themes(panel.grid.major = element_blank(),
                                                                        panel.grid.minor = element_blank()))


cluster_body_upset <- ggarrange(cluster_body_upset,
                                labels = 'Carcass', hjust = -10, vjust = 3)


combined_plots <- ggarrange(head_barplot, carc_barplot,
                            labels = c('    (A) Head', '    (B) Carcass'),
                            nrow = 1,
                            common.legend = TRUE,
                            legend = "bottom")

combined_plots <- ggarrange(combined_plots, cluster_body_upset,
                            labels = c('', '    (C)'),
                            nrow = 2,
                            heights = c(0.8, 1))

ggsave(file = paste0("upset_plot_body_n_barplots_deseq2_results_per_cluster.png"),
       plot = combined_plots,
       width =  150, height = 190, units = c("mm"), dpi = 300,
       bg = 'white')



### Plot X-linked to autosomal sex-biased genes across clusters ###
#### Read in observed/expected genes data per cluster ###
for_barplot_cluster <- read.csv(file = "observed_over_expected_sex_biased_genes_autosomes_and_xchr_per_cluster.csv")

for_barplot_cluster$tissue <- gsub("body", "carcass",
                                   for_barplot_cluster$tissue)

#### Read in FET results per cluster data ###
stats_data <- read.csv(file = "fet_sex_biased_genes_autosomes_and_xchr_per_cluster.csv")

stats_data$tissue <- gsub("body", "carcass",
                          stats_data$tissue)

for_barplot_cluster <- left_join(for_barplot_cluster, stats_data,
                                 by = c('tissue', 'cluster', 'status'),
                                 relationship = 'many-to-one')
for_barplot_cluster <- for_barplot_cluster %>%
  pivot_wider(id_cols = c('tissue', 'cluster', 'status', 'pvalue', 'odds_ratio'),
              names_from = chr_status,
              values_from = obs_exp)

for_barplot_cluster <- for_barplot_cluster %>%
  mutate(signif_label = case_when(pvalue < 0.001 ~ '***',
                                  pvalue < 0.01 ~ '**',
                                  pvalue < 0.05 ~ '*',
                                  TRUE ~ ''))

for_barplot_cluster <- for_barplot_cluster %>%
  mutate(tissue_mod = case_when(tissue == 'head' ~ '       (A) Head',
                                tissue == 'carcass' ~ '       (B) Carcass'))

for_barplot_cluster$cluster <- factor(for_barplot_cluster$cluster,
                                      levels = c("neuron",
                                                 "sensory neuron",
                                                 "glial cell",
                                                 "tracheolar cell",
                                                 "hemocyte",
                                                 "fat cell",
                                                 "epithelial cell",
                                                 "oenocyte",
                                                 "muscle cell",
                                                 "unannotated"))


cluster_barplot_body <- ggplot(data = for_barplot_cluster,
                          aes(x = status,
                              y = odds_ratio,
                              fill = status)) +
  geom_col()  +
  geom_text(aes(label = signif_label),
            vjust = 0,
            colour = "#4C4C4C",
            #position = position_dodge(width = 1),
            size = 3.2)+
  geom_hline(yintercept = 1,
             linetype = 'dashed',
             colour = "#999999") +
  scale_fill_manual(name = 'Status',
                    values = c("#287C8E", "#E69F00", "#999999", "#4C4C4C")) +
  labs(x = "Cell type cluster",
       y = "Odds ratio = X-linked obs/exp:autosomal obs/exp") +
  facet_grid(rows = vars(tissue_mod),
             cols = vars(cluster),
             drop = FALSE,
             scales = "free_y",
             switch = 'both') +
  theme_cowplot(font_size = 11) +
  theme(strip.background = element_rect(colour = NA, fill = "white"), 
        panel.border = element_rect(colour = NA),
        strip.placement = "outside",
        panel.spacing.x = grid::unit(0.26, "lines"),
        panel.spacing.y = grid::unit(2, "lines"), 
        legend.position = "bottom",
        legend.justification = "center",
        plot.margin = unit(c(1,1.7,0.25,0.25), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(hjust = 0.5, size = 8),
        strip.clip = "off",
        strip.text.y.left = element_text(angle = 0, vjust = 1),
        strip.text.y = element_text(margin = margin(t = -15, r = -20),
                                    face = 'bold',
                                    size = 11))

ggsave(file = paste0("barplot_aut_vs_x_head_n_body_deseq2_results_per_cluster.pdf"),
       plot = cluster_barplot_body,
       bg = 'white',
       width = 230, height = 190, units = c("mm"), dpi = 300)


#### Plot log2 F:M number of cells per cluster vs no. of DE genes ####
### Read in file with data on number of cells ###
number_of_cells <- read.csv(file = "./number_of_cells_in_head_and_body_and_clusters",
                            header = TRUE)

number_of_cells <- number_of_cells %>%
  mutate(log2_ratio = log2(ratio))

colnames(number_of_cells)[1] <- 'cluster'

number_of_cells$tissue <- gsub(pattern = 'body',
                               replacement = 'carcass',
                               x = number_of_cells$tissue)

### Summarise tissue and cluster results ###
tissue_results_summary <- tissue_results %>%
  group_by(tissue, status) %>%
  summarise(counts = n()) %>%
  filter(status %in% c('female-biased', 'male-biased'))

cluster_results_summary <- cluster_results %>%
  group_by(tissue, cluster, status) %>%
  summarise(counts = n()) %>%
  filter(status %in% c('female-biased', 'male-biased'))

### Combine cluster and tissue summaries ###
tissue_results_summary <- tissue_results_summary %>%
  mutate(cluster = case_when(tissue %in% c('head', 'carcass') ~
                               'tissue (total)'))

results_summary <- rbind(tissue_results_summary,
                         cluster_results_summary)

results_summary <- full_join(results_summary,
                             number_of_cells)


results_summary <- results_summary[complete.cases(results_summary),]

results_summary <- results_summary %>%
  pivot_wider(id_cols = c('tissue', 'cluster', 'female', 'male', 'ratio',
                           'log2_ratio'),
               names_from = status,
               values_from = counts)

results_summary <- results_summary %>%
  mutate(de_ratio = (`female-biased`/`male-biased`))

results_summary <- results_summary %>%
  mutate(total_de = (`female-biased` + `male-biased`))

### Generate plot ###
## Set plot theme ##
my_theme <- theme(legend.position = "bottom",
                  text = element_text(size = 11),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.spacing = unit(0, "mm"),
                  panel.background = element_blank(),
                  axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                  axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))

## Load required library ##
library(ggrepel)

## Create dataframe for shaded areas ##
rects <- data.frame(xstart = c(-Inf, 0),
                    xend = c(0, Inf),
                    ystart = c(-Inf, 1),
                    yend = c(1, Inf),
                    col = c('more male cells and male-biased genes',
                            'more female cells and female-biased genes'))

## Produce plot ##
de_genes_vs_cell_abundance <- ggplot(data = results_summary,
                                     aes(x = log2_ratio,
                                         y = de_ratio,
                                         shape = str_wrap(tissue),
                                         label = cluster)) +
  geom_rect(data = rects,
            inherit.aes = FALSE,
            aes(xmin = xstart,
                xmax = xend,
                ymin = ystart,
                ymax = yend, 
                fill = col),
            alpha = 0.3) +
  geom_point(aes(size = total_de),
                 na.rm = TRUE) +
  geom_text_repel(show.legend = FALSE,
                  size = 3,
                  box.padding = unit(0.42, "lines"),
                  colour = 'black') +
  geom_vline(xintercept = 0,
             linetype = 'dashed',
             colour = "#999999") +
  geom_hline(yintercept = 1,
             linetype = 'dashed',
             colour = "#999999") +
  annotate(geom = 'text', label = 'More male cells and\nmale-biased genes',
           x = -0.5, y = 0.62, size = 4, hjust = 0,
           colour = "#E69F00",
           fontface = "bold") +
  annotate(geom = 'text', label = 'More female cells and\nfemale-biased genes',
           x = 1.41, y = 2.3, size = 4, hjust = 1,
           colour = "#287C8E",
           fontface = "bold") +
  labs(y = "#female-biased/#male-biased genes",
       x = "log2(female cells:male cells)") +
    scale_shape_manual(name = 'Body part',
                     values = c(19, 1)) +
  scale_fill_manual(guide = 'none',
                    values = c("#287C8E", "#E69F00"),
                    na.translate = FALSE) +
  scale_size(name = 'Total DE genes',
             breaks = c(10, 50, 200, 500, 2000),
             limits = c(10, 2000)) +
  my_theme

## Save plot to file ##
ggsave(file = paste0("fb_mb_genes_vs_female_male_cells_ratios.png"),
       plot = de_genes_vs_cell_abundance,
       width = 210, height = 160, units = c("mm"), dpi = 300)



