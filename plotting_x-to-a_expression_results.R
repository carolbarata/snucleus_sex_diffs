#################################
#      Plot X:A expression      #
#     ratios across tissues     #
#         and cell types        #
#################################

#### Load libraries ####
library(tidyverse)
library(EnvStats)
library(gridExtra)
library(ggridges)
library(see)
library(ggpubr)
library(scales)

## Read in compiled data from csv ##
x_a_ratio_df_plus <- read.csv("x_a_normalised_expression_ratio_n_max_expression_data_relaxed_after_filtering_all_tissues.csv",
                              header = TRUE)

x_a_ratio_df_plus$tissue <- gsub("^body$", 'carcass',
                                 x_a_ratio_df_plus$tissue)

x_a_ratio_df_plus$tissue <- gsub("proboscis and maxillary palps", 'proboscis',
                                 x_a_ratio_df_plus$tissue)

x_a_ratio_df_plus$tissue <- factor(x_a_ratio_df_plus$tissue,
                                   levels = rev(c("gonad", "male reproductive glands",
                                                  "wing", "trachea",
                                                  "proboscis",
                                                  "oenocyte", "malpighian tubule",
                                                  "leg", "heart", "haltere", "gut",
                                                  "fat body", "body wall", "antenna",
                                                  "carcass", "head")))

x_a_ratio_df_plus$sex <- factor(x_a_ratio_df_plus$sex,
                                levels = c('female', 'male'))


x_a_ratio_df_plus$broad_annotation <- gsub("female germline cell",
                                           "germline cell",
                                           x_a_ratio_df_plus$broad_annotation)

x_a_ratio_df_plus$broad_annotation <- gsub("male germline cell",
                                           "germline cell",
                                           x_a_ratio_df_plus$broad_annotation)

x_a_ratio_df_plus$broad_annotation <- gsub("female reproductive system",
                                           "reproductive system",
                                           x_a_ratio_df_plus$broad_annotation)

x_a_ratio_df_plus$broad_annotation <- gsub("male reproductive system",
                                           "reproductive system",
                                           x_a_ratio_df_plus$broad_annotation)


x_a_ratio_df_plus$broad_annotation <- factor(x_a_ratio_df_plus$broad_annotation,
                                             levels = c("neuron",
                                                        "sensory neuron",
                                                        "glial cell",
                                                        "tracheolar cell",
                                                        "salivary gland",
                                                        "cardial cell",
                                                        "hemocyte",
                                                        "excretory system",
                                                        "fat cell",
                                                        "somatic precursor cell",
                                                        "antimicrobial epithelial cell",
                                                        "epithelial cell",
                                                        "oenocyte",
                                                        "muscle cell",
                                                        "unannotated",
                                                        "reproductive system",
                                                        "germline cell",
                                                        "male accessory gland"))


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
                         'oenocyte', 'proboscis and maxillary palps',
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


### Compute sample size per density ###
sample_sizes_per_cluster_tissues <- sample_sizes_per_cluster_per_tissue %>%
  filter(!tissue %in% c('head', 'carcass')) %>%
  ungroup() %>%
  group_by(broad_annotation, sex) %>%
  summarise(sample_size = sum(sample_size))


#### Calculate maximum of each X:A distribution ####
max_per_tissue <- x_a_ratio_df_plus %>%
  group_by(sex,tissue) %>%
  summarise(max_per_tissue = density(normalised_ratio)$x[which.max(density(normalised_ratio)$y)]) %>%
  filter(sex == 'male') %>%
  arrange(max_per_tissue)

max_per_tissue$tissue <- as.character(max_per_tissue$tissue)
max_per_tissue$tissue <- factor(max_per_tissue$tissue,
                                levels = max_per_tissue$tissue)


#### Change order of tissues from least to most compensated ####
x_a_ratio_df_plus$tissue <- factor(x_a_ratio_df_plus$tissue,
                                   levels = max_per_tissue$tissue)

sample_sizes_per_tissue <- x_a_ratio_df_plus %>%
  group_by(tissue, sex) %>%
  summarise(sample_size = n())

##### Plot distributions of X/A expression ratios #####
#### At the tissue level ####
x_a_ratio_df_kstest <- x_a_ratio_df_plus %>%
  select(sex, tissue, normalised_ratio) %>%
  pivot_wider(names_from = c(sex, tissue),
              values_from = normalised_ratio,
              names_sep = "-",
              values_fn = list)

tissues <- max_per_tissue$tissue
tissues <- tissues[!tissues == 'male reproductive glands']

wilcox_results <- data.frame()

for(tissue in tissues){
  new_row <- c(tissue,
               ks.test(x_a_ratio_df_kstest[[paste0('female-', tissue)]][[1]],
                       x_a_ratio_df_kstest[[paste0('male-', tissue)]][[1]])$p.value)
  wilcox_results <- rbind(wilcox_results,
                          new_row)
}

colnames(wilcox_results) <- c('tissue', 'pval')

wilcox_results$pval <- as.numeric(wilcox_results$pval)

wilcox_results$pval_bonf <- p.adjust(wilcox_results$pval,
                                     method = "bonferroni")

wilcox_results <- wilcox_results %>%
  mutate(label = case_when(pval_bonf < 0.001 ~ '***',
                           pval_bonf < 0.01 ~ '**',
                           pval_bonf < 0.05 ~ '*',
                           TRUE ~ 'N.S.'))

wilcox_results <- rbind(c('male reproductive glands',
                          'NA', '', ''),
                        wilcox_results)

wilcox_results$sex <- 'NA'

wilcox_results <- wilcox_results %>%
  mutate(category = case_when(tissue == 'head' ~ 'Whole body parts',
                              tissue == 'carcass' ~ 'Whole body parts',
                              tissue == 'male reproductive glands' ~ 'Reproductive tissues',
                              tissue == 'gonad' ~ 'Reproductive tissues',
                              TRUE ~ 'Somatic tissues'))

wilcox_results$tissue <- as.factor(wilcox_results$tissue)

### Add sample size to tissue as label
sample_sizes_per_tissue_for_label <- sample_sizes_per_tissue %>%
  pivot_wider(names_from = sex, values_from = sample_size)

x_a_ratio_df_plus <- x_a_ratio_df_plus %>%
  left_join(y = sample_sizes_per_tissue_for_label) %>%
  mutate(size_label = case_when(is.na(female) ~ 
                                  paste0('n[M]~`=`~', male),
                                is.na(male) ~ 
                                  paste0('n[F]~`=`~', female),
                                .default = paste0('n[F]~`=`~', female,
                                                  '~~n[M]~`=`~', male))) %>%
  mutate(tissue_label = paste0('bold(', tissue, ')'))


x_a_ratio_df_plus$size_label <- gsub(pattern = ' ',
                                     replacement = '~',
                                     x = x_a_ratio_df_plus$size_label)

x_a_ratio_df_plus$tissue_label <- gsub(pattern = ' ',
                                       replacement = '~',
                                       x = x_a_ratio_df_plus$tissue_label)

# Add tissue_label and size_label to wilcox_results
wilcox_results <- x_a_ratio_df_plus %>%
  group_by(tissue, tissue_label, size_label) %>%
  summarise(count = n()) %>%
  left_join(x = wilcox_results)

# Add tissue_label and size_label to max_per_tissue
max_per_tissue <- x_a_ratio_df_plus %>%
  group_by(tissue, tissue_label, size_label) %>%
  summarise(count = n()) %>%
  left_join(x = max_per_tissue)


# Create groups of tissues for plotting
incomplete_dc <- c('male reproductive glands', 'gonad', 'gut', 'body wall',
                   'fat body', 'heart', 'oenocyte', 'leg')

# Define subsets for max_per_tissue, wilcox_results and x_a_ratio_df_plus
max_per_tissue_inc_dc <- max_per_tissue %>%
  filter(tissue %in% incomplete_dc)

x_a_ratio_df_plus_inc <- x_a_ratio_df_plus %>%
  filter(tissue %in% incomplete_dc) %>%
  mutate(tissue_label =  factor(tissue_label,
                                levels = max_per_tissue_inc_dc$tissue_label),
         size_label =  factor(size_label,
                              levels = max_per_tissue_inc_dc$size_label))

wilcox_results_inc <- wilcox_results %>%
  filter(tissue %in% incomplete_dc) %>%
  mutate(tissue_label =  factor(tissue_label,
                                levels = max_per_tissue_inc_dc$tissue_label),
         size_label =  factor(size_label,
                              levels = max_per_tissue_inc_dc$size_label))



complete_dc <- c('haltere', 'proboscis', 'malpighian tubule', 'trachea',
                 'wing', 'carcass', 'head', 'antenna')

max_per_tissue_comp_dc <- max_per_tissue %>%
  filter(tissue %in% complete_dc)

x_a_ratio_df_plus_comp <- x_a_ratio_df_plus %>%
  filter(tissue %in% complete_dc) %>%
  mutate(tissue_label =  factor(tissue_label,
                                levels = max_per_tissue_comp_dc$tissue_label),
         size_label =  factor(size_label,
                              levels = max_per_tissue_comp_dc$size_label))

wilcox_results_comp <- wilcox_results %>%
  filter(tissue %in% complete_dc) %>%
  mutate(tissue_label =  factor(tissue_label,
                                levels = max_per_tissue_comp_dc$tissue_label),
         size_label =  factor(size_label,
                              levels = max_per_tissue_comp_dc$size_label))

x_a_tissue_plot_inc_dc <- ggplot(data = x_a_ratio_df_plus_inc, # <= mod for inc/comp
                                 aes(y = tissue_label,
                                     x = normalised_ratio,
                                     colour = sex)) +
  ggridges::geom_density_ridges(alpha = 0) +
  scale_x_continuous(breaks = c(0.5, 0.75, 1, 1.25, 1.5),
                     limits = c(0.3, 1.6)) +
  scale_y_discrete(position = 'left') +
  geom_text(data = wilcox_results_inc, # <= mod for inc/comp
            mapping = aes(x = 1.2,
                          y = 5.5,
                          label = label),
            colour = 'black',
            size = 2.3) +
  labs(x = "X:A mean expression",
       y = "Density") +
  scale_colour_manual(name = "Sex",
                      values = c("#287C8E", "#E69F00", "black"),
                      breaks = c('female', 'male')) +
  scale_fill_manual(name = "Sex",
                    values = c("#287C8E", "#E69F00", "black"),
                    breaks = c('female', 'male')) +
  facet_wrap(vars(tissue_label, size_label),
             ncol = 1,
             dir = "v", strip.position = "right",
             scales = 'free',
             drop = TRUE,
             labeller = label_parsed) +
  theme(legend.position = "bottom",
        text = element_text(size = 11),
        strip.text.y = element_text(angle = 0),
        panel.grid.major.x = element_line(color = "#999999",
                                          linetype = 'dashed'),
        panel.background = element_blank(),
        panel.spacing.y = unit(0.5, "lines"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle = 0,
                                         hjust = 0.5,
                                         size = 10),
        strip.text = element_text(size = 8),
        strip.background = element_blank(),
        plot.margin = unit(c(1,0,0.5,0.2), 'cm')) # incomplete DC
        #plot.margin = unit(c(1,0,0,0.2), 'cm')) # complete DC

x_a_tissue_plus_samples <- ggarrange(x_a_tissue_plot_inc_dc,
                                     x_a_tissue_plot_comp_dc,
                                     labels = c('(A) Male X:A < 1',
                                                '(B) Male X:A \u2265 1'),
                                     font.label = list(size = 12),
                                     hjust = -0.2,
                                     nrow = 2,
                                     common.legend = TRUE,
                                     legend = 'bottom')


ggsave(file = "x_to_a_ridges_all_tissues.png",
       plot = x_a_tissue_plus_samples,
       width = 200, height = 297, units = c("mm"), dpi = 300)


#### At the cluster level per tissue ####
sample_sizes_per_cluster_per_tissue <- x_a_ratio_df_plus %>%
  filter(!tissue %in% c('head', 'carcass')) %>%
  group_by(broad_annotation, tissue, sex) %>%
  summarise(sample_size = n())

x_a_ratio_df_plus_subset <- subset(x_a_ratio_df_plus,
                                   !tissue %in% c('head', 'carcass'))

x_a_ratio_df_plus_subset$tissue <- gsub("proboscis and maxillary palps",
                                        "proboscis",
                                        x_a_ratio_df_plus_subset$tissue)

clusters <- as.character(unique(x_a_ratio_df_plus_subset$broad_annotation))

x_a_ratio_df_kstest <- x_a_ratio_df_plus_subset %>%
  select(sex, tissue, broad_annotation, normalised_ratio) %>%
  filter(tissue != 'male reproductive glands') %>%
  pivot_wider(names_from = c(sex, tissue, broad_annotation),
              values_from = normalised_ratio,
              names_sep = "-")

tissues_clusters <- x_a_ratio_df_plus_subset %>%
  select(tissue, broad_annotation) %>%
  filter(tissue != 'male reproductive glands') %>%
  unite(tissue_cluster,
        c("tissue", "broad_annotation"),
        sep = "-") %>%
  distinct(tissue_cluster, .keep_all = TRUE)


for(tissue_cluster in tissues_clusters$tissue_cluster){
  print(tissue_cluster)
  if (paste0('female-', tissue_cluster) %in% names(x_a_ratio_df_kstest) &&
      paste0('male-', tissue_cluster) %in% names(x_a_ratio_df_kstest)) {
    new_row <- c(strsplit(tissue_cluster, '-')[[1]][1],
                 strsplit(tissue_cluster, '-')[[1]][2],
                 ks.test(x_a_ratio_df_kstest[[paste0('female-', tissue_cluster)]][[1]],
                         x_a_ratio_df_kstest[[paste0('male-', tissue_cluster)]][[1]])$p.value)
  } else {
    new_row <- c(strsplit(tissue_cluster, '-')[[1]][1],
                 strsplit(tissue_cluster, '-')[[1]][2],
                 'NA')
  }
  
  if (exists("wilcox_results_c")) {
    wilcox_results_c <- rbind(wilcox_results_c,
                              new_row)
  } else {
    wilcox_results_c <- data.frame(tissue = new_row[1],
                                   broad_annotation = new_row[2],
                                   pval = new_row[3])
    
  }
}


wilcox_results_c <- wilcox_results_c %>%
  filter(!tissue %in% c('head', 'carcass'))

wilcox_results_c$pval <- as.numeric(wilcox_results_c$pval)

wilcox_results_c$pval_bonf <- p.adjust(wilcox_results_c$pval,
                                       method = "bonferroni")

wilcox_results_c <- wilcox_results_c %>%
  mutate(label = case_when(pval_bonf < 0.001 ~ '***',
                           pval_bonf < 0.01 ~ '**',
                           pval_bonf < 0.05 ~ '*',
                           is.na(pval_bonf) ~ '',
                           TRUE ~ 'N.S.'))

wilcox_results_c$tissue <- factor(wilcox_results_c$tissue,
                                  levels = rev(c("gonad", "male reproductive glands",
                                                 "wing", "trachea",
                                                 "proboscis",
                                                 "oenocyte", "malpighian tubule",
                                                 "leg", "heart", "haltere", "gut",
                                                 "fat body", "body wall", "antenna",
                                                 "body", "head")))

clusters <- clusters[!clusters == 'unannotated']

plot_l <- list()
for (cluster in clusters) {
  if (cluster == 'male accessory gland') {
    wilcox_results_cluster <- data.frame(tissue = 'male reproductive glands',
                                         broad_annotation = 'male accessory gland',
                                         pval = '',
                                         pval_bonf = '',
                                         label = '')
    samples_cluster <- sample_sizes_per_cluster_per_tissue %>%
      filter(!tissue %in% c('head', 'carcass')) %>%
      filter(broad_annotation == cluster)
  } else {
    wilcox_results_cluster <- wilcox_results_c %>%
      filter(broad_annotation == cluster) %>%
      select(tissue, broad_annotation, pval, pval_bonf, label)
    
    samples_cluster <- sample_sizes_per_cluster_per_tissue %>%
      filter(!tissue %in% c('head', 'carcass')) %>%
      filter(broad_annotation == cluster)
  }
  
  if (cluster == 'germline cell') {
    cluster_subset <- subset(x_a_ratio_df_plus_subset,
                             broad_annotation == cluster &
                               tissue != 'male reproductive glands')
  } else {
    cluster_subset <- subset(x_a_ratio_df_plus_subset,
                             broad_annotation == cluster)
  }
  
  #### Add info to keep legend only if the cluster is at the bottom ####
  keep_legend_in <- c('excretory system', 'oenocyte', 'sensory neuron',
                      'antimicrobial epithelial cell', 'reproductive system',
                      'unannotated')
  if (cluster %in% keep_legend_in) {
    legend_status <- 'legend'
  } else {
    legend_status <- 'none'
  }
  
  sample_sizes_cluster_for_label <- samples_cluster %>%
    pivot_wider(id_cols = c('broad_annotation', 'tissue'),
                names_from = sex,
                values_from = sample_size)
  
  if (!'female' %in% colnames(sample_sizes_cluster_for_label)) {
    sample_sizes_cluster_for_label$female <- 0
  }
  if (!'male' %in% colnames(sample_sizes_cluster_for_label)) {
    sample_sizes_cluster_for_label$male <- 0
  }
  
  sample_sizes_cluster_for_label[is.na(sample_sizes_cluster_for_label)] <- 0
  
  cluster_subset <- cluster_subset %>%
    left_join(y = sample_sizes_cluster_for_label) %>%
    mutate(size_label = case_when(female == 0 ~ 
                                    paste0('n[M]~`=`~', male),
                                  male == 0 ~ 
                                    paste0('n[F]~`=`~', female),
                                  .default = paste0('n[F]~`=`~', female,
                                                    '~~n[M]~`=`~', male))) %>%
    mutate(tissue_label = paste0('bold(', tissue, ')'))
  
  
  cluster_subset$size_label <- gsub(pattern = ' ',
                                    replacement = '~',
                                    x = cluster_subset$size_label)
  
  cluster_subset$tissue_label <- gsub(pattern = ' ',
                                      replacement = '~',
                                      x = cluster_subset$tissue_label)
  
  cluster_subset <- left_join(x = cluster_subset,
                              y = wilcox_results_cluster,
                              by = c('tissue', 'broad_annotation'),
                              relationship = 'many-to-one')  
  
  title_string <- str_to_sentence(cluster)
  
  cluster_subset <- cluster_subset %>%
    filter(female + male > 20)
  
  x_a_cluster_plot <- ggplot(data = cluster_subset,
                             aes(y = broad_annotation,
                                 x = normalised_ratio,
                                 colour = sex)) +
    ggridges::geom_density_ridges(alpha = 0) +
    scale_x_continuous(breaks = c(0.5, 0.75, 1, 1.25, 1.5),
                       limits = c(0.3, 1.8)) +
    scale_y_discrete(position = 'left') +
    geom_text(mapping = aes(x = 0.70,
                            y = 4,
                            label = label),
              colour = 'black',
              size = 2.3,
              check_overlap = TRUE) +
    labs(x = "X:A mean expression",
         y = "Density",
         title = title_string) +
    scale_colour_manual(guide = legend_status,
                        name = 'Sex',
                        values = c("#287C8E", "#E69F00", "black"),
                        breaks = c('female', 'male')) +
    facet_wrap(vars(tissue_label, size_label),
               ncol = 1,
               dir = "v", strip.position = "right",
               scales = 'free',
               drop = TRUE,
               labeller = label_parsed) +
    theme(legend.position = "bottom",
          text = element_text(size = 10),
          strip.text.y = element_text(angle = 0,
                                      size = 7.5),
          panel.grid.major.x = element_line(color = "#999999",
                                            linetype = 'dashed'),
          panel.background = element_blank(),
          panel.spacing.y = unit(0.25, "lines"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.y.left = element_text(angle = 0,
                                           hjust = 0.5,
                                           size = 9),
          strip.background = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.2, 0.1), "cm"),
          plot.title = element_text(size = 9,
                                    face = 'bold'))
  
  cluster <- gsub(' ', '_', cluster)
  plot_l[[cluster]] <- x_a_cluster_plot
  rm(x_a_cluster_plot)
  rm(cluster_subset)
  rm(wilcox_results_cluster)
  rm(samples_cluster)
  rm(sample_sizes_cluster_for_label)
  gc()
}

## Circulatory & respiratory & excretory system ##
title_plot1 <- as_ggplot(text_grob("Circulatory, respiratory and excretory system",
                                   color = '#999999',
                                   face = "bold",
                                   size = 15,
                                   just = c(0.5,1))) + 
  theme(plot.margin = margin(0, 0, 0.4, 0, "cm"))
circ_resp_exc_plots <- ggarrange(title_plot1,
                                 plot_l[['hemocyte']],
                                 plot_l[['cardial_cell']],
                                 plot_l[['tracheolar_cell']],
                                 plot_l[['excretory_system']],
                                 common.legend = TRUE,
                                 legend = 'bottom',
                                 heights = c(0.1, 1.1, 0.25, 0.4, 0.25),
                                 ncol = 1)

## Endocrine & immune system ##
title_plot2 <- as_ggplot(text_grob("Endocrine and immune system",
                                   color = '#999999',
                                   face = "bold",
                                   size = 15,
                                   just = c(0.5,1))) + 
  theme(plot.margin = margin(0, 0, 0.4, 0, "cm"))
end_immun_plots <- ggarrange(title_plot2,
                             plot_l[['glial_cell']],
                             plot_l[['salivary_gland']],
                             plot_l[['fat_cell']],
                             plot_l[['oenocyte']],
                             common.legend = TRUE,
                             legend = 'bottom',
                             heights = c(0.1, 0.9, 0.3, 1.1, 0.8),
                             ncol = 1)


## Nervous system ##
title_plot3 <- as_ggplot(text_grob("Nervous system",
                                   color = '#999999',
                                   face = "bold",
                                   size = 15,
                                   just = c(0.5,1))) + 
  theme(plot.margin = margin(0, 0, 0.4, 0, "cm"))
nervous_plots <- ggarrange(title_plot3,
                           plot_l[['neuron']],
                           plot_l[['sensory_neuron']],
                           labels = c('', ''),
                           common.legend = TRUE,
                           legend = 'bottom',
                           heights = c(0.1, 1, 1),
                           ncol = 1)

## Other somatic tissues ##
title_plot4 <- as_ggplot(text_grob("Other somatic tissues",
                                   color = '#999999',
                                   face = "bold",
                                   size = 16,
                                   just = c(0.5,1))) + 
  theme(plot.margin = margin(0, 0, 0.4, 0, "cm"))
other_plots <- ggarrange(title_plot4,
                         ggarrange(plot_l[['muscle_cell']],
                                   plot_l[['epithelial_cell']],
                                   plot_l[['somatic_precursor_cell']],
                                   plot_l[['antimicrobial_epithelial_cell']],
                                   ncol = 2,
                                   nrow = 2,
                                   heights = c(1, 0.25),
                                   common.legend = TRUE,
                                   legend = 'bottom'),
                         nrow = 2,
                         heights = c(0.1, 1),
                         common.legend = TRUE,
                         legend = 'bottom')

## Reproductive system ##
title_plot5 <- as_ggplot(text_grob("Reproductive system",
                                   color = '#999999',
                                   face = "bold",
                                   size = 15,
                                   just = c(0.5,1))) + 
  theme(plot.margin = margin(0, 0, 0.4, 0, "cm"))
reprod_plots <- ggarrange(title_plot5,
                          plot_l[['germline_cell']],
                          plot_l[['reproductive_system']],
                          plot_l[['male_accessory_gland']],
                          common.legend = TRUE,
                          legend = 'bottom',
                          heights = c(0.3, 1, 1, 1.1),
                          ncol = 1)

## Save all plots as pdf ##
ggsave(file = "x_to_a_ridges_and_samples_circ_resp_exc_system.png",
       plot = circ_resp_exc_plots,
       bg = 'white',
       width = 160, height = 240, units = c("mm"), dpi = 300)

ggsave(file = "x_to_a_ridges_and_samples_endocrine_immune_system.png",
       plot = end_immun_plots,
       bg = 'white',
       width = 160, height = 320, units = c("mm"), dpi = 300)

ggsave(file = "x_to_a_ridges_and_samples_nervous_system.pdf",
       plot = nervous_plots,
       bg = 'white',
       width = 160, height = 210, units = c("mm"), dpi = 300)

ggsave(file = "x_to_a_ridges_and_samples_other_somatic.pdf",
       plot = other_plots,
       bg = 'white',
       width = 260, height = 320, units = c("mm"), dpi = 300)

ggsave(file = "x_to_a_ridges_and_samples_reproductive_system.pdf",
       plot = reprod_plots,
       bg = 'white',
       width = 160, height = 140, units = c("mm"), dpi = 300)


#### Plot X:A for different male germline cell types ####
### Add UMAP to Seurat object ###
## Read in anndata stringent data ##
testis_anndata <- read_h5ad('./fly_cell_atlas_proj/data/h5ad_initial_files/s_fca_biohub_testis_10x.h5ad',
                            backed = 'r')

## Extract UMAP coordinates ##
testis_umap <- as.data.frame(testis_anndata$obsm[['X_umap']])

rownames(testis_umap) <- rownames(testis_anndata$obs)

testis_umap$cellID <- rownames(testis_umap)

rm(testis_anndata)
gc()

## Save testis UMAP data to file ##
write.csv(testis_umap,
          file = "s_testis_umap_reduction_coordinates.csv",
          row.names = FALSE)

## Read in testis UMAP ##
testis_umap <- read.csv("s_testis_umap_reduction_coordinates.csv",
                        header = TRUE)

rownames(testis_umap) <- testis_umap$cellID

### Read in testis h5Seurat ###
testis_relaxed_seurat <- LoadH5Seurat("r_testis_cell_type_filt_decontx_stringent_filt_sctransformed.h5Seurat")

Idents(testis_relaxed_seurat) <- 's_annotation'

### Add X:A mean expression as metadata ###
testis_relaxed_seurat <- AddMetaData(object = testis_relaxed_seurat,
                                     metadata = subset(x_a_ratio_df_plus,
                                                       tissue == 'gonad' &
                                                         sex == 'male')$normalised_ratio,
                                     col.name = 'normalised_ratio')


# Only keep cells that are present in Seurat file
testis_umap <- subset(testis_umap,
                      cellID %in% testis_relaxed_seurat@meta.data$cellID)

testis_relaxed_seurat <- subset(x = testis_relaxed_seurat,
                                subset = cellID %in% subset(x_a_ratio_df_plus,
                                                            tissue == 'gonad' &
                                                              sex == 'male')$cell_id,
                                invert = FALSE)

testis_relaxed_seurat[['umap']] <- CreateDimReducObject(embeddings = as.matrix(testis_umap[,c(1,2)]),
                                                        key = 'umap_',
                                                        global = TRUE,
                                                        assay = 'SCT')

FeaturePlot(subset(testis_relaxed_seurat,
                   s_broad_annotation == 'male germline cell'),
            features = 'normalised_ratio',
            reduction = "umap",
            label = TRUE,
            repel = TRUE,
            label.size = 2.5,
            pt.size = 1.2,
            order = TRUE,
            min.cutoff = 0,
            max.cutoff = 1.2,
            cols = c('#6b008a', '#E69F00')) +
  ggtitle("mean X:A") +
  labs(x = "UMAP 1",
       y = "UMAP 2") +
  theme_minimal() +
  theme(text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5,
                                  face = 'bold'))

### Plot germline X:A levels ###
## Subset germline Seurat ##
germcyst_seurat <- subset(x = testis_relaxed_seurat,
                          subset = s_broad_annotation %in% c('male germline cell', 
                                                             'male reproductive system'))

germcyst_seurat <- subset(x = germcyst_seurat,
                          subset = s_annotation != 'male gonad associated epithelium')

germcyst_seurat <- subset(x = germcyst_seurat,
                          subset = s_annotation != 'secretory cell of the male reproductive tract')

x_a_ratio_df_germline <- subset(x_a_ratio_df_plus,
                                broad_annotation == 'germline cell' &
                                  sex == 'male' &
                                  tissue == 'gonad' |
                                  broad_annotation == 'reproductive system' &
                                  sex == 'male' &
                                  tissue == 'gonad')

x_a_ratio_df_germline <- x_a_ratio_df_germline %>%
  filter(s_annotation != 'secretory cell of the male reproductive tract') %>%
  filter(s_annotation != 'male gonad associated epithelium')

x_a_ratio_df_germline$s_annotation <- gsub('spermatocyte 0',
                                           'spermatocyte 0-7',
                                           x_a_ratio_df_germline$s_annotation)

x_a_ratio_df_germline$s_annotation <- gsub('spermatocyte 1',
                                           'spermatocyte 0-7',
                                           x_a_ratio_df_germline$s_annotation)

x_a_ratio_df_germline$s_annotation <- gsub('spermatocyte 2',
                                           'spermatocyte 0-7',
                                           x_a_ratio_df_germline$s_annotation)

x_a_ratio_df_germline$s_annotation <- gsub('spermatocyte 3',
                                           'spermatocyte 0-7',
                                           x_a_ratio_df_germline$s_annotation)

x_a_ratio_df_germline$s_annotation <- gsub('spermatocyte 4',
                                           'spermatocyte 0-7',
                                           x_a_ratio_df_germline$s_annotation)

x_a_ratio_df_germline$s_annotation <- gsub('spermatocyte 5',
                                           'spermatocyte 0-7',
                                           x_a_ratio_df_germline$s_annotation)

x_a_ratio_df_germline$s_annotation <- gsub('spermatocyte 6',
                                           'spermatocyte 0-7',
                                           x_a_ratio_df_germline$s_annotation)

x_a_ratio_df_germline$s_annotation <- gsub('spermatocyte 7a',
                                           'spermatocyte 0-7',
                                           x_a_ratio_df_germline$s_annotation)

x_a_ratio_df_germline$s_annotation <- gsub('early cyst cell 1',
                                           'early cyst cell',
                                           x_a_ratio_df_germline$s_annotation)

x_a_ratio_df_germline$s_annotation <- gsub('early cyst cell 2',
                                           'early cyst cell',
                                           x_a_ratio_df_germline$s_annotation)

x_a_ratio_df_germline$s_annotation <- gsub('spermatocyte cyst cell branch a',
                                           'spermatocyte cyst cell',
                                           x_a_ratio_df_germline$s_annotation)

x_a_ratio_df_germline$s_annotation <- gsub('spermatocyte cyst cell branch b',
                                           'spermatocyte cyst cell',
                                           x_a_ratio_df_germline$s_annotation)

x_a_ratio_df_germline$s_annotation <- gsub('cyst cell branch a',
                                           'cyst cell',
                                           x_a_ratio_df_germline$s_annotation)

x_a_ratio_df_germline$s_annotation <- gsub('cyst cell branch b',
                                           'cyst cell',
                                           x_a_ratio_df_germline$s_annotation)

x_a_ratio_df_germline$s_annotation <- gsub('late cyst cell branch a',
                                           'late cyst cell',
                                           x_a_ratio_df_germline$s_annotation)

x_a_ratio_df_germline$s_annotation <- gsub('late cyst cell branch b',
                                           'late cyst cell',
                                           x_a_ratio_df_germline$s_annotation)

x_a_ratio_df_germline$s_annotation <- gsub('spermatogonium-spermatocyte transition',
                                           'spermatogonium to spermatocyte transition',
                                           x_a_ratio_df_germline$s_annotation)


x_a_ratio_df_germline$s_annotation <- factor(x_a_ratio_df_germline$s_annotation,
                                             levels = c('spermatogonium',
                                                        'mid-late proliferating spermatogonia',
                                                        'spermatogonium to spermatocyte transition',
                                                        'spermatocyte 0-7',
                                                        'spermatocyte',
                                                        'late primary spermatocyte',
                                                        'spermatid',
                                                        'early elongation stage spermatid',
                                                        'early-mid elongation-stage spermatid',
                                                        'mid-late elongation-stage spermatid',
                                                        'early cyst cell',
                                                        'cyst cell intermediate',
                                                        'late cyst cell',
                                                        'cyst cell',
                                                        'spermatocyte cyst cell',
                                                        'head cyst cell',
                                                        'male gonad pigment cell'))


germline_boxplot <- ggplot(data = subset(x_a_ratio_df_germline,
                                         broad_annotation == 'germline cell'),
                           aes(x = s_annotation,
                               y = normalised_ratio,
                               fill = broad_annotation)) +
  geom_boxplot(notch = TRUE,
               outlier.alpha = 0.6,
               outlier.shape = 95,
               outlier.size = 5,
               alpha = 0.75) +
  labs(y = "X:A mean expression per cell",
       x = "") +
  scale_fill_manual(values = '#e69f00') +
  scale_y_continuous(breaks = c(0.5, 0.75, 1),
                     limits = c(0.3, 2.5)) +
  scale_x_discrete(labels = label_wrap(10)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 0.5, vjust = 0.65, 
                                   size = 6),
        text = element_text(size = 10),
        panel.grid.major.y = element_line(color = "#999999",
                                          linetype = 'dashed'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))

cystcell_boxplot <- ggplot(data = subset(x_a_ratio_df_germline,
                                         broad_annotation == 'reproductive system'),
                           aes(x = s_annotation,
                               y = normalised_ratio,
                               fill = broad_annotation)) +
  geom_boxplot(notch = TRUE,
               outlier.alpha = 0.6,
               outlier.shape = 95,
               outlier.size = 5,
               alpha = 0.75) +
  labs(y = "X:A mean expression per cell",
       x = "") +
  scale_fill_manual(values = '#288e60') +
  scale_y_continuous(breaks = c(0.5, 0.75, 1),
                     limits = c(0.3, 2.5)) +
  scale_x_discrete(labels = label_wrap(10)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 0.5, vjust = 0.65, 
                                   size = 8),
        text = element_text(size = 10),
        panel.grid.major.y = element_line(color = "#999999",
                                          linetype = 'dashed'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        plot.margin = unit(c(0,2,0,2), "cm"))

testis_combi_plot <- ggarrange(germline_boxplot,
                               cystcell_boxplot,
                               labels = c('(A)', '(B)'),
                               nrow = 2,
                               widths = c(1, 0.7),
                               heights = c(1, 0.9))

ggsave(file = "x_to_a_boxplots_germline_and_cyst_cells_in_testis.pdf",
       plot = testis_combi_plot,
       width = 180, height = 160, units = c("mm"), dpi = 300)

