##############################
#Plot sex bias
# with and w/o
# cell type composition
##################################

#### Load libraries ####
library(tidyverse)
library(ggpubr)
library(broom)

#### Read in data on cells per cluster ####
cells_df <- read.csv(file = 'number_of_cells_per_cluster_all_tissues.csv',
                     header = TRUE)

cells_df <- cells_df %>%
  group_by(tissue, sex) %>%
  mutate(total_cells_tissue = sum(cells))

cells_df <- cells_df %>%
  pivot_wider(id_cols = c('s_broad_annotation', 'tissue'),
              names_from = sex,
              values_from = c(cells, total_cells_tissue))

cells_df[is.na(cells_df)] <- 0

cells_df <- cells_df %>%
  filter(cells_female + cells_male > 20)

cells_df <- cells_df %>%
  mutate(total_cluster_cells = cells_female + cells_male,
         prop_female = cells_female/total_cells_tissue_female,
         prop_male = cells_male/total_cells_tissue_male,
         fminus_m = prop_female - prop_male,
         f2m_ratio = prop_female/prop_male,
         log2_f2m = log2(f2m_ratio))

cells_df <- cells_df %>%
  mutate(bias_status = case_when(fminus_m == 0 ~ 'unbiased',
                                 fminus_m > 0 ~ 'female-biased',
                                 fminus_m < 0 ~ 'male-biased'))

cells_df$tissue <- gsub('^body$', 'carcass',
                        cells_df$tissue)

cells_df$tissue <- gsub('proboscis_and_maxillary_palps',
                        'proboscis',
                        cells_df$tissue)

#### Read in sex bias data ####
sex_bias_df <- read.csv("sex_bias_cell_vs_no_cell_per_tissue_plus.csv",
                        header = TRUE)

#### Replace 'body' with 'carcass' ####
sex_bias_df$tissue <- gsub("^body$", 'carcass', sex_bias_df$tissue)

#### Plot tissue-level bias ####
### Order tissues according to head & carcass distn ###
# sex_bias_df$tissue <- gsub("_", " ",
#                            sex_bias_df$tissue)

### Change 'proboscis and maxillary palps' to 'proboscis' ###
sex_bias_df$tissue <- gsub('proboscis_and_maxillary_palps',
                           'proboscis',
                           sex_bias_df$tissue)

### Set colour status for head & carcass tissues ###
sex_bias_df <- sex_bias_df %>%
  mutate(colour_status = case_when(tissue %in% c('head', 'antenna',
                                                 'proboscis') ~
                                     'head colour',
                                   .default = 'carcass colour'))




### Set sex bias level variable
sex_bias_df <- sex_bias_df %>%
  mutate(log_mean_sb_comp_v1 = log2(mean_sb_cell_comp_v1),
         log_mean_sb_comp_v2 = log2(mean_sb_cell_comp_v2),
         log_mean_sb_no_comp_v1 = log2(mean_sb_no_cell_comp_v1),
         log_mean_sb_no_comp_v2 = log2(mean_sb_no_cell_comp_v2),
         log_sum_sb_comp_v1 = log2(sum_sb_cell_comp_v1),
         log_sum_sb_comp_v2 = log2(sum_sb_cell_comp_v2),
         log_sum_sb_no_comp_v1 = log2(sum_sb_no_cell_comp_v1),
         log_sum_sb_no_comp_v2 = log2(sum_sb_no_cell_comp_v2))


#### Set plot theme ####
my_theme <- theme(legend.position = "bottom",
                  text = element_text(size = 10),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.spacing = unit(2, "mm"),
                  panel.background = element_blank(),
                  axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                  axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))

### Set colour name vector
cols <- c("female-biased" = "#287C8E",
          "male-biased" = "#E69F00",
          "unbiased" = "#999999",
          "compensatory" = "purple",
          "cell type" = "pink",
          "head colour" = '#288e60',
          "carcass colour" = '#905395')



tissues <- unique(cells_df$tissue)

plot_l <- list()
for (tissue_id in tissues) {
  ## Subset per tissue ##
  cells_df_subset <- cells_df %>%
    filter(tissue == tissue_id)
  
  ### Arrange by cluster size ###
  cells_df_subset <- cells_df_subset %>%
    arrange(total_cluster_cells)
  
  cells_df_subset$s_broad_annotation <- as.factor(cells_df_subset$s_broad_annotation)
  
  cells_df_subset <- transform(cells_df_subset,
                               s_broad_annotation = reorder(s_broad_annotation,
                                                            -total_cluster_cells) ) 
  
  
  ### Plot number of female and male cells per tissue ###
  cells_plots <- ggplot(data = cells_df_subset,
                        aes(x = fminus_m,
                            y = s_broad_annotation,
                            fill = bias_status)) +
    geom_bar(position = "stack",
             stat = "identity") +
    geom_vline(xintercept = 0,
               colour = "#4C4C4C") +
    xlim(c(-0.2, 0.2)) +
    labs(y = '',
         x = 'F - M cluster proportion') +
    scale_fill_manual(values = cols,
                      guide = 'none') +
    my_theme +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(1, 0, 0.1, 0),
                             "cm"),
          axis.text.y = element_text(face = "bold", size = 9.5))
  
  cells_total_plot <- ggplot(data = cells_df_subset,
                             aes(x = tissue,
                                 y = s_broad_annotation)) +
    geom_point(aes(size = total_cluster_cells),
               colour = "#4C4C4C") +
    geom_text(aes(label = total_cluster_cells),
              colour = "#4C4C4C",
              size = 3,
              hjust = -0.5) +
    labs(x = '',
         y = '') +
    scale_size(guide = 'none') +
    theme_void() +
    theme(plot.margin = unit(c(1, 0.5, 1, 0), 
                             "cm"))
  
  cells_diffs_plots <- ggarrange(cells_total_plot,
                                 cells_plots,
                                 nrow = 1,
                                 widths = c(0.4, 1))
  
  
  ### Plot sex bias with cell composition vs no cell composition differences
  ## Any tissue
  hnc_tissues <- sex_bias_df %>%
    filter(tissue == tissue_id)
  
  regressions_hnc <- hnc_tissues %>%
    nest(data = -tissue) %>%
    mutate(fit = map(data,
                     ~ lm(log_mean_sb_no_comp_v1 ~ log_mean_sb_comp_v1, 
                          data = (.))),
           augmented = map(fit, augment,
                           interval = "prediction",
                           conf.level = 0.95),
           outliers = map(fit, car::outlierTest,
                          n.max = Inf))
  
  names(regressions_hnc$outliers) <- regressions_hnc$tissue
  
  names(regressions_hnc$data) <- regressions_hnc$tissue
  
  for (tissue_id2 in regressions_hnc$tissue) {
    outlier_index <- data.frame(rownames = rownames_to_column(as.data.frame(regressions_hnc$outliers[[tissue_id2]]$bonf.p)),
                                tissue = tissue_id2)
    colnames(outlier_index) <- c('rownames', 'bonf_p', 'tissue')
    
    outlier_index$gene <- as.data.frame(regressions_hnc$data[[tissue_id2]])[outlier_index$rownames,'gene']
    
    if (!exists("hnc_outlier_index")) {
      hnc_outlier_index <- outlier_index
    } else {
      hnc_outlier_index <- rbind(outlier_index, hnc_outlier_index)
    }
    rm(outlier_index)
    gc()
  }
  
  
  hnc_tissues <- left_join(hnc_tissues,
                           hnc_outlier_index) %>%
    mutate(outlier_status = case_when(!is.na(rownames) ~ 'outlier',
                                      .default = 'non-outlier'))
  
  hnc_tissues_outlier_count <- hnc_tissues %>%
    group_by(tissue, outlier_status) %>%
    summarise(outlier_count = n()) %>%
    filter(outlier_status == 'outlier') %>%
    select(tissue, outlier_count) %>%
    mutate(colour_status_label = case_when(tissue %in% c('head', 'antenna',
                                                         'proboscis') ~
                                             'head colour',
                                           .default = 'carcass colour'))
  
  
  sb_cell_vs_no_cell_hnc <- ggplot(data = hnc_tissues,
                                   aes(x = log_mean_sb_comp_v1,
                                       y = log_mean_sb_no_comp_v1)) +
    geom_hline(yintercept = 1,
               linetype = 'dashed',
               colour = "#999999") +
    geom_hline(yintercept = -1,
               linetype = 'dashed',
               colour = "#999999") +
    geom_vline(xintercept = 1,
               linetype = 'dashed',
               colour = "#999999") +
    geom_vline(xintercept = -1,
               linetype = 'dashed',
               colour = "#999999") +
    geom_point(aes(colour = colour_status),
               alpha = 0.45) +
    geom_point(aes(shape = outlier_status,
                   colour = colour_status),
               alpha = 0.45) +
    geom_line(data = regressions_hnc %>% unnest(augmented),
              aes(y = .fitted),
              colour = "#999999") +
    geom_ribbon(data = regressions_hnc %>% unnest(augmented),
                aes(ymin = .lower, ymax = .upper),
                alpha = 0.3,
                fill = "#999999",
                colour = NA) +
    geom_label(data = hnc_tissues_outlier_count,
               aes(label = paste('outliers = ', outlier_count),
                   colour = colour_status_label),
               x = Inf, y = -Inf,
               vjust = -0.8, hjust = 1.2,
               label.size = NA) +
    labs(x = 'log2(sex bias with cell type composition)',
         y = 'log2(sex bias no cell type composition)') +
    scale_colour_manual(guide = 'none',
                        values = cols) +
    scale_shape_manual(guide = 'none',
                       values = c('outlier' = 4,
                                  'non-outlier' = 19)) +
    my_theme +
    theme(axis.text = element_text(size = 7),
          plot.margin = unit(c(0.6, 0.2, 0.2, 0.2),
                             "cm"))
  
  combined_plots <- ggarrange(cells_diffs_plots,
                              sb_cell_vs_no_cell_hnc,
                              nrow = 1,
                              widths = c(1, 0.8))
  
  plot_l[[tissue_id]] <- combined_plots
  rm(cells_df_subset)
  rm(cells_plots)
  rm(cells_total_plot)
  rm(cells_diffs_plots)
  rm(hnc_tissues)
  rm(hnc_outlier_index)
  rm(hnc_tissues_outlier_count)
  rm(regressions_hnc)
  rm(sb_cell_vs_no_cell_hnc)
  rm(combined_plots)
  gc()
}


### Save plot with most biased tissues ###
sb_comp_vs_no_comp <- ggarrange(plot_l[['carcass']],
                                plot_l[['fat_body']],
                                plot_l[['proboscis']],
                                plot_l[['oenocyte']],
                                nrow = 4,
                                labels = c('(A) Carcass',
                                           '(B) Fat body',
                                           '(C) Proboscis',
                                           '(D) Oenocyte'))

ggsave(file = "tissue_mean_sex_bias_cell_comp_vs_no_comp_all_genes_most_biased.png",
       plot = sb_comp_vs_no_comp,
       width = 240, height = 297, units = c("mm"), dpi = 300,
       bg = 'white')

### Save plot with remaining head tissues ###
sb_comp_vs_no_comp <- ggarrange(plot_l[['head']],
                                plot_l[['antenna']],
                                nrow = 2,
                                labels = c('(A) Head',
                                           '(B) Antenna'))

ggsave(file = "tissue_mean_sex_bias_cell_comp_vs_no_comp_all_genes_other_head.pdf",
       plot = sb_comp_vs_no_comp,
       width = 240, height = 160, units = c("mm"), dpi = 300,
       bg = 'white')

### Save plot with remaining internal organs ###
sb_comp_vs_no_comp <- ggarrange(plot_l[['heart']],
                                plot_l[['trachea']],
                                plot_l[['malpighian_tubule']],
                                plot_l[['gut']],
                                nrow = 4,
                                labels = c('(A) Heart',
                                           '(B) Trachea',
                                           '(C) Malpighian tubule',
                                           '(D) Gut'))

ggsave(file = "tissue_mean_sex_bias_cell_comp_vs_no_comp_all_genes_other_organs.png",
       plot = sb_comp_vs_no_comp,
       width = 240, height = 297, units = c("mm"), dpi = 300,
       bg = 'white')

### Save plot with remaining tissues ###
sb_comp_vs_no_comp <- ggarrange(plot_l[['body_wall']],
                                plot_l[['haltere']],
                                plot_l[['leg']],
                                plot_l[['wing']],
                                nrow = 4,
                                labels = c('(A) Body wall',
                                           '(B) Haltere',
                                           '(C) Leg',
                                           '(D) Wing'))

ggsave(file = "tissue_mean_sex_bias_cell_comp_vs_no_comp_all_genes_other.png",
       plot = sb_comp_vs_no_comp,
       width = 240, height = 297, units = c("mm"), dpi = 300,
       bg = 'white')


