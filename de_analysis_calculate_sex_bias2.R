##############################################
#             Calculate sex bias             #
#         per cluster and per tissue         #
##############################################

#### Load libraries ####
library(Seurat)
library(SeuratDisk)
library(tidyverse)

#### Set working directory ####
setwd("./fca_relaxed_loom/")

#### Loop over all tissue samples ####
tissues <- c('Antenna', 'Malpighian_tubule', 'Body_wall',
             'Oenocyte', 'Fat_body', 'Gut',
             'Proboscis_and_maxillary_palps', 'Haltere',
             'Heart', 'Trachea', 'Leg',
             'Wing', 'Head', 'Body')


#### Get expression data ####
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
  
  ### Get cluster names ###
  tissue_clusters <- unique(tissue_seurat_final@meta.data$s_broad_annotation)
  
  ### Get summary data per cluster ###
  for (cluster in tissue_clusters){
    ## Create subset per cluster ##
    cluster_subset <- subset(tissue_seurat_final,
                             subset = s_broad_annotation == cluster)
    
    ## Get counts data ##
    counts_data <- GetAssayData(object = cluster_subset,
                                layer = 'counts',
                                #slot = 'counts',
                                assay = 'SCT')
    
    metadata_df <- cluster_subset@meta.data[, c('cellID', 
                                                'sex',
                                                's_broad_annotation')]
    
    rm(cluster_subset)
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
    
    
    #### Match cell ID with metadata ####
    counts_data <- left_join(counts_data, metadata_df)
    rm(metadata_df)
    gc()
    
    ### Calculate summary stats for females and males ###
    counts_data_summary <- counts_data %>%
      group_by(gene, sex) %>%
      summarise(median = median(counts, na.rm = TRUE),
                mean = mean(counts, na.rm = TRUE),
                sum = sum(counts, na.rm = TRUE),
                variance = var(counts, na.rm = TRUE),
                cells = n())
    
    rm(counts_data)
    gc()
    
    counts_data_summary <- counts_data_summary %>%
      mutate(std_dev = sqrt(variance),
             cvar = sqrt(variance)/mean)
    
    ### Split female and male data across cols ###
    counts_data_summary_mod <- counts_data_summary %>%
      filter(sex == 'female')
    
    counts_data_summary_m <- counts_data_summary %>%
      filter(sex == 'male')
    
    counts_data_summary_mod <- left_join(counts_data_summary_mod,
                                         counts_data_summary_m,
                                         by = c('gene'))
    
    rm(counts_data_summary)
    rm(counts_data_summary_m)
    gc()

     counts_data_summary_mod$cluster <- cluster
     counts_data_summary_mod$tissue <- tolower(tissue)
    
    if (exists("tissue_bias_df")) {
      tissue_bias_df <- rbind(tissue_bias_df,
                              counts_data_summary_mod)
    } else {
      tissue_bias_df <- counts_data_summary_mod
    }
    
    rm(counts_data_summary_mod)
    gc()
  }
  
  rm(tissue_seurat_final)
  rm(tissue_clusters)
  gc()
  
  ncells_f <- tissue_bias_df %>%
    group_by(cluster) %>%
    summarise(cluster_fcells = unique(cells.x)) %>%
    ungroup() %>%
    summarise(tissue_fcells = sum(cluster_fcells))
  
  ncells_m <- tissue_bias_df %>%
    group_by(cluster) %>%
    summarise(cluster_mcells = unique(cells.y)) %>%
    ungroup() %>%
    summarise(tissue_mcells = sum(cluster_mcells))
    
  tissue_bias_df <- tissue_bias_df %>%
    mutate(prop.x = cells.x/ncells_f$tissue_fcells,
           prop.y = cells.y/ncells_m$tissue_mcells,
           prop.xy = 0.5*(prop.x + prop.y))
  
  sb_cell_comp <- tissue_bias_df %>%
    mutate(mean_prop.x = mean.x*prop.x,
           mean_prop.y = mean.y*prop.y,
           sum_prop.x = sum.x*prop.x,
           sum_prop.y = sum.x*prop.y) %>%
    group_by(gene) %>%
    summarise(mean_ratio_cell_v1 = sum(mean_prop.x + 10e-9)/sum(mean_prop.y + 10e-9),
              mean_ratio_cell_v2 = sum(mean_prop.x + 0.1)/sum(mean_prop.y + 0.1),
              mean_sb_cell_comp_v1 = mean_ratio_cell_v1*(ncells_f$tissue_fcells/ncells_m$tissue_mcells),
              mean_sb_cell_comp_v2 = mean_ratio_cell_v2*(ncells_f$tissue_fcells/ncells_m$tissue_mcells),
              sum_ratio_cell_v1 = sum(sum_prop.x + 10e-9)/sum(sum_prop.y + 10e-9),
              sum_ratio_cell_v2 = sum(sum_prop.x + 0.1)/sum(sum_prop.y + 0.1),
              sum_sb_cell_comp_v1 = sum_ratio_cell_v1*(ncells_f$tissue_fcells/ncells_m$tissue_mcells),
              sum_sb_cell_comp_v2 = sum_ratio_cell_v2*(ncells_f$tissue_fcells/ncells_m$tissue_mcells))
  
  sb_no_cell_comp <- tissue_bias_df %>%
    mutate(mean_prop.x = mean.x*prop.xy,
           mean_prop.y = mean.y*prop.xy,
           sum_prop.x = sum.x*prop.xy,
           sum_prop.y = sum.y*prop.xy) %>%
    group_by(gene) %>%
    summarise(mean_sb_no_cell_comp_v1 = sum(mean_prop.x + 10e-9)/sum(mean_prop.y + 10e-9),
              mean_sb_no_cell_comp_v2 = sum(mean_prop.x + 0.1)/sum(mean_prop.y + 0.1),
              sum_sb_no_cell_comp_v1 = sum(sum_prop.x + 10e-9)/sum(sum_prop.y + 10e-9),
              sum_sb_no_cell_comp_v2 = sum(sum_prop.x + 0.1)/sum(sum_prop.y + 0.1))
  
  rm(tissue_bias_df)
  gc()
  
  sb_df <- full_join(sb_cell_comp, sb_no_cell_comp)
  
  sb_df$tissue <- tolower(tissue)
  
  if (file.exists("sex_bias_cell_vs_no_cell_per_tissue_plus.csv")){
    write.table(sb_df,
                "sex_bias_cell_vs_no_cell_per_tissue_plus.csv",
                append = TRUE,
                sep = ",",
                col.names = FALSE,
                row.names = FALSE)
  } else {
    write.csv(sb_df,
              "sex_bias_cell_vs_no_cell_per_tissue_plus.csv",
              row.names = FALSE)
  }
  
  rm(ncells_f)
  rm(ncells_m)
  rm(sb_cell_comp)
  rm(sb_no_cell_comp)
  rm(sb_df)
  gc()
}


#### Get number of cells per cluster per sex ####
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
  metadata_df <- tissue_seurat_final@meta.data[, c('sex', 
                                                   's_broad_annotation')]
  
  metadata_summary <- metadata_df %>%
    group_by(sex, s_broad_annotation) %>%
    summarise(cells = n())
  
  
  metadata_summary$tissue <- tolower(tissue)
  
  if (file.exists("number_of_cells_per_cluster_all_tissues.csv")){
    write.table(metadata_summary,
                "number_of_cells_per_cluster_all_tissues.csv",
                append = TRUE,
                sep = ",",
                col.names = FALSE,
                row.names = FALSE)
  } else {
    write.csv(metadata_summary,
              "number_of_cells_per_cluster_all_tissues.csv",
              row.names = FALSE)
  }
  
  rm(tissue_seurat_final)
  rm(metadata_df)
  rm(metadata_summary)
  gc()
}
