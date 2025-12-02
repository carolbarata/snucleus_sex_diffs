##############################################
#         Calculating X/A expression         #
#           ratios & plotting their          #
#        distribution across tissues         #
##############################################

#### Load libraries ####
library(ape)
library(anndata)
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(EnvStats)
library(gridExtra)
library(ggridges)
library(see)
library(ggpubr)
library(scales)

#### Loop over all tissue samples ####
tissues <- c('Antenna', 'Malpighian_tubule', 'Body_wall',
             'Oenocyte', 'Fat_body', 'Ovary', 'Gut',
             'Proboscis_and_maxillary_palps', 'Haltere',
             'Testis', 'Heart', 'Trachea', 'Leg',
             'Wing', 'Male_reproductive_glands', 'Head', 'Body')

gff_data <- read.gff(file = "./fly_cell_atlas_proj/data/dmel-genes-r6.31.gff")
gff_data$gene <- str_match(gff_data$attributes, 'Name=(.*?);Alias=')[1,2]
gff_data <- gff_data %>%
  mutate(gene = str_match(attributes, 'Name=(.*?);')[,2])


#### Calculate X/A expression ratio function ####
x_a_ratio_n_max <- function(cell) {
  #### Extract genes that have expression > 0 ####
  expressed_data <- as.data.frame(tissue_seurat@assays$SCT@counts[which(as.data.frame(tissue_seurat@assays$SCT@counts[,cell]) > 0), cell])
  expressed_genes <- data.frame(expression = expressed_data[,1],
                                gene = as.character(rownames(expressed_data)),
                                cell_id = as.character(colnames(tissue_seurat@assays$SCT[, cell])),
                                sex = as.character(tissue_seurat$sex[cell]),
                                cell_type = as.character(tissue_seurat$annotation[cell]))
  rm(expressed_data)
  
  expressed_genes <- merge(expressed_genes,
                           tissue_genes,
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
                                tissue = tolower(tissue),
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
                                tissue = tolower(tissue),
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


for (tissue in tissues) {
  ## Read in Seurat file ##
  tissue_seurat <- LoadH5Seurat(paste0('./r_', tolower(tissue),
                                       "_cell_type_filt_decontx_stringent_filt_sctransformed.h5Seurat"))
  ## Find genes in GFF data ##
  tissue_genes <- as.data.frame(subset(gff_data, gff_data$gene 
                                       %in%
                                         rownames(tissue_seurat@assays$SCT@counts))[,c('seqid', 'gene')])
  colnames(tissue_genes) <- c('chromosome', 'gene')
  ## Set chromosome status ##
  tissue_genes <- tissue_genes %>%
    mutate(status = case_when((chromosome == "X") ~ "X",
                              (chromosome == "Y") ~ "Y",
                              (chromosome == "mitochondrion_genome") ~ "mt",
                              TRUE ~ "autosome"))
  
  for (cell in 1:dim(tissue_seurat@assays$SCT@counts)[2]) {
   cell_data <- c()
   cell_data <- x_a_ratio_n_max(cell)
   if (class(cell_data) == "NULL") {
     next
   }
   
   ## Write output to file ##
   if (cell == 1){
     write.table(cell_data,
                 paste0("x_a_normalised_expression_ratio_n_max_expression_data_",
                        tolower(tissue),
                        "_relaxed_after_filtering.csv"),
                 row.names = FALSE,
                 col.names = TRUE)
   } else {
     write.table(cell_data,
                 paste0("x_a_normalised_expression_ratio_n_max_expression_data_",
                        tolower(tissue),
                        "_relaxed_after_filtering.csv"),
                 row.names = FALSE,
                 col.names = FALSE,
                 append = TRUE)
   }
  }

  rm(tissue_seurat)
  rm(tissue_genes)
  rm(cell_data)
  gc()
}

#### Read in data from csv ####
x_a_ratio_df <- list.files(path = "./",
                           pattern = "x_a_normalised_expression_ratio_n_max_expression_data_(.*)_relaxed_after_filtering[.]csv") %>%
  map_df(~ read.table(., header = TRUE))

x_a_ratio_df <- x_a_ratio_df %>%
  filter(sex != "mix")

x_a_ratio_df$tissue <- gsub("_", " ",
                            x_a_ratio_df$tissue)


#### Add stringent annotation to the data frame ####
## Create empty list for stringent annotation data ##
s_annotation_df <- vector("list", length(tissues))

for (i in 1:length(tissues)) {
  tissue <- tissues[i]
  ## Read in h5Seurat file ##
  tissue_seurat <- LoadH5Seurat(paste0('./r_', tolower(tissue),
                                       "_cell_type_filt_decontx_stringent_filt_sctransformed.h5Seurat"))
  ## Extract s_annotation info ##
  tissue_s_annotation <- as.data.frame(tissue_seurat@meta.data[,c('cellID',
                                                                  'tissue',
                                                                  'sex',
                                                                  's_annotation')])
  ## Add tissue s_annotation to dataframe variable ##
  s_annotation_df[[i]] <- as.data.frame(tissue_s_annotation)
  rm(tissue_seurat)
  rm(tissue_s_annotation)
  gc()
}

s_annotation_df <- dplyr::bind_rows(s_annotation_df)
colnames(s_annotation_df)[1] <- "cell_id"

s_annotation_df$tissue <- gsub("proboscis_and_maxpalp",
                               "proboscis_and_maxillary_palps",
                               s_annotation_df$tissue)

s_annotation_df$tissue <- gsub("testis", "gonad", s_annotation_df$tissue)
s_annotation_df$tissue <- gsub("ovary", "gonad", s_annotation_df$tissue)
s_annotation_df$tissue <- gsub("_", " ", s_annotation_df$tissue)

## Merge with X:A data frame ##
x_a_ratio_df$tissue <- gsub("testis", "gonad", x_a_ratio_df$tissue)
x_a_ratio_df$tissue <- gsub("ovary", "gonad", x_a_ratio_df$tissue)

x_a_ratio_df_plus <- left_join(x = x_a_ratio_df,
                               y = s_annotation_df)

## Manually replace gonadal pigment cell annotation ##
gonad_pigment_cells <- which(x_a_ratio_df_plus$s_annotation == 'pigment cell' &
                               x_a_ratio_df_plus$tissue == 'gonad')

x_a_ratio_df_plus[gonad_pigment_cells, 's_annotation'] <- 'male gonad pigment cell'

### Use broadly annotated cell types ###
annotation <- read.csv("./fly_cell_atlas_proj/data/cell_types_and_broad_annotation2.csv",
                       header = TRUE)

## Remove duplicates
annotation <- annotation[!duplicated(annotation),]

colnames(annotation) <- c("s_annotation", "broad_annotation")

## Merge with stringent annotation ##
x_a_ratio_df_plus <- left_join(x_a_ratio_df_plus,
                               annotation)

tissues_in_df <- unique(x_a_ratio_df_plus$tissue)
x_a_ratio_df_plus <- subset(x_a_ratio_df_plus, tissue %in% tissues_in_df[tissues_in_df != 'body'] |
                              tissue == 'body' &
                              broad_annotation != 'female reproductive system')

#### Save X:A table to file ####
write.csv(x_a_ratio_df_plus,
          file = "x_a_normalised_expression_ratio_n_max_expression_data_relaxed_after_filtering_all_tissues.csv",
          row.names = FALSE)
