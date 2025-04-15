##########################################
#     Removing sex-specific clusters     #
#     of cells from non-reproductive     #
#                 tissues                #
##########################################

#### Load libraries ####
library(reticulate)
use_python("/usr/bin/python3")
library(Seurat)
library(SeuratDisk)
library(loomR)

#### Set working directory ####
setwd("/media/cdecastr/Samsung USB/fly_cell_atlas_proj/data")

#### Loop over all tissue samples ####
tissues <- c('Antenna', 'Malpighian_tubule', 'Body_wall',
             'Oenocyte', 'Fat_body', 'Ovary', 'Gut',
             'Proboscis_and_maxillary_palps', 'Haltere',
             'Testis', 'Heart', 'Trachea', 'Leg',
             'Wing', 'Male_reproductive_glands', 'Head', 'Body')


#### Set clusters to remove ####
to_remove <- c("ovary cell", "testis", "ejaculatory bulb", "female reproductive system",
               "male accessory gland", "male germline differentiating cell",
               "btl-GAL4 positive female cell, cluster 1, likely to be ovary cell",
               "btl-GAL4 positive female cell, cluster 2, likely to be ovary cell",
               "btl-GAL4 positive female cell, cluster 3, likely to be ovary cell",
               "btl-GAL4 positive female cell, likely to be ovary cell, sim+",
               "btl-GAL4 positive female cell, likely to be ovary cell, sim+, H15+",
               "artefact")

##### If file format is loom #####
#### Convert from loom to h5Seurat ####
for (tissue in tissues) {
  # Connect to loom file
  tissue_loom <- Connect(filename = paste0('./r_fca_biohub_', tolower(tissue), '_10x.loom'),
                         mode = "r",
                         force = TRUE)
  # Gene list
  n_genes <- tissue_loom[['row_attrs']][['Gene']][['dims']]
  genes <- tissue_loom[['row_attrs']][['Gene']][1:n_genes]
  # Cell ID list
  n_cells <- tissue_loom[['col_attrs']][['CellID']][['dims']]
  cellids <- tissue_loom[['col_attrs']][['CellID']][1:n_cells]
  # Get raw counts matrix
  raw.counts <- as.data.frame(t(tissue_loom[['matrix']][,]))
  names(raw.counts) <- cellids
  rownames(raw.counts) <- genes
  rm(genes)
  gc()
  # Get metadata
  metadata <- data.frame(cellID = cellids,
                         ClusterID = tissue_loom[['col_attrs']][['ClusterID']][1:n_cells],
                         annotation = tissue_loom[['col_attrs']][['annotation']][1:n_cells],
                         broad_annotation = tissue_loom[['col_attrs']][['annotation_broad']][1:n_cells],
                         #batch = tissue_loom[['col_attrs']][['batch']][1:n_cells],
                         batch = tissue_loom[['col_attrs']][['batch_id']][1:n_cells],
                         sample_id = tissue_loom[['col_attrs']][['sample_id']][1:n_cells],
                         sex = tissue_loom[['col_attrs']][['sex']][1:n_cells],
                         tissue = tissue_loom[['col_attrs']][['tissue']][1:n_cells])
  tissue_loom$close_all()
  rm(tissue_loom)
  gc()
  rownames(metadata) <- cellids
  rm(cellids)
  gc()
  # Create Seurat object
  tissue_seurat <- CreateSeuratObject(counts = raw.counts,
                                      project = "fromLoom",
                                      assay = "RNA",
                                      meta.data = metadata)
  rm(metadata)
  rm(raw.counts)
  gc()
  # Remove sex-specific & artefact cells from non-sex-specific tissues
  if (tissue == 'Ovary' || tissue == 'Testis' || tissue == 'Male_reproductive_glands') {
    tissue_seurat_filtered <- subset(x = tissue_seurat,
                                     subset = annotation == "artefact",
                                     invert = TRUE)
  } else {
    tissue_seurat_filtered <- subset(x = tissue_seurat,
                                     subset = annotation %in% to_remove,
                                     invert = TRUE)
  }
  SaveH5Seurat(tissue_seurat_filtered,
               file = paste0('./r_raw_', tolower(tissue), '_filtered.h5Seurat'))
  rm(tissue_seurat)
  rm(tissue_seurat_filtered)
  gc()
}


