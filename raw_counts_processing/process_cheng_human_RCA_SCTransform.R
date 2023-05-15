library(Seurat)
library(tidyverse)
library(data.table)
library(ggpubr)


################################################
# Take a look at human RCA data (Cheng et al)
cheng_rca_seurat = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Cheng_human_RCA_scRNA_data/RCA.rds")
cheng_rca_seurat
DefaultAssay(cheng_rca_seurat) = "SCT"
DimPlot(cheng_rca_seurat)


# Plot a few genes to get an idea of what cell types we're dealing with
genes_of_interest = c("LMOD1", "TNFRSF11B", "FN1", "LTBP1", "CD68", "IL1B", "APOE", "CDH5", "CD8A", "CSPG4", "CRTAC1", "FBLN1")
plot_list = gene_plot_list(cheng_rca_seurat, genes_of_interest = genes_of_interest)
ggarrange(plotlist = plot_list)


# Take a look at the metadata and sequencing depth
cheng_meta = cheng_rca_seurat@meta.data
table(cheng_meta$Pt)
cheng_meta %>%
  #group_by(Age) %>%
  #summarize(n=n())
  ggplot(aes(x=nCount_RNA, fill=Pt)) + 
  geom_histogram(bins = 50) +
  facet_wrap(~Pt, ncol = 1) +
  theme(aspect.ratio = 0.2) + 
  custom_theme

cheng_meta %>%
  group_by(Pt) %>%
  summarize(mean_seq_depth = mean(nCount_RNA))
  
###############################################################################
# Let's try to split this dataset into 6 individual libraries
cheng_rca_seurat_list = Seurat::SplitObject(cheng_rca_seurat, split.by = "Pt")

# We need to get the raw count matrices from each patient
extract_raw_counts = function(seurat_obj) { 
  return(seurat_obj@assays$RNA@counts)
  }

raw_matrices_list = lapply(cheng_rca_seurat_list, extract_raw_counts)
length(raw_matrices_list)
class(raw_matrices_list[[2]])

cheng_subset_library = raw_matrices_list[[3]]
cheng_processed_outs = doitall(cheng_subset_library, library_id = "Pt1", study_id = "Cheng_et_al", 
                               genes_of_interest = c("MYH11", "LMOD1", "LTBP1", "CRTAC1"))


########################################################################################
# Write a function to convert a Seurat obj to a df with countsfor cell types of interest 

# define target cells
target_clusters = c(1, 2)
target_barcodes = cheng_rca_seurat@meta.data %>%
  filter(orig.ident %in% target_clusters) %>%
  rownames()

# Subset df using base R
meta_df = cheng_rca_seurat@meta.data 
meta_df_subset = meta_df[meta_df$orig.ident %in% target_clusters, ]
target_barcodes = rownames(meta_df_subset)
names(target_barcodes) = meta_df_subset$orig.ident

# Get counts matrix
norm_counts = cheng_rca_seurat@assays$SCT@data
norm_matrix_counts_df = as.data.frame(t(as.matrix(norm_counts)))

# Subset counts matrix
counts_df_subset = norm_matrix_counts_df[target_barcodes, ]


# It'd also be neat to add a column matching the cell barcode to whatever annotation 
counts_df_subset$anno = names(target_barcodes)

# Now wrap everything into a function 
tidy_seurat_counts = function(seurat_obj, assay="SCT", 
                               target_coldata=NULL, 
                               target_anno=NULL) {
  
  # Get counts matrix
  if (assay=="SCT") {
    norm_counts = seurat_obj@assays$SCT@data
  } else if (assay=="RNA") {
    norm_counts = seurat_obj@assays$RNA@data
  } else {
    msg = ("The provided assay should be either SCT or RNA!")
    stop(msg)
  }
  
  norm_matrix_counts_df = as.data.frame(t(as.matrix(norm_counts)))
  
  # What if we want to extract the entire counts matrix
  if (!is.null(target_anno) & !is.null(target_coldata)) {
    meta_df = seurat_obj@meta.data
    
    # Subset df
    meta_df_subset = meta_df[meta_df[[target_coldata]] %in% target_anno, ]
    target_barcodes = rownames(meta_df_subset)
    names(target_barcodes) = meta_df_subset[[target_coldata]]
    
    # Subset counts df to keep only cell types of interest
    counts_df_subset = norm_matrix_counts_df[target_barcodes,]
    
    # Add a column matching the cell barcode to its corresponding annotation 
    counts_df_subset$annotation = names(target_barcodes)
    return(counts_df_subset)
  } else {
    return(norm_matrix_counts_df)
  }
} 



