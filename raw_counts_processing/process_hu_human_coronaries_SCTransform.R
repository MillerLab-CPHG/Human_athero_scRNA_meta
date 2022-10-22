library(Seurat)
library(tidyverse)
library(scDblFinder)
library(celda)
library(sctransform)
library(cluster)


# Set seed for doublet removal reproducibility
set.seed(1)

# Load scRNA analysis utils
source("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/R_scripts/Utils/scRNA_processing_utils.R")

# Run some tests to develop a QC and processing pipeline for the meta-analysis

###############################################################################
# Read sample 1: Patient1CA Coronary Artery (non diseased artery)             #
###############################################################################

coronary_1_p1 = Seurat::Read10X("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Hu_et_al_human_coronaries_aortas_pulmonary/cardiac_arteries_processed_data/patient1_CA/")

# Raw, unfiltered matrix is 33538 genes x 23512 cells
# Why so many genes? 
dim(coronary_1_p1)

###############################################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 19094 genes x 23354 cells
coronary1_p1_seurat_sct = CreateSeuratObject(counts = coronary1_p1, 
                                  project = "hu_p1_coronary1", 
                                  min.cells = 10, 
                                  min.features = 200)

##################################
# STEP 3
# Detect doublets with scDblFinder
# Called 3694 doublets
coronary1_p1_doublet_ids = scDblFinder_clusters(coronary1_p1_seurat_sct, nrep = 3)
length(coronary1_p1_doublet_ids)

# First filtering step: Use the vector of consensus doublet IDs output below
coronary1_p1_singlets = setdiff(colnames(coronary_1_p1), coronary1_p1_doublet_ids)
length(coronary1_p1_singlets)

# New dims after removing doublets: 18509 genes x 19660 cells
coronary1_p1_seurat_sct = CreateSeuratObject(counts = coronary_1_p1[, coronary1_p1_singlets], 
                                      project = "hu_p1_coronary1", 
                                      min.cells = 10, 
                                      min.features = 200)

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
coronary1_p1_seurat_sct =  decontX_remove(coronary1_p1_seurat_sct)
head(coronary1_p1_seurat_sct@assays$RNA@counts)

#################################
# STEP 2 and 4 before/after removing doublets
# SCT normalize
coronary1_p1_seurat_sct = Seurat_SCT_process(coronary1_p1_seurat_sct, seurat_filter = TRUE,
                                             sample_id = "hu_coronary1_p1", 
                                             study_name = "hu_et_al")

# Finding the right resolution is quite important for doublet detection and removal
coronary1_p1_seurat_sct =  FindClusters(coronary1_p1_seurat_sct, resolution = 0.9)   

# Visualize clusters
p1_before_QC = DimPlot(coronary1_p1_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

p2_after_QC = DimPlot(coronary1_p1_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")

FeaturePlot(coronary1_p1_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                       "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)

FeaturePlot(coronary1_p1_seurat_sct, features = c("CD68", "CD14", "CD163", "APOE", 
                                       "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)


###################
# Save .rds object
# Number of cells after processing: 17494 cells
saveRDS(coronary1_p1_seurat_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/hu_et_al_rds_objects/hu_coronary1_p1_processed_seurat_obj.rds")
coronary1_p1_seurat_sct = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/hu_et_al_rds_objects/hu_coronary1_p1_processed_seurat_obj.rds")


###############################################################################
# Read sample 2: Patient2 Coronary Artery 1                                   #
###############################################################################

coronary1_p2 = Seurat::Read10X("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Hu_et_al_human_coronaries_aortas_pulmonary/cardiac_arteries_processed_data/patient2_CA1/")

# Raw, unfiltered matrix is 33538 genes x 15349 cells
# Why so many genes? 
dim(coronary1_p2)

################################################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 18553 genes x 15283 cells
coronary1_p2_seurat_sct = CreateSeuratObject(counts = coronary1_p2, 
                                      project = "hu_p2_coronary1", 
                                      min.cells = 10, 
                                      min.features = 200)

##################################
# STEP 3
# Detect doublets with scDblFinder
# Called 1538 doublets
coronary1_p2_doublet_ids = scDblFinder_clusters(coronary1_p2_seurat_sct, nrep = 3)
length(coronary1_p2_doublet_ids)

# First filtering step: Use the vector of consensus doublet IDs output below
coronary1_p2_singlets = setdiff(colnames(coronary1_p2), coronary1_p2_doublet_ids)
length(coronary1_p2_singlets)

# New dims after removing doublets: 18188 genes x 13745 samples
coronary1_p2_seurat_sct = CreateSeuratObject(counts = coronary1_p2[, coronary1_p2_singlets], 
                                      project = "hu_p2_coronary1", 
                                      min.cells = 10, 
                                      min.features = 200)

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
coronary1_p2_seurat_sct =  decontX_remove(coronary1_p2_seurat_sct)
head(coronary1_p2_seurat_sct@assays$RNA@counts)

###################################
# STEP 2 and 4 before/after removing doublets
# SCT normalize
coronary1_p2_seurat_sct = Seurat_SCT_process(coronary1_p2_seurat_sct, seurat_filter = TRUE,
                                      sample_id = "hu_coronary1_p2", 
                                      study_name = "hu_et_al")

# Finding the right resolution is quite important for doublet detection and removal
coronary1_p2_seurat_sct =  FindClusters(coronary1_p2_seurat_sct, resolution = 0.9)   

# Visualize clusters
p2_before_QC = DimPlot(coronary1_p2_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

p2_after_QC = DimPlot(coronary1_p2_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")

FeaturePlot(coronary1_p2_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                           "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)

FeaturePlot(coronary1_p2_seurat_sct, features = c("CD68", "CD14", "CD163", "APOE", 
                                           "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)

#####################################
# Save .rds object
# Number of cells after processing: 13299 cells
saveRDS(coronary1_p2_seurat_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/hu_et_al_rds_objects/hu_coronary1_p2_processed_seurat_obj.rds")
coronary1_p2_seurat_sct = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/hu_et_al_rds_objects/hu_coronary1_p2_processed_seurat_obj.rds")


###############################################################################
# Read sample 3: Patient2 Coronary Artery 2                                   #
###############################################################################

coronary2_p2 = Seurat::Read10X("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Hu_et_al_human_coronaries_aortas_pulmonary/cardiac_arteries_processed_data/patient2_CA2/")

# Raw, unfiltered matrix is 33538 genes x 14973 cells
dim(coronary2_p2)

#################################################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 18538 genes x 14916 cells
coronary2_p2_seurat_sct = CreateSeuratObject(counts = coronary2_p2, 
                                      project = "hu_p2_coronary2", 
                                      min.cells = 10, 
                                      min.features = 200)


##################################
# STEP 3
# Detect doublets with scDblFinder
# Called 1396 doublets
coronary2_p2_doublet_ids = scDblFinder_clusters(coronary2_p2_seurat_sct, nrep = 3)
length(coronary2_p2_doublet_ids)

# First filtering step: Use the vector of consensus doublet IDs output below
coronary2_p2_singlets = setdiff(colnames(coronary2_p2), coronary2_p2_doublet_ids)
length(coronary2_p2_singlets)

# New dims after removing doublets: 18220 genes x 13520 cells
coronary2_p2_seurat_sct = CreateSeuratObject(counts = coronary2_p2[, coronary2_p2_singlets], 
                                      project = "hu_p2_coronary2", 
                                      min.cells = 10, 
                                      min.features = 200)

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
coronary2_p2_seurat_sct =  decontX_remove(coronary2_p2_seurat_sct)
head(coronary2_p2_seurat_sct@assays$RNA@counts)



##################################
# STEP 2 and 4 before/after removing doublets
# SCT normalize
coronary2_p2_seurat_sct = Seurat_SCT_process(coronary2_p2_seurat_sct, seurat_filter = FALSE,
                                      sample_id = "hu_coronary2_p2", 
                                      study_name = "hu_et_al")

# Finding the right resolution is quite important for doublet detection and removal
coronary2_p2_seurat_sct =  FindClusters(coronary2_p2_seurat_sct, resolution = 0.8)   

# Visualize clusters
p3_before_QC = DimPlot(coronary2_p2_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets") + custom_theme + level2_annotations_scale
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure1/SuppFig1b_library_before_QC.svg",
       plot = p3_before_QC, width = 7, height = 7)

p3_after_QC = DimPlot(coronary2_p2_seurat_sct_new, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX") + custom_theme + level2_annotations_scale

FeaturePlot(coronary2_p2_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                           "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)

FeaturePlot(coronary2_p2_seurat_sct, features = c("CD68", "CD14", "CD163", "APOE", 
                                           "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure1/SuppFig1b_library_after_QC.svg",
       plot = p3_after_QC, width = 7, height = 7)


##########################################
# Save .rds object
# Nunmber of cells after processing: 13167
#saveRDS(coronary2_p2_seurat_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/hu_et_al_rds_objects/hu_coronary2_p2_processed_seurat_obj.rds")
coronary2_p2_seurat_sct_new = readRDS("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/hu_et_al_rds_objects/hu_coronary2_p2_processed_seurat_obj.rds")

###############################################################################
# Read sample 4: GSM4837526 Patient3 Coronary Artery 1    #
###############################################################################

coronary1_p3 = Seurat::Read10X("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Hu_et_al_human_coronaries_aortas_pulmonary/cardiac_arteries_processed_data/patient3_CA1/")

# Raw, unfiltered matrix is 33538 genes x 10123 cells
dim(coronary1_p3)

################################################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 18340 genes x 10098 cells
coronary1_p3_seurat_sct = CreateSeuratObject(counts = coronary1_p3, 
                                      project = "hu_p3_coronary1", 
                                      min.cells = 10, 
                                      min.features = 200)

##################################
# STEP 3
# Detect doublets with scDblFinder
# Called 614 doublets
coronary1_p3_doublet_ids = scDblFinder_clusters(coronary1_p3_seurat_sct, nrep = 3)
length(coronary1_p3_doublet_ids)

# First filtering step: Use the vector of consensus doublet IDs output below
coronary1_p3_singlets = setdiff(colnames(coronary1_p3), coronary1_p3_doublet_ids)
length(coronary1_p3_singlets)

# New dims after removing doublets: 18149 genes x 9484 cells
coronary1_p3_seurat_sct = CreateSeuratObject(counts = coronary1_p3[, coronary1_p3_singlets], 
                                      project = "hu_p3_coronary1", 
                                      min.cells = 10, 
                                      min.features = 200)

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
coronary1_p3_seurat_sct =  decontX_remove(coronary1_p3_seurat_sct)
head(coronary1_p3_seurat_sct@assays$RNA@counts)

###################################
# STEP 2 and 4 before/after removing doublets
# SCT normalize
coronary1_p3_seurat_sct = Seurat_SCT_process(coronary1_p3_seurat_sct, seurat_filter = TRUE,
                                      sample_id = "hu_coronary1_p3", 
                                      study_name = "hu_et_al")

# Finding the right resolution is quite important for doublet detection and removal
coronary1_p3_seurat_sct =  FindClusters(coronary1_p3_seurat_sct, resolution = 0.8)   

# Visualize clusters
p4_before_QC = DimPlot(coronary1_p3_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

p4_after_QC = DimPlot(coronary1_p3_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")

FeaturePlot(coronary1_p3_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                           "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)

FeaturePlot(coronary1_p3_seurat_sct, features = c("CD68", "CD14", "CD163", "APOE", 
                                           "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)

# Save .rds object
# Number of cells after processing: 8600 cells 
saveRDS(coronary1_p3_seurat_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/hu_et_al_rds_objects/hu_coronary1_p3_processed_seurat_obj.rds")



###############################################################################
# Read sample 5: GSM4837527 Patient3 Coronary Artery 2 #
###############################################################################

coronary2_p3 = Seurat::Read10X("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Hu_et_al_human_coronaries_aortas_pulmonary/cardiac_arteries_processed_data/patient3_CA2/")

# Raw, unfiltered matrix is 33538 genes x 8091 cells
dim(coronary2_p3)

################################################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 17861 genes x 7994 cells
coronary2_p3_seurat_sct = CreateSeuratObject(counts = coronary2_p3, 
                                      project = "hu_p3_coronary2", 
                                      min.cells = 10, 
                                      min.features = 200)

##################################
# STEP 3
# Detect doublets with scDblFinder
# Called 365 doublets
coronary2_p3_doublet_ids = scDblFinder_clusters(coronary2_p3_seurat_sct, nrep = 3)
length(coronary2_p3_doublet_ids)

# First filtering step: Use the vector of consensus doublet IDs output below
coronary2_p3_singlets = setdiff(colnames(coronary2_p3), coronary2_p3_doublet_ids)
length(coronary2_p3_singlets)

# New dims after removing doublets: 17716 genes x 7629 cells
coronary2_p3_seurat_sct = CreateSeuratObject(counts = coronary2_p3[, coronary2_p3_singlets], 
                                      project = "hu_p3_coronary2", 
                                      min.cells = 10, 
                                      min.features = 200)

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
coronary2_p3_seurat_sct =  decontX_remove(coronary2_p3_seurat_sct)
head(coronary2_p3_seurat_sct@assays$RNA@counts)

###################################
# STEP 2 and 4 before/after removing doublets
# SCT normalize
coronary2_p3_seurat_sct = Seurat_SCT_process(coronary2_p3_seurat_sct, seurat_filter = TRUE,
                                      sample_id = "hu_coronary2_p3", 
                                      study_name = "hu_et_al")

# Finding the right resolution is quite important for doublet detection and removal
coronary2_p3_seurat_sct =  FindClusters(coronary2_p3_seurat_sct, resolution = 0.8)   

# Visualize clusters
p5_before_QC = DimPlot(coronary2_p3_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

p5_after_QC = DimPlot(coronary2_p3_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")


FeaturePlot(coronary2_p3_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                           "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)

FeaturePlot(coronary2_p3_seurat_sct, features = c("CD68", "CD14", "CD163", "TREM2", 
                                           "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)


# Save .rds objects
# Number of cells after processing: 6933
saveRDS(coronary2_p3_seurat_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/hu_et_al_rds_objects/hu_coronary2_p3_processed_seurat_obj.rds")









