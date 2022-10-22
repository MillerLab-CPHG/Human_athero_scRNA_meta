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
# Read sample 1: GSM4837523 Patient1 Carotid Artery Atherosclerotic Core (AC) #
###############################################################################

carotid_ac_p1 = Seurat::Read10X("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Alsaigh_human_carotids_scRNA_data/GSM4837523_02dat20190515tisCARconDIS_featurebcmatrixfiltered/filtered_feature_bc_matrix/")

# Raw, unfiltered matrix is 33538 genes x 11015 cells
# Why so many genes? 
dim(carotid_ac_p1)

###############################################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 17023 genes x 10982 cells
ac_p1_seurat_sct = CreateSeuratObject(counts = carotid_ac_p1, 
                                  project = "alsaigh_p1_carotid_AC", 
                                  min.cells = 10, 
                                  min.features = 200)

##################################
# STEP 3
# Detect doublets with scDblFinder
# Called 739 doublets
ac_p1_doublet_ids = scDblFinder_clusters(ac_p1_seurat_sct, nrep = 3)
length(ac_p1_doublet_ids)

# First filtering step: Use the vector of consensus doublet IDs output below
ac_p1_singlets = setdiff(colnames(carotid_ac_p1), ac_p1_doublet_ids)
length(ac_p1_singlets)

# New dims after removing doublets: 16805 genes x 10243 cells
ac_p1_seurat_sct = CreateSeuratObject(counts = carotid_ac_p1[, ac_p1_singlets], 
                                      project = "alsaigh_p1_carotid_AC", 
                                      min.cells = 10, 
                                      min.features = 200)

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
ac_p1_seurat_sct =  decontX_remove(ac_p1_seurat_sct)
head(ac_p1_seurat_sct@assays$RNA@counts)

#################################
# STEP 2 and 4 before/after removing doublets
# SCT normalize
ac_p1_seurat_sct = Seurat_SCT_process(ac_p1_seurat_sct, seurat_filter = TRUE,
                                      sample_id = "alsaigh_ac_p1", 
                                      study_name = "alsaigh_et_al")

# Finding the right resolution is quite important for doublet detection and removal
ac_p1_seurat_sct =  FindClusters(ac_p1_seurat_sct, resolution = 1)   

# Visualize clusters
p2_before_QC = DimPlot(ac_p1_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

p2_after_QC = DimPlot(ac_p1_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")

FeaturePlot(ac_p1_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                       "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)

FeaturePlot(ac_p1_seurat_sct, features = c("CD68", "CD14", "CD163", "APOE", 
                                       "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)


###################
# Save .rds object
# Number of cells after processing: 9655 cells
saveRDS(ac_p1_seurat_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_ac_p1_processed_seurat_obj.rds")


###############################################################################
# Read sample 2: GSM4837524 Patient1 Carotid Artery Proximal Adjacent (PA)    #
###############################################################################

carotid_pa_p1 = Seurat::Read10X("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Alsaigh_human_carotids_scRNA_data/GSM4837524_01dat20190515tisCARconHEA_featurebcmatrixfiltered/filtered_feature_bc_matrix/")

# Raw, unfiltered matrix is 33538 genes x 3716 cells
# Why so many genes? 
dim(carotid_pa_p1)

################################################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 15693 genes x 3629 cells
pa_p1_seurat_sct = CreateSeuratObject(counts = carotid_pa_p1, 
                                      project = "alsaigh_p1_carotid_PA", 
                                      min.cells = 10, 
                                      min.features = 200)

##################################
# STEP 3
# Detect doublets with scDblFinder
# Called 72 doublets
pa_p1_doublet_ids = scDblFinder_clusters(pa_p1_seurat_sct, nrep = 3)
length(pa_p1_doublet_ids)

# First filtering step: Use the vector of consensus doublet IDs output below
pa_p1_singlets = setdiff(colnames(carotid_pa_p1), pa_p1_doublet_ids)
length(pa_p1_singlets)

# New dims after removing doublets: 15645 genes x 3557 samples
pa_p1_seurat_sct = CreateSeuratObject(counts = carotid_pa_p1[, pa_p1_singlets], 
                                      project = "alsaigh_p1_carotid_PA", 
                                      min.cells = 10, 
                                      min.features = 200)

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
pa_p1_seurat_sct =  decontX_remove(pa_p1_seurat_sct)
head(pa_p1_seurat_sct@assays$RNA@counts)

###################################
# STEP 2 and 4 before/after removing doublets
# SCT normalize
pa_p1_seurat_sct = Seurat_SCT_process(pa_p1_seurat_sct, seurat_filter = TRUE,
                                      sample_id = "alsaigh_pa_p1", 
                                      study_name = "alsaigh_et_al")

# Finding the right resolution is quite important for doublet detection and removal
pa_p1_seurat_sct =  FindClusters(pa_p1_seurat_sct, resolution = 0.5)   

# Visualize clusters
p2_before_QC = DimPlot(pa_p1_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

p2_after_QC = DimPlot(pa_p1_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")

FeaturePlot(pa_p1_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                           "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)

FeaturePlot(pa_p1_seurat_sct, features = c("CD68", "CD14", "CD163", "APOE", 
                                           "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)

#####################################
# Save .rds object
# Number of cells after processing: 3056
saveRDS(pa_p1_seurat_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_pa_p1_processed_seurat_obj.rds")


# Calculate silhouette width
dist.matrix = dist(x=Embeddings(object = pa_p1_seurat_sct[["pca"]])[, 1:30])
clusters = pa_p1_seurat_sct$seurat_clusters
sil = silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
sil_scores = as.data.frame(sil[, 1:3])
sil_scores %>%
  ggplot(aes(x=sil_width)) + 
  geom_density() +
  facet_wrap(~cluster)


# Sil scores before QC
sil_scores_before_QC = sil_scores %>%
  group_by(cluster) %>%
  summarize(mean_sil_width = mean(sil_width))

sil_scores_beforeQC_mean = mean(sil_scores$sil_width)

# Sil width scores after QC
silQC = sil_scores
sil_scores_QC = sil_scores %>%
  group_by(cluster) %>%
  summarize(mean_sil_width = mean(sil_width))

sil_scoresQC_mean = mean(sil_scores$sil_width)
sil_scoresQC_mean

###############################################################################
# Read sample 3: GSM4837525 Patient2 Carotid Artery Atherosclerotic Core (AC) #
###############################################################################

carotid_ac_p2 = Seurat::Read10X("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Alsaigh_human_carotids_scRNA_data/GSM4837525_02dat20190620tisCARconDIS_featurebcmatrixfiltered/filtered_feature_bc_matrix/")

# Raw, unfiltered matrix is 33538 genes x 15960 cells
dim(carotid_ac_p2)

#################################################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 17901 genes x 15923 cells
ac_p2_seurat_sct = CreateSeuratObject(counts = carotid_ac_p2, 
                                      project = "alsaigh_p2_carotid_AC", 
                                      min.cells = 10, 
                                      min.features = 200)


##################################
# STEP 3
# Detect doublets with scDblFinder
# Called 1693 doublets
ac_p2_doublet_ids = scDblFinder_clusters(ac_p2_seurat_sct, nrep = 3)
length(ac_p2_doublet_ids)

# First filtering step: Use the vector of consensus doublet IDs output below
# Dims of Seurat object: 16840 genes x 10494 cells
ac_p2_singlets = setdiff(colnames(carotid_ac_p2), ac_p2_doublet_ids)
length(ac_p2_singlets)

# New dims after removing doublets: 17585 genes x 14230 samples
ac_p2_seurat_sct = CreateSeuratObject(counts = carotid_ac_p2[, ac_p2_singlets], 
                                      project = "alsaigh_p2_carotid_AC", 
                                      min.cells = 10, 
                                      min.features = 200)

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
ac_p2_seurat_sct =  decontX_remove(ac_p2_seurat_sct)
head(ac_p2_seurat_sct@assays$RNA@counts)



##################################
# STEP 2 and 4 before/after removing doublets
# SCT normalize
ac_p2_seurat_sct = Seurat_SCT_process(ac_p2_seurat_sct, seurat_filter = TRUE,
                                      sample_id = "alsaigh_ac_p2", 
                                      study_name = "alsaigh_et_al")

# Finding the right resolution is quite important for doublet detection and removal
ac_p2_seurat_sct =  FindClusters(ac_p2_seurat_sct, resolution = 1)   

# Visualize clusters
p3_before_QC = DimPlot(ac_p2_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

p3_after_QC = DimPlot(ac_p2_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")

FeaturePlot(ac_p2_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                           "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)

FeaturePlot(ac_p2_seurat_sct, features = c("CD68", "CD14", "CD163", "APOE", 
                                           "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)


##########################################
# Save .rds object
# Nunmber of cells after processing: 12214
saveRDS(ac_p2_seurat_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_ac_p2_processed_seurat_obj.rds")


###############################################################################
# Read sample 4: GSM4837526 Patient2 Carotid Artery Proximal Adjacent (PA)    #
###############################################################################

carotid_pa_p2 = Seurat::Read10X("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Alsaigh_human_carotids_scRNA_data/GSM4837526_01dat20190620tisCARconHEA_featurebcmatrixfiltered//filtered_feature_bc_matrix/")

# Raw, unfiltered matrix is 33538 genes x 5523 cells
dim(carotid_pa_p2)

################################################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 16458 genes x 5504 cells
pa_p2_seurat_sct = CreateSeuratObject(counts = carotid_pa_p2, 
                                      project = "alsaigh_p2_carotid_PA", 
                                      min.cells = 10, 
                                      min.features = 200)

##################################
# STEP 3
# Detect doublets with scDblFinder
# Called 172 doublets
pa_p2_doublet_ids = scDblFinder_clusters(pa_p2_seurat_sct, nrep = 3)
length(pa_p2_doublet_ids)

# First filtering step: Use the vector of consensus doublet IDs output below
pa_p2_singlets = setdiff(colnames(carotid_pa_p2), pa_p2_doublet_ids)
length(pa_p2_singlets)

# New dims after removing doublets: 16347 genes x 5332 cells
pa_p2_seurat_sct = CreateSeuratObject(counts = carotid_pa_p2[, pa_p2_singlets], 
                                      project = "alsaigh_p2_carotid_PA", 
                                      min.cells = 10, 
                                      min.features = 200)

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
pa_p2_seurat_sct =  decontX_remove(pa_p2_seurat_sct)
head(pa_p2_seurat_sct@assays$RNA@counts)

###################################
# STEP 2 and 4 before/after removing doublets
# SCT normalize
pa_p2_seurat_sct = Seurat_SCT_process(pa_p2_seurat_sct, seurat_filter = TRUE,
                                      sample_id = "alsaigh_pa_p2", 
                                      study_name = "alsaigh_et_al")

# Finding the right resolution is quite important for doublet detection and removal
pa_p2_seurat_sct =  FindClusters(pa_p2_seurat_sct, resolution = 0.5)   

# Visualize clusters
p4_before_QC = DimPlot(pa_p2_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

p4_after_QC = DimPlot(pa_p2_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")

FeaturePlot(pa_p2_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                           "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)

FeaturePlot(pa_p2_seurat_sct, features = c("CD68", "CD14", "CD163", "APOE", 
                                           "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)

# Save .rds object
# Number of cells after processing: 4630 cells 
saveRDS(pa_p2_seurat_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_pa_p2_processed_seurat_obj.rds")



###############################################################################
# Read sample 5: GSM4837527 Patient3 Carotid Artery Atherosclerotic Core (AC) #
###############################################################################

carotid_ac_p3 = Seurat::Read10X("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Alsaigh_human_carotids_scRNA_data/GSM4837527_02dat20190717tisCARconDIS_featurebcmatrixfiltered/filtered_feature_bc_matrix/")

# Raw, unfiltered matrix is 33538 genes x 12388 cells
dim(carotid_ac_p3)

################################################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 17876 genes x 12339 cells
ac_p3_seurat_sct = CreateSeuratObject(counts = carotid_ac_p3, 
                                      project = "alsaigh_p3_carotid_AC", 
                                      min.cells = 10, 
                                      min.features = 200)

##################################
# STEP 3
# Detect doublets with scDblFinder
ac_p3_doublet_ids = scDblFinder_clusters(ac_p3_seurat_sct, nrep = 3)
length(ac_p3_doublet_ids)

# First filtering step: Use the vector of consensus doublet IDs output below
ac_p3_singlets = setdiff(colnames(carotid_ac_p3), ac_p3_doublet_ids)
length(ac_p3_singlets)

# New dims after removing doublets: 17671 genes x 11447 cells
ac_p3_seurat_sct = CreateSeuratObject(counts = carotid_ac_p3[, ac_p3_singlets], 
                                      project = "alsaigh_p3_carotid_AC", 
                                      min.cells = 10, 
                                      min.features = 200)

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
ac_p3_seurat_sct =  decontX_remove(ac_p3_seurat_sct)
head(ac_p3_seurat_sct@assays$RNA@counts)

###################################
# STEP 2 and 4 before/after removing doublets
# SCT normalize
ac_p3_seurat_sct = Seurat_SCT_process(ac_p3_seurat_sct, seurat_filter = TRUE,
                                      sample_id = "alsaigh_ac_p3", 
                                      study_name = "alsaigh_et_al")

# Finding the right resolution is quite important for doublet detection and removal
ac_p3_seurat_sct =  FindClusters(ac_p3_seurat_sct, resolution = 1)   

# Visualize clusters
p5_before_QC = DimPlot(ac_p3_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

p5_after_QC = DimPlot(ac_p3_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")


FeaturePlot(ac_p3_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                           "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)

FeaturePlot(ac_p3_seurat_sct, features = c("CD68", "CD14", "CD163", "TREM2", 
                                           "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)

FeaturePlot(ac_p3_seurat_sct, features = c("MYH11", "CNN1", "FHL5", "FOXC1"), order = TRUE,  pt.size = 0.1)

# Save .rds objects
# Number of cells after processing: 9761
saveRDS(ac_p3_seurat_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_ac_p3_processed_seurat_obj.rds")


###############################################################################
# Read sample 6: GSM4837528 Patient3 Carotid Artery Proximal Adjacent (PA) #
###############################################################################

carotid_pa_p3 = Seurat::Read10X("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Alsaigh_human_carotids_scRNA_data/GSM4837528_01dat20190717tisCARconHEA_featurebcmatrixfiltered/filtered_feature_bc_matrix/")

# Raw, unfiltered matrix is 33538 genes x 3379 cells
dim(carotid_pa_p3)

#################################################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 16133 genes x 3344 cells
pa_p3_seurat_sct = CreateSeuratObject(counts = carotid_pa_p3, 
                                      project = "alsaigh_p3_carotid_PA", 
                                      min.cells = 10, 
                                      min.features = 200)

##################################
# STEP 3
# Detect doublets with scDblFinder
pa_p3_doublet_ids = scDblFinder_clusters(pa_p3_seurat_sct, nrep = 3)
length(pa_p3_doublet_ids)

# First filtering step: Use the vector of consensus doublet IDs output below
pa_p3_singlets = setdiff(colnames(carotid_pa_p3), pa_p3_doublet_ids)
length(pa_p3_singlets)

# New dims after removing doublets: 16073 genes x 3285 cells
pa_p3_seurat_sct = CreateSeuratObject(counts = carotid_pa_p3[, pa_p3_singlets], 
                                      project = "alsaigh_p3_carotid_AC", 
                                      min.cells = 10, 
                                      min.features = 200)

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
pa_p3_seurat_sct =  decontX_remove(pa_p3_seurat_sct)
head(pa_p3_seurat_sct@assays$RNA@counts)

###################################
# STEP 2 and 4 before/after removing doublets
# SCT normalize
pa_p3_seurat_sct = Seurat_SCT_process(pa_p3_seurat_sct, seurat_filter = TRUE,
                                      sample_id = "alsaigh_pa_p3", 
                                      study_name = "alsaigh_et_al")

# Finding the right resolution is quite important for doublet detection and removal
pa_p3_seurat_sct =  FindClusters(pa_p3_seurat_sct, resolution = 0.5)   

# Visualize clusters
p6_before_QC = DimPlot(pa_p3_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

p6_after_QC = DimPlot(pa_p3_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")

FeaturePlot(pa_p3_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                           "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)

FeaturePlot(pa_p3_seurat_sct, features = c("CD68", "CD14", "CD163", "APOE", 
                                           "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)

# Save .rds object
# Number of cells after processing: 2791
saveRDS(pa_p3_seurat_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_pa_p3_processed_seurat_obj.rds")





