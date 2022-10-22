library(Seurat)
library(tidyverse)
library(scDblFinder)
library(celda)
library(sctransform)
library(cluster)


# Set seed for doublet removal reproducibility
set.seed(1)


# PROCESS HUMAN CAROTID MATRICES USING SCTransform

#######################################################################
# Read Sample 1: GSM4705589 RPE004 (Caucasian, Male, 83, Symptomatic) #
#######################################################################

# Load Matrix
rpe004_matrix = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Pan_et_al_human_carotids_scRNA_data/GSE155512_RAW/GSM4705589_RPE004_matrix.txt", sep="\t", header = TRUE)

rownames(rpe004_matrix) = rpe004_matrix$gene
head(rownames(rpe004_matrix))
dim(rpe004_matrix)

# Make sparse matrix and remove "gene" column
rpe004_sparse_mtx = as.sparse(rpe004_matrix)
rpe004_sparse_mtx = rpe004_sparse_mtx[, -1]

# 15796 genes x 2614 cells
dim(rpe004_sparse_mtx)

##################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 15179 genes x 2614 cells
rpe004_seurat_sct = CreateSeuratObject(counts = rpe004_sparse_mtx, 
                                      project = "pan_rpe004_carotid", 
                                      min.cells = 10, 
                                      min.features = 200)

##################################
# STEP 3
# Detect doublets with scDblFinder
# Called 8 doublets
rpe004_doublet_ids = scDblFinder_clusters(rpe004_seurat_sct, nrep = 3)
length(rpe004_doublet_ids)

# First filtering step: Use the vector of consensus doublet IDs output below
rpe004_singlets = setdiff(colnames(rpe004_sparse_mtx), rpe004_doublet_ids)
length(rpe004_singlets)

# New dims after removing doublets: 15174 genes x 2606 cells
rpe004_seurat_sct = CreateSeuratObject(counts = rpe004_sparse_mtx[, rpe004_singlets], 
                                      project = "pan_rpe004", 
                                      min.cells = 10, 
                                      min.features = 200)

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
rpe004_seurat_sct =  decontX_remove(rpe004_seurat_sct)
head(rpe004_seurat_sct@assays$RNA@counts)

###################################
# STEP 2 and 4 before/after removing doublets
# SCT normalize
rpe004_seurat_sct = Seurat_SCT_process(rpe004_seurat_sct, seurat_filter = TRUE,
                                      sample_id = "pan_rpe004", 
                                      study_name = "pan_et_al")

# Finding the right resolution is quite important for doublet detection and removal
rpe004_seurat_sct =  FindClusters(rpe004_seurat_sct, resolution = 0.5)   

# Visualize clusters
p1_before_QC = DimPlot(rpe004_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

p1_after_QC = DimPlot(rpe004_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")

FeaturePlot(rpe004_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                           "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)

FeaturePlot(rpe004_seurat_sct, features = c("CD68", "CD14", "CD163", "APOE", 
                                           "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)

FeaturePlot(rpe004_seurat_sct, features = c("MYH11", "CNN1", "FHL5", "FOXC1"), order = TRUE,  pt.size = 0.1)


########################################
# Save RDS object
# NUmber of cells after processing: 2606
saveRDS(rpe004_seurat_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/pan_et_al_rds_objects/pan_rpe004_processed_seurat_obj.rds")


#########################################################################
# Read Sample 2: GSM4705590 RPE005 (Caucasian, Male, 67, Asymptomatic)  #
#########################################################################

# Load Matrix
rpe005_matrix = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Pan_et_al_human_carotids_scRNA_data/GSE155512_RAW/GSM4705590_RPE005_matrix.txt", sep="\t", header = TRUE)

rownames(rpe005_matrix) = rpe005_matrix$gene
head(rownames(rpe005_matrix))
dim(rpe005_matrix)

# Make sparse matrix and remove "gene" column
rpe005_sparse_mtx = as.sparse(rpe005_matrix)
rpe005_sparse_mtx = rpe005_sparse_mtx[, -1]

# 17397 genes x 3486 cells
dim(rpe005_sparse_mtx)

##################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 16670 genes x 3486 cells
rpe005_seurat_sct = CreateSeuratObject(counts = rpe005_sparse_mtx, 
                                       project = "pan_rpe005", 
                                       min.cells = 10, 
                                       min.features = 200)

##################################
# STEP 3
# Detect doublets with scDblFinder
# Called 60 doublets
rpe005_doublet_ids = scDblFinder_clusters(rpe005_seurat_sct, nrep = 3)
length(rpe005_doublet_ids)

# First filtering step: Use the vector of consensus doublet IDs output below
rpe005_singlets = setdiff(colnames(rpe005_sparse_mtx), rpe005_doublet_ids)
length(rpe005_singlets)

# New dims after removing doublets: 16640 genes x 3436 cells
rpe005_seurat_sct = CreateSeuratObject(counts = rpe005_sparse_mtx[, rpe005_singlets], 
                                       project = "pan_rpe005", 
                                       min.cells = 10, 
                                       min.features = 200)

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
rpe005_seurat_sct =  decontX_remove(rpe005_seurat_sct)
head(rpe005_seurat_sct@assays$RNA@counts)

###################################
# STEP 2 and 4 before/after removing doublets
# SCT normalize
rpe005_seurat_sct = Seurat_SCT_process(rpe005_seurat_sct, seurat_filter = TRUE,
                                       sample_id = "pan_rpe005", 
                                       study_name = "pan_et_al")

# Finding the right resolution is quite important for doublet detection and removal
rpe005_seurat_sct =  FindClusters(rpe005_seurat_sct, resolution = 0.5)   

# Visualize clusters
p2_before_QC = DimPlot(rpe005_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

p2_after_QC = DimPlot(rpe005_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")

FeaturePlot(rpe005_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                            "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)

FeaturePlot(rpe005_seurat_sct, features = c("CD68", "CD14", "CD163", "APOE", 
                                            "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)

# Save RDS object
# Number of cells after processing: 3436
saveRDS(rpe005_seurat_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/pan_et_al_rds_objects/pan_rpe005_processed_seurat_obj.rds")



##########################################################################
# Read Sample 3: GSM4705591 RPE006 (Caucasian, Female, 76, Asymptomatic) #
##########################################################################

# Load Matrix
rpe006_matrix = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Pan_et_al_human_carotids_scRNA_data/GSE155512_RAW/GSM4705591_RPE006_matrix.txt", sep="\t", header = TRUE)

rownames(rpe006_matrix) = rpe006_matrix$gene
head(rownames(rpe006_matrix))
dim(rpe006_matrix)

# Make sparse matrix and remove "gene" column
rpe006_sparse_mtx = as.sparse(rpe006_matrix)
rpe006_sparse_mtx = rpe006_sparse_mtx[, -1]

# 15687 genes x 2767 cells
dim(rpe006_sparse_mtx)

##################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 15445 genes x 2767 cells
rpe006_seurat_sct = CreateSeuratObject(counts = rpe006_sparse_mtx, 
                                       project = "pan_rpe006", 
                                       min.cells = 10, 
                                       min.features = 200)

##################################
# STEP 3
# Detect doublets with scDblFinder
rpe006_doublet_ids = scDblFinder_clusters(rpe006_seurat_sct, nrep = 3)
length(rpe006_doublet_ids)

# First filtering step: Use the vector of consensus doublet IDs output below
rpe006_singlets = setdiff(colnames(rpe006_sparse_mtx), rpe006_doublet_ids)
length(rpe006_singlets)

# New dims after removing doublets: 15406 genes x 2730 cells
rpe006_seurat_sct = CreateSeuratObject(counts = rpe006_sparse_mtx[, rpe006_singlets], 
                                       project = "pan_rpe006", 
                                       min.cells = 10, 
                                       min.features = 200)

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
rpe006_seurat_sct =  decontX_remove(rpe006_seurat_sct)
head(rpe006_seurat_sct@assays$RNA@counts)

###################################
# STEP 2 and 4 before/after removing doublets
# SCT normalize
rpe006_seurat_sct = Seurat_SCT_process(rpe006_seurat_sct, seurat_filter = TRUE,
                                       sample_id = "pan_rpe006", 
                                       study_name = "pan_et_al")

# Finding the right resolution is quite important for doublet detection and removal
rpe006_seurat_sct =  FindClusters(rpe006_seurat_sct, resolution = 0.5)   

# Visualize clusters
p3_before_QC = DimPlot(rpe006_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

p3_after_QC = DimPlot(rpe006_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")


FeaturePlot(rpe006_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                            "CRTAC1", "VCAM1"), order = TRUE,  pt.size = 0.1)

FeaturePlot(rpe006_seurat_sct, features = c("CD68", "CD14", "CD163", "APOE", 
                                            "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)

# Save RDS object
# Number of cells after processing: 2730
saveRDS(rpe006_seurat_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/pan_et_al_rds_objects/pan_rpe006_processed_seurat_obj.rds")


#################################################################

##########################################################################################
# INTEGRATION WORKFLOW FOR PAN ET AL HUMAN CAROTID SAMPLES PROCESSED WITH SCTransform    #
##########################################################################################

# Prep datasets for integration
seurat_obj_list = list(rpe004=rpe004_seurat_sct,
                       rpe005=rpe005_seurat_sct,
                       rpe006=rpe006_seurat_sct)
int_features_sct = SelectIntegrationFeatures(seurat_obj_list, nfeatures = 3000)
seurat_obj_list = PrepSCTIntegration(seurat_obj_list, anchor.features = int_features_sct)

# Find integration anchors and integrate data
int_sct_anchors = FindIntegrationAnchors(object.list = seurat_obj_list, 
                                          normalization.method = "SCT",
                                          anchor.features = int_features_sct)
carotid_int_sct = IntegrateData(anchorset = int_sct_anchors, normalization.method = "SCT")

# Run PCA and UMAP on integrated object (default is now integrated object)
DefaultAssay(carotid_int_sct) = "integrated"
carotid_int_sct = RunPCA(carotid_int_sct)
carotid_int_sct = FindNeighbors(carotid_int_sct, reduction = "pca", dims = 1:30, k.param = 20)
carotid_int_sct = RunUMAP(carotid_int_sct, dims = 1:30, n.neighbors = 20)
carotid_int_sct = FindClusters(carotid_int_sct, resolution = 1)

# Visualization of clusters
p1_sct = DimPlot(carotid_int_sct, reduction = "umap", group.by = "sample", pt.size = 0.1)
p2_sct = DimPlot(carotid_int_sct, reduction = "umap", label = TRUE, repel = TRUE)
p1_sct + p2_sct


# Switch back to "RNA" or SCT assay to visualize gene expression or perform DE analyses.
# Integrated assay data is heavily corrected so not appropriate for those analysis.

# The ‘corrected’ UMI counts are stored in [["SCT"]]@counts assuming all cells were sequenced at the
# same depth. 
# You can use the corrected log-normalized counts [["SCT"]]@data for differential expression, visualization 
# and integration. 
# However, in principle, it would be most optimal to perform these calculations directly on the 
# residuals (stored in the scale.data slot) themselves. 

DefaultAssay(carotid_int_sct) = "SCT"
FeaturePlot(carotid_int_sct ,features = c("MYH11", "CNN1", "TNFRSF11B", "VCAM1", 
                                            "CRTAC1", "KRT17"), 
            order = TRUE,  pt.size = 0.1, slot = "data")
FeaturePlot(carotid_int_sct ,features = c("MYH11", "CNN1", "TNFRSF11B", "VCAM1", 
                                          "CRTAC1", "VIM"), 
            order = TRUE,  pt.size = 0.1, slot = "data")

FeaturePlot(carotid_int_sct, features = c("CD68", "CD14", "CD163", "APOE", 
                                            "S100A8", "IL1B"), 
            order = TRUE,  pt.size = 0.1)

VlnPlot(carotid_int_sct, features = c("MYH11"), 
        assay = "SCT", slot = "data")

# Run a DE analysis to get cluster markers and annotate cell types. 
# We could say that we use a Wilcox test to get cluster markers. 
# Could always use a different test for more targeted DE analysis between 2 cell types. 
carotid_markers = FindAllMarkers(carotid_int_sct,
                                 logfc.threshold = 0.5, 
                                 only.pos = TRUE,
                                 min.pct = 0.25)

c1_markers = carotid_markers[carotid_markers$cluster==1, ] # Diff SMC
c1_markers %>%
  arrange(desc(avg_log2FC)) %>%
  head(n=100)

c3_markers = carotid_markers[carotid_markers$cluster==3, ] # Diff SMC
c3_markers %>%
  arrange(desc(avg_log2FC)) %>%
  head(n=100)

# Decreased expression of SMC markers and gradual increase of ECM markers, VCAM1
c4_markers = carotid_markers[carotid_markers$cluster==4, ] # VCAM1 ECM SMC
c4_markers %>%
  arrange(desc(avg_log2FC)) %>%
  head(n=100)

c8_markers = carotid_markers[carotid_markers$cluster==8, ] # VCAM1 ECM SMC
c8_markers %>%
  arrange(desc(avg_log2FC)) %>%
  head(n=100)

# Higher expression of KRT17, OGN, VCAN, IGFBP7, TNFRSF11B.
# There are likely counterparts of the Wirka et al Fibromyo
c12_markers = carotid_markers[carotid_markers$cluster==12, ] # ECM SMC (Fibromyo?)
c12_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

# Higher expression of LUM, IGFBP2, TNFRSF11B, COL1A1, COL1A2, FN1, VCAN
c2_markers = carotid_markers[carotid_markers$cluster==2, ] # ECM SMC (Fibromyo?)
c2_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

# High expression of IGFBP2, TNFAIP6, PRELP, SPARC, FN1, CRTAC1, POSTN, PCOLCE2, LUM, DCN
c13_markers = carotid_markers[carotid_markers$cluster==13, ] # Fibrochondrocytes
c13_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

# This cluster shows almost no expression of SMC markers, high expression of APOE 
# and ECM markers. These might be fibroblasts
c10_markers = carotid_markers[carotid_markers$cluster==10, ] # APOE ECM cells (Fibroblasts?)
c10_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c9_markers = carotid_markers[carotid_markers$cluster==9, ] # Endothelial1
c9_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c5_markers = carotid_markers[carotid_markers$cluster==5, ] # Endothelial1
c5_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c14_markers = carotid_markers[carotid_markers$cluster==14, ] # Endothelial2
c14_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c15_markers = carotid_markers[carotid_markers$cluster==15, ] # T/NK cells
c15_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c7_markers = carotid_markers[carotid_markers$cluster==7, ] # T/NK cells
c7_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c16_markers = carotid_markers[carotid_markers$cluster==16, ] # Mast cell
c16_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

# Higher CD163, APOE, APOC1
c6_markers = carotid_markers[carotid_markers$cluster==6, ] # Macrophage1 (Foamy mac?)
c6_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

# More specific expression of APOE, APOC1, CD74, CD163
c18_markers = carotid_markers[carotid_markers$cluster==18, ] # Macrophage1 (Foamy mac?)
c18_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

# Increasing expression of IL1B
c0_markers = carotid_markers[carotid_markers$cluster==0, ] # Macrophage2 (likely pro-inflammatory mac)
c0_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c11_markers = carotid_markers[carotid_markers$cluster==11, ] # Macrophage3
c11_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

# Specific expression of S100A8 and S100A9.Also high expression of IL1B
c20_markers = carotid_markers[carotid_markers$cluster==20, ] # Macrophage4
c20_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

# Mostly distinguished by expression of RPS/RPL proteins
c17_markers = carotid_markers[carotid_markers$cluster==17, ] # Cluster17
c17_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c19_markers = carotid_markers[carotid_markers$cluster==19, ] # Cluster19
c19_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

##################################
# Assign cells types to clusters
diff_smc = c("1", "3")
vcam1_h_ecm_smc = c("4", "8")
vcam1_l_ecm_smc = c("2", "12")
FC = c("13")
apoe_ecm = c("10")
endothelial1 = c("5", "9")
endothelial2 = c("14")
t_nk_cell = c("7", "15")
macrophage1 = c("6", "18")
macrophage2 = c("0")
macrophage3 = c("11")
macrophage4 = c("20")
mast_cell = c("16")
cluster17 = c("17")
cluster19 = c("19")

carotid_int_sct@meta.data = carotid_int_sct@meta.data %>%
  mutate(manually_annotated_cell_types = case_when(seurat_clusters %in% diff_smc ~ "SMC",
                                                   seurat_clusters %in% vcam1_h_ecm_smc ~ "VCAM1_ECM_SMC",
                                                   seurat_clusters %in% vcam1_l_ecm_smc ~ "ECM_SMC (Fibromyo)",
                                                   seurat_clusters %in% FC ~ "Fibrochondrocyte",
                                                   seurat_clusters %in% apoe_ecm ~ "APOE_ECM_rich",
                                                   seurat_clusters %in% endothelial1 ~ "Endothelial1",
                                                   seurat_clusters %in% endothelial2 ~ "Endothelial2",
                                                   seurat_clusters %in% mast_cell ~ "Mast_cell",
                                                   seurat_clusters %in% macrophage1 ~ "Macrophage1",
                                                   seurat_clusters %in% macrophage2 ~ "Macrophage2",
                                                   seurat_clusters %in% macrophage3 ~ "Macrophage3",
                                                   seurat_clusters %in% macrophage4 ~ "Macrophage4",
                                                   seurat_clusters %in% t_nk_cell ~ "T/NK_cell",
                                                   seurat_clusters %in% cluster17 ~ "Cluster17",
                                                   seurat_clusters %in% cluster19 ~ "Cluster19"))
DimPlot(carotid_int_sct, group.by = "manually_annotated_cell_types", 
        repel = TRUE, label = TRUE)

#################
# Save RDS object
saveRDS(carotid_int_sct, "~/Desktop/SMC_modulation_project/human_scRNA_meta-analysis/rds_objects/pan_et_al_rds_objects/rpe00456_integrated_samples_annotated_seurat_SCTransform_obj.rds")








