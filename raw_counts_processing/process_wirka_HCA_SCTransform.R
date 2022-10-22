library(Seurat)
library(tidyverse)
library(scDblFinder)
library(celda)
library(sctransform)
library(cluster)


# Set seed for doublet removal reproducibility
set.seed(1)


# PROCESS WIRKA ET AL DATA USING THE SCTRANSFORM PIPELINE

# Read the Wirka HCA counts matrix
wirka_sparse_counts = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Wirka_human_coronaries_scRNA_data/scRNA_raw_sparse_counts.rds")
wirka_sparse_counts[1:4, 1:4]

# Dims: 20431 genes x 11756 cells
dim(wirka_sparse_counts)

################################################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 18194 genes x 11756 cells
wirka_seurat_sct = CreateSeuratObject(counts = wirka_sparse_counts, 
                                      project = "wirka_human_coronaries_8patients", 
                                      min.cells = 10, 
                                      min.features = 200)

# Assign names to the 8 samples that were merged. Sample names are at the very end of
# each cell barcode
sample_numbers = stringr::str_sub(rownames(wirka_seurat_sct@meta.data), 
                                  start = -1)
sample_ids = paste("wirka", sample_numbers, sep = "_")
table(sample_ids)

##################################
# STEP 3
# Detect doublets with scDblFinder
wirka_doublet_ids = scDblFinder_clusters(wirka_seurat_sct, nrep = 3)
length(wirka_doublet_ids)

# First filtering step: Use the vector of consensus doublet IDs output below
wirka_singlets = setdiff(colnames(wirka_sparse_counts), wirka_doublet_ids)
length(wirka_singlets)

# New dims after removing doublets: 18059 genes x   11199 cells
wirka_seurat_sct = CreateSeuratObject(counts = wirka_sparse_counts[, wirka_singlets], 
                                      project = "wirka_human_coronaries_8patients", 
                                      min.cells = 10, 
                                      min.features = 200)

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
wirka_seurat_sct =  decontX_remove(wirka_seurat_sct)
head(wirka_seurat_sct@assays$RNA@counts)

# Generate sample names for clean matrix
sample_numbers = stringr::str_sub(rownames(wirka_seurat_sct@meta.data), 
                                  start = -1)
sample_ids = paste("wirka_coronary", sample_numbers, sep = "_")
table(sample_ids)

###################################
# STEP 2 and 4 before/after removing doublets
# SCT normalize
wirka_seurat_sct = Seurat_SCT_process(wirka_seurat_sct, seurat_filter = TRUE,
                                      sample_id = sample_ids, 
                                      study_name = "wirka_et_al")

# Finding the right resolution is quite important for doublet detection and removal
wirka_seurat_sct =  FindClusters(wirka_seurat_sct, resolution = 0.9)   

# Visualize clusters
p1_before_QC = DimPlot(wirka_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

p2_after_QC = DimPlot(wirka_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")

FeaturePlot(wirka_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                           "CRTAC1", "LUM"), order = TRUE, pt.size = 0.1)
FeaturePlot(wirka_seurat_sct, features = c("MYH11", "TNFRSF11B", "CRTAC1", "KRT17", "VCAM1", "FHL5"), order = TRUE,  pt.size = 0.1)

DimPlot(wirka_seurat_sct, group.by = "sample")

FeaturePlot(wirka_seurat_sct, features = c("CD68", "CD14", "LYZ", "APOE", 
                                           "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)

FeaturePlot(wirka_seurat_sct, features = c("CD68", "CD3G", "LYZ", "S100A9", 
                                           "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)


# Save rds object
# Number of cells after processing: 10944
saveRDS(wirka_seurat_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/wirka_et_al_rds_objects/wirka_human_coronaries_processed_seurat_obj.rds")



# Run decontX individually on this dataset
wirka_sce = as.SingleCellExperiment(wirka_seurat_sct)
wirka_sce_decontX = celda::decontX(wirka_sce)
umap = reducedDim(wirka_sce_decontX, "decontX_UMAP")
plotDimReduceCluster(x = wirka_sce_decontX$decontX_clusters,
                     dim1 = umap[, 1], dim2 = umap[, 2])
plotDecontXContamination(wirka_sce_decontX)





# Run a DE analysis to get cluster markers and annotate cell types. 
# We could say that we use a Wilcox test to get cluster markers. 
# Could always use a different test for more targeted DE analysis between 2 cell types. 
coronary_markers = FindAllMarkers(wirka_seurat_sct,
                                 logfc.threshold = 1, 
                                 only.pos = TRUE,
                                 min.pct = 0.25)

c6_markers = coronary_markers[coronary_markers$cluster==6, ] # Diff SMC
c6_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c3_markers = coronary_markers[coronary_markers$cluster==3, ] # Fibromyocytes
c3_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

# Lower expression of VCAN, COL1A1, COL3A1. High expression of DPT, PRELP, PCOLCE2, SERPINE2.
# FCs are embedded within this cluster. There's increased expression of SOX9 compared to
# SMCs and Fibromyocytes
c7_markers = coronary_markers[coronary_markers$cluster==7, ] # Fibrochondrocytes
c7_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c0_markers = coronary_markers[coronary_markers$cluster==0, ] # Fibroblasts
c0_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c5_markers = coronary_markers[coronary_markers$cluster==5, ] # Pericyte1
c5_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c8_markers = coronary_markers[coronary_markers$cluster==8, ] # Pericyte2
c8_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c1_markers = coronary_markers[coronary_markers$cluster==1, ] # Endothelial
c1_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c4_markers = coronary_markers[coronary_markers$cluster==4, ] # T cell 
c4_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c14_markers = coronary_markers[coronary_markers$cluster==14, ] # NK cell 
c14_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c23_markers = coronary_markers[coronary_markers$cluster==23, ] # Cluster23
c23_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c9_markers = coronary_markers[coronary_markers$cluster==9, ] # B cell
c9_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c17_markers = coronary_markers[coronary_markers$cluster==17, ] #  
c17_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c18_markers = coronary_markers[coronary_markers$cluster==18, ] # Plasma cell 1
c18_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c21_markers = coronary_markers[coronary_markers$cluster==21, ] # B cell
c21_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c20_markers = coronary_markers[coronary_markers$cluster==20, ] # Plasma cell 2
c20_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c22_markers = coronary_markers[coronary_markers$cluster==22, ] # Mast cell
c22_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

# Specific expression of APOE and APOC1. These might be foamy macarophages
c10_markers = coronary_markers[coronary_markers$cluster==10, ] # Macrophage1 (APOE)
c10_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c2_markers = coronary_markers[coronary_markers$cluster==2, ] # Macrophage2 (IL1B)
c2_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c12_markers = coronary_markers[coronary_markers$cluster==12, ] # Macrophage3 
c12_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

# Highest expression of S100A8 and S1009
c16_markers = coronary_markers[coronary_markers$cluster==16, ] # Macrophage4 (S100A8)
c16_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c17_markers = coronary_markers[coronary_markers$cluster==17, ] # Macrophage5
c17_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c11_markers = coronary_markers[coronary_markers$cluster==11, ] # Wirka_Cluster11
c11_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c13_markers = coronary_markers[coronary_markers$cluster==13, ] # Endothelial2
c13_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c15_markers = coronary_markers[coronary_markers$cluster==15, ] # Neuron
c15_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c19_markers = coronary_markers[coronary_markers$cluster==19, ] # Cluster19
c19_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

c24_markers = coronary_markers[coronary_markers$cluster==24, ] # Cluster24
c24_markers %>% arrange(desc(avg_log2FC)) %>% head(n=100)

##################################
# Assign cells types to clusters

smc = c("6")
fibromyo = c("3")
fibrochondro = c("7")
fibroblast = c("0")
pericyte1 = c("5")
pericyte2 = c("8")
endothelial1 = c("1")
endothelial2 = c("13")
t_cell = c("4")
nk_cell = c("14")
b_cell = c("9", "21")
plasma1 = c("18")
plasma2 = c("20")
mast = c("22")
macrophage1 = c("10")
macrophage2 = c("2")
macrophage3 = c("12")
macrophage4 = c("16")
macrophage5 = c("17")
neuron = c("15")
c11 = c("11")
c19 = c("19")
c23 = c("23")
c24 = c("24")

wirka_seurat_sct@meta.data = wirka_seurat_sct@meta.data %>%
  mutate(manually_annotated_cell_types = case_when(seurat_clusters %in% smc ~ "SMC",
                                                   seurat_clusters %in% fibromyo ~ "Fibromyocytes",
                                                   seurat_clusters %in% fibrochondro ~ "Fibrochondrocytes",
                                                   seurat_clusters %in% fibroblast ~ "Fibroblasts",
                                                   seurat_clusters %in% pericyte1 ~ "Pericyte1",
                                                   seurat_clusters %in% pericyte2 ~ "Pericyte2",
                                                   seurat_clusters %in% endothelial1 ~ "Endothelial1",
                                                   seurat_clusters %in% endothelial2 ~ "Endothelial2",
                                                   seurat_clusters %in% t_cell ~ "T_cell",
                                                   seurat_clusters %in% nk_cell ~ "NK_cell",
                                                   seurat_clusters %in% b_cell ~ "B_cell",
                                                   seurat_clusters %in% plasma1 ~ "Plasma_cell1",
                                                   seurat_clusters %in% plasma2 ~ "Plasma_cell2",
                                                   seurat_clusters %in% mast ~ "Mast_cells",
                                                   seurat_clusters %in% macrophage1 ~ "Macrophage1(APOE)",
                                                   seurat_clusters %in% macrophage2 ~ "Macrophage2",
                                                   seurat_clusters %in% macrophage3 ~ "Macrophage3",
                                                   seurat_clusters %in% macrophage4 ~ "Macrophage4(S100A8)",
                                                   seurat_clusters %in% macrophage5 ~ "Macrophage5",
                                                   seurat_clusters %in% neuron ~ "Neuron",
                                                   seurat_clusters %in% c11 ~ "Wirka_cluster11",
                                                   seurat_clusters %in% c19 ~ "Wirka_cluster19",
                                                   seurat_clusters %in% c23 ~ "Wirka_cluster23",
                                                   seurat_clusters %in% c24 ~ "Wirka_cluster24"))

DimPlot(wirka_seurat_sct, group.by = "manually_annotated_cell_types", 
        repel = TRUE, label = TRUE)                                                   

#################
# Save RDS object
saveRDS(wirka_seurat_sct, "~/Desktop/SMC_modulation_project/human_scRNA_meta-analysis/rds_objects/wirka_et_al_rs_objects/Wirka_HCA_annotated_seurat_SCTransform_obj.rds")                                                   
                                                   
                                                   
                                                   
                                                   





