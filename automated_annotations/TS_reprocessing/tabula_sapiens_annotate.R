library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(data.table)

####################################################################################################################
# This script contains code to transfer cell type labels from different subsets of the Tabula Sapiens reference    #
# to our meta-analyzed scRNA dataset. Since TS data was normalized using Log normalization, we will have to        #
# reprocess the data using SCTransform. We're mostly interested in immune and vascular cell types so we'll focus   #
# on normalizing the TS vasculature and immune subsets. We'll then use the reprocessed data for label transfer     #
####################################################################################################################

##################################################
# Convert h5ad Tabula Sapiens file to h5Seurat
SeuratDisk::Convert("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/annotation_references/Tabula_sapiens/TabulaSapiens.h5ad",
                    dest = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/annotation_references/Tabula_sapiens/TabulaSapiens.h5seurat",
                    overwrite=TRUE,
                    verbose=TRUE)


# Connect to h5seurat file
ts_file = SeuratDisk::Connect("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/annotation_references/Tabula_sapiens/TabulaSapiens.h5seurat")
ts_file

# Look at summary of the data stored within the h5seurat object. 
ts_file$index()
ts_seurat_obj = LoadH5Seurat("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/annotation_references/Tabula_sapiens/TabulaSapiens.h5seurat",
                             assays="RNA")
saveRDS(ts_seurat_obj, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/reference_annotations_rds_objects/Tabula_sapiens/Tabula_sapiens_whole_seurat_obj.rds")

# Dims of TS entire reference: 58870 genes x 481120 cells
ts_seurat_obj
DimPlot(ts_seurat_obj, label = FALSE, repel = TRUE, reduction = "umap")
FeaturePlot(ts_seurat_obj, features = c("MYH11"), raster = FALSE)



##################################################
# Convert Tabula sapiens immune  h5ad to h5seurat
SeuratDisk::Convert("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/annotation_references/Tabula_sapiens/TS_immune.h5ad",
                    dest="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/annotation_references/Tabula_sapiens/TabulaSapiens.immune.h5seurat")


# Look at summary data and load as seurat object. 
# TS_immune has 264009 cells
ts_immune_seurat = LoadH5Seurat("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/annotation_references/Tabula_sapiens/TabulaSapiens.immune.h5seurat",
                                assays="RNA")
saveRDS(ts_immune_seurat, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/reference_annotations_rds_objects/Tabula_sapiens/Tabula_sapiens_immune_seurat_obj.rds")
TS_immune = readRDS("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/reference_annotations_rds_objects/Tabula_sapiens/Tabula_sapiens_immune_seurat_obj.rds")
DimPlot(ts_immune_seurat, label = TRUE, group.by = "cell_ontology_class", raster = FALSE, repel = TRUE, label.size = 3.5) + custom_theme +
  theme(legend.position = "none")
DimPlot(ts_immune_seurat)
table(TS_immune@meta.data$organ_tissue)

# Since we have 200k cells, it might be better to subset to relevant tissues before SCTransform re-normalizing
# List of relevant tissues: Blood, Bone marrow, Vasculature, potentially lymph node and thymus. This subset includes 142261 cells (This crashes R session)
# List of relevant tissues: Blood, Bone marrow, Vasculature (66542 cells)
relevant_tissues = c("Blood", "Bone_Marrow", "Vasculature", "Thymus")
TS_immune_tissue_subset = subset(TS_immune, subset = organ_tissue %in% relevant_tissues)

# SCTransform normalize prior to running label transfer workflow
TS_immune_tissue_subset = SCTransform(TS_immune_tissue_subset, vst.flavor = "v2")
TS_immune_tissue_subset =  RunPCA(TS_immune_tissue_subset) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30, k.param = 20) %>%
  RunUMAP(dims = 1:30, n.neighbors = 30) 
TS_immune_tissue_subset = FindClusters(TS_immune_tissue_subset, resolution = 1)

# Save TS immune subset for Blood, Bone Marrow and Vasculature
saveRDS(TS_immune_tissue_subset, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/reference_annotations_rds_objects/Tabula_sapiens/TS_immune_seurat_obj_Blood_BM_Vasc_SCT.rds")
TS_immune_tissue_subset = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/reference_annotations_rds_objects/Tabula_sapiens/TS_immune_seurat_obj_Blood_BM_Vasc_SCT.rds")

# Plot SCT re-normalized data 
DimPlot(TS_immune_tissue_subset, label = TRUE, group.by = "cell_ontology_class", raster = FALSE, repel = TRUE) + custom_theme +
  theme(legend.position = "none") +
  ggtitle("TS immune for Blood, Bone marrow and Vasculature")
FeaturePlot(TS_immune_tissue_subset, features = c("IL1B", "NKG7","APOE", "CD3E"), order = TRUE) & new_scale & custom_theme

#####################################################
# Convert Tabula sapiens vasculature h5ad to h5seurat
SeuratDisk::Convert("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/annotation_references/Tabula_sapiens/TS_Vasculature.h5ad",
                    dest = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/annotation_references/Tabula_sapiens/TabulaSapiens.vasculature.h5seurat",
                    overwrite=TRUE,
                    verbose=TRUE)

TS_vasculature = LoadH5Seurat("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/annotation_references/Tabula_sapiens/TabulaSapiens.vasculature.h5seurat",
                              assays="RNA")
saveRDS(TS_vasculature, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/reference_annotations_rds_objects/Tabula_sapiens/Tabula_sapiens_vasculature_seurat_obj.rds")
TS_vasculature = readRDS("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/reference_annotations_rds_objects/Tabula_sapiens/Tabula_sapiens_vasculature_seurat_obj.rds")

table(TS_vasculature@meta.data$cell_ontology_class)
DimPlot(TS_vasculature, label = TRUE, group.by = "cell_ontology_class", raster = FALSE, repel = TRUE, label.size = 3.5) + custom_theme + 
  theme(legend.position = "none")
FeaturePlot(TS_vasculature, features = c("CNN1", "ACTA2","IL1B", "APOE"), order = TRUE) & new_scale & custom_theme

# Normalize data using SCTransform prior to label transfer so that reference and input datasets are normalized in the same manner#
# NOTE: Looks that for FindTransferAnchors to work both reference and query datasets have to be normalized using the same approach.
# Alternative solution: SCTransform tabula sapiens data. Cluster topology should fairly similar and since we already have the
# cell labels, in theory we should be able to run the label transfer workflow. 
TS_vasculature = SCTransform(TS_vasculature, vst.flavor = "v2")
TS_vasculature =  RunPCA(TS_vasculature) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30, k.param = 20) %>%
  RunUMAP(dims = 1:30, n.neighbors = 30) 
TS_vasculature = FindClusters(TS_vasculature, resolution = 1)

# Save TS vasculature subset
saveRDS(TS_vasculature, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/reference_annotations_rds_objects/Tabula_sapiens/TS_vasculature_seurat_obj_SCT.rds")
TS_vasculature = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/reference_annotations_rds_objects/Tabula_sapiens/TS_vasculature_seurat_obj_SCT.rds")

# Plot SCT re-normalized data 
DimPlot(TS_vasculature, label = TRUE, group.by = "cell_ontology_class", raster = FALSE, repel = TRUE) + custom_theme
FeaturePlot(TS_vasculature, features = c("CNN1", "ACTA2","IL1B", "APOE"), order = TRUE) & new_scale & custom_theme


#########################################################
# Transfer labels from TS subsets to meta-analyzed data #
#########################################################

# Before transferring anchors, make sure that the default assay of the integrated reference is set to SCT
DefaultAssay(rpca_int_sct_v3) = "SCT"

############################################
# Find transfer anchors for TS vasculature
# From Seurat issue 3937: As of now, you can only run our transfer workflow with normalization.method = SCT when both query and reference have been SCTransformed. 
# I suggest using the default normalization.method if you are trying to transfer from one integrated object to another.
transfer_anchors_vasc = FindTransferAnchors(reference = TS_vasculature, query = rpca_int_sct_v3, dims = 1:30, 
                                            reference.reduction = "pca", 
                                            normalization.method = "SCT")
vasc_predictions = TransferData(anchorset = transfer_anchors_vasc, 
                           refdata = TS_vasculature$cell_ontology_class, 
                           dims = 1:30)
ts_vasc_predictions = vasc_predictions %>%
  dplyr::select(predicted.id, prediction.score.max) %>%
  rename(TS_vasc_predicted_id = predicted.id,
         TS_vasc_prediction_score_max = prediction.score.max)
  
rpca_int_sct_v3 = AddMetaData(rpca_int_sct_v3, metadata = ts_vasc_predictions)

# Plot transferred labels to Wirka hca data
DimPlot(rpca_int_sct_v3, label = TRUE, group.by = "TS_vasc_predicted_id", raster = FALSE, repel = TRUE) + custom_theme + 
  ggtitle("Predicted ID from TS vasculature") +
  theme(legend.position = "bottom")
ggplot(rpca_int_sct_v3@meta.data, aes(x=TS_vasc_prediction_score_max)) + geom_histogram(bins = 70) + custom_theme + 
  ggtitle("Prediction scores from TS vasculature") +
  xlab("Prediction scores from TS vasculature") +
  ylab("Cell number")


#############################################
# Find transfer anchors for TS immune subset
# From Seurat issue 3937: As of now, you can only run our transfer workflow with normalization.method = SCT when both query and reference have been SCTransformed. 
# I suggest using the default normalization.method if you are trying to transfer from one integrated object to another.
transfer_anchors_immune = FindTransferAnchors(reference = TS_immune_tissue_subset, query = rpca_int_sct_v3, dims = 1:30, reference.reduction = "pca", 
                                             normalization.method = "SCT")
immune_predictions = TransferData(anchorset = transfer_anchors_immune, refdata = TS_immune_tissue_subset$cell_ontology_class, dims = 1:30)

ts_immune_predictions = immune_predictions %>%
  dplyr::select(predicted.id, prediction.score.max) %>%
  rename(TS_immune_predicted_id = predicted.id,
         TS_immune_prediction_score_max = prediction.score.max)

rpca_int_sct_v3 = AddMetaData(rpca_int_sct_v3, metadata = ts_immune_predictions)
table(rpca_int_sct_v3$TS_immune_predicted_id)

# Plot transferred labels to Wirka hca data
DimPlot(rpca_int_sct_v3, label = TRUE, group.by = "TS_immune_predicted_id", raster = FALSE, repel = TRUE) + custom_theme + 
  ggtitle("Predicted ID from TS immune") +
  theme(legend.position = "bottom")
ggplot(rpca_int_sct_v3@meta.data, aes(x=TS_immune_prediction_score_max)) + geom_histogram(bins = 70) + custom_theme + 
  ggtitle("Prediction scores from TS immune") +
  xlab("Prediction scores from TS immune") +
  ylab("Cell number")

rpca_int_sct_v2@meta.data %>%
  filter(predicted.id == "macrophage") %>%
  ggplot(aes(x=prediction.score.max)) + geom_histogram(bins=50) + custom_theme

# Plot transferred labels to Wirka hca data
DimPlot(rpca_int_sct_v2, label = TRUE, group.by = "predicted.id", raster = FALSE, repel = TRUE) + custom_theme + 
  ggtitle("Predicted ID from TS immune of Blood, Bone Marrow and Vasculature") +
  theme(legend.position = "bottom")
ggplot(rpca_int_sct_v2@meta.data, aes(x=prediction.score.max)) + geom_histogram(bins = 70) + custom_theme +
  ggtitle("prediction scores from TS immune subset (Blood, Bone marrow, Vasculature)")




