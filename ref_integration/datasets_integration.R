library(Seurat)
library(harmony)
library(rliger)
library(tidyverse)
library(data.table)
library(cluster)
library(RColorBrewer)
library(parallel)
library(reticulate)
library(UCell)

# Set seed for reproducibility
set.seed(1)

# Source our own scRNA analysis utils functions
source("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/R_scripts/Utils/scRNA_processing_utils.R")


###############################################################
# INTEGRATION SCRIPT FOR ALSAIGH, PAN, WIRKA AND HU ET AL DATA

#######################################
# Load SCTransform processed datasets #
#######################################

# Load Pan et al data
pan_rpe004 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/pan_et_al_rds_objects/pan_rpe004_processed_seurat_obj.rds")
pan_rpe005 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/pan_et_al_rds_objects/pan_rpe005_processed_seurat_obj.rds")
pan_rpe006 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/pan_et_al_rds_objects/pan_rpe006_processed_seurat_obj.rds")

# Load Alsaigh et al data
alsaigh_ac_p1 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_ac_p1_processed_seurat_obj.rds")
alsaigh_pa_p1 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_pa_p1_processed_seurat_obj.rds")

alsaigh_ac_p2 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_ac_p2_processed_seurat_obj.rds")
alsaigh_pa_p2 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_pa_p2_processed_seurat_obj.rds")

alsaigh_ac_p3 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_ac_p3_processed_seurat_obj.rds")
alsaigh_pa_p3 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_pa_p3_processed_seurat_obj.rds")


# Load Wirka et al
wirka_hca = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/wirka_et_al_rds_objects/wirka_human_coronaries_processed_seurat_obj.rds")

# Load Hu et al data
hu_coronary1_p1 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/hu_et_al_rds_objects/hu_coronary1_p1_processed_seurat_obj.rds")
hu_coronary1_p2 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/hu_et_al_rds_objects/hu_coronary1_p2_processed_seurat_obj.rds")
hu_coronary2_p2 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/hu_et_al_rds_objects/hu_coronary2_p2_processed_seurat_obj.rds")
hu_coronary1_p3 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/hu_et_al_rds_objects/hu_coronary1_p3_processed_seurat_obj.rds")
hu_coronary2_p3 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/hu_et_al_rds_objects/hu_coronary2_p3_processed_seurat_obj.rds")


################################################################
# rPCA                                                         #
################################################################

rpca_obj_list = list(rpe004 = pan_rpe004,
                     rpe005 = pan_rpe005,
                     rpe006 = pan_rpe006,
                     ac_p1 = alsaigh_ac_p1,
                     pa_p1 = alsaigh_pa_p1,
                     ac_p2 = alsaigh_ac_p2,
                     pa_p2 = alsaigh_pa_p2,
                     ac_p3 = alsaigh_ac_p3,
                     pa_p3 = alsaigh_pa_p3,
                     wirka_hca = wirka_hca,
                     hu_coronary1_p1 = hu_coronary1_p1,
                     hu_coronary1_p2 = hu_coronary1_p2,
                     hu_coronary2_p2 = hu_coronary2_p2,
                     hu_coronary1_p3 = hu_coronary1_p3,
                     hu_coronary2_p3 = hu_coronary2_p3)

saveRDS(rpca_obj_list, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/alsaigh_pan_wirka_hu_seurat_list_for_int.rds")


# Generate list of highly variable genes across datasets that will be used for integration 
rpca_features_sct = SelectIntegrationFeatures(rpca_obj_list, nfeatures = 3000)
rpca_obj_list = PrepSCTIntegration(rpca_obj_list, anchor.features = rpca_features_sct)
rpca_obj_list = lapply(rpca_obj_list, FUN = function(x) {
  x = RunPCA(x, features=rpca_features_sct)
}) 


# Start timer 
start_time_rpca = Sys.time() 

# Find integration anchors and integrate data
# rPCA integration is tipically more conservative (less likely to align cells in different states
# like naive and memory T cells) across experiments. The k.anchor parameter allows us to modify 
# the strength of the alignment and is by default set to 5. Let's start with 5 and see if that 
# separates SMCs better from pericytes and fibroblasts. 


# How do integrating through set dims=1:50 compare to our current integration (dims=1:30)? 
# Would this be able to provide any additional information to separate cell types better?
# Setting dims=1:50 doesn't really seem to provide a substantial benefit in separating cell types in low dim space

# Since we now have 15 datasets (~100k cells), we might need to use a larger number of dims
# Dims of integrated data: 121316 cells
rpca_sct_anchors = FindIntegrationAnchors(object.list = rpca_obj_list, 
                                          normalization.method = "SCT",
                                          anchor.features = rpca_features_sct, 
                                          dims = 1:30, 
                                          reduction = "rpca",
                                          k.anchor = 5)
rpca_int_sct = IntegrateData(anchorset = rpca_sct_anchors, 
                             normalization.method = "SCT",
                             dims=1:30)

# Measure time for benchmark
end_time_rpca = Sys.time()
measured_time_rpca = end_time_rpca - start_time_rpca # Time difference of 53.3441 mins

###########################################################################
# Now that we identified integration anchors, we'll process the data for
# downstream analyses. 

# Run PCA and UMAP on integrated object (default is now integrated object)
DefaultAssay(rpca_int_sct) = "integrated"

# Run standard workflow for dimensionality reduction and clustering
# 30 dims captures really most of the variation in the data
rpca_int_sct_v3 = RunPCA(rpca_int_sct_v3) 
rpca_int_sct_v3 = FindNeighbors(rpca_int_sct_v3, reduction = "pca", dims = 1:30, k.param = 20)

# In general, this parameter should be 5-50. Higher values result in more global structure 
# being preserved at the loss of detailed local structure
# UMAP looks better with 30 n.neighbors than 20
# Adjusting min.dist doesn't really help
rpca_int_sct_v3 = RunUMAP(rpca_int_sct_v3, dims = 1:30, n.neighbors = 30, min.dist = 0.3) 

# A resolution of 1 seems appropriate based on previous silhouette analysis 
rpca_int_sct_v3 = FindClusters(rpca_int_sct_v3, resolution = 1)
saveRDS(rpca_int_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int_seurat_clustered.rds")

rpca_int_sct_v3 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/seurat_objects/alsaigh_pan_wirka_hu_int_seurat_annotated_v3.rds")
DefaultAssay(rpca_int_sct_v3) = "integrated" 

# Set clusters to latest version of clustering resolutio; res=1
Idents(rpca_int_sct_v3) = "integrated_snn_res.1"

# Visualization of clusters and other variables defined within metadata
p1_rpca = DimPlot(rpca_int_sct_v3, reduction = "umap", group.by = "sample", pt.size = 0.1, raster = FALSE) + ggtitle("rPCA Alsaigh/Pan/Wirka/Hu samples") + custom_theme
p2_rpca = DimPlot(rpca_int_sct_v3, reduction = "umap", group.by = "study", pt.size = 0.1, raster = FALSE) + ggtitle("rPCA Alsaigh/Pan/Wirka/Hu study") + custom_theme + npg_scale
p3_rpca = DimPlot(rpca_int_sct_v3, reduction = "umap", group.by = "level1_annotations", label = FALSE, repel = TRUE, label.size = 5, 
                  pt.size = 0.1, raster = FALSE) + 
  ggtitle("rPCA Alsaigh/Pan/Wirka/Hu clusters") + 
  custom_theme + npg_scale2 +
  theme(legend.position = "right")
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure1/level1_annotations_UMAP.svg",
       plot = p3_rpca, width = 10, height = 8)
p1_rpca + p2_rpca

# Plot level1 annotations 
level1_annotations = DimPlot(rpca_int_sct_v3, reduction = "umap", group.by = "level1_annotations", label = FALSE, repel = TRUE, 
      pt.size = 0.1, raster = FALSE, label.size = 4) + ggtitle("rPCA Alsaigh/Pan/Wirka/Hu level1 annotations") + custom_theme + 
        theme(legend.position = "none") + npg_scale2 + theme(legend.position = "right")
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure2/Fig2a_level1_annotations_UMAP.pdf",
       plot = level1_annotations, width = 9, height = 9)

# Plot level2 annotations
level2_annotations =  DimPlot(rpca_int_sct_v3, reduction = "umap", group.by = "level2_annotations", label = TRUE, repel = TRUE, 
        pt.size = 0.1, raster = FALSE, label.size = 4) + ggtitle("rPCA Alsaigh/Pan/Wirka/Hu clusters") + custom_theme + 
        theme(legend.position = "none") + level2_annotations_scale 
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure3/Fig3a_level2_annotations_UMAP_no_labels.pdf",
       plot = level2_annotations, width = 9, height = 9)

# Plot level2 annotations by disease status
level1_by_disease = DimPlot(rpca_int_sct_v3, reduction = "umap", group.by = "level1_annotations", split.by = "sample_disease_status", label = FALSE,
        repel = TRUE, pt.size = 0.1, raster = FALSE) + npg_scale2 + custom_theme + 
  theme(legend.position = "none",
        strip.text = element_text(size=17))
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure1/Fig1d_level1_UMAP_by_disease_status.svg",
       plot = level1_by_disease, width = 17, height = 8)

# Plot clusters by disease status
DimPlot(rpca_int_sct_v3, reduction = "umap", split.by = "sample_disease_status", group.by = "level1_annotations",
        label = TRUE, repel = TRUE, pt.size = 0.1, raster = FALSE) + custom_theme

# Plot level1 annotations by arterial bed
DimPlot(rpca_int_sct_v3, reduction = "umap", group.by = "level1_annotations", split.by = "arterial_origin", 
         label = TRUE, repel = TRUE, pt.size = 0.1, raster = FALSE) + npg_scale2 + custom_theme + 
  theme(legend.position = "none",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size=14))

###############################################################
# Visualize gene expression using SCTransform-normalized counts
DefaultAssay(rpca_int_sct_v3) = "SCT"
FeaturePlot(rpca_int_sct_v3 ,features = c("MYH11", "ACTA2", "VCAN", 
                                       "FHL5"), raster = FALSE, order = TRUE,  pt.size = 0.1, slot = "data") & new_scale3 & custom_theme

# SMC clusters
FeaturePlot(rpca_int_sct_v3 ,features = c("CNN1", "MYH11", "LMOD1", "COL6A1", "COL6A2", "COL6A3"), order = TRUE,  pt.size = 0.1, slot = "data", raster = FALSE) & new_scale3 & custom_theme

# Fibromyocyte markers
FeaturePlot(rpca_int_sct_v3 ,features = c("CNN1", "TNFRSF11B", "VCAN", "KRT7"), raster = FALSE,
            order = TRUE,  pt.size = 0.1, slot = "data") & new_scale

# Fibrochondro markers
FeaturePlot(rpca_int_sct_v3 ,features = c("CRTAC1", "IBSP", "COMP", "CYTL1", "FN1", "CLU"), 
            order = TRUE,  pt.size = 0.1, slot = "data", raster = FALSE) & custom_theme & new_scale


# Known markers from mural and immune cell subtypes
acta2 = FeaturePlot(rpca_int_sct_v3, features = c("ACTA2"), raster = FALSE, order = TRUE,
                    pt.size = 0.1) & new_scale3 & custom_theme
ggsave(file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure3/Fig3b_subtype_markers_ACTA2.pdf",
       plot = acta2, width = 9, height = 9)

sele = FeaturePlot(rpca_int_sct_v3, features = c("SELE"), raster = FALSE, order = TRUE,
                    pt.size = 0.1) & new_scale3 & custom_theme
ggsave(file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure3/Fig3b_subtype_markers_SELE.pdf",
       plot = sele, width = 9, height = 9)

xcl1 = FeaturePlot(rpca_int_sct_v3, features = c("XCL1"), raster = FALSE, order = TRUE,
                    pt.size = 0.1) & new_scale3 & custom_theme
ggsave(file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure3/Fig3b_subtype_markers_XCL1.pdf",
       plot = xcl1, width = 9, height = 9)

il1b = FeaturePlot(rpca_int_sct_v3, features = c("IL1B"), raster = FALSE, order = TRUE,
                    pt.size = 0.1) & new_scale3 & custom_theme
ggsave(file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure3/Fig3b_subtype_markers_IL1B.pdf",
       plot = il1b, width = 9, height = 9)

s100a8 = FeaturePlot(rpca_int_sct_v3, features = c("S100A8"), raster = FALSE, order = TRUE,
                    pt.size = 0.1) & new_scale3 & custom_theme
ggsave(file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure3/Fig3b_subtype_markers_S100A8.pdf",
       plot = s100a8, width = 9, height = 9)

spp1 = FeaturePlot(rpca_int_sct_v3, features = c("SPP1"), raster = FALSE, order = TRUE,
                    pt.size = 0.1) & new_scale3 & custom_theme
ggsave(file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure3/Fig3b_subtype_markers_SPP1.pdf",
       plot = spp1, width = 9, height = 9)

cd1c = FeaturePlot(rpca_int_sct_v3, features = c("CD1C"), raster = FALSE, order = TRUE,
                    pt.size = 0.1) & new_scale3 & custom_theme
ggsave(file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure3/Fig3b_subtype_markers_CD1C.pdf",
       plot = cd1c, width = 9, height = 9)

cd8a = FeaturePlot(rpca_int_sct_v3, features = c("CD8A"), raster = FALSE, order = TRUE,
                   pt.size = 0.1) & new_scale3 & custom_theme
ggsave(file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure3/Fig3b_subtype_markers_CD8A.pdf",
       plot = cd8a, width = 9, height = 9)

lyve1 = FeaturePlot(rpca_int_sct_v3, features = c("LYVE1"), raster = FALSE, order = TRUE,
                   pt.size = 0.1) & new_scale3 & custom_theme
ggsave(file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure3/Fig3b_subtype_markers_LYVE1.pdf",
       plot = lyve1, width = 9, height = 9)

marker_genes = ggarrange(acta2, sele, xcl1, il1b, s100a8, 
                         spp1, cd1c, cd8a, lyve1, ncol = 3)


ggsave(file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/collabs/Chani_gene_lookups/ARHGAP42_scale2.pdf",
         plot = arhgap42, width = 10, height = 10)

#######################
# Myeloid clusters #
#######################

# Inflammatory
FeaturePlot(rpca_int_sct_v3 ,features = c("IL1B", "CCL3", "CCL2", "TNF", "NFKBIA", "S100A9"), raster=FALSE,
            order = TRUE,  pt.size = 0.1, slot = "data") & new_scale & custom_theme

# Foamy macs
FeaturePlot(rpca_int_sct_v3 ,features = c("SPP1", "APOE", "APOC1", "FTL", "CTSD", "FABP5"), raster=FALSE,
            order = TRUE,  pt.size = 0.1, slot = "data") & new_scale & custom_theme

# Putative neutrophils or classical monocytes
FeaturePlot(rpca_int_sct_v3 ,features = c("S100A8", "S100A9", "LYZ", "CTSS", "FCN1", "S100A12"), raster=FALSE,
            order = TRUE,  pt.size = 0.1, slot = "data") & new_scale & custom_theme

# Nampt neutrophils
FeaturePlot(rpca_int_sct_v3 ,features = c("S100A8", "S100A9", "MNDA", "NAMPT"), raster=FALSE,
            order = TRUE,  pt.size = 0.1, slot = "data") & new_scale & custom_theme

# Intermediate monocytes
FeaturePlot(rpca_int_sct_v3 ,features = c("LYZ", "HLA-DPB1", "HLA-DRA", "HLA-DPA1", "HLA-DQA1", "LYVE1"), raster=FALSE,
            order = TRUE,  pt.size = 0.1, slot = "data") & new_scale & custom_theme

# Non-classical monocytes
FeaturePlot(rpca_int_sct_v3 ,features = c("HLA-DPA1", "HLA-DRA", "HLA-DRB1", "C1QA", "C1QB", "C1QC"), raster=FALSE,
            order = TRUE,  pt.size = 0.1, slot = "data") & new_scale & custom_theme

# Tissue resident macs
FeaturePlot(rpca_int_sct_v3 ,features = c("LYVE1", "MRC1", "F13A1", "FOLR2", "DAB2", "STAB1"), raster=FALSE,
            order = TRUE,  pt.size = 0.1, slot = "data") & new_scale & custom_theme

# Pericytes
FeaturePlot(rpca_int_sct ,features = c("S100A8", "S100A9", "LYVE1", "TPM2", "NET1", "RERGL"), 
            order = TRUE,  pt.size = 0.1, slot = "data") & new_scale

# Neuron
FeaturePlot(rpca_int_sct_v3, features = c("S100B", "MPZ", "PLP1", "GPM6B"),
            order = TRUE, pt.size = 0.1, slot="data", raster = FALSE) & new_scale

#####################
# Lymphoid clusters #
#####################

# XCL1 NK cells 
FeaturePlot(rpca_int_sct_v3, features = c("XCL1", "XCL2", "NKG7", "GZMB"), raster = FALSE,
            order = TRUE, pt.size = 0.1, slot = "data") & custom_theme & new_scale

# NK cells
FeaturePlot(rpca_int_sct_v3, features = c("NKG7", "CD7", "IFITM1", "GZMH"), raster = FALSE,
            order = TRUE, pt.size = 0.1, slot="data") & new_scale & custom_theme

# C7
FeaturePlot(rpca_int_sct_v3, features = c("CCL5", "CCL4", "GZMK", "CD8A"), raster = FALSE,
            order = TRUE, pt.size = 0.1, slot="data") & new_scale & custom_theme

# C10
FeaturePlot(rpca_int_sct_v3, features = c("CXCR4", "CST7", "CD3E", "BTG1"), raster = FALSE,
            order = TRUE, pt.size = 0.1, slot="data") & new_scale & custom_theme

# C13
FeaturePlot(rpca_int_sct_v3, features = c("LTB", "CD52", "CD2", "TRBC2"), raster = FALSE,
            order = TRUE, pt.size = 0.1, slot="data") & new_scale & custom_theme

# C34
FeaturePlot(rpca_int_sct_v3, features = c("DNAJB1", "HSPA1B", "HSPA1A", "CD69"), raster = FALSE,
            order = TRUE, pt.size = 0.1, slot="data") & new_scale & custom_theme

# C5
FeaturePlot(rpca_int_sct_v3, features = c("IL7R", "ZFP36L2", "CXCR4", "KLRB1"), raster = FALSE,
            order = TRUE, pt.size = 0.1, slot="data") & new_scale & custom_theme

# C28
FeaturePlot(rpca_int_sct_v3, features = c("IL7R", "LTB", "EEF1B2", "CCR7"), raster = FALSE,
            order = TRUE, pt.size = 0.1, slot="data") & new_scale & custom_theme

########################
# Endothelial clusters #
########################

# Inflammatory endo
FeaturePlot(rpca_int_sct_v3, features = c("SELE", "CCL2", "VWF", "ICAM1", "NFKBIA", "IFI27"), raster = FALSE,
            order = TRUE, pt.size = 0.1, slot="data") & new_scale & custom_theme

# Inflammatory endo
FeaturePlot(rpca_int_sct_v3, features = c("FN1", "CRTAC1", "ITLN1", "VWF", "ACTA2", "SULF1"), raster = FALSE,
            order = TRUE, pt.size = 0.1, slot="data") & new_scale & custom_theme

###############
# Fibroblasts #
###############

# Inflammatory endo
FeaturePlot(rpca_int_sct_v3, features = c("APOD", "CFD", "DCN", "GSN", "SERPINF1", "CXCL14"), raster = FALSE,
            order = TRUE, pt.size = 0.1, slot="data") & new_scale & custom_theme

# Myofibroblasts
FeaturePlot(rpca_int_sct_v3, features = c("IGFBP5", "C3", "C7", "ACTA2", "PODN", "PDGFRA"), raster = FALSE,
            order = TRUE, pt.size = 0.1, slot="data") & new_scale3 & custom_theme

######################
# Other immune cells #
######################

# Plasma cells
FeaturePlot(rpca_int_sct_v3 ,features = c("CD79B", "IGLC3", "JCHAIN", "MZB1"), 
            order = TRUE,  pt.size = 0.1, slot = "data", raster = FALSE) & custom_theme & new_scale
plasma_cells = WhichCells(rpca_int_sct_v2, idents = c(32), expression = JCHAIN > 2 | IGLC2 > 2 | MZB1 > 2)
length(plasma_cells) # 431 cells
c35 = WhichCells(rpca_int_sct_v3, idents = c(35)) # 925 cells
c15 = WhichCells(rpca_int_sct_v3, idents = c(15))
cells_to_remove = setdiff(c32, plasma_cells)
DimPlot(rpca_int_sct_v3, cells.highlight = c15, sizes.highlight = 0.1, raster = FALSE) + custom_theme

# cells from c31, c23, c4, 

# Make dot plot for level1 annotations
marker_dotplot = DotPlot(rpca_int_sct_v3, group.by = "level1_annotations", features = c("MYH11", "ACTA2", "TPM2", "LMOD1",
                                                                       "ID4", "RGS16", "NET1", "RERGL",
                                                                       "ACKR1", "VWF", "PECAM1", "CLDN5",
                                                                       "CD14", "CD68", "C1QA", "CD74",
                                                                       "TPSAB1", "TPSB2", "MS4A2", "CPA3",
                                                                       "IL32", "NKG7", "CXCR4", "CD3E",
                                                                       "S100B", "MPZ", "PLP1", "GPM6B",
                                                                       "CFD", "APOD", "FBLN1", "C7",
                                                                       "GZMB", "PTGDS", "IRF7", "HERPUD1",
                                                                       "CD79A", "CD79B", "CD37", "MS4A1",
                                                                       "IGLC2", "IGKC", "IGHM", "JCHAIN"), col.min = 0) +
  new_scale3 + custom_theme + theme(axis.text.x = element_text(angle = 45, hjust=1), axis.text = element_text(size=12)) + coord_flip()
ggsave(file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure1/Fig1b_marker_dot_plot.pdf",
       plot = marker_dotplot, width = 12, height = 10)

# Plot cells by disease status
DimPlot(rpca_int_sct_v3, reduction = "umap", group.by = "level1_annotations", split.by = "sample_disease_status", label = FALSE, repel = TRUE, 
        pt.size = 0.1, raster = FALSE) + custom_theme + npg_scale2

############################################################################################################
# Calculate silhouette scores to determine the best clustering res for whole reference. Test 0.8-1.3 range #
############################################################################################################

res = c(0.8, 0.9, 1, 1.1, 1.2, 1.3)
scores_list = list()
dist_matrix = dist(x = Embeddings(object = rpca_int_sct[["pca"]])[, 1:30])
saveRDS(dist_matrix, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/alsaigh_pan_wirka_hu_rpca_dist_matrix.rds")
dist_matrix = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/alsaigh_pan_wirka_hu_rpca_dist_matrix.rds")
dim(dist_matrix)

# Switch back to "integrated" assay to cluster at different resolutions
for (i in seq_along(res)) {
  # Cluster data
  rpca_int_sct = FindClusters(rpca_int_sct, resolution = res[i])
  
  # Calculate silhouette scores
  scores_list[[i]] = calc_sil_scores(rpca_int_sct, dist_matrix)
  names(scores_list)[[i]] = paste("res", as.character(res[i]), sep = "_")
}

# Create a single df with all the sil scores for the tested resolutions
# Max res
resolutions_sil_df = rbindlist(scores_list, idcol = TRUE)
resolutions_sil_df$method = "rPCA"
dim(resolutions_sil_df)
head(resolutions_sil_df)

###########################################
# Update metadata for integrated datasets #
###########################################

metadata_df = rpca_int_sct_v3@meta.data
dim(metadata_df)

# Group cells by artery of origin
carotids = c("alsaigh_et_al", "pan_et_al")
coronaries = c("wirka_et_al", "hu_et_al")

# Group cells by disease status of sample
lesion = c("alsaigh_et_al", "pan_et_al", "wirka_et_al")
non_lesion = c("hu_et_al")

# Add new metadata col for sex (cells from males = 98930; cells from females = 19648) 
# Sex was defined based on metadata from papers and also XIST expression. 
males = c("pan_rpe004", "pan_rpe005", "alsaigh_ac_p1", "alsaigh_pa_p1", "alsaigh_ac_p2", "alsaigh_pa_p2",
          "alsaigh_ac_p3", "alsaigh_pa_p3", "wirka_coronary_1", "wirka_coronary_2", "wirka_coronary_3",
          "wirka_coronary_4", "wirka_coronary_5", "wirka_coronary_6", "wirka_coronary_7", "hu_coronary1_p1", 
          "hu_coronary1_p2", "hu_coronary2_p2")
females = c("pan_rpe006", "wirka_coronary_8", "hu_coronary1_p3", "hu_coronary2_p3")

# Add new metadata cols and cell barcodes as rownames.
# Not sure why cell barcodes are removed when we create the new variables
new_meta_df  = metadata_df %>% 
  mutate(arterial_origin = case_when(study %in% carotids ~ "carotid",
                                     study %in% coronaries ~ "coronary"),
         sample_disease_status = case_when(study %in% lesion ~ "lesion",
                                           study %in% non_lesion ~ "non_lesion"),
         sex = case_when(sample %in% males ~ "males",
                         sample %in% females ~ "females"))
rownames(new_meta_df) = rownames(rpca_int_sct_v3@meta.data)
head(new_meta_df)

# Add cell type labels (level 1 annotations)
pericytes = c(8)
smc = c(0, 2, 3, 14, 31)
fibroblasts = c(6, 11, 12, 23, 25, 33)
neuron = c(35)
endothelial = c(4, 16, 20, 32, 37)
b_cell = c(22)
pDC = c(38)
plasma_cell = c(36)
mast_cell = c(30)
macrophage = c(1, 15, 17, 18, 21, 24, 27, 29, 39, 40)
t_nk = c(5, 7, 10, 13, 19, 26, 28, 34)

new_meta_df = metadata_df %>%
  mutate(level1_annotations = case_when(integrated_snn_res.1 %in% pericytes ~ "Pericyte",
                                        integrated_snn_res.1 %in% smc ~ "SMC",
                                        integrated_snn_res.1 %in% fibroblasts ~ "Fibroblast",
                                        integrated_snn_res.1 %in% neuron ~ "Neuron",
                                        integrated_snn_res.1 %in% endothelial ~ "Endothelial",
                                        integrated_snn_res.1 %in% b_cell ~ "B_cell",
                                        integrated_snn_res.1 %in% pDC ~ "pDC",
                                        integrated_snn_res.1 %in% plasma_cell ~ "Plasma_cell",
                                        integrated_snn_res.1 %in% mast_cell ~ "Mast_cell",
                                        integrated_snn_res.1 %in% macrophage ~ "Macrophage",
                                        integrated_snn_res.1 %in% t_nk ~ "T_NK",
                                        integrated_snn_res.1 %in% c(9) & TS_vasc_predicted_id == "smooth muscle cell" ~ "SMC",
                                        integrated_snn_res.1 %in% c(9) & TS_vasc_predicted_id == "fibroblast" ~ "Fibroblast",
                                        integrated_snn_res.1 %in% c(9) & TS_vasc_predicted_id == "pericyte cell" ~ "SMC",
                                        integrated_snn_res.1 %in% c(9) & TS_vasc_predicted_id == "t cell" ~ "SMC"))
rownames(new_meta_df) = rownames(rpca_int_sct_v3@meta.data)
rpca_int_sct_v3@meta.data = new_meta_df

# Refactor level 1 annotations to make SMCs stand out more 
rpca_int_sct_v3@meta.data$level1_annotations = factor(rpca_int_sct_v3@meta.data$level1_annotations,
                                                                 levels = c("B_cell", "Endothelial", "Fibroblast", "SMC",
                                                                            "Mast_cell", "Neuron", "Pericyte", "Macrophage",
                                                                            "Plasma_cell", "pDC", "T_NK"))
rpca_int_sct_v3@meta.data$study = factor(rpca_int_sct_v3@meta.data$study,
                                         levels = c("alsaigh_et_al", "pan_et_al", "wirka_et_al", "hu_et_al"))

# Check how many variables we currently have within the metadata: 27
# Available patient metadata variables: arterial_origin, sample_disease_status, sex. 
length(names(rpca_smc_fibro_subset_v3@meta.data))


################################################
# Add cell type labels predicted by celltypist #
################################################

##########################################################################################################################
# Add labels from Immune_All_Low model (no majority voting; immune sub-populations combined from 20 tissues of 19 studies)
immune_all_low = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/celltypist_outputs/Immune_All_Low/celltypist_All_Immune_Low_predicted_labels.csv", 
                       header = TRUE, quote = FALSE)
names(immune_all_low)[1] = "cell_barcode"
head(immune_all_low)

immune_all_low_vec = immune_all_low$predicted_labels
names(immune_all_low_vec) = immune_all_low$cell_barcode
head(immune_all_low_vec)

# Add celltypist Immune_All_Low labels into metadata
metadata_df$celltypist_Immune_All_Low = immune_all_low_vec[match(rownames(metadata_df), 
                                                                 names(immune_all_low_vec))]
rpca_int_sct_v3@meta.data = metadata_df
DimPlot(rpca_int_sct_v3, group.by = "celltypist_Immune_All_Low", label = TRUE, raster = FALSE, repel = TRUE) + custom_theme + 
  theme(legend.position = "bottom")

#########################################################################################################################
# Add labels from Immune_All_Low model (majority voting; )
immune_all_low_majority_voting = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/celltypist_outputs/Immune_All_Low/majority_voting_outputs/celltypist_All_Immune_Low_majority_votingpredicted_labels.csv", 
                                       header = TRUE, quote = FALSE)
names(immune_all_low_majority_voting)[1] = "cell_barcode"
head(immune_all_low_majority_voting)

immune_all_low_majority_voting_vec = immune_all_low_majority_voting$majority_voting
names(immune_all_low_majority_voting_vec) = immune_all_low_majority_voting$cell_barcode
head(immune_all_low_majority_voting_vec)

# Add celltypist Immune_All_Low labels (majority voting) into metadata
metadata_df$celltypist_Immune_All_Low_mvoting = immune_all_low_majority_voting_vec[match(rownames(metadata_df),
                                                                                         names(immune_all_low_majority_voting_vec))]
rpca_int_sct_v3@meta.data = metadata_df
DimPlot(rpca_int_sct_v3, group.by = "celltypist_Immune_All_Low_mvoting", label = TRUE, raster = FALSE, repel = TRUE) + custom_theme + 
  theme(legend.position = "bottom")


##########################################################################################
# Add labels from Immune_All_AddPIP model (no majority voting)
immune_all_add_pip = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/celltypist_outputs/Immune_All_AddPIP/celltypist_All_Immune_AddPIP_predicted_labels.csv", 
                           header = TRUE, quote = FALSE)
names(immune_all_add_pip)[1] = "cell_barcode"
head(immune_all_add_pip)

immune_all_add_pip_vec = immune_all_add_pip$predicted_labels
names(immune_all_add_pip_vec) = immune_all_add_pip$cell_barcode

# Add celltypist Immune_All_AddPIP labels into metadata
metadata_df$celltypist_Immune_All_AddPIP = immune_all_add_pip_vec[match(rownames(metadata_df),
                                                                        names(immune_all_add_pip_vec))]
rpca_int_sct_v3@meta.data = metadata_df
DimPlot(rpca_int_sct_v3, group.by = "celltypist_Immune_All_AddPIP", label = TRUE, raster = FALSE, repel = TRUE) + custom_theme + 
  theme(legend.position = "bottom")

############################################################################################
# Add labels from Immune_All_AddPIP model (majority voting)
immune_all_add_pip_mvoting = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/celltypist_outputs/Immune_All_AddPIP/majority_voting_outputs/celltypist_All_Immune_AddPIP_majority_voting_predicted_labels.csv",
                           header = TRUE, quote = FALSE)
names(immune_all_add_pip_mvoting)[1] = "cell_barcode"
head(immune_all_add_pip_mvoting)

immune_all_add_pip_mvoting_vec = immune_all_add_pip_mvoting$majority_voting
names(immune_all_add_pip_mvoting_vec) = immune_all_add_pip_mvoting$cell_barcode
head(immune_all_add_pip_mvoting_vec)

# Add celltypist Immune_All_AddPIP labels (majority voting) into metadata
metadata_df$celltypist_Immune_All_AddPIP_mvoting = immune_all_add_pip_mvoting_vec[match(rownames(metadata_df),
                                                                                            names(immune_all_add_pip_mvoting_vec))]
rpca_int_sct_v3@meta.data = metadata_df
DimPlot(rpca_int_sct_v3, group.by = "celltypist_Immune_All_AddPIP_mvoting", label = TRUE, raster = FALSE, repel = TRUE) + custom_theme + 
  theme(legend.position = "bottom")

############################################
# Extract markers per cluster              #
############################################

####################################
# Find markers per cluster. Need to run PrepSCTFindMarkers first
# IMPORTANT: need to set the proper min.pct and logfc.threshold parameters to identify all relevant markers
DefaultAssay(rpca_int_sct_v3) = "SCT"
rpca_int_sct_v3 = PrepSCTFindMarkers(rpca_int_sct_v3)

# Call DE genes for all clusters. Set min.pct=0.25 and logfc.threshold=0.25. Might have to be
# a bit lenient with the thresholds since there are too many cells. 
cluster_markers = FindAllMarkers(rpca_int_sct_v3,
                                 assay = "SCT",
                                 logfc.threshold = 0.25,
                                 min.pct = 0.25,
                                 only.pos = TRUE)

####################################
# Clean markers df so we can rank genes based on Log2FC
cluster_list = list()
clusters_id = as.character(unique(level2_annotation_markers$cluster))
clusters_n = length(unique(level2_annotation_markers$cluster))
for (i in seq_len(clusters_n)) {
  cluster = level2_annotation_markers %>%
    filter(cluster == clusters_id[i]) %>%
    arrange(desc(avg_log2FC)) %>%
    head(n=100) %>%
    relocate(gene) 
  cluster_list[[i]] = cluster
}

cluster_logfc_ordered_markers_df = data.table::rbindlist(cluster_list)
head(cluster_logfc_ordered_markers_df, n=30)

# Save raw output from FindAllMarkers
saveRDS(cluster_markers, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cluster_markers/rpca_int_sct_v3/cluster_markers_res1_raw.rds")

# Save clean version of markers list with genes arranged in order of descending log2FC
saveRDS(cluster_logfc_ordered_markers_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cluster_markers/rpca_int_sct_v3/cluster_markers_res1_logFC_ordered.rds")
write.table(cluster_logfc_ordered_markers_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cluster_markers/rpca_int_sct_v3/cluster_markers_level1_annotations_logFC_ordered.tsv",
            row.names = FALSE, col.names = TRUE, sep = "\t")


#############################################
# Add cell type labels (level 2 annotations)
pericytes = c(8)
smc = c(0, 2, 3, 14, 31)
fibroblasts = c(11, 12)
myofibroblasts = c(6, 23)
apoe_fibroblasts = c(25)
fibroblast4 = c(33)
neuron = c(35)

# Endo
inflammatory_endo = c(32)
pro_angiogenic_endo = c(4)
intimal_endo = c(20)
lymphatic_endo = c(37)
endomt_endo = c(16)

# Myeloid
foamy_mac1 = c(24)
foamy_mac2 = c(21)
inflammatory_mac = c(1)
monocytes = c(17)
nampt_neutrophils = c(39)
cDCs = c(18) # These could also be dendritic cells, which are crucial in presenting antigens 
tissue_resident_macs = c(29) # Known markers of tissue resident macrophages include F13A, FOLR2, LYVE1, DAB2, STAB1, MRC1, MS4A6A 
monocytes_dc = c(15) # This cluster likely contains monocytes and immature DCs as CD14 and C1Q gene expression is higher than in the cDC cluster. 
phagocytosis_mac = c(27) # Decreased inflammation and increased expression of C1q molecules depicting phagocyte activation; likely M2. 
proliferating_myeloid = c(40)
pDC = c(38)

# Lymphoid cells
recently_activated_nk = c(19) # Higher expression of CD69. Lacks expression of the T cell marker CD3D
cd8_1 = c(7, 10) # These are likely CD8 effector (cytotoxic); There's varying patterns expression of IL7R, CD69, so likely contains effector memory, central memory 
cd8_2 = c(26) # Some cells in this cluster express CD8A and CD8B, as well as the T cell markers CD3E and CD3D. Terminally differentiated cytotoxic CD8 T cell, lacks CD69 expression
cd4_memory_effector_t_cells = c(5, 34) # Based on expression of FOXP3, RORC and GATA3, this cluster contains different T helper types. Also strong expression of IL7R. 
cd4_naive = c(28) # Markers for naive cells include LTB, IL7R, EEF1B2, CCR7, LDHB
t_reg = c(13) # There is subtle but specific expression of FOXP3 and IL2RA in this cluster. 
b_cell = c(22)
plasma_cell = c(36)
mast_cell = c(30)


# Add annotations
new_meta_df = metadata_df %>%
  mutate(level2_annotations = case_when(integrated_snn_res.1 %in% pericytes ~ "Pericyte",
                                        integrated_snn_res.1 %in% smc ~ "SMC",
                                        integrated_snn_res.1 %in% fibroblasts ~ "Fibroblast",
                                        integrated_snn_res.1 %in% myofibroblasts ~ "Myofibroblast",
                                        integrated_snn_res.1 %in% apoe_fibroblasts ~ "APOE_Fibroblast",
                                        integrated_snn_res.1 %in% fibroblast4 ~ "Fibroblast4",
                                        integrated_snn_res.1 %in% neuron ~ "Neuron",
                                        integrated_snn_res.1 %in% inflammatory_endo ~ "Inflammatory_EC",
                                        integrated_snn_res.1 %in% pro_angiogenic_endo ~ "Angiogenic/Vasa_vasorum_EC",
                                        integrated_snn_res.1 %in% intimal_endo ~ "Intimal_EC",
                                        integrated_snn_res.1 %in% lymphatic_endo ~ "Lymphatic_EC",
                                        integrated_snn_res.1 %in% endomt_endo ~ "EndoMT_EC",
                                        integrated_snn_res.1 %in% foamy_mac1 ~ "Foamy_Mac1",
                                        integrated_snn_res.1 %in% foamy_mac2 ~ "Foamy_Mac2",
                                        integrated_snn_res.1 %in% inflammatory_mac ~ "Inflammatory_Mac",
                                        integrated_snn_res.1 %in% monocytes ~ "Monocytes",
                                        integrated_snn_res.1 %in% nampt_neutrophils ~ "NAMPT_Neutrophils",
                                        integrated_snn_res.1 %in% cDCs ~ "cDC",
                                        integrated_snn_res.1 %in% tissue_resident_macs ~ "Tissue_resident_Mac",
                                        integrated_snn_res.1 %in% monocytes_dc ~ "Monocytes/DC",
                                        integrated_snn_res.1 %in% phagocytosis_mac ~ "Phagocytosis_Mac",
                                        integrated_snn_res.1 %in% proliferating_myeloid ~ "Proliferating_myeloid",
                                        integrated_snn_res.1 %in% recently_activated_nk ~ "Activated_NK",
                                        integrated_snn_res.1 %in% cd8_2 ~ "CTL_CD8_terminally_diff",
                                        integrated_snn_res.1 %in% cd8_1 ~ "CTL_CD8_early_activated/mem",
                                        integrated_snn_res.1 %in% cd4_memory_effector_t_cells ~ "CD4_T_effector/mem",
                                        integrated_snn_res.1 %in% t_reg ~ "Treg",
                                        integrated_snn_res.1 %in% cd4_naive ~ "CD4_T_naive",
                                        integrated_snn_res.1 %in% b_cell ~ "B_cell",
                                        integrated_snn_res.1 %in% plasma_cell ~ "Plasma_cell",
                                        integrated_snn_res.1 %in% mast_cell ~ "Mast_cell",
                                        integrated_snn_res.1 %in% pDC ~ "pDC",
                                        integrated_snn_res.1 %in% c(9) & TS_vasc_predicted_id == "smooth muscle cell" ~ "SMC",
                                        integrated_snn_res.1 %in% c(9) & TS_vasc_predicted_id == "fibroblast" ~ "Fibroblast",
                                        integrated_snn_res.1 %in% c(9) & TS_vasc_predicted_id == "pericyte cell" ~ "SMC",
                                        integrated_snn_res.1 %in% c(9) & TS_vasc_predicted_id == "t cell" ~ "SMC"))

rownames(new_meta_df) = rownames(rpca_int_sct_v3@meta.data)
rpca_int_sct_v3@meta.data = new_meta_df


# There are 48806 cells from carotids and 69772 cells from coronary arteries
table(new_meta_df$arterial_origin)

# There are 59691 cells from lesion samples and 58887 cells from non-diseased samples
table(new_meta_df$sample_disease_status)

# Update metadata for integrated object
rpca_int_sct_v3@meta.data = metadata_df
head(rpca_int_sct_v3@meta.data)

# Plot UMAPs according to level 1 and 2 annotations
# Plot level1 annotations 
level1_annotations = DimPlot(rpca_int_sct_v3, reduction = "umap", group.by = "level1_annotations", label = FALSE, repel = TRUE, 
                             pt.size = 0.1, raster = FALSE, label.size = 4) + ggtitle("rPCA Alsaigh/Pan/Wirka/Hu level1 annotations") + custom_theme + 
  theme(legend.position = "none") + npg_scale2 + theme(legend.position = "right")
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure2/Fig2a_level1_annotations_UMAP.pdf",
       plot = level1_annotations, width = 9, height = 9)

# Plot level2 annotations
level2_annotations =  DimPlot(rpca_int_sct_v3, reduction = "umap", group.by = "level2_annotations", label = TRUE, repel = TRUE, 
                              pt.size = 0.1, raster = FALSE, label.size = 4) + ggtitle("rPCA Alsaigh/Pan/Wirka/Hu clusters") + custom_theme + 
  theme(legend.position = "none") + level2_annotations_scale 
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure3/Fig3a_level2_annotations_UMAP_no_labels.pdf",
       plot = level2_annotations, width = 9, height = 9)

# Plot level2 annotations by disease status
level1_by_disease = DimPlot(rpca_int_sct_v3, reduction = "umap", group.by = "level1_annotations", split.by = "sample_disease_status", label = FALSE,
                            repel = TRUE, pt.size = 0.1, raster = FALSE) + npg_scale2 + custom_theme + 
  theme(legend.position = "none",
        strip.text = element_text(size=17))
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure1/Fig1d_level1_UMAP_by_disease_status.svg",
       plot = level1_by_disease, width = 17, height = 8)


#######################
# Plot metadata stats #
#######################

metadata = rpca_int_sct_v3@meta.data


# Plot percentage of cells type annotations by study
level1_study = metadata %>% 
  ggplot(aes(x=study, fill=level1_annotations)) + 
  geom_bar(position = "fill", width = 0.8) +
  custom_theme + 
  ylab("Cell proportion") +
  xlab("Study") + 
  npg_scale2_bars + 
  theme(aspect.ratio = 1.7)

# Get legend
level1_legend = get_legend(level1_study + 
                             theme(legend.box.margin = margin(0, 0, 0, 12)))

# Plot number of cells from lesions and healthy samples
level1_arteries = metadata %>%
  ggplot(aes(x=arterial_origin, fill=level1_annotations)) + 
  geom_bar(position = "fill", width = 0.4) + 
  custom_theme + 
  ylab("Cell proportion") + 
  npg_scale2_bars +
  theme(aspect.ratio = 1.7)

level1_disease = metadata %>%
  ggplot(aes(x=sample_disease_status, fill=level1_annotations)) + 
  geom_bar(position = "fill", width = 0.4) + 
  custom_theme + 
  ylab("Cell proportion") + 
  npg_scale2_bars + 
  theme(aspect.ratio = 1.7)

# Plot metadata figures
prow = cowplot::plot_grid(level1_study + theme(legend.position = "none", 
                                               axis.text.x = element_text(angle = 45, hjust=1)), 
                          level1_arteries + theme(legend.position = "none", 
                                                  axis.text.x = element_text(angle = 45, hjust=1)),
                          level1_disease + theme(legend.position = "none", 
                                                 axis.text.x = element_text(angle = 45, hjust=1)),
                          align = "vh", 
                          hjust = -1, 
                          nrow = 1)
level1_metadata = cowplot::plot_grid(prow, level1_legend, rel_widths = c(3, .4), label_y = "Cell proportion")
ggsave(file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure1/Fig1c_level1_metadata.pdf",
       plot = level1_metadata, width = 14, height = 7)

# Plot number of cells per sequencing library
metadata %>%
  ggplot(aes(x=sample, fill=level1_annotations)) + 
  geom_bar(position = "fill") + 
  custom_theme + 
  ylab("Cell proportion") + 
  npg_scale2_bars + 
  theme(axis.text.x = element_text(angle=45,hjust = 1))

# Plot cell type freq by disease status
level2_disease = metadata %>%
  filter(level1_annotations %in% c("Macrophage", "T_NK")) %>%
  ggplot(aes(x=sample_disease_status, fill=level2_annotations)) + 
  geom_bar(position = "fill", width = 0.4) + 
  custom_theme + 
  ylab("Cell proportion") + 
  level2_annotations_scale_bars + 
  theme(legend.position = "bottom")


# Create vector with cell types of interest
# Normalize cell type freqs to percentage 
cell_types = c("Foamy_Mac1", "Foamy_Mac2", "cDC", "Monocytes", "Monocytes/DC", "B_cell", 
               "Plasma_cell", "pDC", "Mast_cell", "Tissue_resident_Mac")
cell_type_freqs = metadata %>%
  filter(level2_annotations %in% cell_types) %>%
  group_by(level2_annotations, sample_disease_status) %>%
  summarize(n=length(level2_annotations)) %>%
  mutate(percentage = case_when(sample_disease_status == "lesion" ~ n/59691 * 100,
                                TRUE ~ n/58887 * 100))

# Plot cell type percentages 
cell_type_freqs %>%
  ggplot(aes(x=level2_annotations, y=percentage, fill=sample_disease_status)) + 
  geom_col(position = position_dodge(), width = 0.5) + 
  geom_text(aes(label=round(percentage, digits = 2)), position = position_dodge(0.6), hjust=-0.1) + 
  custom_theme + 
  xlab("") + 
  ylab("Cell percentage per disease status") + 
  theme(aspect.ratio = 0.5,
        legend.position = c(0.9, 0.1), 
        axis.text = element_text(size=14)) + 
  coord_flip() +
  scale_fill_lancet()
  
# Subset reference to immune compartments for lesion vs non-lesion comparison
immune_clusters = c("T_NK", "Macrophage", "Plasma_cell", "B_cell", "pDC", "Mast_cell")
rpca_int_sct_v3_immune = subset(rpca_int_sct_v3, subset = level1_annotations %in% immune_clusters)

DimPlot(rpca_int_sct_v3, group.by = "level2_annotations", split.by = "sample_disease_status", raster = FALSE, label = FALSE) + custom_theme + 
  level2_annotations_scale + theme(legend.position = "none")

egg::ggarrange(foam_by_disease, monocytes_by_disease, b_cells_by_disease, plasma_cells_by_disease, pDC_by_disease, neutrophils_by_disease,
               monocytes_dc_by_disease, cdc_by_disease, nrow =  2)

#######################################################
# Remove spurious or artifact clusters for v2 and v3  #
#######################################################

# Remove myeloid/lymphoid cells from whole batch-corrected reference and also cells from cluster 37 since they're likely also artifacts
# Total number of cells before subset: 121316
# Myeloid/lymphoid cells from SMC subclustering: 46
# Total number of cells after subsetting myeloid/lymphoid cells: 121078; REmoved 238 myeloid/lymphoid cells
macs_in_smc_cluster = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/mac_barcodes_in_smc_clusters.rds")
rpca_int_sct_v2 = subset(rpca_int_sct, cells = macs_in_smc_cluster, invert=TRUE)
rpca_int_sct_v2 = subset(rpca_int_sct_v2, idents = c(37), invert=TRUE)
DefaultAssay(rpca_int_sct_v2) = "integrated"
saveRDS(rpca_int_sct_v2, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/alsaigh_pan_wirka_hu_int_seurat_clustered_v2.rds")

# Remove cluster 27 (remaining 119378 cells)
rpca_int_sct_v2 = subset(rpca_int_sct_v2, idents = c(27), invert=TRUE)

# Remove cells clustering along plasma cells (remaining 118884 cells)
plasma_cells = WhichCells(rpca_int_sct_v2, idents = c(32), expression = JCHAIN > 2 | IGLC2 > 2 | MZB1 > 2)
length(plasma_cells) # 431 cells
c32 = WhichCells(rpca_int_sct_v2, idents = c(32)) # 925 cells
cells_to_remove = setdiff(c32, plasma_cells)
rpca_int_sct_v3 = subset(rpca_int_sct_v2, cells = cells_to_remove, invert=TRUE)

# Find and remove T/NK cells in the fibroblast cluster
artifact_t_nk_cells = rpca_int_sct_v3@meta.data %>%
  filter(TS_vasc_predicted_id %in% c("t cell", "nk cell") & integrated_snn_res.1 %in% c(6, 11, 12)) %>%
  rownames()
length(artifact_t_nk_cells)

DimPlot(rpca_int_sct_v3, cells.highlight = artifact_t_nk_cells, pt.size = 0.1, raster = FALSE, sizes.highlight = 0.1) + custom_theme
rpca_int_sct_v3 = subset(rpca_int_sct_v3, cells = artifact_t_nk_cells, invert=TRUE)

# Find and remove B cells from macrophage cluster
artifact_b_cells = rpca_int_sct_v3@meta.data %>%
  filter(TS_vasc_predicted_id %in% c("b cell") & integrated_snn_res.1 %in% c(15, 17, 18, 21, 27)) %>%
  rownames()
length(artifact_b_cells)
DimPlot(rpca_int_sct_v3, cells.highlight = artifact_b_cells, pt.size = 0.1, raster = FALSE, sizes.highlight = 0.1) + custom_theme
rpca_int_sct_v3 = subset(rpca_int_sct_v3, cells = artifact_b_cells, invert=TRUE)

# Find and remove SMCs from T cell cluster
artifact_smcs = rpca_int_sct_v3@meta.data %>%
     filter(TS_vasc_predicted_id %in% c("smooth muscle cell", "pericyte cell") & integrated_snn_res.1 %in% c(28)) %>%
     rownames()
length(artifact_smcs)
rpca_int_sct_v3 = subset(rpca_int_sct_v3, cells = artifact_smcs, invert=TRUE)

####################################################
# Here we'll save the latest updates to the ref v3 #
####################################################
saveRDS(rpca_int_sct_v3, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/seurat_objects/alsaigh_pan_wirka_hu_int_seurat_annotated_v3.rds")


############################################
# Extract markers per cluster or cell type #
############################################

####################################
# Find markers per cluster. Need to run PrepSCTFindMarkers first
# IMPORTANT: need to set the proper min.pct and logfc.threshold parameters to identify all relevant markers
DefaultAssay(rpca_int_sct_v3) = "SCT"
rpca_int_sct_v3 = PrepSCTFindMarkers(rpca_int_sct_v3)

# Call DE genes for level 1 annotations (main cell partitions)
Idents(rpca_int_sct_v3) = "level1_annotations"
level1_annotation_markers = FindAllMarkers(rpca_int_sct_v3,
                                           assay = "SCT",
                                           logfc.threshold = 0.25,
                                           min.pct = 0.25,
                                           only.pos = TRUE)

# Call DE genes for level2 annotations (cell subtypes)
Idents(rpca_int_sct_v3) = "level2_annotations"
level2_annotation_markers = FindAllMarkers(rpca_int_sct_v3,
                                           assay = "SCT",
                                           logfc.threshold = 0.25,
                                           min.pct = 0.25,
                                           only.pos = TRUE)

####################################
# Clean markers df so we can rank genes based on Log2FC
cluster_list = list()
clusters_id = as.character(unique(level2_annotation_markers$cluster))
clusters_n = length(unique(level2_annotation_markers$cluster))
for (i in seq_len(clusters_n)) {
  cluster = level2_annotation_markers %>%
    filter(cluster == clusters_id[i]) %>%
    arrange(desc(avg_log2FC)) %>%
    head(n=100) %>%
    relocate(gene) 
  cluster_list[[i]] = cluster
}

level2_logfc_ordered_markers_df = data.table::rbindlist(cluster_list)
head(cluster_logfc_ordered_markers_df, n=30)


# Save clean version of markers list with genes arranged in order of descending log2FC
#cluster_markers = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cluster_markers/cluster_markers_res1_logFC_ordered.rds")
write.table(cluster_logfc_ordered_markers_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cluster_markers/rpca_int_sct_v3/cluster_markers_level1_annotations_logFC_ordered.tsv",
            row.names = FALSE, col.names = TRUE, sep = "\t")

# Save level 2 annotation markers ordered by Log2FC (complete list)
saveRDS(level2_logfc_ordered_markers_df,
        "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cluster_markers/rpca_int_sct_v3/cluster_markers_level2_annotations_logFC_ordered_full_list.rds")
write.table(level2_logfc_ordered_markers_df,
            "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cluster_markers/rpca_int_sct_v3/cluster_markers_level2_annotations_logFC_ordered_full_list.tsv",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

saveRDS(level2_logfc_ordered_markers_df,
        "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cluster_markers/rpca_int_sct_v3/cluster_markers_level2_annotations_logFC_ordered_top100.rds")

write.table(level2_logfc_ordered_markers_df,
            "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cluster_markers/rpca_int_sct_v3/cluster_markers_level2_annotations_logFC_ordered_top100.tsv",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")


######################################################
# Try different approaches for subclustering of SMCs #
######################################################

# Dims of subset object: 53948 cells
# The following subset is based on a clustering res of 1 for the whole batch-corrected reference
rpca_smc_fibro_subset = subset(rpca_int_sct, idents = c(0, 1, 2, 10, 15, 19, 21, 31))

##########################################################
# APPROACH 1: Subset -> runPCA -> findNeighbors -> runUMAP
# Subcluster SMCs in the most naive way possible 
DefaultAssay(rpca_smc_fibro_subset) = "integrated"

# Run standard workflow for subset seurat obj
rpca_smc_fibro_subset = RunPCA(rpca_smc_fibro_subset)
rpca_smc_fibro_subset = FindNeighbors(rpca_smc_fibro_subset, dims = 1:30, k.param = 20, reduction = "pca")
rpca_smc_fibro_subset = RunUMAP(rpca_smc_fibro_subset, dims = 1:30, n.neighbors = 30, reduction = "pca")

# A res of 0.7 seems appropriate for separating SMCs from pericytes and fibrochondrocytes from fibroblasts
# Best mean sil scores are achieved at a res of 0.7
rpca_smc_fibro_subset = FindClusters(rpca_smc_fibro_subset, resolution = 0.7)


p3_rpca_smc = DimPlot(rpca_smc_fibro_subset, label = TRUE, pt.size = 0.2) + custom_theme
DimPlot(rpca_smc_fibro_subset, group.by = "study")
DimPlot(rpca_smc_fibro_subset, group.by = "sample")
DimPlot(rpca_smc_fibro_subset, split.by = "sample_disease_status")

# How many cells are each study contributing
table(rpca_smc_fibro_subset@meta.data$study)

# Look at markers expression
DefaultAssay(rpca_smc_fibro_subset) = "SCT"
FeaturePlot(rpca_smc_fibro_subset, features = c("MYH11", "VCAN", "CRTAC1", "IBSP", "COL1A2", "VCAM1"), order = TRUE, pt.size = 0.1) & new_scale & custom_theme

# SMC markers
FeaturePlot(rpca_smc_fibro_subset, features = c("MYH11", "CNN1", "ACTA2", "TPM2", "TAGLN", "RGS5"), order = TRUE, pt.size = 0.1) & new_scale & custom_theme

# Fibromyo markers 
FeaturePlot(rpca_smc_fibro_subset, features = c("MYH11", "PLN", "VCAM1", "TNFRSF11B", "KRT17", "IGFBP2"), order = TRUE, pt.size = 0.1) & new_scale & custom_theme

# Fibroblast markers
FeaturePlot(rpca_smc_fibro_subset, features = c("TNFRSF11B", "KRT17", "FN1", "C7", "FBLN1", "SERPINF1"), order = TRUE, pt.size = 0.1) & new_scale & custom_theme



##########################################################
# APPROACH 2: Subset -> SplitBySample -> SCTransform normalize -> runPCA -> FindIntegrationFeatures -> PrepSCTIntegration
DefaultAssay(rpca_smc_fibro_subset) = "RNA"

# It looks like we'll have way too few cells for each Wirka sample so we'll just lump those together
subset_meta = rpca_smc_fibro_subset@meta.data 
subset_meta = subset_meta %>% 
  mutate(sampleID_2 = case_when(study == "wirka_et_al" ~ "wirka_coronary",
                             TRUE ~ sample))
rownames(subset_meta) = rownames(rpca_smc_fibro_subset@meta.data)
table(subset_meta$sampleID_2)
head(subset_meta)

# Update metadata in subset obj (make sure cell barcodes are kept as rownames)
rpca_smc_fibro_subset@meta.data = subset_meta
  
# Split samples in subset obj by "sample2"
rpca_smc_fibro_list = Seurat::SplitObject(rpca_smc_fibro_subset, split.by = "sampleID_2")
class(rpca_smc_fibro_list)
# There are a total of 15 seurat objects
length(rpca_smc_fibro_list)

# Cell cycle score and SCTransform normalize each object (running PCA at this step might not be necessary)
rpca_smc_fibro_list = lapply(rpca_smc_fibro_list, FUN=function(x){
  x = CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes)
  x = SCTransform(x, vst.flavor = "v2", vars.to.regress = c("S.Score","G2M.Score"))})


# Find variable features for integration and prep SCTintegration
subset_int_features = SelectIntegrationFeatures(object.list = rpca_smc_fibro_list, nfeatures = 3000)
saveRDS(subset_int_features, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/pericytes_SMC_fibro_subset_SCT_var_features.rds")
rpca_smc_fibro_list = Seurat::PrepSCTIntegration(object.list = rpca_smc_fibro_list, anchor.features = subset_int_features)
rpca_smc_fibro_list = lapply(rpca_smc_fibro_list, FUN = function(x) {
  x = RunPCA(x, features=subset_int_features)})

# Find integration anchors and integrate data
rpca_sct_smc_fibro_anchors = FindIntegrationAnchors(object.list = rpca_smc_fibro_list, 
                                          normalization.method = "SCT",
                                          anchor.features = subset_int_features, 
                                          dims = 1:30, 
                                          reduction = "rpca",
                                          k.anchor = 5)
rpca_smc_fibro_int_sct = IntegrateData(anchorset = rpca_sct_smc_fibro_anchors, 
                             normalization.method = "SCT",
                             dims=1:30)

# Process integrated data (dim reduction and clustering) like we usually do
DefaultAssay(rpca_smc_fibro_int_sct) = "integrated"
rpca_smc_fibro_int_sct = RunPCA(rpca_smc_fibro_int_sct) # 30 dims captures really most of the variation in the data
rpca_smc_fibro_int_sct = FindNeighbors(rpca_smc_fibro_int_sct, reduction = "pca", dims = 1:30, k.param = 20)
# In general, this parameter should be 5-50. Higher values result in more global structure being preserved at the loss of detailed local structure
# UMAP looks better with 30 n.neighbors than 20
rpca_smc_fibro_int_sct = RunUMAP(rpca_smc_fibro_int_sct, dims = 1:30, n.neighbors = 30, min.dist = 0.3) # Adjusting min.dist doesn't really help

# A resolution of 1 seems appropriate for evaluating different methods for now
# Haven't been able to calulate silhouette scores for this dataset since R keeps crashing. 
# A res of 1 is probably good for now, as long as it distinguishes SMCs from pericytes and fibrochondrocytes from fibroblasts clearly. 
rpca_smc_fibro_int_sct = FindClusters(rpca_smc_fibro_int_sct, resolution = 0.7)

# Visualize clusters 
# Visualization of clusters
p1_smc_rpca = DimPlot(rpca_smc_fibro_int_sct, reduction = "umap", group.by = "sample", pt.size = 0.2, raster = FALSE) + ggtitle("rPCA Alsaigh/Pan/Wirka/Hu samples") + custom_theme
p2_smc_rpca = DimPlot(rpca_smc_fibro_int_sct, reduction = "umap", group.by = "study", pt.size = 0.2, raster = FALSE) + ggtitle("rPCA Alsaigh/Pan/Wirka/Hu study") + custom_theme
p3_smc_rpca = DimPlot(rpca_smc_fibro_int_sct, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.2, raster = FALSE) + ggtitle("rPCA Alsaigh/Pan/Wirka/Hu clusters") + custom_theme
p1_rpca + p2_rpca

# Look at markers expression
DefaultAssay(rpca_smc_fibro_int_sct) = "SCT"
FeaturePlot(rpca_smc_fibro_int_sct, features = c("MYH11", "VCAN", "CRTAC1", "IBSP", "COL1A2", "VCAM1"), order = TRUE, pt.size = 0.3) & new_scale

# SMC markers: CNN1, SMTN, TPM2, RGS5, MYOCD might be good candidates to show in the paper
FeaturePlot(rpca_smc_fibro_int_sct, features = c("MYH11", "CNN1", "ACTA2", "TPM2", "TAGLN", "RGS5"), order = TRUE, pt.size = 0.3) & new_scale

# Fibromyo markers 
FeaturePlot(rpca_smc_fibro_int_sct, features = c("CNN1", "VCAM1", "FN1", "TNFRSF11B", "KRT17", "VCAN"), order = TRUE, pt.size = 0.4) & new_scale

# Osteochondrogenic markers
FeaturePlot(rpca_smc_fibro_int_sct, features = c("CNN1", "COL1A2", "CRTAC1", "COMP", "IBSP", "COL3A1"), order = TRUE, pt.size = 0.4) & new_scale

# Fibroblast markers
FeaturePlot(rpca_smc_fibro_int_sct, features = c("TNFRSF11B", "KRT17", "IGFBP2", "C7", "FBLN1", "SERPINF1"), order = TRUE, pt.size = 0.3) & new_scale

# Pericyte markers
FeaturePlot(rpca_smc_fibro_int_sct, features = c("CNN1", "TPM2", "FABP4", "STEAP4", "RERGL", "NET1"), order = TRUE, pt.size = 0.4) & new_scale



##########################################################
# APPROACH 3: Subset -> SplitBySample -> SCTransform normalize -> runPCA on integrated assay -> Cluster and runUMAP
DefaultAssay(rpca_smc_fibro_subset) = "RNA"

# It looks like we'll have way too few cells for each Wirka sample so we'll just lump those together
subset_meta = rpca_smc_fibro_subset@meta.data 
subset_meta = subset_meta %>% 
  mutate(sampleID_2 = case_when(study == "wirka_et_al" ~ "wirka_coronary",
                                TRUE ~ sample))
rownames(subset_meta) = rownames(rpca_smc_fibro_subset@meta.data)
table(subset_meta$sampleID_2)
head(subset_meta)

# Update metadata in subset obj (make sure cell barcodes are kept as rownames)
rpca_smc_fibro_subset@meta.data = subset_meta

# Split samples in subset obj by "sample2"
rpca_smc_fibro_list = Seurat::SplitObject(rpca_smc_fibro_subset, split.by = "sampleID_2")
class(rpca_smc_fibro_list)
# There are a total of 15 seurat objects
length(rpca_smc_fibro_list)

# Cell cycle score and SCTransform normalize each object (running PCA at this step might not be necessary)
# This is done only for getting a new and more refined set of variable features
rpca_smc_fibro_list = lapply(rpca_smc_fibro_list, FUN=function(x){
  x = CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes)
  x = SCTransform(x, vst.flavor = "v2", vars.to.regress = c("S.Score","G2M.Score"))})


# Find variable features for running PCA on integrated assay of the subset
subset_int_features = SelectIntegrationFeatures(object.list = rpca_smc_fibro_list, nfeatures = 3000)

# Process integrated data (dim reduction and clustering) like we usually do
DefaultAssay(rpca_smc_fibro_subset) = "integrated"
rpca_smc_fibro_subset = RunPCA(rpca_smc_fibro_subset, features = subset_int_features) # 30 dims captures really most of the variation in the data
rpca_smc_fibro_subset = FindNeighbors(rpca_smc_fibro_subset, reduction = "pca", dims = 1:30, k.param = 20)
# In general, this parameter should be 5-50. Higher values result in more global structure being preserved at the loss of detailed local structure
# UMAP looks better with 30 n.neighbors than 20
rpca_smc_fibro_subset = RunUMAP(rpca_smc_fibro_subset, dims = 1:30, n.neighbors = 30, min.dist = 0.3) # Adjusting min.dist doesn't really help

# A resolution of 1 seems appropriate for evaluating different methods for now
# Haven't been able to calulate silhouette scores for this dataset since R keeps crashing. 
# A res of 1 is probably good for now, as long as it distinguishes SMCs from pericytes and fibrochondrocytes from fibroblasts clearly. 
rpca_smc_fibro_subset = FindClusters(rpca_smc_fibro_subset, resolution = 0.7)

# Visualize clusters 
# Visualization of clusters
p1_smc_rpca = DimPlot(rpca_smc_fibro_subset, reduction = "umap", group.by = "sample", pt.size = 0.2, raster = FALSE) + ggtitle("rPCA Alsaigh/Pan/Wirka/Hu samples") + custom_theme
p2_smc_rpca = DimPlot(rpca_smc_fibro_subset, reduction = "umap", group.by = "study", pt.size = 0.2, raster = FALSE) + ggtitle("rPCA Alsaigh/Pan/Wirka/Hu study") + custom_theme
p3_smc_rpca = DimPlot(rpca_smc_fibro_subset, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.2, raster = FALSE) + ggtitle("rPCA Alsaigh/Pan/Wirka/Hu clusters") + custom_theme
p1_rpca + p2_rpca

# Look at markers expression
DefaultAssay(rpca_smc_fibro_subset) = "SCT"
FeaturePlot(rpca_smc_fibro_subset, features = c("MYH11", "VCAN", "CRTAC1", "IBSP", "COL1A2", "VCAM1"), order = TRUE, pt.size = 0.3) & new_scale

# SMC markers: CNN1, SMTN, TPM2, RGS5, MYOCD might be good candidates to show in the paper
FeaturePlot(rpca_smc_fibro_subset, features = c("MYH11", "CNN1", "ACTA2", "TPM2", "TAGLN", "RGS5"), order = TRUE, pt.size = 0.3) & new_scale

# Fibromyo markers
# SMC markers: CNN1, SMTN, TPM2, RGS5, MYOCD might be good candidates to show in the paper
FeaturePlot(rpca_smc_fibro_subset, features = c("TNFRSF11B", "KRT17", "IGFBP2", "C7", "FBLN1", "SERPINF1"), order = TRUE, pt.size = 0.3) & new_scale

# Osteochondrogenic markers
FeaturePlot(rpca_smc_fibro_subset, features = c("CNN1", "COL1A2", "CRTAC1", "COMP", "IBSP", "COL3A1"), order = TRUE, pt.size = 0.3) & new_scale








