library(Seurat)
library(scRNAutils)
library(tidyverse)
library(data.table)
library(cluster)
library(RColorBrewer)
library(parallel)
library(reticulate)
library(UCell)
library(gprofiler2)

# Set seed for reproducibility
set.seed(1)

# Source our own scRNA analysis utils functions
#source("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/R_scripts/Utils/scRNA_processing_utils.R")

################################################################################################
# This script contains code for subclustering of SMCs, pericytes and fibroblasts, in addition  #
# to generation of cluster markers, enrichment of murine gene modules metadata plots           #
# and GSEA for SMC clusters.                                                                   #
################################################################################################

# Load batch corrected reference with cell type annotations. 
# This Seurat object contains the subclustered and annotated data
rpca_smc_fibro_subset_v3 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/seurat_objects/whole_ref_v3/peri_smc_fibro_subset_seurat_obj_res0.9_v2.rds")
DefaultAssay(rpca_smc_fibro_subset_v3) = "SCT"

#########################
# Subclustering of SMCs #
#########################

# Dims of subset object: 37032 cells
# The following subset is based on a clustering res of 1 for the whole batch-corrected reference
# Cluster 37 was excluded from the main reference because it likely represents myeloid/lymphoid cells that got mixed within the main SMC clusters
rpca_smc_fibro_subset_v3 = subset(rpca_int_sct_v3, idents = c(0, 2, 3, 8, 9, 14, 31))


##########################################################
# APPROACH 1: Subset -> runPCA -> findNeighbors -> runUMAP
# Subcluster SMCs in the most naive way possible 
DefaultAssay(rpca_smc_fibro_subset_v3) = "integrated"

# Run standard workflow for subset seurat obj
rpca_smc_fibro_subset_v3 = RunPCA(rpca_smc_fibro_subset_v3)
rpca_smc_fibro_subset_v3 = FindNeighbors(rpca_smc_fibro_subset_v3, dims = 1:30, k.param = 20, reduction = "pca")
rpca_smc_fibro_subset_v3 = RunUMAP(rpca_smc_fibro_subset_v3, dims = 1:30, n.neighbors = 30, reduction = "pca")

# A res of 0.7 seems appropriate for separating SMCs from pericytes and fibrochondrocytes from fibroblasts
# Best mean sil scores are achieved at a res of 0.7
rpca_smc_fibro_subset_v3 = FindClusters(rpca_smc_fibro_subset_v3, resolution = 0.9)

# This is where the latest changes to the SMC subckustered seurat obj are stored
saveRDS(rpca_smc_fibro_subset_v3, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/seurat_objects/whole_ref_v3/peri_smc_fibro_subset_seurat_obj_res0.9_v2.rds")

# Visualize clusters at res0.7 and cell type annotations
p3_rpca_smc = DimPlot(rpca_smc_fibro_subset_v3, label = TRUE, pt.size = 0.2) + custom_theme
prelim_annotations = DimPlot(rpca_smc_fibro_subset_v3, group.by = "prelim_annotations", label = FALSE, repel = TRUE, pt.size = 0.2) + custom_theme + npg_scale +
  theme(legend.position = "bottom")
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure4/Fig4b_SMC_prelim_annotations_UMAP.pdf",
       plot = prelim_annotations, width = 9, height = 9)
DimPlot(rpca_smc_fibro_subset_v3, group.by = "prelim_annotations", split.by = "study") & custom_theme & npg_scale
DimPlot(rpca_smc_fibro_subset_v3, group.by = "sample") & custom_theme & npg_scale 
DimPlot(rpca_smc_fibro_subset_v3, split.by = "sample_disease_status")

# How many cells are each study contributing
# Alsaigh et al: 8734; Hu et al: 36073, Pan et al: 4277; wirka et al: 4864
table(rpca_smc_fibro_subset@meta.data$study)

# Look at markers expression
DefaultAssay(rpca_smc_fibro_subset_v3) = "SCT"
FeaturePlot(rpca_smc_fibro_subset_v3, features = c("MYH11", "VCAN", "COL6A3", "FN1", "COL1A2", "VCAM1"), 
            order = TRUE, pt.size = 0.1) & custom_theme()

# SMC markers
FeaturePlot(rpca_smc_fibro_subset_v3, features = c("MYH11", "CNN1", "ACTA2", "TPM2", "TAGLN", "LMOD1"), 
            order = TRUE, pt.size = 0.1) & custom_theme()

# Fibromyo markers 
FeaturePlot(rpca_smc_fibro_subset_v3, features = c("MYH11", "VCAM1", "TNFRSF11B", "KRT17", "VCAN", "LTBP1"), 
            order = TRUE, pt.size = 0.1) & custom_theme()

# Fibrochondro markers 
FeaturePlot(rpca_smc_fibro_subset_v3, features = c("CNN1", "FN1", "COL1A2", "COL3A1", "CRTAC1", "CLU"), 
            order = TRUE, pt.size = 0.1) & custom_theme()

# Fibroblast markers
FeaturePlot(rpca_smc_fibro_subset_v3, features = c("TNFRSF11B", "KRT17", "FN1", "C7", "FBLN1", "SERPINF1"), 
            order = TRUE, pt.size = 0.1) & custom_theme()

FeaturePlot(rpca_smc_fibro_subset_v3, split.by = "sample_disease_status", features = c("HDAC4", "HDAC5", "PTK2"), raster = FALSE)

# CRTAC1
crtac1 = FeaturePlot(rpca_smc_fibro_subset_v3, features = c("CRTAC1"), pt.size = 0.7) & custom_theme() & miller_continous_scale()
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure6/Fig6a_CRTAC1_UMAP.pdf",
       plot = crtac1, width = 7, height = 7)

# IBSP
ibsp = FeaturePlot(rpca_smc_fibro_subset_v3, features = c("IBSP"), pt.size = 0.7, 
                   order = TRUE) & custom_theme() & miller_continous_scale()
#ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure6/Fig6a_IBSP_UMAP.pdf",
#       plot = ibsp, width = 7, height = 7)

# LTBP1
ltbp1 = FeaturePlot(rpca_smc_fibro_subset_v3, features = c("LTBP1"), pt.size = 0.7, 
                    order = TRUE) & custom_theme() & miller_continous_scale()
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure6/Fig6a_LTBP1_UMAP.pdf",
       plot = ltbp1, width = 7, height = 7)

# VCAN
vcan = FeaturePlot(rpca_smc_fibro_subset_v3, features = c("VCAN"), pt.size = 0.7, 
                   order=TRUE) & custom_theme() & miller_continous_scale()
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure6/Fig6a_VCAN_UMAP.pdf",
       plot = vcan, width = 7, height = 7)

# VCAM1
vcam1 = FeaturePlot(rpca_smc_fibro_subset_v3, features = c("VCAM1"), pt.size = 0.5, 
                    order=TRUE) & custom_theme() & miller_continous_scale()
ggsave("/gpfs/gpfs0/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Rebuttal_Fig1/Reb_Fig1b_VCAM1_UMAP.png",
       plot = vcam1, width = 7, height = 7)

# Fibroblast markers
fibroblast_markers = FeaturePlot(rpca_smc_fibro_subset_v3, features = c("FBLN1", "C7", "SERPINF1", "ACTA2"), pt.size = 0.1, 
            order=TRUE) & custom_theme() & miller_continous_scale()

ggsave("/gpfs/gpfs0/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Rebuttal_Fig4/Reb_Fig4b_fibroblast_markers_UMAP.png",
       plot = fibroblast_markers, width = 10, height = 10)


# Make a dotplot for prelim cell type annotations
smc_marker_dotplot = DotPlot(rpca_smc_fibro_subset_v3, group.by = "prelim_annotations",  features = c("MYH11", "CNN1", "TPM2", "MYL9", "ACTA2", 
                                                                             "VCAM1", "CTGF","TNFRSF11B", "IGFBP2", "BGN", "FN1", "VCAN", "LTBP1", "COL4A2", "COL6A2",
                                                                             "COL1A1", "COL1A2", "COL3A1", "CLU", "CRTAC1", "COMP", "FBLN1", "C7", "SERPINF1", "APOC1", "APOE", 
                                                                             "HLA-A", "HLA-B", "TIMP1", "AGT", "IGFBP3", "CCL19", "CARMN", "MYH10", "NEAT1", "AHNAK",
                                                                             "RERGL", "NET1", "FABP4", "ID4"), col.min = 0) +
   new_scale3 + custom_theme + theme(axis.text.x = element_text(angle = 45, hjust=1), axis.text = element_text(size=12)) + coord_flip()

ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure4/Fig4c_SMC_markers_dot_plot.pdf",
       plot = smc_marker_dotplot, width = 12, height = 10)


###########################################
# Score cell types signatures using UCell #
###########################################

# Create a list to store murine gene modules that will be tested for enrichment in human data
markers = list()
markers$mouse_meta_SMC1 = smc1_human_homologs_vec_top100
markers$mouse_meta_SEM = sem_human_homologs_vec_top50
markers$mouse_meta_Fibrochondrocyte = fc_human_homologs_vec_top50
markers$mouse_meta_Fibroblast1 = fibro1_human_homologs_vec_top50

rpca_smc_fibro_subset =  AddModuleScore_UCell(rpca_smc_fibro_subset_v3, features = markers)
signature.names = paste0(names(markers), "_UCell")

FeaturePlot(rpca_smc_fibro_subset, features = signature.names, 
            order = TRUE, pt.size = 0.3) & miller_continous_scale() & custom_theme()

# Create a list to store murine gene modules from Kim et al to enrich in human data
kim_markers = list()
kim_markers$Kim_Contractile_SMC = smc_kim_homologs_vec
kim_markers$Kim_FMC = fmc_kim_homologs_vec
kim_markers$Kim_FC = fc_kim_homologs_vec

rpca_smc_fibro_subset = AddModuleScore_UCell(rpca_smc_fibro_subset_v3, features = kim_markers)
signature.names = paste0(names(kim_markers), "_UCell")

kim_signatures = FeaturePlot(rpca_smc_fibro_subset, features = signature.names, 
                             order = TRUE, pt.size = 0.1) & miller_continous_scale() & custom_theme()

ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/rebuttal_figures/Reb_Fig1a_Kim_signatures.png",
       plot = kim_signatures, width = 8, height = 8)


#########################################
# Cluster 14 looks like a Mac artifact so it'll be removed to avoid issues during marker generation.
# Cells were reclustered using the above code after removing this artifact
macs_in_smc_cluster = WhichCells(rpca_smc_fibro_subset, idents = c(14))
rpca_smc_fibro_subset_v3 = subset(rpca_smc_fibro_subset, idents = c(14), invert=TRUE) 

##########################
# Find markers per cluster
# Currently using a res of 0.9 
# IMPORTANT: need to set the proper min.pct and logfc.threshold parameters to identify all relevant markers
DefaultAssay(rpca_smc_fibro_subset_v3) = "SCT"
rpca_smc_fibro_subset_v3 = PrepSCTFindMarkers(rpca_smc_fibro_subset_v3)

# Call DE genes for all clusters. FOr now start with default settings and adjust accordingly
# Perhaps a logfc.threshold of 0.4 is too stringent, might be good enough to set min.pct=0.1 and logfc.treshold=0.25 for when we get markers
# for annotated cell types
peri_smc_fibro_cluster_markers = FindAllMarkers(rpca_smc_fibro_subset,
                                                assay = "SCT",
                                                logfc.threshold = 0.4,
                                                min.pct = 0.1,
                                                only.pos = TRUE)
saveRDS(peri_smc_fibro_cluster_markers, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/peri_smc_fibro_markers_res0.7.rds")
peri_smc_fibro_cluster_markers = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/peri_smc_fibro_markers_res0.7.rds")
head(peri_smc_fibro_cluster_markers)

# Call DE genes for preliminary annotations. Set min.pct=0.1 and logfc.treshold=0.25
# First we need to define the cell type annotations as default identities for clusters
Idents(rpca_smc_fibro_subset_v3) = "prelim_annotations"
peri_smc_fibro_annotation_markers = FindAllMarkers(rpca_smc_fibro_subset_v3, 
                                              assay = "SCT",
                                              logfc.threshold = 0.25,
                                              min.pct = 0.1,
                                              only.pos = TRUE)


####################################
# Clean markers df so we can ranked genes based on Log2FC
cluster_list = list()
clusters_id = as.character(unique(peri_smc_fibro_annotation_markers$cluster))
clusters_n = length(unique(peri_smc_fibro_annotation_markers$cluster))
for (i in seq_len(clusters_n)) {
  cluster = peri_smc_fibro_annotation_markers %>%
    filter(cluster == clusters_id[i]) %>%
    arrange(desc(avg_log2FC)) %>%
    relocate(gene) %>%
    head(n=100)
  cluster_list[[i]] = cluster
}
cluster_top_markers_df = data.table::rbindlist(cluster_list)
head(cluster_top_markers_df)

saveRDS(peri_smc_fibro_annotation_markers, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/cell_type_annotations_markers/peri_smc_fibro_prelim_anno_raw_markers.rds")
saveRDS(cluster_top_markers_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/cell_type_annotations_markers/peri_smc_fibro_prelim_anno_top100_markers.rds")
write.table(cluster_top_markers_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/cell_type_annotations_markers/peri_smc_fibro_prelim_anno_top100_markers.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Load cell type markers from current annotations
cell_type_markers = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/cell_type_annotations_markers/peri_smc_fibro_prelim_anno_top100_markers.rds")


###################################################################
# Add preliminary annotations for clusters based on cluster markers (Currently the second set of prelimnary annotations)
smc1 = c(2, 5, 16)
smc2 = c(12)
smc3 = c(3)
transitional = c(0, 4, 8, 9)
pericyte1 = c(10, 14)
pericyte2 = c(7, 11)
foam_cell = c(15)
fibromyocyte = c(6, 17)
fibroblast = c(1)
osteochondrogenic = c(13)

barcodes = rownames(rpca_smc_fibro_subset_v3@meta.data)
rpca_smc_fibro_subset_v3@meta.data = rpca_smc_fibro_subset_v3@meta.data %>%
  mutate(prelim_annotations = case_when(seurat_clusters %in% smc1 ~ "Contractile_SMC",
                   seurat_clusters %in% smc2 ~ "SMC2",
                   seurat_clusters %in% smc3 ~ "SMC3",
                   seurat_clusters %in% transitional ~ "Transitional-ECM-SMC",
                   seurat_clusters %in% pericyte1 ~ "Pericyte1",
                   seurat_clusters %in% pericyte2 ~ "Pericyte2",
                   seurat_clusters %in% foam_cell ~ "Foam-like",
                   seurat_clusters %in% fibromyocyte ~ "Fibromyocyte",
                   seurat_clusters %in% fibroblast ~ "Fibroblast",
                   seurat_clusters %in% osteochondrogenic ~ "Fibrochondrocyte"))
rownames(rpca_smc_fibro_subset_v3@meta.data) = barcodes
rpca_smc_fibro_subset_v3@meta.data$prelim_annotations = factor(rpca_smc_fibro_subset_v3@meta.data$prelim_annotations,
                                                               levels = c("Contractile_SMC", "Transitional-ECM-SMC", "Fibromyocyte",
                                                                          "Fibrochondrocyte", "Fibroblast", "SMC2", "SMC3", "Foam-like", "Pericyte1", "Pericyte2"))
                                                  



# Call DE genes for preliminary annotations. Set min.pct=0.1 and logfc.treshold=0.25 
# First we need to define cell type annotations as the default identities for each cluster
Idents(rpca_smc_fibro_subset_v3) = "prelim_annotations"
peri_smc_fibro_annotation_markers = FindAllMarkers(rpca_smc_fibro_subset_v3, 
                                                   assay = "SCT",
                                                   logfc.threshold = 0.25,
                                                   min.pct = 0.1,
                                                   only.pos = TRUE)


# Clean markers df so we can ranked genes based on Log2FC
cluster_list = list()
clusters_id = as.character(unique(peri_smc_fibro_annotation_markers$cluster))
clusters_n = length(unique(peri_smc_fibro_annotation_markers$cluster))
for (i in seq_len(clusters_n)) {
  cluster = peri_smc_fibro_annotation_markers %>%
    filter(cluster == clusters_id[i]) %>%
    arrange(desc(avg_log2FC)) %>%
    relocate(gene) %>%
    head(n=100)
  cluster_list[[i]] = cluster
}
cluster_top_markers_df = data.table::rbindlist(cluster_list)
head(cluster_top_markers_df)

saveRDS(peri_smc_fibro_annotation_markers, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/cell_type_annotations_markers/peri_smc_fibro_prelim_anno_raw_markers.rds")
saveRDS(cluster_top_markers_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/cell_type_annotations_markers/peri_smc_fibro_prelim_anno_top100_markers.rds")
write.table(cluster_top_markers_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/cell_type_annotations_markers/peri_smc_fibro_prelim_anno_top100_markers.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Load cell type markers from current annotations
cell_type_markers = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/cell_type_annotations_markers/peri_smc_fibro_prelim_anno_top100_markers.rds")
          

+###################################################
# Add metadata for sex

# Sex was defined based on metadata from papers and also XIST expression. 
males = c("pan_rpe004", "pan_rpe005", "alsaigh_ac_p1", "alsaigh_pa_p1", "alsaigh_ac_p2", "alsaigh_pa_p2",
          "alsaigh_ac_p3", "alsaigh_pa_p3", "wirka_coronary_1", "wirka_coronary_2", "wirka_coronary_3",
          "wirka_coronary_4", "wirka_coronary_5", "wirka_coronary_6", "wirka_coronary_7", "hu_coronary1_p1", 
          "hu_coronary1_p2", "hu_coronary2_p2")
females = c("pan_rpe006", "wirka_coronary_8", "hu_coronary1_p3", "hu_coronary2_p3")

  
rpca_smc_fibro_subset_v3@meta.data = rpca_smc_fibro_subset_v3@meta.data %>%
  mutate(sex = case_when(sample %in% males ~ "males",
                         sample %in% females ~ "females"))
rownames(rpca_smc_fibro_subset_v3@meta.data) = barcodes

  
# Check UMAP embeddings according to sex 
DimPlot(rpca_smc_fibro_subset_v3, group.by = "prelim_annotations", split.by = "sex") & custom_theme & npg_scale
    
#######################
# Make metadata plots #
#######################

# Plot cell type proportions
cell_proportions = rpca_smc_fibro_subset_v3@meta.data %>%
  ggplot(aes(x=study, fill=prelim_annotations)) + 
  geom_bar(width = 0.7, position = "fill") +
  theme_bw() + 
  ylab("Cell proportion") +
  xlab("") + 
  custom_theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  npg_scale_bars

cell_proportions

# Plot cell type proportions by disease status
cell_proportions_by_disease = rpca_smc_fibro_subset_v3@meta.data %>%
  ggplot(aes(x=sample_disease_status, fill=prelim_annotations)) + 
  geom_bar(width = 0.5, position = "fill") +
  theme_bw() + 
  ylab("Cell proportion") +
  xlab("") + 
  custom_theme +
  theme(aspect.ratio = 1.5,
        axis.text.x = element_text(angle = 45, hjust=1)) + 
  npg_scale_bars

# Plot cell type proportions by arterial bed
cell_proportions_by_artery = rpca_smc_fibro_subset_v3@meta.data %>%
     ggplot(aes(x=arterial_origin, fill=prelim_annotations)) + 
     geom_bar(width = 0.5, position = "fill") +
     theme_bw() + 
     ylab("Cell proportion") +
     xlab("") + 
     custom_theme +
     theme(aspect.ratio = 1.5,
           axis.text.x = element_text(angle = 45, hjust=1)) + 
     npg_scale_bars

# Get legend
cell_type_legend = get_legend(cell_proportions + 
                             theme(legend.box.margin = margin(0, 0, 0, 12)))

# Plot metadata figures
prow = cowplot::plot_grid(cell_proportions + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1)), 
                          cell_proportions_by_artery + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1)),
                          cell_proportions_by_disease + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1)),
                          align = "vh", 
                          hjust = -1, 
                          nrow = 1)
smc_metadata_plots = cowplot::plot_grid(prow, cell_type_legend, rel_widths = c(4, .7), label_y = "Cell proportion")
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure3/SuppFig3b_SMC_metadata.pdf",
       plot = smc_metadata_plots, width = 11, height = 4.5)



# This is the subset of pericytes, SMCs and fibroblasts. Cluster 14 was removed since it contained myeloid cells with no contractility or ECM related genes expression
saveRDS(macs_in_smc_cluster, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/mac_barcodes_in_smc_clusters.rds")
saveRDS(rpca_smc_fibro_subset_v3, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/peri_smc_fibro_subset_seurat_obj.rds")


############################################
# Make gene set enrichment plots           #
############################################

# Load SEM object ro remove GO terms redundancy if needed. 
hsGO = GOSemSim::godata("org.Hs.eg.db", ont="BP")

############################################
# Fibrochondrogenic cells
fibrochondro_markers = cell_type_markers %>%
  filter(cluster == "Fibrochondrocyte")
fibrochondro_ranks = fibrochondro_markers$avg_log2FC
names(fibrochondro_ranks) = fibrochondro_markers$gene
length(fibrochondro_ranks)

fibrochondro_ranks_top100 = fibrochondro_ranks[1:100]

# Run gene set enrichment
fibrochondro_gost_res = gprofiler2::gost(names(fibrochondro_ranks_top100), organism = "hsapiens",
                                         ordered_query = TRUE, correction_method = "fdr")
fibrochondro_gprofiler_res = fibrochondro_gost_res$result
fibrochondro_go_bp = fibrochondro_gprofiler_res[fibrochondro_gprofiler_res$source == "GO:BP", ]

# Reduce redundancy in GO BP terms
redundant_terms = gprofiler_go_simplify(fibrochondro_go_bp, semData=hsGO, ontology= "BP")
fibrochondro_go_bp_non_redundant = fibrochondro_go_bp[-which(fibrochondro_go_bp$term_id %in% redundant_terms),]

fibrochondro_go_plot = gprofiler_bar_plot(fibrochondro_go_bp_non_redundant, 15, 500, "GO:BP", 10) + new_scale3_bars

##############################################  
# Fibromyocytes
fibromyo_markers = cell_type_markers %>%
  filter(cluster == "Fibromyocyte")
fibromyo_ranks = fibromyo_markers$avg_log2FC
names(fibromyo_ranks) = fibromyo_markers$gene
length(fibromyo_ranks)
  
fibromyo_ranks_top100 = fibromyo_ranks[1:100]

# Run gene set enrichment  
fibromyo_gost_res = gprofiler2::gost(names(fibromyo_ranks_top100), organism = "hsapiens",
                                           ordered_query = TRUE, correction_method = "fdr")  
fibromyo_grofiler_res = fibromyo_gost_res$result
fibromyo_go_bp = fibromyo_grofiler_res[fibromyo_grofiler_res$source == "GO:BP", ]


gprofiler_bar_plot(fibromyo_go_bp, 15, 500, "GO:BP", 10) + new_scale3_bars

##############################################  
# Transitional ECM SMC
pioneer_markers = cell_type_markers %>%
  filter(cluster == "Transitional-ECM-SMC")
pioneer_ranks = pioneer_markers$avg_log2FC
names(pioneer_ranks) = pioneer_markers$gene
length(pioneer_ranks)

pioneer_ranks_top100 = pioneer_ranks[1:100]

# Run gene set enrichment  
pioneer_gost_res = gprofiler2::gost(names(pioneer_ranks_top100), organism = "hsapiens",
                                     ordered_query = TRUE, correction_method = "fdr")  
pioneer_grofiler_res = pioneer_gost_res$result
pioneer_go_bp = pioneer_grofiler_res[pioneer_grofiler_res$source == "GO:BP", ]


gprofiler_bar_plot(pioneer_go_bp, 15, 500, "GO:BP", 10) + new_scale3_bars


################################################
# Contractile SMCs
smc_markers = cell_type_markers %>%
  filter(cluster == "Contractile_SMC")
smc_ranks = smc_markers$avg_log2FC
names(smc_ranks) = smc_markers$gene
length(smc_ranks)

smc_ranks_top100 = smc_ranks[1:100]

# Run gene set enrichment
smc_gost_res = gprofiler2::gost(names(smc_ranks_top100), organism = "hsapiens",
                                     ordered_query = TRUE, correction_method = "fdr")  
smc_grofiler_res = smc_gost_res$result
smc_go_bp = smc_grofiler_res[smc_grofiler_res$source == "GO:BP", ]

smc_go_plot = gprofiler_bar_plot(smc_go_bp, 15, 500, "GO:BP", 10) + new_scale3_bars


##################################################
# SMC2
smc2_markers = cell_type_markers %>%
  filter(cluster == "SMC2")
smc2_ranks = smc2_markers$avg_log2FC
names(smc2_ranks) = smc2_markers$gene
length(smc2_ranks)

smc2_ranks_top100 = smc2_ranks[1:100]

smc2_gost_res = gprofiler2::gost(names(smc2_ranks_top100), organism = "hsapiens",
                                 ordered_query = TRUE, correction_method = "fdr")  

smc2_bar_plot = gprofiler_bar_plot(smc2_gost_res, 10, 500, "GO:BP", 15)

#######################################################################
# Foam cells
foam_like_markers = cell_type_markers %>%
     filter(cluster == "Foam-like")
foam_like_ranks = foam_like_markers$avg_log2FC
names(foam_like_ranks) = foam_like_markers$gene
length(foam_like_ranks)

foam_like_ranks_top100 = foam_like_ranks[1:100]

# Run gene set enrichment
foam_like_gost_res = gprofiler2::gost(names(foam_like_ranks_top100), organism = "hsapiens", 
                                      ordered_query = TRUE, correction_method = "fdr")  
foam_like_gprofiler_res = foam_like_gost_res$result
foam_like_go_bp = foam_like_gprofiler_res[foam_like_gprofiler_res$source == "GO:BP", ]
 
foam_like_smc_go_plot = gprofiler_bar_plot(foam_like_go_bp, 15, 500, "GO:BP", 10)  


#################################################################
# Try a heatmap for SMC annotation GO terms

# We need to write a function that will prep GO BP output dataframes 
# for making the heatmap
# Args go_df A dataframe output by gprofiler
# Args annotation A character denoting the name of the cell type annotation
# value A df where GO terms have been filtered and ordered by p values 
prep_terms = function(go_df, annotation) { 
  filtered_terms = go_df %>%
    filter(term_size > 15 & term_size < 500) %>%
    mutate(log10_pval = -log10(p_value)) %>%
    mutate(annotation = annotation) %>%
    arrange(desc(log10_pval)) %>%
    head(n=9)
  return(filtered_terms)
}

# Apply prep_terms function to gprofiler outputs
contractile_smc_df = prep_terms(smc_go_bp, "SMC")
pioneer_df = prep_terms(pioneer_go_bp, "Transitional-SMC")
fibromyo_df = prep_terms(fibromyo_go_bp, "Fibromyocyte")
fibrochondro_df = prep_terms(fibrochondro_go_bp_non_redundant, "Fibrochondrocyte")
foam_like_smc_df = prep_terms(foam_like_go_bp, "Foam_like_SMC")

# Merge dfs into a single object for plotting the heatmap
merged_df = rbind(contractile_smc_df, pioneer_df, fibromyo_df, fibrochondro_df, foam_like_smc_df)
new_df = merged_df %>%
  dplyr::select(annotation, term_name, log10_pval)
write.table(new_df,
            "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_tables/GSEA_SMC_annotations_top_20_GO_BP_terms.tsv",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

reshaped_df = spread(new_df, key = annotation, value = log10_pval)
reshaped_df[is.na(reshaped_df)] = 0
go_matrix = reshaped_df %>%
  dplyr::select(-term_name) %>% 
  as.matrix()
rownames(go_matrix) = reshaped_df$term_name

# Plot heatmap
palette_length = 100
my_color3 = colorRampPalette(c("white", brewer.pal(n=9, name = "Blues")))(palette_length)

pheatmap(go_matrix, color = my_color3, fontsize = 12, angle_col = 45)

