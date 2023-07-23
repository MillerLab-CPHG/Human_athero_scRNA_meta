library(Seurat)
library(harmony)
library(rliger)
library(tidyverse)
library(data.table)
library(cluster)
library(RColorBrewer)
library(parallel)
library(reticulate)

# Set seed for reproducibility
set.seed(1)


# Source our own scRNA analysis utils functions
#source("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/R_scripts/Utils/scRNA_processing_utils.R")


########################################################################
# BENCHMARK INTEGRATION ALGORITHMS FOR ALSAIGH, PAN AND WIRKA ET AL DATA

#######################################
# Load SCTransform processed datasets #
#######################################

# Load Pan et al data
pan_rpe004 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/pan_et_al_rds_objects/pan_rpe004_processed_seurat_obj.rds")
pan_rpe005 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/pan_et_al_rds_objects/pan_rpe005_processed_seurat_obj.rds")
pan_rpe006 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/pan_et_al_rds_objects/pan_rpe006_processed_seurat_obj.rds")

# Load Alsaigh et al data
ac_p1 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_ac_p1_processed_seurat_obj.rds")
pa_p1 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_pa_p1_processed_seurat_obj.rds")

ac_p2 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_ac_p2_processed_seurat_obj.rds")
pa_p2 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_pa_p2_processed_seurat_obj.rds")

ac_p3 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_ac_p3_processed_seurat_obj.rds")
pa_p3 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/alsaigh_et_al_rds_objects/alsaigh_pa_p3_processed_seurat_obj.rds")


# Load Wirka et al
wirka_hca = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/wirka_et_al_rds_objects/wirka_human_coronaries_processed_seurat_obj.rds")



######################
# Standard CCA + MNN #
######################

# Prep datasets for integration
cca_obj_list = list(rpe004=pan_rpe004,
                    rpe005=pan_rpe005,
                    rpe006=pan_rpe006,
                    ac_p1 = ac_p1,
                    pa_p1 = pa_p1,
                    ac_p2 = ac_p2,
                    pa_p2 = pa_p2,
                    ac_p3 = ac_p3,
                    pa_p3 = pa_p3,
                    wirka_hca=wirka_hca)

# Select variable features for integration
cca_features_sct = SelectIntegrationFeatures(cca_obj_list, nfeatures = 3000)
cca_obj_list = PrepSCTIntegration(cca_obj_list, anchor.features = cca_features_sct)


# Start timer for benchmark
start_time_cca = Sys.time()

# Find integration anchors and integrate data
cca_sct_anchors = FindIntegrationAnchors(object.list = cca_obj_list, 
                                         normalization.method = "SCT",
                                         anchor.features = cca_features_sct)
cca_int_sct = IntegrateData(anchorset = cca_sct_anchors, normalization.method = "SCT")

# Measure time for benchmark
end_time_cca = Sys.time()

# Time difference of 1.161337 hours (69.7 mins)
measured_time_cca = end_time_cca - start_time_cca 
cca_running_time = measured_time_cca

# Save integrated seurat object without running dim reduction and clustering (so we can play with diff parameters 
# wihtoud datasets to integrate
saveRDS(cca_int_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/CCA/cca_integrated_no_dimReduction.rds")
cca_int_sct = readRDS("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/CCA/cca_integrated_no_dimReduction.rds")

# Run PCA and UMAP on integrated object (default is now integrated object)
DefaultAssay(cca_int_sct) = "integrated"
cca_int_sct = RunPCA(cca_int_sct)
cca_int_sct = FindNeighbors(cca_int_sct, reduction = "pca", dims = 1:30, k.param = 20)
cca_int_sct = RunUMAP(cca_int_sct, dims = 1:30, n.neighbors = 30)

# Resolution 1.3 seems to yield an appropriate number of clusters and distinguish fibrochondrocytes from other ECM rich cells
# There doesn't seem to be that much difference between a res parameter of 1.7-1.9
cca_int_sct = FindClusters(cca_int_sct, resolution = 1.3)

# Visualization of clusters
p1_cca = DimPlot(cca_int_sct, reduction = "umap", group.by = "sample", pt.size = 0.1) + ggtitle("CCA+MNN Wirka/Pan samples") + 
  custom_theme() + miller_discrete_scale(option = 2)
p2_cca = DimPlot(cca_int_sct, reduction = "umap", group.by = "study", pt.size = 0.1) +  ggtitle("CCA+MNN Wirka/Pan study") + 
  custom_theme() + miller_discrete_scale(option = 2) + theme(legend.position = "none")
p3_cca = DimPlot(cca_int_sct, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.1) + ggtitle("CCA+MNN Wirka/Pan clusters") + 
  custom_theme() 
p1_sct + p2_sct

# Save plot by study
ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure2/revision_figures/SuppFig2a_CCA_UMAP_embeddings.png",
       plot = p2_cca, width = 10, height = 10)

# Visualize gene expression
DefaultAssay(cca_int_sct) = "SCT"
FeaturePlot(cca_int_sct ,features = c("MYH11", "TPM2", "TNFRSF11B", "TCF21"), order = TRUE,  pt.size = 0.1, slot = "data") & new_scale

FeaturePlot(cca_int_sct ,features = c("CNN1", "KLF4", "CRTAC1", "COMP"), 
            order = TRUE,  pt.size = 0.2, slot = "data") & new_scale

FeaturePlot(cca_int_sct ,features = c("IL1B", "TREM2", "S100A9", "APOE" , "TNF", "LYZ"), 
            order = TRUE,  pt.size = 0.2, slot = "data") & new_scale

FeaturePlot(cca_int_sct ,features = c("IL1B", "NKG7", "S100A9", "APOE" , "CD163", "C7"), 
            order = TRUE,  pt.size = 0.2, slot = "data") & new_scale



################################################################
# rPCA                                                         #
################################################################

rpca_obj_list = list(rpe004=pan_rpe004,
                     rpe005=pan_rpe005,
                     rpe006=pan_rpe006,
                     ac_p1 = ac_p1,
                     pa_p1 = pa_p1,
                     ac_p2 = ac_p2,
                     pa_p2 = pa_p2,
                     ac_p3 = ac_p3,
                     pa_p3 = pa_p3,
                     wirka_hca=wirka_hca)

rpca_features_sct = SelectIntegrationFeatures(rpca_obj_list, nfeatures = 3000)
rpca_obj_list = PrepSCTIntegration(rpca_obj_list, anchor.features = rpca_features_sct)
rpca_obj_list = lapply(rpca_obj_list, FUN = function(x) {
  x = RunPCA(x, features=rpca_features_sct)
}) 
                       

# Start timer for benchmark
start_time_rpca = Sys.time()

# Find integration anchors and integrate data
# rPCA integration is tipically more conservative (less likely to align cells in different states
# like naive and memory T cells) across experiments. The k.anchor parameter allows us to modify 
# the strength of the alignment and is by default set to 5. Let's start with 5 and see if that 
# separates SMCs better from pericytes and fibroblasts. 

# How do integrating through set dims=1:50 compare to our current integration (dims=1:30)? 
# Would this be able to provide any additional information to separate cell types better?
# Setting dims=1:50 doesn't really seem to provide a substantial benefit in separating cell types in low dim space
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

# Time difference of 23.13806 mins
measured_time_rpca = end_time_rpca - start_time_rpca 
rpca_running_time = measured_time_rpca

# Save integrated seurat object without running dim reduction and clustering (so we can play with diff parameters 
# without waiting for datasets to integrate). 
saveRDS(rpca_int_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/rpca_integrated_no_dimReduction.rds")
saveRDS(rpca_int_sct, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/rpca_integrated50dims_no_dimReduction.rds")

rpca_int_sct = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_int/rpca_integrated_no_dimReduction.rds")


# Run PCA and UMAP on integrated object (default is now integrated object)
DefaultAssay(rpca_int_sct) = "integrated"
rpca_int_sct = RunPCA(rpca_int_sct) # 30 dims captures really most of the variation in the data
rpca_int_sct = FindNeighbors(rpca_int_sct, reduction = "pca", dims = 1:30, k.param = 20)
# In general, this parameter should be 5-50. Higher values result in more global structure being preserved at the loss of detailed local structure
# UMAP looks better with 30 n.neighbors than 20
rpca_int_sct = RunUMAP(rpca_int_sct, dims = 1:30, n.neighbors = 30, min.dist = 0.3) # Adjusting min.dist doesn't really help

# A resolution of 1 seems appropriate for evaluating different methods for now
# Update: res of 1.2-1.4 provide the highest mean silhouette scores and ~ 30 clusters
rpca_int_sct = FindClusters(rpca_int_sct, resolution = 1.3)

# Plot cell embeddings
p1_rpca = DimPlot(rpca_int_sct, reduction = "umap", group.by = "sample", pt.size = 0.1) + ggtitle("rPCA Alsaigh/Pan/Wirka samples") + 
  custom_theme() + miller_discrete_scale(option = 2)
p2_rpca = DimPlot(rpca_int_sct, reduction = "umap", group.by = "study", pt.size = 0.1) + ggtitle("rPCA Alsaigh/Pan/Wirka study") + 
  custom_theme() + miller_discrete_scale(option = 2) + theme(legend.position = "none")
p3_rpca = DimPlot(rpca_int_sct, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.1) + ggtitle("rPCA Alsaigh/Pan/Wirka clusters") + 
  custom_theme() 
p1_rpca + p2_rpca

# Save plot by study
ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure2/revision_figures/SuppFig2a_rPCA_UMAP_embeddings.png",
       plot = p2_rpca, width = 10, height = 10)

# Visualize gene expression
DefaultAssay(rpca_int_sct) = "SCT"
FeaturePlot(rpca_int_sct ,features = c("CNN1", "ACTA2", "VCAN", "CRTAC1", "FBLN1", 
                                      "ACTA2"), order = TRUE,  pt.size = 0.2, slot = "data") & miller_continous_scale()


# NOTES: Using rPCA yields a better separation between SMC-derived cells and fibroblasts than CCA and Harmony (mostly from Wirka et al)
# Separation between SMCs and pericytes seems pretty comparable to Harmony. Running time is fairly efficient but substantially slower 
# than Harmony. Also some gene expression patters (CRTAC1 and COMP) seem better delineated than in CCA and Harmony.
# This approach seems like a pretty viable option given the number of cells we have. 
# Better separation of expression patterns between SMCs, fibromyocytes, fibrochondrocytes and fibroblasts (don't get a lot of ACTA2,
# or CNN1 in this cluster). 
# In general, rPCA seems to conserve cell type specific experssion patterns better introducing less artifacts and with a lower
# degree of alignment of cells in different biological states. 

################################################################
# INTEGRATE WIRKA ET AL TO INTEGRATED PAN ET AL USING HARMONY  #
################################################################

h_test_list = list(rpe004=pan_rpe004,
                   rpe005=pan_rpe005,
                   rpe006=pan_rpe006,
                   ac_p1 = ac_p1,
                   pa_p1 = pa_p1,
                   ac_p2 = ac_p2,
                   pa_p2 = pa_p2,
                   ac_p3 = ac_p3,
                   pa_p3 = pa_p3,
                   wirka_hca=wirka_hca)

h_test_features_sct = SelectIntegrationFeatures(h_test_list, nfeatures = 3000)

# The following merging step combines seurat objects but cleans var features and reduced dim objects
# Harmony embeds cells in PCA space and aligns them according to cell type to remove batch effects. 
# We use SelectIntegrationFeatures() in order to find var features across the datsasets in the list to be merged
# These var features will be added once the seurat objects are merged, since merging does not conserve var features 
harmony_seurat_merge = merge(pan_rpe004, y = c(pan_rpe005, pan_rpe006, ac_p1, pa_p1, ac_p2, pa_p2, ac_p3, pa_p3, wirka_hca),
                             add.cell.ids = c("pan_rpe004", "pan_rpe005", "pan_rpe006", 
                                              "alsaigh_ac_p1", "alsaigh_pa_p1", "alsaigh_ac_p2", "alsaigh_pa_p2", "alsaigh_ac_p3", "alsaigh_pa_p3",
                                              "wirka_hca"),
                             merge.data = TRUE)
VariableFeatures(harmony_seurat_merge) = h_test_features_sct


# Start timer for benchmark
start_time_harmony = Sys.time()

harmony_seurat_int = RunPCA(harmony_seurat_merge) %>%
  RunHarmony(assay.use = "SCT", group.by.vars = "sample", dims.use=1:30) %>%
  RunUMAP(dims = 1:30, n.neighbors = 30, reduction="harmony") %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, k.param = 20)

# Cluster data 
harmony_seurat_int =  FindClusters(harmony_seurat_int ,resolution = 1.2) # Resolution of 1 yields a similar number of clusters as in the above methods

# Measure time for benchmark
end_time_harmony = Sys.time()

# Time difference of 5.433346 mins
measured_time_harmony = end_time_harmony - start_time_harmony 
harmony_running_time = measured_time_harmony

# Save rds object
saveRDS(harmony_seurat_int, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Harmony/harmony_integrated_30dims_clustered_seurat_obj.rds")
harmony_seurat_int = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Harmony/harmony_integrated_30dims_clustered_seurat_obj.rds")

# Plot embeddings produced by Harmony
p1_harmony = DimPlot(harmony_seurat_int, reduction = "umap", group.by = "sample", pt.size = 0.1) + ggtitle("Harmony Alsaigh/Pan/Wirka samples") + 
  custom_theme() + miller_discrete_scale(option = 2)
p2_harmony = DimPlot(harmony_seurat_int, reduction = "umap", group.by = "study", pt.size = 0.1) + ggtitle("Harmony Alsaigh/Pan/Wirka study") + 
  custom_theme() + miller_discrete_scale(option = 2) + theme(legend.position = "none")
p3_harmony = DimPlot(harmony_seurat_int, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.1) + ggtitle("Harmony Alsaigh/Pan/Wirka clusters") + 
  custom_theme() 

# Save plot by study
ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure2/revision_figures/SuppFig2a_Harmony_UMAP_embeddings.png",
       plot = p2_harmony, width = 10, height = 10)


# Plot expression of some representative markers
FeaturePlot(harmony_seurat_int, features = c("CNN1", "ACTA2", "VCAN", "CRTAC1", "FBLN1", "ACTA2"), 
            order = TRUE,  pt.size = 0.2, slot = "data") & miller_continous_scale()


# NOTES:Harmony doesn't provide a very well defined separation between putative fibromyocytes/fibrochondrocytes and fibroblasts.
# CRTAC1 and COMP expression seem more spread out and overlapping fibroblasts markers at some extent. Separation between
# SMCs and pericytes seems pretty comparable to that of rPCA. 
# There are also some cells in the T/NK clusters expressing SMC markers strongly but they're not very well clustered. Number of 
# clusters is fairly comparable to rPCA (43-44 clusters)


######################
# Scanorama          #
######################

# Import scanorama with reticulate (make sure the proper python version is loaded)
scanorama = import("scanorama")

# Scanorama requires a list of matrices and a list of list of genes
seurat_objects = list(ac_p1 = ac_p1,
                      pa_p1 = pa_p1,
                      ac_p2 = ac_p2,
                      pa_p2 = pa_p2,
                      ac_p3 = ac_p3,
                      pa_p3 = pa_p3,
                      rpe004 = pan_rpe004,
                      rpe005 = pan_rpe005,
                      rpe006 = pan_rpe006,
                      wirka_hca = wirka_hca)

matrix_list = list()
genes_list = list()

# Get normalized data from SCT assays
for (i in 1:length(seurat_objects)) {
  matrix_list[[i]] = t(as.matrix(seurat_objects[[i]]@assays$SCT@data))
  genes_list[[i]] = rownames(seurat_objects[[i]])
}

matrix_list[[1]]

# Start timer for benchmark
start_time_scanorama = Sys.time()

# Run scanorama 
integrated_corrected_data = scanorama$correct(datasets_full = matrix_list,
                                              genes_list = genes_list,
                                              return_dimred = TRUE,
                                              return_dense = TRUE)

# Measure time for benchmark
end_time_scanorama = Sys.time()

# Time difference of 82.242552 mins
measured_time_scanorama = end_time_scanorama - start_time_scanorama 
scanorama_running_time = measured_time_scanorama

# Save integrated/batch corrected data
saveRDS(integrated_corrected_data, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Scanorama/Scanorama_integrated_corrected_data.rds")
integrated_corrected_data = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Scanorama/Scanorama_integrated_corrected_data.rds")


# Add back integrated data into seurat obj
int_data = lapply(integrated_corrected_data[[2]], t)
panorama = do.call(cbind, int_data)

# Use common genes to assign rownames to the integrated, batch-corrected expression matrix
rownames(panorama) = as.character(integrated_corrected_data[[3]])
colnames(panorama) = unlist(sapply(matrix_list, rownames))
head(panorama)

# Dim 12681 genes x 61823 cells
dim(panorama)

# Process dim reduction embeddings
int_dim_red = do.call(rbind, integrated_corrected_data[[1]])
colnames(int_dim_red) = paste0("PC_", 1:100)

# Add SD to draw wlbow plots in Seurat
stdevs =  apply(int_dim_red, MARGIN = 2, FUN = sd)

# Create Seurat obj. Skip normalization and variable gene selection since we already have the scanorama PCA embeddings
scanorama_seurat = CreateSeuratObject(counts = panorama, 
                                      assay = "scanorama_integrated",
                                      project = "pan_wirka_alsaigh_integration")

# Add metadata from previous objects
metadata = lapply(seurat_objects, function(x){x@meta.data})
metadata_df = data.table::rbindlist(metadata, use.names = FALSE)
rownames(metadata_df) = colnames(scanorama_seurat)
rownames(int_dim_red) = colnames(scanorama_seurat)

# Create dim reduced object in Seurat
scanorama_seurat@meta.data = metadata_df
scanorama_seurat[["pca"]] = CreateDimReducObject(embeddings = int_dim_red, stdev = stdevs,
                                                 key = "PC_", assay = "scanorama_integrated")

# Save rds object: scanorama integrated seurat object with PCA dim reduced embeddings, but no clusters or UMAP embeddings
saveRDS(scanorama_seurat, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Scanorama/scanorama_int_seurat_obj_dims_reduced_no_clusters.rds")
scanorama_seurat = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Scanorama/scanorama_int_seurat_obj_dims_reduced_no_clusters.rds")

# Proceed with Seurat standard workflow (findneighbors, and UMAP embeddings)
scanorama_seurat = FindNeighbors(scanorama_seurat, reduction = "pca", dims = 1:30, k.param = 20)
# In general, this parameter should be 5-50. Higher values result in more global structure being preserved at the loss of detailed local structure
# UMAP looks better with 30 n.neighbors than 20
scanorama_seurat = RunUMAP(scanorama_seurat, dims = 1:30, n.neighbors = 30, min.dist = 0.3) # Adjusting min.dist doesn't really help

# A resolution of 1 seems appropriate for evaluating different methods for now
# Update: res of 1.1 provide the highest mean silhouette scores and ~ 30 clusters
scanorama_seurat = FindClusters(scanorama_seurat, resolution = 1.1)

# Add study metadata for benchmarking
scanorama_seurat$study = rpca_int_sct@meta.data$study

# Visualization of clusters
p1_scanorama = DimPlot(scanorama_seurat, reduction = "umap", group.by = "sample", pt.size = 0.1) + ggtitle("Scanorama Alsaigh/Pan/Wirka samples") + 
  custom_theme() + miller_discrete_scale(option = 2)
p2_scanorama = DimPlot(scanorama_seurat, reduction = "umap", group.by = "study", pt.size = 0.1) + ggtitle("Scanorama Alsaigh/Pan/Wirka study") + 
  custom_theme() + miller_discrete_scale(option = 2) + theme(legend.position = "none")
p3_scanorama = DimPlot(scanorama_seurat, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.1) + ggtitle("Scanorama Alsaigh/Pan/Wirka clusters") + 
  custom_theme() 
p1_scanorama + p2_scanorama

# Save plot by study
ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure2/revision_figures/SuppFig2a_Scanorama_UMAP_embeddings.png",
       plot = p2_scanorama, width = 10, height = 10)

# Visualize gene expression
FeaturePlot(scanorama_seurat ,features = c("CNN1", "ACTA2", "VCAN", "CRTAC1", "FBLN1", 
                                       "NET1"), order = TRUE,  pt.size = 0.2, slot = "data") & miller_continous_scale()

# NOTES: separation between cell types seems decent but clustering quality overall poor. There's also a lack
# of a substantial number of meaningufl genes like KRT17 or CD79A. 

###########################################
# Plotting running time for each approach #
###########################################

int_running_times_df = data.frame(integration_method = c("CCA_MNN", "rPCA", "Harmony", "Scanorama"),
                               running_time_minutes = c(cca_running_time, rpca_running_time,
                                                        harmony_running_time, scanorama_running_time))


int_running_times_df$integration_method = factor(int_running_times_df$integration_method, 
                                                 levels = c("Scanorama", "CCA_MNN", "rPCA", "Harmony"))


# Make running time plot
running_time_plot = int_running_times_df %>%
  #fct_reorder(integration_method, running_time_minutes) %>%
  ggplot(aes(x=reorder(integration_method, running_time_minutes), 
             y=running_time_minutes, fill=integration_method)) + 
  geom_col(width = 0.6) + 
  xlab("Method") +
  ylab("Running time (mins)") + 
  ggtitle("Running time benchmark for integration of Alsaigh/Pan/Wirka datasets") +  
  custom_theme() +
  scale_fill_manual(values = c("#4DBBD5FF", "#E64B35FF", "#3C5488FF", "#00A087FF")) +  
  theme(aspect.ratio = 1.3,
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none")

ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure2/revision_figures/SuppFig2b_integration_running_time_benchmark.pdf",
       plot = running_time_plot, width = 3, height = 4)









