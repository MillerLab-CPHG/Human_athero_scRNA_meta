library(kBET)
library(lisi)
library(Seurat)
library(data.table)
library(tidyverse)
library(cluster)


# Source scRNA custom utils script
source("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/R_scripts/Utils/scRNA_processing_utils.R")

###############################################
# Calculate Sihouette scores for each method  #
###############################################


# Silhouette scores measure clustering purity by taking into account cohesion and separation of cells
# in clusters. In other words, how similar a cells is to its own cluster compared to other clusters.
# This gives a score in the range of −1 to +1, where a higher score indicates higher performance.
# The silhouette coefficient captures elements of both sample mixing and local structure.
# Here we'll calculate sil scores across a range of clustering resolutions to control for the 
# Seurat clustering granularity parameter. For this calculations, we'll source the wrapper from the utils script

###############################################
# CCA

# Test a few resolution parameters and see which ones produce the higher silhouette scores
# Figure out how many clusters each res yields and what the overall sil scores are
cca_res = c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8)
cca_scores_list = list()
cca_dist_matrix = dist(x = Embeddings(object = cca_int_sct[["pca"]])[, 1:30])
saveRDS(cca_dist_matrix, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/CCA/cca_distance_matrix.rds")
cca_dist_matrix = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/CCA/cca_distance_matrix.rds")

for (i in seq_along(cca_res)) {
  # Cluster data
  cca_int_sct = FindClusters(cca_int_sct, 
                             resolution = cca_res[i])
  
  # Calculate silhouette scores
  cca_scores_list[[i]] = calc_sil_scores(cca_int_sct, cca_dist_matrix)
  names(cca_scores_list)[[i]] = paste("res", as.character(cca_res[i]), sep = "_")
}

# Create a single df with all the sil scores for the tested resolutions
cca_resolutions_sil_df = rbindlist(cca_scores_list, idcol = TRUE)
cca_resolutions_sil_df$method = "CCA"
dim(cca_resolutions_sil_df)
saveRDS(cca_resolutions_sil_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/CCA/CCA_sil_scores_per_clustering_resolution.rds")
cca_resolutions_sil_df =  read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/CCA/CCA_sil_scores_per_clustering_resolution.rds")

# Plot silhouette scores
resolutions_sil_df %>%
  ggplot(aes(x=sil_width, color=.id)) + geom_density() +
  ggtitle("CCA sil scores")

# Plot mean silhouette score per resolution.
# Max sil score 0.15
cca_resolutions_sil_df %>%
  group_by(.id) %>%
  summarize(mean_sil = mean(sil_width)) %>%
  ggplot(aes(x=.id, y=mean_sil)) + geom_col() +
  ggtitle("CCA sil scores") + 
  custom_theme +
  xlab("Clustering resolution") + 
  ylab("Mean Silhouette score") + 
  theme(axis.text.x = element_text(angle=70))



################################
# rPCA

# Figure out how many clusters each res yields and what the overall sil scores are
rpca_res = c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8)
rpca_scores_list = list()
dist_matrix = dist(x = Embeddings(object = rpca_int_sct[["pca"]])[, 1:30])
saveRDS(dist_matrix, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/distance_matrix.rds")
dist_matrix = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/distance_matrix.rds")

for (i in seq_along(res)) {
  # Cluster data
  rpca_int_sct = FindClusters(rpca_int_sct, 
                              resolution = rpca_res[i])
  
  # Calculate silhouette scores
  rpca_scores_list[[i]] = calc_sil_scores(rpca_int_sct, dist_matrix)
  names(rpca_scores_list)[[i]] = paste("res", as.character(res[i]), sep = "_")
}

# Create a single df with all the sil scores for the tested resolutions
rpca_resolutions_sil_df = rbindlist(rpca_scores_list, idcol = TRUE)
rpca_resolutions_sil_df$method = "rPCA"
saveRDS(resolutions_sil_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/sil_scores_per_clustering_resolution.rds")
rpca_resolutions_sil_df = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_int/sil_scores_per_clustering_resolution.rds")

# Plot silhouette scores
rpca_resolutions_sil_df %>%
  ggplot(aes(x=sil_width, color=.id)) + geom_density() +
  ggtitle("rPCA sil scores")

# Plot mean silhouette score per clustering resolution
rpca_resolutions_sil_df %>%
  group_by(.id) %>%
  summarize(mean_sil = mean(sil_width)) %>%
  ggplot(aes(x=.id, y=mean_sil)) + geom_col() +
  ggtitle("rPCA sil scores") + 
  custom_theme +
  xlab("Clustering resolution") + 
  ylab("Mean Silhouette score") + 
  theme(axis.text.x = element_text(angle=70))

############################
# Harmony

# Test a few resolution parameters and see which ones produce the higher silhouette scores
# Figure out how many clusters each res yields and what the overall sil scores are
h_res = c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8)
h_scores_list = list()
h_dist_matrix = dist(x = Embeddings(object = harmony_seurat_int[["harmony"]])[, 1:30])
saveRDS(h_dist_matrix, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Harmony/distance_matrix.rds")
h_dist_matrix = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Harmony/distance_matrix.rds")

for (i in seq_along(res)) {
  # Cluster data
  harmony_seurat_int = FindClusters(harmony_seurat_int, 
                                    resolution = h_res[i])
  
  # Calculate silhouette scores
  h_scores_list[[i]] = calc_sil_scores(harmony_seurat_int, h_dist_matrix)
  names(h_scores_list)[[i]] = paste("res", as.character(res[i]), sep = "_")
}


# Create a single df with all the sil scores for the tested resolutions
# Max res
h_resolutions_sil_df = rbindlist(h_scores_list, idcol = TRUE)
h_resolutions_sil_df$method = "harmony"
#res_col_names = paste("res", res, sep="_")
#colnames(h_resolutions_sil_df) = res_col_names
saveRDS(h_resolutions_sil_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Harmony/harmony_sil_scores_per_clustering_resolution.rds")
h_resolutions_sil_df = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Harmony/harmony_sil_scores_per_clustering_resolution.rds")

# Plot silhouette scores
h_resolutions_sil_df %>%
  ggplot(aes(x=sil_width, color=.id)) + geom_density() +
  ggtitle("Harmony sil scores")

# Plot mean silhouette score per resolution. Res 1.1 and 1.2 yield the highest avg scores
# Max sil score is 0.115
h_resolutions_sil_df %>%
  group_by(.id) %>%
  summarize(mean_sil = mean(sil_width)) %>%
  ggplot(aes(x=.id, y=mean_sil)) + geom_col() +
  ggtitle("Harmony sil scores") + 
  custom_theme +
  scale_x_discrete(labels=res_col_names) +
  xlab("Clustering resolution") + 
  ylab("Mean Silhouette score") + 
  theme(axis.text.x = element_text(angle=70))



################################################
# Scanorama

# Test a few resolution parameters and see which ones produce the higher silhouette scores
# Figure out how many clusters each res yields and what the overall sil scores are
scanorama_res = c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8)
scanorama_scores_list = list()
scanorama_dist_matrix = dist(x = Embeddings(object = scanorama_seurat[["pca"]])[, 1:30])
saveRDS(scanorama_dist_matrix, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Scanorama/scanorama_distance_matrix.rds")
scanorama_dist_matrix = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Scanorama/scanorama_distance_matrix.rds")

for (i in seq_along(scanorama_res)) {
  # Cluster data
  scanorama_seurat = FindClusters(scanorama_seurat, 
                                  resolution = scanorama_res[i])
  
  # Calculate silhouette scores
  scanorama_scores_list[[i]] = calc_sil_scores(scanorama_seurat, scanorama_dist_matrix)
  names(scanorama_scores_list)[[i]] = paste("res", as.character(scanorama_res[i]), sep = "_")
}

# Create a single df with all the sil scores for the tested resolutions
scanorama_resolutions_sil_df = rbindlist(scanorama_scores_list, idcol = TRUE)
scanorama_resolutions_sil_df$method = "Scanorama"
dim(scanorama_resolutions_sil_df)
saveRDS(scanorama_resolutions_sil_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Scanorama/Scanorama_sil_scores_per_clustering_resolution.rds")
rpca_res_sil_df = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_int/sil_scores_per_clustering_resolution.rds")

# Plot silhouette scores
scanorama_resolutions_sil_df %>%
  ggplot(aes(x=sil_width, color=.id)) + geom_density() +
  ggtitle("rPCA sil scores")

# Plot mean silhouette score per resolution. 
scanorama_resolutions_sil_df %>%
  group_by(.id) %>%
  summarize(mean_sil = mean(sil_width)) %>%
  ggplot(aes(x=.id, y=mean_sil)) + geom_col() +
  ggtitle("Scanorama sil scores") + 
  custom_theme +
  xlab("Clustering resolution") + 
  ylab("Mean Silhouette score") + 
  theme(axis.text.x = element_text(angle=70))

########################################################
# Compare rPCA vs CCA vs Harmony vs Scanorama sil scores 
res_sil_scores_df = rbind.data.frame(rpca_res_sil_df, h_resolutions_sil_df, cca_resolutions_sil_df, scanorama_resolutions_sil_df)
saveRDS(res_sil_scores_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA_Harmony_CCA_Scanorama_silhouette_scores_per_clustering_res.rds")
res_sil_scores_df = readRDS("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA_Harmony_CCA_Scanorama_silhouette_scores_per_clustering_res.rds")

dim(res_sil_scores_df)
head(res_sil_scores_df)
table(res_sil_scores_df$method)

res_sil_scores_df$method = factor(res_sil_scores_df$method, 
                                  levels = c("rPCA", "CCA", "Scanorama", "harmony"))

res_sil_scores_plot = res_sil_scores_df %>%
  group_by(.id, method) %>%
  summarize(mean_sil = mean(sil_width)) %>%
  ggplot(aes(x=.id, y=mean_sil, fill=method)) + geom_col(position="dodge", width = 0.7) +
  xlab("Clustering resolution") +
  ylab("Mean Silhouette score") +
  ggtitle("rPCA vs Harmony vs CCA vs Scanorama mean sil scores per clustering resolution") + 
  custom_theme + 
  npg_scale_bars + 
  scale_fill_manual(values = c("#3C5488FF", "#E64B35FF", "#4DBBD5FF", "#00A087FF")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure1/SuppFig1d_rPCA_Harmony_CCA_Scanorama_sil_coeffs.pdf",
       plot = res_sil_scores_plot, width = 4.5, height = 4.5)



################################
# kBET                         #
################################

# # Run a few tests with kBET
# cca_int_sct = readRDS("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/CCA/cca_integrated_no_dimReduction.rds")
# 
# # Get the matrix
# cca_sct_counts = cca_int_sct@assays$SCT@data
# cca_sct_counts_t = t(cca_sct_counts)
# 
# 
# # Get PCA embeddings for downstream calculations
# cca_pca_embeddings = Embeddings(object=cca_int_sct[["pca"]])[, 1:30]
# 
# # Get batch labels 
# batch_labels = as.factor(cca_int_sct@meta.data$study)
# cell_barcodes = rownames(cca_int_sct@meta.data)
# names(cell_barcodes) = batch_labels
# 
# # Calculate kBET
# batch_estimate = kBET(Embeddings(object = cca_int_sct[["pca"]])[, 1:30], batch=batch_labels, verbose = TRUE)

#########################
# Calculate LISI scores #
#########################

#########################
# CCA
cca_pca_embeddings = Embeddings(object=cca_int_sct[["pca"]])[, 1:30]
cca_batch_metadata = data.frame(batch_label=cca_int_sct@meta.data$study, 
                                cluster=cca_int_sct@meta.data$integrated_snn_res.1.3)
rownames(cca_batch_metadata) = rownames(cca_int_sct@meta.data)
cca_batch_metadata$method = rep("CCA_MNN", 
                                length(rownames(cca_batch_metadata)))


# Calculate iLISI
# Mean iLISI = 1.796; Median iLISI = 1.788
cca_ilisi_res = lisi::compute_lisi(cca_pca_embeddings, cca_batch_metadata, c("batch_label"))
cca_batch_metadata$ilisi = cca_ilisi_res$batch_label


# Plot lisi scores distribution by cluster
cca_ilisi_means = cca_batch_metadata %>%
  group_by(cluster) %>%
  summarize(mean_ilisi = mean(ilisi)) %>%
  ggplot(aes(x=cluster, y=mean_ilisi)) + 
  geom_col() + 
  facet_wrap(~cluster, scales = "free_y") + 
  custom_theme

# Try a boxplot instead
# cca_batch_metadata %>%
#   ggplot(aes(x=method, y=ilisi)) + 
#   geom_boxplot(aes(x=method, y=ilisi)) + 
#   geom_jitter()

# Calculate cLISI
# mean cLISI = 1.22; median clISI = 1.02
cca_clisi_res = lisi::compute_lisi(cca_pca_embeddings, 
                                   cca_batch_metadata, c("cluster"))
cca_batch_metadata$clisi = cca_clisi_res$cluster

# plot cLISI scores
cca_batch_metadata %>%
  ggplot(aes(x=clisi)) + 
  geom_density() + 
  xlab("cLISI") + 
  custom_theme



##########################
# rPCA
rpca_pca_embeddings = Embeddings(object = rpca_int_sct[["pca"]])[, 1:30]
rpca_batch_metadata = data.frame(batch_label = rpca_int_sct@meta.data$study,
                                 cluster=rpca_int_sct@meta.data$integrated_snn_res.1.3)
rownames(rpca_batch_metadata) = rownames(cca_int_sct@meta.data)
rpca_batch_metadata$method = rep("rPCA", 
                                 length(rownames(rpca_batch_metadata)))

# Calculate iLISI
# Mean ilisi = 1.61; Median ilisi = 1.484
rpca_ilisi_res = lisi::compute_lisi(rpca_pca_embeddings, rpca_batch_metadata, 
                                    c("batch_label"))
rpca_batch_metadata$ilisi = rpca_ilisi_res$batch_label


# Plot results
rpca_batch_metadata %>%
  #group_by(cluster) %>%
  #summarize(mean_ilisi = mean(ilisi)) %>%
  ggplot(aes(x=ilisi)) + 
  geom_density() + 
  facet_wrap(~cluster, scales = "free_y") + 
  custom_theme

# rpca_batch_metadata %>%
#   ggplot(aes(x=method, y=ilisi)) + 
#   geom_boxplot(aes(x=method, y=ilisi)) 

# Calculate cLISI
# mean cLISI = 1.159; median cLISI = 1.002
rpca_clisi_res = lisi::compute_lisi(rpca_pca_embeddings, rpca_batch_metadata, 
                                    c("cluster"))
rpca_batch_metadata$clisi = rpca_clisi_res$cluster
mean(rpca_clisi_res$cluster)
median(rpca_clisi_res$cluster)


# plot cLISI scores
rpca_batch_metadata %>%
  ggplot(aes(x=clisi)) + 
  geom_density() + 
  custom_theme


############################
# Harmony
harmony_pca_embeddings = Embeddings(object = harmony_seurat_int[["pca"]])[, 1:30]
harmony_batch_metadata = data.frame(batch_label = harmony_seurat_int@meta.data$study,
                                    cluster=harmony_seurat_int@meta.data$seurat_clusters)
rownames(harmony_batch_metadata) = rownames(harmony_seurat_int@meta.data)
harmony_batch_metadata$method = rep("Harmony", 
                                    length(rownames(harmony_batch_metadata)))

# Calculate iLISI scores
# Mean ilisi = 1.08; Median ilisi = 1
harmony_ilisi_res = lisi::compute_lisi(harmony_pca_embeddings, harmony_batch_metadata, c("batch_label"))
harmony_batch_metadata$ilisi = harmony_ilisi_res$batch_label

# plot iLISI scores
harmony_batch_metadata %>%
    group_by(cluster) %>%
    summarize(mean_ilisi = mean(ilisi)) %>%
    ggplot(aes(x=cluster, y=mean_ilisi)) + 
    geom_col() + 
    #facet_wrap(~cluster, scales = "free_y") + 
    custom_theme

# harmony_batch_metadata %>%
#   ggplot(aes(x=method, y=ilisi)) + 
#   geom_boxplot(aes(x=method, y=ilisi))

# Calculate cLISI scores
# mean cLISI = 1.385; median clISI = 1.134
harmony_clisi_res = lisi::compute_lisi(harmony_pca_embeddings, harmony_batch_metadata, 
                                       c("cluster"))
harmony_batch_metadata$clisi = harmony_clisi_res$cluster
mean(harmony_clisi_res$cluster)
median(harmony_clisi_res$cluster)

# plot cLISI scores
harmony_batch_metadata %>%
  ggplot(aes(x=clisi)) + 
  geom_density() + 
  custom_theme


##############################
# Scanorama

scanorama_pca_embeddings = Embeddings(object = scanorama_seurat[["pca"]])[, 1:30]
scanorama_batch_metadata = data.frame(batch_label = scanorama_seurat_int@meta.data$study,
                                      cluster = scanorama_seurat$seurat_clusters)
rownames(scanorama_batch_metadata) = rownames(scanorama_seurat@meta.data)
scanorama_batch_metadata$method = rep("Scanorama", length(rownames(scanorama_batch_metadata)))

# Calculate iLISI scores
# Mean ilisi = 1.311; Median ilisi = 1.137
scanorama_ilisi_res = lisi::compute_lisi(scanorama_pca_embeddings, scanorama_batch_metadata, c("batch_label"))
scanorama_batch_metadata$ilisi = scanorama_ilisi_res$batch_label

# plot iLISI scores
scanorama_batch_metadata %>%
  group_by(cluster) %>%
  summarize(mean_ilisi = mean(ilisi)) %>%
  ggplot(aes(x=cluster, y=mean_ilisi)) + 
  geom_col() + 
  #facet_wrap(~cluster, scales = "free_y") + 
  custom_theme

scanorama_batch_metadata %>%
  ggplot(aes(x=ilisi)) + 
  geom_density() + 
  facet_wrap(~cluster, scales = "free_y")

scanorama_batch_metadata %>%
  ggplot(aes(x=method, y=ilisi)) + 
  geom_boxplot(aes(x=method, y=ilisi))

# Calculate cLISI scores
# mean cLISI = 1.192; median cLISI = 1.013
scanorama_clisi_res = lisi::compute_lisi(scanorama_pca_embeddings, scanorama_batch_metadata,
                                         c("cluster"))
scanorama_batch_metadata$clisi = scanorama_clisi_res$cluster
mean(scanorama_clisi_res$cluster)
median(scanorama_clisi_res$cluster)

# plot cLISI scores
scanorama_batch_metadata %>%
  ggplot(aes(x=clisi)) + 
  geom_density() + 
  custom_theme




#######################################################
# Plot of mean iLISI and cLISI values for each approach
lisi_scores = data.frame(method = c("CCA_MNN", "rPCA", "Scanorama", "Harmony"),
                         mean_lisi = c(mean(cca_batch_metadata$clisi),
                                       mean(rpca_batch_metadata$clisi),
                                       mean(scanorama_batch_metadata$clisi),
                                       mean(harmony_batch_metadata$clisi)))
#lisi_scores$method = factor(lisi_scores$method, levels = c("CCA_MNN", "rPCA", "Scanorama", "Harmony"))
lisi_scores$method = factor(lisi_scores$method, levels = c("rPCA", "Scanorama", "CCA_MNN", "Harmony"))

mean_lisi_scores = lisi_scores %>%
  ggplot(aes(x=method, y=mean_lisi, fill=method)) + 
  geom_col(width = 0.8) + 
  custom_theme + 
  ylab("Mean cLISI") + 
  xlab("Method") + 
  #scale_fill_manual(values = c("#E64B35FF", "#3C5488FF", "#4DBBD5FF", "#00A087FF")) + 
  scale_fill_manual(values = c("#3C5488FF", "#4DBBD5FF", "#E64B35FF", "#00A087FF")) + 
  theme(aspect.ratio = 1.3)
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure1/SuppFig1d_integration_mean_iLISI_scores.pdf",
       plot = mean_lisi_scores, width = 8, height = 8)

# Merge dataframes with lisi scores for all methods
merged_lisi_df = data.table::rbindlist(list(cca_batch_metadata,
                                            rpca_batch_metadata,
                                            scanorama_batch_metadata,
                                            harmony_batch_metadata), use.names = TRUE)


saveRDS(merged_lisi_df,
        "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA_Harmony_CCA_Scanorama_clISI_scores_per_cell.rds")
merged_lisi_df = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA_Harmony_CCA_Scanorama_clISI_scores_per_cell.rds")

mean_clisi_scores = merged_lisi_df %>%
  group_by(method) %>%
  summarize(mean_cLISI = mean(clisi)) %>%
  ggplot(aes(x=factor(method, 
                      levels = c("rPCA", "Scanorama", "CCA_MNN", "Harmony")), 
             y=mean_cLISI, fill=method)) + 
  geom_col(width = 0.8) + 
  xlab("Method") + 
  ylab("Mean cLISI") + 
  custom_theme + 
  scale_fill_manual(values = c("#3C5488FF", "#4DBBD5FF", "#E64B35FF", "#00A087FF")) + 
  theme(aspect.ratio = 1.3)
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure1/SuppFig1e_integration_mean_cLISI_scores.pdf",
       plot = mean_clisi_scores, width = 8, height = 8)

# Try a boxplot instead
# merged_lisi_scores = merged_lisi_df %>%
#   ggplot(aes(x=method, y=clisi, fill=method)) + 
#   geom_boxplot(aes(x=method, y=clisi), outlier.size = 0.5) +
#   #stat_summary(fun = mean, geom = "point", shape=20, size=10, color="black", fill="black") + 
#   scale_fill_manual(values = c("#E64B35FF", "#3C5488FF", "#4DBBD5FF", "#00A087FF")) + 
#   ylab("iLISI") + 
#   xlab("Method") + 
#   ylim(0.7, 3) + 
#   custom_theme + 
#   theme(aspect.ratio = 1.2)

# Try a density plot instead 
# merged_lisi_df %>%
#   ggplot(aes(x=clisi, color=method)) + 
#   geom_density(alpha=0.5) +
#   xlab("cLISI") + 
#   xlim(1, 3) + 
#   ylim(0, 5) + 
#   custom_theme + 
#   scale_color_manual(values = c("#E64B35FF", "#3C5488FF", "#4DBBD5FF", "#00A087FF")) + 
#   theme(aspect.ratio = 1.2,
#         legend.position = "bottom")


