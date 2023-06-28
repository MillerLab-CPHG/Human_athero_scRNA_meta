library(kBET)
library(FNN)
library(lisi)
library(cluster)
library(Seurat)
library(scRNAutils)
library(stats)
library(DescTools)
library(data.table)
library(tidyverse)
library(ggsci)


# Source scRNA custom utils script
#source("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/R_scripts/Utils/scRNA_processing_utils.R")


############################################################################
# We should load PCA embeddings at once so things don't get too repetitive #
############################################################################

embeddings_dir = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/int_benchmark/dim_red_embeddings/"
metadata_dir = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/int_benchmark/metadata/"

#######################################################
# Get embeddings and metadata from integrated objects #
#######################################################

#########
# CCA
cca_pca_list = cca_int_sct@reductions$pca
saveRDS(cca_pca_list,
        paste(embeddings_dir, "cca_pca_dim_red_list.rds", sep=""))

cca_metadata = cca_int_sct@meta.data
saveRDS(cca_metadata,
        paste(metadata_dir, "cca_int_metadata.rds", sep=""))

# Load data
cca_pca_list = read_rds(paste(embeddings_dir, "cca_pca_dim_red_list.rds", sep=""))
cca_metadata = read_rds(paste(metadata_dir, "cca_int_metadata.rds", sep=""))

#########
# rPCA
rpca_pca_list = rpca_int_sct@reductions$pca
saveRDS(rpca_pca_list,
        paste(embeddings_dir, "rpca_pca_dim_red_list.rds", sep=""))

rpca_metadata = rpca_int_sct@meta.data
saveRDS(rpca_metadata,
        paste(metadata_dir, "rpca_int_metadata.rds", sep=""))

# Load data
rpca_pca_list = read_rds(paste(embeddings_dir, "rpca_pca_dim_red_list.rds", sep=""))
rpca_metadata = read_rds(paste(metadata_dir, "rpca_int_metadata.rds", sep=""))

##########
# Harmony
harmony_list = harmony_seurat_int@reductions$harmony
saveRDS(harmony_list,
        paste(embeddings_dir, "harmony_dim_red_list.rds", sep=""))

harmony_metadata = harmony_seurat_int@meta.data
saveRDS(harmony_metadata,
        paste(metadata_dir, "harmony_int_metadata.rds", sep=""))

# Load data
harmony_list = read_rds(paste(embeddings_dir, "harmony_dim_red_list.rds", sep=""))
harmony_metadata = read_rds(paste(metadata_dir, "harmony_int_metadata.rds", sep=""))

############
# Scanorama
scanorama_pca_list = scanorama_seurat@reductions$pca
saveRDS(scanorama_pca_list,
        paste(embeddings_dir, "scanorama_pca_dim_red_list.rds", sep=""))

scanorama_metadata = scanorama_seurat@meta.data
saveRDS(scanorama_metadata,
        paste(metadata_dir, "scanorama_int_metadata.rds", sep=""))

# Load data
scanorama_pca_list = read_rds(paste(embeddings_dir, "scanorama_pca_dim_red_list.rds", sep=""))
scanorama_metadata = read_rds(paste(metadata_dir, "scanorama_int_metadata.rds", sep=""))

#######################################################################
# Load PCA embeddings (First 30 PCs to speed up computations)
cca_pca_embeddings = cca_pca_list@cell.embeddings[, 1:30]
rpca_pca_embeddings = rpca_pca_list@cell.embeddings[, 1:30]
harmony_embeddings = harmony_list@cell.embeddings[, 1:30]
scanorama_pca_embeddings = scanorama_pca_list@cell.embeddings[, 1:30]

###############################################
# Calculate Sihouette scores for each method  #
###############################################

# Silhouette scores measure clustering purity by taking into account cohesion and separation of cells
# in clusters. In other words, how similar a cells is to its own cluster compared to other clusters.
# This gives a score in the range of −1 to +1, where a higher score indicates higher performance.
# The silhouette coefficient captures elements of both sample mixing and local structure.
# Here we'll calculate sil scores across a range of clustering resolutions to control for the 
# Seurat clustering granularity parameter. For this calculations, we'll source the wrapper from the utils script

# Define a range of resolutions to test
res = seq(0.8, 1.8, by=0.1)

###############################################
# CCA

# Test a few resolution parameters and see which ones produce the higher silhouette scores
# Figure out how many clusters each res yields and what the overall sil scores are
cca_scores_list = list()
cca_dist_matrix = dist(x = cca_pca_embeddings)
saveRDS(cca_dist_matrix, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/CCA/cca_distance_matrix.rds")
cca_dist_matrix = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/CCA/cca_distance_matrix.rds")

for (i in seq_along(res)) {
  # Cluster data
  cca_int_sct = FindClusters(cca_int_sct, 
                             resolution = res[i])
  
  # Calculate silhouette scores
  cca_scores_list[[i]] = calc_sil_scores(cca_int_sct, 
                                         cca_dist_matrix)
  names(cca_scores_list)[[i]] = paste("res", 
                                      as.character(res[i]), 
                                      sep = "_")
}

# Create a single df with all the sil scores for the tested resolutions
cca_resolutions_sil_df = rbindlist(cca_scores_list, idcol = TRUE)
cca_resolutions_sil_df$method = "CCA"
saveRDS(cca_resolutions_sil_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/CCA/CCA_sil_scores_per_clustering_resolution.rds")
cca_resolutions_sil_df =  read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/CCA/CCA_sil_scores_per_clustering_resolution.rds")


################################
# rPCA

# Figure out how many clusters each res yields and what the overall sil scores are
rpca_scores_list = list()
dist_matrix = dist(x = rpca_pca_embeddings)
saveRDS(dist_matrix, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/distance_matrix.rds")
dist_matrix = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/distance_matrix.rds")

for (i in seq_along(res)) {
  # Cluster data
  rpca_int_sct = FindClusters(rpca_int_sct, 
                              resolution = res[i])
  
  # Calculate silhouette scores
  rpca_scores_list[[i]] = calc_sil_scores(rpca_int_sct, 
                                          dist_matrix)
  names(rpca_scores_list)[[i]] = paste("res", 
                                       as.character(res[i]), 
                                       sep = "_")
}

# Create a single df with all the sil scores for the tested resolutions
rpca_resolutions_sil_df = rbindlist(rpca_scores_list, idcol = TRUE)
rpca_resolutions_sil_df$method = "rPCA"
saveRDS(resolutions_sil_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/sil_scores_per_clustering_resolution.rds")
rpca_resolutions_sil_df = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_int/sil_scores_per_clustering_resolution.rds")


############################
# Harmony

# Test a few resolution parameters and see which ones produce the higher silhouette scores
# Figure out how many clusters each res yields and what the overall sil scores are
h_scores_list = list()
h_dist_matrix = dist(x = harmony_embeddings)
saveRDS(h_dist_matrix, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Harmony/distance_matrix.rds")
h_dist_matrix = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Harmony/distance_matrix.rds")

for (i in seq_along(res)) {
  # Cluster data
  harmony_seurat_int = FindClusters(harmony_seurat_int, 
                                    resolution = res[i])
  
  # Calculate silhouette scores
  h_scores_list[[i]] = calc_sil_scores(harmony_seurat_int, 
                                       h_dist_matrix)
  names(h_scores_list)[[i]] = paste("res", 
                                    as.character(res[i]), 
                                    sep = "_")
}

# Create a single df with all the sil scores for the tested resolutions
# Max res
h_resolutions_sil_df = rbindlist(h_scores_list, idcol = TRUE)
h_resolutions_sil_df$method = "harmony"
#res_col_names = paste("res", res, sep="_")
#colnames(h_resolutions_sil_df) = res_col_names
saveRDS(h_resolutions_sil_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Harmony/harmony_sil_scores_per_clustering_resolution.rds")
h_resolutions_sil_df = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Harmony/harmony_sil_scores_per_clustering_resolution.rds")


################################################
# Scanorama

# Test a few resolution parameters and see which ones produce the higher silhouette scores
# Figure out how many clusters each res yields and what the overall sil scores are
scanorama_scores_list = list()
scanorama_dist_matrix = dist(x = scanorama_pca_embeddings)
saveRDS(scanorama_dist_matrix, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Scanorama/scanorama_distance_matrix.rds")
scanorama_dist_matrix = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Scanorama/scanorama_distance_matrix.rds")

for (i in seq_along(res)) {
  # Cluster data
  scanorama_seurat = FindClusters(scanorama_seurat, 
                                  resolution = res[i])
  
  # Calculate silhouette scores
  scanorama_scores_list[[i]] = calc_sil_scores(scanorama_seurat, 
                                               scanorama_dist_matrix)
  names(scanorama_scores_list)[[i]] = paste("res", 
                                            as.character(scanorama_res[i]), 
                                            sep = "_")
}

# Create a single df with all the sil scores for the tested resolutions
scanorama_resolutions_sil_df = rbindlist(scanorama_scores_list, idcol = TRUE)
scanorama_resolutions_sil_df$method = "Scanorama"
dim(scanorama_resolutions_sil_df)
saveRDS(scanorama_resolutions_sil_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/Scanorama/Scanorama_sil_scores_per_clustering_resolution.rds")
rpca_res_sil_df = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_int/sil_scores_per_clustering_resolution.rds")


########################################################
# Compare rPCA vs CCA vs Harmony vs Scanorama sil scores 
res_sil_scores_df = rbind.data.frame(rpca_res_sil_df, 
                                     h_resolutions_sil_df, 
                                     cca_resolutions_sil_df, 
                                     scanorama_resolutions_sil_df)
saveRDS(res_sil_scores_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA_Harmony_CCA_Scanorama_silhouette_scores_per_clustering_res.rds")
res_sil_scores_df = readRDS("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA_Harmony_CCA_Scanorama_silhouette_scores_per_clustering_res.rds")

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

# Wrapper function for kBET calculations
calc_kbet = function(dim_red_embeddings, batch_labels,
                     k_range, method_id) {
  kbet_list = list()
  for (i in seq_along(k_range)) { 
    msg = paste("Calculating kBET for", k_range[i], 
                "neighbors, please bear with us...")
    message(msg)
    batch_est = kBET(df=dim_red_embeddings, 
                     batch=batch_labels, 
                     k0 = k_range[i], 
                     knn = NULL,
                     do.pca=FALSE, 
                     heuristic = FALSE,
                     plot = FALSE,
                     adapt = FALSE, 
                     verbose = TRUE)
    
    plotting_df = data.frame(class = rep(c("observed", "expected"), 
                                         each=length(batch_est$stats$kBET.observed)),
                             data = c(batch_est$stats$kBET.observed,
                                      batch_est$stats$kBET.expected))
    plotting_df$k = paste("k_", k_range[i], sep="")
    plotting_df$method = method_id
    kbet_list[[i]] = plotting_df
    names(kbet_list)[i] = paste("kBET_", k_range[i], "_neighbors",
                                sep = "")
    
  }
  return(kbet_list)
}

# Define a k range
k_range = c(10, 25, 50, 100, 500, 1000)

###############################
# CCA

# CCA  Get batch labels 
cca_batch_labels = as.factor(cca_metadata$study)

# Calculate kBET
# Standard implementation of kBET performs a k-nearest neighbors search (if knn=NULL)
# Let's try to run the k-nearest neighbors search in PCA space separately.

# Neighborhood size suggested by kBET authors
k0 = floor(mean(table(harmony_batch_labels))) / 4
k_range = c(k_range, k0)

# Let's try 10 neighbors just for fun. Avg obs rejection rate=0.3 
# Observed rejection rate for k0=5105 is 1. Finding knn takes ~25 mins for this k. This value doesn't make any sense. 
# Maybe let's try a range of k like 10, 25, 50, 100
# It takes 39 secs to get the knn for k=10
#k_range = c(10, 25, 50, 100)

# Calculate kBET results for CCA
cca_kbet_list = calc_kbet(cca_pca_embeddings, 
                          batch_labels=cca_batch_labels,
                          k_range=k_range,
                          method_id="CCA_MNN")

cca_kbet_merged_df = rbindlist(cca_kbet_list)


###############
# rPCA

# Get batch labels 
rpca_batch_labels = as.factor(rpca_metadata$study)

# Calculate kBET metric across the defined range of k (10, 25, 50 and 100)
rpca_kbet_list = calc_kbet(rpca_pca_embeddings, 
                           batch_labels=rpca_batch_labels,
                           k_range=k_range, 
                           method_id="rPCA")

rpca_kbet_merged_df = rbindlist(rpca_kbet_list)

###############
# Harmony

# Get batch labels 
harmony_batch_labels = as.factor(harmony_metadata$study)

# Calculate kBET metric across the defined range of k (10, 25, 50 and 100)
harmony_kbet_list = calc_kbet(harmony_embeddings, 
                              batch_labels=harmony_batch_labels,
                              k_range=k_range, 
                              method_id="Harmony")

harmony_kbet_merged_df = rbindlist(harmony_kbet_list)

###############
# Scanorama

# Get batch labels 
scanorama_batch_labels = as.factor(scanorama_metadata$study)

# Calculate kBET metric across the defined range of k (10, 25, 50 and 100)
scanorama_kbet_list = calc_kbet(scanorama_pca_embeddings, 
                              batch_labels=scanorama_batch_labels,
                              k_range=k_range, 
                              method_id="Scanorama")

scanorama_kbet_merged_df = rbindlist(scanorama_kbet_list)

###################################################################
# Compare CCA, rPCA, Harmony and Scanorama observed rejection rates
merged_res = rbind.data.frame(cca_kbet_merged_df, 
                              rpca_kbet_merged_df,
                              harmony_kbet_merged_df,
                              scanorama_kbet_merged_df)
saveRDS(merged_res,
        "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/int_benchmark/kBET_results/kBET_rejection_rates_df.rds")

merged_res = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/int_benchmark/kBET_results/kBET_rejection_rates_df.rds")

merged_res %>%
  filter(class == "observed") %>%
  group_by(method, k) %>%
  summarize(mean_obs_rejection_rate = mean(data)) %>%
  ggplot(aes(x=reorder(k, mean_obs_rejection_rate), 
             mean_obs_rejection_rate, 
             color=method)) + 
  geom_line(aes(x=reorder(k, mean_obs_rejection_rate),
                y=mean_obs_rejection_rate, group=method), 
            linetype="dashed") + 
  geom_point(size=1.5) +
  #geom_col(position = "dodge") + 
  labs(x="Neighborhood size (k)",
       y="Mean kBET observed rejection rate") + 
  custom_theme() + 
  theme(aspect.ratio = 1.3,
        legend.position = "right",
        axis.text.x = element_text(angle=45, hjust=1)) + 
  miller_discrete_scale()


# Calculate the area under the curve (AUC) for each method
methods = c("CCA_MNN", "rPCA", "Harmony", "Scanorama")
auc_list = list()

for (i in seq_along(methods)) { 
  
  sorted_rej_rates = merged_res %>% 
    filter(class == "observed" & method == methods[i]) %>%
    group_by(method, k) %>%
    summarize(mean_obs_rejection_rate = mean(data)) %>%
    arrange(mean_obs_rejection_rate) %>%
    mutate(k_index = seq(1, 7, by=1))
  
  auc = DescTools::AUC(x=sorted_rej_rates$k_index,
                       y=sorted_rej_rates$mean_obs_rejection_rate,
                       method="spline")
  auc_list[[i]] = auc
  names(auc_list)[i] = methods[i]

  }







##################################################
# Principal Components Regression (PCRegression) #
##################################################


###############################################
# Write a loop to do all of the methods at once
cca = list(dim_red=cca_pca_list, metadata=cca_metadata)
rpca = list(dim_red=rpca_pca_list, metadata=rpca_metadata)
harmony = list(dim_red=harmony_list, metadata=harmony_metadata)
scanorama = list(dim_red=scanorama_pca_list, metadata=scanorama_metadata)

# Define a vector with the names of the methods to test
methods = c("cca", "rpca", "harmony", "scanorama")

# Create a list to store results
pcr_list = list()

# Looks like doing both PC regression and calculating the batch sil score
# might be too intensive in a single loop.
# Calculate metric for methods defined in the vector
for (i in methods) { 
  
  # Prep inputs for linear regression
  batch_cor_method = get(i)

  method_reduc_data = list(x=batch_cor_method$dim_red@cell.embeddings,
                         sdev=batch_cor_method$dim_red@stdev)
  
  # Get batch labels
  batch_labels = as.factor(batch_cor_method$metadata$study)
  
  # Carry out PC regression using the first 30 PCs
  method_pc_reg = kBET::pcRegression(pca.data=method_reduc_data,
                                     batch=batch_labels,
                                     n_top=30)
  
  # Produce index to store PC regression results
  idx = which(methods == i)
  pcr_list[[idx]] =  method_pc_reg
  names(pcr_list)[idx] = i
  
}

# Prep a df to plot results
r2var_df = data.frame(method=c("CCA", "rPCA", "Harmony", "Scanorama"),
                      R2Var=c(pcr_list$cca$R2Var,
                              pcr_list$rpca$R2Var,
                              pcr_list$harmony$R2Var,
                              pcr_list$scanorama$R2Var))
r2var_df %>%
  mutate(method=fct_reorder(method, 
                            desc(-log10(R2Var)))) %>%
  ggplot(aes(x=method, y=-log10(R2Var), 
             fill=method)) + 
  geom_col(width = 0.6) + 
  xlab("Method") + 
  ylab("-Log10(R2Var") +
  custom_theme() + 
  miller_discrete_scale(style = "bars") + 
  theme(aspect.ratio = 1.3,
        legend.position = "none",
        axis.text.x = element_text(angle = 45,
                                   hjust=1))

############################################
# Run a quick test for the batch sil score #
############################################

# # To speed up calculations we will only carry out the calculation in 3 PCs
# batch_sil_new = function(pca_data, batch, nPCs=2) {
#   dd = stats::dist(pca_data$x[, seq_len(nPCs)])
#   abs(summary(cluster::silhouette(as.numeric(batch), 
#                               dist=dd))$avg.width)
# }
# 
# method_reduc_obj = cca_int_sct@reductions$pca
# method_reduc_data = list(x=method_reduc_obj@cell.embeddings,
#                          sdev=method_reduc_obj@stdev)
# 
# # Get batch labels
# batch_labels = as.factor(cca_int_sct@meta.data$study)
# 
# # Calculate Avg Silhouette Width
# batch_asw = batch_sil_new(pca_data=method_reduc_data,
#                           batch = batch_labels,
#                           nPCs=3)
# 
# # Create a list to store results
# batch_asw_list = list()

# # Calculate metric for methods defined in the vector
# for (i in methods_test) { 
#   
#   # Prep inputs for linear regression
#   batch_cor_method = get(i)
#   if (startsWith(i, "harmony")) {
#     method_reduc_obj = batch_cor_method@reductions$harmony
#   } else {
#     method_reduc_obj = batch_cor_method@reductions$pca
#   }
#   
#   method_reduc_data = list(x=method_reduc_obj@cell.embeddings,
#                            sdev=method_reduc_obj@stdev)
#   
#   # Get batch labels
#   batch_labels = as.factor(batch_cor_method@meta.data$study)
# 
#   # Calculate Avg Silhouette Width
#   batch_asw = batch_sil_new(pca_data=method_reduc_data,
#                             batch = batch_labels,
#                             nPCs=3)
#   
#   # Produce index to store PC regression results
#   idx = which(methods_test == i)
# 
#   # Store Batch Avg Sil width scores
#   batch_asw_list[[idx]] = batch_asw
#   names(batch_asw_list)[idx] = i
# }


#########################
# Calculate LISI scores #
#########################

#########################
# CCA
cca_batch_metadata = data.frame(batch_label=cca_metadata$study, 
                                cluster=cca_metadata$integrated_snn_res.1.3)
rownames(cca_batch_metadata) = rownames(cca_metadata)
cca_batch_metadata$method = rep("CCA_MNN", 
                                length(rownames(cca_batch_metadata)))


# Calculate iLISI
# Mean iLISI = 1.796; Median iLISI = 1.788
cca_ilisi_res = lisi::compute_lisi(cca_pca_embeddings, 
                                   cca_batch_metadata, 
                                   c("batch_label"))
cca_batch_metadata$ilisi = cca_ilisi_res$batch_label

# Calculate cLISI
# mean cLISI = 1.22; median clISI = 1.02
cca_clisi_res = lisi::compute_lisi(cca_pca_embeddings, 
                                   cca_batch_metadata, 
                                   c("cluster"))
cca_batch_metadata$clisi = cca_clisi_res$cluster


##########################
# rPCA
rpca_batch_metadata = data.frame(batch_label = rpca_metadata$study,
                                 cluster=rpca_metadata$integrated_snn_res.1.3)
rownames(rpca_batch_metadata) = rownames(rpca_metadata)
rpca_batch_metadata$method = rep("rPCA", 
                                 length(rownames(rpca_batch_metadata)))

# Calculate iLISI
# Mean ilisi = 1.61; Median ilisi = 1.484
rpca_ilisi_res = lisi::compute_lisi(rpca_pca_embeddings, 
                                    rpca_batch_metadata, 
                                    c("batch_label"))
rpca_batch_metadata$ilisi = rpca_ilisi_res$batch_label

# Calculate cLISI
# mean cLISI = 1.159; median cLISI = 1.002
rpca_clisi_res = lisi::compute_lisi(rpca_pca_embeddings, 
                                    rpca_batch_metadata, 
                                    c("cluster"))
rpca_batch_metadata$clisi = rpca_clisi_res$cluster


############################
# Harmony

harmony_batch_metadata = data.frame(batch_label=harmony_metadata$study,
                                    cluster=harmony_metadata$seurat_clusters)
rownames(harmony_batch_metadata) = rownames(harmony_metadata)
harmony_batch_metadata$method = rep("Harmony", 
                                    length(rownames(harmony_batch_metadata)))

# Calculate iLISI scores
# Mean ilisi = 1.55;  median ilisi = 1.402
harmony_ilisi_res = lisi::compute_lisi(harmony_embeddings, 
                                       harmony_batch_metadata, 
                                       c("batch_label"))
harmony_batch_metadata$ilisi = harmony_ilisi_res$batch_label


# Calculate cLISI scores
# mean cLISI = 1.286; median cLISI = 1.062
harmony_clisi_res = lisi::compute_lisi(harmony_embeddings, 
                                       harmony_batch_metadata, 
                                       c("cluster"))
harmony_batch_metadata$clisi = harmony_clisi_res$cluster


##############################
# Scanorama

scanorama_batch_metadata = data.frame(batch_label = scanorama_metadata$study,
                                      cluster = scanorama_metadata$seurat_clusters)
rownames(scanorama_batch_metadata) = rownames(scanorama_metadata)
scanorama_batch_metadata$method = rep("Scanorama", length(rownames(scanorama_batch_metadata)))

# Calculate iLISI scores
# Mean ilisi = 1.311; Median ilisi = 1.137
scanorama_ilisi_res = lisi::compute_lisi(scanorama_pca_embeddings, 
                                         scanorama_batch_metadata, 
                                         c("batch_label"))
scanorama_batch_metadata$ilisi = scanorama_ilisi_res$batch_label


# Calculate cLISI scores
# mean cLISI = 1.192; median cLISI = 1.013
scanorama_clisi_res = lisi::compute_lisi(scanorama_pca_embeddings, 
                                         scanorama_batch_metadata,
                                         c("cluster"))
scanorama_batch_metadata$clisi = scanorama_clisi_res$cluster


#######################################################
# Plot of mean iLISI and cLISI values for each approach

# Merge dataframes with lisi scores for all methods
merged_lisi_df = data.table::rbindlist(list(cca_batch_metadata,
                                            rpca_batch_metadata,
                                            scanorama_batch_metadata,
                                            harmony_batch_metadata), use.names = TRUE)

# Save df with LISI scores (iLISI and cLISI)
saveRDS(merged_lisi_df,
        "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA_Harmony_CCA_Scanorama_LISI_scores_per_cell.rds")
merged_lisi_df = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA_Harmony_CCA_Scanorama_clISI_scores_per_cell.rds")

# Plot iLISI scores
mean_ilisi_scores = merged_lisi_df %>%
  group_by(method) %>%
  summarize(mean_iLISI = mean(ilisi)) %>%
  mutate(method = fct_reorder(method, desc(mean_iLISI))) %>%
  ggplot(aes(x=method, y=mean_iLISI, fill=method)) + 
  geom_col(width = 0.6) + 
  xlab("Method") + 
  ylab("Mean iLISI") + 
  custom_theme() + 
  scale_fill_manual(values= c("#E64B35FF", "#00A087FF", 
                              "#4DBBD5FF", "#3C5488FF")) + 
  theme(aspect.ratio = 1.3,
        legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1))

ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure1/SuppFig1e_integration_mean_cLISI_scores.pdf",
       plot = mean_clisi_scores, width = 8, height = 8)


# Plot cLISI scores
mean_clisi_scores = merged_lisi_df %>%
  group_by(method) %>%
  summarize(mean_cLISI = mean(clisi)) %>%
  mutate(method = fct_reorder(method, mean_cLISI)) %>%
  ggplot(aes(x=method, y=mean_cLISI, fill=method)) + 
  geom_col(width = 0.6) + 
  xlab("Method") + 
  ylab("Mean cLISI") + 
  custom_theme() + 
  scale_fill_manual(values= c("#E64B35FF", "#00A087FF", 
                              "#4DBBD5FF", "#3C5488FF")) + 
  theme(aspect.ratio = 1.3,
        legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1))

ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure1/SuppFig1e_integration_mean_cLISI_scores.pdf",
       plot = mean_clisi_scores, width = 8, height = 8)






