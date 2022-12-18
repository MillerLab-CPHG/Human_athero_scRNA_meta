library(Seurat)
library(CellChat)
library(liana)
library(tidyverse)
library(data.table)



# This script will run cell communication analyses

####################################################
# Set parallel run                                 #
future::plan("multiprocess", workers=4)            #
options(future.globals.maxSize = 8000 * 1024^2)    #
####################################################

# We need to load the normalized scRNA matrix and metadata
sct_norm_counts = rpca_int_sct_v3@assays$SCT@data

# Load the metadata for level2 annotations
rpca_v3_meta = rpca_int_sct_v3@meta.data
meta_df = data.frame(labels=rpca_v3_meta$level2_annotations, row.names = rownames(rpca_v3_meta))
unique(meta_df$labels)

# Create cellchat object
rpca_v3_cellchat = createCellChat(object = sct_norm_counts, meta = meta_df, group.by = "labels")

###################################################
# Set the ligand-receptor interaction database    #
cellchatdb_hs = CellChatDB.human                  #
showDatabaseCategory(cellchatdb_hs)               #
###################################################


# Glance at structure of the db
glimpse(cellchatdb_hs)

# Set the used db withiin the object
rpca_v3_cellchat@DB = cellchatdb_hs


# Pre-process the expression data for downstream analyses
rpca_v3_cellchat = subsetData(rpca_v3_cellchat)

#rpca_v3_cellchat = identifyOverExpressedGenes(rpca_v3_cellchat)
rpca_v3_cellchat = identifyOverExpressedInteractions(rpca_v3_cellchat)

# Compute communication probabilities
rpca_v3_cellchat = computeCommunProb(rpca_v3_cellchat, population.size = TRUE)
rpca_v3_cellchat = computeCommunProbPathway(rpca_v3_cellchat)

# Calculate aggregated cell commnication network
rpca_v3_cellchat = aggregateNet(rpca_v3_cellchat)

# Save cell communication object with level2 annotations
saveRDS(rpca_v3_cellchat, 
        "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cellchat_objects/cellchat_obj_level2_annotations.rds")

# Load cellchat object for level2 annotations
rpca_v3_cellchat_level2 = readRDS("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cellchat_objects/cellchat_obj_level2_annotations.rds")

# Visualize aggregated cell communication network 
group_size = as.numeric(table(rpca_v3_cellchat@idents))
netVisual_circle(rpca_v3_cellchat_level2@net$weight, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weight/strength")

# Compare edge weights 
mat = rpca_v3_cellchat@net$weight
par(mfrow=c(8,4))
for (i in 1:nrow(mat)) {
  mat2 = matrix(0, nrow=nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] = mat[i, ]
  netVisual_circle(mat2, vertex.weight = group_size, weight.scale = TRUE, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Subset cell communication

# Pathways to follow up from level1 annotations: TNF, GRN, TWEAK, PDGF, GALECTIN, CD46
# Pathways specific to level2 annotations: SPP1, RESISTIN, CD6
df_net_smc_level2 = subsetCommunication(rpca_v3_cellchat_level2, sources.use = c("Foamy_Mac1", "Foamy_Mac2", "Monocytes", "Inflammatory_Mac", "Fibroblast", "Myofibroblast",
                                                                      "Activated_NK", "CTL_CD8_early_activated/mem", "CTL_CD8_terminally_diff",
                                                                      "Monocytes/DC", "NAMPT_Neutrophils", "Treg", "cDC", "Tissue_resident_Mac", "Phagocytosis_Mac"),
                             targets.use = c("SMC"))
netVisual_aggregate(rpca_v3_cellchat_level2, signaling = c("GALECTIN"), layout = "circle", thresh = 0.01, 
                    sources.use = c("Inflammatory_Mac", "Foamy_Mac1", "Foamy_Mac2", 'Monocytes', "Monocytes/DC",
                                    "NAMPT_Neutrophils", "Tissue_resident_Mac", "cDC", "Phagocytosis_Mac", "Proliferating_myeloid") ,targets.use = c("SMC"))



df_net_endo_level2 = subsetCommunication(rpca_v3_cellchat_level2, sources.use = c("Foamy_Mac2", "Monocytes", "Inflammatory_Mac", "Fibroblast", "Myofibroblast",
                                                                           "Activated_NK", "CTL_CD8_early_activated/mem", "CTL_CD8_terminally_diff",
                                                                           "Monocytes/DC", "NAMPT_Neutrophils", "Treg", "cDC", "Proliferating_myeloid"),
                                  targets.use = c("Inflammatory_EC", "Intimal_EC", "EndoMT_EC"))
netVisual_aggregate(rpca_v3_cellchat_level2, signaling = c("VCAM"), layout = "circle", thresh = 0.01,
                    sources.use = c("Inflammatory_Mac", "Foamy_Mac1", "Foamy_Mac2", 'Monocytes', "Monocytes/DC",
                                    "NAMPT_Neutrophils", "Tissue_resident_Mac", "cDC", "Phagocytosis_Mac", "Proliferating_myeloid"),
                    targets.use = c("Inflammatory_EC", "Intimal_EC", "EndoMT_EC"))


# Extract significant interactions
extractEnrichedLR(rpca_v3_cellchat_level2, signaling = c("ITGB2"), geneLR.return = TRUE)
netAnalysis_contribution(rpca_v3_cellchat_level2, signaling = c("RESISTIN"))



###################################################################
# Run communication analyses for level1 annotations
meta_df_level1 = data.frame(labels=rpca_v3_meta$level1_annotations, row.names = rownames(rpca_v3_meta))
unique(meta_df_level1$labels)

#Create cellchat object
rpca_v3_cellchat_level1 = createCellChat(object = sct_norm_counts, meta = meta_df_level1, group.by = "labels")

# Set the used db withiin the object
rpca_v3_cellchat_level1@DB = cellchatdb_hs
 
  
# Pre-process the expression data for downstream analyses
rpca_v3_cellchat_level1 = subsetData(rpca_v3_cellchat_level1)
rpca_v3_cellchat_level1 = identifyOverExpressedGenes(rpca_v3_cellchat_level1)
rpca_v3_cellchat_level1 = identifyOverExpressedInteractions(rpca_v3_cellchat_level1)

saveRDS(rpca_v3_cellchat_level1,
        "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cellchat_objects/cellchat_obj_level1_annotations.rds")

# Load cellchat object for level 1 annnotations   
rpca_v3_cellchat_level1 = readRDS("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cellchat_objects/cellchat_obj_level1_annotations.rds")

# Still need to compute communication probabilities for level1 annotations
rpca_v3_cellchat_level1 = computeCommunProb(rpca_v3_cellchat_level1, population.size = TRUE)
rpca_v3_cellchat_level1 = computeCommunProbPathway(rpca_v3_cellchat_level1)

#Calculate aggregated cell commnication network
rpca_v3_cellchat_level1 = aggregateNet(rpca_v3_cellchat_level1)
  
# Visualize aggregated cell communication network 
group_size = as.numeric(table(rpca_v3_cellchat_level1@idents))
netVisual_circle(rpca_v3_cellchat_level1@net$count, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")

netVisual_circle(rpca_v3_cellchat_level1@net$weight, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, title.name = "Interactions weight/strength")

# Subset cell communication network to focus on incoming signals to SMCs
# Interesting pathways to follow up for SMCs: PDGF, TNF, GRN, GALECTIN, FGF, MK, CD46 (COMPLEMENT)
# TWEAK (TNF-related)
df_net_smc_level1 = subsetCommunication(rpca_v3_cellchat_level1, sources.use = c("Macrophage", "T_NK", "Fibroblast"),
                             targets.use = c("SMC"))

# Interesting pathways to follow up for endothelial cells: VEGF, CXCL, CCL, TNF, VISFATIN (induces endothelial dysfunction), ANGPT, MK
# ITGB2 (could look at the weights of the interactions here), VCAM
df_net_endo_level1 = subsetCommunication(rpca_v3_cellchat_level1, sources.use = c("Macrophage", "T_NK", "Fibroblast"),
                             targets.use = c("Endothelial"))

netVisual_aggregate(rpca_v3_cellchat_level1, signaling = c("TWEAK"), layout = "circle" )

# Check weights of ligand-receptor interactions
netAnalysis_contribution(rpca_v3_cellchat_level1, signaling = "PDGF")

# Extract significant interactions
extractEnrichedLR(rpca_v3_cellchat_level1, signaling = c("GALECTIN"), geneLR.return = TRUE)
netVisual_bubble(rpca_v3_cellchat_level1, sources.use = c("Macrophage", "T_NK", "Fibroblast"), targets.use = c("SMC", "Endothelial"),
                 remove.isolate = TRUE)

###########################################################
# Run Cellchat comparing lesion vs non-lesion libraries   #
###########################################################

# To run this analysis, we need to create a separate cellchat objects for each condition.

# Load metadata and get cell barcodes from lesion and non-lesion libraries to subset counts matrix
metadata = rpca_int_sct_v3@meta.data

# There's a total of 59691 cell barcodes from lesion libraries
lesion_barcodes = metadata %>% 
  filter(sample_disease_status == "lesion") %>%
  rownames()
length(lesion_barcodes)

# There's a total of 58887 cell barcodes from non-lesion libraries
non_lesion_barcodes = metadata %>%
  filter(sample_disease_status == "non_lesion") %>%
  rownames()
length(non_lesion_barcodes)

# Extract SCT normalized matrix and subset
sct_counts = rpca_int_sct_v3@assays$SCT@data

# Subset to lesion barcodes and get metadata for each condition
lesion_sct_counts = sct_counts[, lesion_barcodes]
dim(lesion_sct_counts)

lesion_meta = metadata %>%
  filter(sample_disease_status == "lesion")
lesion_meta_df = data.frame(labels=lesion_meta$level1_annotations, 
                            row.names = rownames(lesion_meta))


# Subset to non-lesion barcodes 
non_lesion_sct_counts = sct_counts[, non_lesion_barcodes]
dim(non_lesion_sct_counts)

non_lesion_meta = metadata %>%
  filter(sample_disease_status == "non_lesion")
non_lesion_meta_df = data.frame(labels=non_lesion_meta$level1_annotations,
                                row.names = rownames(non_lesion_meta))

# Create Cellchat objects for lesion and non-lesion libraries
# Looks like we need to run the regular workflow for each cellchat object. 
lesion_level1_labels_cellchat = createCellChat(object = lesion_sct_counts, 
                                               meta = lesion_meta_df,
                                               group.by = "labels")

# Set the used db withiin the object
lesion_level1_labels_cellchat@DB = cellchatdb_hs
   
   
# Pre-process the expression data for downstream analyses
lesion_level1_labels_cellchat = subsetData(lesion_level1_labels_cellchat)
   
lesion_level1_labels_cellchat = identifyOverExpressedGenes(lesion_level1_labels_cellchat)
lesion_level1_labels_cellchat = identifyOverExpressedInteractions(lesion_level1_labels_cellchat)

# Compute communication probabilities
lesion_level1_labels_cellchat = computeCommunProb(lesion_level1_labels_cellchat, population.size = TRUE)
lesion_level1_labels_cellchat = computeCommunProbPathway(lesion_level1_labels_cellchat)
   
# Calculate aggregated cell commnication network
lesion_level1_labels_cellchat = aggregateNet(lesion_level1_labels_cellchat)

# Visualize aggregated cell communication network 
group_size_lesions = as.numeric(table(lesion_level1_labels_cellchat@idents))
netVisual_circle(lesion_level1_labels_cellchat@net$count, vertex.weight = group_size_lesions, weight.scale = TRUE, label.edge = FALSE, title.name = "Lesion number of interactions")
netVisual_circle(lesion_level1_labels_cellchat@net$weight, vertex.weight = group_size_lesions, weight.scale = TRUE, label.edge = FALSE, title.name = "Lesion interactions weight/strength")

# Save cellchat obj 
saveRDS(lesion_level1_labels_cellchat,
        "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cellchat_objects/cellchat_disease_diff_analysis/lesion_level1_cellchat_obj.rds")
lesion_level1_labels_cellchat = readRDS("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cellchat_objects/cellchat_disease_diff_analysis/lesion_level1_cellchat_obj.rds")

# Check signals from ECs to SMCs
endo_smc_level1_lesion = subsetCommunication(lesion_level1_labels_cellchat, 
                                             sources.use = c("Endothelial"),
                                             targets.use = c("SMC"))
endo_smc_level1_lesion

mac_smc_level1_lesion = subsetCommunication(lesion_level1_labels_cellchat, 
                                            sources.use = c("Macrophage"),
                                            targets.use = c("SMC"))
mac_smc_level1_lesion


##############################################################
# Calculate communication probabilities for non-lesion samples
non_lesion_level1_labels_cellchat = createCellChat(object = non_lesion_sct_counts,
                                                   meta = non_lesion_meta_df,
                                                   group.by = "labels")

#Set the used db withiin the object
non_lesion_level1_labels_cellchat@DB = cellchatdb_hs

# Pre-process the expression data for downstream analyses
non_lesion_level1_labels_cellchat = subsetData(non_lesion_level1_labels_cellchat)
    
non_lesion_level1_labels_cellchat = identifyOverExpressedGenes(non_lesion_level1_labels_cellchat)
non_lesion_level1_labels_cellchat = identifyOverExpressedInteractions(non_lesion_level1_labels_cellchat)
  
# Compute communication probabilities
non_lesion_level1_labels_cellchat = computeCommunProb(non_lesion_level1_labels_cellchat, population.size = TRUE)
non_lesion_level1_labels_cellchat = computeCommunProbPathway(non_lesion_level1_labels_cellchat)
  
# Calculate aggregated cell commnication network
non_lesion_level1_labels_cellchat = aggregateNet(non_lesion_level1_labels_cellchat)

# Visualize aggregated cell communication network 
group_size_non_lesions = as.numeric(table(non_lesion_level1_labels_cellchat@idents))
netVisual_circle(non_lesion_level1_labels_cellchat@net$count, vertex.weight = group_size_non_lesions, weight.scale = TRUE, label.edge = FALSE, title.name = "Non-lesion number of interactions")
netVisual_circle(non_lesion_level1_labels_cellchat@net$weight, vertex.weight = group_size_non_lesions, weight.scale = TRUE, label.edge = FALSE, title.name = "Non-lesion interactions weight/strength")

# Save cellchat obj
saveRDS(non_lesion_level1_labels_cellchat,
        "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cellchat_objects/cellchat_disease_diff_analysis/non_lesion_level1_cellchat_obj.rds")
non_lesion_level1_labels_cellchat = readRDS("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cellchat_objects/cellchat_disease_diff_analysis/non_lesion_level1_cellchat_obj.rds")

# Calculate centrality scores
non_lesion_level1_labels_cellchat = netAnalysis_computeCentrality(non_lesion_level1_labels_cellchat)

# Check for signaling coming from ECs to SMCs
endo_smc_level1_non_lesion = subsetCommunication(non_lesion_level1_labels_cellchat, 
                                                 sources.use = c("Endothelial"),
                                                 targets.use = c("SMC"))
endo_smc_level1_non_lesion

fibro_smc_level1_non_lesion = subsetCommunication(non_lesion_level1_labels_cellchat, 
                                                 sources.use = c("Fibroblast"),
                                                 targets.use = c("SMC"))
fibro_smc_level1_non_lesion

mac_smc_level1_non_lesion = subsetCommunication(non_lesion_level1_labels_cellchat, 
                                                  sources.use = c("Macrophage"),
                                                  targets.use = c("SMC"))
mac_smc_level1_non_lesion


########################################################
# Merge level 1 lesion and non-lesion cellchat objects.
level1_obj_list = list(lev1_non_lesion=non_lesion_level1_labels_cellchat, lev1_lesion=lesion_level1_labels_cellchat)
lev1_merged_obj = mergeCellChat(level1_obj_list, add.names = names(level1_obj_list))

# Compare the total number of interactions and interaction strength
compareInteractions(lev1_merged_obj, show.legend = FALSE, group = c(1, 2)) + 
  custom_theme + 
  npg_scale2_bars +
  theme(aspect.ratio = 1.5,
        legend.position = "none")

# Plot differential number of interactions or interaction strength among different cell populations
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(lev1_merged_obj, weight.scale = TRUE, top = 0.5)
netVisual_diffInteraction(lev1_merged_obj, weight.scale = TRUE, measure = "weight", top = 0.3, comparison = c(1,2))

netVisual_heatmap(lev1_merged_obj, targets.use = c("SMC", "Macrophage"))
netVisual_heatmap(lev1_merged_obj, measure = "weight", targets.use = c("SMC", "Macrophage"))

# Control for node size and edge weights of inferred networks across different datasets
weight.max = getMaxWeight(level1_obj_list, attribute = c("idents","count"))
par(mfrow = c(1,1), xpd=TRUE)
plot_list = list()
for (i in 1:length(level1_obj_list)) {
  plot_list[[i]] = netVisual_circle(level1_obj_list[[i]]@net$weight, weight.scale = T, label.edge= F,  top = 0.3,
                   edge.width.max = 12, title.name = paste0("Interaction strength - ", names(level1_obj_list)[i]))
}

# Compare Mac-SMC pathways enriched in lesion vs non-lesion libraries
smc_mac_diff_pathways = rankNet(lev1_merged_obj, mode = "comparison", stacked = TRUE, do.stat = TRUE, sources.use = c("Macrophage"), 
        targets.use = c("SMC"), return.data = TRUE) + 
  custom_theme + 
  scale_fill_lancet() + 
  ggtitle("Mac-Endothelial signaling pathway enrichment") + 
  theme(legend.position = "bottom")

write.table(mac_smc_diff_pathways, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_tables/mac_smc_diff_pathways.tsv",
            row.names = FALSE, col.names = TRUE, sep = "\t")

ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure2/SuppFig2d_lesion_non_lesion_SMC_Mac_diff_pathways.pdf",
       plot = smc_mac_diff_pathways, width = 7, height = 7)


####################################################################
# Run cell communication analysis for macrophage subtypes and SMCs #
####################################################################

# Load metadata and get cell barcodes from lesion and non-lesion libraries to subset counts matrix
metadata = rpca_int_sct_v3@meta.data

# Extract barcodes from lesion samples within the Macrophage compartment. There's a total of 24321 barcodes.  
cell_type_vec = c("Macrophage", "SMC")
mac_smc_lesion_barcodes = metadata %>% 
  filter(sample_disease_status == "lesion" & level1_annotations %in% cell_type_vec) %>%
  rownames()
length(mac_smc_lesion_barcodes)

# Extract SCT normalized matrix and subset
sct_counts = rpca_int_sct_v3@assays$SCT@data

# Subset to mac lesion barcodes and get metadata for each condition
mac_smc_lesion_sct_counts = sct_counts[, mac_smc_lesion_barcodes]
dim(mac_smc_lesion_sct_counts)

# Create metadata df for SMC and Mac lesion barcodes
mac_smc_lesion_meta = metadata %>%
  filter(sample_disease_status == "lesion" & level1_annotations %in% cell_type_vec)
#mac_smc_lesion_meta_df = data.frame(labels=mac_smc_lesion_meta$level2_annotations,
#                                row.names = rownames(mac_smc_lesion_meta))
#table(mac_smc_lesion_meta_df$labels)

# Create Cellchat objects for lesion and non-lesion libraries
# Looks like we need to run the regular workflow for each cellchat object. 
#mac_smc_lesion_level2_cellchat = createCellChat(object = mac_smc_lesion_sct_counts, 
#                                               meta = mac_smc_lesion_meta_df,
#                                               group.by = "labels")

#Set the used db withiin the object
#mac_smc_lesion_level2_cellchat@DB = cellchatdb_hs

# Pre-process the expression data for downstream analyses
#mac_smc_lesion_level2_cellchat = subsetData(mac_smc_lesion_level2_cellchat)

#mac_smc_lesion_level2_cellchat = identifyOverExpressedGenes(mac_smc_lesion_level2_cellchat)
#mac_smc_lesion_level2_cellchat = identifyOverExpressedInteractions(mac_smc_lesion_level2_cellchat)

# Compute communication probabilities
#mac_smc_lesion_level2_cellchat = computeCommunProb(mac_smc_lesion_level2_cellchat, population.size = TRUE)
#mac_smc_lesion_level2_cellchat = computeCommunProbPathway(mac_smc_lesion_level2_cellchat)

# Calculate aggregated cell commnication network
mac_smc_lesion_level2_cellchat = aggregateNet(mac_smc_lesion_level2_cellchat)

# Save cellchat obj
#saveRDS(mac_smc_lesion_level2_cellchat,
#        "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cellchat_objects/cellchat_disease_diff_analysis/mac_smc_lesion_level2_cellchat_obj.rds")
mac_smc_lesion_level2_cellchat = readRDS("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cellchat_objects/cellchat_disease_diff_analysis/mac_smc_lesion_level2_cellchat_obj.rds")

# # Visualize aggregated cell communication network 
group_size_mac_smc_lesion_level2 = as.numeric(table(mac_smc_lesion_level2_cellchat@idents))
netVisual_circle(mac_smc_lesion_level2_cellchat@net$count, vertex.weight = group_size_mac_smc_lesion_level2, weight.scale = TRUE, label.edge = FALSE, top = 0.3, title.name = "Mac-SMC lesion number of interactions")
netVisual_circle(mac_smc_lesion_level2_cellchat@net$weight, vertex.weight = group_size_mac_smc_lesion_level2, weight.scale = TRUE, label.edge = FALSE, top=0.3, title.name = "Mac-SMC lesion interactions weight/strength")
 
# Inspect pathways enriched in lesions
smc_incoming_pathways = subsetCommunication(mac_smc_lesion_level2_cellchat, targets.use = c("SMC"))
netVisual_aggregate(mac_smc_lesion_level2_cellchat, signaling = c("SPP1"), layout = "circle" , targets.use = c("SMC"))
netVisual_heatmap(mac_smc_lesion_level2_cellchat, signaling = c("SPP1"), color.heatmap = "Reds", targets.use = c("SMC"))
 
tweak_signaling = netVisual_aggregate(mac_smc_lesion_level2_cellchat, signaling = c("TWEAK"), layout = "circle" , targets.use = c("SMC"))
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure2/SuppFig2e_TWEAK_signaling.pdf",
        plot=tweak_signaling, width = 6, height = 6)

resistin = netVisual_aggregate(mac_smc_lesion_level2_cellchat, signaling = c("RESISTIN"), layout = "circle" , targets.use = c("SMC"))
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure2/SuppFig2e_RESISTIN_signaling.pdf",
       plot=tweak_signaling, width = 6, height = 6)

spp1_LR_contributions = netAnalysis_contribution(mac_smc_lesion_level2_cellchat, signaling = "SPP1", return.data = TRUE)
spp1_lr_plot = spp1_LR_contributions$LR.contribution %>%
   ggplot(aes(x=reorder(name, contribution), y=contribution)) +
   geom_col(width=0.7, fill="#053061") + 
   xlab("") + 
   ylab("Relative contribution") + 
   coord_flip() + 
   custom_theme + 
   theme(aspect.ratio = 0.4,
         axis.ticks.y = element_blank())
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure3/Fig3d_SPP1_signaling_LR_contributions.pdf",
       plot=spp1_lr_plot, width = 7, height = 4)


#################################################################
# Run cell communication analyses for myeloid and SMC subtypes  #
#################################################################

# Load entire reference
rpca_int_sct_v3 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/seurat_objects/alsaigh_pan_wirka_hu_int_seurat_annotated_v3.rds")
DefaultAssay(rpca_int_sct_v3) = "integrated" 
 
# Set clusters to latest version of clustering resolutio; res=1
Idents(rpca_int_sct_v3) = "integrated_snn_res.1"
table(rpca_int_sct_v3@meta.data$level1_annotations)

# Subset reference to myeloid and SMCs
mac_smc_rpca_int_sct_v3 = subset(rpca_int_sct_v3, subset = level1_annotations %in% c("SMC", "Macrophage"))
table(mac_smc_rpca_int_sct_v3@meta.data$level1_annotations)
dim(mac_smc_rpca_int_sct_v3@meta.data)

mac_smc_metadata = mac_smc_rpca_int_sct_v3@meta.data
mac_smc_barcodes = rownames(mac_smc_metadata)
names(mac_smc_barcodes) = mac_smc_metadata$level2_annotations


# There's an overlap of 27228 SMCs between the main reference and the smc peri fibro subset. 
# These are the cells that don't overlap (1796) 

vec = mac_smc_barcodes[names(mac_smc_barcodes) == "SMC"]
not_ol = setdiff(vec, smc_barcodes)


# Remove those from the main subset
mac_smc_rpca_int_sct_v3_2 = subset(mac_smc_rpca_int_sct_v3, cells = not_ol, invert=TRUE)
table(mac_smc_rpca_int_sct_v3_2@meta.data$level1_annotations)

# Get vector with Mac and SMC barcodes
mac_smc_barcodes_2 = rownames(mac_smc_rpca_int_sct_v3_2@meta.data)
names(mac_smc_barcodes_2) = mac_smc_rpca_int_sct_v3_2@meta.data$level2_annotations 


# Create a new vector where we will update the SMC annotations
mac_smc_barcodes_new = mac_smc_barcodes_2
for (i in seq_len(length(mac_smc_barcodes_new))) { 
  if (names(mac_smc_barcodes_new[i]) == "SMC") { 
      idx = match(mac_smc_barcodes_new[i], smc_barcodes)
      names(mac_smc_barcodes_new)[i] = names(smc_barcodes[idx])
      }
  }

# Add additional col into metadata
mac_smc_rpca_int_sct_v3_2@meta.data$level2_annotations_SMC_updated = names(mac_smc_barcodes_new)

###############################################################################
# Load SMC, peri and fibro batch corrected reference with cell type annotations
rpca_smc_fibro_subset_v3 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/seurat_objects/whole_ref_v3/peri_smc_fibro_subset_seurat_obj_res0.9_v2.rds")
DefaultAssay(rpca_smc_fibro_subset_v3) = "integrated"
table(rpca_smc_fibro_subset_v3@meta.data$prelim_annotations)

# Get SMC barcodes with their annotations
smc_barcodes = c("Contractile_SMC", "Transitional-ECM-SMC", "Fibromyocyte", "Fibrochondrocyte", "SMC2", "SMC3", "Foam-like")
smc_metadata = rpca_smc_fibro_subset_v3@meta.data %>%
  filter(prelim_annotations %in% smc_barcodes) %>%
  dplyr::select(prelim_annotations)

smc_barcodes = rownames(smc_metadata)
names(smc_barcodes) = smc_metadata$prelim_annotations 

##################################################################################
# Get matrix and metada with annotations for Myeloid and SMC subtypes from lesion

mac_smc_lesion_meta = mac_smc_rpca_int_sct_v3_2@meta.data %>%
  filter(sample_disease_status == "lesion") 


# Get sct counts
mac_smc_sct_counts = mac_smc_rpca_int_sct_v3_2@assays$SCT@data

# Filter sct counts matrix to contain cells only from lesions
mac_smc_lesion_sct_counts = mac_smc_sct_counts[, rownames(mac_smc_lesion_meta)]

# Create metadata df for cellchat
mac_smc_lesion_meta_df = data.frame(labels=mac_smc_lesion_meta$level2_annotations_SMC_updated,
                                    row.names = rownames(mac_smc_lesion_meta))
table(mac_smc_lesion_meta_df$labels)
                                      
# Create Cellchat objects for lesion and non-lesion libraries
# Looks like we need to run the regular workflow for each cellchat object. 
mac_smc_lesion_level2_cellchat = createCellChat(object = mac_smc_lesion_sct_counts, 
                                                meta = mac_smc_lesion_meta_df,
                                                group.by = "labels")
                                      
#Set the used db withiin the object
mac_smc_lesion_level2_cellchat@DB = cellchatdb_hs
              
# Pre-process the expression data for downstream analyses
mac_smc_lesion_level2_cellchat = subsetData(mac_smc_lesion_level2_cellchat)
                                      
mac_smc_lesion_level2_cellchat = identifyOverExpressedGenes(mac_smc_lesion_level2_cellchat)
mac_smc_lesion_level2_cellchat = identifyOverExpressedInteractions(mac_smc_lesion_level2_cellchat)

# Compute communication probabilities
mac_smc_lesion_level2_cellchat = computeCommunProb(mac_smc_lesion_level2_cellchat, population.size = TRUE)
mac_smc_lesion_level2_cellchat = computeCommunProbPathway(mac_smc_lesion_level2_cellchat)
 
# Calculate aggregated cell commnication network
mac_smc_lesion_level2_cellchat = aggregateNet(mac_smc_lesion_level2_cellchat)

# Save cellchat object
saveRDS(mac_smc_lesion_level2_cellchat,
        "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cellchat_objects/cellchat_disease_diff_analysis/mac_smc_lesion_subtypes_level2_cellchat_obj.rds")
mac_smc_lesion_level2_cellchat = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cellchat_objects/cellchat_disease_diff_analysis/mac_smc_lesion_subtypes_level2_cellchat_obj.rds")

# Visualize aggregated cell communication network 
group_size_mac_smc_lesion_level2 = as.numeric(table(mac_smc_lesion_level2_cellchat@idents))
int_number = netVisual_circle(mac_smc_lesion_level2_cellchat@net$count, vertex.weight = group_size_mac_smc_lesion_level2, weight.scale = TRUE, label.edge = FALSE, top = 0.15, title.name = "Mac-SMC lesion number of interactions")
int_weight = netVisual_circle(mac_smc_lesion_level2_cellchat@net$weight, vertex.weight = group_size_mac_smc_lesion_level2, weight.scale = TRUE, label.edge = FALSE, top=0.15, title.name = "Mac-SMC lesion interactions weight/strength")
 


# Inspect pathways enriched in lesions
smc_incoming_pathways = subsetCommunication(mac_smc_lesion_level2_cellchat, 
                                            sources.use = c("Foamy_Mac1", "Foamy_Mac2", "Inflammatory_Mac", "Monocytes", "Monocytes/DC",
                                                            "NAMPT_Neutrophils", "Phagocytosis_Mac"),
                                            targets.use = c("Contractile_SMC", "Transitional-ECM-SMC", "Fibromyocyte", "Fibrochondrocyte", "Foam-like"))
write.table(smc_incoming_pathways,
            "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_tables/Myeloid_SMC_level2_LR_interactions.tsv",
            row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

# Make a summary figure
smc_incoming_pathways_plot = smc_incoming_pathways %>%
  mutate(cell_type_pair = paste(source,"-",target, sep = "")) %>%
  ggplot(aes(x=cell_type_pair, y=interaction_name, color=-log2(prob))) + 
  geom_point(size=4) +
  ylab("Interaction") + 
  custom_theme +
  new_scale3 +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 10) 
        #axis.text = element_text(size=8),
        #axis.text.y = element_text(size = 9))
  )

ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/remaining_cell_communication_figures/SMC_incoming_pathways_summary_new.pdf",
       plot = smc_incoming_pathways_plot, width = 14, height = 8)

netVisual_aggregate(mac_smc_lesion_level2_cellchat, signaling = c("RESISTIN"), layout = "circle" , 
                    targets.use = c("Contractile_SMC", "Transitional-ECM-SMC", "Fibromyocyte", "Fibrochondrocyte", "Foam-like"), top = 0.15)
netVisual_heatmap(mac_smc_lesion_level2_cellchat, signaling = c("SPP1"), color.heatmap = "Reds", targets.use = c("Contractile_SMC", "Transitional-ECM-SMC", "Fibromyocyte", "Fibrochondrocyte", "Foam-like"))
 
tweak_signaling = netVisual_aggregate(mac_smc_lesion_level2_cellchat, signaling = c("SPP1"), layout = "circle" , targets.use = c("SMC"))
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure2/SuppFig2e_TWEAK_signaling.pdf",
       plot=tweak_signaling, width = 6, height = 6)

spp1_LR_contributions = netAnalysis_contribution(mac_smc_lesion_level2_cellchat, signaling = "SPP1", return.data = TRUE)
spp1_lr_plot = spp1_LR_contributions$LR.contribution %>%
  ggplot(aes(x=reorder(name, contribution), y=contribution)) +
  geom_col(width=0.7, fill="#053061") + 
  xlab("") + 
  ylab("Relative contribution") + 
  coord_flip() + 
  custom_theme + 
  theme(aspect.ratio = 0.4,
  axis.ticks.y = element_blank())

# Save plots
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure3/Fig3d_SPP1_signaling_LR_contributions.pdf",
       plot=spp1_lr_plot, width = 7, height = 4)





