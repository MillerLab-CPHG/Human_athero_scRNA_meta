library(Seurat)
library(scRNAutils)
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(parallel)
library(monocle3)
library(pheatmap)
library(ggrepel)
library(ggsci)
library(ggridges)

# Set seed for reproducibility
set.seed(1)

# Source our own scRNA analysis utils functions
#source("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/R_scripts/Utils/scRNA_processing_utils.R")


#############################################################################################################
# This script contains code for pseudotime analysis of the subclustered SMCs                                #
# The general strategy will be to create the monocle3 cds object using SCT normalized counts                #
# and then importing PCA, UMAP embeddings and clusters information from the seurat integrated object into   #
# the new monocle object                                                                                    #
#############################################################################################################


# Read subclustered seurat object
rpca_smc_fibro_subset_v3 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/seurat_objects/whole_ref_v3/peri_smc_fibro_subset_seurat_obj_res0.9_v1.rds")
DefaultAssay(rpca_smc_fibro_subset_v3) = "SCT"

# Further subset to keep only SMCs
rpca_smc_subset_v3 = subset(rpca_smc_fibro_subset_v3, idents = c(5, 2, 16, 4, 0 ,9, 6, 13, 17))

# Extract SCT counts sparse matrix
sct_counts = rpca_smc_subset_v3@assays$SCT@counts
class(sct_counts)

# Extract cell metadata as a data.frame object 
cell_metadata = as.data.frame(rpca_smc_subset_v3@meta.data)
class(cell_metadata)

# Extract feature metadata for each gene, should simply be the name of the gene
gene_metadata = data.frame(gene_short_name = row.names(sct_counts),
                           row.names = row.names(sct_counts))

# Create cell data set object for monocle3. 
# We shouldn't need to normalize expression data since we're importing data from the SCT assay
cds = monocle3::new_cell_data_set(expression_data = sct_counts,
                                  cell_metadata = cell_metadata,
                                  gene_metadata = gene_metadata)
class(cds)

# Create toy dim reduction object without normalizing data
cds_ordered = monocle3::preprocess_cds(cds = cds_ordered)


# Now that gene expression was entered, we should insert dim reduced embeddings information from seurat object
# VERY IMPORTANT: if data is not to be normalized, we need to manually input cds@preprocess_aux$prop_var_expl and
# cds@preprocess_aux$gene_loadings from the corresponding Seurat object. 
reducedDim(cds, type = "PCA") = rpca_smc_subset_v3@reductions$pca@cell.embeddings
cds@preprocess_aux$prop_var_expl = rpca_smc_subset_v3@reductions$pca@stdev
cds@preprocess_aux$gene_loadings = rpca_smc_fibro_subset_v3@reductions$pca@feature.loadings
plot_pc_variance_explained(cds)

# Create toy UMAP embeddings
#cds = monocle3::reduce_dimension(cds = cds, preprocess_method = "PCA")

# Transfer UMAP embeddings from seurat object
cds@int_colData@listData$reducedDims$UMAP = rpca_smc_subset_v3@reductions$umap@cell.embeddings

# Plot cells using seurat extracted PCA and UMAP embeddings
monocle3::plot_cells(cds = cds, color_cells_by = "prelim_annotations", cell_size = 0.7, 
                     show_trajectory_graph = FALSE) + custom_theme() +
  theme(legend.position = "bottom")
monocle3::plot_cells(cds = cds, genes = c("CNN1", "FN1", "CLU", "TNFRSF11B"), cell_size = 0.5, 
                     show_trajectory_graph = FALSE, 
                     norm_method = "log") + custom_theme() & new_scale

# We already have clusters information from the subset Seurat object
# Extract cluster IDs from seurat object (metadata from seurat object is stored in cds@colData)
seurat_clusters = rpca_smc_subset_v3@meta.data$seurat_clusters
head(seurat_clusters)
names(seurat_clusters) = rownames(rpca_smc_subset_v3@meta.data)
cds@clusters$UMAP_so$clusters = seurat_clusters


# Run cluster_cells so that we can get the partitions required for building the pseudotime trajectory.
# Monocle3 clusters are located at cds@clusters$UMAP. Might be handy to store both monocle and Seurat cluster IDs. 
cds = monocle3::cluster_cells(cds, reduction_method = "UMAP", 
                              cluster_method = "leiden", 
                              resolution = 1e-4)

# Replace monocle clusters by Seurat clusters
cds@clusters$UMAP$clusters = seurat_clusters

# Fix from https://gitmemory.com/cole-trapnell-lab
rownames(cds@principal_graph_aux$UMAP$dp_mst) = NULL
colnames(cds@int_colData@listData$reducedDims$UMAP) = NULL

# Visualize cells according to Monocle clusters/partitions and Seurat clusters.
# There should be only 1 partition since one trajectory is drawn for each partition. 
# When you are learning trajectories, each partition will eventually become a separate trajectory.
plot_cells(cds, color_cells_by = "partition", 
           group_cells_by = "partition", 
           show_trajectory_graph = FALSE, 
           label_groups_by_cluster = TRUE, 
           labels_per_group = 1, 
           cell_size = 0.7) + custom_theme()
plot_cells(cds, color_cells_by = "cluster", 
           show_trajectory_graph = FALSE, 
           cell_size = 0.7) + custom_theme
plot_cells(cds, color_cells_by = "seurat_clusters", 
           show_trajectory_graph = FALSE, 
           label_groups_by_cluster = TRUE, labels_per_group = 1, 
           label_cell_groups = FALSE, 
           cell_size = 0.7) + custom_theme()

##############################
# Learn the trajectory graph #
##############################

# Keeping pericytes and fibroblasts might be tricky since MOnocle3 forces a trajectory across all cell types
# even if not biologically meaningful. Might have to further subset cells to keep only SMC-derived cells. 
cds = monocle3::learn_graph(cds = cds)
plot_cells(cds = cds_ordered, color_cells_by = "cluster", label_branch_points = FALSE, 
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_principal_points = TRUE,
           label_groups_by_cluster = FALSE,
           cell_size = 0.7) + custom_theme()

# Select root cells and order according to pseudotime
# We're defining cells with the higghest MYH11 expression as the root of
# the trajectory
cds_ordered = order_cells(cds_ordered, root_pr_nodes = c("Y_86"))

# Plot cells according to calculated pseudotime 
plot_cells(cds = cds_ordered, color_cells_by = "pseudotime", 
           label_branch_points = FALSE, 
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_groups_by_cluster = FALSE,
           cell_size = 0.7, graph_label_size = 5,
           alpha = 1) + miller_continous_scale()
  custom_theme() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12))

# Plot cells according to sample disease status
plot_cells(cds = cds_ordered, 
           color_cells_by = "sample_disease_status", 
           label_branch_points = FALSE, 
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_groups_by_cluster = FALSE,
           cell_size = 0.5, graph_label_size = 5,
           alpha = 1) +  
  custom_theme() + 
  miller_discrete_scale() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "bottom",
        legend.text = element_text(size = 14)) 

# Plot cells according to arterial origin
arterial_origin_pseudotime = plot_cells(cds = cds_ordered, 
           color_cells_by = "arterial_origin", 
           label_branch_points = FALSE, 
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_groups_by_cluster = FALSE,
           cell_size = 0.5, graph_label_size = 5,
           alpha = 1) +  
  custom_theme() + 
  miller_discrete_scale(option = 2) + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "bottom",
        legend.text = element_text(size = 14))

ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Revisions_SuppFig7/SuppFig7a_arterial_origin_pseudotime.pdf",
       arterial_origin_pseudotime, width = 7, height = 7)



#########################################################################
# Plot distribution of cells from lesion status along pseudotime
# Create df of subset SMCs for gene expression across pseudotime plotting
metadata = colData(cds_ordered)

# Add pseudotime values for each cell
pseudotime = pseudotime(cds_ordered)
metadata$pseudotime = pseudotime[match(rownames(metadata), 
                                        names(pseudotime))]
metadata_df = as.data.frame(metadata)
lesion_density = metadata_df %>%
  ggplot(aes(x=pseudotime, y=sample_disease_status, fill=sample_disease_status
             )) + 
  geom_density_ridges() + 
  xlim(0, 34) + 
  ylab("Cell density") + 
  custom_theme() +
  theme(aspect.ratio = 0.6,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom") + 
  miller_discrete_scale(style="bars") 

ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Revisions_SuppFig7/SuppFig7a_lesion_pseudotime_density.pdf",
       plot = arterial_pseudotime_density, width = 7, height = 4)

#####################
# Save monocle object 
saveRDS(cds_ordered, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/smcs_ordered_cds_monocle3_obj.rds")
cds_ordered = readRDS("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/monocle_objects/smcs_ordered_cds_monocle3_v3_obj.rds")

#####################################################
# Distribution of gene expression across pseudotime #
#####################################################

# In this section we'll find genes and modules that change as a function of pseudotime
# Before running this function, need to run trace('calculateLW', edit = T, where = asNamespace("monocle3"))
# and change Matrix::rBind to rbind
smc_cds_pr_test_res = monocle3::graph_test(cds_ordered, neighbor_graph = "principal_graph", cores = 4)
saveRDS(smc_cds_pr_test_res, "~/Desktop/human_athero_scRNA_meta-analysis/rds_objects/monocle_objects/smc_monocle3_cds_v3_graph_test.rds")
smc_cds_pr_test_res = readRDS("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/monocle_objects/smc_monocle3_cds_v3_graph_test.rds")

# Filter df to keep only significant genes at FDR <= 0.05
significant_genes_df = smc_cds_pr_test_res %>% 
  filter(q_value <= 0.05) %>%
  arrange(q_value)

# Add numeric indexes to genes
dim(significant_genes_df)
significant_genes_df$gene_index = 1:nrow(significant_genes_df)

# Create a plot showing significant genes
significant_genes_df %>%
  head(n=1000) %>%
  ggplot(aes(x=gene_index ,y=-log10(q_value),
             label= ifelse(-log10(q_value) > 100, 
                           gene_short_name, ""))) + 
  geom_point() + 
  geom_text_repel() + 
  theme_bw()

##################################################################
# Need to add pseudotime values into seurat object so can create a 
# matrix for the heatmap

# set default assay to SCT, just in case
DefaultAssay(rpca_smc_subset_v3) = "SCT"

# Add pseudotime values into Seurat object metadata
rpca_smc_subset_v3@meta.data$pseudotime = pseudotime(cds_ordered)

# Round pseudotime values
rpca_smc_subset_v3$pseudotime_value_round = round(rpca_smc_subset_v3$pseudotime,
                                                  digits = 1)
rpca_smc_subset_v3$pseudotime_value_round = factor(rpca_smc_subset_v3$pseudotime_value_round,
                                                   levels = unique(rpca_smc_subset_v3$pseudotime_value_round)[order(unique(rpca_smc_subset_v3$pseudotime_value_round))])

# Define pseudotime points as idents to get avga gene expression
Idents(rpca_smc_subset_v3) = "pseudotime_value_round"

# Get Avg expression of genes at each point in pseudotime
avgexp = AverageExpression(rpca_smc_subset_v3, 
                           assays = c("SCT"), 
                           slot = "data",
                           return.seurat = TRUE)

# Define top genes (check significant genes based on q value)
# Get the top 500 significant genes based on q values. 
top_de_genes = head(rownames(significant_genes_df), n=500)

# Create matrix of average gene expression by pseudotime values
expression.values = FetchData(avgexp, vars = top_de_genes, 
                              slot = "scale.data")

# Use cluster.pseudotime.genes function from Peisker et al., 2022
cluster_pseudotime_genes = function(expression.values, k.treecut=5, keep.hclust=F){
  #expression.values: expression values from FetchData
  #k.treecut: cutoff value for the clustertree, based on this there will be more or less groups of genes
  expression.values = as.data.frame(t(expression.values))
  d = dist(expression.values, method = "euclidean")
  clust = hclust(d, method = "complete")
  plot(clust, labels = FALSE)
  clust = cutree(clust, k = k.treecut) %>%  data.frame()
  names(clust)="hcluster"
  expression.values = cbind(expression.values,clust)
  expression.values = arrange(expression.values,hcluster)
  
  order.df=NULL
  for (k in levels(factor(expression.values$hcluster))) {
    k_sub = expression.values[expression.values$hcluster==k,]
    k_sub$hcluster=NULL
    k_sub.means = colMeans(k_sub)
    df = data.frame("ps"=names(k_sub.means[k_sub.means == max(k_sub.means)]),"k"=k)
    order.df=rbind(order.df,df)
  }
  order.df = arrange(order.df,ps)
  expression.values$hcluster = factor(expression.values$hcluster,levels=order.df$k )
  expression.values = arrange(expression.values,hcluster)
  if (keep.hclust) {
    return(expression.values) 
  }else{
    expression.values$hcluster = sapply(expression.values$hcluster,function(k) strrep("_",k) )
    rownames(expression.values) = paste0(expression.values$hcluster,rownames(expression.values))
    expression.values$hcluster=NULL
    return(expression.values)}
}

# Cluster genes in pseudotime
expression.values = cluster_pseudotime_genes(expression.values)

# Highlight both DE genes across pseudotime and also hits from the MAGMA gene effector analysis
genes_to_highlight = c("_LMOD1", "_MYH11", "__ACTA2", "_FHL5", "_CARMN", "_ACTA2", "_TNS1", "_MYOCD",  
                       "__CRIP1", "__MFAP4", "__PFN1", "__BTF3", 
                       "___AEBP1", "___MMP2", "_PALLD", "MRV1", "_RBPMS", "_PALLD", "___PDGFD", "___GEM",
                       "____LUM", "____COMP", "____TIMP1", "___COL8A1", "____CRTAC1", 
                       "___FN1", "____COL6A2", "___LTBP1", 
                       "___PDGFRB", "____COL1A2", "___TGFB1", 
                       "___COL4A1", "___CDH13", "____AGT", "____LOXL1", "____SPRY1",
                       "____TIMP1", "____CYTL1", "____PCOLCE2")

newnames=NULL
for (row.nr in 1:nrow(expression.values)){
  if(row.names(expression.values[row.nr,]) %in% genes_to_highlight) 
    {newnames=c(newnames,bquote(italic(.(row.names(expression.values[row.nr,])))))
  } else {
    newnames=c(newnames,"")
  }
}

# Plot heatmap of gene expression across pseudotime
expression.heatmap = pheatmap(expression.values,
                              labels_row = as.expression(newnames), labels_col = "",
                              cluster_cols = FALSE, cluster_rows = FALSE,
                              angle_col = 0,border_color = 0,
                              treeheight_row = 20,
                              fontsize = 9, 
                              scale = "column",
                              breaks = seq(-4, 4, length.out = 100),
                              color = colorRampPalette(c("#053061", "#2171B5", "white", "#FFD92F", "#A50F15"))(100)
                              ) 

ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure6/revision_figures/Fig6b_pseudotime_expression_heatmap.pdf",
       plot = expression.heatmap, width = 6, height = 8)

# Find which CAD effector genes are expressed across pseudotime
# Overlapping genes: CLEC18A, CDH13, FHL5, MRVI1, LMOD1, AGT, MYH11, COL4A2, COL4A1,
# RBPMS, RBPMS2, TNS1, PALLD, PDGFD, GEM
idx = unlist(lapply(cad_smc_genes_vec, grep, x=rownames(expression.values)))
cad_pseudotime_genes = rownames(expression.values)[idx]

#############################
# Group DE genes into modules

# Keep only genes that pass the FDR < 0.05 threshold
smc_deg_ids = row.names(subset(smc_cds_pr_test_res, q_value < 0.05))

# Collect trajectory-variable genes into modules
gene_modules_df = find_gene_modules(cds_ordered[smc_deg_ids,], resolution=c(10^seq(-6,-1)))
saveRDS(gene_modules_df,
        "~/Desktop/human_athero_scRNA_meta-analysis/rds_objects/monocle_objects/smc_monocle3_cds_v3_gene_modules_df.rds")

# Plot module scores within each cell type annotation
cell_group_df = tibble::tibble(cell=row.names(colData(cds_ordered)), 
                                cell_group=colData(cds_ordered)$prelim_annotations)
saveRDS(cell_group_df,
        "~/Desktop/human_athero_scRNA_meta-analysis/rds_objects/monocle_objects/smc_monocle3_cds_v3_cell_group_df.rds")
agg_mat = aggregate_gene_expression(cds_ordered, gene_modules_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))

# Make a heatmap of modules enrichment per cell type 
palette_length = 100
#my_color = colorRampPalette(c("#053061", "white", "darkred"))(palette_length)
my_color3 = colorRampPalette(c("#053061", "#2171B5", "white", "#FFD92F", "#A50F15"))(palette_length)
pseudotime_heatmap = pheatmap::pheatmap(agg_mat, scale="column", fontsize = 14,
                                        color = my_color3,
                                        angle_col = 45,
                                        clustering_method="ward.D2") + custom_theme
ggsave(file="~/Desktop/human_athero_scRNA_meta-analysis/manuscript_figures/Supplemental_Fig4/pseudotime_DE_modules_heatmap.pdf", plot = pseudotime_heatmap, width = 8, height = 8)

# Plot gene modules in UMAP space
plot_cells(cds_ordered,
           genes=gene_modules_df %>% filter(module %in% c(9, 4, 5, 10)),
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 5) & new_scale & custom_theme

# Inspect interesting modules
module1 = gene_modules_df %>% filter(module == 1)
module1$id


#######################
# Fibromyo module
# Interesting genes: COL4A2 (CAC GWAS), LRP1, AEBP1, PDGFRB, LGALS3BP, ITGBL1, VCAN, FN1, AEBP1
# COL6A1, COL6A2, COL6A3, TIMP3, COL4A1, COL15A1, LTBP1, LTBP2, JAG1, TIMP2, OGN, FMOD
module4 = gene_modules_df %>% filter(module == 4)
module4$id

#######################
# Transitional SMC
# Module 5: 
module5 = gene_modules_df %>% filter(module == 5)
module5$id

# Module 10: ARID5B, LGALS3, BGN, KRT8, KRT18, SPARC, RGS5, NOTCH3, TUB1X
module10 = gene_modules_df %>% filter(module == 10)
module10$id

######################
# Fibrochondro modules
# Module 9 interesting genes:  BMP4, MMP2, GSN, FBLN2, PCOLCE2, SERPINE2, WISP2, MFAP5, DPT, CCDC80, LUM, DCN, 
# PCOLCE2, PODN
# SPRY genes are expressed in pre-hypertrophic and hypertrophic chondrocytes. 
module9 = gene_modules_df %>% filter(module == 9)
module9$id

module11 = gene_modules_df %>% filter(module == 11)
module11$id

# Module 13 interesting genes: IGFBP3, LAMA2, C3, NRP1, SFRP2, CXCL12
module13 = gene_modules_df %>% filter(module == 13)
module13$id

#######################################################################################
# Here we'll use the custom script we wrote for plotting gene expression as a function 
# of pseudotime (sourced from our own scRNA utils script)

#cds_ordered_v2 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/monocle_objects/smcs_ordered_cds_monocle3_v2_obj.rds")

# First we need to define the cell types and genes if interest we want to model
cells_types_of_interest = c("Contractile_SMC", "Transitional-ECM-SMC", "Fibromyocyte")
cells_types_of_interest2 = c("Contractile_SMC", "Transitional-ECM-SMC", "Fibrochondrocyte")
smc_fibromyo_genes = c("CNN1", "FN1", "LGALS3", "AEBP1", "LTBP1", "PDGFRB")
contractile_markers = c("LMOD1", "MYH11", "ACTA2", "TAGLN")
smc_fibrochondro_genes = c("CNN1", "COL1A2", "PCOLCE2", "CRTAC1", "COMP", "MMP2")


# Input cell type and genes of interest args to function. 
fibromyo_genes_pseudotime = plot_expression_on_pseudotime(cds_ordered = cds_ordered_v2, 
                                                          genes = smc_fibromyo_genes, 
                              cell_annotations = cells_types_of_interest,
                              facet_wrap_plot = TRUE) + 
  miller_discrete_scale(style = "bars")

ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure6/revision_figures/Fig6c_SMC_fibromyo_markers_pseudotime.pdf",
       plot = fibromyo_genes_pseudotime, width = 8, height = 4)

contractile_genes_pseudotime = plot_expression_on_pseudotime(cds_ordered = cds_ordered_v2, 
                                                          genes = contractile_markers, 
                                                          cell_annotations = cells_types_of_interest,
                                                          facet_wrap_plot = TRUE)

ggsave("~/Desktop/human_athero_scRNA_meta-analysis/manuscript_figures/Supplementary_Figure5/SuppFig5b_contractile_genes_across_pseudotime.pdf",
       plot = contractile_genes_pseudotime, width = 4.5, height = 4.5)

fibrochondro_genes_pseudotime = plot_expression_on_pseudotime(cds_ordered = cds_ordered_v2, 
                                                          genes = smc_fibrochondro_genes, 
                                                          cell_annotations = cells_types_of_interest2,
                                                          facet_wrap_plot = TRUE) + 
  miller_discrete_scale(style = "bars")

ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure6/revision_figures/Fig6d_SMC_fibrochondro_markers_pseudotime.pdf",
       plot = fibrochondro_genes_pseudotime, width = 8, height = 4)












