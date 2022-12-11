library(Seurat)
library(dorothea)
library(decoupleR)
library(tidyverse)
library(data.table)
library(cluster)
library(RColorBrewer)


# Set seed for reproducibility
set.seed(1)

# Source our own scRNA analysis utils functions
source("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/R_scripts/Utils/scRNA_processing_utils.R")


###################################################################################################
# This script will be used to infer TF activity in SMC clusters using DoRothEA regulons and VIPER #
###################################################################################################


# Load batch corrected reference with cell type annotations
rpca_smc_fibro_subset_v3 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/seurat_objects/whole_ref_v3/peri_smc_fibro_subset_seurat_obj_res0.9_v2.rds")
DefaultAssay(rpca_smc_fibro_subset_v3) = "integrated"

# Subset seurat object to keep only SMCs (Contractile, transitional, fibromyocytes and fibrochondrocytes; 20127 cells)
rpca_smc_subset_v3 = subset(rpca_smc_fibro_subset_v3, idents = c(5, 2, 16, 4, 0 ,9, 6, 13, 17))
DimPlot(rpca_smc_subset_v3, group.by = "prelim_annotations", pt.size = 0.4) + npg_scale + custom_theme
rpca_smc_subset_v3

# Read Dorothea regulons from human
dorothea_regs_human = get(data("dorothea_hs", package = "dorothea"))
dorothea_regs_human = as.data.frame(dorothea_regs_human)
dim(dorothea_regs_human)
head(dorothea_regs_human)
dorothea_regs_human[dorothea_regs_human$tf == "SOX9", ]
dorothea_regs_human[dorothea_regs_human$target == "LGALS3", ]

# Let's get regulons based on interactions with high confidence levels  
high_conf_regulons = dorothea_regs_human %>%
  dplyr::filter(confidence %in% c("A", "B", "C") | tf %in% c("TCF21", "MEF2D"))

# There are only 13223 TF-target interactions from high confidence levels
dim(high_conf_regulons)

# Compute Viper scores using SCT normalized data (include A-C scores to see what results look like; this might take a while)
DefaultAssay(rpca_smc_subset_v3) = "SCT"
rpca_smc_subset_v3 = dorothea::run_viper(rpca_smc_subset_v3, 
                                         regulons=high_conf_regulons,
                                         assay_key = "SCT",
                                         options = list(method="scale", 
                                                        minsize=4,
                                                        eset.filter=FALSE, 
                                                        cores=5))

# Scale tf activity scores
DefaultAssay(rpca_smc_subset_v3) = "dorothea"
rpca_smc_subset_v3 = ScaleData(rpca_smc_subset_v3)


# Load seurat object with viper scores
rpca_smc_subset_v3 = readRDS("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/seurat_objects/whole_ref_v3/rpca_smc_subset_res0.9_v1_dorothea.rds")

# Transform scaled viper scores into a df to better handle results
viper_scores_df = GetAssayData(rpca_smc_subset_v3, slot = "scale.data", assay = "dorothea") %>%
  data.frame() %>%
  t()
dim(viper_scores_df)

# Create a df containing cells and clusters
cells_clusters = data.frame(cell = names(Idents(rpca_smc_subset_v3)),
                            cell_type = rpca_smc_subset_v3$prelim_annotations)

# Create a df containing viper scores per cell and their respeective cell type annotations
viper_scores_clusters = viper_scores_df %>%
  data.frame() %>%
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(cells_clusters)
head(viper_scores_clusters)
table(viper_scores_clusters$tf)

# Summarize viper scores per cell type
summarized_viper_scores = viper_scores_clusters %>%
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))
head(summarized_viper_scores)

# For visualization purposes, select the 20 most variable TFs across clusters according to viper scores
highly_var_tfs = summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg)) %>%
  ungroup() %>%
  top_n(180, var) %>%
  distinct(tf)
highly_var_tfs

# Prepare data for plot
summarized_viper_scores_df = summarized_viper_scores %>%
  semi_join(highly_var_tfs, by="tf") %>%
  dplyr::select(-std) %>%
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE)
palette_length = 100

my_breaks = c(seq(min(summarized_viper_scores_df), 0,
                  length.out = ceiling(palette_length/2) + 1),
              seq(max(summarized_viper_scores_df)/palette_length,
                  max(summarized_viper_scores_df),
                  length.out = floor(palette_length/2)))

my_color1 = colorRampPalette(c("Darkblue", "white", "yellow", "red"))(palette_length)
my_color2 = colorRampPalette(c("Darkblue", "skyblue", "yellow", "red"))(palette_length)
my_color3 = colorRampPalette(c("#053061", "#2171B5", "white", "#FFD92F", "#A50F15"))(palette_length) #A50F15. The other orange alternative is #FEB24C

rownames(summarized_viper_scores_df)
viper_hmap = pheatmap::pheatmap(t(summarized_viper_scores_df), fontsize = 14,
                                fontsize_row = 10,
                                color = my_color3 , breaks = my_breaks,
                                main = "DoRothEA (ABC)", angle_col = 45,
                                treeheight_col = 0, border_color = NA, cluster_cols = FALSE)
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure5/Fig5d_SMC_dorothea_regs_heatmap.pdf",
       plot = viper_hmap, width = 7, height = 9)

# Visualize TF activity scores in a different way
VlnPlot(rpca_smc_subset_v3, features = c("SRF", "TCF21", "SOX9", "KLF5", "SPI1", "SMAD3"), 
        pt.size = 0, group.by = "prelim_annotations", y.max = 5) & npg_scale_bars & theme(legend.position = "none")
FeaturePlot(rpca_smc_subset_v3, features = c("KLF4", "KLF5", "PAX6", "TFAP2C"), raster = FALSE) & custom_theme & new_scale3

saveRDS(rpca_smc_subset_v3, 
        "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/seurat_objects/whole_ref_v3/rpca_smc_subset_res0.9_v1_dorothea.rds")












