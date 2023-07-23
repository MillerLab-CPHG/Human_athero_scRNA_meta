library(Seurat)
library(scRNAutils)
library(tidyverse)
library(data.table)
library(ggsci)

#########################################################################################
# This script will plot results from the pilot scDRS analysis using the whole reference #
#########################################################################################

# Load meta-analyzed reference
rpca_int_sct_v3_1 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/seurat_objects/alsaigh_pan_wirka_hu_int_seurat_clustered_v3_1.rds")
DimPlot(rpca_int_sct_v3_1, group.by = "updated_level2_annotations", raster = FALSE, label = TRUE, label.size = 3) + 
  custom_theme() + theme(legend.position = "none") + level2_annotations_scale_new

# Plot by study
int_by_study = DimPlot(rpca_int_sct_v3_1, group.by = "study", split.by = "study", raster = FALSE, ncol = 2) + 
  custom_theme() + theme(legend.position = "none") + 
  theme(strip.background = element_rect(fill="white")) + 
  miller_discrete_scale()

# Save plot by study
ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure2/revision_figures/SuppFig2h_ref_by_study.png",
       int_by_study, width = 7, height = 7)


# Load scDRS outputs. In this case we have 250 control gene sets
cad_mvp_eur_meta = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_out_scores/CAD_MVP_EUR_meta/CAD_MVP_EUR_meta.full_score.gz")
wbc_count = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_out_scores/WBC_count/White_Blood_Cell_count.full_score.gz")
alzheimer_disease = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_out_scores/Alzheimer_disease/Alzheimer_disease.full_score.gz")

# Create list with traits
traits_vec = c("cad_mvp_eur_meta", "wbc_count", "alzheimer_disease")

# Add normscDRS scores into metadata
for (i in traits_vec) {
  trait = get(i)
  rpca_int_sct_v3_1@meta.data[[i]] = trait$norm_score
  
}

# Explore some of the stats
cad_pvals = cad_mvp_eur_meta %>%
  ggplot(aes(x=pval)) + 
  geom_histogram(bins = 100) + 
  #ggtitle("CAD MVP EUR meta") + 
  geom_vline(xintercept = 0.05, color="darkred") + 
  custom_theme()

cad_norm_scores = cad_mvp_eur_meta %>%
  ggplot(aes(x=norm_score)) + 
  geom_histogram(bins = 100) + 
  custom_theme()

wbc_pvals = wbc_count %>%
  ggplot(aes(x=pval)) + 
  #geom_density() + 
  geom_histogram(bins = 100) +
  #ggtitle("WBC count") + 
  geom_vline(xintercept = 0.05, color="darkred") + 
  custom_theme()

wbc_norm_scores = wbc_count %>%
  ggplot(aes(x=norm_score)) + 
  #geom_density() + 
  geom_histogram(bins = 100) + 
  custom_theme()

plot_list = list(wbc_pvals, wbc_norm_scores,
                 cad_pvals, cad_norm_scores)
scdrs_p_vals_norm_scores_distribution = cowplot::plot_grid(plotlist = plot_list)

# Save plot
ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Rebuttal_Fig9/Reb_Fig9b_scDRS_p_val_distribution.pdf",
       plot = scdrs_p_vals_norm_scores_distribution, width = 8, height = 8)


# Plot scDRS scores in UMAP embeddings
scdrs_plot = FeaturePlot(rpca_int_sct_v3_1, 
            features = c("cad_mvp_eur_meta",
                         "wbc_count",
                         "alzheimer_disease"), 
            order = TRUE,
            raster = FALSE, ncol = 3,
            pt.size = 0.1) & custom_theme() & scale_color_gradientn(colours = c("#053061", "white", "darkred"))

ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure4/revision_figures/Fig4e_scDRS_UMAP_norm_scores_AD.pdf",
       plot = scdrs_plot, width = 7, height = 7)




