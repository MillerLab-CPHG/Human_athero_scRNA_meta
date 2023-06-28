library(Seurat)
library(scRNAutils)
library(tidyverse)
library(data.table)
library(ggsci)

####################################################################################
# This script will plot results from the pilot scDRS analysis using the Wirka data #
####################################################################################

# Load Wirka data
rpca_int_sct_v3_1 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/seurat_objects/alsaigh_pan_wirka_hu_int_seurat_clustered_v3_1.rds")

# Load scDRS outputs. In this case we have 250 control gene sets
cad_van_der_harst = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_out_scores/CAD_van_der_Harst/CAD_van_der_Harst.full_score.gz")
cad_mvp_eur_meta = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_out_scores/CAD_MVP_EUR_meta/CAD_MVP_EUR_meta.full_score.gz")
wbc_count = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_out_scores/WBC_count/White_Blood_Cell_count.full_score.gz")
alzheimer_disease = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_out_scores/Alzheimer_disease/Alzheimer_disease.full_score.gz")
carotid_plaque = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_out_scores/Carotid_plaque/Carotid_Plaque.full_score.gz")

# Create list with traits
traits_vec = c("cad_van_der_harst" , "cad_mvp_eur_meta", "wbc_count",
               "alzheimer_disease", "carotid_plaque")

# Add normscDRS scores into Wirka metadata
for (i in traits_vec) {
  trait = get(i)
  rpca_int_sct_v3_1@meta.data[[i]] = trait$norm_score
  
}

# Explore some of the stats
head(cad_van_der_harst)

cad_van_der_harst %>%
  ggplot(aes(x=pval)) + 
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = 0.05) + 
  custom_theme()

wbc_count %>%
  ggplot(aes(x=norm_score)) + 
  #geom_density() + 
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = 0.05) + 
  custom_theme()
  

# Plot manual annotations
DimPlot(rpca_int_sct_v3_1, group.by = "manually_annotated_labels", 
        label = TRUE, repel = TRUE) + custom_theme()  

# Plot scDRS scores 
FeaturePlot(rpca_int_sct_v3_1, features = c("cad_van_der_harst", 
                                            "cad_mvp_eur_meta",
                                            "carotid_plaque"), 
            order = TRUE,
            raster = FALSE,
            pt.size = 0.1) & custom_theme() &  scale_color_gradientn(colours = c("#053061", "white", "darkred"))

FeaturePlot(rpca_int_sct_v3_1, features = c("wbc_count",
                                            "alzheimer_disease"), 
             order = TRUE,
             raster = FALSE, ncol = 1,
             pt.size = 0.1) & custom_theme() &  scale_color_gradientn(colours = c("#053061", "white", "darkred"))

plot_list = scRNAutils::plot_gene_UMAP_list(wirka_data, genes_of_interest = traits_vec, pt.size = 0.3, target_assay = "RNA") & scale_color_gradient2(low = "navy", mid="white", high = "red")
cowplot::plot_grid(plotlist = plot_list)


