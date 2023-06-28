library(Seurat)
library(scRNAutils)
library(tidyverse)
library(data.table)
library(ggsci)

####################################################################################
# This script will plot results from the pilot scDRS analysis using the Wirka data #
####################################################################################

# Load Wirka data
wirka_data = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scRNA_raw_matrices/Wirka_human_coronaries_scRNA_data/past_analyses/rds_objects/scRNA_athero_data_seurat_manual_labels1_res1.7.rds")

# Load scDRS outputs
cad_van_der_harst = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_pilot/scDRS_out_scores/CAD_van_der_Harst/CAD_van_der_Harst.full_score.gz")
cad_mvp_eur_meta = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_pilot/scDRS_out_scores/CAD_MVP_EUR_meta/CAD_MVP_EUR_meta.full_score.gz")
wbc_count = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_pilot/scDRS_out_scores/WBC_count/White_Blood_Cell_count.full_score.gz")
alzheimer_disease = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_pilot/scDRS_out_scores/Alzheimer_disease/Alzheimer_disease.full_score.gz")
carotid_plaque = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_pilot/scDRS_out_scores/Carotid_plaque/Carotid_Plaque.full_score.gz")

# Create list with traits
traits_vec = c("cad_van_der_harst", "cad_mvp_eur_meta",
                   "wbc_count", "alzheimer_disease", "carotid_plaque")

# Add normscDRS scores into Wirka metadata
for (i in traits_vec) {
  trait = get(i)
  wirka_data@meta.data[[i]] = trait$norm_score
}

# Plot manual annotations
wirka_clusters = DimPlot(wirka_data, group.by = "manually_annotated_labels", label = TRUE, repel = TRUE) + custom_theme()  

# Plot scDRS scores 
FeaturePlot(wirka_data, features = c("cad_mvp_eur_meta", "carotid_plaque", "wbc_count", "alzheimer_disease"), order = TRUE,
             pt.size = 0.9) & custom_theme() & scale_color_gradientn(colours = c("#053061", "white", "darkred")) 

plot_list = scRNAutils::plot_gene_UMAP_list(wirka_data, genes_of_interest = traits_vec, pt.size = 0.3, target_assay = "RNA") & scale_color_gradient2(low = "navy", mid="white", high = "red")
cowplot::plot_grid(plotlist = plot_list)


