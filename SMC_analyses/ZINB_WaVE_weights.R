library(Seurat)
library(zinbwave)
library(DESeq2)
library(scRNAutils)
library(scran)
library(tidyverse)
library(data.table)
library(ggrepel)
library(BiocParallel)


####################################################################################
# The goal of this script is to generate the observation weights from ZINB-WaVE to #  
# run DESeq2 analyses
####################################################################################


# Turn off the parallelization
register(SerialParam())


##########################################
# Run DE analysis for meta-analyzed SMCs #


# Load meta-analyzed SMC seurat object
# Read subclustered seurat object
rpca_smc_fibro_subset_v3 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/seurat_objects/whole_ref_v3/peri_smc_fibro_subset_seurat_obj_res0.9_v1.rds")
DefaultAssay(rpca_smc_fibro_subset_v3) = "SCT"

# Further subset to keep only SMCs
rpca_smc_subset_v3 = subset(rpca_smc_fibro_subset_v3, idents = c(5, 2, 16, 4, 0 ,9, 6, 13, 15, 17))
DimPlot(rpca_smc_subset_v3)

# Add metadata for sex
# Sex was defined based on metadata from papers and also XIST expression. 
males = c("pan_rpe004", "pan_rpe005", "alsaigh_ac_p1", "alsaigh_pa_p1", "alsaigh_ac_p2", "alsaigh_pa_p2",
          "alsaigh_ac_p3", "alsaigh_pa_p3", "wirka_coronary_1", "wirka_coronary_2", "wirka_coronary_3",
          "wirka_coronary_4", "wirka_coronary_5", "wirka_coronary_6", "wirka_coronary_7", "hu_coronary1_p1", 
          "hu_coronary1_p2", "hu_coronary2_p2")
females = c("pan_rpe006", "wirka_coronary_8", "hu_coronary1_p3", "hu_coronary2_p3")

# Update metadata
rpca_smc_subset_v3@meta.data = rpca_smc_subset_v3@meta.data %>%
  mutate(sex = case_when(sample %in% males ~ "males",
                         sample %in% females ~ "females"),
         sample_disease_status = case_when(sample_disease_status == "non_diseased" ~ "non_lesion",
                                           TRUE ~ "lesion"))
rownames(rpca_smc_subset_v3@meta.data) = barcodes
names(rpca_smc_subset_v3@meta.data)

# Re factor lesion levels
rpca_smc_subset_v3$sample_disease_status = factor(rpca_smc_subset_v3$sample_disease_status,
                                                  levels = c("non_lesion", "lesion"))

# Find variable features
DefaultAssay(rpca_smc_subset_v3) = "RNA"
rpca_smc_subset_v3 = NormalizeData(rpca_smc_subset_v3,
                                   normalization.method = "LogNormalize",
                                   scale.factor = 10000)
rpca_smc_subset_v3 = FindVariableFeatures(rpca_smc_subset_v3, 
                                          selection.method = "vst",
                                          nfeatures = 3000, assay = "RNA")
rpca_smc_var_features = rpca_smc_subset_v3@assays$RNA@var.features

# Make SE object
rpca_smc_se = Seurat::as.SingleCellExperiment(rpca_smc_subset_v3)
class(rpca_smc_se)

# First subset to keep only variable features
# This SE object has 21148 cells x 3000 highly variable genes
rpca_smc_se_var_feats = rpca_smc_se[rpca_smc_var_features,]

###########################################################
# Model zero component using ZINB-WaVE

# Filter out lowly expressed genes by removing genes that
# do not have at least 5 reads in at least 5 samples
keep = rowSums(assay(rpca_smc_se_var_feats) >= 5) >= 5

# By filtering, we're keeping 1231 genes
table(keep)
rpca_smc_se_filtered = rpca_smc_se_var_feats[keep,]

# Make sure that input counts matrix is of class matrix
assay(rpca_smc_se_filtered) = round(as.matrix(assay(rpca_smc_se_filtered)))
rpca_smc_se_filtered

# Fit model and get observation weights (Takes around 82 mins setting BPPARAM=SerialParam())
# Get ZINB-WaVE weights using the same design matrix as for DESEq2
# Add sex and arterial origin as covariates
system.time({
  rpca_smc_zinb = zinbwave(rpca_smc_se_filtered, K=0, 
                           observationalWeights=TRUE,
                           X = "~ sex + arterial_origin + sample_disease_status",
                           BPPARAM=SerialParam())
})

# Save SE object with ZINB-WaVE weights
saveRDS(rpca_smc_zinb,
        "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/DE_analyses/rPCA_SMC_ZINB_WaVE_weights_sex_artery_lesion_covars.rds")

rpca_smc_zinb = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/DE_analyses/rPCA_SMC_ZINB_WaVE_weights_sex_artery_lesion_covars.rds")

# Check that we have all cells
table(rpca_smc_zinb@colData$sample_disease_status)

#########################################
# Run DE analysis
# Adjust for same covariates in scDRS analyses: sex, arterial origin and disease status
dds = DESeqDataSet(rpca_smc_zinb,
                   design = ~ sex + arterial_origin + sample_disease_status)

# Re factor lesion categories
dds@colData$sample_disease_status = factor(dds@colData$sample_disease_status, 
                                           levels = c("non_lesion", "lesion"))

# Estimate size factors
dds = estimateSizeFactors(dds, type="poscounts")

# Use scran's sum factors
#scr = scran::computeSumFactors(dds)
#sizeFactors(dds) = sizeFactors(scr)

# Need to set the right parameters for single cell analysis
dds = DESeq(dds, test="LRT", minmu = 1e-6, 
            useT=FALSE,
            reduced=~1, 
            minReplicatesForReplace = Inf)

# Get results
#res = results(dds, contrast=c("sample_disease_status", "lesion", "non_lesion"))
res = DESeq2::results(dds, contrast = c("sample_disease_status",
                                        "lesion", "non_lesion"))
res_df = as.data.frame(res)
res_df$gene = rownames(res_df)

# Set threshold for filtering (FDR < 0.05 and log2FC > 1)
res_df_new = res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 0.9 & baseMean > 0.3) %>%
  arrange(desc(log2FoldChange))
head(res_df_new, n=100)
tail(res_df_new, n=100)
dim(res_df_new)

# Do a quick volcano plot
res_df_new %>%
  filter(abs(log2FoldChange) < 5) %>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj), label=gene)) + 
  geom_point() + 
  geom_vline(xintercept = c(-1, 0, 1), color="darkred", linetype="dashed") + 
  geom_hline(yintercept = -log10(1e-10), color="darkred", linetype="dashed") + 
  geom_text_repel() + 
  custom_theme()

# Save results table
write.table(res_df_new,
            "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_tables/DESe2_lesion_status_SMCs.tsv",
            quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

# Plot dispersion estimates
disp_estimates = plotDispEsts(dds)
class(disp_estimates)

# Make MA plot
plotMA(dds, ylim=c(-5, 5))

############################################################################
# GWAS genes

# Load effector genes 
effector_genes = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/whole_ref/CELLECT-GENES/results/effector_genes.csv")

# Filter for CAD effecor genes in SMCs
cad_smc_genes_df = effector_genes %>%
  filter(gwas == "CAD_MVP_EUR_meta" & annotation == "SMC")

cad_smc_genes_vec = cad_smc_genes_df$gene_symbol
length(cad_smc_genes_vec)

# Check how many CAD genes in SMCs overlap DE genes
# Genes fulfilling this criteria: MYH11, PDE5A, FHL5, MYH11, COL4A1/2
idx = which(cad_smc_genes_vec %in% res_df_new$gene)
cad_smc_genes_vec[idx]

# Define vector with genes of interest
genes_vec = c("CDH13", "PDE5A", "FHL5", "MYH11", "COL4A1", "COL4A2")
markers_vec = c("LTBP1", "CRTAC1")

# Write a function to speed this up
plot_deseq_counts_df = function(dds, gene_id) { 
  d = plotCounts(dds, gene=gene_id, intgroup = "sample_disease_status",
                 returnData = TRUE)
  d$gene_name = gene_id
  return(d)
  }

# Plot counts for all genes of interest 
genes_list = lapply(markers_vec, plot_deseq_counts_df, 
                    dds=dds)

genes_df = data.table::rbindlist(genes_list)

gwas_smc_de_genes = genes_df %>%
  ggplot(aes(x=sample_disease_status, y=count)) + 
  geom_violin() + 
  #geom_point() + 
  geom_jitter(aes(color=sample_disease_status), height = 0, width = 0.05, size=0.3, alpha=0.1) +
  facet_wrap(~ gene_name, scales = c("free_y")) + 
  labs(x="Disease category", y="Normalized expression") + 
  scale_y_log10() + 
  custom_theme() + 
  theme(strip.background = element_rect(fill="white"),
        legend.position = "none") + 
  miller_discrete_scale()

ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/SMC_DE_figure/ltbp1_crtac1_lesion_DE_plot.pdf",
       plot = gwas_smc_de_genes, width = 6, height = 4)


#######################################################################
# Shrink LFCs
res_lfc_shrunken = lfcShrink(dds, contrast = c("sample_disease_status", 
                                               "lesion", "non_lesion"),
                             type = "normal", lfcThreshold = 2)

# Plot before and after shrinking
unshrunk_ma = plotMA(dds, ylim=c(-5, 5))
shrunken_ma = plotMA(res_lfc_shrunken, ylim=c(-5, 5))

# Inspect results after LFC shrinking
res_lfc_shrunk = as.data.frame(res_lfc_shrunken)
head(res_lfc_shrunk, n=40)

res_new = res_lfc_shrunk %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange))
head(res_new, n=100)

# Do a quick volcano plot 
res_lfc_shrunk %>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point() + 
  custom_theme()
















