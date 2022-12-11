library(tidyverse)
library(data.table)
library(ggsci)



# Helper functions

# Write a function to reshape RNA seq data and prep for plotting
# Arg rna_df A matrix where rows are gene names and columns are samples with norm. expression data 
# Arg gene_vec A vector with ensembl IDs with the respective gene symbols as names
# Arg metadata A df with samples metadata
# Value A df with normalized expression for genes of interest that can be used for plotting
reshape_rna_data = function(rna_df, gene_vec, metadata) {
  
  # Take the main matrix and filter to keep only genes of interest
  rna_df = rna_df[rna_df$Name %in% gene_vec,]
  rna_df$Gene_symbol = names(gene_vec)[match(rna_df$Name, gene_vec)]
  rownames(rna_df) = rna_df$Gene_symbol
  rna_df = rna_df %>%
    dplyr::select(-Name, -Gene_symbol)
  rna_df_t = t(rna_df)
  rna_df_t = as.data.frame(rna_df_t)
  rna_df_t$Sample = rownames(rna_df_t)
  target_genes_df = reshape2::melt(rna_df_t)
  names(target_genes_df) = c("Sample", "Gene", "TPM")
  
  # Add details about samples disease status
  disease_cat = metadata$Diseased[match(target_genes_df$Sample, 
                                        filtered_metadata$SUBJID)]
  
  target_genes_df$disease_cat = disease_cat
  target_genes_df = target_genes_df %>%
    mutate(disease = case_when(disease_cat == 1 ~ "Diseased",
                               disease_cat == 0 ~ "Healthy"))
  target_genes_df$disease = factor(target_genes_df$disease, 
                                   levels = c("Healthy", "Diseased"))
  return(target_genes_df)
}

# This is a function to generate qqplots using the output df produced by the function above
# Arg target_genes_df A df output by the reshape_rna_data() function  
# Arg gene A character with the target gene for the qqplots
# Value A ggplot object sshowing qqplots for genes of interest

make_qqplot = function(target_genes_df, gene) {
  qq_plots = target_genes_df %>%
    filter(Gene == gene) %>%
    ggplot(aes(sample=TPM, color=disease)) +
    geom_qq(size=0.5) +
    geom_qq_line(color="black") +
    facet_wrap(~disease, scales = "free_y") + 
    custom_theme + 
    scale_color_manual(values = c("#ED0000FF", "#00468BFF")) +
    theme(aspect.ratio = 1,
          legend.position = "none",
          strip.text = element_text(size=14),
          strip.background = element_rect(fill="azure3"))
  return(qq_plots)
  }



#####################################
# Read metadata for RNA-seq samples 
metadata = fread("~/Desktop/human_athero_scRNA_meta-analysis/coronary_bulk_RNA_seq/QTL_covariates_disease_status.txt")
names(metadata)
table(metadata$Lesion)
table(metadata$Classification)
table(metadata$Diseased)
table(metadata$Tissue)
glimpse(metadata)

disease_status = metadata %>%
  select(SUBJID, Diseased)


# Look at samples of interest 
diseased_samples = c("UVA069", "UVA048", "UVA108", "UVA176", "UVA098", "UVA068")
healthy_samples = c("UVA157", "UVA005", "UVA154", "UVA156", "UVA026", "UVA007", "UVA138")
metadata[metadata$SUBJID %in% diseased_samples, ]
metadata[metadata$SUBJID %in% healthy_samples, ]

################################################################################
# Filter metadata file to keep only samples that have an annotated lesion status
filtered_metadata = metadata %>%
  filter(Lesion %in% c("N", "D"))
dim(filtered_metadata)

# There are 27 healthy samples and and 21 with lesions
table(filtered_metadata$Diseased)

#############################################################
# Load TPM normalized counts from lesion vs healthy patients
patient_tpms = fread("~/Desktop/human_athero_scRNA_meta-analysis/coronary_bulk_RNA_seq/TPM_filtered_invn_UVAcoronary_ENSGs_July2020.txt")
patient_tpms = as.data.frame(patient_tpms)
rownames(patient_tpms) = patient_tpms$Name
head(patient_tpms)
dim(patient_tpms)

##################################################################################
# Check which samples from the bulk RNA-seq are included within that metadata file 
smpl_idx = which(names(patient_tpms) %in% filtered_metadata$SUBJID)
length(names(patient_tpms)[smpl_idx])

patient_tpms_filtered = patient_tpms[, c(1,smpl_idx)]
dim(patient_tpms_filtered)
head(patient_tpms_filtered)

#################################################################
# Use helper function to compare expression across disease status

# Prep inputs for plotting function

ensembl_ids = c( "ENSG00000049323", "ENSG00000164761", "ENSG00000038427", "ENSG00000118526")
gene_symbols = c("LTBP1", "TNFRSF11B", "VCAN", "TCF21")
names(ensembl_ids) = gene_symbols

ensembl_ids2 = c("ENSG00000095713", "ENSG00000029559", "ENSG00000111341", "ENSG00000118785", 
                 "ENSG00000162551", "ENSG00000105664")
gene_symbols2 = c("CRTAC1", "IBSP", "MGP", "SPP1", "ALPL", "COMP")
names(ensembl_ids2) = gene_symbols2

ensembl_ids3 = c("ENSG00000095713", "ENSG00000029559", "ENSG00000049323", "ENSG00000118526")
gene_symbols3 = c("CRTAC1", "IBSP", "LTBP1", "TCF21")
names(ensembl_ids3) = gene_symbols3

# Process rna exp matrix with helper function
target_genes_df = reshape_rna_data(patient_tpms_filtered, ensembl_ids2, filtered_metadata)
target_genes_df2 = reshape_rna_data(patient_tpms_filtered, ensembl_ids, filtered_metadata)
target_genes_df3 = reshape_rna_data(patient_tpms_filtered, ensembl_ids3, filtered_metadata)


# Plot TPMs across disease status
target_genes_df3$Gene = factor(target_genes_df3$Gene, levels = c("CRTAC1", "IBSP", "LTBP1", "TCF21"))
crtac1_ibsp_tpms = target_genes_df3 %>%
  ggplot(aes(x=disease, y=TPM)) + 
  geom_boxplot(aes(x=disease, y=TPM), outlier.shape = NA) + 
  geom_point(aes(x=disease, y=TPM, color=disease), size=3,
             position = position_dodge(width = 0.77)) + 
  ylab("Norm expression (TPM)") + 
  xlab("Disease category") + 
  facet_wrap(~Gene, scales = "free_y", ncol = 2) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + 
  #scale_color_manual(values = c("#ED0000FF", "#00468BFF")) +
  custom_theme + 
  scale_fill_manual(values = c(rep("white", 4))) +
  theme(aspect.ratio = 1,
        legend.position = "none",
        strip.text = element_text(size=14, face = "italic"),
        strip.background = element_rect(fill="white")) +
  npg_scale
        

ggsave(file="~/Desktop/human_athero_scRNA_meta-analysis/manuscript_figures/Figure7/Fig7c_CRTAC1_IBSP_RNA-seq_tpms.pdf",
       plot = crtac1_ibsp_tpms, width = 8, height = 7)

##################################################
# Calculate statistics for genes of interest
# Run t test for CRTAC1 (Check assumtpions first)
make_qqplot(target_genes_df3, "IBSP")

crtac1_var_plot = target_genes_df %>%
  filter(Gene == "CRTAC1") %>%
  ggplot(aes(x=TPM, fill=disease)) + 
  geom_density(alpha=0.5) + 
  custom_theme + 
  scale_fill_manual(values = c("#ED0000FF", "#00468BFF"))

# Since TPMs don't seem normally distributed, run a Wilcoxon Mann U Whitney test
# CRTAC1 Wilcox t test p value = 5.429E-05
wilcox.test(TPM ~ disease, 
            data = target_genes_df3[target_genes_df3$Gene == "CRTAC1",])

# IBSP Wilcox t test p value = 7.661E-05
wilcox.test(TPM ~ disease, 
            data = target_genes_df3[target_genes_df3$Gene == "IBSP",])

# LTBP1 Wilcox t test p value = 9.754E-05
wilcox.test(TPM ~ disease, 
            data = target_genes_df3[target_genes_df3$Gene == "LTBP1",])

# TCF21 Wilcox t test p value = 0.05536 
wilcox.test(TPM ~ disease, 
            data = target_genes_df3[target_genes_df3$Gene == "TCF21",])

# Calculate fold change for CRTAC1
target_genes_df3 %>%
  group_by(disease) %>%
  summarise(median(TPM), min(TPM), max(TPM))




