library(tidyverse)
library(data.table)
library(ggsci)
library(ggpubr)
library(ggrepel)

#############################################################################
# This script  will correlate rna and protein levels for genes of interest  #
#############################################################################


# Write helper function to correlate RNA and proteomics data
# Args rna_df A dataframe of bulk RNA-seq (150 samples)
# Args proteo_df A dataframe output by the reshape_proteo_data() function with the proteins of interest to correlate
# Args gene_vec A character of ensembl IDs of interest that'll be used to filter the RNA-seq data. They should be named with the correpsnding gene symbols
# Value A dataframe where there is matching RNA and proteomics data for each patient
cor_rna_protein =  function(rna_df, proteo_df, gene_vec) { 
  
  # There are 47 samples from the proteomics disease category 1 vs 3 samples within 
  # the RNA-seq df.
  proteo_idx = which(proteo_df$UVA_ID %in% colnames(rna_df))
  disease_proteo_uva_ids = proteo_df$UVA_ID
  disease_proteo_uva_ids = disease_proteo_uva_ids[proteo_idx]
  
  # Filter RNA-seq and proteo df to contain only matching samples from proteomics data
  rna_df = rna_df[, c("Name", disease_proteo_uva_ids)]
  proteo_df_new = proteo_df[proteo_idx, ]
  
  # Filter RNA-seq matrix to contain only genes of interest 
  rna_df = rna_df[rna_df$Name %in% gene_vec,]
  rna_df$Gene_symbol = names(gene_vec)[match(rna_df$Name, gene_vec)]
  rownames(rna_df) = rna_df$Gene_symbol
  rna_df = rna_df %>%
    dplyr::select(-Name, -Gene_symbol)
  rna_df_t = t(rna_df)
  rna_df_t = as.data.frame(rna_df_t)
  rna_df_t$Sample = rownames(rna_df_t)
  target_genes_df = reshape2::melt(rna_df_t)
  names(target_genes_df) = c("UVA_ID", "Gene", "TPM")
  
  # Merge matched RNA-seq and proteomics samples
  #rna_proteo_merged_df = merge(proteo_df_new, target_genes_df, by = c("UVA_ID"))
  rna_proteo_merged_df = proteo_df_new %>%
    inner_join(target_genes_df, by = c("UVA_ID", "Gene"))
  dim(rna_proteo_merged_df)
  
  return(rna_proteo_merged_df)
  
}


##########################################################################
# Check for matching samples in RNA-seq data
# Check that we have matching samples in the bulk RNA-seq data (42 samples)
length(which(unique(proteo_interesting_hits_df$UVA_ID) %in% metadata$SUBJID))

# There are 47 samples from the proteomics disease category 1 vs 3 samples within 
# the RNA-seq df.
proteo_idx = which(unique(proteo_interesting_hits_df$UVA_ID) %in% colnames(patient_tpms))

# Define genes of interest to correlate 
ensembl_ids = c("ENSG00000095713", "ENSG00000049323")
gene_symbols = c("CRTAC1", "LTBP1" )
names(ensembl_ids) = gene_symbols

# Use helper function to create df with matching RNA and proteomics data
rna_proteo_cor = cor_rna_protein(patient_tpms, proteo_interesting_hits_df, ensembl_ids )

# Plot correlation results
rna_proteo_cor %>%
  #filter(UVA_ID != "UVA047") %>%
  ggplot(aes(x=TPM, y=Normalized_expression, label=UVA_ID)) +
  geom_smooth(method = "lm", color="black", se = TRUE) + 
  geom_point(aes(color=disease_category), size=2) +
  ylab("Norm. protein exp.") + 
  xlab("RNA TPMs") + 
  facet_wrap(~Gene, scales = c("free")) + 
  #geom_text_repel(size=3) + 
  custom_theme + 
  #scale_color_manual(values = c("#ED0000FF", "#00468BFF")) + 
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        strip.text = element_text(size=14, face = "italic"),
        strip.background = element_rect(fill="white")) + 
  npg_scale

# Calculate correlation between CRTAC1 protein and RNA expression
# cor r = 0.4938183;  p-value = 0.0004874
new_merged = rna_proteo_cor %>%
  filter(UVA_ID != "UVA047") %>%
  filter(Gene == "CRTAC1")
cor.test(new_merged$Normalized_expression, new_merged$TPM,
         method = "pearson")

# Calculate correlation between LTBP1 protein and RNA expression
# cor r = 0.3599226;  p-value = 0.01401
new_merged = rna_proteo_cor %>%
  filter(UVA_ID != "UVA047") %>%
  filter(Gene == "LTBP1")
cor.test(new_merged$Normalized_expression, new_merged$TPM,
         method = "pearson")


