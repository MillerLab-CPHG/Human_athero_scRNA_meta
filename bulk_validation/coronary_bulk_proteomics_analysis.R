library(tidyverse)
library(data.table)
library(ggsci)
library(RColorBrewer)


#######################
# Helper functions

# Write a function that will take gene names as input and return a df ready for plotting
# Args proteo_df A dataframe where rows are proteins and columns sample names
# Args genes A vector with genes of interest for plotting

reshape_proteo_data = function(proteo_df, genes) { 
  proteo_data_subset = proteo_df[proteo_df$GeneSymbol %in% genes, ]
  head(proteo_data_subset)
  rownames(proteo_data_subset) = proteo_data_subset$GeneSymbol
  proteo_data_subset = proteo_data_subset %>% dplyr::select(-GeneSymbol)
  proteo_data_subset_t = as.data.frame(t(proteo_data_subset))
  proteo_data_subset_t$Sample_ID = rownames(proteo_data_subset_t)
  proteo_subset_reshaped = suppressWarnings(melt(proteo_data_subset_t)) 
  names(proteo_subset_reshaped) = c("Sample_ID", "Gene", "Normalized_expression")
  
  # Introduce desired metadata into reshaped df
  proteo_subset_reshaped = proteo_subset_reshaped %>%
    mutate(disease_category = case_when(Sample_ID %in% cat1_samples_ids ~ "1",
                                        Sample_ID %in% cat3_samples_ids ~ "3")) 
  
  proteo_subset_reshaped$disease_category = factor(proteo_subset_reshaped$disease_category, 
                                                   levels = c("1", "3"))
  
  # Match UVA IDs to proteomics sample IDs
  idx = match(proteo_subset_reshaped$Sample_ID, protein_meta$KCL_ID2)
  uva_ids = protein_meta$UVA_ID
  uva_ids_disease_status = uva_ids[idx]
  proteo_subset_reshaped$UVA_ID = uva_ids_disease_status
  
  
  return(proteo_subset_reshaped)
}


#####################################
# Load metadata for disease status
protein_meta = fread("~/Desktop/human_athero_scRNA_meta-analysis/coronary_proteomics_data/protein_coronary_metadata_092122.csv", quote = FALSE)
head(protein_meta)

# Only 43 samples are classified with disease status (Diseased=29, Healthy=14)
table(protein_meta$`Healthy/Diseased`)

#############################################################
# Plot CRTAC1 expression according to Healthy/Diseased status
# healthy_samples = protein_meta[protein_meta$`Healthy/Diseased` == "Healthy"]
# healthy_samples_ids = healthy_samples$KCL_ID2
# 
# diseased_samples = protein_meta[protein_meta$`Healthy/Diseased` == "Diseased"]
# diseased_samples_ids = diseased_samples$KCL_ID2

#############################################################
# Plot CRTAC1 expression according to Category 1 vs 3 status
# We need KCL_ID2 because these are the col names of the proteomics df
# A total of 56 samples are classified as either category 1 or 3
# Category 1 n= 27
cat1_samples = protein_meta[protein_meta$Category == 1]
cat1_samples_ids = cat1_samples$KCL_ID2
length(cat1_samples_ids)

# Category 3 n=29
cat3_samples = protein_meta[protein_meta$Category == 3]
cat3_samples_ids = cat3_samples$KCL_ID2
length(cat3_samples_ids)

#############################################################
# Load log2normalized protein expression
proteomics_data = fread("~/Desktop/human_athero_scRNA_meta-analysis/coronary_proteomics_data/Coronary_proteomics_data.txt")
proteomics_data = as.data.frame(proteomics_data)
head(proteomics_data)

# Filter proteomics data to keep samples with healthy/diseased labels
sample_names = names(proteomics_data)
# proteomics_data_filtered = proteomics_data[, c("GeneSymbol", healthy_samples_ids, 
#                                                diseased_samples_ids)]

# Filter proteomics data to keep samples with category 1 and 3 labels
proteomics_data_filtered = proteomics_data[, c("GeneSymbol", 
                                               cat1_samples_ids, 
                                               cat3_samples_ids)]

# 56 samples (cat 1 and 3)
dim(proteomics_data_filtered)

####################################
# Prep inputs for plotting function
interesting_genes = c("LTBP1", "TNFRSF11B", "VCAN")
proteo_interesting_hits_df = reshape_proteo_data(proteomics_data_filtered, interesting_genes)
proteo_interesting_hits_df$Gene = factor(proteo_interesting_hits_df$Gene, levels = c("LTBP1",
                                                                                     "TNFRSF11B",
                                                                                     "VCAN"))
interesting_genes = c("CRTAC1", "SPP1", "LTBP1", "VCAN")
proteo_interesting_hits_df = reshape_proteo_data(proteomics_data_filtered, interesting_genes)
proteo_interesting_hits_df$Gene = factor(proteo_interesting_hits_df$Gene, levels = c("CRTAC1", "SPP1", "LTBP1", "VCAN"))

interesting_genes = c("CRTAC1", "LTBP1")
proteo_interesting_hits_df = reshape_proteo_data(proteomics_data_filtered, interesting_genes)
proteo_interesting_hits_df$Gene = factor(proteo_interesting_hits_df$Gene, levels = c("CRTAC1", "LTBP1"))

# Compare protein expression across disease status
proteo_interesting_hits_df %>%
  ggplot(aes(x=disease_category, y=Normalized_expression)) + 
  geom_boxplot(aes(x=disease_category, y=Normalized_expression, fill=Gene), 
               outlier.shape = NA) + 
  geom_point(aes(x=disease_category, y=Normalized_expression, color=Gene), size=3,
             position = position_dodge(width = 0.77)) + 
  ylab("Protein norm. expression") + 
  xlab("Disease category") + 
  #geom_dotplot(aes(fill=Gene, color=Gene), binaxis='y', 
  #             stackdir='center', dotsize=1) +
  #scale_fill_lancet() +
  #scale_color_lancet() + 
  #scale_color_manual(values = c("#ED0000FF", "#00468BFF")) +
  scale_fill_manual(values = c(rep("white", 4))) + 
  #facet_wrap(~Gene, scales = "free_y", ncol = 3) + 
  theme_bw() + 
  custom_theme + 
  theme(aspect.ratio = 1.3,
        legend.position = "right",
        strip.text = element_text(size=14),
        strip.background = element_rect(fill="white")) + 
  npg_scale 

ggsave("~/Desktop/human_athero_scRNA_meta-analysis/manuscript_figures/Figure5/Fig5e_CRTAC1_protein_exp_cat_1vs3.pdf",
       plot = crtac1_proteo_plot, width = 7, height = 7)



########################
# Statistical analyses #
########################

# Check assumptions for calculating the p value
crtac1_proteo_df %>%
  ggplot(aes(sample=CRTAC1_norm_expression, color=Category)) +
  geom_qq(size=0.5) +
  geom_qq_line(color="black") +
  facet_wrap(~Category, scales = "free_y") + 
  custom_theme + 
  scale_color_manual(values = c("#ED0000FF", "#00468BFF")) +
  theme(aspect.ratio = 1,
        legend.position = "none",
        strip.text = element_text(size=14),
        strip.background = element_rect(fill="azure3"))

# Check for variance
crtac1_proteo_df %>%
  ggplot(aes(x=CRTAC1_norm_expression, fill=Category)) + 
  geom_density(alpha=0.5) +
  custom_theme + 
  scale_fill_manual(values = c("#ED0000FF", "#00468BFF")) 


# CRTAC1 p value = 4.29e-07
t.test(Normalized_expression ~ disease_category, 
       data = proteo_interesting_hits_df[proteo_interesting_hits_df$Gene == "CRTAC1",])

# SPP1 p  value = 7.731e-11
t.test(Normalized_expression ~ disease_category, 
       data = proteo_interesting_hits_df[proteo_interesting_hits_df$Gene == "SPP1",])

# LTBP1 p value = 0.009616
t.test(Normalized_expression ~ disease_category, 
       data = proteo_interesting_hits_df[proteo_interesting_hits_df$Gene == "LTBP1",])

# VCAN p value = 2.124e-07
t.test(Normalized_expression ~ disease_category, 
       data = proteo_interesting_hits_df[proteo_interesting_hits_df$Gene == "VCAN",])


# Calculate fold changes
crtac1_proteo_means = crtac1_proteo_df %>%
  group_by(Disease_status) %>%
  summarize(CRTAC1_mean_norm_expression = mean(CRTAC1_norm_expression))

crtac1_proteo_log2fc = log2(crtac1_proteo_means$CRTAC1_mean_norm_expression[2] / crtac1_proteo_means$CRTAC1_mean_norm_expression[1])






