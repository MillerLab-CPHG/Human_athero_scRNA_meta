library(tidyverse)
library(data.table)
library(ggsci)
library(scRNAutils)

################################################################################
# The goal of this script is to plot results from the LDSC and MAGMA analyses  #
################################################################################

# Reorder GWAS traits
ordered_gwas_traits = c("CAD_MVP_EUR_meta", "CAD_GWAS_van_der_Harst", "CAC_GWAS", "Myocardial_infarction", "Carotid_plaque", "Carotid_IMT",
                        "Pulse_pressure_MVP", "Diastolic_blood_pressure_MVP", "Systolic_blood_pressure_MVP", "Alzheimer_disease_Jansen",
                        "White_blood_cell_count", "Body_mass_index", "Type_2_diabetes")
cardiovascular_traits = c("CAD_MVP_EUR_meta", "Myocardial_infarction", "CAC_GWAS", "Carotid_plaque", 
                          "Alzheimer_disease_Jansen", "White_blood_cell_count")
cardiovascular_traits_supp = c("CAD_GWAS_van_der_Harst", "Pulse_pressure_MVP", "Body_mass_index", "Type_2_diabetes")

cardiovascular_traits_3 = c("CAD_MVP_EUR_meta", "Myocardial_infarction", "CAC_GWAS", "Carotid_plaque", 
                          "Alzheimer_disease_Jansen", "Type_2_diabetes")

# Select representative CAD and control traits to include in manuscript
gwas_traits_of_interest = c("CAD_MVP_EUR_meta", "CAD_GWAS_van_der_Harst", "CAC_GWAS", "Myocardial_infarction", 
                            "Carotid_plaque", "Pulse_pressure_MVP",
                            "Alzheimer_disease_Jansen", "White_blood_cell_count", 
                            "Body_mass_index", "Type_2_diabetes")

# Reorder cell annotations
ordered_level1_anno = c("B_cell", "Endothelial", "Fibroblast", "SMC",
                        "Mast_cell", "Neuron", "Pericyte", "Macrophage",
                        "Plasma_cell", "pDC", "T_NK")
ordered_level2_smc = c("Contractile_SMC", "Transitional-SMC", "Fibromyocyte",
                       "Fibrochondrocyte", "Fibroblast", "SMC2", "SMC3", "Foam-like", 
                       "Pericyte1", "Pericyte2")

#########################################
# LDSC from combined level1 annotations #
#########################################

# Load results for level1 LDSC and MAGMA analyses
ldsc_level1_prioritization_df = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/whole_ref/CELLECT-LDSC/results/prioritization.csv")
head(ldsc_level1_prioritization_df)

magma_level1_prioritization_df = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/whole_ref/CELLECT-MAGMA/results/prioritization.csv")
head(magma_level1_prioritization_df)

# Re factor GWAS traits and annotations
methods_vec = c("ldsc_level1_prioritization_df", 
                "magma_level1_prioritization_df")
for (i in methods_vec) {
  df = get(i)
  df$gwas = factor(df$gwas, 
                   levels = ordered_gwas_traits)
  df$annotation = factor(df$annotation, 
                         levels = ordered_level1_anno)
}

# Rename cols and prep df for joint LDSC and MAGMA plotting 
ldsc_level1_df = ldsc_level1_prioritization_df %>%
  dplyr::rename(LDSC = pvalue) %>%
  select(gwas, annotation, LDSC)

magma_level1_df = magma_level1_prioritization_df %>%
  dplyr::rename(MAGMA = pvalue) %>%
  select(gwas, annotation, MAGMA)

level1_merged_df = merge(ldsc_level1_df, magma_level1_df,
                         by.x = c("gwas", "annotation"), 
                         by.y = c("gwas", "annotation"))
level1_merged_df = reshape2::melt(level1_merged_df, 
                                   id=c("gwas", "annotation"))
names(level1_merged_df)[c(3, 4)] = c("method", "pvalue")
level1_merged_df$gwas = factor(level1_merged_df$gwas,
                               levels = ordered_gwas_traits)


# Plot results as dot plot (deprecated in manuscript)
ldsc_level1_results = ldsc_level1_prioritization_df %>%
  filter(gwas %in% cardiovascular_traits_supp) %>%
  ggplot(aes(x=annotation, y=-log10(pvalue), color=annotation,
             label=ifelse(-log10(pvalue) >= 1.3, annotation, ""))) + 
  geom_point(aes(size=-log10(pvalue) >= 1.3)) + 
  #geom_text_repel(color="black") + 
  geom_hline(yintercept = 1.3, color="darkred", linetype="dashed") + 
  scale_size_manual(values = c(2,4)) + 
  ggtitle("LDSC results from level 1 annotations") + 
  facet_wrap(~gwas, scales = "free_y", ncol = 4) + 
  custom_theme() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size=14)) + 
  miller_discrete_scale() +
  guides(size=FALSE,
         color = guide_legend(override.aes = list(size=4)))

ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure1/Fig1e_level1_ldsc_new.pdf",
       plot=level1_ldsc, width = 12, height = 5)


# Plot ldsc and magma results jointly as a bar plot 
ldsc_magma_level1 = level1_merged_df %>%
  filter(gwas %in% cardiovascular_traits_supp) %>%
  ggplot(aes(x=annotation, y=-log10(pvalue), fill=method)) + 
  geom_col(position = "dodge", width = 0.7) + 
  geom_hline(yintercept = 1.3, color="black", linetype="dashed") + 
  ggtitle(" LDSC and MAGMA level1 results") + 
  coord_flip() + 
  facet_wrap(~gwas, scales = c("free_x"), ncol = 4) + 
  custom_theme() + 
  theme(legend.position = "bottom",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size=12)) + 
  miller_discrete_scale(style = "bars", option = 2) 

ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/new_Supp_Figure3/SuppFig3_ldsc_magma_level1_res.pdf",
       plot = ldsc_magma_level1, width = 9, height = 4)



#################################################################
# LDSC results according to lesion status (level 1 annotations) #
#################################################################

# Plot results according to lesion status
lesion_level1_df = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/lesion_status/lesion/CELLECT-LDSC/results/prioritization.csv")
head(lesion_level1_df)

non_lesion_level1_df = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/lesion_status/non_lesion/CELLECT-LDSC/results/prioritization.csv")
head(non_lesion_level1_df)

# Refactor traits and cell annotations
lesion_dfs = c("lesion_level1_df", "non_lesion_level1_df")
for (i in lesion_dfs) {
  df = get(i)
  df$gwas = factor(df$gwas, 
                   levels = ordered_gwas_traits)
  df$annotation = factor(df$annotation, 
                         levels = ordered_level1_anno)
}

# Get p vals from each df
lesion_level1_df = lesion_level1_df %>%
  rename(lesion = pvalue) %>%
  select(gwas, annotation, lesion)
non_lesion_level1_df = non_lesion_level1_df %>%
  rename(non_lesion = pvalue) %>%
  select(gwas, annotation, non_lesion)

# Merge p values from each lesion category into a single df
lesion_status_merged_df = merge(lesion_level1_df, non_lesion_level1_df,
                                by.x = c("gwas", "annotation"),
                                by.y = c("gwas", "annotation"))

reshaped_df = reshape2::melt(lesion_status_merged_df, id = c("gwas", "annotation"))
reshaped_df = reshaped_df %>%
  rename(lesion_status = variable,
         p_value = value)
reshaped_df$gwas = factor(reshaped_df$gwas, 
                          levels = ordered_gwas_traits)
head(reshaped_df)

# Plot results
ldsc_level1_lesion_status = reshaped_df %>%
  filter(gwas %in% gwas_traits_of_interest) %>%
  ggplot(aes(x=annotation, y=-log10(p_value), fill=lesion_status)) + 
  geom_col(position = "dodge", width = 0.7) + 
  geom_hline(yintercept = 1.3, color="black", linetype="dashed") + 
  scale_size_manual(values=c(2,4)) + 
  ggtitle("LDSC results according to lesion status") + 
  coord_flip() + 
  facet_wrap(~gwas, ncol = 5) + 
  custom_theme() + 
  theme(legend.position = "bottom",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size=12)) + 
  miller_discrete_scale(style = "bars", option = 2)

ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/new_Supp_Figure3/new_SuppFig3b_ldsc_lesion_status.pdf",
       plot = ldsc_level1_lesion_status, width = 14, height = 6)



####################################################
# Plot results for SMCs, pericytes and fibroblasts #
####################################################

level2_ldsc_smc_df = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/SMC_Peri_Fibro_CELLECT_outputs/CELLECT-LDSC/results/prioritization.csv")
head(level2_ldsc_smc_df)

level2_magma_smc_df = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/SMC_Peri_Fibro_CELLECT_outputs/CELLECT-MAGMA/results/prioritization.csv")
head(level2_magma_smc_df)

# Refactor SMC annotations
smc_vec = c("level2_ldsc_smc_df", "level2_magma_smc_df")
for (i in smc_vec) {
  df = get(i)
  df$gwas = factor(df$gwas, levels = ordered_gwas_traits)
  df$annotation = factor(df$annotation, levels = ordered_level2_smc)
}

# Get p vals from each df
ldsc_smc_df = level2_ldsc_smc_df %>%
  rename(LDSC = pvalue) %>%
  select(gwas, annotation, LDSC)
magma_smc_df = level2_magma_smc_df %>%
  rename(MAGMA = pvalue) %>%
  select(gwas, annotation, MAGMA)
level2_smc_merged_df = merge(ldsc_smc_df, magma_smc_df,
                             by.x = c("gwas", "annotation"),
                             by.y = c("gwas", "annotation"))
level2_smc_merged_df = reshape2::melt(level2_smc_merged_df,
                                      id=c("gwas", "annotation"))
names(level2_smc_merged_df)[c(3, 4)] = c("method", "pvalue")
level2_smc_merged_df$gwas = factor(level2_smc_merged_df$gwas, 
                                   levels = ordered_gwas_traits)


# Plot results as a dot plot (deprecated in current version of manuscript)
smc_ldsc_plot = smc_peri_fibro_prioritization_df %>%
  filter(gwas %in% cardiovascular_traits2) %>%
  ggplot(aes(x=annotation, y=-log10(pvalue), color=annotation)) + 
  geom_point(aes(size=-log10(pvalue) >= 1.3)) +
  scale_size_manual(values = c(2, 4)) + 
  xlab("Cell type") +
  ylab("-log10(LDSC p-value)") + 
  geom_hline(yintercept = 1.3, linetype="dashed", color="grey57") + 
  custom_theme + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size=14),
        legend.position = "right") + 
  npg_scale + 
  guides(size=FALSE,
         color = guide_legend(override.aes = list(size = 4))) + 
  facet_wrap(~gwas, scales = "free_y", ncol = 2)

ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure4/Fig4c_SMC_subtypes_LDSC_long_new2.pdf",
        plot=smc_ldsc_plot, width = 7, height = 14)


smc_peri_fibro_prioritization_df = smc_peri_fibro_prioritization_df %>%
  filter(gwas %in% cardiovascular_traits2)

write.table(smc_peri_fibro_prioritization_df,
            "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_tables/LDSC_level2_SMCs_annotations_results.txt",
            row.names = FALSE, quote = FALSE)


# Plot ldsc and magma results jointly as a bar plot 
ldsc_magma_level2_smc = level2_smc_merged_df %>%
  filter(gwas %in% cardiovascular_traits_3) %>%
  ggplot(aes(x=annotation, y=-log10(pvalue), fill=method)) + 
  geom_col(position = "dodge", width = 0.7) + 
  geom_hline(yintercept = 1.3, color="black", linetype="dashed") + 
  ggtitle(" LDSC and MAGMA level2 SMC results") + 
  coord_flip() + 
  facet_wrap(~gwas, scales = c("free_x"), ncol = 3) + 
  custom_theme() + 
  theme(legend.position = "bottom",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size=12)) + 
  miller_discrete_scale(style = "bars", option = 2) 

ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure4/revision_figures/Fig4d_ldsc_magma_level2_SMC.pdf",
       plot = ldsc_magma_level2_smc, width = 12, height = 7)



#################################################
# Plot results for Myeloid level 2 annotations  #
#################################################

ldsc_myeloid_prioritization_df = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/Myeloid/CELLECT-LDSC/results/prioritization.csv")
head(ldsc_myeloid_prioritization_df)

magma_myeloid_prioritization_df = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/Myeloid/CELLECT-MAGMA/results/prioritization.csv")
head(magma_myeloid_prioritization_df)

# Re factor GWAS traits for both LDSC and MAGMA results
ldsc_magma_vec = c("ldsc_myeloid_prioritization_df", 
                   "magma_myeloid_prioritization_df")
for (i in ldsc_magma_vec) {
  df = get(i)
  df$gwas = factor(df$gwas,
                   levels = ordered_gwas_traits)
}

# Rename cols
ldsc_myeloid_df = ldsc_myeloid_prioritization_df %>%
  rename(LDSC = pvalue) %>%
  select(gwas, annotation, LDSC)

magma_myeloid_df = magma_myeloid_prioritization_df %>%
  rename(MAGMA = pvalue) %>%
  select(gwas, annotation, MAGMA)

myeloid_merged_df = merge(ldsc_myeloid_df, magma_myeloid_df,
                  by.x = c("gwas", "annotation"),
                  by.y = c("gwas", "annotation"))
myeloid_merged_df = reshape2::melt(myeloid_merged_df, 
                           id=c("gwas", "annotation"))
names(myeloid_merged_df)[c(3, 4)] = c("method", "pvalue")
myeloid_merged_df$gwas = factor(myeloid_merged_df$gwas, 
                                levels = ordered_gwas_traits)
  

# Plot ldsc and magma results jointly
myeloid_merged_df %>%
  filter(gwas %in% c("CAD_MVP_EUR_meta", "Alzheimer_disease_Jansen", "White_blood_cell_count")) %>%
  #filter(gwas %in% gwas_traits_of_interest) %>%
  ggplot(aes(x=annotation, y=-log10(pvalue), fill=method)) + 
  geom_col(position = "dodge", width = 0.7) + 
  geom_hline(yintercept = 1.3, color="black", linetype="dashed") + 
  scale_size_manual(values=c(2,4)) + 
  ggtitle("Myeloid LDSC and MAGMA results") + 
  coord_flip() + 
  facet_wrap(~gwas) + 
  custom_theme() + 
  theme(legend.position = "right",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size=14)) + 
  miller_discrete_scale(style = "bars", option = 2)
  




