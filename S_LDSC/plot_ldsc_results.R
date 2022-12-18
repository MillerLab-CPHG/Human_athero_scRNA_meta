library(tidyverse)
library(data.table)
library(ggsci)

######################################################################
# The goal of this script is to plot results from the LDSC analysis  #
######################################################################

prioritization_df = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/CELLECT-LDSC/results/prioritization.csv")
head(prioritization_df)

# Add multiple testing correction p-values
prioritization_df$FDR = p.adjust(prioritization_df$pvalue, method = "fdr")

prioritization_df$gwas = factor(prioritization_df$gwas,
                                levels = c("CAD_GWAS_van_der_Harst", "CAD_GWAS_Koyama", "Myocardial_infarction", "Carotid_plaque", "Carotid_IMT",
                                           "Pulse_pressure_MVP", "Diastolic_blood_pressure_MVP", "Systolic_blood_pressure_MVP", "Alzheimer_disease_Jansen",
                                           "White_blood_cell_count", "Body_mass_index", "Type_2_diabetes"))

gwas_traits = c("Alzheimer_disease_Jansen", "Body_mass_index",  "Carotid_plaque",
                "Pulse_pressure_MVP", "Diastolic_blood_pressure_MVP", "Systolic_blood_pressure_MVP",
                "Type_2_diabetes", "White_blood_cell_count")

# Refactor level 1 annotations
prioritization_df$annotation = factor(prioritization_df$annotation,
                                      levels = c("B_cell", "Endothelial", "Fibroblast", "SMC",
                                                 "Mast_cell", "Neuron", "Pericyte", "Macrophage",
                                                 "Plasma_cell", "pDC", "T_NK"))

level1_ldsc = prioritization_df %>%
  filter(gwas %in% c("CAD_GWAS_van_der_Harst", "Myocardial_infarction", "Pulse_pressure_MVP", "Carotid_plaque",
                    "Alzheimer_disease_Jansen", "Body_mass_index", "Type_2_diabetes", "White_blood_cell_count")) %>%
  #filter(gwas %in% c(gwas_traits)) %>%
  ggplot(aes(x=annotation, y=-log10(pvalue), color=annotation)) + 
  #geom_point(aes(size=ifelse(-log10(pvalue) >= 1.3, 1.5, 1))) + 
  geom_point(aes(size=-log10(pvalue) >= 1.3)) + 
  #scale_size_manual(values=c(2.7, 6)) + 
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
  npg_scale2 + 
  guides(size=FALSE,
         color = guide_legend(override.aes = list(size=4))) + 
  facet_wrap(~gwas, scales = "free_y", nrow = 2)


ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Figure1/Fig1e_level1_ldsc_new.pdf",
       plot=level1_ldsc, width = 12, height = 5)

prioritization_df = prioritization_df %>%
  filter(gwas %in% c("CAD_GWAS_van_der_Harst", "Myocardial_infarction", "Pulse_pressure_MVP", "Carotid_plaque",
                     "Alzheimer_disease_Jansen", "Body_mass_index", "Type_2_diabetes", "White_blood_cell_count"))
write.table(prioritization_df, 
            "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_tables/LDSC_level1_annotations_results.txt",
            row.names = FALSE, quote = FALSE)


# Plot results for SMCs, pericytes and fibroblasts
smc_peri_fibro_prioritization_df = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/SMC_Peri_Fibro_CELLECT_outputs/CELLECT-LDSC/results/prioritization.csv")
head(smc_peri_fibro_prioritization_df)

smc_peri_fibro_prioritization_df$gwas = factor(smc_peri_fibro_prioritization_df$gwas,
                                                  levels = c("CAD_GWAS_van_der_Harst", "CAD_GWAS_Koyama", "Myocardial_infarction", "CAC_GWAS", "Carotid_plaque", "Carotid_IMT",
                                                  "Pulse_pressure_MVP", "Diastolic_blood_pressure_MVP", "Systolic_blood_pressure_MVP", "Alzheimer_disease_Jansen",
                                                  "White_blood_cell_count", "Type_2_diabetes", "Body_mass_index"))

cardiovascular_traits = c("CAD_GWAS_van_der_Harst", "Myocardial_infarction", "CAC_GWAS", "Pulse_pressure_MVP", "Diastolic_blood_pressure_MVP",
                          "Carotid_plaque", "Alzheimer_disease_Jansen" ,"Type_2_diabetes" , "White_blood_cell_count" ,"Systolic_blood_pressure_MVP")

cardiovascular_traits2 = c("CAD_GWAS_van_der_Harst", "Myocardial_infarction", "CAC_GWAS", "Pulse_pressure_MVP",
                           "Carotid_plaque", "Type_2_diabetes")

# Refactor SMC annotations

smc_peri_fibro_prioritization_df %>%
  mutate(annotation2 = case_when(annotation == "Unknown" ~ "SMC3",
                                TRUE ~ annotation))

idx = which(smc_peri_fibro_prioritization_df$annotation == "Unknown")
smc_annotations = smc_peri_fibro_prioritization_df$annotation
smc_annotations[idx] = "SMC3"

smc_peri_fibro_prioritization_df$annotation = smc_annotations
smc_peri_fibro_prioritization_df$annotation = factor(smc_peri_fibro_prioritization_df$annotation,
                                                               levels = c("Contractile_SMC", "Transitional-ECM-SMC", "Fibromyocyte",
                                                                          "Fibrochondrocyte", "Fibroblast", "SMC2", "SMC3", "Foam-like", 
                                                                          "Pericyte1", "Pericyte2"))

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


