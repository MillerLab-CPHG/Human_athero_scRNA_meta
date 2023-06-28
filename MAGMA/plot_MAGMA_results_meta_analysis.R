library(Seurat)
library(tidyverse)
library(scRNAutils)
library(data.table)
library(ggsci)
library(ggrepel)

##################################################################
# Load prioritization results from MAGMA for level 1 annotations #
##################################################################

level1_prioritization_df = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/whole_ref/CELLECT-MAGMA/results/prioritization.csv")
head(level1_prioritization_df)

# Reorder GWAS traits
level1_prioritization_df$gwas = factor(level1_prioritization_df$gwas,
                                levels = c("CAD_MVP_EUR_meta", "CAD_GWAS_van_der_Harst", "CAC_GWAS", "Myocardial_infarction", "Carotid_plaque", "Carotid_IMT",
                                "Pulse_pressure_MVP", "Diastolic_blood_pressure_MVP", "Systolic_blood_pressure_MVP", "Alzheimer_disease_Jansen",
                                "White_blood_cell_count", "Body_mass_index", "Type_2_diabetes"))

gwas_traits_of_interest = c("CAD_MVP_EUR_meta", "CAD_GWAS_van_der_Harst", "CAC_GWAS", "Myocardial_infarction", "Carotid_plaque", "Pulse_pressure_MVP",
                            "Alzheimer_disease_Jansen", "White_blood_cell_count", "Body_mass_index", "Type_2_diabetes")

# Refactor level 1 annotations
level1_prioritization_df$annotation = factor(level1_prioritization_df$annotation,
                                             levels = c("B_cell", "Endothelial", "Fibroblast", "SMC",
                                                        "Mast_cell", "Neuron", "Pericyte", "Macrophage",
                                                        "Plasma_cell", "pDC", "T_NK"))


# Plot results 
magma_level1_results = level1_prioritization_df %>%
  filter(gwas %in% gwas_traits_of_interest) %>%
  ggplot(aes(x=annotation, y=-log10(pvalue), color=annotation,
         label= ifelse(-log10(pvalue) >= 1.3, annotation, ""))) + 
  geom_point(aes(size=-log10(pvalue) >= 1.3)) + 
  #geom_text_repel(color="black") + 
  geom_hline(yintercept = 1.3, color="darkred", linetype="dashed") + 
  scale_size_manual(values=c(2,4)) + 
  ggtitle("MAGMA results from level 1 annotations") + 
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
  
ggsave(filename = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/MAGMA_figures/MAGMA_results_level1_annotations.PDF",
       plot = magma_level1_results, width = 15, height = 7)


################################################################## 
# Load prioritization results from MAGMA for level 2 annotations #
##################################################################

smc_prioritization_df = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/SMC_Peri_Fibro_CELLECT_outputs/CELLECT-MAGMA/results/prioritization.csv")
head(smc_prioritization_df)

# Reorder GWAS traits
smc_prioritization_df$gwas = factor(smc_prioritization_df$gwas,
                                    levels = c("CAD_MVP_EUR_meta", "CAD_GWAS_van_der_Harst", "CAC_GWAS", "Myocardial_infarction", "Carotid_plaque", "Carotid_IMT",
                                                  "Pulse_pressure_MVP", "Diastolic_blood_pressure_MVP", "Systolic_blood_pressure_MVP", "Alzheimer_disease_Jansen",
                                                  "White_blood_cell_count", "Body_mass_index", "Type_2_diabetes"))

gwas_traits_of_interest2 = c("CAD_MVP_EUR_meta", "CAD_GWAS_van_der_Harst", "CAC_GWAS", "Myocardial_infarction", "Carotid_plaque", "Pulse_pressure_MVP",
                             "Type_2_diabetes")

# Refactor level 1 annotations
smc_prioritization_df$annotation = factor(smc_prioritization_df$annotation,
                                                     levels = c("Contractile_SMC", "Transitional-SMC", "Fibromyocyte",
                                                                "Fibrochondrocyte", "Fibroblast", "SMC2", "SMC3", "Foam-like", 
                                                                "Pericyte1", "Pericyte2"))

# Plot results 
magma_level2_smc_plot = smc_prioritization_df %>%
  filter(gwas %in% gwas_traits_of_interest2) %>%
  ggplot(aes(x=annotation, y=-log10(pvalue), color=annotation,
             label= ifelse(-log10(pvalue) >= 1.3, annotation, ""))) + 
  geom_point(aes(size=-log10(pvalue) >= 1.3)) + 
  #geom_text_repel(color="black") + 
  geom_hline(yintercept = 1.3, color="darkred", linetype="dashed") + 
  scale_size_manual(values=c(2,4)) + 
  ggtitle("MAGMA results from level 2 SMC annotations") + 
  facet_wrap(~gwas, scales = "free_y", ncol = 3) + 
  custom_theme() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size=14)) + 
  miller_discrete_scale() +
  guides(size=FALSE,
         color = guide_legend(override.aes = list(size=4)))

ggsave(filename = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/MAGMA_figures/MAGMA_results_smc_level2_annotations.pdf",
       plot = magma_level2_smc_plot, width = 14, height = 7)


# Plot results from effector gene analyses for SMCs level1 annotations
effector_genes = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/whole_ref/CELLECT-GENES/results/effector_genes.csv")

eff_genes_level1_smc = effector_genes %>%
  filter(gwas %in% c("CAD_MVP_EUR_meta", "Myocardial_infarction")) %>% 
  filter(annotation %in% c("SMC")) %>%
  ggplot(aes(x=-log10(magma_gene_pval), y=esmu, label=ifelse(esmu >= 0.7 & -log10(magma_gene_pval) >= 4.5, gene_symbol, ""))) + 
  geom_point() + 
  facet_wrap(~ gwas, scales = "free_x") + 
  ggtitle("Meta-analysis SMC effector genes for CAD and MI") + 
  geom_text_repel(color="darkred", fontface="bold.italic") + 
  ylab("Cell type expression specificity (Esmu)") + 
  xlab("-log10(MAGMA genes p val)") + 
  custom_theme() + 
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size=14))

ggsave(filename = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/MAGMA_figures/MAGMA_CAD_MI_level1_SMC_effector_genes.pdf",
       plot = eff_genes_level1_smc, width = 15, height = 8)

# Plot results from effector gene analyses for ECs
eff_genes_level1_endo = effector_genes %>%
  filter(gwas %in% c("Carotid_plaque")) %>%
  filter(annotation %in% c("Endothelial")) %>%
  ggplot(aes(x=-log10(magma_gene_pval), y=esmu, label=gene_symbol)) + 
  geom_point() + 
  facet_wrap(~ gwas, scales = "free_x") + 
  ggtitle("Meta-analysis Endothelial effector genes for Carotid plaque") + 
  geom_text_repel(color="darkgrey") + 
  xlab("-log10(Magma genes p val)") +
  ylab("Cell type expression specificity (Esmu)") +
  custom_theme() + 
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size=14))

ggsave(filename = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/MAGMA_figures/MAGMA_Carotid_plaque_level1_EC_effector_genes.pdf",
       plot = eff_genes_level1_endo, width = 8, height = 8)
  
  
# Plot effector genes for SMC subtypes
smc_subtypes_eff_genes = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/SMC_Peri_Fibro_CELLECT_outputs/CELLECT-GENES/results/effector_genes.csv")

eff_genes_fibromyo = smc_subtypes_eff_genes %>%
  filter(gwas %in% c("CAD_MVP_EUR_meta")) %>%
  filter(annotation %in% c("Fibrochondrocyte")) %>%
  ggplot(aes(x=-log10(magma_gene_pval), y=esmu, label=ifelse(esmu >= 0.8 & -log10(magma_gene_pval) >= 4, gene_symbol, ""))) + 
  geom_point() + 
  facet_wrap(~ gwas, scales = "free_x") + 
  ggtitle("Meta-analysis Fibromyocyte effector genes for CAD, MI and Pulse pressure") + 
  geom_text_repel(color="darkred", fontface = "bold.italic") + 
  xlab("-log10(MAGMA genes p val)") +
  ylab("Cell type expression specificity (Esmu)") +
  custom_theme() + 
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size=14))
  
ggsave(filename = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/MAGMA_figures/MAGMA_CAD_MI_Pulse_pressure_level2_Fibromyo_effector_genes.pdf",
       plot = eff_genes_fibromyo, width = 16, height = 8)
  
  











