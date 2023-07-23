library(Seurat)
library(scRNAutils)
library(tidyverse)
library(scRNAutils)
library(data.table)
library(ggsci)
library(ggrepel)


#########################################################
# Plot results from MAGMA-based effector genes analyses #
#########################################################

# Plot results from effector gene analyses for SMCs level1 annotations
effector_genes = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/whole_ref/CELLECT-GENES/results/effector_genes.csv")

eff_genes_level1_smc = effector_genes %>%
  filter(gwas %in% c("CAD_MVP_EUR_meta", "Myocardial_infarction")) %>% 
  filter(annotation %in% c("SMC")) %>%
  ggplot(aes(x=-log10(magma_gene_pval), y=esmu, 
             label=ifelse(esmu >= 0.7 & -log10(magma_gene_pval) >= 4.5, 
                          gene_symbol, ""))) + 
  geom_point() + 
  facet_wrap(~ gwas, scales = "free_x") + 
  ggtitle("Meta-analysis SMC effector genes for CAD and MI") + 
  geom_text_repel(color="darkred", fontface="bold.italic") + 
  ylab("Cell type expression specificity (Esmu)") + 
  xlab("-log10(MAGMA genes p val)") + 
  custom_theme() + 
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size=14))

ggsave(filename = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/new_Supp_Figure3/new_SuppFig3c_MAGMA_CAD_MI_level1_SMC_effector_genes.pdf",
       plot = eff_genes_level1_smc, width = 10, height = 5)

##################################################
# Plot results from effector gene analyses for ECs
eff_genes_level1_mac = effector_genes %>%
  filter(gwas %in% c("White_blood_cell_count", "Alzheimer_disease_Jansen")) %>%
  filter(annotation %in% c("Macrophage")) %>%
  ggplot(aes(x=-log10(magma_gene_pval), y=esmu, 
             label=ifelse(esmu >= 0.7 & -log10(magma_gene_pval) >= 4.5, 
                          gene_symbol, ""))) + 
  geom_point() + 
  facet_wrap(~ gwas, scales = "free_x") + 
  ggtitle("Meta-analysis Macrophage effector genes for WBC count and AD") + 
  geom_text_repel(color="darkred") + 
  xlab("-log10(Magma genes p val)") +
  ylab("Cell type expression specificity (Esmu)") +
  custom_theme() + 
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size=14))

ggsave(filename = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/new_Supp_Figure3/new_SuppFig3c_MAGMA_WBC_AD_level1_Mac_effector_genes.pdf",
       plot = eff_genes_level1_mac, width = 10, height = 5)
  
####################################################  
# Plot effector genes for SMC subtypes
smc_subtypes_eff_genes = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/SMC_Peri_Fibro_CELLECT_outputs/CELLECT-GENES/results/effector_genes.csv")

eff_genes_fc = smc_subtypes_eff_genes %>%
  filter(gwas %in% c("CAD_MVP_EUR_meta")) %>%
  filter(annotation %in% c("Fibrochondrocyte")) %>%
  ggplot(aes(x=-log10(magma_gene_pval), y=esmu, 
             label=ifelse(esmu >= 0.8 & -log10(magma_gene_pval) >= 4, 
                          gene_symbol, ""))) + 
  geom_point() + 
  facet_wrap(~ gwas, scales = "free_x") + 
  ggtitle("Meta-analysis FC effector genes for CAD") + 
  geom_text_repel(color="darkred", fontface = "bold.italic") + 
  xlab("-log10(MAGMA genes p val)") +
  ylab("Cell type expression specificity (Esmu)") +
  custom_theme() + 
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size=14))
  
ggsave(filename = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Rebuttal_Fig13/CAD_CRTAC1_effector_gene.pdf",
       plot = eff_genes_fc, width = 4, height = 4)
  
  











