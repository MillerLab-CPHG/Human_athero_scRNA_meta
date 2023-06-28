library(Seurat)
library(tidyverse)
library(scRNAutils)
library(data.table)
library(ggsci)
library(ggrepel)


# This script will plot results from MAGMA analyses

#########################################
# Load MAGMA test run using window=100kb
parent_dir = c("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_test/MAGMA_outs/")
magma_100kb = fread(paste(parent_dir, "/window_100kb/CELLECT_magma_test_outs/CELLECT-MAGMA/results/prioritization.csv", sep = ""))
head(magma_100kb)

# Load effector genes
magma_eff_genes_100kb = fread(paste(parent_dir, "window_100kb/CELLECT_magma_test_outs/CELLECT-GENES/results/effector_genes.csv", sep = ""))
head(magma_eff_genes_100kb)
magma_eff_genes_100kb_smc_fibromyo = magma_eff_genes_100kb %>%
  filter(annotation %in% c("SMC", "Fibromyocyte"))

# Load MAGMA test run using window=10kb
magma_10kb = fread(paste(parent_dir, "/window_10kb/CELLECT_magma_test_outs/CELLECT-MAGMA/results/prioritization.csv", sep = ""))
head(magma_10kb)

# Load MAGMA test run using window=5kb
magma_5kb = fread(paste(parent_dir, "/window_5kb/CELLECT_magma_test_outs/CELLECT-MAGMA/results/prioritization.csv", sep = ""))
head(magma_5kb)


##########################################################################
# Visualize results 
magma_5kb %>%
  ggplot(aes(x=annotation, y=-log10(pvalue), color=annotation, 
             label= ifelse(-log10(pvalue) >= 1.3, annotation, ""))) + 
  geom_point(aes(size=-log10(pvalue) >= 1.3)) + 
  geom_hline(yintercept = 1.3, color="darkred", linetype="dashed") + 
  geom_text_repel(color="black") + 
  ggtitle("MAGMA window 5kb") + 
  custom_theme() + 
  scale_size_manual(values=c(2, 4)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


# Visualize effector genes
magma_eff_genes_100kb_smc_fibromyo %>%
  ggplot(aes(x=-log10(magma_gene_pval), y=esmu, label=gene_symbol)) + 
  geom_point() + 
  geom_text_repel() + 
  ggtitle("MAGMA effector genes window 100kb") + 
  ylab("Cell type expression specificity (Esmu)") + 
  facet_wrap(~annotation) + 
  custom_theme()










