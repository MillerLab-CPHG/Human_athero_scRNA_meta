library(Seurat)
library(scRNAutils)
library(tidyverse)
library(data.table)


#################################################################################
# The goal of this script is to prep the magma gene sets for munging with scDRS #
#################################################################################

# Load hg19 gene coords
# Load gene coordinates directly from CELLECT files
gene_coords_hg19 = fread("/project/cphg-millerlab/Jose/CELLECT/data/shared/gene_coordinates.GRCh37.ensembl_v91.txt")
names(gene_coords_hg19) = c("Ensembl_id", "Chrom", "Start", "End", 
                            "Strand", "Gene_symbol")
# This file contains 19430 genes 
dim(gene_coords_hg19)
head(gene_coords_hg19)

#########################################################
# Prep CAD GWAS MVP EUR meta MAGMA GS for scDRS munging #
#########################################################

# Load MAGMA GS for CAD MVP EUR meta sum stats
cad_mvp_magma_gs = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/CELLECT-MAGMA/precomputation/CAD_MVP_EUR_meta/CAD_MVP_EUR_meta.genes.out")
head(cad_mvp_magma_gs)

# We have 18776 genes x 9 vars
dim(cad_mvp_magma_gs)

# We simply need to match the ensembl ids on the MAGMA gs to the ensembl ids in the gene coords file used for SNP mapping 
cad_mvp_magma_gs$GENE_SYMBOL = gene_coords_hg19$Gene_symbol[match(cad_mvp_magma_gs$GENE, 
                                                            gene_coords_hg19$Ensembl_id)]

# Format MAGMA GS as pval file for munging in scDRS
cad_mvp_scdrs_pval_file = cad_mvp_magma_gs %>%
  select(GENE_SYMBOL, P) %>%
  group_by(GENE_SYMBOL) %>%
  filter(P == min(P)) %>%
  distinct()
names(cad_mvp_scdrs_pval_file) = c("GENE", "CAD_MVP_EUR_meta")
write.table(cad_mvp_scdrs_pval_file,
            file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_pilot/magma_gene_sets/CAD_MVP_EUR_meta/CAD_GWAS_MVP_EUR_meta_pval.tsv",
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

# Format MAGMA GS as zstat file for munging in scDRS
cad_mvp_scdrs_zscore_file = cad_mvp_magma_gs %>%
  select(GENE_SYMBOL, ZSTAT) %>%
  group_by(GENE_SYMBOL) %>%
  filter(ZSTAT == max(ZSTAT)) %>%
  distinct()
names(cad_mvp_scdrs_zscore_file) = c("GENE", "CAD_MVP_EUR_meta")
write.table(cad_mvp_scdrs_zscore_file,
            file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_pilot/magma_gene_sets/CAD_MVP_EUR_meta/CAD_GWAS_MVP_EUR_meta_zscores.tsv",
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

##########################################################
# Prep White Blood Cell count MAGMA GS for scDRS munging #
##########################################################

# Load MAGMA GS from White Blood Cell sum stats 
wbc_count_magma_gs = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/CELLECT-MAGMA/precomputation/White_blood_cell_count/White_blood_cell_count.genes.out")
head(wbc_count_magma_gs)

# We have 18876 genes x 9 vars
dim(wbc_count_magma_gs)

# We simply need to match the ensembl ids on the MAGMA gs to the ensembl ids in the gene coords file used for SNP mapping 
wbc_count_magma_gs$GENE_SYMBOL = gene_coords_hg19$Gene_symbol[match(wbc_count_magma_gs$GENE, 
                                                              gene_coords_hg19$Ensembl_id)]

# Format MAGMA GS as pval file for munging in scDRS
wbc_count_scdrs_pval_file = wbc_count_magma_gs %>%
  select(GENE_SYMBOL, P) %>%
  group_by(GENE_SYMBOL) %>%
  filter(P == min(P)) %>%
  distinct()
names(wbc_count_scdrs_pval_file) = c("GENE", "White_Blood_Cell_count")
write.table(wbc_count_scdrs_pval_file,
            file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_pilot/magma_gene_sets/White_Blood_Cell_count/WBC_count_pval.tsv",
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

# Format MAGMA GS as zscores file for munging in scDRS 
wbc_count_scdrs_zcores_file = wbc_count_magma_gs %>%
  select(GENE_SYMBOL, ZSTAT) %>%
  group_by(GENE_SYMBOL) %>%
  filter(ZSTAT == max(ZSTAT)) %>%
  distinct()
names(wbc_count_scdrs_zcores_file) = c("GENE", "White_Blood_Cell_count")
write.table(wbc_count_scdrs_zcores_file,
            file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_pilot/magma_gene_sets/White_Blood_Cell_count/WBC_count_zscores.tsv",
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

################################################
# Prep Alzheimer MAGMA GS for munging in scDRS #
################################################

# Load MAGMA GS from White Blood Cell sum stats 
alzheimer_magma_gs = fread("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/CELLECT_outputs/whole_ref/CELLECT-MAGMA/precomputation/Alzheimer_disease_Jansen/Alzheimer_disease_Jansen.genes.out")
head(alzheimer_magma_gs)

# We have 18876 genes x 9 vars 
dim(alzheimer_magma_gs)

# We simply need to match the ensembl ids on the MAGMA gs to the ensembl ids in the gene coords file used for SNP mapping 
alzheimer_magma_gs$GENE_SYMBOL = gene_coords_hg19$Gene_symbol[match(alzheimer_magma_gs$GENE, 
                                                                    gene_coords_hg19$Ensembl_id)]
# Format MAGMA GS as pval file for munging in scDRS
alzheimer_scdrs_pval_file = alzheimer_magma_gs %>%
  select(GENE_SYMBOL, P) %>%
  group_by(GENE_SYMBOL) %>%
  filter(P == min(P)) %>%
  distinct()
names(alzheimer_scdrs_pval_file) = c("GENE", "Alzheimer_disease")

write.table(alzheimer_scdrs_pval_file,
            file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_pilot/magma_gene_sets/Alzheimer_disease/Alzheimer_disease_MAGMA_pval.tsv",
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

# Format MAGMA GS as zscores file for munging in scDRS
alzheimer_scdrs_zscores_file = alzheimer_magma_gs %>%
  select(GENE_SYMBOL, ZSTAT) %>%
  group_by(GENE_SYMBOL) %>%
  filter(ZSTAT == max(ZSTAT)) %>%
  distinct()
names(alzheimer_scdrs_zscores_file) = c("GENE", "Alzheimer_disease")

write.table(alzheimer_scdrs_zscores_file,
            file = "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/scDRS_analyses/scDRS_pilot/magma_gene_sets/Alzheimer_disease/Alzheimer_disease_MAGMA_zscores.tsv",
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)








