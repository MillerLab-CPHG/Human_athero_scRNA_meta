library(biomaRt)
library(tidyverse)
library(data.table)
library(scRNAutils)

# Source our own scRNA analysis utils functions
#source("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/R_scripts/Utils/scRNA_processing_utils.R")

######################################################################################################################
# This script has code to convert a list of mice genes to human homologs (only one-to-one orthology). This           #
# also involves creating mice SMC gene modules to test for their enrichment in the meta-analyzed human scRNA data.   #  
######################################################################################################################


# For this, we will source our custom functions (convert_mouse_to_human)

# Load mouse mart
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

##########################################
# Load df with mouse meta-analysis markers
mouse_meta_markers = read.csv("/project/cphg-millerlab/Jose/mouse_scRNA_meta_analyses/SMC_lin_tracing_meta_analysis/Supplementary_Tables/IJ-00992.R1_Supplementary_Table_S2_cluster_markers.csv", header = TRUE)
de_col_names = c("Gene_symbol", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "cluster")
mouse_meta_markers = mouse_meta_markers %>%
  dplyr::select(all_of(de_col_names)) %>%
  dplyr::filter(avg_logFC > 0)
dim(mouse_meta_markers)
head(mouse_meta_markers, n=20)

##################
# Get SEM homologs
# Looks like genes are already ordered by Log2FC
sem_df = mouse_meta_markers %>%
  filter(cluster == "SEM")
dim(sem_df)
sem_markers = sem_df$Gene_symbol
length(sem_markers)

# Convert SEM markers to human nomenclature
sem_human_homologs_list = convert_mouse_to_human(sem_markers, mouse_biomart = mouse)
sem_human_homologs_df = sem_human_homologs_list[[1]]

# Add human homologs back to markers df for SEM cells
sem_df_filtered = add_human_homologs_to_mice_df(sem_df, sem_human_homologs_df)
head(sem_df_filtered)

# Get vector of human homologs for mice SEM markers
# There are 244 human homologs for SEM markers from mouce meta analysis of SMC lineage-tracing scRNA
sem_human_homologs_vec = sem_df_filtered$human_homolog
sem_human_homologs_vec_top50 = sem_human_homologs_vec[1:50]
sem_human_homologs_vec_top100 = sem_human_homologs_vec[1:100]
length(sem_human_homologs_vec)

# Prep table for export
sem_df_filtered_new = sem_df_filtered %>%
  dplyr::select(Gene_symbol, human_homolog, cluster) %>%
  head(n=50)



##################
# Get FC homologs
# Looks like genes are already ordered by Log2FC
fc_df = mouse_meta_markers %>%
  filter(cluster == "FC")
dim(fc_df)
fc_markers = fc_df$Gene_symbol
length(fc_markers)

# Convert FC markers to human nomenclature
fc_human_homologs_list = convert_mouse_to_human(fc_markers, mouse_biomart = mouse)
fc_human_homologs_df = fc_human_homologs_list[[1]]

# Add human homologs back to markers df for FC cells
fc_df_filtered = add_human_homologs_to_mice_df(fc_df, fc_human_homologs_df)
head(fc_df_filtered)
tail(fc_df_filtered)

# Get vector of human homologs for mice FC markers
# There are 324 human homologs for FC markers from mouce meta analysis of SMC lineage-tracing scRNA
fc_human_homologs_vec = fc_df_filtered$human_homolog
fc_human_homologs_vec_top50 = fc_human_homologs_vec[1:50]
fc_human_homologs_vec_top100 = fc_human_homologs_vec[1:100]
length(fc_human_homologs_vec)

# Prep table for export 
fc_df_filtered_new = fc_df_filtered %>%
  dplyr::select(Gene_symbol, human_homolog, cluster) %>%
  head(n=50)

######################################
# Get SMC1 (Contractile SMCs) homologs
# Looks like genes are already ordered by Log2FC
smc1_df = mouse_meta_markers %>%
  filter(cluster == "SMC1")
head(smc1_df)
dim(smc1_df)
smc1_markers = smc1_df$Gene_symbol
length(smc1_markers)

# Convert FC markers to human nomenclature
smc1_human_homologs_list = convert_mouse_to_human(smc1_markers, mouse_biomart = mouse)
smc1_human_homologs_df = smc1_human_homologs_list[[1]]

# Add human homologs back to markers df for SMC1 cells
smc1_df_filtered = add_human_homologs_to_mice_df(smc1_df, smc1_human_homologs_df)
head(smc1_df_filtered)

# Get vector of human homologs for mice FC markers
# There are 324 human homologs for FC markers from mouce meta analysis of SMC lineage-tracing scRNA
smc1_human_homologs_vec = smc1_df_filtered$human_homolog
smc1_human_homologs_vec_top50 = smc1_human_homologs_vec[1:50]
smc1_human_homologs_vec_top100 = smc1_human_homologs_vec[1:100]
length(smc1_human_homologs_vec)

# Prep table for export 
smc1_df_filtered_new = smc1_df_filtered %>%
  dplyr::select(Gene_symbol, human_homolog, cluster) %>%
  head(n=50)

###################################################
# Get Fibro1 (Non-SMC-derived fibroblasts) homologs
# Looks like genes are already ordered by Log2FC
fibro1_df = mouse_meta_markers %>%
  filter(cluster == "Fibro1")
head(fibro1_df)
dim(fibro1_df)
fibro1_markers = fibro1_df$Gene_symbol
length(fibro1_markers)

# Convert FC markers to human nomenclature
fibro1_human_homologs_list = convert_mouse_to_human(fibro1_markers, mouse_biomart = mouse)
fibro1_human_homologs_df = fibro1_human_homologs_list[[1]]

# Add human homologs back to markers df for Fibro1 cells
fibro1_df_filtered = add_human_homologs_to_mice_df(fibro1_df, fibro1_human_homologs_df)
head(fibro1_df_filtered)

# Get vector of human homologs for mice Fibro1 markers
# There are 620 human homologs for Fibro1 markers from mouce meta analysis of SMC lineage-tracing scRNA
fibro1_human_homologs_vec = fibro1_df_filtered$human_homolog
fibro1_human_homologs_vec_top50 = fibro1_human_homologs_vec[1:50]
fibro1_human_homologs_vec_top100 = fibro1_human_homologs_vec[1:100]
length(fibro1_human_homologs_vec)

# Prep table for export
fibro1_df_filtered_new = fibro1_df_filtered %>%
  dplyr::select(Gene_symbol, human_homolog, cluster) %>%
  head(n=50)

# Merge dfs and write .tsv file
df_list = list(sem_df_filtered_new,
               fc_df_filtered_new,
               smc1_df_filtered_new,
               fibro1_df_filtered_new)
merged_gene_sets = data.table::rbindlist(df_list)
dim(merged_gene_sets)

write.table(merged_gene_sets,
            "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_tables/meta_analyzed_murine_SMC_gene_sets.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


#######################################################
# Test for signature enrichments from Kim et al. 2020 #
#######################################################

# Load markers from Kim et al 
kim_smc_markers = fread("/project/cphg-millerlab/Jose/mouse_scRNA_data/Kim_et_al/SMC_markers.csv")

# Get Contractile SMC homologs 
# Looks like genes are already sorted by avg Log2FC (Only top50 are available)
smc_kim_df = kim_smc_markers %>%
  filter(cluster == "SMC")
head(smc_kim_df)
dim(smc_kim_df)

names(smc_kim_df)[8] = "Gene_symbol"

smc_kim_markers = smc_kim_df$Gene_symbol

# Convert SMC markers to human nomenclature
smc_kim_human_homologs_list = convert_mouse_to_human(smc_kim_markers, 
                                                     mouse_biomart=mouse)
smc_kim_homologs_df = smc_kim_human_homologs_list[[1]]

# Add human homologs back to markers df for SEM cells
# Kept 45 homologs
smc_kim_filtered = add_human_homologs_to_mice_df(smc_kim_df, smc_kim_homologs_df)
head(smc_kim_filtered)
dim(smc_kim_filtered)

# Get vector of human homologs for Kim et al mice contractile SMC markers
smc_kim_homologs_vec = smc_kim_filtered$human_homolog
smc_kim_homologs_vec

##############################################################################
# Get Fibromyocyte homologs 
# Looks like genes are already sorted by avg Log2FC (Only top50 are available)
fmc_kim_df = kim_smc_markers %>%
  filter(cluster == "FMC")
head(fmc_kim_df)
dim(fmc_kim_df)

names(fmc_kim_df)[8] = "Gene_symbol"

fmc_kim_markers = fmc_kim_df$Gene_symbol

# Convert SMC markers to human nomenclature
fmc_kim_human_homologs_list = convert_mouse_to_human(fmc_kim_markers, 
                                                     mouse_biomart=mouse)
fmc_kim_homologs_df = fmc_kim_human_homologs_list[[1]]

# Add human homologs back to markers df for SEM cells
fmc_kim_filtered = add_human_homologs_to_mice_df(fmc_kim_df, fmc_kim_homologs_df)
head(fmc_kim_filtered)
dim(fmc_kim_filtered)

# Get vector of human homologs for Kim et al mice contractile SMC markers
fmc_kim_homologs_vec = fmc_kim_filtered$human_homolog
fmc_kim_homologs_vec

################################################################################
# Get FC homologs 
# Looks like genes are already sorted by avg Log2FC (Only top50 are available)
fc_kim_df = kim_smc_markers %>%
  filter(cluster == "CMC")
head(fc_kim_df)
dim(fc_kim_df)

names(fc_kim_df)[8] = "Gene_symbol"

fc_kim_markers = fc_kim_df$Gene_symbol

# Convert FC markers to human nomenclature
fc_kim_human_homologs_list = convert_mouse_to_human(fc_kim_markers, 
                                                     mouse_biomart=mouse)
fc_kim_homologs_df = fc_kim_human_homologs_list[[1]]

# Add human homologs back to markers df for FCs
# Kept 45 homologs
fc_kim_filtered = add_human_homologs_to_mice_df(fc_kim_df, fc_kim_homologs_df)
head(fc_kim_filtered)
dim(fc_kim_filtered)

# Get vector of human homologs for Kim et al mice contractile SMC markers
fc_kim_homologs_vec = fc_kim_filtered$human_homolog
fc_kim_homologs_vec













