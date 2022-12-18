library(Seurat)
library(tidyverse)
library(pheatmap)




# This script will run gene set enrichments for level 1 annotations


# Load list of markers for level 1 annotations
level1_markers = readRDS("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/cluster_markers/rpca_int_sct_v3/cluster_markers_level1_annotations_logFC_ordered.rds")

# Load SEM object ro remove GO terms redundancy if needed. 
#hsGO = GOSemSim::godata("org.Hs.eg.db", ont="BP")



############################################
# Macrophage cells
mac_markers = level1_markers %>%
  filter(cluster == "Macrophage")
mac_ranks = mac_markers$avg_log2FC
names(mac_ranks) = mac_markers$gene
length(mac_ranks)

mac_ranks_top100 = mac_ranks[1:100]
mac_ranks_top50 = mac_ranks[1:50]

# Run gene set enrichment
mac_gost_res = gprofiler2::gost(names(mac_ranks_top100), organism = "hsapiens",
                                         ordered_query = TRUE, correction_method = "fdr")
mac_gprofiler_res = mac_gost_res$result
mac_go_bp = mac_gprofiler_res[mac_gprofiler_res$source == "GO:BP", ]

# Plot GO BP terms
mac_go_plot = gprofiler_bar_plot(mac_go_bp, 15, 500, "GO:BP", 10) + new_scale3_bars


##############################################
# T/NK

t_nk_markers = level1_markers %>%
  filter(cluster == "T_NK")
t_nk_ranks = t_nk_markers$avg_log2FC
names(t_nk_ranks) = t_nk_markers$gene
length(t_nk_ranks)

t_nk_ranks_top100 = t_nk_ranks[1:100]
t_nk_ranks_top50 = t_nk_ranks[1:50]

# Run gene set enrichment
t_nk_gost_res = gprofiler2::gost(names(t_nk_ranks_top100), organism = "hsapiens",
                                         ordered_query = TRUE, correction_method = "fdr")
t_nk_gprofiler_res = t_nk_gost_res$result
t_nk_go_bp = t_nk_gprofiler_res[t_nk_gprofiler_res$source == "GO:BP", ]

# Plot GO BP terms
t_nk_go_plot = gprofiler_bar_plot(t_nk_go_bp, 15, 500, "GO:BP", 10) + new_scale3_bars


############################################
# Endothelial cells
endo_markers = level1_markers %>%
  filter(cluster == "Endothelial")
endo_ranks = endo_markers$avg_log2FC
names(endo_ranks) = endo_markers$gene
length(endo_ranks)

endo_ranks_top100 = endo_ranks[1:100]
endo_ranks_top50 = endo_ranks[1:50]

# Run gene set enrichment
endo_gost_res = gprofiler2::gost(names(endo_ranks_top100), organism = "hsapiens",
                                         ordered_query = TRUE, correction_method = "fdr")
endo_gprofiler_res = endo_gost_res$result
endo_go_bp = endo_gprofiler_res[endo_gprofiler_res$source == "GO:BP", ]

# Plot GO BP terms
endo_go_plot = gprofiler_bar_plot(endo_go_bp, 15, 500, "GO:BP", 10) + new_scale3_bars

#################################################
# SMCs
smc_markers = level1_markers %>%
  filter(cluster == "SMC")
smc_ranks = smc_markers$avg_log2FC
names(smc_ranks) = smc_markers$gene
length(smc_ranks)

smc_ranks_top100 = smc_ranks[1:100]
smc_ranks_top50 = smc_ranks[1:50]

# Run gene set enrichment
smc_gost_res = gprofiler2::gost(names(smc_ranks_top100), organism = "hsapiens",
                                         ordered_query = TRUE, correction_method = "fdr")
smc_gprofiler_res = smc_gost_res$result
smc_go_bp = smc_gprofiler_res[smc_gprofiler_res$source == "GO:BP", ]

# Plot GO BP terms
smc_go_plot = gprofiler_bar_plot(smc_go_bp, 15, 500, "GO:BP", 10) + new_scale3_bars


################################################################
# B cells
b_cell_markers = level1_markers %>%
  filter(cluster == "B_cell")
b_cell_ranks = b_cell_markers$avg_log2FC
names(b_cell_ranks) = b_cell_markers$gene
length(b_cell_ranks)

b_cell_ranks_top100 = b_cell_ranks[1:100]
b_cell_ranks_top50 = b_cell_ranks[1:50]

# Run gene set enrichment
b_cell_gost_res = gprofiler2::gost(names(b_cell_ranks_top100), organism = "hsapiens",
                                         ordered_query = TRUE, correction_method = "fdr")
b_cell_gprofiler_res = b_cell_gost_res$result
b_cell_go_bp = b_cell_gprofiler_res[b_cell_gprofiler_res$source == "GO:BP", ]

# Plot GO BP terms
b_cell_go_plot = gprofiler_bar_plot(b_cell_go_bp, 15, 500, "GO:BP", 15) + new_scale3_bars

#################################################################
# Plasma cells
plasma_cell_markers = level1_markers %>%
     filter(cluster == "Plasma_cell") 
plasma_cell_ranks = plasma_cell_markers$avg_log2FC
names(plasma_cell_ranks) = plasma_cell_markers$gene
length(plasma_cell_ranks)
 
plasma_cell_ranks_top100 = plasma_cell_ranks[1:100]
plasma_cell_ranks_top50 = plasma_cell_ranks[1:50]

# Run gene set enrichment
plasma_cell_gost_res = gprofiler2::gost(names(plasma_cell_ranks_top100), organism = "hsapiens",
                                       ordered_query = TRUE, correction_method = "fdr")
plasma_cell_gprofiler_res = plasma_cell_gost_res$result
plasma_cell_go_bp = plasma_cell_gprofiler_res[plasma_cell_gprofiler_res$source == "GO:BP", ]

# Plot GO BP terms
plasma_cell_go_plot = gprofiler_bar_plot(plasma_cell_go_bp, 15, 500, "GO:BP", 10) + new_scale3_bars

###############################################################
# Fibroblasts
fib_markers = level1_markers %>%
        filter(cluster == "Fibroblast") 
fib_ranks = fib_markers$avg_log2FC
names(fib_ranks) = fib_markers$gene
length(fib_ranks)

fib_ranks_top100 = fib_ranks[1:100]

# Run gene set enrichment
fib_gost_res = gprofiler2::gost(names(fib_ranks_top100), organism = "hsapiens",
                                ordered_query = TRUE, correction_method = "fdr")
fib_gprofiler_res = fib_gost_res$result
fib_go_bp = fib_gprofiler_res[fib_gprofiler_res$source == "GO:BP", ]

# Plot GO BP terms
fib_go_plot = gprofiler_bar_plot(fib_go_bp, 15, 500, "GO:BP", 10) + new_scale3_bars

##################################################################
# Pericytes
peri_markers = level1_markers %>%
           filter(cluster == "Pericyte") 
peri_ranks = peri_markers$avg_log2FC
names(peri_ranks) = peri_markers$gene


peri_ranks_top100 = peri_ranks[1:100]

# Run gene set enrichment
peri_gost_res = gprofiler2::gost(names(peri_ranks_top100), organism = "hsapiens",
                                 ordered_query = TRUE, correction_method = "fdr")
peri_gprofiler_res = peri_gost_res$result
peri_go_bp = peri_gprofiler_res[peri_gprofiler_res$source == "GO:BP", ]

# Plot GO BP terms
peri_go_plot = gprofiler_bar_plot(peri_go_bp, 15, 500, "GO:BP", 10) + new_scale3_bars

####################################################################
# Mast cells
mast_markers = level1_markers %>%
              filter(cluster == "Mast_cell") 
mast_ranks = mast_markers$avg_log2FC
names(mast_ranks) = mast_markers$gene
length(mast_ranks)

mast_ranks_top100 = mast_ranks[1:100]
mast_ranks_top50 = mast_ranks[1:50]

# Run gene set enrichment
mast_gost_res = gprofiler2::gost(names(mast_ranks_top50), organism = "hsapiens",
                                ordered_query = TRUE, correction_method = "fdr")
mast_gprofiler_res = mast_gost_res$result
mast_go_bp = mast_gprofiler_res[mast_gprofiler_res$source == "GO:BP", ]
 
# Plot GO BP terms
mast_go_plot = gprofiler_bar_plot(mast_go_bp, 15, 500, "GO:BP", 10) + new_scale3_bars


######################################################################
# pDC
pDC_markers = level1_markers %>%
                filter(cluster == "pDC") 
pDC_ranks = pDC_markers$avg_log2FC
names(pDC_ranks) = pDC_markers$gene
length(pDC_ranks)

 
pDC_ranks_top100 = pDC_ranks[1:100]
pDC_ranks_top50 = pDC_ranks[1:50]

# Run gene set enrichment
pDC_gost_res = gprofiler2::gost(names(pDC_ranks_top50), organism = "hsapiens",
                                ordered_query = TRUE, correction_method = "fdr")
pDC_gprofiler_res = pDC_gost_res$result
pDC_go_bp = pDC_gprofiler_res[pDC_gprofiler_res$source == "GO:BP", ]
  
# Plot GO BP terms
pDC_go_plot = gprofiler_bar_plot(pDC_go_bp, 15, 500, "GO:BP", 10) + new_scale3_bars

#################################################################
# Try a heatmap for level 1 annotation GO terms
prep_terms = function(go_df, annotation, nterms) { 
  filtered_terms = go_df %>%
    filter(term_size > 15 & term_size < 500) %>%
    mutate(log10_pval = -log10(p_value)) %>%
    mutate(annotation = annotation) %>%
    arrange(desc(log10_pval)) %>%
    head(n=nterms)
  return(filtered_terms)
}


mac_df = prep_terms(mac_go_bp, "Macrophage", n=7)
t_nk_df = prep_terms(t_nk_go_bp, "T_NK", n=7)
endo_df = prep_terms(endo_go_bp, "Endothelial", n=7)
smc_df = prep_terms(smc_go_bp, "SMC", n=7)
peri_df = prep_terms(peri_go_bp, "Pericyte", n=7)
b_cell_df = prep_terms(b_cell_go_bp, "B_cell", n=7)
plasma_df = prep_terms(plasma_cell_go_bp, "Plasma/B_cell", n=7)
mast_df = prep_terms(mast_go_bp, "Mast_cell", n=7)
pDC_df = prep_terms(pDC_go_bp, "pDC", n=7)
fib_df = prep_terms(fib_go_bp, "Fibroblasts", n=7)

merged_df = rbind(mac_df, 
                  t_nk_df, 
                  endo_df, 
                  smc_df,
                  peri_df,
                  plasma_df,
                  mast_df,
                  pDC_df,
                  fib_df)
new_df = merged_df %>%
  dplyr::select(annotation, term_name, log10_pval, term_id, source)
write.table(new_df, "/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_tables/GSEA_level1_annotations_top30_GO_BP_terms.tsv",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


reshaped_df = spread(new_df, key = annotation, value = log10_pval)
reshaped_df[is.na(reshaped_df)] = 0
go_matrix = reshaped_df %>%
  dplyr::select(-term_name, -term_id, -source) %>% 
  as.matrix()
rownames(go_matrix) = reshaped_df$term_name

# Plot heatmap
palette_length = 100
my_color3 = colorRampPalette(c("white", brewer.pal(n=9, name = "Blues")))(palette_length)

level1_annotations_heatmap = pheatmap(go_matrix, color = my_color3, fontsize = 12, angle_col = 45)
ggsave(file="/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Supplementary_Figure2/SuppFig2a_level1_GSEA_heatmap.pdf",
       plot = level1_annotations_heatmap, width = 10, height = 10.5)





