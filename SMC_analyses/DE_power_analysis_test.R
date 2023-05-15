library(Seurat)
library(tidyverse)
library(data.table)
library(scPower)
library(cowplot)



# Load latest version of the reference 
rpca_int_sct_v3_1 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/seurat_objects/alsaigh_pan_wirka_hu_int_seurat_clustered_v3_1.rds")
rpca_smc_fibro_subset_v3 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/smc_pericytes_fibroblasts_subclustering_method1/seurat_objects/whole_ref_v3/peri_smc_fibro_subset_seurat_obj_res0.9_v2.rds")



# Extract df with UMI counts
raw_umi_counts = rpca_smc_fibro_subset_v3@assays$RNA@data
raw_umi_counts[1:4, 1:4]

# Look at the read depth (UMI counts for the first cell) and also the count distribution
sum(raw_umi_counts[, 1])
hist(raw_umi_counts[, 1])

# How many samples do we have? 22 libraries (sample size) 
rpca_smc_subset_meta = rpca_smc_fibro_subset_v3@meta.data
rpca_smc_fibro_subset_v3@meta.data %>%
  count(sample)

# What are the cell type frequencies in the subset
# Total number of cells in the subset = 37032 cells
cell_freqs = rpca_smc_subset_meta %>%
  group_by(prelim_annotations) %>%
  summarize(n=n(), 
            cell_freq = n/nrow(rpca_smc_subset_meta))
sum(cell_freqs$cell_freq)

# Calculate n cells per sample
cells_per_sample = rpca_smc_subset_meta %>%
  count(sample)


# Calculate the read depth for each cell in the reference
cells_umi_counts = apply(raw_umi_counts, MARGIN = 2, FUN = sum)
cells_umi_counts_df = as.data.frame(cells_umi_counts)
cells_umi_counts_df %>%
  ggplot(aes(x=cells_umi_counts)) + 
  geom_histogram() + 
  geom_freqpoly(color="darkred") +
  ggtitle("SMC read depth (UMI counts)") + 
  custom_theme

# Do a quick run of scPower
# We have a mean UMI cell count of 7120
smc_pwr = power.general.restrictedDoublets(nSamples = 22, nCells = 37000, readDepth = 7120,
                                     ct.freq = 0.193, ref.study = de.ref.study, 
                                     ref.study.name = "Blueprint (CLL) iCLL-mCLL", type = "de",
                                     read.umi.fit = read.umi.fit[read.umi.fit$type == "10X_PBMC_1",],
                                     gamma.mixed.fits = gamma.mixed.fits, ct="CD14+ Monocytes",
                                     cellsPerLane = 20000,
                                     disp.fun.param = disp.fun.param,
                                     min.UMI.counts = 200,
                                     multipletRateGrowth = "constant",
                                     multipletRate = 0.1)

################################################
# Run a test to define priors with our own data 
toy_matrix_list = scPower::count.matrix.example
length(toy_matrix_list)
sapply(toy_matrix_list, dim)

count_matrix1 = toy_matrix_list[[1]]
count_matrix1[1:4, 1:4]


# Count the number of expressed genes 
expressed_genes_df = NULL

# Iterate over each count matrix
for (name in names(toy_matrix_list)) { 
  count_matrix = toy_matrix_list[[name]]
  
  # Create an annotation file (here containing only one cell type, but can have more)
  annot_df = data.frame(individual=paste0("S", rep(1:14, length.out=ncol(count_matrix))),
                                          cell.type=rep("default_ct", ncol(count_matrix)))
  
  # Reformat count matrix into pseudobulk matrix
  pseudo_bulk = create.pseudobulk(count_matrix, annot_df)
  
  # Calculate expressed genes in the pseudobulk matrix
  expressed_genes = calculate.gene.counts(pseudo_bulk, min.counts = 3, perc.indiv.expr = 0.5)
  
  # Get the number of expressed genes
  num_expressed_genes = nrow(expressed_genes)
  
  # Save expressed genes
  expressed_genes_df = rbind(expressed_genes_df,
                             data.frame(matrix=name,
                                        num_cells=ncol(count_matrix),
                                        expressed_genes=num_expressed_genes))
  
}

print(expressed_genes_df)


# Estimate the negative binomial parameters for each gene
norm_mean_values =  NULL
disp_param = NULL
for (name in names(toy_matrix_list)) {
  temp = nbinom.estimation(toy_matrix_list[[name]])
  
  # Save normalized mean values
  norm_mean_values_temp = temp[[1]]
  norm_mean_values_temp$matrix = name
  norm_mean_values = rbind(norm_mean_values, norm_mean_values_temp)
  
  # Save the parameter of the mean-dspersion function
  disp_param_temp = temp[[3]]
  disp_param_temp$matrix = name
  disp_param = rbind(disp_param, disp_param_temp)
  
}

# First rows with normalized mean values
head(norm_mean_values)

# First rows with parameter of the mean-dispersion function
head(disp_param)


############################################################
# Estimation of a gamma mixed distribution over all means

gamma_fits = NULL

for (name in names(toy_matrix_list)) { 
  
  # Number of cells per cell types as censoring point
  censored_point = 1 / ncol(toy_matrix_list[[name]])
  
  norm_mean_values_temp = norm_mean_values[norm_mean_values$matrix==name,]
  gamma_fit_temp = mixed.gamma.estimation(norm_mean_values_temp$mean,
                                          num.genes.kept = 21000,
                                          censoredPoint = censored_point)
  gamma_fit_temp$matrix = name
  gamma_fits = rbind(gamma_fits, gamma_fit_temp)
  
}
print(gamma_fits)

# Comparison of gamma mixed fits with original means
g = visualize.gamma.fits(norm_mean_values$mean[norm_mean_values$matrix],
                         gamma_fits[gamma_fits$matrix=="subsampled25",],
                         nGenes = 21000)

print(g)


######################################################################################
# Parameterization of the parameters of the gamma fits by the mean UMI counts per cell

# Estimate the mean umi values per cell for each matrix
umi_values = NULL
for(name in names(toy_matrix_list)) {
  mean.umi = meanUMI.calculation(toy_matrix_list[[name]])
  umi_values = rbind(umi_values, data.frame(mean.umi, matrix=name))
}
print(umi_values)

gamma_fits = merge(gamma_fits, umi_values, by = "matrix")

# Convert the gamma fits from the shap rate parameterization to
# the mean-sd parmeterization 
gamma_fits = convert.gamma.parameters(gamma.fits = gamma_fits)

# Visualize the linear relationship between gamma parameters and UMI values in plots
plot_values = melt(gamma_fits, id.vars = c("matrix", "mean_umi"))
plot_values = plot_values[plot_values$variable %in% c("mean1", "mean2", "sd1", "sd2", "p1", "p2"), ]
ggplot(plot_values, aes(x=mean_umi, y=value)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~variable, ncol=2, scales="free")

# Fit relationship between gamma parameters and UMI values
gamma_linear_fit_new = umi.gamma.relation(gamma_fits)
print(gamma_linear_fit_new)


##############################################################
# Estimation of median dispersion function for each cell type
disp_fun_general_new = dispersion.function.estimation(disp_param)
print(disp_fun_general_new)


###############################################################
# Annotation of cell type for all fitted data frames
gamma_linear_fit_new$ct = "New_ct"
disp_fun_general_new$ct = "New_ct"

###############################################################
# Fitting a function for UMI counts dependent on read depth

# Number of mapped reads taken from cellranger
mapped_reads = data.frame(matrix=c("complete", "subsampled75", "subsampled50", "subsampled25"),
                          transcriptome.mapped.reads=c(23130, 17422, 11666, 5859))

# Plot relationship between mean reads per cell and mean UMI per cell
read_umis = merge(umi_values, mapped_reads, by="matrix")
print(read_umis)

ggplot(read_umis, aes(x=transcriptome.mapped.reads, y=mean.umi)) + 
  geom_point() + 
  geom_line()

# Fit relationship between mean reads per cell and mean UMI per cell
read_umi_fit_new = umi.read.relation(read_umis)
print(read_umi_fit_new)

###############################################################
# Validation of expression probability model 

# Merge the observed numbers of expressed genes with the read depth
expressed_genes_df = merge(expressed_genes_df, mapped_reads, by="matrix")

# Get the number of cells per cell type and individual 
expressed_genes_df$cells_indiv = expressed_genes_df$num_cells/14
expressed_genes_df$estimated.genes = NA

for (i in 1:nrow(expressed_genes_df)) { 
  
  # Vector with the expression probability for each gene
  expr.prob = estimate.exp.prob.param(nSamples = 14,
                                      readDepth = expressed_genes_df$transcriptome.mapped.reads[i], 
                                      nCellsCt = expressed_genes_df$cells_indiv[i],
                                      read.umi.fit = read_umi_fit_new,
                                      gamma.mixed.fits = gamma_linear_fit_new,
                                      ct="New_ct",
                                      disp.fun.param = disp_fun_general_new,
                                      min.counts = 3,
                                      perc.indiv = 0.5)
  
  # Expected number of expressed genes
  expressed_genes_df$estimated.genes[i] = round(sum(expr.prob))
  
  }
print(expressed_genes_df)

plot_expressed_genes_df = reshape2::melt(expressed_genes_df,
                                         id.vars=c("matrix", "num_cells", "cells_indiv",
                                                   "transcriptome.mapped.reads"))

ggplot(plot_expressed_genes_df, aes(x=transcriptome.mapped.reads, y=value,
                                    color=variable)) + 
  geom_point() + 
  geom_line()





























