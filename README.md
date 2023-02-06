# Human_athero_scRNA_meta

This repository contains scripts used for generation of the meta-analyzed reference and other downstream analyses currently included in the manuscript. The repository can be broadly divided into 3 main sections: 1) QC and processing of each sequencing library 2) Integration benchmark, building of reference and annotations 3) Downstream analyses

All analyses were carried using R (v.4.0.3) and python (v.3.8)

### 0) Packages and libraries used

The following packages are required to run the scripts:

*Single cell data processing and integration:*

- Seurat (v.4.1.0)
- celda (v.1.6.1)
- scDblFinder (v.1.4.0)
- tidyverse (v.1.3.1)
- data.table (v.1.14.2)
- reticulate (v.1.18)
- RColorBrewer (v.1.1.2)
- ggsci (v.2.9)
- scanorama (v.1.7.1)
- lisi (v.1.0)
- cluster (v.2.1.0)
- SeuratDisk (v.0.0.0.9019)
- celltypist (v.1.0.0)

*LDSC workflow:*
- numpy (v.1.19.2)
- scipy (v.1.5.2)
- pandas (v.1.4.2)
- cellex (v.1.2.2)
- matplotlib (v.3.5.2)
- snakemake (v.6.0.5)

*Downstream analyses:*
- monocle3 (v.1.0.0)
- CellChat (v.1.5.0)
- UCell (v.1.3.1)
- dorothea (v.1.8.0)
- biomaRt (v.2.46.0)
- viper (v.1.24.0)
- gprofiler2 (v.0.2.1)
- ggridges (v.0.5.3)
- stats (v.4.0.3)



### 1) QC and processing of each individual sequencing library

Scripts used for QC and normalization of libraries can be found in the following dirs:

- `raw_counts_processing`

As detailed in the manuscript, the main criteria for processing of each of the libraries includes:

- removal of doublets 
- removal of ambient RNA
- Number of unique UMIs and unique genes expressed. % of reads mapped to the mitochondrial genome, Hb genes.
- dimensionality reduction


### 2) Benchmark of integration methods

Scripts used for the benchmark and generation of the integrated reference were included in the following dirs:

- `integration_benchmarking`

Main criteria/metrics calculated in the benchmark include:
- Running time
- Local Inverse Simpson Index (LISI)
- Silhouette scores calculated across 0.8-1.8 resolutions

### 3) Generation of harmonized scRNA reference 

Scripts used for building the reference were included in the following dirs:

- `ref_integration`
- `automated_annotations`

Level 1 annotations were done through Transfer Learning using the [Tabula Sapiens](https://tabula-sapiens-portal.ds.czbiohub.org/) cardiovascular subset as reference. T. Sapiens data was redownloaded and reprocessed using SCTransform normalization to match processing and normalization parameters of the built reference.

Level 2 annotations were done through a combination of automated annotation and manual curation. For the automated portion we used [celltypist](https://github.com/Teichlab/celltypist) with 4 different models as follows:

```
#!/bin/bash

# Annotate cell types based on Immune_All_AddPIP model (no majority voting)
celltypist --indata /project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/CELLECT_inputs/rpca_int_sct_v3_norm_counts.csv \
--model Immune_All_AddPIP.pkl \
--outdir /project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/celltypist_outputs/Immune_All_AddPIP \
--prefix celltypist_All_Immune_AddPIP_ \
--transpose-input \
--plot-results
```

### 4) Downstream analyses (including SMC modulation)

Donwstream analyses including annotation and characterization of SMC phenotypes are included in the following dirs:

- `SMC_analyses`
- `S-LDSC`

Characterization of SMC modulated phenotypes includes:

- Enrichment of murine SMC gene signatures
- Pseudotime analyses
- Transcription Factor (TF activity inference) using human [DoRothEA](https://saezlab.github.io/dorothea/) regulons

For LD score regression applied to specifically expressed genes (LDSC-SEG) analyses, we first generated a gene expression specificity matrix using SCTransformed-normalized gene expression for level 1 and 2 annotations using CELLEX. GWAS summary statistics were munged as follows: 

```
#!/bin/bash

# Run the mtag_munge script to make sure that GWAS summary stats are correctly formatted for LDSC
python /project/cphg-millerlab/Jose/CELLECT/ldsc/mtag_munge.py \
--sumstats /project/cphg-millerlab/Jose/GWAS_summary_stats/GWAS_pulse_pressure_MVP/pp_MVP_summary_for_ldsc.txt.gz \
--merge-alleles /project/cphg-millerlab/Jose/CELLECT/data/ldsc/w_hm3.snplist \
--out /project/cphg-millerlab/Jose/human_scRNA_meta_analysis/CELLECT_analyses/munged_GWAS_summary_stats/pulse_pressure_MVP_GWAS/munged_pulse_pressure_MVP \
--n-value 300000 \
--keep-pval \
--p PVAL \
--chunksize 500000
```
Prior integration with GWAS summary statistics, we defined a `cellect_config_multi_GWAS.yml` file with attributes pointing to munged GWAS summary statistics and the gene expression specificity matrix. We then ran the LDSC snakemake workflow as follows: 

```
#!/bin/bash

# Load modules. Looks like python3.8 is needed for the script to finish running but need to do a few more tests. 
module load anaconda/2020.11-py3.8
module load bzip2/1.0.6
module load snakemake/6.0.5

# Run ldsc analysis using CELLECT wrapper.
# Before running this command, it might be necessary to gunzip the gene expression specificity file. Somehow it seems to help make things run more smoothly. 
snakemake --use-conda -j -s /project/cphg-millerlab/Jose/CELLECT/cellect-ldsc.snakefile --configfile cellect_config_multi_GWAS.yml
