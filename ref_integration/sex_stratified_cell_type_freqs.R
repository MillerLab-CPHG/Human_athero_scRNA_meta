library(Seurat)
library(tidyverse)
library(scRNAutils)
library(ggsci)
library(ggpubr)


# Load meta-analyzed data
rpca_int_sct_v3_1 = read_rds("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/rds_objects/integration_rds_objects/rPCA/alsaigh_pan_wirka_hu_int/seurat_objects/alsaigh_pan_wirka_hu_int_seurat_clustered_v3_1.rds")
DimPlot(rpca_int_sct_v3_1, group.by = "updated_level2_annotations", raster = FALSE) + 
  theme(legend.position = "none")

# Plot XIST expression by library
xist_by_library = VlnPlot(rpca_int_sct_v3_1, features = c("XIST"), group.by = "sample") + 
  theme(legend.position = "none")
xist_by_library


# Plot cell type frequencies by level 1 annotation
meta_df = rpca_int_sct_v3_1@meta.data

# How many cells do we have per sex?
# We have 19648 cells from females and 98930 from males 
table(meta_df$sex)

# Make stacked bar plot
meta_df %>%
  ggplot(aes(x=level1_annotations, fill=sex)) + 
  geom_bar(position = "fill") + 
  custom_theme() + 
  miller_discrete_scale(style = "bars") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###########################################################################
# Make bar plot (normalizing for total number of cells in each sex category)
level1_freq_by_sex = meta_df %>%
  group_by(level1_annotations, sex) %>%
  summarize(n=length(level1_annotations)) %>%
  mutate(percentage = case_when(sex=="males" ~ n/98930 * 100,
                                TRUE ~ n/19648 * 100))

level1_freq_by_sex_plot = level1_freq_by_sex %>%
  ggplot(aes(x=level1_annotations, y=percentage, fill=sex)) + 
  geom_col(position = "dodge", width = 0.5) + 
  ylab("Cell type percentage") + 
  xlab("Level 1 annotations") + 
  custom_theme() + 
  ggtitle("Cell type frequency normalized by total number of cells by sex") + 
  miller_discrete_scale(style = "bars") +
  coord_flip()

ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Level1_anno_by_sex_bar_plot.pdf",
       plot = level1_freq_by_sex_plot, width = 8, height = 8)


#########################################################################################################
# Make bar plot for level 2 Endo annotations (normalizing for total number of cells in the Endo category)
endo_meta = meta_df %>%
  filter(level1_annotations %in% c("Endothelial")) 

# How many SMCs from males and how many from females?
# We have 11292 cells from males and 1387 from females
table(endo_meta$sex)

endo_freq_by_sex = endo_meta %>%
  group_by(updated_level2_annotations, sex) %>%
  summarize(n=length(updated_level2_annotations)) %>%
  mutate(percentage = case_when(sex=="males" ~ n/11292 * 100,
                                TRUE ~ n/1387 * 100))

endo_freq_by_sex_plot = endo_freq_by_sex %>%
  ggplot(aes(x=updated_level2_annotations, y=percentage, fill=sex)) + 
  geom_col(position = "dodge", width = 0.5) + 
  ylab("Cell type percentage") + 
  xlab("Level 2 Endo annotations") + 
  custom_theme() + 
  ggtitle("Endothelial cell type frequency normalized by total number of cells by sex") + 
  miller_discrete_scale(style = "bars") + 
  coord_flip()

ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Level2_EC_by_sex_bar_plot.pdf",
       plot = endo_freq_by_sex_plot, width = 8, height = 8)

#######################################################################################################
# Make bar plot for level 2 Mac annotations (normalizing for total number of cells in the Mac category)
mac_meta = meta_df %>%
  filter(level1_annotations %in% c("Macrophage"))

# How many Macs from males and females
# We have 17966 Macs from males and 4948 from females 
table(mac_meta$sex)

mac_freq_by_sex = mac_meta %>%
  group_by(updated_level2_annotations, sex) %>%
  summarize(n=length(updated_level2_annotations)) %>%
  mutate(percentage=case_when(sex=="males" ~ n/17966 * 100,
                              TRUE ~ n/4948 * 100))

mac_freq_by_sex_plot = mac_freq_by_sex %>%
  ggplot(aes(x=updated_level2_annotations, y=percentage, fill=sex)) + 
  geom_col(position="dodge", width=0.5) + 
  ylab("Cell type percentage") + 
  xlab("Level 2 Myeloid annotations") + 
  custom_theme() + 
  ggtitle("Macs cell type frequency normalized by total number of cells by sex") + 
  miller_discrete_scale(style = "bars") + 
  coord_flip()

ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/sex_tratification_figures/Level2_Mac_by_sex_bar_plot.pdf",
       plot = mac_freq_by_sex_plot, width = 8, height = 8)

#########################################################################################################
# Make bar plot for level 2 SMC annotations (normalizing for total number of mural cells)
smc_annotations = c("Contractile_SMC", "Transitional-ECM-SMC", "SMC2", "SMC3", "Fibromyocyte", "Fibrochondrocyte")
smc_meta = meta_df %>%
  filter(updated_level2_annotations %in% smc_annotations)

# How many SMCs from males and females
# We have 21622 cells from males and 4767 from females
table(smc_meta$sex)

smc_freq_by_sex = smc_meta %>%
  group_by(updated_level2_annotations, sex) %>%
  summarize(n=length(updated_level2_annotations)) %>%
  mutate(percentage=case_when(sex=="males" ~ n/21622 * 100,
                              TRUE ~ n/4767 * 100))

smc_freq_by_sex_plot = smc_freq_by_sex %>%
  ggplot(aes(x=updated_level2_annotations, y=percentage, fill=sex)) + 
  geom_col(position = "dodge", width = 0.5) + 
  ylab("Cell type percentage") + 
  xlab("Level 2 SMC annotations") + 
  custom_theme() + 
  ggtitle("SMC cell type frequency normalized by total number of cells by sex") + 
  miller_discrete_scale(style = "bars") + 
  coord_flip()
  
ggsave("/project/cphg-millerlab/Jose/human_scRNA_meta_analysis/manuscript_figures/Level2_SMC_by_sex_bar_plot.pdf",
       plot = smc_freq_by_sex_plot, width = 8, height = 8)
  
  
  






