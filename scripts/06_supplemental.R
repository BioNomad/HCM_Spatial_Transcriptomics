# --- Load Libraries -----------------------------------------------------------
LIB1='/cluster/tufts/hpc/tools/R/4.0.0'
LIB2='/cluster/home/jlaird01/R/x86_64-pc-linux-gnu-library/4.0'
.libPaths(c(LIB1,LIB2))
library(Seurat)
library(SpatialDecon)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(clusterProfiler)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(ComplexUpset)
library(openxlsx)

# --- Load Data ----------------------------------------------------------------
# loads data from previous script "05_deconvolution.R"
# if you did not run the previous script the same day
# swap out Sys.Date() with the date you ran the previous step
load(paste0("./results/deconvolution.RData"))

# --- All DEG Comparison UpsetR plot -------------------------------------------

# isolate genes with an FDR < 0.05
sig_list <- disarray_pheno_combined_results %>%
  filter(FDR<0.05)

# split into a list of degs per contrast level
sig_list <- split(sig_list$Gene,sig_list$Contrast)

# convert list into a boolean matrix
sig_gene_int_mat <- reshape2::melt(
  sig_list
) %>% reshape2::acast(L1~value) %>%
  t()
sig_gene_int_mat = ifelse(is.na(sig_gene_int_mat),
                          FALSE,
                          TRUE) %>%
  as.data.frame()

# extract all contrast levels
all_comparison_names <- names(sig_gene_int_mat)

# set column coloring 
set.cols <- data.frame(
  cols=c(rep('#c795a8',8), rep('#9b95c7',6)),
  set=all_comparison_names,
  Comparison=ifelse(grepl("CONTROL_Normal|CONTROL_Mild|CONTROL", all_comparison_names) &
                      grepl("- HCM", all_comparison_names),
                    "Between HCM Status",
                    "Between Disarray Levels")
)

# plot the overlap between deg groups
upset(
  sig_gene_int_mat,
  colnames(sig_gene_int_mat),
  name="",
  min_size = 2,
  stripes=upset_stripes(
    mapping=aes(color=Comparison,size=14),
    colors=c(
      "Between HCM Status"='#9b95c7',
      "Between Disarray Levels"='#c795a8'
    ),
    data=set.cols
  )
)

# --- DotPlot Of All Comparisons -----------------------------------------------
dotplot(
  object = disarray_pheno_combined_enrichment_cd$enriched %>%
    mutate(Contrast=gsub("_",":",Contrast)) %>%
    mutate(Cluster=gsub("_",":",Cluster)) ,
  showCategory=5,
  by='count'
)+
  theme_pubr(legend="right")+
  facet_grid(.~Contrast,
             space="free",
             scales="free",
             labeller = label_wrap_gen(width = 22))+
  rotate_x_text(angle = 45)+
  scale_color_gradient(low = "midnightblue",
                       high="thistle")+
  scale_x_discrete(labels=function(x)gsub(".*\\.","",x))+
  labs(
    color="Adjusted P-Value"
  )
# --- Significant DEG Table ----------------------------------------------------
# set header style for result tables
hs <- createStyle(fontColour = "#ffffff", fgFill = "#232436",
                  halign = "center", valign = "center", 
                  textDecoration = "Bold",
                  border = "TopBottomLeftRight",
                  borderColour = "#ffffff",
                  fontSize = 14)

# write deg table to a file 
write.xlsx(
  disarray_pheno_combined_results %>%
             filter(FDR<0.05) %>%
    rename("P-Value"=`Pr(>|t|)`) %>%
    rename("FDR Adjusted P-Value"="FDR") %>%
    rename("Change Direction" = "change_direction"), 
  file = "./results/phenotype_disarray_differential_expression_results.xlsx",
  headerStyle = hs,
  colWidths="auto")

# --- Significant LR Table -----------------------------------------------------

# write significant LR gene pairs to a file
write.xlsx(
  diff_comb_ints, 
  file = "./results/significant_ligand_receptor_pairs_results.xlsx",
  headerStyle = hs,
  colWidths="auto")

# --- DotPlot of Overlapping Cell Marker Genes ---------------------------------

# create a dotplot overlapping marker genes in single cell data
DotPlot(matched_single_cell,features = unique(unlist(overlapping_markers)))+
  coord_flip()+
  RotatedAxis()+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))+
  labs(
    title = "Expression of Overlapping Single-Cell/Geomx Cell Marker Genes in Single-Cell Data"
  )

# --- Overlapping Cell Marker Genes Used In Deconvolution ----------------------

# write cell marker genes to a file
write.xlsx(
  data.frame(
    lapply(
      overlapping_markers, 
      "length<-", 
      max(lengths(overlapping_markers))
      )
    ), 
  file = "./results/overlapping_cell_marker_genes.xlsx",
  headerStyle = hs,
  colWidths="auto")

# --- Deconvolution Broken Down By Patient -------------------------------------

avg_prop_plot = ggpubr::ggscatter(
  data = per_patient_deconed_overlapping %>%
    mutate(Cell=factor(Cell,
                       levels=c(
                         "T Lymphocyte",
                         "Neuron",
                         "Pericyte",
                         "Lymphatic Endothelial",
                         "Fibroblast",
                         "Smooth Muscle",
                         "Macrophage",
                         "Endothelial",
                         "Dendritic cell",
                         "Cardiomyocyte"
                       )))   %>%
    group_by(composite,slide.name,Cell) %>%
    summarise(avg_cell_prop=mean(value)*100) %>%
    filter(composite != "CONTROL: Mild"),
  x="composite",
  y="Cell",
  size="avg_cell_prop",
  color="composite",
  palette = colors[c(1,2,6,4,3,5)],
  xlab = "",
  ylab = "",
  title = ""
)+
  theme_pubr(legend = "right")+
  scale_size(range = c(0.1,6))+
  guides(color = FALSE)+
  labs(
    size="Average Cell Proportion (%)",
  )+
  rotate_x_text(angle = 55)+
  ggh4x::facet_grid2(~composite+slide.name,
                     drop=TRUE,
                     scales='free_x',
                     space = "free_x",
                     #labeller = label_wrap_gen(width = 10),
                     strip = ggh4x::strip_nested(bleed = TRUE))


div_plot = ggpubr::ggscatter(
  data = meta %>%
    filter(composite != "CONTROL: Mild"),
  x="composite",
  y="per_patient_overlapping_diversity_statistic",
  color="composite",
  palette = colors[c(1,2,6,4,3,5)]
)+
  theme_pubr(legend = "right")+
  rotate_x_text(angle = 55)+
  ggh4x::facet_grid2(~composite+slide.name,
                     drop=TRUE,
                     scales='free_x',
                     space = "free_x",
                     strip = ggh4x::strip_nested(bleed = TRUE))+
  labs(
    x = "",
    y = "Cell Type Diversity Statistic",
    color = ""
  )+
avg_prop_plot/div_plot
