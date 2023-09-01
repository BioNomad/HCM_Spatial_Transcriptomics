# --- Load Libraries -----------------------------------------------------------
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
load(paste0("./results/deconvolution.RData"))


# --- DEG Violin Plots by Patient ----------------------------------------------
merged = meta %>%
  inner_join(.,
             q3_counts %>% 
               t() %>%
               as.data.frame() %>%
               mutate(id=rownames(.)),
             by="id")

comparisons_of_interest = c(
  "CONTROL: Normal - HCM: Normal",
  "HCM: Normal - HCM: Severe"
)
top_deg_plots <- list()
for(comparison in comparisons_of_interest) {
  
  colors = get_palette("npg",k=10)[1:length(unique(merged$slide.name))]
  names(colors) = unique(merged$slide.name)
  
  top_deg = disarray_pheno_combined_results %>% 
    mutate(Contrast = gsub("_",": ",Contrast)) %>%
    filter(Contrast == comparison) %>%
    arrange(FDR) %>% slice_head(n=1) %>%
    dplyr::select(Gene) %>%
    unlist()
  
  plot = ggviolin(
    data=merged %>%
      filter(grepl(
        gsub(" - ","|",comparison),
        composite
      )),
    x="slide.name",
    y=top_deg,
    color="slide.name",
    fill = "slide.name",
    palette = get_palette("npg",k=10)
  )+
    facet_grid(.~composite,
               scales = "free",
               space = "free")+
    theme_pubr(legend="right")+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    scale_fill_manual(
      values=colors
    )+
    scale_color_manual(
      values=colors
    )+
    labs(
      x="Patient",
      y=top_deg,
      fill="Patient",
      color="Patient"
    )
  top_deg_plots[[comparison]] = plot
}
do.call(wrap_plots,top_deg_plots)
  
# --- DEG Scatter Plots by Patient ---------------------------------------------
merged = meta %>%
  inner_join(.,
             q3_counts %>% 
               t() %>%
               as.data.frame() %>%
               mutate(id=rownames(.)),
             by="id")

comparisons_of_interest = c(
  "CONTROL: Normal - HCM: Normal",
  "HCM: Normal - HCM: Severe"
)
top_deg_plots <- list()
for(comparison in comparisons_of_interest) {
  
  colors = get_palette("npg",k=10)[1:length(unique(merged$slide.name))]
  names(colors) = unique(merged$slide.name)
  
  top_deg = disarray_pheno_combined_results %>% 
    mutate(Contrast = gsub("_",": ",Contrast)) %>%
    filter(Contrast == comparison) %>%
    arrange(FDR) %>% slice_head(n=1) %>%
    dplyr::select(Gene) %>%
    unlist()
  
  plot = ggscatter(
    data=merged %>%
      filter(grepl(
        gsub(" - ","|",comparison),
        composite
      )),
    x="id",
    y=top_deg,
    color="slide.name",
    palette = get_palette("npg",k=10)
  )+
    facet_grid(.~composite,
               scales = "free",
               space = "free")+
    theme_pubr(legend="right")+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    scale_color_manual(
      values=colors
    )+
    labs(
      x="ROI",
      y=top_deg,
      color="Patient"
    )
  top_deg_plots[[comparison]] = plot
}
do.call(wrap_plots,top_deg_plots)

# --- ZBTB Investigation -------------------------------------------------------
zbtb_genes = c(
  "ZBTB16",
  "ZBTB20",
  "HDAC5",
  "NCOR1", 
  "SIN3A"
)

zbtb_plots = list()
for(gene in zbtb_genes){
  plot=ggviolin(
    data=merged %>%
      filter(!(composite == "CONTROL: Mild")),
    x="slide.name",
    y=gene,
    color="slide.name",
    fill = "slide.name",
    palette = get_palette("npg",k=10)
  )+
    facet_grid(.~composite,
               scales = "free",
               space = "free",
               labeller = label_wrap_gen(width = 10))+
    theme_pubr(legend="right")+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    scale_fill_manual(
      values=colors
    )+
    scale_color_manual(
      values=colors
    )+
    labs(
      x="Patient",
      y=gene,
      fill="Patient",
      color="Patient"
    )
  zbtb_plots[[gene]] = plot
}
zbtb_plots$ZBTB16/zbtb_plots$ZBTB20/zbtb_plots$HDAC5/zbtb_plots$NCOR1/zbtb_plots$SIN3A+plot_layout(guides = "collect") & theme(legend.position = "bottom")

# --- Dotplot ------------------------------------------------------------------
obj =disarray_pheno_combined_enrichment_cd$enriched %>%
  filter(qvalue<0.1) %>%
  filter(grepl("CONTROL_Normal - HCM_Normal|HCM_Normal - HCM_Severe",Cluster)) %>%
  filter(!is.na(Cluster)) 
obj@compareClusterResult$p.adjust=obj@compareClusterResult$qvalue
obj@compareClusterResult$Contrast=gsub("_",": ",obj@compareClusterResult$Contrast)
dotplot(
 object = obj,
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

# --- Volcano Plot -------------------------------------------------------------

volcanoPlot(results = disarray_pheno_combined_results %>%
              mutate(Contrast=gsub("_",":",Contrast)) %>%
              filter(Contrast %in% c(
                "CONTROL:Normal - HCM:Normal",
                "HCM:Normal - HCM:Severe"
              )),
            title = "")+
  theme_pubr(legend = "bottom")

# --- Module Plots -------------------------------------------------------------

colors <- paletteer::paletteer_d("colorBlindness::paletteMartin")
colors <- colors[-1]

# atp module expression
atp_row <- obj %>%
  filter(Contrast == "CONTROL_Normal - HCM_Normal" & 
           Description == "ATP synthesis coupled electron transport") %>%
  slice_head(n=1)
atp_geneset <- strsplit(atp_row$geneID,"/")[[1]]

moduleGenePlot(
  target_combineData[,pData(target_combineData)$pheno_disarray != "CONTROL_Mild"],
  atp_geneset,
  meta_levels = c("CONTROL:Normal","CONTROL:Mild","HCM:Normal","HCM:Mild","HCM:Moderate","HCM:Severe"),
  meta_data_col="pheno_disarray",
  title="ATP synthesis coupled electron transport",
  legend_title="",
  contrast_col="",
  assay="log_q"
)

# interferon regulation module expression
interfereon_row <- obj %>%
  filter(Contrast == "CONTROL_Normal - HCM_Normal" & 
           Description == "regulation of type I interferon production") %>%
  slice_head(n=1)
interferon_geneset <- strsplit(interfereon_row$geneID,"/")[[1]]

moduleGenePlot(
  target_combineData[,pData(target_combineData)$pheno_disarray != "CONTROL_Mild"],
  interferon_geneset,
  meta_levels = c("CONTROL:Normal","CONTROL:Mild","HCM:Normal","HCM:Mild","HCM:Moderate","HCM:Severe"),
  meta_data_col="pheno_disarray",
  title="Regulation of type I interferon production",
  legend_title="",
  contrast_col="",
  assay="log_q"
)

# 	extracellular matrix module expression
ext_cell_mat_row <- obj %>%
  filter(Contrast == "HCM: Normal - HCM: Severe" & 
           Description == "extracellular matrix") %>%
  slice_head(n=1)
ext_cell_mat_geneset <- strsplit(ext_cell_mat_row$geneID,"/")[[1]]

moduleGenePlot(
  target_combineData[,pData(target_combineData)$pheno_disarray != "CONTROL_Mild"],
  ext_cell_mat_geneset,
  meta_levels = c("CONTROL:Normal","CONTROL:Mild","HCM:Normal","HCM:Mild","HCM:Moderate","HCM:Severe"),
  meta_data_col="pheno_disarray",
  title="Extracellular Matrix",
  legend_title="",
  contrast_col="",
  assay="log_q"
)

# --- Gene Expression of Genes of Interest -------------------------------------
# https://www.ahajournals.org/doi/10.1161/JAHA.119.015473
mazzarotto_hcm_genes=openxlsx::read.xlsx("./data/mazzarotto_hcm_genes.xlsx",colNames = T) %>%
  filter(!is.na(CLINGEN) | CLINGEN != "No")
hcm_genes_present = mazzarotto_hcm_genes$GENE[mazzarotto_hcm_genes$GENE %in% rownames(target_combineData)]
assayDataElement(target_combineData[hcm_genes_present,],
                 elt = "log_q") %>%
  t() %>% 
  as.data.frame() %>%
  mutate(slide_id=rownames(.)) %>% 
  pivot_longer(!slide_id) %>%
  inner_join(.,
             pData(target_combineData) %>%
               as.data.frame() %>%
               mutate(slide_id=rownames(.)),
             by="slide_id") %>%
  filter(composite != "CONTROL: Mild") %>%
  ggviolin(
    data=.,
    x="composite",
    y="value",
    fill="composite",
    color="composite",
    palette = colors[c(9,8,6,4,3,5)]
  )+
  theme_pubr(legend = "right")+
  rotate_x_text(angle = 55)+
  facet_grid(.~name)+
  labs(
    x="",
    y="",
    fill="",
    color=""
  )
  
# --- Avg proportions ----------------------------------------------------------

patient_2828 = per_patient_deconed_overlapping %>%
  filter(`slide name` == "2828")

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
                     labeller = label_wrap_gen(width = 10),
                     strip = ggh4x::strip_nested(bleed = TRUE))

# --- Diversity Statistic ------------------------------------------------------
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
  #stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = -.2) +
  ggh4x::facet_grid2(~composite+slide.name,
                     drop=TRUE,
                     scales='free_x',
                     space = "free_x",
                     strip = ggh4x::strip_nested(bleed = TRUE))+
  labs(
    x = "",
    y = "Cell Type Diversity Statistic",
    color = ""
  )
# --- Combined Supplemental ----------------------------------------------------

div_plot/avg_prop_plot

# --- Decon Stat try -----------------------------------------------------------
  
mod_lst =list()
for (i in unique(per_patient_deconed_overlapping$Cell)){
  a=per_patient_deconed_overlapping %>%
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
    filter(Cell == i) 
  # %>%
  #   mutate(value=100*value+1)
  b=emmeans::emmeans(nlme::lme(value ~ composite,random=~1|slide.name, a),pairwise~composite)
  mod_lst[[i]] = b$contrasts %>% as.data.frame()
}

c=do.call(rbind, mod_lst)

# --- Diversity Statisitc Mixed Model ------------------------------------------

b= emmeans::emmeans(nlme::lme(per_patient_overlapping_diversity_statistic ~ composite,random=~1|slide.name, meta),pairwise~composite) 
c = b$contrasts %>% as.data.frame()
c

# define our model
model = nlme::lme(per_patient_overlapping_diversity_statistic ~ composite,random=~1|slide.name, meta)


# model qc, code taken from
# https://dfzljdn9uc3pi.cloudfront.net/2020/9522/1/MixedModelDiagnostics.html
# Check if Residuals are normally distributed
qqnorm(resid(model))
qqline(resid(model))

plot(model)

#Calculate leverage
lev<-hat(model.matrix(model,data = meta))

#Plot leverage against standardized residuals
plot(resid(model,type="pearson")~lev,
     las=1,
     ylab="Standardised residuals",
     xlab="Leverage")

#Inspect the random effect distribution
ggdensity(data=nlme::ranef(model, condVar = TRUE),
          x="(Intercept)",
          fill = "midnightblue")+
  labs(
    title="random effect distribution"
  )+
  annotate(
    "text",
    x = .04, y=7, 
    label =paste("Shapiro-Wilks normality test:",
                 as.character(
                   round(shapiro.test(
                     as.numeric(
                       unlist(
                         nlme::ranef(model, condVar = TRUE)
                         )))$p.value,digits = 5) )))
