# --- Load Libraries -----------------------------------------------------------
library(Seurat)
library(SpatialDecon)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(tidyverse)
library(ggpubr)

# --- Load Data ----------------------------------------------------------------
# loads data from previous script "04_ligand_receptor.R"
load(paste0("./results/ligand-receptor.RData"))

# --- Prep for Deconvolution ---------------------------------------------------
# load the seurat object
matched_single_cell<-readRDS("./data/6_CellAssignment_1.rds")

matched_single_cell$patient <- gsub("-.*","",gsub(".*-h","",matched_single_cell$cells))

#Cell Type Assignments
matched_single_cell$cleaned_assignments <- str_replace_all(
  matched_single_cell$final_assignments,
  c(
    "\\bCM\\b" = "Cardiomyocyte",
    "\\bFB\\b" = "Fibroblast",
    "\\bEC\\b" = "Endothelial",
    "\\bPC\\b" = "Pericyte",
    "\\bSM\\b" = "Smooth Muscle",
    "\\bMAC\\b" = "Macrophage",
    "\\bLEC\\b" = "Lymphatic Endothelial",
    "\\bNEURO\\b" = "Neuron",
    "\\bDC\\b" = "Dendritic cell",
    "\\bTLYMPH\\b" = "T Lymphocyte"
  )
) 
# set the identities as the cleaned assignments
Idents(matched_single_cell) <- matched_single_cell$cleaned_assignments

#Noramlize Data
matched_single_cell <- NormalizeData(object = matched_single_cell)

# prepare nanostring counts and meta data
meta <- pData(target_combineData) 
q3_counts <- assayDataElement(target_combineData,
                              elt = "log_q")
q3_counts <- q3_counts[,rownames(meta)]
colnames(q3_counts) <- meta$id
raw_counts <- assayDataElement(target_combineData,
                               elt = "exprs")
raw_counts <- raw_counts[,rownames(meta)]
colnames(raw_counts) <- meta$id
rownames(meta)<- meta$id
meta$Phenotype <- gsub("NORMAL","Control",meta$Phenotype)
meta$composite<-  factor(paste(meta$Phenotype,meta$Segment.tags.pathologist_1,sep=": "),
                         levels = c(
                           "CONTROL: Normal",
                           "CONTROL: Mild",
                           "HCM: Normal",
                           "HCM: Mild",
                           "HCM: Moderate",
                           "HCM: Severe"
                         )
)
# --- UMAP plot ----------------------------------------------------------------
# umap plot of single cell data split by HCM status
DimPlot(
  matched_single_cell,
  split.by = "genotype",
  cols = ggpubr::get_palette("npg",k=10),
  ncol = 1)

# --- Investigate Available Marker Genes ---------------------------------------

# find all markers
sn_markers <- FindAllMarkers(matched_single_cell, 
                             only.pos = TRUE,
                             min.pct = 0.25, 
                             logfc.threshold = 0.25)
# filter for significant markers
sn_sig_markers <- sn_markers %>%
  filter(p_val_adj < 0.05 &
           abs(avg_log2FC) > 1.5)

# how many of these markers are available in the geomx data
# isolate geomx genes
geomx_genes <- rownames(q3_counts)

# create a list of marker genes by cell type
sn_sig_marker_genes <- split(sn_sig_markers$gene,sn_sig_markers$cluster)

# find out how many marker genes are in our geomx data
overlapping_markers <- lapply(
  sn_sig_marker_genes,
  function(x){
    return(x[x %in% geomx_genes])
  }
)

# create a dotplot overlapping marker genes in single cell data
DotPlot(matched_single_cell,features = unique(unlist(overlapping_markers)))+
  coord_flip()+
  RotatedAxis()+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))+
  labs(
    #title = "Expression of Overlapping Single-Cell/Geomx Cell Marker Genes in Single-Cell Data",
    x="",
    y=""
  )

# isolate geomx expression of these shared markers
geomx_marker_exp <- as.data.frame(q3_counts[rownames(q3_counts) %in%
                                unique(unlist(overlapping_markers)),]) %>%
  mutate(gene=rownames(.)) %>%
  reshape2::melt() %>%
  group_by(variable) %>%
  mutate(total_per_roi=sum(value)) %>%
  ungroup() %>%
  mutate(percent_exp = value/total_per_roi*100)

# create a dotplot overlapping marker genes in geomx data
ggplot(geomx_marker_exp,
       aes(
         x=variable,
         y=gene,
         size=percent_exp,
         fill=value
       ))+
  geom_point()+
  ggpubr::theme_pubr(legend = "right") +
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))+
  RotatedAxis()+
  scale_size(range = c(0.01,3))+
  labs(
    x="",
    y="",
    fill="Log2(Exp+1)",
    size="Percent of Total ROI Expression",
    title="Expression of Overlapping Single-Cell/Geomx Markers in Geomx Data"
  )

# --- Overlapping Marker Gene Deconvolution Per Patient ------------------------
# define a function to deconvolute data on a per paitent basis
# where each patient's GeoMx data is deconvoluted using the snRNA-seq 
# data from that particular patient
per_patient_decon_overlapping <- function(seurat,
                              patient_col,
                              id,
                              nanostring_counts,
                              raw_nanostring_counts,
                              overlapping_genes,
                              nanostring_meta){
  # grab log2 averaged expression values per patient
  avg <- AverageExpression(seurat[,seurat[[patient_col]]==id])$RNA
  avg <- log2(avg+1)
  avg <- avg[rownames(avg) %in% overlapping_genes,]
  # derive a background using the negative probes for one patient
  bg = derive_GeoMx_background(
    norm = nanostring_counts[,grepl(id,colnames(nanostring_counts))],
    probepool = rep(1, nrow(nanostring_counts[,grepl(id,colnames(nanostring_counts))])),
    negnames = "Negative Probe")
  # run deconvolution
  spd <- spatialdecon(norm = as.matrix(
    nanostring_counts[,grepl(
      id,colnames(nanostring_counts))]),
                      bg = bg,
                      X = avg,
                      align_genes = TRUE)
  #prepare a df for plotting
  df=reshape2::melt(spd$prop_of_all)
  colnames(df) <- c("Cell","id","value")
  df <- left_join(df,
                  nanostring_meta,
                  by=c("id"="id"))
  df$idSegment <- paste(df$id,df$Segment.tags.pathologist_1,sep = "_")
  return(df)
}

# apply deconvolution function to all patient ids
patients <- unique(meta$`slide name`)

# edit the spatial decon function to remove arbitrary cutoff:
# trace("spatialdecon",edit=T)
# change 100 shared genes to 50

# apply our per-patient deconvolution function
per_patient_deconed_overlapping <- as.data.frame(
  do.call(rbind,
          lapply(patients,
                 function(x){
                   tmp=per_patient_decon_overlapping(seurat = matched_single_cell,
                                         patient_col = "patient",
                                         id=x,
                                         nanostring_counts = q3_counts,
                                         overlapping_genes = unique(unlist(overlapping_markers)),
                                         raw_nanostring_counts = raw_counts,
                                         nanostring_meta = meta);
                   return(tmp)
                 })
  )) 

# cell diversity statistic described here:
# https://github.com/tanya-karagiannis/Cell-Type-Diversity-Statistic
# samples by row, cell types by column
per_patient_overlapping_deconed_mat <- reshape2::dcast(
  per_patient_deconed_overlapping, 
  id~Cell,
  value.var="value") %>%
  as.data.frame() %>%
  column_to_rownames("id")

div_res_per_patient_overlapping <- apply(
  per_patient_overlapping_deconed_mat,
  1,
  function(x){(-sum(x*log(x), na.rm = T)/log(ncol(per_patient_overlapping_deconed_mat))-1)})

# merge diversity statistic into the meta data
meta$per_patient_overlapping_diversity_statistic <- div_res_per_patient_overlapping


# --- Plotting Cell Types Over Phenotype/Disarray Levels -----------------------
# proceeding with Deconvolution Performed Using Averaged Per-Patient Single Cell Expression Data
# Filtered By Overlapping Cell Type Markers

final_avg_cell_prop_plot =  ggpubr::ggscatter(
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

# --- Final Deconvolution Bar Plot ---------------------------------------------
final_decon_bar = ggpubr::ggbarplot(
  data = per_patient_deconed_overlapping %>%
    filter(!is.na(value) & 
             composite != "CONTROL: Mild"	) %>%
    mutate(Phenotype=gsub("NORMAL","Control",Phenotype)),
  x="idSegment",
  y="value",
  color="Cell",
  fill = "Cell",
  xlab = "ROIs",
  ylab="Cell Proportions"
) +
  ggpubr::theme_pubr(legend = "top")+
  ggh4x::facet_grid2(~Phenotype+slide.name+Segment.tags.pathologist_1,
             drop=TRUE,
             scales='free_x',
             space = "free_x",
             strip = ggh4x::strip_nested(bleed = TRUE))+
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top")+
scale_fill_manual(values=colors,
                  name="Cell Types")+
scale_color_manual(values=colors,
                   name="Cell Types")

final_decon_bar
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
# --- Diversity Statisitc Mixed Model ------------------------------------------

emmeans_res= emmeans::emmeans(nlme::lme(per_patient_overlapping_diversity_statistic ~ composite,random=~1|slide.name, meta),pairwise~composite) 
emmeans_res_df = emmeans_res$contrasts %>% as.data.frame()
emmeans_res_df

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
# --- Save Data ----------------------------------------------------------------
save.image(
  file = paste0("./results/deconvolution.RData")
)
