# --- Load Libraries -----------------------------------------------------------
LIB1='/cluster/tufts/hpc/tools/R/4.0.0'
LIB2='/cluster/home/jlaird01/R/x86_64-pc-linux-gnu-library/4.0'
.libPaths(c(LIB1,LIB2))
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(paletteer)
library(reshape2) 
library(scales)
library(ggrepel)
library(clusterProfiler)
library(reshape2) 
library(enrichplot)

# --- Load in Data From Previous Step ------------------------------------------
# loads data from previous script "differential_expression.R"
load(paste0("./results/differential-expression.RData"))

# --- Two Group DEG Enrichment Analysis ----------------------------------------

applyEnrichment <- function(
  results,
  ont = "ALL",
  pvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "none",
  qvalueCutoff = 1,
  OrgDb='org.Hs.eg.db',
  keyType='SYMBOL',
  pval_col="Pr(>|t|)",
  lfc_col="Estimate",
  fdr_col="FDR",
  gene_col="Gene",
  contrast_col="Contrast"
){
  # change names to expected format
  colnames(results)[colnames(results) == pval_col] = "Pr(>|t|)"
  colnames(results)[colnames(results) == fdr_col] = "FDR"
  colnames(results)[colnames(results) == lfc_col] = "Estimate"
  colnames(results)[colnames(results) == gene_col] = "Gene"
  colnames(results)[colnames(results) == contrast_col] = "Contrast"
  
  # change contrast column to character
  results$Contrast <- as.character(results$Contrast)
  
  # filter the results object to enrich
  res_filt <- results %>% 
    filter(FDR <0.05)
  
  # run over-representation test on results
  enriched <- enrichGO(res_filt$Gene,
                       universe = unique(results$Gene),
                       ont = ont,
                       pvalueCutoff = pvalueCutoff,
                       minGSSize = minGSSize,
                       maxGSSize = maxGSSize,
                       pAdjustMethod = pAdjustMethod,
                       qvalueCutoff = qvalueCutoff,
                       OrgDb=OrgDb,
                       keyType=keyType)
  
  # run pairwise termism to determine overlap in genes between annotations
  enriched <- pairwise_termsim(enriched)
  
  # create a dotplot
  dot <- dotplot(enriched,showCategory=10, orderBy="Count")
  
  return(list(
    enriched=enriched,
    dotplot=dot
  ))
}

# --- Cluster DEG Enrichment Analysis ------------------------------------------

applyClusterEnrichment <- function(
  results,
  formula,
  enrich_fun='enrichGO',
  ont = "ALL",
  pvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "none",
  qvalueCutoff = 1,
  OrgDb='org.Hs.eg.db',
  keyType='SYMBOL',
  pval_col="Pr(>|t|)",
  lfc_col="Estimate",
  fdr_col="FDR",
  gene_col="Gene",
  contrast_col="Contrast"
){
  # change names to expected format
  colnames(results)[colnames(results) == pval_col] = "Pr(>|t|)"
  colnames(results)[colnames(results) == fdr_col] = "FDR"
  colnames(results)[colnames(results) == lfc_col] = "Estimate"
  colnames(results)[colnames(results) == gene_col] = "Gene"
  colnames(results)[colnames(results) == contrast_col] = "Contrast"
  
  # add in a change direction column
  results$change_direction <- ifelse(results$Estimate>0,"Upregulated","Downregulated")
  
  # change contrast column to character
  results$Contrast <- as.character(results$Contrast)
  
  # filter the results object to enrich
  res_filt <- results %>% 
    filter(FDR <0.05)
  
  # run over-representation test on results
  enriched <- compareCluster(formula, 
                             data=res_filt,
                             fun=enrich_fun,
                             universe = unique(results$Gene),
                             ont = ont,
                             pvalueCutoff = pvalueCutoff,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             pAdjustMethod = pAdjustMethod,
                             qvalueCutoff = qvalueCutoff,
                             OrgDb=OrgDb,
                             keyType=keyType)
  
  # run pairwise termism to determine overlap in genes between annotations
  enriched <- pairwise_termsim(enriched)
  
  # create a dotplot
  dot <- dotplot(enriched)
  
  # ensure the p.adjusted column is an adjusted p-value
  enriched@compareClusterResult$p.adjust=enriched@compareClusterResult$qvalue
  
  return(list(
    enriched=enriched,
    dotplot=dot
  ))
}

disarray_pheno_combined_enrichment_cd <- applyClusterEnrichment(
  results = disarray_pheno_combined_results,
  formula = as.formula("Gene~Contrast+change_direction"))

# isolate only significant go terms
disarray_pheno_combined_enrichment_cd_sig <- disarray_pheno_combined_enrichment_cd$enriched %>%
  mutate(Contrast=gsub("_",":",Contrast)) %>%
  mutate(Cluster=gsub("_",":",Cluster)) %>%
  filter(qvalue<0.1)

# --- DotPlot of Enriched Terms ------------------------------------------------
dotplot(
  object = disarray_pheno_combined_enrichment_cd_sig %>%
    filter(grepl("CONTROL:Normal - HCM:Normal|HCM:Normal - HCM:Severe",Contrast)),
  showCategory=10,
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
# --- Expression of Module Genes -----------------------------------------------

moduleGenePlot <- function(
  NanoStringGeoMxSet,
  geneset,
  meta_data_col,
  meta_levels=NULL,
  title="",
  legend_title="",
  contrast_col="",
  assay="log_q"
){
  meta_exp <- pData(NanoStringGeoMxSet)[meta_data_col] %>%
    rownames_to_column("Var2")
  gene_exp <- assayDataElement(NanoStringGeoMxSet[rownames(assayDataElement(NanoStringGeoMxSet,elt = assay)) %in% geneset,],elt = assay) %>%
    reshape2::melt() %>%
    inner_join(meta_exp,
               by=c("Var2"="Var2"))
  gene_exp[[meta_data_col]] <- gsub(
    "_",
    ":",
    gene_exp[[meta_data_col]]
  )
  colnames(gene_exp)[colnames(gene_exp) == meta_data_col] = "Contrast"
  if(!is.null(meta_levels)){
    gene_exp$Contrast <- factor(gene_exp$Contrast,
                                levels = meta_levels)
  }else{
    gene_exp <- gene_exp
  }
  
  
  mod_plot <- ggplot(gene_exp,
                     aes(x = Contrast, 
                         y = value,
                         fill = Contrast,
                         color = Contrast
                     )) +
    geom_violin(alpha=.8) +
    ggpubr::theme_pubr(legend = "right")+
    scale_fill_manual(values=colors[c(
      9,8,6,4,3,5
    )])+
    scale_color_manual(values=colors[c(
      9,8,6,4,3,5
    )],guide = "none")+
    theme(axis.text.x = element_text(angle = 45,hjust=1))+
    labs(
      x="",
      y = "Module Expression",
      fill="",
      title=title
    )
  return(mod_plot)
}

# atp module expression
atp_row <- disarray_pheno_combined_enrichment_cd_sig@compareClusterResult %>%
  filter(Contrast == "CONTROL:Normal - HCM:Normal" & 
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
interfereon_row <- disarray_pheno_combined_enrichment_cd_sig@compareClusterResult %>%
  filter(Contrast == "CONTROL:Normal - HCM:Normal" & 
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
ext_cell_mat_row <- disarray_pheno_combined_enrichment_cd_sig@compareClusterResult %>%
  filter(Contrast == "HCM:Normal - HCM:Severe" & 
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
# --- Save Data ----------------------------------------------------------------
save.image(
  file = paste0("./results/deg-enrichment.RData")
)
