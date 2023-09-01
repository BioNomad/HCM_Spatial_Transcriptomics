# --- Load Libraries -----------------------------------------------------------
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
library(ComplexUpset)



# --- Load in Data From Previous Step ------------------------------------------
# loads data from previous script "geomx_preprocessing.R"

load(paste0("./results/geomx-preporcessing.RData"))

# --- Differential Expression --------------------------------------------------

# convert test variables to factors
pData(target_combineData)$batch <- 
  factor(pData(target_combineData)$batch, c("batch_1", "batch_2"))

pData(target_combineData)$Segment.tags.pathologist_1 <- 
  factor(pData(target_combineData)$Segment.tags.pathologist_1,
         c("Severe", "Moderate", "Mild", "Normal"))

pData(target_combineData)$slide.name<- 
  factor(pData(target_combineData)[["slide name"]])

pData(target_combineData)$Phenotype <- 
  factor(pData(target_combineData)$Phenotype, c("HCM","CONTROL"))

# add in column for combined phenotype-disarray
pData(target_combineData)$pheno_disarray <- 
  factor(paste(
    pData(target_combineData)$Phenotype,
    pData(target_combineData)$Segment.tags.pathologist_1,
    sep="_"
  ),
  levels = c(
    "HCM_Severe",
    "HCM_Moderate",
    "HCM_Mild",
    "HCM_Normal",
    "CONTROL_Mild",
    "CONTROL_Normal"
  ))

assayDataElement(object = target_combineData, elt = "log_q") <-
  assayDataApply(target_combineData, 2, FUN = log, base = 2, elt = "q_norm")

# define a function to perform mixed model differential expression
# and add in additional data columns
mixedModelApply <- function(
  NanoStringGeoMxSet,
  assay="log_q",
  formula,
  groupVar
){
  results <- c()
  mixedOutmc <-
    GeomxTools::mixedModelDE(NanoStringGeoMxSet,
                             elt = assay,
                             modelFormula = formula,
                             groupVar = groupVar,
                             multiCore = FALSE)
  
  # format results as data.frame
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  # use lapply in case you have multiple levels of your test factor to
  # correctly associate gene name with it's row in the results table
  r_test$Gene <-
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Contrast", "Estimate",
                       "Pr(>|t|)", "FDR")]
  results <- rbind(results, r_test)
  rownames(results) <- NULL
  results$change_direction <- ifelse(results$Estimate>0,"Upregulated","Downregulated")
  
  # swap names of results for interpretability
  results$Contrast <- sapply(results$Contrast,
                             function(x){
                               return(
                                 paste(
                                   strsplit(x," - ")[[1]][2],
                                   strsplit(x," - ")[[1]][1],
                                   sep=" - "
                                 )
                               )
                             })
  
  return(results)
}

# test for degs between phenotype+disarray
disarray_pheno_combined_results <- mixedModelApply(
  NanoStringGeoMxSet=target_combineData,
  assay="log_q",
  formula=as.formula(
    " ~ pheno_disarray + (1 | slide.name)"
  ),
  groupVar="pheno_disarray")

# --- Plot Where ROI's Are Coming From -----------------------------------------

roi_patient_df <- pData(target_combineData) %>%
  mutate(count=1) %>%
  mutate(Phenotype=gsub("NORMAL","CONTROL",Phenotype))

# plot patient v. ROI count, colored in by Disarray Level

patient_roi_plot <- ggplot(
  roi_patient_df,
  aes(
    y=`slide name`,
    x=count,
    fill=Segment.tags.pathologist_1
  )
)+
  geom_bar(stat = "identity")+
  theme_bw()+
  scale_fill_manual(values= c(
    "#6b2a3a",
    "#ba8693",
    "#869dba",
    "#2a4d6b"
  ))+
  theme(axis.text.x = element_text(angle = 45,hjust=1),
        text = element_text(size = 12))+
  facet_grid(Phenotype~.,
             drop=TRUE,
             scales='free_y', 
             space = "free")+
  labs(
    x="ROI Count",
    y="Patient",
    fill="Disarray Level",
    title = "ROI Count Per Patient"
  )

patient_roi_plot

# --- Volcano Plot -------------------------------------------------------------

volcanoPlot <- function(
  results,
  title,
  pval_thresh = 0.05,
  fdr_thresh = 0.05,
  fdr_strict_thresh = 0.001,
  pval_col="Pr(>|t|)",
  lfc_col="Estimate",
  fdr_col="FDR",
  gene_col="Gene",
  contrast_col="Contrast"
){
  require(ggrepel) 
  # change names to expected format
  colnames(results)[colnames(results) == pval_col] = "Pr(>|t|)"
  colnames(results)[colnames(results) == fdr_col] = "FDR"
  colnames(results)[colnames(results) == lfc_col] = "Estimate"
  colnames(results)[colnames(results) == gene_col] = "Gene"
  colnames(results)[colnames(results) == contrast_col] = "Contrast"
  
  # Categorize Results based on P-value & FDR for plotting
  results$Color <- "Not Significant"
  results$Color[results$`Pr(>|t|)` < pval_thresh] <- paste(
    "P < ",as.character(pval_thresh),sep = "")
  results$Color[results$FDR < fdr_thresh] <- paste(
    "FDR < ",as.character(fdr_thresh),sep = "")
  results$Color[results$FDR < fdr_strict_thresh] <- paste(
    "FDR < ",as.character(fdr_strict_thresh),sep = "")
  results$Color <- factor(results$Color,
                          levels = c(
                            "Not Significant", 
                            paste("P < ",
                                  as.character(pval_thresh),
                                  sep = ""),
                            paste("FDR < ",
                                  as.character(fdr_thresh),
                                  sep = ""),
                            paste("FDR < ",
                                  as.character(fdr_strict_thresh),
                                  sep = "")
                          ))
  
  # pick top genes for either side of volcano to label
  # order genes for convenience:
  results$invert_P <- (-log10(results[[pval_col]])) * sign(results[[lfc_col]])
  top_g <- c(
    results[, gene_col][
      order(results[, 'invert_P'], decreasing = TRUE)[1:15]],
    results[, gene_col][
      order(results[, 'invert_P'], decreasing = FALSE)[1:15]])
  top_g <- unique(top_g)
  results <- results[, -1*ncol(results)] # remove invert_P from matrix
  
  # Graph results
  volcano <- ggplot(results,
                    aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                        color = Color, label = Gene)) +
    geom_hline(yintercept = -log10(pval_thresh), lty = "dashed") +
    geom_point() +
    labs(x = "log2(FC)",
         y = "Significance, -log10(P)",
         color = "Significance",
         title = title) +
    scale_color_manual(values = c("grey",
                                  "#006DDBFF",
                                  "#FF6DB6FF",
                                  "#000000FF"),
                       guide = guide_legend(override.aes = list(size = 4))) +
    ggrepel::geom_label_repel(data = subset(results, Gene %in% top_g & FDR < fdr_strict_thresh),
                              color="darkslategrey",max.overlaps = 50,size=3,
                              min.segment.length = .1) +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom") +
    facet_wrap(~Contrast, scales = "free_y")
  return(volcano)
}

volcanoPlot(results = disarray_pheno_combined_results %>%
              mutate(Contrast=gsub("_",":",Contrast)) %>%
              filter(Contrast %in% c(
                "CONTROL:Normal - HCM:Normal",
                "HCM:Mild - HCM:Moderate",
                "HCM:Moderate - HCM:Severe"
              )),
            title = "Gene Expression Changes Between Phenotypes and HCM Disarray Levels")
# --- Plot DEG Barplot ---------------------------------------------------------
# define a function to plot DEGs per comparison group

degSumPlot <- function(
  results,
  title="",
  pval_thresh = 0.05,
  fdr_thresh = 0.05,
  fdr_strict_thresh = 0.001,
  pval_col="Pr(>|t|)",
  lfc_col="Estimate",
  fdr_col="FDR",
  gene_col="Gene",
  contrast_col="Contrast"
){
  require(ggrepel) 
  # change names to expected format
  colnames(results)[colnames(results) == pval_col] = "Pr(>|t|)"
  colnames(results)[colnames(results) == fdr_col] = "FDR"
  colnames(results)[colnames(results) == lfc_col] = "Estimate"
  colnames(results)[colnames(results) == gene_col] = "Gene"
  colnames(results)[colnames(results) == contrast_col] = "Contrast"
  
  # Categorize Results based on P-value & FDR for plotting
  results$Color <- "Not Significant"
  results$Color[results$`Pr(>|t|)` < pval_thresh] <- paste(
    "P < ",as.character(pval_thresh),sep = "")
  results$Color[results$FDR < fdr_thresh] <- paste(
    "FDR < ",as.character(fdr_thresh),sep = "")
  results$Color[results$FDR < fdr_strict_thresh] <- paste(
    "FDR < ",as.character(fdr_strict_thresh),sep = "")
  results$Color <- factor(results$Color,
                          levels = c(
                            "Not Significant", 
                            paste("P < ",
                                  as.character(pval_thresh),
                                  sep = ""),
                            paste("FDR < ",
                                  as.character(fdr_thresh),
                                  sep = ""),
                            paste("FDR < ",
                                  as.character(fdr_strict_thresh),
                                  sep = "")
                          ))
  deg_sum <- results %>%
    filter(Color !="Not Significant") %>%
    count(Color,Contrast) 
  
  # plot deg bar plot
  degs<- ggplot(deg_sum, 
                aes(x=n, 
                    y=reorder(Contrast,n),
                    fill=Color)) +
    geom_bar(stat="identity")+
    theme_bw() +
    scale_fill_manual(
      values=c("#006DDBFF",
               "#FF6DB6FF",
               "#000000FF"
      )
    )+
    theme(axis.text.x = element_text(angle = 45,hjust=1),
          text = element_text(size = 12)
    )+
    labs(
      x="DEG Count",
      y="",
      title = title,
      fill=""
    )
  return(degs)
}

degSumPlot(results = disarray_pheno_combined_results,
           title = "~ Phenotype_Disarray + (1|Patient)")

# --- Plot Genes of Interest ---------------------------------------------------

genePlot <- function(
  NanoStringGeoMxSet,
  gene,
  meta_data_col,
  title="",
  contrast_col="",
  assay="log_q"
){
  if(contrast_col==""){
    plotting_df <- data.frame(
      meta=as.character(pData(target_combineData)[[meta_data_col]]),
      gene=as.numeric(assayDataElement(target_combineData[gene,],
                                       elt = assay))
    )
    gene_plot <- ggplot(plotting_df,
                        aes(x = meta, 
                            color = meta,
                            y = gene
                        )) +
      geom_boxplot() +
      geom_jitter() +
      labs(y = paste(gene,"Expression")) +
      theme_bw() +
      labs(
        x="",
        color="",
        title=title
      )
  }else{
    plotting_df <- data.frame(
      meta=as.character(pData(target_combineData)[[meta_data_col]]),
      gene=as.numeric(assayDataElement(target_combineData[gene,],
                                       elt = assay)),
      contrast=as.character(pData(target_combineData)[[contrast_col]])
    )
    gene_plot <- ggplot(plotting_df,
                        aes(x = meta, 
                            color = contrast,
                            y = gene
                        )) +
      geom_boxplot() +
      geom_jitter() +
      labs(y = paste(gene,"Expression")) +
      theme_bw() +
      labs(
        x="",
        color="",
        title=title
      )
  }
  return(gene_plot)
}

# confirm direction of change
genePlot(NanoStringGeoMxSet=target_combineData,
         gene="GOT2",
         meta_data_col="Segment.tags.pathologist_1",
         assay="log_q",
         contrast_col = "Phenotype",
         title = "")


# --- Save Data ----------------------------------------------------------------
save.image(
  file = paste0("./results/differential-expression.RData")
)
