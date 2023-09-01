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

# set colors

colors <- paletteer_d("colorBlindness::paletteMartin")
colors <- colors[-1]

# following this tutorial:
# https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html#4_QC__Pre-processing


# --- Seattle Data -------------------------------------------------------------

dcc_seattle = "/cluster/tufts/chinlab/jlaird01/nanostring/data/seattle_raw_files/DCC_newPKC"
ann_seattle = "/cluster/tufts/chinlab/jlaird01/nanostring/data/seattle_raw_files/Annotation_seattle_newpkc_rebecca.xlsx"
pkc_seattle = "/cluster/tufts/chinlab/jlaird01/nanostring/data/seattle_raw_files/Alpha_CTA_v2.1.pkc"

# --- Grafton data -------------------------------------------------------------

dcc_grafton = ./data/grafton_raw_files/DCC_Files"
ann_grafton = "./data/grafton_raw_files/Annotation_grafton_rebecca.xlsx"
pkc_grafton = "./data/grafton_raw_files/GeoMx_Hs_CTA_v1.0.pkc"

# --- combined_annotation ------------------------------------------------------

# load annotation file
ann_combined = "./data/annotation_cleaned/annotation_cleaned_with_two_pathologists.xlsx"


# read in dcc files
DCCFiles_grafton <- dir(dcc_grafton,
                        pattern = ".dcc$",
                        full.names = TRUE, 
                        recursive = TRUE)

# read in dcc files
DCCFiles_seattle <- dir(dcc_seattle, 
                        pattern = ".dcc$",
                        full.names = TRUE, 
                        recursive = TRUE)

# combined dcc, pkc and annotation files into a NanoString object
# using the grafton probeset
combineData <-
  GeomxTools::readNanoStringGeoMxSet(dccFiles = c(DCCFiles_grafton, DCCFiles_seattle),
                                     pkcFiles = pkc_grafton,
                                     phenoDataFile = ann_combined,
                                     phenoDataSheet = "Sheet 1", 
                                     phenoDataDccColName = "Sample_ID",      
                                     protocolDataColNames = c("roi"))   

# check the feature/sample count
dim(combineData)

# --- Segment QC ---------------------------------------------------------------
pkcs <- annotation(combineData)
modules <- gsub(".pkc", "", pkcs)

combineData <- GeomxTools::shiftCountsOne(combineData, useDALogic = TRUE)
combineData <- combineData[,!(pData(combineData)$`slide name` %in% c("2884","2900"))]
pData(combineData)$Phenotype=gsub("NORMAL","CONTROL",pData(combineData)$Phenotype)

QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 80,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 10,   # Minimum negative control counts (10) 
       maxNTCCount = 1000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 100,         # Minimum # of nuclei estimated (100)
       minArea = 5000)         # Minimum segment area (5000)
combineData <-
  GeomxTools::setSegmentQCFlags(combineData, 
                                qcCutoffs = QC_params)        

# Collate QC Results
QCResults <- protocolData(combineData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

# set qc variable to color by
col_by <- "batch"

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

QC_histogram(GeomxTools::sData(combineData), "Trimmed (%)", col_by, 80)


QC_histogram(GeomxTools::sData(combineData), "Saturated (%)", col_by, 50) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")



# calculate the negative geometric means for each module
negativeGeoMeans <- 
  NanoStringNCTools::esBy(
    NanoStringNCTools::negativeControlSubset(combineData), 
    GROUP = "Module", 
    FUN = function(x) { 
      NanoStringNCTools::assayDataApply(
        x, MARGIN = 2, FUN = GeomxTools::ngeoMean, elt = "exprs") 
    }) 
protocolData(combineData)[["NegGeoMean"]] <- negativeGeoMeans


# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(combineData)[, negCols] <- GeomxTools::sData(combineData)[["NegGeoMean"]]
for(ann in negCols) {
  plt <- QC_histogram(pData(combineData), ann, col_by, 2, scale_trans = "log10")
  print(plt)
}

combineData <- combineData[, QCResults$QCStatus == "PASS"]

# Subsetting our dataset has removed samples which did not pass QC
dim(combineData)

# --- Probe QC -----------------------------------------------------------------

# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
combineData <- GeomxTools::setBioProbeQCFlags(combineData, 
                                              qcCutoffs = list(minProbeRatio = 0.1,
                                                               percentFailGrubbs = 20), 
                                              removeLocalOutliers = TRUE)

ProbeQCResults <- fData(combineData)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
#Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
  subset(combineData, 
         fData(combineData)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(combineData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)

# check feature/sample count
dim(ProbeQCPassed)

# save probes that passed qc metrics
combineData <- ProbeQCPassed 

# --- Gene Level Counts --------------------------------------------------------

# number of unique genes
length(unique(featureData(combineData)[["TargetName"]]))

# collapse to targets
target_combineData <- GeomxTools::aggregateCounts(combineData)

# check feature/sample counts
dim(target_combineData)

# check expression data frame
exprs(target_combineData)[1:5, 1:2]


# --- Limit of Quantification --------------------------------------------------

# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_combineData))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_combineData)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_combineData)[, vars[1]] * 
             pData(target_combineData)[, vars[2]] ^ cutoff)
  }
}

# add LOQ value to meta data
pData(target_combineData)$LOQ <- LOQ


# 4.4 Filtering
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_combineData)$Module == module
  Mat_i <- t(esApply(target_combineData[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}

# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_combineData)$TargetName, ]

# Save detection rate information to pheno data
pData(target_combineData)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_combineData)$GeneDetectionRate <-
  pData(target_combineData)$GenesDetected / nrow(target_combineData)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_combineData)$DetectionThreshold <- 
  cut(pData(target_combineData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(target_combineData),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = batch)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")


### Filter segments
target_combineData <-
  target_combineData[, pData(target_combineData)$GeneDetectionRate >= .01]

# check feature/sample count
dim(target_combineData)

# --- Gene Detection Rate/Filtering --------------------------------------------
# gene detection rate
library(scales) # for percent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_combineData)]
fData(target_combineData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_combineData)$DetectionRate <-
  fData(target_combineData)$DetectedSegments / nrow(pData(target_combineData))

# Gene of interest detection table
goi <- c("PDCD1", "CD274", "IFNG", "CD8A", "CD68", "EPCAM",
         "KRT18", "NPHS1", "NPHS2", "CALB1", "CLDN8")
goi_df <- data.frame(
  Gene = goi,
  Number = fData(target_combineData)[goi, "DetectedSegments"],
  DetectionRate = fData(target_combineData)[goi, "DetectionRate"])

## gene filtering
# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_combineData)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_combineData))
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")

# Subset to target genes detected in at least 10% of the samples.
#   Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_combineData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_combineData <- 
  target_combineData[fData(target_combineData)$DetectionRate >= 0.01 |
                       fData(target_combineData)$TargetName %in% neg_probes, ]

# check feature/sample count
dim(target_combineData)

# retain only detected genes of interest
goi <- goi[goi %in% rownames(target_combineData)]

# --- Normalization ------------------------------------------------------------

# Graph Q3 value vs negGeoMean of Negatives

## here change to 
ann_of_interest <- "Segment.tags.pathologist_1"

Stat_data <- 
  data.frame(row.names = colnames(exprs(target_combineData)),
             Segment = colnames(exprs(target_combineData)),
             Annotation = pData(target_combineData)[, ann_of_interest],
             Q3 = unlist(apply(exprs(target_combineData), 2,
                               quantile, 0.75, na.rm = TRUE)),
             NegProbe = exprs(target_combineData)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) + theme_bw() +
  # scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) + 
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- cowplot::plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43,0.57))
cowplot::plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))


# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
target_combineData <- GeomxTools::normalize(target_combineData ,
                                            norm_method = "quant", 
                                            desiredQuantile = .75,
                                            toElt = "q_norm")

# Background normalization for WTA/CTA without custom spike-in
target_combineData <- GeomxTools::normalize(target_combineData ,
                                            norm_method = "neg", 
                                            fromElt = "exprs",
                                            toElt = "neg_norm")

# visualize the first 10 segments with each normalization method
boxplot(exprs(target_combineData)[,1:10],
        col = "#9EDAE5", main = "Raw Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Raw")

boxplot(assayDataElement(target_combineData[,1:10], elt = "q_norm"),
        col = "#2CA02C", main = "Q3 Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Q3 Normalized")


boxplot(assayDataElement(target_combineData[,1:10], elt = "neg_norm"),
        col = "#FF7F0E", main = "Neg Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Neg. Normalized")

# --- Unsupervised Analysis ----------------------------------------------------

# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42

# run UMAP
umap_out <-
  umap::umap(t(log2(assayDataElement(target_combineData , elt = "q_norm"))),  
       config = custom_umap)

# add umap data to geomx meta data
pData(target_combineData)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]

# clean meta data variables
pData(target_combineData)$cleaned_batch <- gsub("_"," ",pData(target_combineData)$batch)
pData(target_combineData)$Disarray <- pData(target_combineData)$Segment.tags.pathologist_1
pData(target_combineData)$Patient <- pData(target_combineData)$`slide name`

# plot umap with patient labels
ggpubr::ggscatter(pData(target_combineData),
                  x="UMAP1",
                  y="UMAP2",
                  color="Patient",
                  shape = "Phenotype",
                  legend="right",
                  label = "id",
                  size = 1
)+
  scale_color_manual(values=colors)

# plot umap with just shapes
ggpubr::ggscatter(pData(target_combineData),
                  x="UMAP1",
                  y="UMAP2",
                  color="Disarray",
                  shape = "Phenotype",
                  legend="right",
                  size = 4
)+
  scale_color_manual(values=colors[c(5,3,4,8)]) 

# check that umap data was added to geomx meta data
head(target_combineData@phenoData@data)

# --- Save Data ----------------------------------------------------------------
save.image(
  file = paste0("./results/geomx-preporcessing.RData")
)
