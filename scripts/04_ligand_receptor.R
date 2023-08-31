# --- Load Libraries -----------------------------------------------------------
LIB1='/cluster/tufts/hpc/tools/R/4.0.0'
LIB2='/cluster/home/jlaird01/R/x86_64-pc-linux-gnu-library/4.0'
.libPaths(c(LIB1,LIB2))
library(tidymodels)
library(broom.mixed)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(tidyverse)
library(patchwork)
library(ggnewscale) 
library(paletteer)
library(clusterProfiler)
library(OmnipathR)
library(MASS)


# --- Load in Data From Previous Step ------------------------------------------
# loads data from previous script "deg_enrichment.R"
load(paste0("./results/deg-enrichment.RData"))

# --- Load CellChat Data -------------------------------------------------------
load("/cluster/tufts/chinlab/jlaird01/nanostring/data/cellchat/CellChatDB.human.rda")


# --- Ligand Receptor Modelling ------------------------------------------------
# clean up composite variable for plotting
pData(target_combineData)$composite <- factor(paste(
  pData(target_combineData)$Phenotype,
  pData(target_combineData)$Segment.tags.pathologist_1,
  sep=": "
),
levels = c("CONTROL: Normal",
           "CONTROL: Mild",
           "HCM: Normal",
           "HCM: Mild",
           "HCM: Moderate",
           "HCM: Severe"))

pData(target_combineData)$composite_slide <- paste(
  pData(target_combineData)$composite,
  pData(target_combineData)$slide.name,
  sep="_"
)

# --- Differential Combination Pathways ----------------------------------------
# upregulated lr genes
diff_comb_hits_up <- disarray_pheno_combined_results %>%
  filter(FDR < 0.05 & Estimate > 0 &
           !grepl("CONTROL_Mild",Contrast))
diff_comb_hits_up=split(diff_comb_hits_up$Gene, diff_comb_hits_up$Contrast)

# downregulated lr genes
diff_comb_hits_down <- disarray_pheno_combined_results %>%
  filter(FDR < 0.05 & Estimate < 0 &
           !grepl("CONTROL_Mild",Contrast))
diff_comb_hits_down=split(diff_comb_hits_down$Gene, diff_comb_hits_down$Contrast)

# identify lr interactions in upregulated degs
diff_comb_ints_up = lapply(
  diff_comb_hits_up,
  function(x){
    ints=strsplit(CellChatDB.human$interaction$interaction_name,"_");
    res=c();
    for (int in ints){
      present_lr=all(int %in% x)
      res=c(res,present_lr)
    };
    lr_hits=CellChatDB.human$interaction[res,];
    return(lr_hits)
  }
)

# identify lr interactions in downregulated degs
diff_comb_ints_down = lapply(
  diff_comb_hits_down,
  function(x){
    ints=strsplit(CellChatDB.human$interaction$interaction_name,"_");
    res=c();
    for (int in ints){
      present_lr=all(int %in% x)
      res=c(res,present_lr)
    };
    lr_hits=CellChatDB.human$interaction[res,];
    return(lr_hits)
  }
)

# remove 0 length lists
diff_comb_ints_up = diff_comb_ints_up[lapply(diff_comb_ints_up,nrow) !=0]
diff_comb_ints_down = diff_comb_ints_down[lapply(diff_comb_ints_down,nrow) !=0]

# add in up/downregulated column
diff_comb_ints_up = lapply(diff_comb_ints_up,function(x){return(x %>% mutate(group="Upregulated"))})
diff_comb_ints_down = lapply(diff_comb_ints_down,function(x){return(x %>% mutate(group="Downregulated"))})

# collapse data frames
diff_comb_ints_up=map_df(diff_comb_ints_up, ~as.data.frame(.x), .id="id")
diff_comb_ints_down=map_df(diff_comb_ints_down, ~as.data.frame(.x), .id="id")

# combine into one data frame
diff_comb_ints = rbind(
  diff_comb_ints_up,
  diff_comb_ints_down
) %>%
  mutate(id=gsub("_",":",id)) %>%
  mutate(id=factor(id,
                   levels = c(
                     "CONTROL:Normal - HCM:Normal",
                     "CONTROL:Normal - HCM:Mild",
                     "CONTROL:Normal - HCM:Moderate",
                     "CONTROL:Normal - HCM:Severe",
                     "HCM:Normal - HCM:Moderate",
                     "HCM:Normal - HCM:Severe"
                   ))) %>%
  mutate(pathway=str_replace_all(pathway_name,c("PDGF" = "Platelet-Derived Growth Factor Signaling", 
                                   "NT" = "Neurotrophin Signaling", 
                                   "NOTCH" = "NOTCH Signaling",
                                   "JAM"="Junctional Adhesion Molecules",
                                   "CD99" = "CD99 Signaling",
                                   "CD46" = "CD46 Signaling",
                                   "FN1"="Fibronectin Signaling",
                                   "CDH" ="Cadherin Signaling",
                                   "APP"="Amyloid Precursor Protein Signaling")))

diff_comb_ints_sum = diff_comb_ints %>%
  group_by(id,group,pathway) %>%
  summarise(count=length(interaction_name)) %>%
  ungroup()%>%
  mutate(upregulated_count=ifelse(group=="Upregulated",count,0),
         downregulated_count=ifelse(group=="Downregulated",count,0))

ggplot() +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1,
            data=diff_comb_ints_sum[diff_comb_ints_sum$upregulated_count != 0,], 
            aes(x=factor(id,levels = levels(diff_comb_ints_sum$id)), 
                y=pathway, 
                fill=upregulated_count)) +
  scale_fill_gradient("Number of Ligand-Receptors",
                      low="brown1", 
                      high="brown4") +
  new_scale("fill") +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1,
            data=diff_comb_ints_sum[diff_comb_ints_sum$downregulated_count != 0,],
            aes(x=id, 
                y=pathway,
                fill=downregulated_count)) +
  scale_fill_gradient("Number of Ligand-Receptors",low="dodgerblue", high="midnightblue")+
  ggpubr::theme_pubr()+
  ggpubr::rotate_x_text(angle = 45)+
  facet_grid(.~group,scales = "free",space = "free")+
  labs(
    x="",
    y=""
  )
    
# --- Save Data ----------------------------------------------------------------

save.image(
  file = paste0("./results/ligand-receptor.RData")
)
