# Spatial Transcriptomic Analysis of Focal and Normal Areas of Myocyte Disarray in Human Hypertrophic Cardiomyopathy

by Jason Laird 1, Gayani Perera 2, Rebecca Batorsky 3, Hongjie Wang 3, Knarik Arkun 4 and Michael T. Chin 2*
1. Research Technology, Tufts University, Medford, MA 02144, USA
2. Molecular Cardiology Research Institute, Tufts Medical Center, Boston, MA 02111, USA
3. Data Intensive Studies Center, Tufts University, Medford, MA 02155, USA
4. Department of Pathology, Tufts Medical Center, Boston, MA 02111, USA

\* Author to whom correspondence should be addressed.

Code used in the analysis of Hypertrophic Cardiomyopathy Spatial Transcriptomic Data (https://doi.org/10.3390/ijms241612625)

## Pipeline 
1. `scripts/01_geomx_preprocessing.R`: GeoMx pre-processing (probe, segment, Gene QC, etc.)
2. `scripts/02_differential_expression.R`: Mixed model differential dxpression, volcano plots of DEGs, top DEGs for comparisons of interest, ROI distribution
3. `scripts/03_deg_enrichment.R`: Gene ontology enrichment of significant DEGs, dot plot of enriched terms, module expression plots
4. `scripts/04_ligand_receptor.R`: Ligand receptor analysis, CellChat DB pathways
5. `scripts/05_deconvolution.R`: Deconvolution using matched snRNA-seq data and overlapping cell type markers, bar plot of cell type proportions, dot plot of cell type proportions, scatterplot of cell type diversity statistic
6. `scripts/06_supplemental.R`: Supplemental figures, overlap of significant DEGs between all comparisons, overall gene ontology enrichment dotplot, UMAP plot of snRNA-seq cell types, expresison of overlapping marker genes, supplemental table creation
