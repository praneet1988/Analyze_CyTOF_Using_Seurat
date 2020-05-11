# Analyze CyTOF Using Seurat
Cytometry by time-of-flight(CyTOF) data is very useful in studying the presence/absence of antigens/surface markers at single cell level. There are multiple tools to analyze CyTOF data but here I am presenting a tutorial of how one can quickly use Seurat (R package for scRNA-Seq analysis) [https://satijalab.org/seurat/] for analyzing CyTOF data and understanding the cellular and phenotypic diversity at cellular level. 

# There are four major steps:
Step1: Reading .fcs files using read.flowSet function from flowCore R package

Step2: Using Arcsinh transformation to normalize the fcs files. Read more [https://support.cytobank.org/hc/en-us/articles/206148057-About-the-Arcsinh-transform]

Step3: Create a Seurat Object using normalized counts from .fcs files

Step4: Run Dimensionality reduction, clustering and then Visualize cells

# Step1: Reading .fcs files using read.flowSet function from flowCore R package

Let's read the .fcs files using read.flowSet

fcs_raw1 <- read.flowSet('fcs_raw1.fcs', path = getwd(), transformation = FALSE,  truncate_max_range = FALSE)

fcs_raw2 <- read.flowSet('fcs_raw2.fcs', path = getwd(), transformation = FALSE,  truncate_max_range = FALSE)

Let's read the panel file which provides information on the markers

panel_filename <- "CyTOF_panel.xlsx"

panel <- read_excel(panel_filename)

panel$Antigen <- gsub("/", "_", panel$Antigen)

panel$Antigen <- gsub("-", "_", panel$Antigen)

panel_fcs_raw1 <- pData(parameters(fcs_raw1[[1]]))

panel_fcs_raw2 <- pData(parameters(fcs_raw2[[1]]))

panel_fcs_raw1$desc <- gsub("/", "_", panel_fcs_raw1$desc)

panel_fcs_raw1$desc <- gsub("-", "_", panel_fcs_raw1$desc)

panel_fcs_raw2$desc <- gsub("/", "_", panel_fcs_raw2$desc)

panel_fcs_raw2$desc <- gsub("-", "_", panel_fcs_raw2$desc)

(lineage_markers <- panel$Antigen[panel$Lineage == 1])

(functional_markers <- panel$Antigen[panel$SurfaceMarkers == 1])

# Step2: Using Arcsinh transformation to normalize the fcs files

Let's perform the Arcsinh transformation
fcs_1 <- fsApply(fcs_raw1, function(x, cofactor = 5){
  colnames(x) <- panel_fcs_raw1$desc
  expr <- exprs(x)
  expr <- asinh(expr[, c(panel_fcs_raw1$desc)] / cofactor)
  exprs(x) <- expr
  x
})
fcs_1

fcs_2 <- fsApply(fcs_raw2, function(x, cofactor = 5){
  colnames(x) <- panel_fcs_raw2$desc
  expr <- exprs(x)
  expr <- asinh(expr[, c(panel_fcs_raw2$desc)] / cofactor)
  exprs(x) <- expr
  x
})
fcs_2
