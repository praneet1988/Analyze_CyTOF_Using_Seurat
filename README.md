# Analyze CyTOF Using Seurat
Cytometry by time-of-flight(CyTOF) data is very useful in studying the presence/absence of antigens/surface markers at single cell level. There are multiple tools to analyze CyTOF data but here I am presenting a tutorial of how one can quickly use Seurat (R package for scRNA-Seq analysis) [https://satijalab.org/seurat/] for analyzing CyTOF data and understanding the cellular and phenotypic diversity at cellular level. 

# There are five major steps:
Step1: Reading .fcs files using read.flowSet function from flowCore R package

Step2: Using Arcsinh transformation to normalize the fcs files. Read more [https://support.cytobank.org/hc/en-us/articles/206148057-About-the-Arcsinh-transform]

Step3: Create a Seurat Object using normalized counts from .fcs files

Step4: Run Dimensionality reduction, clustering and then visualize cells

Step5: Integrate Seurat objects to understand simmilarities and differences among cell-types learnt from each .fcs file

# R packages required

flowCore: https://www.bioconductor.org/packages/release/bioc/html/flowCore.html

Seurat: https://satijalab.org/seurat/install.html

readxl: https://readxl.tidyverse.org/

cowplot: https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html

ggplot2: https://cran.r-project.org/web/packages/ggplot2/index.html

# Step1: Reading .fcs files using read.flowSet function from flowCore R package

## Let's read the .fcs files using read.flowSet. [Here I am showing analysis tutorial for 2 fcs samples]

fcs_raw1 <- read.flowSet('fcs_raw1.fcs', path = getwd(), transformation = FALSE,  truncate_max_range = FALSE)

fcs_raw2 <- read.flowSet('fcs_raw2.fcs', path = getwd(), transformation = FALSE,  truncate_max_range = FALSE)

## Let's read the panel file which provides information on the markers

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

## Let's perform the Arcsinh transformation for fcs_raw_1 and fcs_raw_2

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

# Create Seurat Object using normalized expression from fcs files

## Let's create Seurat objects for two fcs samples:

### Seurat Object for first fcs sample

expr_fcs_1 <- fsApply(fcs_1, exprs)

expr_fcs_1 <- as.matrix(t(expr_fcs_1))

Cells <- c()

SampleName <- c()

TimePoint <- c()

for (f in 1:ncol(expr_fcs_1))

{

	a <- paste("sample_",f, sep="")
	
	Cells <- c(Cells,a)
	
	a <- paste("sample")
	
	SampleName <- c(SampleName, a)
	
	a <- paste("time")
	
	TimePoint <- c(TimePoint, a)
	
}

colnames(expr_fcs_1) <- Cells

metadata <- cbind(Cells, SampleName, TimePoint)

metadata <- data.frame(metadata)

pd <- new("AnnotatedDataFrame", data = metadata)

rownames(pd) <- pd$Cells

Obj_fcs_1 <- CreateSeuratObject(expr_fcs_1)

CellsMeta = Obj_fcs_1@meta.data

CellsMeta["SampleName"] <- pd$SampleName

CellsMetaTrim <- subset(CellsMeta, select = c("SampleName"))

Obj_fcs_1 <- AddMetaData(Obj_fcs_1, CellsMetaTrim)

CellsMeta = Obj_fcs_1@meta.data

CellsMeta["TimePoint"] <- pd$TimePoint

CellsMetaTrim <- subset(CellsMeta, select = c("TimePoint"))

Obj_fcs_1 <- AddMetaData(Obj_fcs_1, CellsMetaTrim)

### Seurat Object for second fcs sample

expr_fcs_2 <- fsApply(fcs_2, exprs)

expr_fcs_2 <- as.matrix(t(expr_fcs_2))

Cells <- c()

SampleName <- c()

TimePoint <- c()

for (f in 1:ncol(expr_fcs_2))

{

	a <- paste("sample_",f, sep="")
	
	Cells <- c(Cells,a)
	
	a <- paste("sample")
	
	SampleName <- c(SampleName, a)
	
	a <- paste("time")
	
	TimePoint <- c(TimePoint, a)
	
}

colnames(expr_fcs_2) <- Cells

metadata <- cbind(Cells, SampleName, TimePoint)

metadata <- data.frame(metadata)

pd <- new("AnnotatedDataFrame", data = metadata)

rownames(pd) <- pd$Cells

Obj_fcs_2 <- CreateSeuratObject(expr_fcs_2)

CellsMeta = Obj_fcs_2@meta.data

CellsMeta["SampleName"] <- pd$SampleName

CellsMetaTrim <- subset(CellsMeta, select = c("SampleName"))

Obj_fcs_2 <- AddMetaData(Obj_fcs_2, CellsMetaTrim)

CellsMeta = Obj_fcs_2@meta.data

CellsMeta["TimePoint"] <- pd$TimePoint

CellsMetaTrim <- subset(CellsMeta, select = c("TimePoint"))

Obj_fcs_2 <- AddMetaData(Obj_fcs_2, CellsMetaTrim)

# Step4: Variable feature selection, dimension reduction, clustering and visualization

## For Obj_fcs_1
VariableFeatures(Obj_fcs_1) <- rownames(Obj_fcs_1)

Obj_fcs_1 <- ScaleData(Obj_fcs_1)

Obj_fcs_1 <- RunPCA(Obj_fcs_1, verbose = TRUE)

Obj_fcs_1 <- FindNeighbors(Obj_fcs_1, dims = 1:10, verbose = TRUE)

Obj_fcs_1 <- FindClusters(Obj_fcs_1, verbose = TRUE)

Obj_fcs_1 <- RunUMAP(Obj_fcs_1, dims = 1:10, verbose = TRUE)

Obj_fcs_1 <- RunTSNE(Obj_fcs_1, dims = 1:10, verbose = TRUE)

DimPlot(Obj_fcs_1, reduction="umap", pt.size=1)

DimPlot(Obj_fcs_1, reduction="tsne", pt.size=1)

## For Obj_fcs_2

VariableFeatures(Obj_fcs_2) <- rownames(Obj_fcs_2)

Obj_fcs_2 <- ScaleData(Obj_fcs_2)

Obj_fcs_2 <- RunPCA(Obj_fcs_2, verbose = TRUE)

Obj_fcs_2 <- FindNeighbors(Obj_fcs_2, dims = 1:10, verbose = TRUE)

Obj_fcs_2 <- FindClusters(Obj_fcs_2, verbose = TRUE)

Obj_fcs_2 <- RunUMAP(Obj_fcs_2, dims = 1:10, verbose = TRUE)

Obj_fcs_2 <- RunTSNE(Obj_fcs_2, dims = 1:10, verbose = TRUE)

DimPlot(Obj_fcs_2, reduction="umap", pt.size=1)

DimPlot(Obj_fcs_2, reduction="tsne", pt.size=1)

# Step5: Integration of Obj_fcs_1 and Obj_fcs_2 to uncover similarities and differences among cell types from each sample

Integration.anchors <- FindIntegrationAnchors(object.list = list(Obj_fcs_1, Obj_fcs_2), dims = 1:20)

Integration.combined <- IntegrateData(anchorset = Integration.anchors, dims = 1:20)

DefaultAssay(Integration.combined) <- "integrated"

Integration.combined <- ScaleData(Integration.combined)

Integration.combined <- RunPCA(Integration.combined, verbose = TRUE)

Integration.combined <- RunUMAP(Integration.combined, reduction = "pca", dims = 1:10)

Integration.combined <- RunTSNE(Integration.combined, reduction = "pca", dims = 1:10)

Integration.combined <- FindNeighbors(Integration.combined, reduction = "pca", dims = 1:10)

Integration.combined <- FindClusters(Integration.combined, resolution = 0.5)

DimPlot(Integration.combined, pt.size=1)

# UMAP Plot after integration (Please Note this is example plot)
![Graph](IntegrationUMAPPlot.png)

# All the steps here are based on the analysis for determining cell types from CyTOF data and if you are using these steps then please change parameters accordingly
