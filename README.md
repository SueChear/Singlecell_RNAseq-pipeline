# Single cell RNA-Seq data analysis pipeline with an example of reprocessing data from a study investigating susceptibility of multiple types of pancreatic islet cells to SARS-CoV-2.

Download cellranger

```
wget -O cellranger-6.1.2.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.1.2.tar.gz?Expires=1641366506&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjEuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NDEzNjY1MDZ9fX1dfQ__&Signature=kaV8~ZabHhyDykUhbN~F78PDQfNZ64IamgsGc1nOSghFKPr0fbZ3WJk-2eWYh7IEt-KupenYP89W1zHi4lrxF~ZBbuP4NTaKEAa-G6ILJoX-VdyFnktkXFYDHgzEJ8ABq-NM6RWn20WD3a9BITNHTIWPtxjM-NaXAuR5uc5PuAEgjSDaQ2QBAQr~1q4aSM-~vJt~ia5e8acTz9RlM24EluLqfO59VCtAorP-5iJRwvLw9DjfrTlDtWfy3M2LSXp5OGmVJH1WUQReLK~0iZX2e8~vrHlAYpuxMa0Lgil6oHQ5s6vc~Dod3Aqpjb9sM~wuVo80zi4EqJ5nq0LU8SNbiQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
```

Download human reference genome data file

```
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
```

Unzip cellranger

```
Tar -xf cellranger-7.1.0.tar.gz
```

Unzip human reference genome file
```
Tar -xf refdata-gex-GRCh38-2020-A.tar.gz
```

Download SRR12831415, SRR12831416,SRR12831417 and SRR12831418 files from GEO database with accession # GSE159556.
Place Hg38 reference genome folder, cellranger folder in the same folder, hereby named "scrnaseq". Within this folder, create a folder to store all downloaded fastq files.

Run cellranger count pipeline to align sequencing reads in FASTQ files to the reference transcriptome.
Set up the command to run cellranger count.

Specify an --id, which is the folder name for output files, here this directory is called "run_count". The --fastqs is the path to the directory containing the FASTQ files. Use the --sample argument to specify which samples to use. The --transcriptome argument specifies the path to the transcriptome reference package.
```
cd scrnaseq
Cellranger-7.1.0/bin/cellranger count --id=run_count --fastqs=fastqs --sample=SRR12831418 --transcriptome=refdata-gex-GRCh38-2020-A
```
When the output of the cellranger count command says, “Pipestance completed successfully!”, this means the job is done.

Run all of the following code in R.

Install and load required libraries.

```
library(Seurat)
library(hdf5r)
library(ggplot2)
library(gridExtra)
```
Run the Read10X() function to read in the output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix. The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column).

```
h5_files<-list.files(pattern="*.h5", recursive = F, full.names=F)
h5_read<-lapply(h5_files, Read10X_h5)
```
We next use the count matrix to create a Seurat object. The object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results) for a single-cell dataset. 
```
h5_seurat<-lapply(h5_read, CreateSeuratObject, min.cells=3, min.features=200)

merged_h5_seurat<-merge(h5_seurat[[1]], y=h5_seurat[2:length(h5_seurat)],
                        add.cell.ids=c("Mock_1_M1","Cov_1_C1","Mock_2_M2","Cov_2_C2"),
                        project="Covid")
merged_h5_seurat                   
```       
```
##An object of class Seurat 
##26784 features across 40096 samples within 1 assay 
##Active assay: RNA (26784 features, 0 variable features)
```
##QC and selecting cells for further analysis

Every row shows a cell with its barcode, number of transcripts and number of genes
```
head(merged_h5_seurat@meta.data)
```
```
orig.ident nCount_RNA nFeature_RNA
Mock_1_M1_AAACCTGAGACTGGGT-1 SeuratProject       1978         1101
Mock_1_M1_AAACCTGAGAGGTAGA-1 SeuratProject       1363          806
Mock_1_M1_AAACCTGAGGGATGGG-1 SeuratProject       8555         3473
Mock_1_M1_AAACCTGAGTAGCGGT-1 SeuratProject       1329          828
Mock_1_M1_AAACCTGCAACTTGAC-1 SeuratProject       8675         3257
Mock_1_M1_AAACCTGCACTTACGA-1 SeuratProject       1706         1086
```

Create a sample column
```
merged_h5_seurat$sample<-rownames(merged_h5_seurat@meta.data)
```

Split sample column
```
merged_h5_seurat@meta.data<-separate(merged_h5_seurat@meta.data, col='sample',
                                     into=c('Type','Batch','ID','Barcode'),
                                     sep='_')
head(merged_h5_seurat@meta.data)
```
```
orig.ident nCount_RNA nFeature_RNA Type Batch ID            Barcode
Mock_1_M1_AAACCTGAGACTGGGT-1 SeuratProject       1978         1101 Mock     1 M1 AAACCTGAGACTGGGT-1
Mock_1_M1_AAACCTGAGAGGTAGA-1 SeuratProject       1363          806 Mock     1 M1 AAACCTGAGAGGTAGA-1
Mock_1_M1_AAACCTGAGGGATGGG-1 SeuratProject       8555         3473 Mock     1 M1 AAACCTGAGGGATGGG-1
Mock_1_M1_AAACCTGAGTAGCGGT-1 SeuratProject       1329          828 Mock     1 M1 AAACCTGAGTAGCGGT-1
Mock_1_M1_AAACCTGCAACTTGAC-1 SeuratProject       8675         3257 Mock     1 M1 AAACCTGCAACTTGAC-1
Mock_1_M1_AAACCTGCACTTACGA-1 SeuratProject       1706         1086 Mock     1 M1 AAACCTGCACTTACGA-1
```

Calculate percentage of reads that map to the mitochondrial genome
```
merged_h5_seurat$mitoPercent<-PercentageFeatureSet(merged_h5_seurat, pattern='^MT-')
```

Calculate amount of ribosomal genes
```
merged_h5_seurat$riboPercent<-PercentageFeatureSet(merged_h5_seurat,pattern='^RP[SL]')
```

Visualize QC metrics as violin plots
```
VlnPlot(merged_h5_seurat, group.by="ID", features=c("nFeature_RNA","nCount_RNA","mitoPercent","riboPercent"),ncol=3,pt.size=0.1)+NoLegend()
```

Image 1

FeatureScatter is typically used to visualize feature-feature relationships.Quality data should follow the straight line trend, the scattered points on lower
right quadrant indicate some low number of genes that have been sequenced again and again

```
FeatureScatter(merged_h5_seurat, feature1="nCount_RNA", feature2="nFeature_RNA")+
  geom_smooth(method='lm')
```

Image 2


We filter cells that have unique feature counts <500 or >6000 and unique RNA counts <1000 or >60000 and cells that have >15% mitochondrial counts. These thresholds are taken from the published paper.

```
merged_h5_seurat_filtered<-subset(merged_h5_seurat,subset=nFeature_RNA >500 & nFeature_RNA<=6000 & nCount_RNA >1000 & nCount_RNA <= 60000 &
                                    mitoPercent< 15)

```

After removing unwanted cells, we normalize data in order to compare gene expression across multiple cells. Wwe divide gene expression measurement in each cell by the total expression, then multiply it by a scaling factor followed by a  log transformation. Default normalization simply scales the counts by total counts in each cell, multiplies by 10,000 and then log transforms.

```
merged_h5_seurat_filtered<-NormalizeData(merged_h5_seurat_filtered)
```

##Identification of highly variable features (feature selection)
Get a list of the genes which have high cell-to-cell variation in the dataset where they are highly expressed in some cells and lowly expressed in others to be used for downstream analysis. This command returns 2000 features per dataset.

```
merged_h5_seurat_filtered<-FindVariableFeatures(merged_h5_seurat_filtered, selection.method = "vst", nfeatures=2000)
```

Identify the 10 most highly variable genes
```
top10<-head(VariableFeatures(merged_h5_seurat_filtered),10)
top10
```
```
[1] "MMP1"     "IL1B"     "HSPA6"    "COL1A1"   "CSF3"     "IL11"     "MMP3"     "IGFBP5"   "REG3G"    "SERPINB2"
```

Plot variable features with and without labels
```
plot1<-VariableFeaturePlot(merged_h5_seurat_filtered)
LabelPoints(plot=plot1, points=top10, repel=T)
```

Image3



##Scaling the data
We scale the data prior to PCA.
```
merged_h5_seurat_filtered<-ScaleData(object=merged_h5_seurat_filtered)
```

##Linear dimension reduction
```
merged_h5_seurat_filtered<-RunPCA(object=merged_h5_seurat_filtered)
```

Find dimension of dataset using elbowplot. This ranks principle components based on the percentage of variance explained by each one. From the plot below we see ~ first 15 principal components captured majority of the variation. 

```
ElbowPlot(merged_h5_seurat_filtered)
```

Image4


We use all 20 dimensions for cell clustering. 

```
merged_h5_seurat_filtered<-FindNeighbors(merged_h5_seurat_filtered, dims=1:20)
merged_h5_seurat_filtered<-FindClusters( merged_h5_seurat_filtered, resolution=0.5)
```

##Non-linear dimensional reduction (UMAP)
Use the same PCs as input to the clustering analysis.

merged_h5_seurat_filtered<-RunUMAP(merged_h5_seurat_filtered, dims=1:20)

DimPlot(merged_h5_seurat_filtered, reduction = "umap")
```

Image 5


Perform Seurat integration to correct for batch effects

Split the dataset into a list of four seurat objects (M1,M2, C1,C2), then normalize and identify variable features for each dataset independently.
```
obj.list<-SplitObject(merged_h5_seurat_filtered, split.by='ID')
obj.list
for(i in 1:length(obj.list)){
  obj.list[[i]]<-NormalizeData(object=obj.list[[i]])
  obj.list[[i]]<-FindVariableFeatures(object=obj.list[[i]])
}
```

Select features that are repeatedly variable across datasets for integration
```
features<-SelectIntegrationFeatures(object.list = obj.list)
```


Find integration anchors (CCA) to correct technical differences
```
anchors<-FindIntegrationAnchors(object.list=obj.list, anchor.features=features)
```

Integrate data
```
seurat.integrated<-IntegrateData(anchorset=anchors)
```

#specify that we will perform downstream analysis on integrated data, note that the original unmodified data still resides in the 'RNA' assay
```
DefaultAssay(seurat.integrated)<-"integrated"
```

Run the standard workflow for visualization and clustering. RunPCA by default compute and stores 50 PCs

```
seurat.integrated<-ScaleData(seurat.integrated)
seurat.integrated<-RunPCA(seurat.integrated)
seurat.integrated<-RunUMAP(seurat.integrated, dims=1:50)
seurat.integrated<-FindNeighbors(seurat.integrated,reduction="pca",dims=1:50)
seurat.integrated<-FindClusters(seurat.integrated, resolution=0.1)

```

p3<-DimPlot(seurat.integrated, reduction='umap', group.by = 'Type')

p4<-DimPlot(seurat.integrated, reduction='umap', group.by='ID',
            cols=c('red','green','blue','black'))


grid.arrange(p1,p2,p3,p4, ncol=2, nrow=2)


#change default assay to 'rna'
DefaultAssay(seurat.integrated)<-"RNA"


#find markers for cluster 6
cluster6.markers<-FindMarkers(seurat.integrated, ident.1=6, min.pct=0.5, only.pos=TRUE, logfc.threshold=0.25)

#find markers for every cluster compared to all remaining cells
ALL.markers<-FindAllMarkers(seurat.integrated, min.pct=0.5, only.pos=TRUE, logfc.threshold=0.25)

View(ALL.markers)

#Differential gene expression
#find markers distinguishing cluster 3 from cluster  6
cluster3_6.markers<-FindMarkers(seurat.integrated, ident.1=3, ident.2=6, min.pct=0.25)

#find all markers distinguishing cluster 3 from clusters 1 and 6
cluster3_1and6.markers<-FindMarkers(seurat.integrated, ident.1=3, ident.2=c(1,6), min.pct=0.25)

str(seurat.integrated)  

View(seurat.integrated@meta.data)

clusters<-DimPlot(seurat.integrated, reduction='umap', group.by='seurat_clusters', label=TRUE)
condition<-DimPlot(seurat.integrated, reduction='umap', group.by = 'Type')

clusters|condition


#this separates cells by conditions in this case treatment vs control, then comparing gene expressing
#in each cell type in this cluster vs all other clusters
markers_cluster6<-FindConservedMarkers(seurat.integrated,
                                       ident.1=1,
                                       grouping.var='Type')
#the pct.1 refers to percentage of cells in the cluster expressing that gene and pct.2 is the percentage
#of cells combined in all other clusters expressing that gene
head(markers_cluster6)

#let's visualize top features
FeaturePlot(seurat.integrated, features=c('ZNF385D'), min.cutoff='q10')

#rename cluster 6 ident
Idents(seurat.integrated)
seurat.integrated<-RenameIdents(seurat.integrated, '6'='X Cells')

DimPlot(seurat.integrated, reduction='umap', label=T)

#findMarkers between conditions
seurat.integrated$celltype.cnd<-paste0(seurat.integrated$seurat_clusters,'_', seurat.integrated$Type)
View(seurat.integrated@meta.data)
Idents(seurat.integrated)<-seurat.integrated$celltype.cnd

DimPlot(seurat.integrated, reduction='umap', label=T)

#find markers
response<-FindMarkers(seurat.integrated, ident.1='1_Mock', ident.2 = '1_Cov')


#genes up or down regulated in ident.1 vs ident.2
head(response)

#plotting conserved features vs De features between conditions
head(markers_cluster6)

FeaturePlot(seurat.integrated, features=c('CXCL2',"FGB",'MX1'), split.by='Type', min.cutoff = 'q10')


#show on umap the cell clusters by genes
FeaturePlot(seurat.integrated, features=c('INS',"GCG",'SST','VIM'), split.by='Type', cols=c('light gray','blue'))


