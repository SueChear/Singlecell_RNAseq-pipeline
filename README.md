# Single cell RNA-Seq data analysis pipeline: Reprocessing data from a study investigating susceptibility of multiple types of pancreatic islet cells to SARS-CoV-2.
Tang, Xuming et al. “SARS-CoV-2 Infection Induces Beta Cell Transdifferentiation.” Cell Metabolism, 19 May. 2021, [doi:10.1016/j.cmet.2021.05.015
]([url](https://www.sciencedirect.com/science/article/pii/S1550413121002321?via%3Dihub))

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
## QC and selecting cells for further analysis

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

<img width="950" alt="1_qcviolinplots" src="https://user-images.githubusercontent.com/117556524/215320903-d68970e4-0031-4dc8-87b3-daaa1f4e3bdc.PNG">

FeatureScatter is typically used to visualize feature-feature relationships.Quality data should follow the straight line trend, the scattered points on lower
right quadrant indicate some low number of genes that have been sequenced again and again

```
FeatureScatter(merged_h5_seurat, feature1="nCount_RNA", feature2="nFeature_RNA")+
  geom_smooth(method='lm')
```


<img width="951" alt="2_featurefeatureplot" src="https://user-images.githubusercontent.com/117556524/215320967-82bdcd47-9008-431d-b4ca-d0e1afa6390a.PNG">


We filter cells that have unique feature counts <500 or >6000 and unique RNA counts <1000 or >60000 and cells that have >15% mitochondrial counts. These thresholds are taken from the published paper.

```
merged_h5_seurat_filtered<-subset(merged_h5_seurat,subset=nFeature_RNA >500 & nFeature_RNA<=6000 & nCount_RNA >1000 & nCount_RNA <= 60000 &
                                    mitoPercent< 15)

```

After removing unwanted cells, we normalize data in order to compare gene expression across multiple cells. Wwe divide gene expression measurement in each cell by the total expression, then multiply it by a scaling factor followed by a  log transformation. Default normalization simply scales the counts by total counts in each cell, multiplies by 10,000 and then log transforms.

```
merged_h5_seurat_filtered<-NormalizeData(merged_h5_seurat_filtered)
```

## Identification of highly variable features (feature selection)
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

<img width="953" alt="3_highlyvariablefeatures" src="https://user-images.githubusercontent.com/117556524/215320994-6572cac7-e6dd-43b6-beae-6d9f8aa4f1fd.PNG">



## Scaling the data
We scale the data prior to PCA.
```
merged_h5_seurat_filtered<-ScaleData(object=merged_h5_seurat_filtered)
```

## Linear dimension reduction
```
merged_h5_seurat_filtered<-RunPCA(object=merged_h5_seurat_filtered)
```

Find dimension of dataset using elbowplot. This ranks principle components based on the percentage of variance explained by each one. From the plot below we see ~ first 15 principal components captured majority of the variation. 

```
ElbowPlot(merged_h5_seurat_filtered)
```

<img width="957" alt="4_elbowplot" src="https://user-images.githubusercontent.com/117556524/215321017-93d6dd71-444b-4421-b036-efc16132aa8e.PNG">


We use all 20 dimensions for cell clustering. 

```
merged_h5_seurat_filtered<-FindNeighbors(merged_h5_seurat_filtered, dims=1:20)
merged_h5_seurat_filtered<-FindClusters( merged_h5_seurat_filtered, resolution=0.5)
```

## Non-linear dimensional reduction (UMAP)
Use the same PCs as input to the clustering analysis.
```
merged_h5_seurat_filtered<-RunUMAP(merged_h5_seurat_filtered, dims=1:20)

DimPlot(merged_h5_seurat_filtered, reduction = "umap")
```

<img width="952" alt="5_umap" src="https://user-images.githubusercontent.com/117556524/215321039-2b4b2eb0-3e30-4fb3-afa9-6863ec87eda5.PNG">



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

Sspecify that we will perform downstream analysis on integrated data, note that the original unmodified data still resides in the 'RNA' assay
```
DefaultAssay(seurat.integrated)<-"integrated"
```

Run the standard workflow for visualization and clustering. RunPCA by default computes and stores 50 PCs.

```
seurat.integrated<-ScaleData(seurat.integrated)
seurat.integrated<-RunPCA(seurat.integrated)
seurat.integrated<-RunUMAP(seurat.integrated, dims=1:50)
seurat.integrated<-FindNeighbors(seurat.integrated,reduction="pca",dims=1:50)
seurat.integrated<-FindClusters(seurat.integrated, resolution=0.1)
```
Comparison of plots before and after integration. Plots in the second row showed cells from all samples integrated with each other compared to plots in the first row (before integration).
```
p1<-DimPlot(merged_h5_seurat_filtered, reduction='umap', group.by = 'Type')

p2<-DimPlot(merged_h5_seurat_filtered, reduction='umap', group.by='ID',
            cols=c('red','green','blue','black'))
p3<-DimPlot(seurat.integrated, reduction='umap', group.by = 'Type')

p4<-DimPlot(seurat.integrated, reduction='umap', group.by='ID',
            cols=c('red','green','blue','black'))

p5 <- DimPlot(seurat.integrated, reduction = "umap", label = TRUE, repel = TRUE)


grid.arrange(p1,p2,p3,p4,p5, ncol=2, nrow=3)
```

<img width="948" alt="6_integration" src="https://user-images.githubusercontent.com/117556524/215321072-1fa4fa8d-2224-41a9-9730-b46943a52047.PNG">

To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
```
DimPlot(seurat.integrated, reduction = "umap", split.by = "Type")
```

<img width="952" alt="7_bycondition" src="https://user-images.githubusercontent.com/117556524/215321083-fdd185fd-a804-4e93-90fa-82a83f6dc998.PNG">

For performing differential expression after integration, we switch back to the original data
Change default assay to 'rna'
```
DefaultAssay(seurat.integrated)<-"RNA"
```

Here, marker genes from multiple cell types were refered from the published paper. For a quick screening of cell type in various clusters, we compare the plots of markers genes and clusters.
Marker genes: PRSS1+ acinar cells, GCG+ alpha cells, INS+ beta cells, KRT19+ ductal cells, COL1A1+ fibroblasts, PPY+ PP cells, SST+ delta cells, PECAM1+ endothelial cells, and LAPTM5+ immune cells. Plot shows acinar cells corresponding to cluster 3, alpha cells to cluster 2, beta cells to cluster 1, ductal cells to cluster 4,delta cells to cluster 6 and etc.
```
a1<-FeaturePlot(seurat.integrated, features=c('PRSS1','GCG','INS','KRT19','COL1A1','PPY','SST','PECAM1','LAPTM5'), cols=c('light gray', 'blue'))
a2<-DimPlot(immune.combined, reduction = "umap", label=TRUE)

a1|a2
```

<img width="950" alt="8_clusterandmarkergenes" src="https://user-images.githubusercontent.com/117556524/215321099-c4e37080-708f-4642-bf1c-407286ea1bdb.PNG">


To find marker genes for cluster 1, we use FindMarkers function.
```
cluster1.markers<-FindMarkers(seurat.integrated, ident.1=1, min.pct=0.5, only.pos=TRUE, logfc.threshold=0.25)
head(cluster1.markers)
```
```
p_val avg_log2FC pct.1 pct.2 p_val_adj
ERO1B       0  1.2805018 0.744 0.347         0
PTPRN       0  0.6095381 0.808 0.443         0
KIF1A       0  0.6230766 0.649 0.331         0
ZNF385D     0  2.0605839 0.540 0.047         0
UCHL1       0  1.5207275 0.750 0.371         0
HADH        0  1.1885251 0.657 0.206         0
```
This analysis compares each cluster against all others and outputs the genes that are differentially expressed/present using the FindAllMarkers() function.
```
ALL.markers<-FindAllMarkers(seurat.integrated, min.pct=0.5, only.pos=TRUE, logfc.threshold=0.25)
```

## Differential gene expression
To find markers distinguishing cluster 3 from cluster  1.

```
cluster3_1.markers<-FindMarkers(seurat.integrated, ident.1=3, ident.2=1, min.pct=0.25)
head(cluster3_1.markers)
```
```
     p_val avg_log2FC pct.1 pct.2 p_val_adj
DHRS3      0  0.6663316 0.392 0.035         0
CTRC       0  3.4762435 0.935 0.437         0
CELA2A     0  3.5000203 0.889 0.496         0
CELA2B     0  1.6969351 0.523 0.065         0
CELA3B     0  3.5588546 0.917 0.546         0
CELA3A     0  3.7938090 0.967 0.781         0
```


FindConservedMarkers separates cells by conditions in this case Mock vs Cov, then comparing gene expressing in each cell type in this cluster vs all other clusters
Pct.1 refers to percentage of cells in the cluster expressing that gene and pct.2 is the percentage of cells combined in all other clusters expressing that gene
```
markers_cluster6<-FindConservedMarkers(seurat.integrated,ident.1=1,grouping.var='Type')
head(markers_cluster6)
```
```
Cov_p_val Cov_avg_log2FC Cov_pct.1 Cov_pct.2 Cov_p_val_adj    Mock_p_val Mock_avg_log2FC Mock_pct.1 Mock_pct.2 Mock_p_val_adj      max_pval
SAMD11           0      0.3434806     0.249     0.017             0  0.000000e+00       0.5562374      0.323      0.022   0.000000e+00  0.000000e+00
C1orf127         0      0.4813396     0.289     0.016             0  0.000000e+00       0.3988611      0.252      0.013   0.000000e+00  0.000000e+00
S100A11          0     -1.6016597     0.736     0.821             0 7.507749e-262      -1.3442169      0.779      0.840  2.010876e-257 7.507749e-262
RGS16            0      0.5597420     0.289     0.042             0  0.000000e+00       0.5932222      0.282      0.041   0.000000e+00  0.000000e+00
SUSD4            0      0.3681491     0.333     0.059             0  0.000000e+00       0.3983245      0.323      0.063   0.000000e+00  0.000000e+00
ERO1B            0      1.2688898     0.736     0.340             0  0.000000e+00       1.2876819      0.753      0.353   0.000000e+00  0.000000e+00
         minimump_p_val
SAMD11                0
C1orf127              0
S100A11               0
RGS16                 0
SUSD4                 0
ERO1B                 0
```

Rename cluster to identified cell types based on marker genes expression.

```
seurat.labeled <- RenameIdents(seurat.integrated, `3` = "acinar cells", `2` = "alpha cells", `1` = "beta cells",
                               `4` = "ductal cells", `6` = "delta cells")
DimPlot(seurat.labeled, label = TRUE)

```

<img width="951" alt="9_renameclusters" src="https://user-images.githubusercontent.com/117556524/215321112-88d45ecf-edbc-4512-9d8b-d8daabafb52e.PNG">


## Identify differential expressed genes across conditions

Create a column named celltype.cnd which merges cluster name and condition, and assign celltype.cnd as new cell identity.
```
seurat.integrated$celltype.cnd<-paste0(seurat.integrated$seurat_clusters,'_', seurat.integrated$Type)

Idents(seurat.integrated)<-seurat.integrated$celltype.cnd

DimPlot(seurat.integrated, reduction='umap', label=T)
```

Find markers between mock and covid infected cells in cluster 1 beta cells. Genes showed below are up/down regulated in Mock samples relative to Cov samples.

```
response<-FindMarkers(seurat.integrated, ident.1='1_Cov', ident.2 = '1_Mock')
head(response)
```

```
 p_val avg_log2FC pct.1 pct.2    p_val_adj
FGB     9.350848e-83 -0.2690794 0.661 0.129 1.870170e-79
CXCL2   8.210034e-64 -0.2576961 0.774 0.316 1.642007e-60
MX1     4.194869e-61 -0.2588754 0.991 0.556 8.389739e-58
SPON2   7.968229e-56  0.3101120 0.656 0.205 1.593646e-52
PRG4    9.478011e-43 -0.3641632 0.674 0.231 1.895602e-39
GADD45B 2.679879e-12 -0.2705537 0.835 0.513 5.359758e-09
```

Plot differentially expressed genes between conditions
```
FeaturePlot(seurat.integrated, features=c('SPP1',"PRG4",'HMOX1'), split.by='Type', min.cutoff = 'q10', cols=c("light gray","red"))
```
<img width="946" alt="10_DEG" src="https://user-images.githubusercontent.com/117556524/215321138-78ab5d50-4c68-4a24-8624-3fc3ded67260.PNG">




