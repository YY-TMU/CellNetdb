###  CCA  batch correction Demo 

library(stringr)
library(Seurat)
setwd("/data1/lizekun/sc_cancer/0.rawdata/neuroblastoma_GSE137804/h5")           #######  cellranger output 
fs=list.files(pattern = '.h5')

sceList = lapply(fs, function(x){
  #x=fs[1]
  print(x)
  a=Read10X_h5( x )
  
  p=str_split(x,'_',simplify = T)[,1]
  sce <- CreateSeuratObject( a ,project = p , min.cells = 3, min.features =200)
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce <- subset(sce, subset =nFeature_RNA > 200  & percent.mt < 25 & nCount_RNA > 800)
  
})
for (i in 1:length(sceList)) {
  sceList[[i]] <- NormalizeData(sceList[[i]], verbose = FALSE,normalization.method = "LogNormalize", scale.factor = 10000)
  sceList[[i]] <- FindVariableFeatures(sceList[[i]], selection.method = "vst", 
                                       nfeatures = 2000, verbose = FALSE)
}

#CCA
sce.anchors <- FindIntegrationAnchors(object.list = sceList, dims = 1:30,reduction  = "cca")
sce.CCA <- IntegrateData(anchorset = sce.anchors, dims = 1:30,new.assay.name = "CCA")
sce.CCA <- ScaleData(sce.CCA,vars.to.regress = c("nCount_RNA","percent.mt"))
sce.CCA <- RunPCA(sce.CCA,features = VariableFeatures(sce.CCA))

ElbowPlot(sce.CCA,ndims = 30)
sce.CCA <- FindNeighbors(sce.CCA, dims = 1:18)
sce.CCA <- FindClusters(sce.CCA, resolution = 0.8)
sce.CCA <- RunUMAP(sce.CCA, dims = 1:18)
DimPlot(sce.CCA,reduction = "umap", pt.size = 0.3, seed = 1234, label = TRUE)





