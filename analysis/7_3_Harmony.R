####   Different  batch correction methonds 

library(Seurat)
library(tidyverse)
library(harmony)
library(anndata)
library(reticulate)
library(SeuratDisk)
use_python("~/apps/anaconda3/envs/scvi2/bin/python")

library(sceasy)
sc <- import("scanpy", convert = FALSE)

args <- commandArgs(trailingOnly = TRUE)
cancer_type <- args[1]
print(cancer_type)

tumor_path <- "/data1/lizekun/sc_cancer/2.re_cluster/total_cells/RDS/"


cancer_data <- readRDS(paste0(solid_tumor_path,cancer_type,".rds"))
DefaultAssay(cancer_data) <- "RNA"
cancer_data[["percent.mt"]] <- PercentageFeatureSet(cancer_data, pattern = "^MT-")

over_metadata <-  read.csv(paste0("/data1/lizekun/CellNetdb/add_analysis/scvi/00.adj_rds/cancer_metadata/",cancer_type,".csv"),row.names = 1)
identical(colnames(cancer_data),rownames(over_metadata)) %>% message()

cancer_data <- AddMetaData(cancer_data,metadata = over_metadata$Cell_TypeID %>% as.character(),col.name = "Cell_TypeID")
cancer_data <- AddMetaData(cancer_data,metadata = over_metadata$BatchID %>% as.character(),col.name = "BatchID")

if("DonorID" %in% colnames(over_metadata)){
  cancer_data <- AddMetaData(cancer_data,metadata = over_metadata$DonorID %>% as.character(),col.name = "DonorID")
}

# metadata <- cancer_data@meta.data

cancer_data <- NormalizeData(cancer_data, verbose = FALSE,normalization.method = "LogNormalize", scale.factor = 10000)
cancer_data <- FindVariableFeatures(cancer_data,selection.method = "vst",verbose = F,nfeatures = 2000)
cancer_data <- ScaleData(cancer_data,vars.to.regress = c("nCount_RNA","percent.mt"))

cancer_data <- RunPCA(cancer_data,features = VariableFeatures(cancer_data),verbose = F)

# DimPlot(cancer_data,reduction = "pca", group.by = "BatchID",pt.size = 0.3, seed = 1234, label = TRUE)
# DimPlot(cancer_data,reduction = "pca", group.by = "Cell_TypeID",pt.size = 0.3, seed = 1234, label = TRUE)

options(repr.plot.height=2.5,repr.plot.width = 6)
cancer_data <- RunHarmony(object = cancer_data,
                          group.by.vars = "BatchID",
                          reduction = "pca",
                          plot_convergence = TRUE)


cancer_harmony_emb_df <- cancer_data@reductions$harmony@cell.embeddings

# DimPlot(cancer_data,reduction = "harmony", group.by = "BatchID",pt.size = 0.3, seed = 1234, label = TRUE)
# DimPlot(cancer_data,reduction = "harmony", group.by = "Cell_TypeID",pt.size = 0.3, seed = 1234, label = TRUE)

cancer_data <- FindNeighbors(cancer_data,reduction = "harmony",dims = 1:30)
cancer_data <- FindClusters(cancer_data)
cancer_data <- RunUMAP(cancer_data,reduction = "harmony",dims = 1:30)
# DimPlot(cancer_data,reduction = "umap", group.by = "Cell_TypeID",pt.size = 0.3, seed = 1234, label = TRUE)
# DimPlot(cancer_data,reduction = "umap", group.by = "BatchID",pt.size = 0.3, seed = 1234, label = TRUE)

cancer_umap_emb_df <- cancer_data@reductions$umap@cell.embeddings

write.csv(cancer_harmony_emb_df,paste0("/data1/lizekun/CellNetdb/add_analysis/harmony/BatchID/01.cancer_emd_data/01.harmony/harmony/",cancer_type,".csv"),quote = F)
write.csv(cancer_umap_emb_df,paste0("/data1/lizekun/CellNetdb/add_analysis/harmony/BatchID/01.cancer_emd_data/01.harmony/umap/",cancer_type,".csv"),quote = F)


adata <- convertFormat(cancer_data, from="seurat", to="anndata", drop_single_values=FALSE)
print(adata)
write_h5ad(adata,paste0("/data1/lizekun/CellNetdb/add_analysis/harmony/BatchID/02.h5ad/",cancer_type,".h5ad"))



