### Generate cell-type-specific networks for each cancer type or Pan-cancer immune cell

library(Seurat)
library(ACTIONet)
library(SCINET)
library(igraph)
library(dplyr)
library(tidyverse)

file_paths = read.table("./pathway_for_SinglecellDATA/Celltype_data_frame.txt",header = T,sep = "\t")    ### Seurat data path 
all_reference_networks <- c("ConsensusPathDB","HumanNet","Reactome","STRING_PPI")

for (i in 1:nrow(file_paths)) {
  
  cancer = file_paths$cancer[i]
  if(file_paths$processed[i]){
    next
  }
  
  seurat_object = readRDS(file_paths$path[i])
  DefaultAssay(seurat_object) <- "CCA"
  celltype_data.frame = data.frame(celltype = Idents(seurat_object),barcode = colnames(seurat_object))
  
  celltype_table = celltype_data.frame$celltype %>% table()
  celltype_data.frame <- celltype_data.frame[celltype_data.frame$celltype %in% names(celltype_table)[which(celltype_table>100)] ,]

  seurat_object <- subset(seurat_object,cells = rownames(celltype_data.frame))

  ace <- import.ace.from.Seurat(seurat_object)
  
  ace = reduce.ace(ace)
  ace = run.ACTIONet(ace)

  ace$clusters = celltype_data.frame$celltype[match(colnames(ace),celltype_data.frame$barcode)]
  
  ace = compute.cluster.feature.specificity(ace, ace$clusters, 'cluster_specificity_scores')
  
  for(j in 1:length(all_reference_networks)){
    
    ref_net <- readRDS(paste0("./reference/",all_reference_networks[j],".RDS"))
    network.list = run.SCINET.clusters(ace, 'cluster_specificity_scores_feature_specificity',G = ref_net)
    saveRDS(network.list,file =paste0("./SCINET_Network/",cancer,"_",all_reference_networks[j],".RDS"))

  }
  
  message(paste0(cancer,"   Finished"))
  
}






