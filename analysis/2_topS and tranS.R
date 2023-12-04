##  Topological and transcriptional specificity score 


############
####  topS 
############

library(igraph)
library(tidyverse)
library(tidygraph)
library(Matrix)


topo.spec <- function(G, sample_no = 10000) {
  d = igraph::degree(G)
  w = igraph::strength(G)  
  
  edge_weights = E(G)$weight
  N = length(edge_weights)
  
  UD = setdiff(sort(unique(d)), 0)
  mean_rand_stat = array(0, max(UD))  
  sd_rand_stat = array(0, max(UD))  
  for(deg in UD){
    rand_samples = sample(edge_weights, sample_no*deg, replace=TRUE);
    start_idx = seq(from=1, to = length(rand_samples), by = deg);
    rand_sums = sapply(start_idx, function(i) {sum(rand_samples[i:(i+deg-1)])})
    
    mean_rand_stat[deg] = mean(rand_sums);
    sd_rand_stat[deg] = sd(rand_sums);    
  }
  
  specificity = sapply( 1:length(w), function(i) {if(w[i] == 0) {x = -Inf;} else {x = (w[i] - mean_rand_stat[d[i]]) / sd_rand_stat[d[i]];}; return(x) } )
  
  specificity[is.na(specificity)] = -Inf
  
  names(specificity) = V(G)$name
  return(specificity)
}


used_cancer <- c("B_cells",
                 "CD4_Tcells",
                 "CD8_Tcells",
                 "Myeloid_cells")


cancer_net_path <- "/home/liugerui/pro/FGN/SCINET_Network(Imm)/"


for(i in 1:length(used_cancer)){
  
  this.cancer <- used_cancer[i]
  cancer_net_list <- readRDS(paste0(cancer_net_path,this.cancer,"_STRING.RDS"))
  
  cell_types <- names(cancer_net_list)
  
  topS_list <- list()
  
  for(j in 1:length(cell_types)){
    this.cell.type <- cell_types[j]
    this.net <- cancer_net_list[[this.cell.type]]
    
    test_out <- topo.spec(this.net,sample_no = 10000)
    
    test_out_df <- data.frame(spec=test_out,gene=names(test_out))
    
    topS_list[[this.cell.type]] <- test_out_df
    
  }
  
  saveRDS(topS_list,paste0("/data1/lizekun/CellNetdb/add_analysis/06.topS/TIMEs/",this.cancer,".rds"))
  
  message(i)
}



############
####  tranS
############
####  use the exp profile after  imputation and batch correction

library(Seurat)
library(ACTIONet)
library(SCINET)
library(igraph)
library(dplyr)
library(tidyverse)


transc.spec.celltype <- function(Arch.imputed.profile, Labels,this.celltype) {
  
  gene_names = rownames(Arch.imputed.profile)
  
  if(!is.factor(Labels)){
    Labels = factor(Labels, levels = sort(unique(Labels)))
  }
  
  UL = levels(Labels)
  
  celltype <- this.celltype
  
  transcriptional.gene.specificity.Z = array(0, dim = c(length(gene_names), 1))
  rownames(transcriptional.gene.specificity.Z) = gene_names
  colnames(transcriptional.gene.specificity.Z) = celltype
  
  
  R.utils::printf('%s\n', celltype)
  
  Group = which(Labels %in% celltype)
  Null = which(!(Labels %in% celltype))
  
  A_group = Arch.imputed.profile[, Group] %>% as.matrix()
  A_null = Arch.imputed.profile[, Null]   %>% as.matrix()
  
  delta_mean = rowMeans(A_group) - rowMeans(A_null);
  sigma1_sq = apply(A_group, 1, var)
  sigma2_sq = apply(A_null, 1, var)
  sigma_pooled = sqrt( (sigma1_sq / length(Group)) + (sigma2_sq / length(Null)) )
  z.stat = delta_mean / sigma_pooled
  
  transcriptional.gene.specificity.Z[, celltype] = z.stat
  
  
  return(transcriptional.gene.specificity.Z)
}


impute.genes.using.archetypes <- function(ace, genes, features_use = NULL) {
  
  Z <- rowMaps(ace)[["unified_feature_profile"]][features_use, ]
  
  H <- Matrix::t(colMaps(ace)[["H_unified"]])
  
  expression_imputed <- Z %*% H
  
  return(expression_imputed)
}




this.cancer <- "Myeloid_cells"

seurat_object <- readRDS(paste0("/data1/lizekun/CellNetdb/Imm/",this.cancer,"/",this.cancer,".rds"))
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
##### saveRDS(ace,"/data1/lizekun/CellNetdb/add_analysis/07.tranS/TIME/Myeloid_cells/Myeloid_cells.rds")

metadata=colData(ace) %>%data.frame()

NETgene_list <- readRDS("/data1/lizekun/CellNetdb/Imm/data/immnet_geneList.rds")
cancer_net <- names(NETgene_list)[str_detect(string = names(NETgene_list),pattern = this.cancer)]
cancer_net <- cancer_net[grep("STRING",cancer_net)]

cancer_net_list <- NETgene_list[cancer_net]

lapply(cancer_net_list,length) %>% unlist()


for(i in 1:length(cancer_net_list)){
  
  this.cell_type <- names(cancer_net_list)[i]  %>% gsub(pattern = "Myeloid_cells_STRING_",replacement = "")
  
  net_gene <- cancer_net_list[[i]] %>% names()
  
  expression_imputed <- impute.genes.using.archetypes(ace = ace,features_use = net_gene)
  
  tranS <- transc.spec.celltype(Arch.imputed.profile = expression_imputed ,Labels = as.character(metadata$cluster),this.celltype=this.cell_type  )
  
  saveRDS(tranS,paste0("/data1/lizekun/CellNetdb/add_analysis/07.tranS/2.tranS_output/TIME/Myeloid_cells/tranS/",
                       this.cell_type,".rds"))
  message(i)
  
}





