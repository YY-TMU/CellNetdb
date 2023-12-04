###  Generate CellChat data

library(tidyverse)
library(CellChat)
library(Seurat)
library(stringi)
library(stringr)
library(dplyr)

file_paths = read.table("./Cellchat_dataFrame.txt",header = T,sep = "\t")

for(i in 1:nrow(file_paths)){
  
  if(file_paths$processed[i]){
    next
  }
  
  seurat_object = readRDS(file_paths$path[i])
  data.input  <- seurat_object@assays$RNA@data
  identity = data.frame( group = Idents(seurat_object) , row.names = rownames(seurat_object@meta.data) )
  
  identity$group = identity$group %>% as.character() %>% as.factor()
  
  cellchat <- createCellChat(data.input)
  CellChatDB <- CellChat::CellChatDB.human
  cellchat@DB = CellChatDB
  
  cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
  cellchat <- setIdent(cellchat, ident.use = "labels")
  
  cellchat <- subsetData(cellchat)
  # future::plan("multiprocess", workers = 4)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, CellChat::PPI.human)
  
  cellchat <- computeCommunProb(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  saveRDS(cellchat,file = paste("./CellChat_result/",file_paths$cancer[i],".RDS",sep = ""))
  # cellchat <- updateCellChat(cellchat)
  cellchat@net$adj_prob <- cellchat@net$prob
  LR_list = list()
  for(j in 1:length(cellchat@netP$pathways)){
    
    the.path = cellchat@netP$pathways[j]
    
    tab = searchPair(signaling = the.path, pairLR.use = cellchat@LR$LRsig,key = "pathway_name", matching.exact = T, pair.only = T) #get L-R
    tab$score = 0
    for(net_name in tab$interaction_name){
      net = cellchat@net$prob[,, net_name]
      pvalue_net = cellchat@net$pval[,, net_name]
      adj_net <- net
      adj_net[pvalue_net > 0.05] <- 0
      
      tab$score[which(tab$interaction_name == net_name)] = adj_net %>% sum()
      cellchat@net$adj_prob[,, net_name] <- adj_net
    }
    tab = tab[which(tab$score != 0),]
    LR_list[[1 + length(LR_list)]] = tab
    
  }
  LR_list = do.call("rbind",LR_list)
  net_list = cellchat@net$adj_prob[,,LR_list$interaction_name]
  
  result_list = list(LR = LR_list, net = net_list)
  
  saveRDS(result_list , file = paste0("./data/new_cellchat/",file_paths$cancer[i],".RDS"))
  
  file_paths$processed[i] = TRUE
  file_paths$state[i] = "*"
  write.table(file_paths,file = "./Cellchat_dataFrame.txt",col.names = T,row.names = F,sep = "\t",quote = F)
  
  print( paste(file_paths$cancer[i] , "is finished") )
  

}




