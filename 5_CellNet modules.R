###  CellNet demo example

library(tidyverse)
library(networkD3)
library(CellChat)
input <- list()
input$ctnet_cancer <- "Adrenal neuroblastoma"
input$ctnet_reference_net <- "STRING"
input$ctnet_celltype <- "Neuroendocrine cells (malignant)"
input$ctnet_gene <- "MDK"


ct_namechange = read.table("./data/ct_namechange.txt",sep = "\t",header = T)

text_Optimize <- function(st_set , the.type,input = NA){

  
  if(the.type == "input_cancer"){
    st_set = ct_namechange$old_cancer[match(st_set,ct_namechange$new_cancer)]
  }else if(the.type == "input_celltype"){
    cancer = input$ctnet_cancer
    tab = ct_namechange[which(ct_namechange$new_cancer == cancer),]
    st_set = tab$old_celltype[match(st_set,tab$new_celltype)]
    
  }else if(the.type == "output_cancer"){
    st_set = ct_namechange$new_cancer[match(st_set,ct_namechange$old_cancer)]
  }else if(the.type == "output_celltype"){
    
    
    cancer = input$ctnet_cancer
    tab = ct_namechange[which(ct_namechange$new_cancer == cancer),]
    st_set = tab$new_celltype[match(st_set,tab$old_celltype)]
    
  }else if(the.type == "change_raw"){
    namechange =read.table("./data/name_change.txt",sep = "\t",header = T)
    st_set = namechange$cell_name[match(st_set,namechange$cell_subset)]
  }else if(the.type == "change_old"){
    namechange =read.table("./data/name_change.txt",sep = "\t",header = T)
    st_set = namechange$cell_subset[match(st_set,namechange$cell_name)]
  }

  return(st_set)
}

generate_Net <- function(input){
  
  cancer = input$ctnet_cancer %>% text_Optimize(the.type = "input_cancer")
  celltype = input$ctnet_celltype %>% text_Optimize(the.type = "input_celltype",input)
  the.gene = input$ctnet_gene
  Reference_net = input$ctnet_reference_net
  
  net_path = paste("./data/new_Network/",cancer,"_",Reference_net,"_",celltype,".RDS",sep = "")
  net = readRDS(net_path)
  net = net$net
  
  geneset = c(
    rownames(net)[which(net[,the.gene] != 0)],
    colnames(net)[which(net[the.gene,] != 0)],the.gene
  ) %>% unique()
  
  bg = rownames(net)
  net = net[geneset,geneset]
  adj_net = net + t(net)
  result_net=adj_net
  
  Node_size = rowSums(adj_net)
  
  net = reshape2::melt(net)
  net = net[which(net$value!=0),]
  
  src <- net$Var1 %>% as.character()
  target <- net$Var2 %>% as.character()
  networkData <- data.frame(src, target)
  
  MisNodes <- data.frame(name=unique(c(src, target)),
                         group= 1,size = 1)
  MisLinks <- data.frame(source = (match(networkData$src, MisNodes$name)-1),
                         target = (match(networkData$target, MisNodes$name)-1),
                         value = 1)
  MisNodes$group[which(MisNodes$name %in% the.gene)] = 2
  MisNodes$size = Node_size[MisNodes$name] %>% as.integer()
  
  result = list(MisNodes = MisNodes,
                MisLinks = MisLinks,
                bg = bg,
                net =result_net)
  
  
}

the.Net = generate_Net(input)


## plot the "MDK" network
ColourScale <- 'd3.scaleOrdinal()
            .domain(["1", "2"])
           .range(["#ABC2E1","#0868A3"]);'

forceNetwork(Links = the.Net$MisLinks, Nodes = the.Net$MisNodes,
             Source = "source", Target = "target",
             Value = "value", NodeID = "name",Nodesize = "size",fontSize = 20,
             Group = "group",legend = F,arrows = F,opacity = 1,colourScale = JS(ColourScale),zoom = TRUE,fontFamily = "Arial")


###   CellNet Mutation module
tab_mut = readRDS(paste0("./data/final_mutation_Cancer/",input$ctnet_cancer,"_",input$ctnet_reference_net,"_",input$ctnet_celltype,".RDS"))
tab_mut <- unique(tab_mut)
tab_mut = tab_mut[which(tab_mut$Symble %in% the.Net$MisNodes$name),]
tab_mut = tab_mut[order(tab_mut$Symble),]
tab_mut$HGVSG = paste("chr",tab_mut$HGVSG,sep = "")
tab_mut$`Mutation genome position` = paste("chr",tab_mut$`Mutation genome position`,sep = "")


###   CellNet GO module
GO.list = readRDS("./data/GO_reference.RDS")                             ##  each GO Term genelist
GO_gene_set = readRDS("./data/GO_gene_set.RDS")                          ##  BP,CC,MF genelist

fun_for_GO <- function(total_geneset,GO.list,total_universe){
  
  overlap = data.frame(
    
    GO_id = names(GO.list),
    GO_name = "",
    overlap_for_inter = 0,
    overlap_for_uniANDref = 0,
    overlap_for_uniDIFFref = 0,
    genesetLEN = 0,
    pvalue = 1,
    genes = "",
    GeneRatio = "",
    BgRatio = "",
    fold_enrichment = 0,
    type = ""
    
  )
  
  geneset_List = list(
    `GO:BP` = intersect(GO_gene_set$`GO:BP` , total_geneset) ,
    `GO:CC` = intersect(GO_gene_set$`GO:CC` , total_geneset) ,
    `GO:MF` = intersect(GO_gene_set$`GO:MF` , total_geneset)
  )
  
  universe_List = list(
    `GO:BP` = intersect(GO_gene_set$`GO:BP` , total_universe) ,
    `GO:CC` = intersect(GO_gene_set$`GO:CC` , total_universe) ,
    `GO:MF` = intersect(GO_gene_set$`GO:MF` , total_universe)
  )
  
  for(i in 1:nrow(overlap)){
    
    GO_id = overlap$GO_id[i]
    the.type = GO.list[[GO_id]]$type
    
    universe = universe_List[[the.type]]
    geneset = geneset_List[[the.type]]
    
    
    overlap_for_inter = intersect(GO.list[[GO_id]]$geneset , geneset)
    if(length(overlap_for_inter)<2 ){
      next
    }
    
    overlap_for_uniANDref = intersect(GO.list[[GO_id]]$geneset,universe )
    overlap_for_uniDIFFref = setdiff(universe , GO.list[[GO_id]]$geneset)
    
    overlap$overlap_for_inter[i] = length(overlap_for_inter)
    overlap$overlap_for_uniANDref[i] = overlap_for_uniANDref %>% length()
    overlap$overlap_for_uniDIFFref[i] = overlap_for_uniDIFFref %>% length()
    overlap$genesetLEN[i] = length(geneset)
    
    overlap$genes[i] = paste0(overlap_for_inter,collapse = ",")
    overlap$GeneRatio[i] = paste(length(overlap_for_inter),"/", length(geneset) ,sep = " ")
    overlap$BgRatio[i] = paste(length(overlap_for_uniANDref),"/", length(universe) ,sep = " ")
    
    overlap$fold_enrichment[i] = ( length(overlap_for_inter)/length(geneset) ) / ( length(overlap_for_uniANDref)/length(universe) )
    
    overlap$pvalue[i] = phyper(
      
      overlap$overlap_for_inter[i] - 1,
      overlap$overlap_for_uniANDref[i],
      overlap$overlap_for_uniDIFFref[i],
      overlap$genesetLEN[i],
      lower.tail = F
      
    )
    
    overlap$GO_name[i] = GO.list[[GO_id]]$name
    overlap$type[i] = GO.list[[GO_id]]$type
    
  }
  
  overlap <- overlap[overlap$overlap_for_inter >0,]
  
  BP_overlap <- overlap[overlap$type %in% "GO:BP" ,]
  CC_overlap <- overlap[overlap$type %in% "GO:CC" ,]
  MF_overlap <- overlap[overlap$type %in% "GO:MF" ,]
  
  
  BP_overlap$FDR <- p.adjust(BP_overlap$pvalue,method = "fdr")
  CC_overlap$FDR <- p.adjust(CC_overlap$pvalue,method = "fdr")
  MF_overlap$FDR <- p.adjust(MF_overlap$pvalue,method = "fdr")
  
  overlap_adj <- rbind(BP_overlap,CC_overlap,MF_overlap)
  
  
  overlap_adj = overlap_adj[which(overlap_adj$pvalue<0.01),]
  
  return(overlap_adj)
  
}

Integration_GO <- function(input_tab){

  input_tab = input_tab[,c(1,2,9,10,11,12,7,13)]
  colnames(input_tab) = c(
    "GO term ID",
    "GO term name",
    "GeneRatio",
    "BgRatio",
    "Fold enrichment",
    "Type",
    "P-value",
    "FDR"
  )
  input_tab$`GO term name` = gsub(input_tab$`GO term name`,pattern = "GOBP_|GOMF|GOCC_",replacement = "")
  input_tab$Type = factor(input_tab$Type,levels = c("GO:BP","GO:MF","GO:CC"))
  
  input_tab = input_tab[order(input_tab$Type,input_tab$`P-value`),]
  input_tab$`P-value` = sprintf("%0.2e", input_tab$`P-value`)
  input_tab$FDR = sprintf("%0.2e", input_tab$FDR)
  input_tab$`Fold enrichment` = round(input_tab$`Fold enrichment`,digits = 1)
  
  
  return(input_tab)
  
}

tab_go = fun_for_GO(the.Net$MisNodes$name,GO.list,the.Net$bg)
tab_go = Integration_GO(tab_go)


###   CellNet Disease module 
disease.ls = readRDS("./data/disease.RDS")                             ###  DisGeNET disease genelist
Disease_gene_set = readRDS("./data/Disease_gene_set.RDS")              ### all DisGeNET disease genes

fun_for_Disease <- function(total_geneset,disease.ls,total_universe){
  
  
  universe = intersect( total_universe , Disease_gene_set )
  geneset = intersect( total_geneset , Disease_gene_set )
  
  overlap = data.frame(
    
    disease = names(disease.ls),
    overlap_for_inter = 0,
    overlap_for_uniANDref = 0,
    overlap_for_uniDIFFref = 0,
    genesetLEN = 0,
    pvalue = 1,
    genes = "",
    GeneRatio = "",
    BgRatio = "",
    fold_enrichment = 0
    
  )
  
  for(i in 1:nrow(overlap)){
    
    disease_name = overlap$disease[i]
    
    overlap_for_inter = intersect(disease.ls[[disease_name]] , geneset)
    if(length(overlap_for_inter)<2 ){
      next
    }
    
    overlap_for_uniANDref = intersect(disease.ls[[disease_name]],universe )
    overlap_for_uniDIFFref = setdiff(universe , disease.ls[[disease_name]])
    
    overlap$overlap_for_inter[i] = length(overlap_for_inter)
    overlap$overlap_for_uniANDref[i] = overlap_for_uniANDref %>% length()
    overlap$overlap_for_uniDIFFref[i] = overlap_for_uniDIFFref %>% length()
    overlap$genesetLEN[i] = length(geneset)
    
    overlap$genes[i] = paste0(overlap_for_inter,collapse = ",")
    overlap$GeneRatio[i] = paste(length(overlap_for_inter),"/", length(geneset) ,sep = " ")
    overlap$BgRatio[i] = paste(length(overlap_for_uniANDref),"/", length(universe) ,sep = " ")
    
    overlap$fold_enrichment[i] = ( length(overlap_for_inter)/length(geneset) ) / ( length(overlap_for_uniANDref)/length(universe) )
    
    overlap$pvalue[i] = phyper(
      
      overlap$overlap_for_inter[i] - 1,
      overlap$overlap_for_uniANDref[i],
      overlap$overlap_for_uniDIFFref[i],
      overlap$genesetLEN[i],
      lower.tail = F
      
    )
    
  }
  
  overlap <- overlap[overlap$overlap_for_inter >0,]
  overlap$FDR <- p.adjust(overlap$pvalue,method = "fdr")
  overlap = overlap[which(overlap$pvalue<0.01),]
  
  return(overlap)
  
}
Integration_Disease <- function(input_tab){

  input_tab = input_tab[,c(1,8,9,10,6,11)]
  colnames(input_tab) = c(
    "Disease",
    "GeneRatio",
    "BgRatio",
    "Fold enrichment",
    "P-value",
    "FDR"
  )
  
  options(scipen = 1)
  
  input_tab = input_tab[order(input_tab$`P-value`),]
  input_tab$`P-value` = sprintf("%0.2e", input_tab$`P-value`)
  input_tab$FDR = sprintf("%0.2e", input_tab$FDR)
  input_tab$`Fold enrichment` = round(input_tab$`Fold enrichment`,digits = 1)
  
  return(input_tab)
  
}

tab_disease = fun_for_Disease(the.Net$MisNodes$name,disease.ls,the.Net$bg)
tab_disease = tab_disease %>% Integration_Disease()


###   CellNet Communication module 
generate_Cellchat_data <- function(input,geneset){
  
  cellchat_data = readRDS(paste("./data/new_cellchat/",input$ctnet_cancer %>% text_Optimize(the.type = "input_cancer"),".RDS",sep = ""))
  LR = cellchat_data$LR
  net = cellchat_data$net
  
  used_LR = c()
  for(i in 1:nrow(LR)){
    
    intersect_gene = LR$interaction_name[i] %>% str_split(pattern = "_") %>% unlist()
    
    if(intersect(intersect_gene,geneset) %>% length() >0){
      used_LR = c( used_LR , LR$interaction_name[i] )
    }
  }
  
  used_LR = LR[which(LR$interaction_name %in% used_LR),]
  used_net = net[,,used_LR$interaction_name]
  
  result = list(LR = used_LR,
                net = used_net)
  return(result)
  
}
Integration_Cellchat <- function(input_tab,geneset){
  
  input_tab = input_tab$LR
  colnames(input_tab) = c(
    "Ligand-receptor pair",
    "Gene in network",
    "Ligand",
    "Receptor",
    "Communication score"
  )
  if(nrow(input_tab) == 0){
    return(input_tab)
  }
  
  input_tab = input_tab[order(input_tab$`Communication score`,decreasing = T),]
  for(i in 1:nrow(input_tab)){
    
    genes = paste(input_tab$Ligand[i],input_tab$Receptor[i],sep = "_") %>% str_split(pattern = "_") %>% unlist()
    input_tab$`Gene in network`[i] = intersect(genes,geneset) %>% paste0(collapse = ",")
    
  }
  input_tab$`Communication score` = round(input_tab$`Communication score`,digits = 3)
  # input_tab = input_tab[,-1]
  return(input_tab)
  
}


tab_cellchat = generate_Cellchat_data(input,geneset=the.Net$MisNodes$name)
tab_cellchat_end =  Integration_Cellchat(input_tab = tab_cellchat ,geneset = the.Net$MisNodes$name)

#  "MDK-NCL"  example
net = tab_cellchat$net[,, "MDK_NCL"]
rownames(net) = rownames(net) %>% text_Optimize(the.type = "output_celltype",input = input)
colnames(net) = colnames(net) %>% text_Optimize(the.type = "output_celltype",input = input)
demo_cellchat_plot = netVisual_circle(net,weight.scale = T,arrow.size = 0.2)






