###  Infomap  |Community Detection 

library(tidygraph)
library(tidyverse)
library(ggraph)
library(igraph)
library(ggsci)
library(viridis)


net_df_list <- readRDS("/data1/lizekun/CellNetdb/add_analysis/05.CD4T_similarity/0.pre/32_solid_CD4Tnet_edge_df_list.rds")

ct_namechange = read.table("/data1/lizekun/CellNetdb/data/ct_namechange.txt",sep = "\t",header = T)

net_file = dir("/data1/lizekun/CellNetdb/data/new_Network/")
net_file = net_file[which(str_detect(net_file,pattern = "STRING") & str_detect(net_file,pattern = "CD4[+] T"))]

net_df = data.frame(
  cancer_type = "",
  cell_type = "",
  file = net_file
)
for(i in 1:nrow(net_df)){
  
  the_path = net_df$file[i]
  the_path = str_split(the_path,pattern = "_STRING_") %>% unlist()
  the.cancer = ct_namechange$new_cancer[which(ct_namechange$old_cancer == the_path[1])] %>% unique()
  the.celltype = ct_namechange$new_celltype[which(ct_namechange$old_celltype == gsub(the_path[2],pattern = "[.]RDS",replacement = ""))] %>% unique()
  
  net_df$cancer_type[i] = the.cancer
  net_df$cell_type[i] = the.celltype
  
}


netList = list()
for(i in 1:nrow(net_df)){

  the.path = net_df$file[i]
  net = readRDS( paste("/data1/lizekun/CellNetdb/data/new_Network/",the.path,sep = "")  )
  netList[[net_df$cancer_type[i]]] = net$net
  message(i)

}

lapply(netList, nrow) %>% unlist()

# saveRDS(netList,"/data1/lizekun/CellNetdb/add_analysis/08.infomap/CD4T/0.pre/net_df_list.rds")


################
# con - network 
################


replenish_Net <- function(net, Sum_geneset){
  
  # for(i in 1:nrow(net)){
  #   for(j in i:ncol(net)){
  #     net[i,j] = 0
  #   }
  # }
  
  added_gene = setdiff(Sum_geneset,rownames(net))
  
  added_mat_r = matrix( 0,nrow = nrow(net),ncol = length(added_gene) )
  rownames(added_mat_r) = rownames(net)
  colnames(added_mat_r) = added_gene
  net = cbind(net,added_mat_r)
  
  added_mat_d = matrix( 0,nrow = length(added_gene),ncol = length(Sum_geneset))
  rownames(added_mat_d) = added_gene  
  colnames(added_mat_d) = c(colnames(net))
  net = rbind(net,added_mat_d)
  
  net = net[Sum_geneset,Sum_geneset]
  
  return(net)
  
}


fun_for_GO <- function(total_geneset,total_universe){
  
  GO_gene_set = readRDS("/data1/lizekun/CellNetdb/data/GO_gene_set.RDS")
  GO.list = readRDS("/data1/lizekun/CellNetdb/data/GO_reference.RDS")

  
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
  
  overlap = overlap[which(overlap$type == "GO:BP"),]
  overlap$adj_p = overlap$pvalue %>% p.adjust(method = "fdr")
  overlap = overlap[which(overlap$pvalue<0.01),]
  overlap = overlap[order(overlap$adj_p,decreasing = F),]
  
  return(overlap)
  
}

fun_for_KEGG <- function(total_geneset,total_universe){
  
  kegg_list = readRDS("/data1/lizekun/CellNetdb/data/KEGG_list.RDS")
  Disease_gene_set = readRDS("/data1/lizekun/CellNetdb/data/KEGG_geneset.RDS")
  
  
  universe = intersect( total_universe , Disease_gene_set )
  geneset = intersect( total_geneset , Disease_gene_set )
  
  overlap = data.frame(
    
    disease = names(kegg_list),
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
    
    overlap_for_inter = intersect(kegg_list[[disease_name]] , geneset)
    if(length(overlap_for_inter)<2 ){
      next
    }
    
    overlap_for_uniANDref = intersect(kegg_list[[disease_name]],universe )
    overlap_for_uniDIFFref = setdiff(universe , kegg_list[[disease_name]])
    
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
  
  colnames(overlap)[1] = "KEGG"
  overlap = overlap[which(overlap$pvalue!=1),]
  overlap$adj_p = overlap$pvalue %>% p.adjust(method = "fdr")
  overlap = overlap[order(overlap$adj_p,decreasing = F),]
  
  return(overlap)
  
}


con_network <- function(input_netlist,weight = 3){
  
  sum_nodes = c()
  for(i in 1:length(input_netlist)){
    
    the_net = input_netlist[[i]]
    sum_nodes = c(sum_nodes,rownames(the_net))
    
  }
  sum_nodes = unique(sum_nodes)
  
  for(i in 1:length(input_netlist)){
    input_netlist[[i]] = replenish_Net(input_netlist[[i]],sum_nodes)
    print(i)
  }
  
  con_net = matrix(1,nrow = length(sum_nodes),ncol = length(sum_nodes))
  rownames(con_net) = sum_nodes
  colnames(con_net) = sum_nodes
  
  for(i in 1:length(input_netlist)){
    
    the_net = input_netlist[[i]]
    the_net[which(the_net < weight)] = 0
    
    con_net[which(the_net == 0)] = 0
    
  }
  
  return(con_net)
  
}


con_netSet = c(
  "Bladder urothelial carcinoma",
  "Prostate adenocarcinoma"
)


con_netList <- netList[con_netSet]

con_net = con_network(con_netList,weight = 2.5)
aaa = rowSums(con_net)
bbb = colSums(con_net)

con_net = con_net[which(aaa != 0 & bbb !=0) , which(aaa != 0 & bbb !=0)]

con.graph <- as_tbl_graph(con_net, directed = FALSE)
con.graph = con.graph %>%
  activate(nodes) %>%
  mutate(cluster = as.factor(group_infomap()))

nodes = data.frame(
  genes = V(con.graph)$name,
  cluster = V(con.graph)$cluster
)


write.table(nodes,"/data1/lizekun/CellNetdb/add_analysis/08.infomap/CD4T/1.cluster_out/BLCA_PRAD_cluster.txt",
            sep = "\t",col.names = T,row.names = F,quote = F)


####  show top15 centrality genes in top  cluster 


out_list = list()
for(cluster in 1:9){
  
  used_gene = nodes$genes[which(nodes$cluster == cluster)]
  sub_net = con_net[used_gene,used_gene]
  sub.graph <- as_tbl_graph(sub_net)
  
  sub.graph <- sub.graph %>%
    activate(nodes) %>%
    mutate(centrality = centrality_authority())

  V(sub.graph)$label = ""
  out_gene <- rank(V(sub.graph)$centrality) %>% sort(decreasing = T) %>% head(15) %>%names()
  V(sub.graph)$label[match(out_gene,V(sub.graph)$name)] = out_gene
  
  
  the_nodes = data.frame(
    gene = V(sub.graph)$name,
    centrality = V(sub.graph)$centrality
  )
  the_nodes = the_nodes[order(the_nodes$centrality,decreasing = T),]
  
  p = sub.graph %>%
    ggraph(layout = "graphopt") +
    geom_edge_link(width = 0.5, colour = "lightgray") +
    geom_node_point(size = 4,aes(color = centrality) )+
    theme_graph()+
    geom_node_text(aes(label = label), repel = TRUE)+
    coord_fixed(ratio = 1)
  
  
  
  out_list[[cluster]] = list(
    go = fun_for_GO(total_geneset = nodes$genes[which(nodes$cluster == cluster)],total_universe = nodes$genes),
    kegg = fun_for_KEGG(total_geneset = nodes$genes[which(nodes$cluster == cluster)],total_universe = nodes$genes),
    plot = p,
    node = the_nodes,
    graph=sub.graph
    
  )
  
  # print(p)
  
}

# saveRDS(out_list,"/data1/lizekun/CellNetdb/add_analysis/08.infomap/CD4T/1.fun_list_out/BLCA_PRAD.rds")



