####  network Shared-edge  Similarity   

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(corrplot)
library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)
library(lubridate)
library(scales)
library(stringr)

my.cor.test_spearman <- function(...) {
  obj<- try(cor.test(...,method="spearman"), silent=TRUE)
  if (is(obj, "try-error"))  return(c(NA,NA)) else return(obj[c("estimate","p.value")])
}


used_solid_network_dir_df <- data.table::fread("/data1/lizekun/CellNetdb/add_analysis/txt/all_solid_cancer_plus.txt",sep = "\t",header = T) %>% data.frame()
used_solid_network_dir_df$cell_type %>% unique() %>% sort()

select_cell_type <- c("B cells","CD4+ T cells","CD8+ T cells","Endothelial cells","Epithelial cells (malignant)","Fibroblasts","Myeloid cells","NK cells")

select_used_solid_network_dir_df <- used_solid_network_dir_df[used_solid_network_dir_df$cell_type %in% select_cell_type,]
tab <- table(select_used_solid_network_dir_df$adj_cancer_type) %>% data.frame()
select_cancer <- tab$Var1[tab$Freq==8] %>% as.character()              ###  17 cancer types 
select_used_solid_network_dir_df <- select_used_solid_network_dir_df[select_used_solid_network_dir_df$adj_cancer_type %in% select_cancer,]  


pair_net_df <- combn(select_used_solid_network_dir_df$net_name,2) %>% t() %>% data.frame()
colnames(pair_net_df) <- c("net1","net2")
pair_net_df$shared_node_num <- NA
pair_net_df$shared_node_edge_num <- NA

pair_net_df$Spearman <- NA 



######  net_edge_df List

net_df_list = list()

for(i in 1:nrow(select_used_solid_network_dir_df)){

  the.net_name = select_used_solid_network_dir_df$net_name[i]

  net_ig = readRDS( paste0("/data1/lizekun/CellNetdb/data/new_Network/",the.net_name,".RDS"))

  the.net = net_ig$net
  the.node =net_ig$node


  the.net_df = reshape2::melt(the.net)
  the.net_df = the.net_df[which(the.net_df$value != 0),]

  the.net_df$Var1 = the.net_df$Var1 %>% as.character()
  the.net_df$Var2 = the.net_df$Var2 %>% as.character()

  the.net_df$gene1ID = the.net_df$Var1 > the.net_df$Var2
  the.net_df$gene1ID = the.net_df$gene1ID + 1

  the.net_df$gene2ID = the.net_df$Var1 < the.net_df$Var2
  the.net_df$gene2ID = the.net_df$gene2ID + 1

  the.net_df$edge = ""
  for(j in 1:nrow(the.net_df)){

    gene1 = the.net_df[ j , the.net_df$gene1ID[j] ]
    gene2 = the.net_df[ j , the.net_df$gene2ID[j] ]
    the.net_df$edge[j] = paste(gene1,gene2,sep = "_")

  }

  net_df_list[[the.net_name]] = the.net_df
  message(i)


}

lapply(net_df_list,nrow) %>% unlist()



for(i in 1:nrow(pair_net_df)){
  
  net1_name <- pair_net_df$net1[i]
  net2_name <- pair_net_df$net2[i]
  
  net1 <- net_df_list[[net1_name]]
  net2 <- net_df_list[[net2_name]]
  
  net1_gene <- c(net1$Var1,net1$Var2) %>% unique()
  net2_gene <- c(net2$Var1,net2$Var2) %>% unique()
  
  int_gene <- intersect(net1_gene,net2_gene)
  
  sub_net1 <- net1[net1$Var1 %in% int_gene & net1$Var2 %in% int_gene,]
  sub_net2 <- net2[net2$Var1 %in% int_gene & net2$Var2 %in% int_gene,]
  
  shared_edge <- intersect(sub_net1$edge,sub_net2$edge)
  
  pair_net_df$shared_node_num[i] <- length(int_gene)
  pair_net_df$shared_node_edge_num[i] <- length(shared_edge)
  
  if(length(shared_edge) ==0){
    i=i+1
  }else{
    
    adj_sub_net1 <- sub_net1[match(shared_edge,sub_net1$edge),]
    adj_sub_net2 <- sub_net2[match(shared_edge,sub_net2$edge),]
    
    corr_out <- my.cor.test_spearman(adj_sub_net1$value %>%as.numeric(),adj_sub_net2$value %>%as.numeric())

    
    pair_net_df$Spearman[i] <- as.numeric(corr_out[1])
    
  }
  
  message(i)
}






pair_net_df$Spearman[which(pair_net_df$shared_node_edge_num < 10)] <- NA
pair_net_df$Spearman[which(pair_net_df$Spearman < 0)] <- NA


####################
####  turn to matrix 
####################


corr_mat <- matrix(ncol = 136,nrow = 136)
rownames(corr_mat) <- c(pair_net_df$net1,pair_net_df$net2) %>% unique() %>% sort()
colnames(corr_mat) <- rownames(corr_mat)


Spearman_corr_mat <- corr_mat


for(j in 1:nrow(pair_net_df)){
  net1_name <- pair_net_df$net1[j]
  net2_name <- pair_net_df$net2[j]
  
  Spearman_corr_mat[net1_name,net2_name] <- pair_net_df$Spearman[j]
  Spearman_corr_mat[net2_name,net1_name] <- pair_net_df$Spearman[j]
  
}

diag(Spearman_corr_mat) <- 1

adj_Spearman_corr_mat <- Spearman_corr_mat
rownames(adj_Spearman_corr_mat) <- gsub(rownames(adj_Spearman_corr_mat),pattern = "_STRING_",replacement = " ")
rownames(adj_Spearman_corr_mat) <- gsub(rownames(adj_Spearman_corr_mat),pattern = "_",replacement = " ")
rownames(adj_Spearman_corr_mat) <- gsub(rownames(adj_Spearman_corr_mat),pattern = "[(]malignant[)]",replacement = "")
colnames(adj_Spearman_corr_mat) <- rownames(adj_Spearman_corr_mat)




####  ComplexHeatmap plot 

col1=colorRampPalette(colors =c("red","white","darkgreen"),space="Lab") 

adj_plot_Spearman_corr_mat <- adj_Spearman_corr_mat[adj_net_name,adj_net_name]

plot_path <- "/data1/lizekun/CellNetdb/add_analysis/net_similarity/Shared_edge/plot/"


pdf(paste0(plot_path,"shared_edge_Spearman_corr.pdf"),width = 22,height = 18.5)
Hplot <- ComplexHeatmap::Heatmap(adj_plot_Spearman_corr_mat,
                                 show_row_names = T,cluster_rows = F,
                                 show_column_names = F,cluster_columns =F,
                                 row_names_gp = gpar(fontsize= 11) ,
                                 show_row_dend = F,
                                 show_column_dend = F,
                                 na_col = "#EAEAEA",
                                 # na_col = "#E2E1DC",
                                 column_title_side = "top",
                                 row_names_side = "left",
                                 # heatmap_legend_param = list(
                                 #   title = "Shared-edge Similarity", title_gp = gpar(fontsize = 16), at = c(-0.5,0,0.5,1), title_position= "topcenter", 
                                 #   title_gp = gpar(col = "black"),
                                 #   labels = c("-0.5","0","0.5","1"),legend_direction = "horizontal",legend_width = unit(8, "cm"),labels_gp = gpar(fontsize = 18)
                                 # ),
                                 # col=c(col1(25)[1:12],col1(50)[26:50]),
                                 heatmap_legend_param = list(
                                   title = "Shared-edge Similarity", title_gp = gpar(fontsize = 16), at = c(0,1), title_position= "topcenter", 
                                   title_gp = gpar(col = "black"),
                                   labels = c("0","1"),legend_direction = "horizontal",legend_width = unit(8, "cm"),labels_gp = gpar(fontsize = 18)
                                 ),
                                 col= col1(50)[26:50],
                                 row_names_max_width = unit(12, 'cm'))


draw(Hplot,heatmap_legend_side = "bottom")
dev.off()







