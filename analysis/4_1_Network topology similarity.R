####  network Topology  Similarity   

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


net_data_path <- "/data1/lizekun/CellNetdb/data/new_Network/"

net_df_list = list()
for(i in 1:nrow(select_used_solid_network_dir_df)){
  
  the.net_name = select_used_solid_network_dir_df$rds[i]
  
  net_ig = readRDS( paste(net_data_path,the.net_name,sep = "")  )
  
  the.node =net_ig$node
  
  adj_the.net_name <- the.net_name %>% gsub(pattern = ".RDS",replacement = "")
  net_df_list[[adj_the.net_name]] = the.node
  message(i)
}



get_node_specificity_sim <- function(net1_name,net2_name,cor_method){
  
  net1 = net_df_list[[net1_name]]
  net2 = net_df_list[[net2_name]]
  
  shared_node = intersect(net1$gene,net2$gene)
  
  net1_f = net1[match(shared_node,net1$gene),]
  
  net2_f = net2[match(shared_node,net2$gene),]
  
  return( cor(net1_f$specificity,net2_f$specificity,method = cor_method) )
  
}


corr_mat = matrix(0,nrow = nrow(select_used_solid_network_dir_df),
                  ncol = nrow(select_used_solid_network_dir_df))
rownames(corr_mat) = select_used_solid_network_dir_df$net_name
colnames(corr_mat) = rownames(corr_mat)


for(i in 1:nrow(corr_mat)){
  for(j in 1:nrow(corr_mat)){
    
    net1 = rownames(corr_mat)[i]
    net2 = rownames(corr_mat)[j]
    corr_mat[i,j] = get_node_specificity_sim(net1,net2,cor_method="spearman")

  }
  message(i)
}



adj_corr_mat <- corr_mat
rownames(adj_corr_mat) <- gsub(rownames(adj_corr_mat),pattern = "_STRING_",replacement = " ")
rownames(adj_corr_mat) <- gsub(rownames(adj_corr_mat),pattern = "_",replacement = " ")
rownames(adj_corr_mat) <- gsub(rownames(adj_corr_mat),pattern = "[(]malignant[)]",replacement = "")
colnames(adj_corr_mat) <- rownames(adj_corr_mat)

data_df <- data.frame(raw=rownames(corr_mat),new=colnames(adj_corr_mat))

plot_corr_mat <- adj_corr_mat[adj_net_name,adj_net_name]    
plot_corr_mat[plot_corr_mat<0 ] <- NA


plot_path <- "/data1/lizekun/CellNetdb/add_analysis/net_similarity/node_specificity/plot/"
pdf(paste0(plot_path,"spearman_corr.pdf"),width = 22,height = 18.5)
Hplot <- ComplexHeatmap::Heatmap(plot_corr_mat,
                                 show_row_names = T,cluster_rows = F,
                                 show_column_names = F,cluster_columns =F,
                                 row_names_gp = gpar(fontsize= 11) ,
                                 show_row_dend = F,
                                 show_column_dend = F,
                                 na_col = "#EAEAEA",
                                 column_title_side = "top",
                                 row_names_side = "left",
                                 heatmap_legend_param = list(
                                   title = "Topology Similarity", title_gp = gpar(fontsize = 16), at = c(0,1), title_position= "topcenter", 
                                   title_gp = gpar(col = "black"),
                                   labels = c("0","1"),legend_direction = "horizontal",legend_width = unit(8, "cm"),labels_gp = gpar(fontsize = 18)
                                 ),
                                 col=brewer.pal(11,"PuOr")[6:10],
                                 row_names_max_width = unit(12, 'cm'))


draw(Hplot,heatmap_legend_side = "bottom")
dev.off()



