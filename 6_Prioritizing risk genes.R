### Prioritize risk genes based on random walk with restart(RWR)

library(stringi)
library(stringr)
library(dplyr)
library(tidyverse)
library(igraph)



#' Title. Convert the upper triangular matrix to a symmetric matrix
#' 
#' @param net a upper triangular matrix
#'
#' @return Symmetric matrices

get_W <- function(net) {
  
  net = net + t(net) - diag(diag(net))
  return(net)
  
}

#' Title. Calculate the transition probability matrix, which is the 1/degrees
#' 
#' @param net a upper triangular matrix
#' @return transition probability matrix

get_Q <- function(net){
  
  if(!isSymmetric(net[1:100,1:100])){
    net = get_W(net)
  }
  
  WIj = rowSums(net)
  
  for(i in 1:nrow(net)){
    net[i,] = net[i,] / WIj
  }
  
  return(net)
}

#' Title. get the vector of initial probability distribution for all genes
#' 
#' @param net a upper triangular matrix
#' @param seeds seed genes
#' 
#' @return the vector of initial probability distribution

get_R0 <- function(net,seeds){
  
  is_seed = rep(0,nrow(net))
  is_seed[which(rownames(net) %in% seeds)] = 1
  is_seed = is_seed/sum(is_seed)
  
  return(is_seed)
  
}

#' Title get the final gene scores
#'
#' @param net  a parameter to measure the importance of genes and interactions
#' @param seeds the vector of initial disease risk scores for all genes
#' @param threshold the threshold for ending the iteration
#' 
#' @return final gene probability distribution (vector)
#' @export

get_R <- function(net, seeds , threshold = 10^(-10) , bet = 0.5) {
  
  QQ = get_Q(net)
  R_0 = get_R0(net,seeds)
  
  R_dise <- 1
  R_old <- R_0 
  kk <- 1
  
  while (R_dise > threshold) {
    
    R_new <- (1 - bet) * QQ %*% R_old + bet * R_0
    Po_old = norm(as.matrix(R_old), "F")
    Po_new = norm(as.matrix(R_new), "F")
    R_dise <- abs(Po_new - Po_old)
    R_old <- R_new
    kk <- kk + 1
    
    print(kk)
    
  }
  
  print(kk)
  names(R_new) = rownames(QQ)
  R_new = sort(R_new,decreasing = T)
  
  return(R_new)
  
}


###  demo : CD8T Tex network 
CGC = read_tsv("./CancerDriver/CGC/CGC.tsv")                                                      ####  https://cancer.sanger.ac.uk/census
NCG = read_tsv("./CancerDriver/NCG/NCG_cancerdrivers_annotation_supporting_evidence.tsv")         ####  http://ncg.kcl.ac.uk/download.php

seeds_in_CGC = c("PDCD1",
                 "CD274",
                 "BTLA",
                 "HAVCR2",
                 "LAG3",
                 "TIGIT")
genes_in_NCG = NCG$symbol %>% unique()

cp_net = readRDS("./data/new_Network(Imm)/CD8_Tcells_STRING_Tex.RDS") 
cp_net <- cp_net$net
cp_net[cp_net>0] <- 1
cp_Gene_prio = get_R(cp_net ,seeds = seeds_in_CGC)

cp_Gene_prio_df <- data.frame(gene=names(cp_Gene_prio),score=as.numeric(cp_Gene_prio))



###  plot top30 genes

top30_gene_plot_df <- cp_Gene_prio_df[!cp_Gene_prio_df$gene %in% seeds_in_CGC,]
top30_gene_plot_df <- top30_gene_plot_df[1:30,]
top30_gene_plot_df$color = "black"
top30_gene_plot_df$color[which(top30_gene_plot_df$gene %in% genes_in_NCG)] = "red"

ggplot(top30_gene_plot_df,aes(x = gene,y = score,fill = score ))+
  geom_bar(stat = "identity")+
  scale_x_discrete(limit = top30_gene_plot_df$gene %>% rev())+coord_flip()+
  theme_classic()+
  scale_fill_viridis_c()+
  theme(axis.text.y = element_text(color = top30_gene_plot_df$color %>% rev()))








