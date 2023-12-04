####  Deconvolution of cancer prognostic signatures

library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(maxstat)
library(survival)
library(survminer)
library(patchwork)
library(igraph)


##################
#########  part1 |   top 500 prognostic genes(FPKM)  filter by gene expression 
##################

cancerTP_TCGA = readRDS("/data1/lizekun/CellNetdb/data/TCGA_project.rds")
cancerTP_TCGA$sub_project <- gsub(cancerTP_TCGA$project,pattern = "TCGA-",replacement = "")
cancerTP_TCGA$sub_project[grep("BRCA",cancerTP_TCGA$sub_project)] <- c("BRCA-ER","BRCA-HER2","BRCA-TNBC")

all_tcga_cancer <- cancerTP_TCGA$sub_project %>% unique() %>% sort() ## 24 

cancer_gene_list <- list()

for(i in 1:length(all_tcga_cancer)){
  
  used_cancer_type <- all_tcga_cancer[i]
  
  
  clinical_data <- data.table::fread(paste0("/data1/lizekun/TCGA/clinical_info/all_survival/",used_cancer_type,".txt"),
                                     sep = "\t",header = T)
  cancer_fpkm_data <- readRDS(paste0("/data1/lizekun/TCGA/FPKM/0.all_samples/",used_cancer_type,".rds"))
  
  cancer_fpkm_sample_df <- data.frame(raw_sample_id=colnames(cancer_fpkm_data))
  cancer_fpkm_sample_df$sub_id <- substr(cancer_fpkm_sample_df$raw_sample_id,1,12)
  cancer_fpkm_sample_df$label1 <- substr(cancer_fpkm_sample_df$raw_sample_id,14,15)
  table(cancer_fpkm_sample_df$label1)
  cancer_fpkm_sample_df <- cancer_fpkm_sample_df[!cancer_fpkm_sample_df$label1 %in% "11",]
  
  int_used_sample <- intersect(clinical_data$submitter_id,cancer_fpkm_sample_df$sub_id)
  clinical_data <- clinical_data[clinical_data$submitter_id %in% int_used_sample,]
  cancer_fpkm_sample_df <- cancer_fpkm_sample_df[match(int_used_sample,cancer_fpkm_sample_df$sub_id),]
  identical(clinical_data$submitter_id,cancer_fpkm_sample_df$sub_id)
  
  cancer_fpkm_data <- cancer_fpkm_data[,cancer_fpkm_sample_df$raw_sample_id]
  
  
  used_cancer_index <- c()
  for(j in 1:nrow(cancer_fpkm_data)){
    gene_used <- rownames(cancer_fpkm_data)[j]
    
    gene_exp <- cancer_fpkm_data[j,] %>%as.numeric()
    
    half_sample_size <- nrow(clinical_data)/2
    if(length(which(gene_exp > 1)) > half_sample_size){
      used_cancer_index <- c(used_cancer_index,j)
    }
    message(j)
  }
  
  
  used_cancer_gene_id <- rownames(cancer_fpkm_data)[used_cancer_index]
  
  raw_gene_path <- "/home/lizekun/projects/CellNetdb_analysis/Connectivity/0.prognostic_gene/1.fpkm/1.raw/"
  id_tran_df <- read.table("/data1/lizekun/txt/FPKM_gene_tran_df.txt",sep = "\t",header = T)
  
  all_gene_df <- read.table(paste0(raw_gene_path,used_cancer_type,".txt"),sep = "\t",header = T)
  all_gene_df$gene_symbol <- id_tran_df$gene_name[match(all_gene_df$gene,id_tran_df$gene_id)]
  
  filter_gene_df <- all_gene_df[all_gene_df$gene %in% used_cancer_gene_id ,]
  
  filter_gene_df = filter_gene_df[which(filter_gene_df$p_value_median< 0.05),]
  filter_gene_df$cancer=filter_gene_df[i] 
  filter_gene_df <- filter_gene_df[order(filter_gene_df$p_value_median,decreasing = F),]
  filter_gene_df <- filter_gene_df[1:500,]
  
  cancer_gene_list[[used_cancer_type]] <- filter_gene_df
  
  message(paste(used_cancer_type,"OVER"))
  
}



##################
#########  part2 |   top 500 prognostic genes connectivity 
##################

network_edge_list_STRING <- readRDS("/data1/lizekun/CellNetdb/add_analysis/00.pre/network_edge_list.rds")

solid_network_dir_df <- data.table::fread("/data1/lizekun/CellNetdb/add_analysis/txt/all_solid_cancer.txt",sep = "\t",header = T) %>% data.frame()
solid_network_dir_df$adj_cancer_type <- solid_network_dir_df$cancer_type %>% gsub(pattern = " ",replacement = '_')
solid_network_dir_df$net_name <- solid_network_dir_df$rds %>% gsub(pattern = ".RDS",replacement = '')



cal.connectivity <- function(network.edge.list = NULL, geneset = NULL){
  
  node.list <- lapply(network.edge.list, function(net){
    genes.all <- unique(c(as.character(net[,1]),as.character(net[,2])))
  })
  names(node.list) <- names(network.edge.list)
  
  connectivity.sig <- lapply(network.edge.list,function(net){
    nrow(net[(net[,1] %in% geneset & net[,2] %in% geneset), ])
  })
  
  connectivity.sig.all <- as.data.frame(dplyr::bind_rows(connectivity.sig))
  df <- t(connectivity.sig.all)
  data <- as.data.frame(df)
  colnames(data) <- 'connectivity'
  
  data$signature_gene_num <- rep(length(geneset), nrow(data))
  
  sig.f <- lapply(node.list, function(node){length(geneset[geneset %in% node])})
  
  data$detected.sig.num <- unlist(sig.f)
  
  return(data)
  
}


int_output_df <- data.frame()

for(i in 1:nrow(cancerTP_TCGA)){
  
  this.cancer <- cancerTP_TCGA$cancer_type[i]
  this.tcga.project <- cancerTP_TCGA$sub_project[i]
  
  prognostic_gene_df <- cancer_gene_list[[this.tcga.project]]
  this.cancer_df <- solid_network_dir_df[solid_network_dir_df$adj_cancer_type %in% this.cancer,]
  output_df <- cal.connectivity(network_edge_list_STRING[this.cancer_df$net_name],prognostic_gene_df$gene_symbol)
  output_df$connectivity.normalized <- output_df$connectivity / output_df$detected.sig.num  
  output_df$cancer <- this.cancer_df$cancer_type
  output_df$cell_type <- this.cancer_df$cell_type   
  
  int_output_df <- rbind2(int_output_df ,output_df)
  
  message(i)
  
}



write.table(int_output_df,"/home/lizekun/projects/CellNetdb_analysis/Connectivity/deconvolute_prognostic_gene/deconvolute_top500_out.txt",
            sep = "\t",col.names = T,row.names = T,quote = F)



##################
#########  part3 |   evaluate significance
##################

sig.Connectivity <- function(network = NULL, geneset = NULL, simulate.num = 10000){
  
  network.graph <- igraph::graph_from_data_frame(network)
  degree.centrality <- igraph::degree(network.graph)
  
  detected.genes <- geneset[geneset %in% names(degree.centrality)]

  pb = txtProgressBar(min = 0, max = simulate.num, style = 3)
  
  connectivity.random <- vector()
  for (n in 1:simulate.num){
    Sys.sleep(0.05)
    setTxtProgressBar(pb,n)
    
    random.geneset <- vector()
    selected <- vector()
    for (i in seq_along(detected.genes)){
      gene <- detected.genes[i]
      gene.degree <- degree.centrality[gene]
      degree.centrality.f <- degree.centrality[!(names(degree.centrality) %in% selected)]
      degree.range <- c(floor(0.8 * gene.degree), ceiling(1.2 * gene.degree))
      node.pool <- degree.centrality.f[degree.range[1] <= degree.centrality.f & degree.centrality.f <= degree.range[2]]
      
      if (length(node.pool) == 0){
        break
      }
      
      gene.random <- sample(names(node.pool), 1)
      
      random.geneset[i] <- gene.random
      selected[i] <- gene.random
    }
    
    connectivity.random.n <- nrow(network[(network[,1] %in% random.geneset & network[,2] %in% random.geneset), ])
    connectivity.random[n] <- connectivity.random.n
  }
  close(pb)
  
  
  connectivity <- nrow(network[(network[,1] %in% geneset & network[,2] %in% geneset), ])
  connectivity.final <- c(connectivity.random, connectivity)
  pvalue <- rank(-(connectivity.final),ties.method = 'last')[simulate.num + 1] / simulate.num
  
  output.list <- list(null.distribution = connectivity.random, p.value = pvalue, detected.geneset = detected.genes, observed = connectivity)
  
  return(output.list)
}

for(i in 1:nrow(cancerTP_TCGA)){
  
  this_cancer_type <- cancerTP_TCGA$cancer_type[i]
  this_tcga_project <- cancerTP_TCGA$sub_project[which(cancerTP_TCGA$cancer_type %in% this_cancer_type)]
  
  prognostic_gene_df <- cancer_gene_list[[this_tcga_project]]
  
  this.cancer_df <- solid_network_dir_df[solid_network_dir_df$adj_cancer_type %in% this_cancer_type,]
  
  cancer_Connectivity_list <- list()
  
  
  for(j in 1:nrow(this.cancer_df)){
    this.net.name <- this.cancer_df$net_name[j]
    
    this.net_out <- sig.Connectivity(network = network_edge_list_STRING[[this.net.name]],geneset = prognostic_gene_df$gene_symbol)
    
    cancer_Connectivity_list[[this.net.name]] <- this.net_out
    message(paste0(this.net.name,"      OVER"))
  }
  
  saveRDS(cancer_Connectivity_list,paste0("/home/lizekun/projects/CellNetdb_analysis/Connectivity/deconvolute_prognostic_gene/connectivity_pvalue/top500/",
                                          this_cancer_type,".rds"))
  
  
}

