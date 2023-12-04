###  Taxonomy plot 
library(tidyverse)
library(Seurat)
library(patchwork)
library(dplyr)
library(cowplot)
library(stringi)
library(stringr)
library(plotly)
library(readr)
library(hexbin)
library(ggplot2)


############
## UMAP plot
############

data <- readRDS("/data1/lizekun/sc_cancer/2.re_cluster/total_cells/RDS/Solid_tumor/Adrenal_neuroblastoma.rds")   ### Seurat object
plot.data = as.data.frame(data@reductions$umap@cell.embeddings)
plot.data$celltype = as.character(data@active.ident[match(rownames(plot.data),rownames(data@meta.data))])


theme1 <-   theme(panel.background = element_rect(fill='white', colour='black'),
                  panel.grid=element_blank(), axis.title = element_text(color='black',size=35),
                  legend.title=element_blank(),legend.text=element_text(size=35),
                  legend.key=element_blank(),plot.margin = unit(c(10,10,10,5),'lines'),
                  legend.position = 'none')
theme2 <- theme(panel.background = element_rect(fill='white', colour='black'),
                panel.grid=element_blank(), axis.title = element_text(color='black',size=35),
                legend.title=element_blank(),legend.text=element_text(size=35),legend.position = c(0.22,0.5),
                legend.key=element_blank(),legend.key.width = unit(50,"pt"),legend.key.height = unit(50,"pt"))


plot <- ggplot(plot.data,aes(x = UMAP_1,y= UMAP_2))+
  geom_point(aes(color=celltype),size = 0.5)+
  scale_color_manual(values=c("#FFD100","#DEB887","#FFB5C5","#82589F","#827717","#00e676","#7FBC41","#bf360c","#35978F"))+
  theme(text = element_text(size=35))+
  ggtitle("Adrenal neuroblastoma(160,566 cells)")+
  guides(colour = guide_legend(override.aes = list(size=8)))+
  theme(plot.title = element_text(size=35,colour = "black",hjust = 0.5))+
  coord_fixed(ratio = (max(plot.data$UMAP_1)-min(plot.data$UMAP_1))/(max(plot.data$UMAP_2)-min(plot.data$UMAP_2)))

plot1 <- plot+theme1
plot2 <- plot+theme2

legend = cowplot::get_legend(plot2)
legend.plot <- plot_grid(legend)
plot1+legend.plot

pdf(file = "/data1/lizekun/sc_cancer/2.re_cluster/total_cells/UMAP_pdf/Adrenal_neuroblastoma.pdf",width = 25,height = 15)
print(plot1+legend.plot)
dev.off()


########################
## Hexbin marker plot 
########################

.make_hexbin_helper <- function(dr, nbins = 80, use_dims) {
  if (dim(dr)[2] < max(use_dims)) {
    stop("Please specify use_dims that are calculated.")
  }
  
  xbnds <- range(c(dr[, use_dims[1]]))
  ybnds <- range(c(dr[, use_dims[2]]))
  
  drhex <- hexbin(dr[, use_dims[1]],
                  dr[, use_dims[2]],
                  nbins,
                  xbnds = xbnds,
                  ybnds = ybnds,
                  IDs = TRUE
  )
  cID <- drhex@cID
  drhex <- cbind(
    as.numeric(hcell2xy(drhex)$x),
    as.numeric(hcell2xy(drhex)$y),
    as.numeric(drhex@count)
  )
  
  colnames(drhex) <- c("x", "y", "number_of_cells")
  
  res <- list(cID = cID, hexbin.matrix = drhex)
  
  return(res)
}

.make_hexbin_function <- function(x, action, cID, na.rm=FALSE, no) {
  if (action == "majority") {
    func_if <- !(is.factor(x) | is.character(x))
    
    if (func_if) {
      stop("For action 'majority' x needs to be a factor or character.")
    } else {
      res <- tapply(x, cID, FUN = function(z) {
        names(sort(table(z),
                   decreasing = TRUE
        )[1])
      })
      res <- as.factor(res)
      return(res)
    }
  }
  
  if (action == "prop") {
    func_if <- !(is.factor(x) | is.character(x))
    
    if (func_if) {
      stop("For action 'prop' x needs to be a factor or character.")
    } else {
      res <- tapply(x, cID, FUN = function(z) {
        sum(z == unique(x)[no], na.rm = na.rm) / sum(!is.na(z))
      })
      res <- as.numeric(res)
      return(res)
    }
  }
  
  if (action == "median") {
    func_if <- !is.numeric(x)
    
    if (func_if) {
      stop("For action 'median' x needs to be numeric")
    } else {
      res <- tapply(x, cID, FUN = function(z) median(z, na.rm = na.rm))
      res <- as.numeric(res)
      return(res)
    }
  }
  
  if (action == "mode") {
    func_if <- !is.numeric(x)
    
    if (func_if) {
      stop("For action 'mode' x needs to be numeric")
    } else {
      res <- tapply(x, cID, FUN = function(z) .get_mode(z))
      res <- as.numeric(res)
      return(res)
    }
  }
  
  if (action == "prop_0") {
    func_if <- !is.numeric(x)
    
    if (func_if) {
      stop("For action 'prop_0' x needs to be numeric")
    } else {
      res <- tapply(x, cID, FUN = function(z) sum(z > 0, na.rm = na.rm) / 
                      sum(!is.na(z)))
      res <- as.numeric(res)
      return(res)
    }
  }
  
  if (action == "sd") {
    func_if <- !is.numeric(x)
    
    if (func_if) {
      stop("For action 'prop_0' x needs to be numeric")
    } else {
      res <- tapply(x, cID, FUN = function(z) sd(z, na.rm = na.rm))
      res <- as.numeric(res)
      return(res)
    }
  }
  
  if (action == "mean") {
    func_if <- !is.numeric(x)
    
    if (func_if) {
      stop("For action 'mean' x needs to be numeric")
    } else {
      res <- tapply(x, cID, FUN = function(z) mean(z, na.rm = na.rm))
      res <- as.numeric(res)
      return(res)
    }
  } else {
    stop("Specify valid action!")
  }
}

.plot_hexbin <- function(drhex,
                         colour_by = "Cluster_majority",
                         action,
                         colors = NULL,
                         title = NULL,
                         xlab = NULL,
                         ylab = NULL) {
  if (any(!c("x", "y", colour_by) %in% colnames(drhex))) {
    stop("The dataframe must contain columns named 'x', 'y' and label.")
  }
  
  if (is.null(xlab)) {
    xlab <- "x"
  }
  
  if (is.null(ylab)) {
    ylab <- "y"
  }
  
  if (action=="majority") {
    if (is.null(colors)) {
      ggplot(drhex, aes_string("x", "y", fill = colour_by)) +
        geom_hex(stat = "identity") +
        theme_classic() + theme(legend.position = "bottom") +
        ggtitle(title) +
        labs(x = xlab, y = ylab) + theme(legend.title = element_blank())
    } else {
      ggplot(drhex, aes_string("x", "y", fill = colour_by)) +
        geom_hex(stat = "identity") + scale_fill_manual(values = colors) +
        theme_classic() + theme(legend.position = "bottom") +
        ggtitle(title) +
        labs(x = xlab, y = ylab) + theme(legend.title = element_blank())
    }
  } else {
    ggplot(drhex, aes_string("x", "y", fill = colour_by)) +
      geom_hex(stat = "identity") +
      theme_classic() + scale_fill_viridis_c() + ggtitle(title) +
      labs(x = xlab, y = ylab) + theme(legend.title = element_blank())
  }
}



load("./data/Taxonomy_data.RData")


###CT
name_change = data.frame( final_cancer = unique(Taxonomy_cancer_data$`Cancer type`), raw_cancer = "" )
name_change$raw_cancer = gsub(name_change$final_cancer,pattern = " ",replacement = "_")
name_change$type = Taxonomy_cancer_data$Projects[match(name_change$final_cancer,Taxonomy_cancer_data$`Cancer type`)]
name_change$raw_cancer[which(name_change$raw_cancer == "T_cell_acute_lymphoblastic_leukemia")] = "T_cALL"
name_change$raw_cancer[which(name_change$raw_cancer == "Mixed_phenotype_acute_leukemia")] = "Mixed_phenotype_acute_leukemias"
name_change$markers = ""

for(i in 1:nrow(name_change)){
  cancer_final = name_change$final_cancer[i]
  tab = Taxonomy_cancer_data[which(Taxonomy_cancer_data$`Cancer type` == cancer_final),]
  markers = paste0(tab$`Marker gene`,collapse = ",") %>% str_split(pattern = ",") %>% unlist() %>% unique()
  markers = paste0(markers,collapse = ",")
  name_change$markers[i] = markers
}

for(i in 1:nrow(name_change)){
  
  type = name_change$type[i]
  if(type == "Solid tumor"){type = "Solid_tumor"}else{type = "Hematological_tumor/"}
  cancer_raw = name_change$raw_cancer[i]
  cancer_final = name_change$final_cancer[i]
  markergenes = name_change$markers[i] %>% str_split(pattern = ",") %>% unlist()
  
  message(cancer_final)
  
  seurat_Data = readRDS(paste0("/data1/lizekun/sc_cancer/2.re_cluster/total_cells/RDS/",type,"/",cancer_raw,".rds"))
  
  if(unique(seurat_Data$orig.ident) %>% length() == 1){
    
    counts = seurat_Data@assays$RNA@data
    markergenes = intersect(markergenes,rownames(counts))
    
    tab = counts[markergenes,]%>% as.matrix() %>% t()  %>% as.data.frame()
    
  }else{
    
    counts = seurat_Data@assays$CCA@data
    the.count = seurat_Data@assays$RNA@data
    
    tab = counts[intersect(markergenes,rownames(counts)),] %>% as.matrix() %>% t()  %>% as.data.frame()
    
    if( length( setdiff(markergenes,rownames(counts)) ) == 1){
      tab_ = the.count[setdiff(markergenes,rownames(counts)),]%>% as.matrix()  %>% as.data.frame()
      colnames(tab_) = setdiff(markergenes,rownames(counts))
    }else if( length( setdiff(markergenes,rownames(counts)) ) == 0 ){
      tab_ = matrix(0,nrow = nrow(tab),ncol = 1) %>% as.data.frame()
      rownames(tab_) = rownames(tab)
      colnames(tab_) = "the_null"
    }else{
      tab_ = the.count[setdiff(markergenes,rownames(counts)),]%>% as.matrix() %>% t()  %>% as.data.frame()
    }
    
    final_tab = cbind(tab,tab_[rownames(tab),])
    colnames(final_tab) = c(colnames(tab),colnames(tab_)) 
    tab = final_tab
    
  }
  
  # markergenes = intersect(markergenes,rownames(counts))
  
  UMAP_data = seurat_Data@reductions$umap@cell.embeddings %>% as.data.frame()
  UMAP_data = cbind(UMAP_data,tab[rownames(UMAP_data),])
  
  res <- .make_hexbin_helper(UMAP_data[,c(1:2)], 50, use_dims = c(1, 2))
  out <- res$hexbin.matrix
  cID <- res$cID
  out <- as_tibble(out)
  
  for(the.col in 3:ncol(UMAP_data)){
    
    the.gene = colnames(UMAP_data)[the.col]
    x <- as.numeric(UMAP_data[,the.col])
    hh <- .make_hexbin_function(x, "mean", cID)
    library(tibble)
    out[,the.gene]<- hh
    
  }
  
  saveRDS(out,file = paste0("./data/Tax_expOpzData/Cancer/",cancer_final,".RDS"))
  message(cancer_final)
  
}



####hexbin  plot
demo_markergene <- "CD79A"

for(i in 4:ncol(out)){
  
  tab = out[,i]
  tab = scale(tab)
  out[,i] = tab
  
}

a <- .plot_hexbin(out, colour_by=demo_markergene, action="mean",
                  title="title", xlab="xlab", ylab="ylab")+
  ggplot2::theme_void()+
  ggplot2::theme(panel.background = element_blank(),
                 panel.border =element_blank(),
                 axis.line = element_blank(),
                 plot.title = element_blank(),
                 plot.background = element_blank(),
                 panel.grid = element_blank(),
                 panel.grid.major = element_blank(),
                 plot.margin = unit(c(3.5,3.5,3.5,3.5),'lines'),
                 text = element_text(size = 13))+
  ggtitle("")+
  coord_fixed(ratio = (max(out$x)-min(out$x))/(max(out$y)-min(out$y)))









