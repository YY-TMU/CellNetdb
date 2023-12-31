####   Compare different  batch correction methonds 


####   We calculate three metrics (ARI, NMI, ASW_celltype) were adopted to evaluate the clustering performance and 
####   three metrics (iLISI, BatchKL, ASW_batch) were used to evaluate the ability to remove batch effect.


library(tidyverse)
library(lisi)



##BatchKL    The lower value denotes the better mixing performance.
BatchKL=function(df,dimensionData=NULL,replicates=200,n_neighbors=100,n_cells=100,batch=NULL){
  #entropy of batch mixiing
  #replicates is the number of boostrap times
  #df is metadata which includes BatchID
  #dimensionData is multidimensional coordinates of all cells 
  #n_neighbors is the number of nearest neighbours of cell(from all batchs)
  #n_cells is the number of randomly picked cells
  set.seed(1)
  if (is.null(dimensionData)){
    umapdata=as.matrix(df[,c("UMAP_1","UMAP_2")])
  }else{
    umapdata=as.matrix(dimensionData)
  }
  batchdata=factor(as.vector(df[,batch]))
  table.batchdata=as.matrix(table(batchdata))[,1]
  tmp00=table.batchdata/sum(table.batchdata)#proportation of population
  n=dim(df)[1]
  KL=sapply(1:replicates,function(x){
    bootsamples=sample(1:n,n_cells)
    #nearest=nn2(umapdata,umapdata[bootsamples,],k=n_neighbors)
    nearest=nabor::knn(umapdata,umapdata[bootsamples,],k=min(5*length(tmp00),n_neighbors))
    KL_x=sapply(1:length(bootsamples),function(y){
      id=nearest$nn.idx[y,]
      tmp=as.matrix(table(batchdata[id]))[,1]
      tmp=tmp/sum(tmp)
      return(sum(tmp*log2(tmp/tmp00),na.rm = T))
    })
    return(mean(KL_x,na.rm = T))
  })
  return(mean(KL))#
}


########higher iLISI indicates better performance for batch mixing.
CalLISI=function(emb,meta){
  if('DonorID' %in% colnames(meta)){
    lisi_index <- lisi::compute_lisi(emb, meta, c('Cell_TypeID', 'BatchID',"DonorID"))
    clisi = median(lisi_index$Cell_TypeID)
    ilisi_Batch = median(lisi_index$BatchID)
    ilisi_Donor = median(lisi_index$DonorID)
    return(c(clisi,ilisi_Batch,ilisi_Donor))
  }else{
    lisi_index <- lisi::compute_lisi(emb, meta, c('Cell_TypeID', 'BatchID'))
    clisi = median(lisi_index$Cell_TypeID)
    ilisi_Batch = median(lisi_index$BatchID)
    return(c(clisi,ilisi_Batch))
  }
}







