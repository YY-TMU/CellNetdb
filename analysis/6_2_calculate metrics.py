####   Compare different  batch correction methonds 


####   We calculate three metrics (ARI, NMI, ASW_celltype) were adopted to evaluate the clustering performance and 
####   three metrics (iLISI, BatchKL, ASW_batch) were used to evaluate the ability to remove batch effect.


import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import sklearn
from sklearn.metrics.cluster import pair_confusion_matrix
from sklearn.metrics import silhouette_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score

import sys
cancer_used=sys.argv[1]
print(cancer_used)


path_to_anndata=["/data1/lizekun/CellNetdb/add_analysis/01.h5ad/",cancer_used,".h5ad"]           #####  save raw Seurat count and metadata to h5ad data
adata = anndata.read_h5ad("".join(path_to_anndata))

ASW_celltype = sklearn.metrics.silhouette_score(labels = adata.obs["Cell_TypeID"],X = adata.obsm["X_umap"])
ASW_celltype = (ASW_celltype+1)/2

ASW_batch_Batch = sklearn.metrics.silhouette_score(labels = adata.obs["BatchID"],X = adata.obsm["X_umap"])
ASW_batch_Batch = (ASW_batch_Batch+1)/2

if "DonorID" in adata.obs.columns :
    ASW_batch_Donor = sklearn.metrics.silhouette_score(labels = adata.obs["DonorID"],X = adata.obsm["X_umap"])
    ASW_batch_Donor = (ASW_batch_Donor+1)/2  
else:
        ASW_batch_Donor = None

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)



# sklearn ari bug
def ari(labels_true,labels_pred): 
    '''safer implementation of ari score calculation'''
    (tn, fp), (fn, tp) = pair_confusion_matrix(labels_true, labels_pred)
    tn=int(tn)
    tp=int(tp)
    fp=int(fp)
    fn=int(fn)

    # Special cases: empty data or full agreement
    if fn == 0 and fp == 0:
        return 1.0

    return 2. * (tp * tn - fn * fp) / ((tp + fn) * (fn + tn) +
                                       (tp + fp) * (fp + tn))


# Find optimal resolution given ncluster
def find_resolution(adata_, n_clusters, random):
    adata = adata_.copy()
    obtained_clusters = -1
    iteration = 0
    resolutions = [0., 1000.]
    while obtained_clusters != n_clusters and iteration < 50:
        current_res = sum(resolutions)/2
        sc.tl.louvain(adata, resolution = current_res, random_state = random)
        labels = adata.obs['louvain']
        obtained_clusters = len(np.unique(labels))

        if obtained_clusters < n_clusters:
            resolutions[0] = current_res
        else:
            resolutions[1] = current_res
        iteration = iteration + 1
    return current_res


def calulate_ari_nmi(adata,n_clusters):
    sc.pp.neighbors(adata,random_state=0,n_pcs=cancer_nPC)
    reso=find_resolution(adata,n_clusters,0)
    sc.tl.louvain(adata,reso,random_state=0)
    sc.tl.umap(adata)
    if(adata.X.shape[1]==2):
        adata.obsm["X_emb"]=adata.X
#         sc.pl.embedding(adata, basis='emb', color = ['louvain'], wspace = 0.5)
#     else:
#         sc.pl.umap(adata,color=["louvain"])

    cancer_ARI= ari(adata.obs["Cell_TypeID"].astype(str), adata.obs["louvain"])
    cancer_NMI= normalized_mutual_info_score(adata.obs["Cell_TypeID"].astype(str), adata.obs["louvain"])
    print("louvain clustering result(resolution={}):n_clusters={}".format(reso,n_clusters))
    print("ARI:",cancer_ARI)
    print("NMI:",cancer_NMI)
    return cancer_ARI,cancer_NMI

n_celltype=adata.obs["Cell_TypeID"].nunique()
cancer_ARI,cancer_NMI=calulate_ari_nmi(adata,n_celltype)




data = [['ARI',cancer_ARI],['NMI',cancer_NMI],['ASW_celltype',ASW_celltype],['ASW_batch_Batch',ASW_batch_Batch],['ASW_batch_Donor',ASW_batch_Donor]]
df = pd.DataFrame(data,columns=['metrics','value']) 

df.to_csv("/home/lizekun/projects/CellNetdb_analysis/metrics/01.CCA_out/four_metrics/"+cancer_used+".csv",index=False,header=True,quoting=False)



    
   
   

