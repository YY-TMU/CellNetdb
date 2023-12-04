####   Different  batch correction methonds 

import numpy as np
import anndata
import scanpy as sc
import scvi
import scipy
from scipy.sparse import csr_matrix
from rich import print
from scvi.model.utils import mde
from scvi_colab import install
from scib_metrics.benchmark import Benchmarker
from scvi.model.utils import mde
import sklearn
from sklearn.metrics.cluster import pair_confusion_matrix
from sklearn.metrics import silhouette_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score


import sys

cancer_used=sys.argv[1]

path_to_anndata=["/data1/lizekun/CellNetdb/add_analysis/scvi/01.h5ad/",cancer_used,".h5ad"]
adata = anndata.read_h5ad("".join(path_to_anndata))
adata.X=scipy.sparse.csr_matrix(adata.X)
adata.layers["counts"] = adata.X.copy()
adata.raw = adata
sc.pp.highly_variable_genes(
  adata,
  flavor="seurat_v3",
  n_top_genes=2000,
  layer="counts",
  batch_key="BatchID",
  subset=True,
)

TF_CPP_MIN_LOG_LEVEL=0
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="BatchID")
# scvi.model.SCVI.setup_anndata(adata, layer="counts", categorical_covariate_keys=["BatchID","DonorID"])
vae = scvi.model.SCVI(adata)

vae.train()

# model_path = ["/data1/lizekun/CellNetdb/add_analysis/scvi/02.model_output/" ,cancer_used]
# vae.save("".join(model_path))
# 
# new_h5ad_path=["/data1/lizekun/CellNetdb/add_analysis/scvi/02.model_output/" ,cancer_used,"/",cancer_used,".h5ad"]
# anndata.AnnData.write_h5ad(adata,"".join(new_h5ad_path))


adata.obsm["X_scVI"] = vae.get_latent_representation()
adata.obsm["X_normalized_scVI"] = vae.get_normalized_expression()

##save 
addmodel_h5ad_path=["/data1/lizekun/CellNetdb/add_analysis/scvi/03.add_model/" ,cancer_used,".h5ad"]
anndata.AnnData.write_h5ad(adata,"".join(addmodel_h5ad_path))





#scVI_emd_data = pd.DataFrame(adata.obsm['X_scVI'])
#scVI_emd_data.index = adata.obs.index
#scVI_emd_data.to_csv("/data1/lizekun/CellNetdb/add_analysis/scvi/00.adj_rds/cancer_scVI_emd_data/"+cancer_used+".csv",index=True,header=True,quoting=False)




ASW_celltype = sklearn.metrics.silhouette_score(labels = adata.obs["Cell_TypeID"],X = adata.obsm["X_scVI"])
ASW_celltype = (ASW_celltype+1)/2

ASW_batch_Batch = sklearn.metrics.silhouette_score(labels = adata.obs["BatchID"],X = adata.obsm["X_scVI"])
ASW_batch_Batch = (ASW_batch_Batch+1)/2

if "DonorID" in adata.obs.columns :
    ASW_batch_Donor = sklearn.metrics.silhouette_score(labels = adata.obs["DonorID"],X = adata.obsm["X_scVI"])
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
    sc.pp.neighbors(adata,random_state=0,use_rep="X_scVI")
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


## save 
df.to_csv("/home/lizekun/projects/CellNetdb_analysis/scVI/metrics/02.scVI_out/four_metrics/"+cancer_used+".csv",index=False,header=True,quoting=False)




