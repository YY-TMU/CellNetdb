import patsy
import pandas as pd
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
import matplotlib.pyplot as plt
sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80)
# %matplotlib inline

import os
import torch
from multiprocessing import cpu_count

cpu_num = 10
os.environ ['OMP_NUM_THREADS'] = str(cpu_num)
os.environ ['OPENBLAS_NUM_THREADS'] = str(cpu_num)
os.environ ['MKL_NUM_THREADS'] = str(cpu_num)
os.environ ['VECLIB_MAXIMUM_THREADS'] = str(cpu_num)
os.environ ['NUMEXPR_NUM_THREADS'] = str(cpu_num)
torch.set_num_threads(cpu_num)


import sys
cancer_used=sys.argv[1]



adata_raw = anndata.read_h5ad("/data1/lizekun/CellNetdb/add_analysis/01.h5ad/" + cancer_used + ".h5ad")
print(adata_raw)


adata_combat=adata_raw.copy()
sc.pp.normalize_total(adata_combat, target_sum=1e4)
sc.pp.log1p(adata_combat)

# run combat
sc.pp.combat(adata_combat, key='BatchID',covariates=['nCount_RNA',"percent.mt"])

# ## save 
# anndata.AnnData.write_h5ad(adata_combat,"/data1/lizekun/CellNetdb/add_analysis/Combat/BatchID/03.h5ad/" + cancer_used +"_combat_out.h5ad")


sc.pp.highly_variable_genes(adata_combat)
print("Highly variable genes: %d"%sum(adata_combat.var.highly_variable))
sc.pl.highly_variable_genes(adata_combat)


sc.pp.pca(adata_combat, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata_combat)
sc.tl.umap(adata_combat)

# fig, axs = plt.subplots(1, 2, figsize=(10,4),constrained_layout=True)
# sc.pl.umap(adata_combat, color="BatchID", title="Combat Corrected umap", ax=axs[0] , show=False)
# sc.pl.umap(adata_raw, color="BatchID", title="CCA umap", ax=axs[1], show=False)
# 
# sc.pl.umap(adata_combat, color="Cell_TypeID", title="Combat Corrected umap" , show=False)
# sc.pl.umap(adata_raw, color="Cell_TypeID", title="CCA umap", show=False)



emd_data = pd.DataFrame(adata_combat.obsm['X_umap'])
emd_data.index = adata_combat.obs.index
emd_data.to_csv("/data1/lizekun/CellNetdb/add_analysis/Combat/BatchID/01.cancer_emd_data/"+cancer_used+".csv",index=True,header=True,quoting=False)



ASW_celltype = sklearn.metrics.silhouette_score(labels = adata_combat.obs["Cell_TypeID"],X = adata_combat.obsm["X_umap"])
ASW_celltype = (ASW_celltype+1)/2

ASW_batch_Batch = sklearn.metrics.silhouette_score(labels = adata_combat.obs["BatchID"],X = adata_combat.obsm["X_umap"])
ASW_batch_Batch = (ASW_batch_Batch+1)/2

if "DonorID" in adata_combat.obs.columns :
    ASW_batch_Donor = sklearn.metrics.silhouette_score(labels = adata_combat.obs["DonorID"],X = adata_combat.obsm["X_umap"])
    ASW_batch_Donor = (ASW_batch_Donor+1)/2  
else:
        ASW_batch_Donor = None



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
    sc.pp.neighbors(adata,random_state=0)
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


n_celltype=adata_combat.obs["Cell_TypeID"].nunique()
cancer_ARI,cancer_NMI=calulate_ari_nmi(adata_combat,n_celltype)


data = [['ARI',cancer_ARI],['NMI',cancer_NMI],['ASW_celltype',ASW_celltype],['ASW_batch_Batch',ASW_batch_Batch],['ASW_batch_Donor',ASW_batch_Donor]]
df = pd.DataFrame(data,columns=['metrics','value']) 




## save 
df.to_csv("/data1/lizekun/CellNetdb/add_analysis/Combat/BatchID/02.metric/four_metrics/"+cancer_used+".csv",index=False,header=True,quoting=False)
anndata.AnnData.write_h5ad(adata_combat,"/data1/lizekun/CellNetdb/add_analysis/Combat/BatchID/03.h5ad/" + cancer_used +".h5ad")













