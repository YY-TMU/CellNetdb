####   Different  batch correction methonds 

import os
import anndata
import scanpy as sc 
import matplotlib
import numpy as np
import pandas as pd 
import seaborn as sns
# %matplotlib inline
import matplotlib.pyplot as plt
import plotly
import scDML 
print(scDML.__version__)
import sys
print(sys.version)
from scDML import scDMLModel
from scDML.utils import print_dataset_information

import sys
cancer_used=sys.argv[1]


adata_raw = anndata.read_h5ad("/data1/lizekun/CellNetdb/add_analysis/01.h5ad/" + cancer_used + ".h5ad")
print(adata_raw)
print_dataset_information(adata_raw,batch_key="BatchID",celltype_key="Cell_TypeID")


scdml=scDMLModel(save_dir="/data1/lizekun/CellNetdb/add_analysis/scDML/01.model/"+cancer_used+"/",verbose=True)

##  preprocessing 
adata=adata_raw.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,n_top_genes=2000,subset=True)
sc.pp.regress_out(adata, ['nCount_RNA', 'percent.mt'])
sc.pp.scale(adata,max_value=10)
sc.tl.pca(adata,n_comps=100)
sc.pp.neighbors(adata,random_state=0)
sc.tl.louvain(adata,resolution=3.0,key_added="init_cluster")

print(adata)

from scDML.utils import plotDendrogram,plotHeatMap
# convert adata to training data for neural network
scdml.convertInput(adata,batch_key="BatchID")
_,_,cor,_ = scdml.calculate_similarity()

fig = plt.figure(figsize=(15,10))
Z=plotDendrogram(scdml.cor_matrix.copy(),scdml.nn_matrix.copy(),
                           scdml.merge_df["init_cluster"].value_counts().values.copy())
hf=plotHeatMap(scdml.cor_matrix.copy(),Z)
plt.show()

from scDML.utils import plotSankey
ncelltype=adata.obs["Cell_TypeID"].nunique()

merge_df=scdml.merge_cluster(ncluster_list=[ncelltype],merge_rule="rule2")
merge_df["celltype"]=adata.obs["Cell_TypeID"].values

cols=["init_cluster"]+[str(ncelltype)]+["celltype"]
fig=plotSankey(merge_df,cat_cols=cols,value_cols='value',title="Sanky plot")
plotly.offline.iplot(fig)


scdml.build_net()


# delete trained model to retrain the model
if os.path.isfile(os.path.join(scdml.save_dir,"scDML_model.pkl")):
    os.remove(os.path.join(scdml.save_dir,"scDML_model.pkl"))
    
embedding=scdml.train(expect_num_cluster=ncelltype)

adata.obsm["X_emb"]=embedding
adata.obs["reassign_cluster"]=scdml.train_label.astype(int).astype(str)
adata.obs["reassign_cluster"]=adata.obs["reassign_cluster"].astype("category")



# plot loss
fig=plt.figure()
plt.title("scDML loss")
plt.plot(range(1,len(scdml.loss)+1),scdml.loss,c="r")
plt.xlabel("training epoch")
plt.ylabel("number of mined hard triplet")
plt.show()

## or  plt.plot(range(1,len(scdml.loss)+1),scdml.loss)

# scdml.integrate(adata,batch_key="BatchID",ncluster_list=[ncelltype],
#                expect_num_cluster=ncelltype,merge_rule="rule2")
               

sc.pp.neighbors(adata,use_rep="X_emb",random_state=0)
sc.tl.umap(adata)
sc.pl.umap(adata,color=["Cell_TypeID","BatchID"])
sc.pl.umap(adata,color=["Cell_TypeID","reassign_cluster"],legend_loc="on data",legend_fontsize="xx-large")


## step1 : 提取scDML坐标

scDML_emd_data = pd.DataFrame(adata.obsm['X_emb'])
scDML_emd_data.index = adata.obs.index
scDML_emd_data.to_csv("/data1/lizekun/CellNetdb/add_analysis/scDML/02.cancer_scDML_emd_data/"+cancer_used+".csv",index=True,header=True,quoting=False)



## step2 :   ##evaluation
import sklearn
from sklearn.metrics import silhouette_score
from sklearn.metrics.cluster import adjusted_rand_score,normalized_mutual_info_score

cancer_ARI= adjusted_rand_score(adata.obs["Cell_TypeID"].astype(str), adata.obs["reassign_cluster"])
cancer_NMI= normalized_mutual_info_score(adata.obs["Cell_TypeID"].astype(str), adata.obs["reassign_cluster"])
print("ARI={}".format(cancer_ARI))
print("NMI={}".format(cancer_NMI))


ASW_celltype = sklearn.metrics.silhouette_score(labels = adata.obs["Cell_TypeID"],X = adata.obsm["X_emb"])
ASW_celltype = (ASW_celltype+1)/2

ASW_batch_Batch = sklearn.metrics.silhouette_score(labels = adata.obs["BatchID"],X = adata.obsm["X_emb"])
ASW_batch_Batch = (ASW_batch_Batch+1)/2

if "DonorID" in adata.obs.columns :
    ASW_batch_Donor = sklearn.metrics.silhouette_score(labels = adata.obs["DonorID"],X = adata.obsm["X_emb"])
    ASW_batch_Donor = (ASW_batch_Donor+1)/2  
else:
        ASW_batch_Donor = None


data = [['ARI',cancer_ARI],['NMI',cancer_NMI],['ASW_celltype',ASW_celltype],['ASW_batch_Batch',ASW_batch_Batch],['ASW_batch_Donor',ASW_batch_Donor]]
df = pd.DataFrame(data,columns=['metrics','value']) 




##save h5ad
anndata.AnnData.write_h5ad(adata,"/data1/lizekun/CellNetdb/add_analysis/scDML/02.h5ad/" +cancer_used+".h5ad")
df.to_csv("/home/lizekun/projects/CellNetdb_analysis/scDML/metrics/four_metrics/"+cancer_used+".csv",index=False,header=True,quoting=False)




