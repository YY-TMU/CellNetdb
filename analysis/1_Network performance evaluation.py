###   evaluate network performance
###   Malignant cell networks from four reference networks  recover DisGeNet genes


from network_evaluation_tools import data_import_tools as dit
from network_evaluation_tools import network_evaluation_functions as nef
from network_evaluation_tools import network_propagation as prop
import pandas as pd
import numpy as np

import sys
cancer_type=sys.argv[1]
network_id=sys.argv[2]


print cancer_type
print network_id

wd="/home/lizekun/projects/CellNetdb_analysis/Network_Evaluation/01.malignant_net_edge/"
network = dit.load_network_file(wd+cancer_type+"_"+network_id+".txt", verbose=True)


gene_set_path="/home/lizekun/projects/CellNetdb_analysis/Network_Evaluation/00.DisGeNET_cancer_paired_geneset/"

f = open(gene_set_path+cancer_type+".txt")
node_set_lines = f.read().splitlines()
f.close()
genesets = {node_set_lines[0]:set(node_set_lines[1:])}

# Calculate geneset sub-sample rate
genesets_p = nef.calculate_p(network, genesets)

# Determine optimal alpha for network (can also be done automatically by next step)
alpha = prop.calculate_alpha(network)
print 'alpha',alpha

# Calculate network kernel for propagation
kernel = nef.construct_prop_kernel(network, alpha=alpha, verbose=True)


kernel_out_wd="/home/lizekun/projects/CellNetdb_analysis/Network_Evaluation/02.propagation_score/"
kernel.to_csv(kernel_out_wd+cancer_type+"_"+network_id+".csv",index=True,header=True,quoting=False)


if len(network.edges) < 250000 :
    # Calculate the AUPRC values for each gene set
    AUPRC_values = nef.small_network_AUPRC_wrapper(kernel, genesets, genesets_p, n=50, cores=5, verbose=True)
    
    # Construct null networks and calculate the AUPRC of the gene sets of the null networks
    # We can use the AUPRC wrapper function for this
    null_AUPRCs = []
    for i in range(50):
        shuffNet = nef.shuffle_network(network, max_tries_n=50, verbose=True)
        shuffNet_kernel = nef.construct_prop_kernel(shuffNet, alpha=alpha, verbose=False)
        shuffNet_AUPRCs = nef.small_network_AUPRC_wrapper(shuffNet_kernel, genesets, genesets_p, n=50, cores=5, verbose=False)
        null_AUPRCs.append(shuffNet_AUPRCs)
        print 'shuffNet', repr(i+1), 'AUPRCs calculated'
else :
    # Calculate the AUPRC values for each gene set
    AUPRC_values = nef.large_network_AUPRC_wrapper(kernel, genesets, genesets_p, n=50, cores=5, verbose=True)
    
    # Construct null networks and calculate the AUPRC of the gene sets of the null networks
    # We can use the AUPRC wrapper function for this
    null_AUPRCs = []
    for i in range(50):
        shuffNet = nef.shuffle_network(network, max_tries_n=50, verbose=True)
        shuffNet_kernel = nef.construct_prop_kernel(shuffNet, alpha=alpha, verbose=False)
        shuffNet_AUPRCs = nef.large_network_AUPRC_wrapper(shuffNet_kernel, genesets, genesets_p, n=50, cores=5, verbose=False)
        null_AUPRCs.append(shuffNet_AUPRCs)
        print 'shuffNet', repr(i+1), 'AUPRCs calculated'


# Construct table of null AUPRCs
null_AUPRCs_table = pd.concat(null_AUPRCs, axis=1)
null_AUPRCs_table.columns = ['shuffNet'+repr(i+1) for i in range(len(null_AUPRCs))]


# Calculate performance metric of gene sets
network_performance = nef.calculate_network_performance_score(AUPRC_values, null_AUPRCs_table, verbose=True)
network_performance.name = 'Test Network'


# Calculate network performance gain over median null AUPRC
network_perf_gain = nef.calculate_network_performance_gain(AUPRC_values, null_AUPRCs_table, verbose=True)
network_perf_gain.name = 'Test Network'



# Network Performance
network_performance_metric_data = pd.concat([network_performance, network_perf_gain], axis=1)
network_performance_metric_data.columns = ['Network Performance', 'Network Performance Gain']

out_wd="/home/lizekun/projects/CellNetdb_analysis/Network_Evaluation/03.Network_Performance_output/"

network_performance_metric_data.to_csv(out_wd+cancer_type+"_"+network_id+".csv",index=True,header=True,quoting=False)







