# CellNetdb
CellNetdb is a database of a large-scale atlas of cell-type-specific interactome networks within tumor microenvironments. A total of 55 tumor scRNA-seq datasets, corresponding to 563 different patients and covering over two million cells in 44 tumor types, were incorporated for constructing these networks. We have provided all the main codes and parameters for constructing the cell-type-specific interactome networks. Moreover, detailed elucidations have also been provided regarding the execution of statistical calculations within each CellNet function module. The main processes of network function analysis were described below.
## Reconstruction of cell-type-specific interactome networks

We utilized the workflow of [SCINET](https://github.com/shmohammadi86/SCINET) implemented in [ACTIONet](https://github.com/shmohammadi86/ACTIONet) along with four widely-used reference interactome networks, namely STRING, HumanNet, ConsensusPathDB, and Reactome, to construct cell-type-specific interactome networks for each cell type.

## CellNet function modules

The data in Mutation module of each solid tumor type was collected from [COSMIC](https://cancer.sanger.ac.uk/cosmic). Then, users could browse or search somatic mutations in gene which are involved in the queried subnetwork.

To facilitate users in acquiring functional insights into the network, we conducted enrichment analysis on Gene Ontology (GO) and disease-associated gene sets. The GO gene sets were obtained from the [Gene Ontology database](https://release.geneontology.org/2022-06-15) (release 2022-06-15), while the disease-associated gene sets were sourced from the [DisGeNET database](https://www.disgenet.org/) (v7.0). The statistical significance of the enrichment of genes in GO terms or disease-associated gene sets within each queried local network was determined using the hypergeometric test.

In the Survival module, we gathered clinical data from several large-scale cohorts, including TCGA, MMRF and TARGET. Then, users could browse or search genes involved in the queried subnetwork whose expression level are associated with patients' overall survival.

In the Communication module, users could browse or search ligand-receptor pairs involved in the queried subnetwork. We only showed the significant interactions (*P* < 0.05). The Communication score in web table was the overall strength of cell-cell communication.

## Gene prioritization

We have implemented the random walk with restart (RWR) algorithm to prioritize interested genes based on the cell-type-specific interactome networks. Specifically, the random walk with restart is mathematically defined as follows:

<p align="center">${ p^{t+1}=(1-γ)Wp^t+γp^0  }$</p>

$W$ represents the column-normalized adjacency matrix of the network. The vector $p^t$ denotes the probability for the random walk to be at node $v$ at time $t$, while $p^0$ is the initial probability vector where only the seed genes have non-zero values. The restart probability, $γ$, is set to 0.5. By iteratively repeating the process until the difference between $p^t$ and $p^{t+1}$ falls below $10^{-10}$, we can numerically approximate the steady-state probability vector. Ultimately, this allows for the ranking of all genes in the network.

## Network performance evaluation

We evaluated the performance of various malignant cell networks generated by four reference networks in terms of their capacity to recover DisGeNET disease genes associated with specific tumor type. In this part, we utilized two network performance metrics, which were previously defined by [Huang et al.](https://doi.org/10.1016/j.cels.2018.03.001
), as follows:

<p align="center">${ Performance \ Score =\frac{AUPRC_{avg\ }-Median{(AUPRC_{null\ })}}{k \times MAD{(AUPRC_{null\ })}}}$</p>

<p align="center">${ Performance \ Gain  =\frac{AUPRC_{avg\ }-Median{(AUPRC_{null\ })}}{Median{(AUPRC_{null\ })}} \times 100 \%  }$</p>

The area under the precision-recall curve ($AUPRC$) used to evaluate the performance of a recovery task. $AUPRC$<sub>$avg$</sub> was calculated by averaging the $AUPRC$ values obtained from 50 sub-sampling iterations. Then, to establish a reference, a null distribution was created by repeating the above steps for 50 degree-preserving network shuffles, resulting in the creation of $AUPRC$<sub>$null$</sub>. The $MAD$ stands for "Median Absolute Deviation". $k$ is the reciprocal of the value corresponding to a probability of 0.75 in the cumulative distribution function, which is equal to 1.48.

## Topological specificity and transcriptional specificity analysis

To assess the functional application of cell-type-specific networks in understanding the context-specific role of genes, we utilized a metric called topological specificity ($topS$), previously introduced by [Mohammadi et al.](https://doi.org/10.1016/j.cels.2019.10.007). This metric allows for the direct quantification of a gene's influence in a network, beyond what is captured by connectivity and strength alone. We initially calculated the total strength of their local neighbors, denoted as $w^{(celltype)}(i)$, for each gene $i$. Subsequently, a random model was constructed to preserve the underlying network topology while uniformly reshuffling the edge weights. This ensemble of random networks allowed us to recomputed the strength of interactions, thereby enabling the generation of a distribution of gene neighborhood strengths for each gene. By utilizing the mean and standard deviation of each distribution ${μ^{(celltype)}_{R}(i)}$ and 

${σ^{(celltype)}_{R}(i)}$ , respectively; the topological specificity of each gene in a given cell-type-specific network can be defined as follows:


<p align="center">${ topS(i) =\frac{w^{(celltype)}(i)-μ^{(celltype)}_{R}(i)}{σ^{(celltype)}_{R}(i)}  }$</p>

Transcriptional specificity of genes pertains to their degree of specificity in expression within a particular cell type. To determine this, we employed the gene expression profile to calculate the average expression of various genes in a given cell type $x_{celltype}(i)$ and other cell types $x_{else}(i)$. By considering the variance of each group ${s^2_{celltype}(i)}$ and ${s^2_{else}(i)}$, the transcriptional specificity of each gene in a given cell type can be defined as:

$$
\normalsize    tranS(i) =\frac{x_{celltype}(i)-x_{else}(i)}{ \sqrt{\frac{s^2_{celltype}(i)}{n_{celltype}}  + \frac{s^2_{else}(i)}{n_{else}}  }  }
$$


## Network connectivity analysis

We ranked the top 500 prognostic genes for each tumor type based on their statistical significance using clinical data from the [TCGA project](https://www.cancer.gov/ccg/research/genome-sequencing/tcga). The normalized within-group connectivity of each cancer prognostic signature in all cell-type-specific networks was calculated. Furthermore, $10,000$ random gene sets, each containing the same number of genes as the test gene set, were selected to generate a null model. The statistical significance was assessed by measuring the rank of the observed within-group connectivity within the null distribution.


## Network similarity evaluation

We employed two distinct metrics, namely shared-edge similarity and topology similarity, to assess the degree of similarity between networks. Initially, shared nodes were identified between any pair of networks. To quantify the shared-edge similarity, the edges connecting these nodes were extracted from both networks, resulting in the creation of subgraphs for each network. The shared-edge similarity was subsequently determined by calculating the Spearman correlation coefficient between the weights assigned to the shared edges in the respective subgraphs of both networks. To evaluate the topology similarity, the Spearman correlation coefficient was computed for the transformed topological specificity ($topS_{transf}$) across all shared nodes. The transformation function used for $topS_{transf}$ was defined as:

$$
\normalsize   topS_{transf}(i) =\frac{1}{1+e^{-topS(i)}}
$$





