# ILoReg

## Introduction

ILoReg (abbr. for **I**terative **Lo**gistic **Reg**ression) is an R package for high-precision cell type identification from single-cell RNA-seq (scRNA-seq) data. 
High-precision refers to the ability of ILoReg to identify subpopulations with subtle transcriptomic differences. 
In our study ([Article title](https://gitlab.utu.fi/pajosm/iloreg)), we showed that ILoReg was able to identify, by both unsupervised clustering and visually, immune cell types that other scRNA-seq data analysis pipelines were not able to find. 
Moreover, ILoReg can identify subpopulations that are differentiable by only a single gene.

The figure below depicts the workflows of ILoReg and a feature selection -based pipeline, which is commonly used by many cell type identification methods, e.g. Seurat.

![*Figure: Analysis workflows of ILoReg and a feature-selection based approach*](vignettes/figure.png)


Unlike most scRNA-seq data analysis pipelines, ILoReg does not reduce the dimensionality of the gene expression matrix by feature selection. 
Instead, it performs probabilistic feature extraction, in which the Iterative Clustering Projection (ICP) clustering algorithm is run *L* times, which yields 
*L* *k*-dimensional probability matrices that contain the new features. ICP is a novel clustering algorithm that iteratively seeks a clustering of size *k* 
that maximizes the adjusted Rand index (ARI) between the clustering and its projection by L1-regularized logistic regression. 
The *L* probability matrices are then merged and transformed by the principal component analysis (PCA) into a lower dimension *p*. 
The second and final clustering step is performed using hierarhical clustering by the Ward's method, from which the user can efficiently (~1 s with 3,000 cells) 
select a clustering of size *K*. Two-dimensional visualization is supported using two popular nonlinear dimensionality reduction methods: 
*t*-Distributed Stochastic Neighbor Embedding (t-SNE) and Uniform Manifold Approximation and Projection (UMAP).

## Installation

The latest version of ILoReg can be downloaded from GitHub using devtools R package.

```R
library(devtools)

creds = git2r::cred_ssh_key("~/../.ssh/id_rsa.pub",
                            "~/../.ssh/id_rsa")
devtools::install_git("gitlab@gitlab.utu.fi:pajosm/iloreg.git",
                      credentials = creds)

```

## Contact information

If you find bugs from ILoReg or have suggestions on how to improve our pipeline, please contact us in [GitHub](https://gitlab.utu.fi/pajosm/iloreg). 

