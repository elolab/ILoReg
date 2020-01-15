# ILoReg

## Introduction



`ILoReg` is a tool for high-resolution *de novo* cell population identification from single-cell RNA-seq (scRNA-seq) data. High-resolution refers to the ability of `ILoReg` to identify subpopulations with subtle transcriptomic differences. In our study [1](https://gitlab.utu.fi/pajosm/iloreg), we showed that `ILoReg` identitied, by both unsupervised clustering and visually, immune cell populations that other scRNA-seq data analysis pipelines had difficulties to find.

The figure below illustrates the workflows of `ILoReg` and a conventional pipeline that applies feature selection prior to dimensionality reduction by principal component analysis (PCA), e.g. Seurat.

![*Figure: Analysis workflows of ILoReg and a feature-selection based approach*](vignettes/figure.png)



In contrast to most scRNA-seq data analysis pipelines, `ILoReg` does not reduce the dimensionality of the gene expression matrix by feature selection. Instead, it performs probabilistic feature extraction, in which the **iterative clustering projection (ICP)** clustering algorithm is run *L* times, yielding *L* probability matrices that contain probabilities of each of the *N* cells belonging to the *k* clusters. ICP is a new clustering algorithm that iteratively seeks a clustering with *k* clusters that maximizes the adjusted Rand index (ARI) between the clustering and its projection by L1-regularized logistic regression. The *L* probability matrices are then merged into a joint probability matrix and subsequently transformed by PCA into a lower dimensional matrix (consensus matrix). The final clustering step is performed using hierarhical clustering by the Ward's method, after which the user can efficiently (~1 s with 3,000 cells) extract a clustering with *K* consensus clusters. Two-dimensional visualization is supported using two popular nonlinear dimensionality reduction methods: *t*-distributed stochastic neighbor embedding (t-SNE) and uniform manifold approximation and projection (UMAP). Additionally, `ILoReg` provides user-friendly functions that enable identification of differentially expressed genes and visualization of gene expression.

## Installation

The latest version of `ILoReg` can be downloaded from GitHub using devtools R package.

```R
library(devtools)

creds = git2r::cred_ssh_key("~/../.ssh/id_rsa.pub",
                            "~/../.ssh/id_rsa")
devtools::install_git("gitlab@gitlab.utu.fi:pajosm/iloreg.git",
                      credentials = creds)

```

## Example

Please follow this [link](https://gitlab.utu.fi/pajosm/iloreg) to an example, where a peripheral blood mononuclear cell (PBMC) dataset was analyzed using `ILoReg`.

## Contact information

If you have questions related to `ILoReg`, please contact us in [GitHub](https://gitlab.utu.fi/pajosm/iloreg). 

## References

1. Johannes Smolander, Sini Junttila, Mikko S Venäläinen, Laura L Elo. "ILoReg enables high-resolution cell population identification from single-cell RNA-seq data". BioRxiv (2020).
