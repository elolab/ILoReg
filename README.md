# ILoReg

## Introduction

`ILoReg` is a novel tool for cell population identification from single-cell RNA-seq (scRNA-seq) data. In our study [[1]](https://doi.org/10.1093/bioinformatics/btaa919), we showed that `ILoReg` was able to identify, by both unsupervised clustering and visually, rare cell populations that other scRNA-seq data analysis pipelines were unable to identify. 


The figure below illustrates the workflows of `ILoReg` and a typical pipeline that applies feature selection prior to dimensionality reduction by principal component analysis (PCA).

![*Figure: Analysis workflows of ILoReg and a feature-selection based approach*](vignettes/figure.png)


In contrast to most scRNA-seq data analysis pipelines, `ILoReg` does not reduce the dimensionality of the gene expression matrix by feature selection. Instead, it performs probabilistic feature extraction using **iterative clustering projection (ICP)**, yielding a probability matrix, which contains probabilities of each of the *N* cells belonging to the *k* clusters. ICP is a novel machine learning algorithm that iteratively seeks a clustering with *k* clusters that maximizes the adjusted Rand index (ARI) between the clustering and its projection by L1-regularized logistic regression. In the ILoReg consensus approach, ICP is run *L* times and the *L* probability matrices are merged into a joint probability matrix and subsequently transformed by principal component analysis (PCA) into a lower dimensional matrix (consensus matrix). The final clustering step is performed using hierarhical clustering by the Ward's method, after which the user can extract a clustering with *K* consensus clusters. Two-dimensional visualization is supported using two popular nonlinear dimensionality reduction methods: *t*-distributed stochastic neighbor embedding (t-SNE) and uniform manifold approximation and projection (UMAP). Additionally, ILoReg provides user-friendly functions that enable identification of differentially expressed (DE) genes and visualization of gene expression.

## Installation

The latest version of `ILoReg` can be downloaded from GitHub using the devtools R package.

```R

devtools::install_github("elolab/ILoReg")

```

## Example

Please follow this [link](https://github.com/elolab/ILoReg/blob/master/vignettes/ILoReg.Rmd) to an example, in which a peripheral blood mononuclear cell (PBMC) dataset is analyzed using `ILoReg`. In [Bioconductor](https://bioconductor.org/packages/release/bioc/html/ILoReg.html) the vignette can be accessed in a readable format. 

## Contact information

If you have questions related to `ILoReg`, please contact us [here](https://github.com/elolab/ILoReg/issues). 

## References

1. Johannes Smolander, Sini Junttila, Mikko S Venäläinen, Laura L Elo. " ILoReg: a tool for high-resolution cell population identification from single-cell RNA-seq data". Bioinformatics, Volume 37, Issue 8, 15 April 2021, Pages 1107–1114, [https://doi.org/10.1093/bioinformatics/btaa919](https://doi.org/10.1093/bioinformatics/btaa919).
