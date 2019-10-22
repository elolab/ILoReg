#' @title Iterative Logistic Regression (ILoReg) consensus clustering
#'
#' @description R package that enables supervised learning -based clustering through L1-regularized (LASSO) logistic regression.
#' Consensus clustering is performed by running ILoReg L times. PCA is used to reduce dimensionality of the
#' consensus probability matrix. The Ward's method is used to perform the final clustering.
#' The silhouette method can be used to choose the optimal number of clusters (K). t-SNE or UMAP can be performed
#' onto the PCA-transformed data to visualize the results. The package provides functions for differential expression analysis
#' and visualization of gene.
#'
#' @author Johannes Smolander <pajosm@utu.fi>, Mikko Ven?l?inen, Sini Junttila, Laura L. Elo
#'
#' @docType package
#' @name ILoReg-package
#' @aliases ILoReg
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' The iloreg Class
#'
#' iloreg class description...
#'
#' @name iloreg-class
#' @rdname iloreg-class
#' @importClassesFrom Matrix dgCMatrix
#' @exportClass iloreg
#'

setOldClass("hclust")
setClass("iloreg", representation(normalized.data = "dgCMatrix",
                                  k = "numeric",
                                  d = "numeric",
                                  L = "numeric",
                                  r = "numeric",
                                  C = "numeric",
                                  type = "character",
                                  max.number.of.iterations = "numeric",
                                  consensus.probability = "list",
                                  metrics = "list",
                                  rotated.consensus.probability = "matrix",
                                  hc = "hclust",
                                  umap.embeddings = "matrix",
                                  tsne.embeddings = "matrix",
                                  threads = "numeric",
                                  number.of.pcs = "numeric",
                                  scale.pca = "logical",
                                  silhouette.information = "numeric",
                                  K.optimal = "numeric",
                                  clustering.optimal = "factor",
                                  K.manual = "numeric",
                                  clustering.manual = "factor"),prototype=prototype(hc=structure(list(), class="hclust")))


#' Create an ILoReg object
#' @param normalized.data An object of class dgCMatrix, matrix or data.frame. A normalized gene expression matrix with genes in rows and cells in columns. If the object is of class matrix or data.frame, it is transformed into a LogNormalized expression matrix. The non-expressing genes are removed when the CreateILoRegObject function is called.
#'
#' @importFrom Matrix rowSums
#' @export
#'
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#'
CreateILoRegObject <- function(normalized.data)
{

  normalized.data_class <- class(normalized.data)

  if (normalized.data_class=="matrix")
  {
    normalized.data <- Matrix(normalized.data,sparse = TRUE)
    cat("Transforming object of class 'matrix' to 'dgCMatrix'\n")
  }

  if (normalized.data_class=="data.frame")
  {
    normalized.data <- Matrix(as.matrix(normalized.data),sparse = TRUE)
    cat("Transforming object of class 'data.frame' to 'dgCMatrix'\n")
  }

  if (normalized.data_class=="dgCMatrix")
  {
    cat("Data already of class 'dgCMatrix'\n")
  }

  cat("Removing non-expressing genes\n")
  genes_before_filtering <- nrow(normalized.data)
  normalized.data <- normalized.data[rowSums(normalized.data)!=0,]
  genes_after_filtering <- nrow(normalized.data)
  cat(genes_after_filtering,"/",genes_before_filtering," left\n")

  iloreg_object <- new(Class="iloreg",normalized.data = normalized.data)

  return(iloreg_object)
}

setGeneric("RunParallelICP", function(iloreg.object=NULL,k=15,d=0.3,L=200,r=5,C=0.3,type="L1",max.number.of.iterations=200,threads=0,seed=1){
  standardGeneric("RunParallelICP")
})

#' @title ICP consensus method
#'
#' @rdname RunParallelICP
#' @name RunParallelICP
#'
#' @description
#' Enables running L ICP runs in parallel.
#'
#' @details
#' @param iloreg.object object of class 'iloreg'
#' @param k A positive integer greater or equal to 2, which denotes the number of clusters in iterative logistic regression. (default 15)
#' @param d A numeric greater than 0 and smaller than 1 that determines how many cells 'n' are down- or oversampled from each cluster into the training data. (d in n=N/k*d), where N is the total number of cells, k is the number of clusters in ICP. (default 0.3)
#' @param L A positive integer greater than 1. Number of ICP runs. Contraining L has not been studied.  (Default 200)
#' @param r A positive integer that denotes the maximum number of reiterations performed until the ICP algorithm stops. (default 5). Increasing recommended with a significantly larger sample size (tens of thousands).
#' @param C A positive real number, the cost of constraints violation in the L1-regularized logistic regression model from the LIBLINEAR library. (default 0.3). Decreasing leads to more stringent feature selection, i.e. less genes are selected that are used to build the classifier.
#' @param type "L1" or "L2". L2-regularization was not investigated in the manuscript. (default "L1")
#' @param threads A positive integer that specifies how many logical processors (threads) should be used. Use threads=1 to disable parallelism and threads=0 to use all available core minus one. (default 0)
#' @param max.number.of.iterations A positive integer that denotes the maximum number of iterations performed until the algorithm ends. (default 200)
#' @param seed A positive integer that specifies the random seed to be set before randomly generating random seeds for each parallel task. (default 1)
#' @return iloreg object
#' @keywords iterative logistic regression ILoReg consensus clustering
#' @importFrom parallel makeCluster
#' @importFrom parallel detectCores
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @import Matrix
#' @import tictoc
#' @import aricode
#' @import LiblineaR
#' @import SparseM
#' @export
#' @examples
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
setMethod("RunParallelICP", "iloreg", function(iloreg.object, k,d,L,r,C,type,max.number.of.iterations,threads,seed){

  if (class(iloreg.object) != "iloreg") {
    stop("iloreg.object must of of class 'iloreg'")
  }

  if (!is.numeric(k) | k < 2 | k%%1!=0)
  {
    stop("k must be positive integer and greater than 1")
  } else {
    iloreg.object@k <- k
  }

  if (!is.numeric(d) | d >= 1 | d <= 0)
  {
    stop("d must be numeric and in range of (0,1)")
  } else {
    iloreg.object@d <- d
  }

  if (!is.numeric(L) | L <= 0 | L%%1!=0)
  {
    stop("L must be positive integer and greater than 0")
  } else {
    iloreg.object@L <- L
  }

  if (!is.numeric(r) | r <= 0 | r%%1!=0)
  {
    stop("r must be positive integer and greater than 0")
  } else {
    iloreg.object@r <- r
  }

  if (!is.numeric(C) | C <= 0)
  {
    stop("C must be numeric and greater than 0")
  } else {
    iloreg.object@C <- C
  }

  if (!is.character(type) | (type!="L1" & type!="L2"))
  {
    stop("type must be 'L1' or 'L2'")
  } else {
    iloreg.object@type <- type
  }

  if (!is.numeric(max.number.of.iterations) | max.number.of.iterations <= 0 | max.number.of.iterations%%1!=0)
  {
    stop("max.number.of.iterations must be positive integer and greater than 0")
  } else {
    iloreg.object@max.number.of.iterations <- max.number.of.iterations
  }

  if (!is.numeric(threads) | threads < 0 | threads%%1!=0)
  {
    stop("threads must be positive integer (0=detect automatically)")
  } else {
    iloreg.object@threads <- threads
  }

  set.seed(seed)

  parallelism <- TRUE

  if (threads==0) {
    cl <- makeCluster(detectCores(logical=TRUE)-1)
    registerDoParallel(cl)
  } else if(threads==1) {
    cat("Parallelism disabled, because threads=1\n")
    parallelism <- FALSE
  } else {
    cl<-makeCluster(threads)
    registerDoParallel(cl)
  }

  if (parallelism) {

    seeds <- sample(1:10^9,L,replace=FALSE)
    iloreg_out <- foreach(task = 1:L,.verbose = FALSE,
                          .combine = list,
                          .maxcombine = 1000,
                          .inorder = FALSE,
                          .export = c('RunICP','LogisticRegression'),
                          .packages=c("tictoc","Matrix","aricode","LiblineaR","SparseM"),
                          .multicombine = TRUE)  %dopar% {
                            set.seed(seeds[task])
                            RunICP(normalized.data = iloreg.object@normalized.data,
                                   k = iloreg.object@k,
                                   d = iloreg.object@d,
                                   r = iloreg.object@r,
                                   C = iloreg.object@C,
                                   type = type,
                                   max.number.of.iterations = max.number.of.iterations)
                          }
    stopCluster(cl)

  } else {
    iloreg_out <- list()
    for (l in 1:L) {
      res <- RunICP(normalized.data = iloreg.object@normalized.data,
                    k = iloreg.object@k,
                    d = iloreg.object@d,
                    r = iloreg.object@r,
                    C = iloreg.object@C,
                    type = type,
                    max.number.of.iterations = max.number.of.iterations)
      consensus_probability[[l]] <- res
    }
  }

  iloreg.object@consensus.probability <- lapply(iloreg_out,function(x) x$probabilities)
  iloreg.object@metrics <- lapply(iloreg_out,function(x) x$metrics)

  return(iloreg.object)
})



setGeneric("VisualizeQC", function(iloreg.object=NULL,return.plot=FALSE){
  standardGeneric("VisualizeQC")
})

#' @title Visualize the terminal projection accuracy, number of epochs/iterations and the average pairwise ARI values.
#'
#' @rdname VisualizeQC
#' @name VisualizeQC
#'
#' @description
#' Draw violin plots of the terminal projection accuracy, the number of epochs/iterations and the average pairwise ARI values.
#'
#'
#' @param iloreg.object object of class 'iloreg'
#' @param return.plot logical indicating if the ggplot2 object should be returned (default FALSE)
#' @return ggplot2 object if return.plot=TRUE
#' @keywords terminal projection accuracy average pairwise ARI quality control QC
#' @import ggplot2
#' @import cowplot
#' @importFrom aricode ARI
#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
setMethod("VisualizeQC", "iloreg", function(iloreg.object,return.plot){

  final_aris <- unlist(lapply(iloreg.object@metrics,function(x) x["ARI",ncol(x)]))

  df <- data.frame(matrix(final_aris,ncol = 1,dimnames = list(1:length(final_aris),c("value"))))
  df$Measure <- "CPA"
  df1 <- df

  number_of_runs <- unlist(lapply(iloreg.object@metrics,function(x) ncol(x)))

  df <- data.frame(matrix(number_of_runs,ncol = 1,dimnames = list(1:length(number_of_runs),c("value"))))
  df$Measure <- "Epochs"
  df2 <- df

  cluster_list <- lapply(iloreg.object@consensus.probability,function(x) apply(x,1,which.max))
  ARI_matrix <- matrix(NA,nrow = iloreg.object@L,ncol = iloreg.object@L)
  for (i in 1:iloreg.object@L)
  {
    for (j in 1:iloreg.object@L)
    {
      if (i < j)
      {
        next
      }
      ari <- ARI(cluster_list[[i]],cluster_list[[j]])
      ARI_matrix[i,j] <- ari
      ARI_matrix[j,i] <- ari

    }
  }

  pairwise_aris <- apply(ARI_matrix,1,mean)
  df <- data.frame(matrix(pairwise_aris,ncol = 1,dimnames = list(1:length(pairwise_aris),c("value"))))
  df$Measure <- "Average Pairwise ARI"
  df3 <- df

  # df <- rbind(df1,df2,df3)


  p1 <- ggplot(df1, aes(x=Measure,y=value))+
    geom_violin(trim=TRUE,fill="#F8766D") +
    theme_bw() +
    ylab("Clustering prediction accuracy") +
    xlab("") +
    # geom_boxplot(width=0.2)+
    # geom_jitter(shape=16, position=position_jitter(0.2)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ylim(0,1)


  p2 <- ggplot(df2, aes(x=Measure,y=value))+
    geom_violin(trim=TRUE,fill="#F8766D") +
    theme_bw() +
    ylab("Number of Epochs") +
    xlab("") +
    # geom_boxplot(width=0.2)+
    # geom_jitter(shape=16, position=position_jitter(0.2)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())



  p3 <- ggplot(df3, aes(x=Measure,y=value))+
    geom_violin(trim=TRUE,fill="#F8766D") +
    theme_bw() +
    ylab("Average pairwise ARI") +
    xlab("") +
    # geom_boxplot(width=0.2)+
    # geom_jitter(shape=16, position=position_jitter(0.2)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ylim(0,1)

  p <- plot_grid(p1,p2,p3,nrow=1)

  if (return.plot) {
    return(p)
  } else {
    print(p)
  }
})



setGeneric("RunPCA", function(iloreg.object=NULL,p=50,scale=FALSE){
  standardGeneric("RunPCA")
})

#' @title PCA transformation of the joint probability matrix
#'
#' @rdname RunPCA
#' @name RunPCA
#'
#' @description
#' Perform the PCA transformation of the joint probability matrix, which reduces the dimensionality from k*L to p
#'
#' @param iloreg.object object of class 'iloreg'
#' @param p a positive integer denoting the number of principal components to calculate and select (default 50)
#' @param scale a logical specifying whether the probabilities should be standardized to unit-variance before running PCA (default FALSE)
#' @return iloreg object
#' @keywords PCA eigendecomposition
#' @importFrom RSpectra eigs_sym
#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
setMethod("RunPCA", "iloreg", function(iloreg.object,p,scale){

  iloreg.object@number.of.pcs <- p
  iloreg.object@scale.pca <- scale

  X <- do.call(cbind,iloreg.object@consensus.probability)

  X <- scale(X,scale = scale,center = TRUE)

  # X^T %*% X
  A = crossprod(X)

  # Perform eigendecomposition
  eigs_sym_out <- eigs_sym(A, p, which = "LM")

  rotated <- X %*% eigs_sym_out$vectors
  colnames(rotated) <- paste0("PC",1:ncol(rotated))

  iloreg.object@rotated.consensus.probability <- rotated

  return(iloreg.object)

})


setGeneric("PCAElbowPlot", function(iloreg.object=NULL,return.plot=FALSE){
  standardGeneric("PCAElbowPlot")
})

#' @title Elbow plot of the standard deviations of the principal components
#'
#' @rdname PCAElbowPlot
#' @name PCAElbowPlot
#'
#' @description
#' Draw an elbow plot of the standard deviations of the principal components to deduce an appropriate value for p.
#'
#'
#' @param iloreg.object object of class 'iloreg'
#' @param return.plot logical indicating if the ggplot2 object should be returned (default FALSE)
#' @return ggplot2 object if return.plot=TRUE
#' @keywords PCA elbow plot
#' @import ggplot2
#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)

setMethod("PCAElbowPlot", "iloreg", function(iloreg.object,return.plot){


  df <- matrix(apply(iloreg.object@rotated.consensus.probability,2,sd),nrow = ncol(iloreg.object@rotated.consensus.probability),ncol = 1,dimnames = list(1:ncol(iloreg.object@rotated.consensus.probability),"SD"))
  df <- reshape2::melt(df)

  p<-ggplot(df, aes(x=Var1, y=value)) +
    geom_line(color="blue")+
    geom_point(color="black") +
    theme_bw()+
    ylab("Standard Deviation")+
    xlab("PC")


  if (return.plot) {
    return(p)
  } else {
    print(p)
  }
})



setGeneric("RunUMAP", function(iloreg.object=NULL){
  standardGeneric("RunUMAP")
})

#' @title Uniform Manifold Approximation and Projection (UMAP)
#'
#' @rdname RunUMAP
#' @name RunUMAP
#'
#' @description
#' Run nonlinear dimensionality reduction using UMAP with the PCA-transformed consensus matrix as input.
#'
#'
#' @param iloreg.object object of class 'iloreg'
#' @return iloreg object
#' @keywords Uniform Manifold Approximation and Projection UMAP
#' @importFrom umap umap
#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- RunUMAP(iloreg_object)
setMethod("RunUMAP", "iloreg", function(iloreg.object){

  umap_out <- umap(iloreg.object@rotated.consensus.probability)

  iloreg.object@umap.embeddings <- umap_out$layout

  return(iloreg.object)

})



setGeneric("RunTSNE", function(iloreg.object=NULL,perplexity=30){
  standardGeneric("RunTSNE")
})

#' @title Barnes-Hut implementation of t-Distributed Stochastic Neighbor Embedding (t-SNE)
#'
#' @rdname RunTSNE
#' @name RunTSNE
#'
#' @description
#' Run nonlinear dimensionality reduction using t-SNE with the PCA-transformed consensus matrix as input.
#'
#'
#' @param iloreg.object object of class 'iloreg'
#' @param perplexity perplexity of t-SNE
#' @return iloreg object
#' @keywords  Barnes-Hut implementation of t-Distributed Stochastic Neighbor Embedding t-SNE
#' @importFrom Rtsne Rtsne
#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- RunTSNE(iloreg_object)
setMethod("RunTSNE", "iloreg", function(iloreg.object,perplexity){

  Rtsne_out <- Rtsne(iloreg.object@rotated.consensus.probability,
                     is_distance=FALSE,
                     perplexity=perplexity,
                     pca=FALSE)

  iloreg.object@tsne.embeddings <- Rtsne_out$Y

  return(iloreg.object)

})




setGeneric("HierarchicalClustering", function(iloreg.object=NULL){
  standardGeneric("HierarchicalClustering")
})

#' @title Hierarchical clustering using the Ward's method
#'
#' @rdname HierarchicalClustering
#' @name HierarchicalClustering
#'
#' @description
#' Perform Hierarchical clustering using the Ward's method.
#'
#' @param iloreg.object object of class 'iloreg'
#' @return iloreg object
#' @keywords ward hierarchical clustering
#' @importFrom fastcluster hclust.vector
#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- HierarchicalClustering(iloreg_object)
setMethod("HierarchicalClustering", "iloreg", function(iloreg.object){


  hc <- hclust.vector(iloreg.object@rotated.consensus.probability,
                      method = "ward")

  iloreg.object@hc <- hc

  return(iloreg.object)

})






setGeneric("CalculateSilhouetteInformation", function(iloreg.object=NULL,K.range=2:50){
  standardGeneric("CalculateSilhouetteInformation")
})

#' @title Silhouette method for estimating K
#'
#' @rdname CalculateSilhouetteInformation
#' @name CalculateSilhouetteInformation
#'
#' @description
#' Estimate the optimal number of clusters K from the dendrogram of the hierarhical clustering using the silhouette method.
#'
#'
#' @param iloreg.object object of class 'iloreg'
#' @param K.range a numeric vector for the different K values to be tested (default 2:50)
#' @return iloreg Object
#' @keywords silhouette
#' @importFrom parallelDist parDist
#' @importFrom cluster silhouette
#' @importFrom dendextend cutree

#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- HierarchicalClustering(iloreg_object)
#' iloreg_object <- CalculateSilhouetteInformation(iloreg_object,K.range = 2:50)
setMethod("CalculateSilhouetteInformation", "iloreg", function(iloreg.object,K.range){

  distance_matrix <- parDist(iloreg.object@rotated.consensus.probability, method = "euclidean", threads = 1)
  distance_matrix <- as.matrix(distance_matrix)
  sis <- c()
  for (k in k.range)
  {
    clustering <- cutree(iloreg.object@hc,k=k)

    si <- silhouette(clustering,dmatrix = distance_matrix)
    avgsi <- summary(si)$avg.width
    sis <- c(sis,avgsi)
  }
  # Select optimal k and cluster the data
  k_optimal <- which.max(sis)+1

  cat(paste0("optimal k: ",k_optimal, ", average silhouette score: ",sis[which.max(sis)],"\n"))

  clustering <- factor(cutree(as.dendrogram(iloreg.object@hc),k=k_optimal))
  names(clustering) <- colnames(iloreg.object@normalized.data)

  iloreg.object@clustering.optimal <- clustering
  iloreg.object@K.optimal <- k_optimal

  names(sis) <- k.range
  iloreg.object@silhouette.information <- sis

  return(iloreg.object)

})




setGeneric("SilhouetteCurve", function(iloreg.object=NULL,return.plot=FALSE){
  standardGeneric("SilhouetteCurve")
})

#' @title Silhouette curve
#'
#' @rdname SilhouetteCurve
#' @name SilhouetteCurve
#'
#' @description
#' Draw the silhouette curve: the average silhouette value across the cells for different values of K.
#'
#'
#' @param iloreg.object object of class 'iloreg'
#' @param return.plot logical indicating if the ggplot2 object should be returned (default FALSE)
#' @return ggplot2 object if return.plot=TRUE
#' @keywords silhouette curve
#' @import ggplot2
#' @importFrom DescTools AUC

#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- HierarchicalClustering(iloreg_object)
#' iloreg_object <- CalculateSilhouetteInformation(iloreg_object,K.range = 2:50)
#' SilhouetteCurve(iloreg_object)
setMethod("SilhouetteCurve", "iloreg", function(iloreg.object,return.plot){

  df <- data.frame(cbind(names(iloreg.object@silhouette.information),iloreg.object@silhouette.information),
                   stringsAsFactors = FALSE)
  colnames(df) <- c("K","AvgSilhouette")
  df$AvgSilhouette <- as.numeric(df$AvgSilhouette)
  df$K <- as.numeric(df$K)

  auc <- round(AUC(df$K,df$AvgSilhouette),3)

  p<-ggplot(df, aes(x=K, y=AvgSilhouette)) +
    geom_line(color="red")+
    geom_point(color="black")+
    ylab("Average silhouette")+
    theme_bw() +
    ggtitle(paste0("AUSC=",auc))


  if (return.plot)
  {
    return(p)
  } else {
    print(p)
  }




})




setGeneric("SelectKClusters", function(iloreg.object=NULL,K=10){
  standardGeneric("SelectKClusters")
})

#' @title Selecting K clusters from hierarchical clustering
#'
#' @rdname SelectKClusters
#' @name SelectKClusters
#'
#' @description
#' Selects K clusters from the dendrogram using the cutree function from the dendextend R package.
#'
#'
#' @param iloreg.object object of class 'iloreg'
#' @param K a positive integer denoting how many clusters to select
#' @return iloreg Object
#' @keywords iterative logistic regression ILoReg consensus clustering
#' @importFrom dendextend cutree

#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- HierarchicalClustering(iloreg_object)
#' iloreg_object <- SelectKClusters(iloreg_object,K=25)
setMethod("SelectKClusters", "iloreg", function(iloreg.object,K){

  clustering <- factor(cutree(as.dendrogram(iloreg.object@hc),k=K))
  names(clustering) <- colnames(iloreg.object@normalized.data)

  iloreg.object@clustering.manual <- clustering
  iloreg.object@K.manual <- K

  return(iloreg.object)

})


setGeneric("MergeClusters", function(iloreg.object=NULL,clusters.to.merge="",new.name=""){
  standardGeneric("MergeClusters")
})

#' @title Merge clusters
#'
#' @rdname MergeClusters
#' @name MergeClusters
#'
#' @description
#' MergeClusters function enables merging clusters into one.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param clusters.to.merge a character or numeric vector for the names of the clusters to merge
#' @param new.name a character for the new name of the merged cluster. If empty, the new cluster name is formed by separating the cluster names by "_".
#' @return iloreg Object
#' @keywords iterative logistic regression ILoReg consensus clustering

#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- HierarchicalClustering(iloreg_object)
#' iloreg_object <- SelectKClusters(iloreg_object,K=25)
#' iloreg_object <- MergeClusters(iloreg_object,clusters.to.merge = c(23,6),"NewCluster")
setMethod("MergeClusters", "iloreg", function(iloreg.object,clusters.to.merge,new.name){

  clusters.to.merge <- as.character(clusters.to.merge)

  clustering_old <- iloreg.object@clustering.manual
  clusters_old <- levels(clustering_old)

  if (sum(clusters.to.merge %in% clusters_old)!=length(clusters.to.merge))
  {
    stop("invalid clusters.to.merge argument")
  }

  if (new.name=="")
  {
    new_cluster_name <- paste(clusters.to.merge,collapse = ",")
  } else {
    new_cluster_name <- new.name
  }

  clustering_new <- as.character(clustering_old)
  clustering_new[clustering_new %in% clusters.to.merge] <- new_cluster_name

  clustering_new <- factor(clustering_new)
  names(clustering_new) <- names(clustering_old)

  iloreg.object@clustering.manual <- clustering_new
  iloreg.object@K.manual <- length(levels(clustering_new))

  return(iloreg.object)

})



setGeneric("RenameAllClusters", function(iloreg.object=NULL,new.cluster.names=""){
  standardGeneric("RenameAllClusters")
})

#' @title Renaming all clusters at once
#'
#' @rdname RenameAllClusters
#' @name RenameAllClusters
#'
#' @description
#' RenameAllClusters function enables renaming all cluster at once.
#'
#'
#' @param iloreg.object object of class 'iloreg'
#' @param new.cluster.names object of class 'iloreg'
#' @return iloreg object
#' @keywords rename all clusters
#' @importFrom plyr mapvalues

#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- HierarchicalClustering(iloreg_object)
#' iloreg_object <- SelectKClusters(iloreg_object,K=15)
#' ## Rename all 15 clusters to letters
#' iloreg_object <- RenameAllClusters(iloreg_object,new.cluster.names = LETTERS[1:iloreg_object@K.manual])

setMethod("RenameAllClusters", "iloreg", function(iloreg.object,new.cluster.names){

  new.cluster.names <- as.character(new.cluster.names)

  clustering_old <- iloreg.object@clustering.manual
  clusters_old <- levels(clustering_old)

  if (length(clusters_old) != length(new.cluster.names))
  {
    stop("number of elements in clusters.to.merge is not equal to the current number of clusters")
  }


  clustering_new <- mapvalues(clustering_old,clusters_old,new.cluster.names)

  iloreg.object@clustering.manual <- clustering_new

  return(iloreg.object)

})




setGeneric("RenameCluster", function(iloreg.object=NULL,old.cluster.name="",new.cluster.name=""){
  standardGeneric("RenameCluster")
})

#' @title Renaming one cluster
#'
#' @rdname RenameCluster
#' @name RenameCluster
#'
#' @description
#' RenameCluster function enables renaming a cluster
#'
#' @param iloreg.object object of class 'iloreg'
#' @param old.cluster.name object of class 'iloreg'
#' @param new.cluster.name object of class 'iloreg'
#' @return iloreg object
#' @keywords rename one cluster

#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- HierarchicalClustering(iloreg_object)
#' iloreg_object <- SelectKClusters(iloreg_object,K=15)
#' ## Rename cluster 1 to "A"
#' iloreg_object <- RenameCluster(iloreg_object,old.cluster.name=1,new.cluster.name="A")
setMethod("RenameCluster", "iloreg", function(iloreg.object,old.cluster.name,new.cluster.name){


  old.cluster.name <- as.character(old.cluster.name)
  new.cluster.name <- as.character(new.cluster.name)

  if (old.cluster.name=="" | new.cluster.name=="")
  {
    stop("'old.cluster.name' or 'new.cluster.name' empty\n")
  }

  clustering_old <- iloreg.object@clustering.manual
  clusters_old <- levels(clustering_old)
  clustering_old <- as.character(clustering_old)

  if (!(old.cluster.name %in% clusters_old))
  {
    stop("'old.cluster.name' unvalid cluster name\n")
  }

  clustering_new <- clustering_old
  clustering_new[clustering_new==old.cluster.name] <- new.cluster.name

  clustering_new <- factor(clustering_new)
  names(clustering_new) <- names(clustering_old)

  iloreg.object@clustering.manual <- clustering_new

  return(iloreg.object)

})



setGeneric("GeneScatterPlot", function(iloreg.object=NULL,genes="",return.plot=FALSE,dim.reduction.type="tsne",point.size=0.7,title=""){
  standardGeneric("GeneScatterPlot")
})

#' @title Visualize gene expression over nonlinear dimensionality reduction
#'
#' @rdname GeneScatterPlot
#' @name GeneScatterPlot
#'
#' @description
#' GeneScatterPlot enables visualizing gene expression of a gene over nonlinear dimensionality reduction with t-SNE or UMAP.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param genes a character vector of the genes to be visualized
#' @param return.plot whether to return the ggplot2 object or just draw it (default \code{FALSE})
#' @param dim.reduction.type "tsne" or "umap" (default "tsne")
#' @param point.size point size (default 0.7)
#' @param title text to write above the plot
#' @return ggplot2 object if return.plot=TRUE
#' @keywords gene scatter plot visualization
#' @import ggplot2
#' @importFrom scales muted
#' @importFrom cowplot plot_grid

#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- RunTSNE(iloreg_object,perplexity=30)
#' GeneScatterPlot(iloreg_object,c("CCR7","CD3D","S100A4"),dim.reduction.type = "tsne",point.size=0.7)
#' iloreg_object <- RunUMAP(iloreg_object,perplexity=30)
#' GeneScatterPlot(iloreg_object,c("CCR7","CD3D","S100A4"),dim.reduction.type = "umap",point.size=0.7)

setMethod("GeneScatterPlot", "iloreg", function(iloreg.object,genes,return.plot,dim.reduction.type,point.size,title){

  if (dim.reduction.type=="umap")
  {
    two.dim.data <- iloreg.object@umap.embeddings
    xlab <- "UMAP_1"
    ylab <- "UMAP_2"
  } else if (dim.reduction.type=="tsne"){
    two.dim.data <- iloreg.object@tsne.embeddings
    xlab <- "tSNE_1"
    ylab <- "tSNE_2"
  } else {
    stop("dim.reduction.type must be either 'tsne' or 'umap'")
  }

  if (length(genes)==1)
  {

    df <- as.data.frame(two.dim.data)

    if (!(genes %in% rownames(iloreg.object@normalized.data)))
    {
      stop("invalid gene name")
    }

    color.by <- iloreg.object@normalized.data[genes,]
    df$group <- color.by
    colnames(df) <- c("dim1","dim2","group")

    if (title!="")
    {
      p<-ggplot(df, aes(x=dim1, y=dim2)) +
        geom_point(size=point.size,aes(color=group)) +
        scale_colour_gradient2(low = muted("red"), mid = "lightgrey",
                               high = "blue",name = genes) +
        xlab(xlab) +
        ylab(ylab) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))

    } else {
      p<-ggplot(df, aes(x=dim1, y=dim2)) +
        geom_point(size=point.size,aes(color=group)) +
        scale_colour_gradient2(low = muted("red"), mid = "lightgrey",
                               high = "blue",name = genes) +
        xlab(xlab) +
        ylab(ylab) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5))

    }
    if (return.plot) {
      return(p)
    } else {
      print(p)
    }

  } else {

    plot_list <- list()
    for (gene in genes)
    {
      df <- as.data.frame(two.dim.data)

      if (!(gene %in% rownames(iloreg.object@normalized.data)))
      {
        stop(paste0("invalid gene name: ",gene))
      }

      color.by <- iloreg.object@normalized.data[gene,]
      df$group <- color.by
      colnames(df) <- c("dim1","dim2","group")

      if (title!="") {
        p<-ggplot(df, aes(x=dim1, y=dim2)) +
          geom_point(size=point.size,aes(color=group)) +
          scale_colour_gradient2(low = muted("red"), mid = "lightgrey",
                                 high = "blue",name = gene) +
          xlab(xlab) +
          ylab(ylab) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
      } else {
        p<-ggplot(df, aes(x=dim1, y=dim2)) +
          geom_point(size=point.size,aes(color=group)) +
          scale_colour_gradient2(low = muted("red"), mid = "lightgrey",
                                 high = "blue",name = gene) +
          xlab(xlab) +
          ylab(ylab) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          ggtitle(title) +
          theme(plot.title = element_text(hjust = 0.5))

      }

      plot_list[[gene]] <- p

    }

    p <- plot_grid(plotlist = plot_list,align = "hv")

    if (return.plot) {
      return(p)
    } else {
      print(p)
    }

  }
})


setGeneric("ClusteringScatterPlot", function(iloreg.object=NULL,clustering.type="manual",return.plot=FALSE,dim.reduction.type="",point.size=0.7,title=""){
  standardGeneric("ClusteringScatterPlot")
})

#' @title Visualize the clustering over nonliner dimensionality reduction
#'
#' @rdname ClusteringScatterPlot
#' @name ClusteringScatterPlot
#'
#' @description
#' ClusteringScatterPlot function enables visualizing the clustering over nonliner dimensionality reduction (t-SNE or UMAP)
#'
#' @param iloreg.object object of class 'iloreg'
#' @param clustering.type "manual" or "optimal". "manual" refers to the clustering formed using the "SelectKClusters" function and "optimal" to the clustering formed using the "CalculateSilhouetteInformation" function. (default "manual")
#' @param return.plot object of class 'iloreg' (default \code{FALSE})
#' @param dim.reduction.type "tsne" or "umap" (default "tsne")
#' @param point.size point size (default 0.7)
#' @param title text to write above the plot
#' @return ggplot2 object if return.plot=TRUE
#' @keywords clustering scatter plot nonlinear dimensionality reduction
#' @import ggplot2

#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- RunTSNE(iloreg_object,perplexity=30)
#' iloreg_object <- HierarchicalClustering(iloreg_object)
#' iloreg_object <- SelectKClusters(iloreg_object,K=15)
#' ClusteringScatterPlot(iloreg_object,clustering.type="manual",dim.reduction.type = "tsne",point.size=0.7)
#' iloreg_object <- CalculateSilhouetteInformation(iloreg_object,k.range = 2:50)
#' ClusteringScatterPlot(iloreg_object,clustering.type="optimal",dim.reduction.type = "tsne",point.size=0.7)
#' iloreg_object <- RunUMAP(iloreg_object,perplexity=30)
#' ClusteringScatterPlot(iloreg_object,clustering.type="manual",dim.reduction.type = "umap",point.size=0.7)
#' ClusteringScatterPlot(iloreg_object,clustering.type="optimal",dim.reduction.type = "umap",point.size=0.7)
setMethod("ClusteringScatterPlot", "iloreg", function(iloreg.object,clustering.type,return.plot,dim.reduction.type,point.size,title){

  if (dim.reduction.type=="umap")
  {
    two.dim.data <- iloreg.object@umap.embeddings
    xlab <- "UMAP_1"
    ylab <- "UMAP_2"
  } else if (dim.reduction.type=="tsne"){
    two.dim.data <- iloreg.object@tsne.embeddings
    xlab <- "tSNE_1"
    ylab <- "tSNE_2"
  } else {
    stop("dim.reduction.type must be either 'tsne' or 'umap'")
  }

  if (clustering.type=="manual")
  {
    color.by <- iloreg.object@clustering.manual
  } else if (clustering.type=="optimal")
  {
    color.by <- iloreg.object@clustering.optimal
  } else {
    clustering <- iloreg.object@clustering.manual
    cat("clustering.type='manual'")
  }

  df <- as.data.frame(two.dim.data)

  df$cluster <- color.by
  colnames(df) <- c("dim1","dim2","cluster")

  two.dim.data_ <- two.dim.data
  rownames(two.dim.data_) <- names(color.by)
  cluster_centers <- lapply(levels(color.by),function(x) apply(two.dim.data_[names(color.by)[color.by==x],,drop=FALSE],2,median))
  cluster_centers <- do.call(rbind,cluster_centers)

  if (title=="")
  {
    p<-ggplot(df, aes(x=dim1, y=dim2)) +
      geom_point(size=point.size,aes(color=cluster)) +
      xlab(xlab) +
      ylab(ylab) +
      theme_classic() +
      annotate("text", x = cluster_centers[,1], y = cluster_centers[,2], label = levels(color.by))

  } else {

    p<-ggplot(df, aes(x=dim1, y=dim2)) +
      geom_point(size=point.size,aes(color=cluster)) +
      xlab(xlab) +
      ylab(ylab) +
      theme_classic() +
      annotate("text", x = cluster_centers[,1], y = cluster_centers[,2], label = levels(color.by)) +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))

  }

  if (return.plot) {
    return(p)
  } else {
    plot(p)
  }

})

setGeneric("FindAllGeneMarkers", function(iloreg.object=NULL,clustering.type="",test="wilcox",
                                          logfc.threshold = 0.25,
                                          min.pct = 0.1,
                                          min.diff.pct = NULL,
                                          min.cells.group = 3,
                                          max.cells.per.cluster = NULL,
                                          random.seed = 1,
                                          pseudocount.use = 1,
                                          return.thresh = 0.01,
                                          only.pos=FALSE){
  standardGeneric("FindAllGeneMarkers")
})

#' @title identification of gene markers for all clusters
#'
#' @rdname FindAllGeneMarkers
#' @name FindAllGeneMarkers
#'
#' @description
#' FindAllGeneMarkers enables identifying gene markers for all clusters at once. This is done by differential expresission analysis
#' where cells from one cluster are compared against the cells from the rest of the clusters. Gene and cell filters can be applied to accelerate
#' the analysis, but this might lead to missing weak signals.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param clustering.type "manual" or "optimal". "manual" refers to the clustering formed using the "SelectKClusters" function and "optimal" to the clustering formed using the "CalculateSilhouetteInformation" function. (default "manual")
#' @param test Which test to use. Only "wilcoxon" (the Wilcoxon rank-sum test, AKA Mann-Whitney U test) is supported at the moment.
#' @param logfc.threshold Filters out genes that have log2 fold-change of the averaged gene expression values (with the pseudo-count value added to the averaged values before division if pseudocount.use > 0) below this threshold. (default 0.25)
#' @param min.pct Filters out genes that have dropout rate (fraction of cells expressing a gene) below this threshold in both comparison groups (default 0.1).
#' @param min.diff.pct Filters out genes that do not have this minimum difference in the dropout rates (fraction of cells expressing a gene) between the two comparison groups  (default NULL).
#' @param min.cells.group The minimum number of cells in the two comparison groups to perform the DE analysis. If the number of cells is below the threshold (default 3), then the DE analysis of this cluster is skipped.
#' @param max.cells.per.cluster The maximun number of cells per cluster if downsampling is performed to speed up the DE analysis. (default NULL, i.e. no downsampling).
#' @param random.seed random seed for for downsampling if max.cells.per.cluster is not NULL (default 1)
#' @param pseudocount.use A positive integer (default 1), which is added to the average gene expression values before calculating the fold-change. This makes sure that no divisions by zero occur.
#' @param return.thresh If only.pos=TRUE, then return only genes that have the adjusted p-value (adjusted by the Bonferroni method) below or equal to this threshold (default 0.01).
#' @param only.pos Whether to return only genes that have an adjusted p-value (adjusted by the Bonferroni method) below or equal to the threshold (default FALSE).
#' @return a data frame of the results if positive results were found, else NULL
#' @keywords differential expression DE analysis gene markers

#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- RunTSNE(iloreg_object,perplexity=30)
#' iloreg_object <- HierarchicalClustering(iloreg_object)
#' iloreg_object <- SelectKClusters(iloreg_object,K=15)
#' gene_markers <- FindAllGeneMarkers(iloreg_object,
#'                                   clustering.type = "manual",
#'                                   test = "wilcox",
#'                                   logfc.threshold = 0.25,
#'                                   min.pct = 0.1,
#'                                  min.diff.pct = NULL,
#'                                  pseudocount.use = 1,
#'                                  min.cells.group = 3,
#'                                 return.thresh = 0.01,
#'                                 only.pos = FALSE,max.cells.per.cluster = NULL)
#' library(dplyr)
#' ## Select top 10 genes per cluster by log2 fold-change
#' gene_markers %>% group_by(cluster) %>% top_n(10, log2FC) -> top10
#' gene_markers %>% group_by(cluster) %>% top_n(1, log2FC) -> top1


setMethod("FindAllGeneMarkers", "iloreg", function(iloreg.object,
                                                   clustering.type,
                                                   test,
                                                   logfc.threshold,
                                                   min.pct,
                                                   min.diff.pct,
                                                   min.cells.group,
                                                   max.cells.per.cluster,
                                                   random.seed,
                                                   pseudocount.use,
                                                   return.thresh,
                                                   only.pos)
{

  number.of.expressed.genes <- nrow(iloreg.object@normalized.data)

  if (clustering.type=="manual")
  {
    clustering <- iloreg.object@clustering.manual
  } else if (clustering.type=="optimal")
  {
    clustering <- iloreg.object@clustering.optimal
  } else {
    clustering <- iloreg.object@clustering.manual
    cat("clustering.type='manual'")
  }

  data <- iloreg.object@normalized.data

  clusters <- levels(clustering)

  # Downsample each cluster
  if (!is.null(max.cells.per.cluster))
  {
    cells_downsampled <- c()
    for (cluster in clusters)
    {
      cells_in_cluster <- table(clustering)[cluster]
      if (max.cells.per.cluster < cells_in_cluster)
      {
        set.seed(random.seed)
        inds <- sample(1:cells_in_cluster,size = max.cells.per.cluster,replace = FALSE)
        cells_downsampled <- c(cells_downsampled,names(clustering[clustering==cluster])[inds])
      }
    }
    data <- data[,cells_downsampled]
    clustering <- clustering[cells_downsampled]
  }

  # Compare cells from each cluster against all other clusters
  results_list <- list()

  for (cluster in clusters)
  {
    cat("-----------------------------------\n")
    cat(paste0("testing cluster ",cluster,"\n"))
    # Extract data

    data_cluster <- data[,clustering==cluster]
    data_other <- data[,clustering!=cluster]

    # Skip if the number of cells in the test or the reference set is lower than min.cells.group
    if (ncol(data_cluster) < min.cells.group | ncol(data_other) < min.cells.group)
    {
      cat("-----------------------------------\n")
      next
    }

    # min.pct filter
    genes.pct_cluster <- apply(data_cluster,1,function(x) sum(x!=0))/ncol(data_cluster)
    genes.pct_other <- apply(data_other,1,function(x) sum(x!=0))/ncol(data_other)

    genes_to_include <- rownames(data_cluster)[genes.pct_cluster>=min.pct | genes.pct_other >= min.pct]

    data_cluster <- data_cluster[genes_to_include,,drop=FALSE]
    data_other <- data_other[genes_to_include,,drop=FALSE]

    cat(paste0(nrow(data_cluster)," genes left after min.pct filtering\n"))
    if (nrow(data_cluster)==0)
    {
      cat("-----------------------------------\n")
      next
    }

    # min.diff.pct filter
    if (!is.null(min.diff.pct))
    {
      genes.pct_cluster <- genes.pct_cluster[genes_to_include]
      genes.pct_other <- genes.pct_other[genes_to_include]

      genes_to_include <- rownames(data_cluster)[abs(genes.pct_cluster-genes.pct_other) >= min.diff.pct]

      data_cluster <- data_cluster[genes_to_include,,drop=FALSE]
      data_other <- data_other[genes_to_include,,drop=FALSE]

    }

    cat(paste0(nrow(data_cluster)," genes left after min.diff.pct filtering\n"))
    if (nrow(data_cluster)==0)
    {
      cat("-----------------------------------\n")
      next
    }

    # logfc.threshold filter
    # Calculate log2 fold changes
    cluster_aves <- apply(data_cluster,1,mean)
    other_aves <- apply(data_other,1,mean)

    log2FC <- log2((cluster_aves+pseudocount.use)/(other_aves+pseudocount.use))

    genes_to_include <- rownames(data_cluster)[log2FC >= logfc.threshold | log2FC <= -logfc.threshold]

    data_cluster <- data_cluster[genes_to_include,,drop=FALSE]
    data_other <- data_other[genes_to_include,,drop=FALSE]


    cat(paste0(nrow(data_cluster)," genes left after logfc.threshold filtering\n"))
    if (nrow(data_cluster)==0)
    {
      cat("-----------------------------------\n")
      next
    }

    # Run DE test

    if (test=="wilcox")
    {
      wilcox.res <- lapply(rownames(data_cluster),function(x) wilcox.test(x=data_cluster[x,],y=data_other[x,]))
      p_values <- unlist(lapply(wilcox.res,function(x) x$p.value))
      names(p_values) <- rownames(data_cluster)

      # Adjust p-values
      adj_p_values <- p.adjust(p_values, method = "bonferroni", n = number.of.expressed.genes)

      res <- cbind(p_values,adj_p_values,log2FC[names(p_values)],genes.pct_cluster[names(p_values)],genes.pct_other[names(p_values)],abs(genes.pct_cluster-genes.pct_other)[names(p_values)])
      colnames(res) <- c("p.value","adj.p.value","log2FC","pct.1","pct.2","diff.pct")
      res <- as.data.frame(res)
      res$cluster <- cluster
      res$gene <- names(p_values)
    }

    results_list[[cluster]] <- res

    cat("-----------------------------------\n")

  }

  results_df <- do.call(rbind,results_list)
  rownames(results_df) <- make.unique(unlist(lapply(results_list,rownames)))

  if(only.pos) {
    results_df <- results_df[results_df$adj.p.value <= return.thresh,]
    return(results_df)
  }
  return(results_df)
})





setGeneric("FindGeneMarkers", function(iloreg.object=NULL,
                                       clusters.1 = NULL,
                                       clusters.2 = NULL,
                                       clustering.type="",
                                       test="wilcox",
                                       logfc.threshold = 0.25,
                                       min.pct = 0.1,
                                       min.diff.pct = NULL,
                                       min.cells.group = 3,
                                       max.cells.per.cluster = NULL,
                                       random.seed = 1,
                                       pseudocount.use = 1,
                                       return.thresh = 0.01,
                                       only.pos=FALSE){
  standardGeneric("FindGeneMarkers")
})

#' @title Identification of gene markers for a cluster or two arbitrary combinations of clusters
#'
#' @rdname FindGeneMarkers
#' @name FindGeneMarkers
#'
#' @description
#' FindGeneMarkers enables identifying gene markers for one cluster or two arbitrary combinations of clusters, e.g. 1_2 vs. 3_4_5. Gene and cell filters can be applied to accelerate
#' the analysis, but this might lead to missing weak signals.
#'
#'
#' @param iloreg.object object of class 'iloreg'
#' @param clusters.1 a character or numeric vector denoting which clusters to use in the first group (named group.1 in the results)
#' @param clusters.2 a character or numeric vector denoting which clusters to use in the second group (named group.2 in the results)
#' @param clustering.type "manual" or "optimal". "manual" refers to the clustering formed using the "SelectKClusters" function and "optimal" to the clustering formed using the "CalculateSilhouetteInformation" function. (default "manual")
#' @param test Which test to use. Only "wilcoxon" (the Wilcoxon rank-sum test, AKA Mann-Whitney U test) is supported at the moment.
#' @param logfc.threshold Filters out genes that have log2 fold-change of the averaged gene expression values (with the pseudo-count value added to the averaged values before division if pseudocount.use > 0) below this threshold.
#' @param min.pct Filters out genes that have dropout rate (fraction of cells expressing a gene) below this threshold in both comparison groups (default 0.1).
#' @param min.diff.pct Filters out genes that do not have this minimum difference in the dropout rates (fraction of cells expressing a gene) between the two comparison groups  (default NULL).
#' @param min.cells.group The minimum number of cells in the two comparison groups to perform the DE analysis. If the number of cells is below the threshold (default 3), then the DE analysis of this cluster is skipped.
#' @param max.cells.per.cluster The minimum number of cells in the two comparison groups to perform the DE analysis. If the number of cells is below the threshold (default 3), then the DE analysis of this cluster is skipped.
#' @param random.seed random seed for for downsampling if max.cells.per.cluster is not NULL (default 1)
#' @param pseudocount.use A positive integer (default 1), which is added to the average gene expression values before calculating the fold-change. This makes sure that no divisions by zero occur.
#' @param return.thresh If only.pos=TRUE, then return only genes that have the adjusted p-value (adjusted by the Bonferroni method) below or equal to this threshold (default 0.01).
#' @param only.pos Whether to return only genes that have an adjusted p-value (adjusted by the Bonferroni method) below or equal to the threshold (default FALSE).
#' @return a data frame of the results if positive results were found, else NULL
#' @keywords differential expression DE analysis gene markers

#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- RunTSNE(iloreg_object,perplexity=30)
#' iloreg_object <- HierarchicalClustering(iloreg_object)
#' iloreg_object <- SelectKClusters(iloreg_object,K=15)
#' gene_markers_1_vs_2 <- FindGeneMarkers(iloreg_object,clusters.1 = c(1),clusters.2 = c(2),
#'                                   clustering.type = "manual",
#'                                   test = "wilcox",
#'                                   logfc.threshold = 0.25,
#'                                   min.pct = 0.1,
#'                                  min.diff.pct = NULL,
#'                                  pseudocount.use = 1,
#'                                  min.cells.group = 3,
#'                                 return.thresh = 0.01,
#'                                 only.pos = FALSE,max.cells.per.cluster = NULL)
#' library(dplyr)
#' ## Select top 10 genes per cluster by log2 fold-change
#' gene_markers  %>% top_n(10, log2FC) -> top10
#' gene_markers  %>% top_n(1, log2FC) -> top1
setMethod("FindGeneMarkers", "iloreg", function(iloreg.object,
                                                clusters.1,
                                                clusters.2,
                                                clustering.type,
                                                test,
                                                logfc.threshold,
                                                min.pct,
                                                min.diff.pct,
                                                min.cells.group,
                                                max.cells.per.cluster,
                                                random.seed,
                                                pseudocount.use,
                                                return.thresh,
                                                only.pos)
{

  if (clustering.type=="manual")
  {
    clustering <- iloreg.object@clustering.manual
  } else if (clustering.type=="optimal")
  {
    clustering <- iloreg.object@clustering.optimal
  } else {
    clustering <- iloreg.object@clustering.manual
    cat("clustering.type='manual'")
  }

  data <- iloreg.object@normalized.data

  cells_to_include_1 <- names(clustering)[clustering %in% clusters.1]
  clustering_1 <- factor(rep("group.1",length(cells_to_include_1)))
  names(clustering_1) <- cells_to_include_1

  cells_to_include_2 <- names(clustering)[clustering %in% clusters.2]
  clustering_2 <- factor(rep("group.2",length(cells_to_include_2)))
  names(clustering_2) <- cells_to_include_2

  data <- data[,c(cells_to_include_1,cells_to_include_2)]

  clustering <- factor(c(as.character(clustering_1),as.character(clustering_2)))
  names(clustering) <- c(cells_to_include_1,cells_to_include_2)

  clusters <- levels(clustering)

  # Remove genes that are not expressed in any of the cells
  data <- data[Matrix::rowSums(data)!=0,]

  clusters <- levels(clustering)

  # Downsample each cluster
  if (!is.null(max.cells.per.cluster))
  {
    cells_downsampled <- c()
    for (cluster in clusters)
    {
      cells_in_cluster <- table(clustering)[cluster]
      if (max.cells.per.cluster < cells_in_cluster)
      {
        set.seed(random.seed)
        inds <- sample(1:cells_in_cluster,size = max.cells.per.cluster,replace = FALSE)
        cells_downsampled <- c(cells_downsampled,names(clustering[clustering==cluster])[inds])
      }
    }
    data <- data[,cells_downsampled]
    clustering <- clustering[cells_downsampled]
  }

  # Compare cells from each cluster against all other clusters
  results_list <- list()

  # for (cluster in clusters)
  # {

  cluster <- "group.1"

  cat(paste0("testing cluster ",cluster,"\n"))
  # Extract data

  data_cluster <- data[,clustering==cluster]
  data_other <- data[,clustering!=cluster]

  # Skip if the number of cells in the test or the reference set is lower than min.cells.group
  if (ncol(data_cluster) < min.cells.group | ncol(data_other) < min.cells.group)
  {
    next
  }

  # min.pct filter
  genes.pct_cluster <- apply(data_cluster,1,function(x) sum(x!=0))/ncol(data_cluster)
  genes.pct_other <- apply(data_other,1,function(x) sum(x!=0))/ncol(data_other)

  genes_to_include <- rownames(data_cluster)[genes.pct_cluster>=min.pct | genes.pct_other >= min.pct]


  data_cluster <- data_cluster[genes_to_include,,drop=FALSE]
  data_other <- data_other[genes_to_include,,drop=FALSE]

  cat(paste0(nrow(data_cluster)," genes left after min.pct filtering\n"))
  if (nrow(data_cluster)==0)
  {
    cat("-----------------------------------\n")
    next
  }

  # min.diff.pct filter
  if (!is.null(min.diff.pct))
  {
    genes.pct_cluster <- genes.pct_cluster[genes_to_include]
    genes.pct_other <- genes.pct_other[genes_to_include]

    genes_to_include <- rownames(data_cluster)[abs(genes.pct_cluster-genes.pct_other) >= min.diff.pct]

    data_cluster <- data_cluster[genes_to_include,,drop=FALSE]
    data_other <- data_other[genes_to_include,,drop=FALSE]
  }

  cat(paste0(nrow(data_cluster)," genes left after min.diff.pct filtering\n"))
  if (nrow(data_cluster)==0)
  {
    cat("-----------------------------------\n")
    next
  }

  # logfc.threshold filter
  # Calculate log2 fold changes
  cluster_aves <- apply(data_cluster,1,mean)
  other_aves <- apply(data_other,1,mean)

  log2FC <- log2((cluster_aves+pseudocount.use)/(other_aves+pseudocount.use))

  genes_to_include <- rownames(data_cluster)[log2FC >= logfc.threshold | log2FC <= -logfc.threshold]

  data_cluster <- data_cluster[genes_to_include,,drop=FALSE]
  data_other <- data_other[genes_to_include,,drop=FALSE]

  cat(paste0(nrow(data_cluster)," genes left after logfc.threshold filtering\n"))
  if (nrow(data_cluster)==0)
  {
    cat("-----------------------------------\n")
    next
  }

  # Run DE test

  if (test=="wilcox")
  {
    wilcox.res <- lapply(rownames(data_cluster),function(x) wilcox.test(x=data_cluster[x,],y=data_other[x,]))
    p_values <- unlist(lapply(wilcox.res,function(x) x$p.value))
    names(p_values) <- rownames(data_cluster)

    # Adjust p-values
    adj_p_values <- p.adjust(p_values, method = "bonferroni", n = nrow(iloreg.object@normalized.data))

    res <- cbind(p_values,adj_p_values,log2FC[names(p_values)],genes.pct_cluster[names(p_values)],genes.pct_other[names(p_values)],abs(genes.pct_cluster-genes.pct_other)[names(p_values)])
    colnames(res) <- c("p.value","adj.p.value","log2FC","pct.1","pct.2","diff.pct")
    res <- as.data.frame(res)
    res$cluster <- cluster
    res$gene <- names(p_values)
  }

  results_list[[cluster]] <- res


  # }

  results_df <- do.call(rbind,results_list)
  rownames(results_df) <- make.unique(unlist(lapply(results_list,rownames)))

  results_df$cluster <- NULL

  if(only.pos) {
    results_df <- results_df[results_df$adj.p.value <= return.thresh,]
    return(results_df)
  }
  return(results_df)
})









setGeneric("VlnPlot", function(iloreg.object=NULL,
                               clustering.type="",
                               genes="",
                               return.plot=FALSE){
  standardGeneric("VlnPlot")
})

#' @title Gene expression visualization using violin plots
#'
#' @rdname VlnPlot
#' @name VlnPlot
#'
#' @description
#' The VlnPlot function enables visualizing expression levels of a gene, or multiple genes, across clusters using Violin plots.
#'
#'
#' @param iloreg.object object of class 'iloreg'
#' @param clustering.type "manual" or "optimal". "manual" refers to the clustering formed using the "SelectKClusters" function and "optimal" to the clustering formed using the "CalculateSilhouetteInformation" function. (default "manual")
#' @param genes a character vector denoting the gene names that are visualized
#' @param return.plot return.plot whether to return the ggplot2 object or just draw it (default \code{FALSE})
#' @return ggplot2 object if return.plot=TRUE
#' @keywords violin plot
#' @import ggplot2
#' @importFrom cowplot plot_grid

#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- RunTSNE(iloreg_object,perplexity=30)
#' iloreg_object <- HierarchicalClustering(iloreg_object)
#' iloreg_object <- SelectKClusters(iloreg_object,K=15)
#' VlnPlot(iloreg_object,genes = c("CD3D","CD8B","CD79A"))

setMethod("VlnPlot", "iloreg", function(iloreg.object,
                                        clustering.type,
                                        genes,
                                        return.plot)
{

  if (clustering.type=="manual")
  {
    clustering <- iloreg.object@clustering.manual
  } else if (clustering.type=="optimal")
  {
    clustering <- iloreg.object@clustering.optimal
  } else {
    clustering <- iloreg.object@clustering.manual
    cat("clustering.type='manual'")
  }

  data <- iloreg.object@normalized.data

  df <- as.numeric(t(data[genes,]))
  df <- data.frame(matrix(df,ncol = 1,dimnames = list(1:length(df),"Expression")))
  df$gene  <- unlist(lapply(genes,function(x) rep(x,ncol(data))))
  df$gene <- factor(df$gene)
  df$Cluster <- rep(as.character(clustering),length(genes))
  df$Cluster <- factor(df$Cluster)

  plotlist <- lapply(genes,function(x) ggplot(df[df$gene==x,], aes(x=Cluster, y=Expression, fill=Cluster))+geom_violin(trim=TRUE)+geom_jitter(height = 0, width = 0.1)+theme_classic()+ggtitle(x)+theme(plot.title = element_text(hjust = 0.5),legend.position = "none"))

  p <- plot_grid(plotlist = plotlist)
  print(p)

  if (return.plot)
  {
    return(p)
  } else {
    print(p)
  }



})




setGeneric("GeneHeatmap", function(iloreg.object=NULL,
                                   clustering.type="",
                                   gene.marker.data.frame=data.frame()){
  standardGeneric("GeneHeatmap")
})

#' @title Heatmap visualization of the gene markers identified by FindAllGeneMarkers
#'
#' @rdname GeneHeatmap
#' @name GeneHeatmap
#'
#' @description
#' The GeneHeatmap function enables drawing a heatmap of the gene markers identified by FindAllGeneMarkers, where the cell are grouped
#' by the clustering.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param clustering.type "manual" or "optimal". "manual" refers to the clustering formed using the "SelectKClusters" function and "optimal" to the clustering using the "CalculateSilhouetteInformation" function. (default "manual")
#' @param gene.marker.data.frame a data frame of the gene markers generated by FindAllGeneMarkers function. To accelerate the drawing, filtering the dataframe by selecting e.g. top 10 genes is recommended.
#' @return iloreg Object
#' @keywords iterative logistic regression ILoReg consensus clustering
#' @import pheatmap

#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- HierarchicalClustering(iloreg_object)
#' iloreg_object <- SelectKClusters(iloreg_object,K=15)
#' gene_markers <- FindAllGeneMarkers(iloreg_object,
#'                                   clustering.type = "manual",
#'                                   test = "wilcox",
#'                                   logfc.threshold = 0.25,
#'                                   min.pct = 0.1,
#'                                  min.diff.pct = NULL,
#'                                  pseudocount.use = 1,
#'                                  min.cells.group = 3,
#'                                 return.thresh = 0.01,
#'                                 only.pos = FALSE,max.cells.per.cluster = NULL)
#' library(dplyr)
#' ## Select top 10 genes per cluster by log2 fold-change
#' gene_markers %>% group_by(cluster) %>% top_n(10, log2FC) -> top10
#' gene_markers %>% group_by(cluster) %>% top_n(1, log2FC) -> top1
#'
#' GeneHeatmap(iloreg_object,clustering.type = "manual",gene.marker.data.frame = top10)
setMethod("GeneHeatmap", "iloreg", function(iloreg.object,
                                            clustering.type,
                                            gene.marker.data.frame)
{

  if (clustering.type=="manual")
  {
    clustering <- iloreg.object@clustering.manual
  } else if (clustering.type=="optimal")
  {
    clustering <- iloreg.object@clustering.optimal
  } else {
    clustering <- iloreg.object@clustering.manual
    cat("clustering.type='manual'")
  }

  data <- iloreg.object@normalized.data
  data <- data[unique(gene.marker.data.frame$gene),]
  # data <- scale(data,center = TRUE,scale = TRUE)
  data <- data[,order(clustering)]


  # Generate column annotations
  annotation = data.frame(cluster=sort(clustering))

  pheatmap(data,show_colnames = FALSE,
           gaps_col = cumsum(table(clustering[order(clustering)])),
           gaps_row = cumsum(table(gene.marker.data.frame[!duplicated(gene.marker.data.frame$gene),"cluster"])),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           annotation_col = annotation)



})








setGeneric("AnnotationScatterPlot", function(iloreg.object=NULL,
                                             annotation=c(),
                                             return.plot=FALSE,
                                             dim.reduction.type="",
                                             point.size=0.7){
  standardGeneric("AnnotationScatterPlot")
})

#' @title Visualiation of a custom annotation over nonlinear dimensionality reduction
#'
#' @rdname AnnotationScatterPlot
#' @name AnnotationScatterPlot
#'
#' @description
#' The AnnotationScatterPlot enables visualizing arbitrary class labels over the nonliner dimensionality reduction, e.g. t-SNE or UMAP.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param annotation a character vector, factor or numeric for the class labels.
#' @param return.plot return.plot whether to return the ggplot2 object or just draw it (default \code{FALSE})
#' @param dim.reduction.type "tsne" or "umap" (default "tsne")
#' @param point.size point size (default 0.7)
#' @return ggplot2 object if return.plot=TRUE
#' @keywords annotation custom visualization t-sne umap nonlinear dimensionality reduction
#' @import ggplot2

#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- RunTSNE(iloreg_object,perplexity=30)
#' AnnotationScatterPlot(iloreg_object,annotation=c(),point.size=0.5,dim.reduction.type="tsne")
setMethod("AnnotationScatterPlot", "iloreg", function(iloreg.object,
                                                      annotation,
                                                      return.plot,
                                                      dim.reduction.type,
                                                      point.size)
{

  if (dim.reduction.type=="umap")
  {
    two.dim.data <- iloreg.object@umap.embeddings
    xlab <- "UMAP_1"
    ylab <- "UMAP_2"
  } else if (dim.reduction.type=="tsne"){
    two.dim.data <- iloreg.object@tsne.embeddings
    xlab <- "tSNE_1"
    ylab <- "tSNE_2"
  } else {
    stop("dim.reduction.type must be either 'tsne' or 'umap'")
  }

  annotation <- factor(as.character(annotation))
  names(annotation) <- colnames(iloreg.object@normalized.data)

  df <- as.data.frame(two.dim.data)
  df$group <- annotation
  colnames(df) <- c("dim1","dim2","group")

  two.dim.data_ <- two.dim.data
  rownames(two.dim.data_) <- names(annotation)
  cluster_centers <- lapply(levels(annotation),function(x) apply(two.dim.data_[names(annotation)[annotation==x],,drop=FALSE],2,median))
  cluster_centers <- do.call(rbind,cluster_centers)


  p<-ggplot(df, aes(x=dim1, y=dim2)) +
    geom_point(size=point.size,aes(color=group)) +
    xlab(xlab) +
    ylab(ylab) +
    theme_classic() +
    annotate("text", x = cluster_centers[,1], y = cluster_centers[,2], label = levels(annotation)) +
    guides(colour = guide_legend(override.aes = list(size=2)))

  if (return.plot)
  {
    return(p)
  } else {
    print(p)
  }




})






setGeneric("GeneDropoutRatePlot", function(iloreg.object=NULL,
                                           genes="",
                                           return.plot=FALSE,use.clusters=NULL,clustering.type="manual"){
  standardGeneric("GeneDropoutRatePlot")
})

#' @title Dropout rate versus average non-zero expression
#'
#'
#' @rdname GeneDropoutRatePlot
#' @name GeneDropoutRatePlot
#'
#' @description
#' The GeneDropoutRatePlot function enables visualizion of dropout rates (fraction of cells expressing a gene) against
#' the average non-zero expression values. This function can aid the user to determinine if a gene is differentially expressed or not,
#' since differentially expressed genes are likely to be above the expected curve. Cells can be filtered based on clustering before
#' calculating the values.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param genes a character vector of the genes to visualize over the dropout rate curve
#' @param return.plot return.plot whether to return the ggplot2 object or just draw it (default \code{FALSE})
#' @param use.clusters use data from these clusters only
#' @param clustering.type "manual" or "optimal". "manual" refers to the clustering formed using the "SelectKClusters" function and "optimal" to the clustering using the "CalculateSilhouetteInformation" function. (default "manual")
#' @return ggplot2 object if return.plot=TRUE
#' @keywords dropout rate curve non-zero average expression gene
#' @import ggplot2
#' @importFrom reshape2 melt

#' @export
#' @examples
#' ## Load 10X Chromium data and LogNormalized them.
#' raw_data <- Seurat::Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- Seurat::LogNormalize(raw_data)
#' ## Initialize the iloreg object with a normalized gene expression dataset.
#' iloreg_object <- CreateILoRegObject(normalized.data=data)
#' iloreg_object <- RunParallelICP(iloreg_object,threads = 12,seed=1)
#' VisualizeQC(iloreg_object)
#' iloreg_object <- RunPCA(iloreg_object,p = 50,scale = FALSE)
#' PCAElbowPlot(iloreg_object,return.plot=FALSE)
#' iloreg_object <- HierarchicalClustering(iloreg_object)
#' iloreg_object <- SelectKClusters(iloreg_object,K=15)
#' GeneDropoutRatePlot(iloreg_object,genes = c("CD79A","TCL1A","IGLL5","VPREB3"),use.clusters = c(24,15,14,2),clustering.type = "manual")
setMethod("GeneDropoutRatePlot", "iloreg", function(iloreg.object,
                                                    genes,
                                                    return.plot,use.clusters,clustering.type)
{

  if (is.null(use.clusters))
  {
    normalized_data <- iloreg.object@normalized.data

  } else {

    if (clustering.type=="manual")
    {
      clustering <- iloreg.object@clustering.manual
    } else if (clustering.type=="optimal") {
      clustering <- iloreg.object@clustering.optimal
    } else {
      clustering <- iloreg.object@clustering.manual
      cat("clustering.type='manual'")
    }

    cells.to.use <- names(clustering)[clustering %in% use.clusters]
    print(length(cells.to.use))

    normalized_data <- iloreg.object@normalized.data[,cells.to.use]

  }

  dropout_rates <- apply(normalized_data,1,function(x) sum(x==0))/ncol(normalized_data)

  average_nonzero_expression <- apply(normalized_data,1,function(x) mean(x[x!=0]))

  df <- melt(dropout_rates)
  df$average_nonzero_expression <- average_nonzero_expression
  df$logical <- rownames(df) %in% genes


  p <- ggplot(data=df, aes(x=average_nonzero_expression, y=value,color=logical)) +
    geom_point()+
    theme_bw()+
    scale_colour_discrete(name="gene",labels=c("FALSE","TRUE")) +
    ylab("Dropout rate")+
    xlab("Average non-zero expression")+
    theme(legend.position = "none")

  genes <- genes[order(df[genes,"value"])]

  dropout_previous_gene <- NA
  for (gene in genes)
  {
    text_annotate_x_space <- 0.25
    if (nchar(gene) > 7)
    {
      text_annotate_x_space <- 0.5
    }
    seg_length <- runif(1,0.1,2)
    if (is.na(dropout_previous_gene))
    {
      p <- p + annotate("segment", x = df[gene,"average_nonzero_expression"], xend = df[gene,"average_nonzero_expression"]+seg_length, y = df[gene,"value"], yend = df[gene,"value"], colour = "black") +
        annotate("text", x = df[gene,"average_nonzero_expression"]+seg_length+text_annotate_x_space, y = df[gene,"value"], label = gene)
    } else {
      if ((dropout_previous_gene - df[gene,"value"]) < 0.05)
      {
        random_y_space <- runif(1,-0.25,0.25)
        p <- p + annotate("segment", x = df[gene,"average_nonzero_expression"], xend = df[gene,"average_nonzero_expression"]+seg_length, y = df[gene,"value"], yend = df[gene,"value"]+random_y_space, colour = "black") +
          annotate("text", x = df[gene,"average_nonzero_expression"]+seg_length+text_annotate_x_space, y = df[gene,"value"]+random_y_space, label = gene)

      } else {
        p <- p + annotate("segment", x = df[gene,"average_nonzero_expression"], xend = df[gene,"average_nonzero_expression"]+seg_length, y = df[gene,"value"], yend = df[gene,"value"], colour = "black") +
          annotate("text", x = df[gene,"average_nonzero_expression"]+seg_length+text_annotate_x_space, y = df[gene,"value"], label = gene)
      }
    }
    dropout_previous_gene <- df[gene,"value"]
  }


  if (return.plot)
  {
    return(p)
  } else {
    print(p)
  }


})



