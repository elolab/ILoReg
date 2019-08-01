#' @title Iterative Logistic Regression (ILoReg) consensus clustering
#'
#' @description R package that enables supervised learning -based clustering through L1-regularized (LASSO) logistic regression.
#' Consensus clustering is performed by running ILoReg multiple times. PCA is used to reduce dimensionality of the
#' consensus probability matrix. The Ward's method is used to perform the final clustering.
#' The silhouette method can be used to choose the optimal number of clusters (K). Moreover, t-SNE or UMAP can be performed
#' onto the PCA-transformed data to visualize the results.
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
#' @param normalized.data a gene expression matrix with genes in rows and cells in columns...
#'
#' @importFrom Matrix rowSums
#' @export
#'
#' @examples
#' ## Generate simulated single-cell RNA-Seq tags.
#' library(Seurat)
#' raw_data <- Read10x("~/../Downloads/10x_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' data <- LogNormalize(raw_data)
#' iloreg <- CreateILoRegObject(normalized.data=data)
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

setGeneric("RunILoRegConsensus", function(iloreg.object=NULL,k=15,d=0.3,L=100,r=5,C=0.5,type="L1",max.number.of.iterations=200,threads=0){
  standardGeneric("RunILoRegConsensus")
})

#' @title Iterative logistic regression (ILoReg) consensus clustering
#'
#' @rdname RunILoRegConsensus
#' @name RunILoRegConsensus
#'
#' @description
#' Enables iterative logistic regression (ILoReg) consensus clustering with parallel computation.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#' @param iloreg.object object of class 'iloreg'
#' @param k A positive integer greater or equal to 2, which denotes the number of clusters in iterative logistic regression. If stratified.downsampling=TRUE, the final k is likely to be smaller. (default 15)
#' @param d A numeric that defines how many cells per cluster should be down- and oversampled (d in N/k*d), when stratified.downsampling=FALSE,  or what fraction should be downsampled in the stratified approach ,stratified.downsampling=TRUE. (default 0.3)
#' @param L A positive integer that specifies how many ILoReg runs in the consensus method are ran. (default 50)
#' @param r A positive integer that denotes the maximum number of reiterations performed until the algorithm ends. (default 5)
#' @param C Cost of constraints (C). Recommendation is in range [0.2,0.5] for LASSO.  (default 0.3)
#' @param type "L1" for LASSO and "L2" for Ridge.  (default "L1")
#' @param threads A positive integer that specifies how many logical processors (threads) should be used. Use threads=1 to disable parallelism and threads=0 to use all available core minus one. (default 0)
#' @param max.number.of.iterations A positive integer that denotes the maximum number of iterations performed until the algorithm ends. (default 200)
#' @return iloreg object
#' @keywords iterative logistic regression ILoReg consensus clustering
#' @import foreach
#' @import doParallel
#' @import ParallelLogger
#' @import parallel
#' @import Matrix
#' @import tictoc
#' @import aricode
#' @import LiblineaR
#' @import SparseM
#' @export
#' @examples
#' a <- c(0,1,2)
setMethod("RunILoRegConsensus", "iloreg", function(iloreg.object, k,d,L,r,C,type,max.number.of.iterations,threads){

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

    iloreg_out <- foreach(clusters = rep(k,L),.verbose = FALSE,
                                     .combine = list,
                                     .maxcombine = 1000,
                                     .inorder = FALSE,
                                     .export = c('RunSingleILoRegClustering','LogisticRegression'),
                                     .packages=c("tictoc","Matrix","aricode","LiblineaR","SparseM"),
                                     .multicombine = TRUE)  %dopar%
      RunSingleILoRegClustering(normalized.data = iloreg.object@normalized.data,
                                k = iloreg.object@k,
                                d = iloreg.object@d,
                                r = iloreg.object@r,
                                C = iloreg.object@C,
                                type = type,
                                max.number.of.iterations = max.number.of.iterations)
    stopCluster(cl)

  } else {
    iloreg_out <- list()
    for (l in 1:L) {
      res <- RunSingleILoRegClustering(normalized.data = iloreg.object@normalized.data,
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



setGeneric("VisualizeConsensusInformation", function(iloreg.object=NULL,return.plot=FALSE){
  standardGeneric("VisualizeConsensusInformation")
})

#' @title Iterative logistic regression (ILoReg) consensus clustering
#'
#' @rdname VisualizeConsensusInformation
#' @name VisualizeConsensusInformation
#'
#' @description
#' Plot histogram of the ARIs.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param return.plot logical indicating if the ggplot2 object should be returned (default FALSE)
#' @return ggplot2 object if return.plot=TRUE
#' @keywords iterative logistic regression ILoReg consensus clustering
#' @import ggplot2
#' @import cowplot
#' @export
#' @examples
#' a <- c(0,1,2)
setMethod("VisualizeConsensusInformation", "iloreg", function(iloreg.object,return.plot){

  final_aris <- unlist(lapply(iloreg.object@metrics,function(x) x["ARI",ncol(x)]))

  df <- data.frame(matrix(final_aris,ncol = 1,dimnames = list(1:length(final_aris),c("ARI"))))

  p1 <- ggplot(df, aes(x=ARI))+
    geom_histogram(color="darkblue", fill="lightblue") +
    xlim(0,1)

  number_of_runs <- unlist(lapply(iloreg.object@metrics,function(x) ncol(x)))

  df <- data.frame(matrix(number_of_runs,ncol = 1,dimnames = list(1:length(number_of_runs),c("Runs"))))

  p2 <- ggplot(df, aes(x=Runs))+
    geom_histogram(color="darkblue", fill="lightblue") +
    xlab("Number of iterations")

  p <- plot_grid(p1,p2,nrow=1)

  print(p)

  if (return.plot) {
    return(p)
  }
})



setGeneric("RunPCA", function(iloreg.object=NULL,number.of.pcs=50,scale=FALSE){
  standardGeneric("RunPCA")
})

#' @title Iterative logistic regression (ILoReg) consensus clustering
#'
#' @rdname RunPCA
#' @name RunPCA
#'
#' @description
#' Plot histogram of the ARIs.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param number.of.pcs logical indicating if the ggplot2 object should be returned (default FALSE)
#' @param scale logical indicating if the ggplot2 object should be returned (default FALSE)
#' @return ggplot2 object if return.plot=TRUE
#' @keywords iterative logistic regression ILoReg consensus clustering
#' @importFrom RSpectra eigs_sym
#' @export
#' @examples
#' a <- c(0,1,2)
setMethod("RunPCA", "iloreg", function(iloreg.object,number.of.pcs,scale){

  iloreg.object@number.of.pcs <- number.of.pcs
  iloreg.object@scale.pca <- scale

  X <- do.call(cbind,iloreg.object@consensus.probability)

  X <- scale(X,scale = scale,center = TRUE)

  # X^T %*% X
  A = crossprod(X)

  # Perform eigendecomposition
  eigs_sym_out <- eigs_sym(A, number.of.pcs, which = "LM")

  rotated <- X %*% eigs_sym_out$vectors
  colnames(rotated) <- paste0("PC",1:ncol(rotated))

  iloreg.object@rotated.consensus.probability <- rotated

  return(iloreg.object)

})



setGeneric("RunUMAP", function(iloreg.object=NULL){
  standardGeneric("RunUMAP")
})

#' @title Iterative logistic regression (ILoReg) consensus clustering
#'
#' @rdname RunUMAP
#' @name RunUMAP
#'
#' @description
#' Plot histogram of the ARIs.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#'
#' @param iloreg.object object of class 'iloreg'
#' @return ggplot2 object if return.plot=TRUE
#' @keywords iterative logistic regression ILoReg consensus clustering
#' @importFrom umap umap
#' @export
#' @examples
#' a <- c(0,1,2)
setMethod("RunUMAP", "iloreg", function(iloreg.object){

  umap_out <- umap(iloreg.object@rotated.consensus.probability)

  iloreg.object@umap.embeddings <- umap_out$layout

  return(iloreg.object)

})



setGeneric("RunTSNE", function(iloreg.object=NULL,perplexity=30){
  standardGeneric("RunTSNE")
})

#' @title Iterative logistic regression (ILoReg) consensus clustering
#'
#' @rdname RunTSNE
#' @name RunTSNE
#'
#' @description
#' Plot histogram of the ARIs.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param perplexity perplexity
#' @return iloreg Object
#' @keywords iterative logistic regression ILoReg consensus clustering
#' @importFrom Rtsne Rtsne
#' @export
#' @examples
#' a <- c(0,1,2)
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

#' @title Iterative logistic regression (ILoReg) consensus clustering
#'
#' @rdname HierarchicalClustering
#' @name HierarchicalClustering
#'
#' @description
#' Plot histogram of the ARIs.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#'
#' @param iloreg.object object of class 'iloreg'
#' @return iloreg Object
#' @keywords iterative logistic regression ILoReg consensus clustering
#' @importFrom fastcluster hclust.vector
#' @export
#' @examples
#' a <- c(0,1,2)
setMethod("HierarchicalClustering", "iloreg", function(iloreg.object){


  hc <- hclust.vector(iloreg.object@rotated.consensus.probability,
                      method = "ward")

  iloreg.object@hc <- hc

  return(iloreg.object)

})






setGeneric("CalculateSilhouetteInformation", function(iloreg.object=NULL,k.range=2:50){
  standardGeneric("CalculateSilhouetteInformation")
})

#' @title Iterative logistic regression (ILoReg) consensus clustering
#'
#' @rdname CalculateSilhouetteInformation
#' @name CalculateSilhouetteInformation
#'
#' @description
#' Plot histogram of the ARIs.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#'
#' @param iloreg.object object of class 'iloreg'
#' @return iloreg Object
#' @keywords iterative logistic regression ILoReg consensus clustering
#' @importFrom parallelDist parDist
#' @importFrom cluster silhouette
#' @importFrom dendextend cutree

#' @export
#' @examples
#' a <- c(0,1,2)
setMethod("CalculateSilhouetteInformation", "iloreg", function(iloreg.object,k.range){

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

#' @title Iterative logistic regression (ILoReg) consensus clustering
#'
#' @rdname SilhouetteCurve
#' @name SilhouetteCurve
#'
#' @description
#' Plot histogram of the ARIs.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param return.plot object of class 'iloreg'
#' @return iloreg Object
#' @keywords iterative logistic regression ILoReg consensus clustering
#' @import ggplot2
#' @importFrom DescTools AUC

#' @export
#' @examples
#' a <- c(0,1,2)
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

  print(p)

  if (return.plot)
  {
    return(p)
  }




})




setGeneric("ManualClustering", function(iloreg.object=NULL,K=10){
  standardGeneric("ManualClustering")
})

#' @title Iterative logistic regression (ILoReg) consensus clustering
#'
#' @rdname ManualClustering
#' @name ManualClustering
#'
#' @description
#' Plot histogram of the ARIs.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#'
#' @param iloreg.object object of class 'iloreg'
#' @return iloreg Object
#' @keywords iterative logistic regression ILoReg consensus clustering
#' @importFrom dendextend cutree

#' @export
#' @examples
#' a <- c(0,1,2)
setMethod("ManualClustering", "iloreg", function(iloreg.object,K){

  clustering <- factor(cutree(as.dendrogram(iloreg.object@hc),k=K))
  names(clustering) <- colnames(iloreg.object@normalized.data)

  iloreg.object@clustering.manual <- clustering
  iloreg.object@K.manual <- K

  return(iloreg.object)

})


setGeneric("MergeClusters", function(iloreg.object=NULL,clusters.to.merge=""){
  standardGeneric("MergeClusters")
})

#' @title Iterative logistic regression (ILoReg) consensus clustering
#'
#' @rdname MergeClusters
#' @name MergeClusters
#'
#' @description
#' Plot histogram of the ARIs.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param clusters.to.merge object of class 'iloreg'
#' @return iloreg Object
#' @keywords iterative logistic regression ILoReg consensus clustering

#' @export
#' @examples
#' a <- c(0,1,2)
setMethod("MergeClusters", "iloreg", function(iloreg.object,clusters.to.merge){

  clusters.to.merge <- as.character(clusters.to.merge)

  clustering_old <- iloreg.object@clustering.manual
  clusters_old <- levels(clustering_old)

  if (sum(clusters.to.merge %in% clusters_old)!=length(clusters.to.merge))
  {
    stop("invalid clusters.to.merge argument")
  }

  new_cluster_name <- paste(clusters.to.merge,collapse = ",")

  clustering_new <- as.character(clustering_old)
  clustering_new[clustering_new %in% clusters.to.merge] <- new_cluster_name

  clustering_new <- factor(clustering_new)
  names(clustering_new) <- names(clustering_old)

  iloreg.object@clustering.manual <- clustering_new
  iloreg.object@K.manual <- length(levels(clustering_new))

  return(iloreg.object)

})



setGeneric("RenameClusters", function(iloreg.object=NULL,new.cluster.names=""){
  standardGeneric("RenameClusters")
})

#' @title Iterative logistic regression (ILoReg) consensus clustering
#'
#' @rdname RenameClusters
#' @name RenameClusters
#'
#' @description
#' Plot histogram of the ARIs.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param new.cluster.names object of class 'iloreg'
#' @return iloreg Object
#' @keywords iterative logistic regression ILoReg consensus clustering
#' @importFrom plyr mapvalues

#' @export
#' @examples
#' a <- c(0,1,2)
setMethod("RenameClusters", "iloreg", function(iloreg.object,new.cluster.names){

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




setGeneric("GeneScatterPlot", function(iloreg.object=NULL,genes="",return.plot=FALSE,dim.reduction.type="",point.size=0.7){
  standardGeneric("GeneScatterPlot")
})

#' @title Iterative logistic regression (ILoReg) consensus clustering
#'
#' @rdname GeneScatterPlot
#' @name GeneScatterPlot
#'
#' @description
#' Plot histogram of the ARIs.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param genes object of class 'iloreg'
#' @param return.plot object of class 'iloreg'
#' @param dim.reduction.type object of class 'iloreg'
#' @param point.size object of class 'iloreg'
#' @return iloreg Object
#' @keywords iterative logistic regression ILoReg consensus clustering
#' @import ggplot2
#' @importFrom scales muted

#' @export
#' @examples
#' a <- c(0,1,2)
setMethod("GeneScatterPlot", "iloreg", function(iloreg.object,genes,return.plot,dim.reduction.type,point.size){

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

    p<-ggplot(df, aes(x=dim1, y=dim2)) +
      geom_point(size=point.size,aes(color=group)) +
      scale_colour_gradient2(low = muted("red"), mid = "lightgrey",
                             high = muted("blue"),name = genes) +
      xlab(xlab) +
      ylab(ylab) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))


    if (return.plot) {
      return(p)
    }

    print(p)

  }
})





setGeneric("ClusterScatterPlot", function(iloreg.object=NULL,clustering.type="",return.plot=FALSE,dim.reduction.type="",point.size=0.7){
  standardGeneric("ClusterScatterPlot")
})

#' @title Iterative logistic regression (ILoReg) consensus clustering
#'
#' @rdname ClusterScatterPlot
#' @name ClusterScatterPlot
#'
#' @description
#' Plot histogram of the ARIs.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param clustering.type object of class 'iloreg'
#' @param return.plot object of class 'iloreg'
#' @param dim.reduction.type object of class 'iloreg'
#' @param point.size object of class 'iloreg'
#' @return iloreg Object
#' @keywords iterative logistic regression ILoReg consensus clustering
#' @import ggplot2

#' @export
#' @examples
#' a <- c(0,1,2)
setMethod("ClusterScatterPlot", "iloreg", function(iloreg.object,clustering.type,return.plot,dim.reduction.type,point.size){

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
    stop("clustering.type must be 'optimal' or 'manual'")
  }

  df <- as.data.frame(two.dim.data)

  df$cluster <- color.by
  colnames(df) <- c("dim1","dim2","cluster")

  two.dim.data_ <- two.dim.data
  rownames(two.dim.data_) <- names(color.by)
  cluster_centers <- lapply(levels(color.by),function(x) apply(two.dim.data_[names(color.by)[color.by==x],,drop=FALSE],2,median))
  cluster_centers <- do.call(rbind,cluster_centers)

  p<-ggplot(df, aes(x=dim1, y=dim2)) +
    geom_point(size=point.size,aes(color=cluster)) +
    xlab(xlab) +
    ylab(ylab) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    annotate("text", x = cluster_centers[,1], y = cluster_centers[,2], label = levels(color.by))
  plot(p)

  if (return.plot) {
    return(p)
  }

})





setGeneric("FindAllGeneMarkers", function(iloreg.object=NULL,clustering.type="",test="wilcox",
                                          logfc.threshold = 0.25,
                                          min.pct = 0.1,
                                          min.diff.pct = NULL,
                                          min.cells.cluster = 3,
                                          max.cells.per.cluster = NULL,
                                          random.seed = 1,
                                          pseudocount.use = 1,
                                          return.thresh = 0.01,
                                          only.pos=FALSE){
  standardGeneric("FindAllGeneMarkers")
})

#' @title Iterative logistic regression (ILoReg) consensus clustering
#'
#' @rdname FindAllGeneMarkers
#' @name FindAllGeneMarkers
#'
#' @description
#' Plot histogram of the ARIs.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param clustering.type object of class 'iloreg'
#' @param test object of class 'iloreg'
#' @param logfc.threshold object of class 'iloreg'
#' @param min.pct object of class 'iloreg'
#' @param min.diff.pct object of class 'iloreg'
#' @param min.cells.cluster object of class 'iloreg'
#' @param max.cells.per.cluster object of class 'iloreg'
#' @param random.seed object of class 'iloreg'
#' @param pseudocount.use object of class 'iloreg'
#' @param return.thresh object of class 'iloreg'
#' @param only.pos object of class 'iloreg'
#' @return iloreg Object
#' @keywords iterative logistic regression ILoReg consensus clustering

#' @export
#' @examples
#' a <- c(0,1,2)
setMethod("FindAllGeneMarkers", "iloreg", function(iloreg.object,
                                                   clustering.type,
                                                   test,
                                                   logfc.threshold,
                                                   min.pct,
                                                   min.diff.pct,
                                                   min.cells.cluster,
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
    stop("clustering.type must be 'optimal' or 'manual'")
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
    cat(paste0("testing cluster ",cluster,"\n"))
    # Extract data

    data_cluster <- data[,clustering==cluster]
    data_other <- data[,clustering!=cluster]

    # Skip if the number of cells in the test or the reference set is lower than min.cells.group
    if (ncol(data_cluster) < min.cells.cluster | ncol(data_other) < min.cells.cluster)
    {
      next
    }

    # min.pct filter
    genes.pct_cluster <- apply(data_cluster,1,function(x) sum(x!=0))/ncol(data_cluster)
    genes.pct_other <- apply(data_other,1,function(x) sum(x!=0))/ncol(data_other)

    genes_to_include <- rownames(data_cluster)[genes.pct_cluster>=min.pct | genes.pct_other >= min.pct]

    data_cluster <- data_cluster[genes_to_include,]
    data_other <- data_other[genes_to_include,]

    # min.diff.pct filter
    if (!is.null(min.diff.pct))
    {
      genes.pct_cluster <- genes.pct_cluster[genes_to_include]
      genes.pct_other <- genes.pct_other[genes_to_include]

      genes_to_include <- rownames(data_cluster)[abs(genes.pct_cluster-genes.pct_other) >= min.diff.pct]

      data_cluster <- data_cluster[genes_to_include,]
      data_other <- data_other[genes_to_include,]
    }

    # logfc.threshold filter
    # Calculate log2 fold changes
    cluster_aves <- apply(data_cluster,1,mean)
    other_aves <- apply(data_other,1,mean)

    log2FC <- log2((cluster_aves+pseudocount.use)/(other_aves+pseudocount.use))

    genes_to_include <- rownames(data_cluster)[log2FC >= logfc.threshold]

    data_cluster <- data_cluster[genes_to_include,]
    data_other <- data_other[genes_to_include,]

    # Run DE test

    if (test=="wilcox")
    {
      wilcox.res <- lapply(rownames(data_cluster),function(x) wilcox.test(x=data_cluster[x,],y=data_other[x,]))
      p_values <- unlist(lapply(wilcox.res,function(x) x$p.value))
      names(p_values) <- rownames(data_cluster)

      # Adjust p-values
      adj_p_values <- p.adjust(p_values, method = "BH", n = length(p_values))

      res <- cbind(p_values,adj_p_values,log2FC[names(p_values)],genes.pct_cluster[names(p_values)],genes.pct_other[names(p_values)],abs(genes.pct_cluster-genes.pct_other)[names(p_values)])
      colnames(res) <- c("p.value","adj.p.value","log2FC","pct.1","pct.2","diff.pct")
      res <- as.data.frame(res)
      res$cluster <- cluster
      res$gene <- names(p_values)
    }

    results_list[[cluster]] <- res


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
                                   min.cells.cluster = 3,
                                   max.cells.per.cluster = NULL,
                                   random.seed = 1,
                                   pseudocount.use = 1,
                                   return.thresh = 0.01,
                                   only.pos=FALSE){
  standardGeneric("FindGeneMarkers")
})

#' @title Iterative logistic regression (ILoReg) consensus clustering
#'
#' @rdname FindGeneMarkers
#' @name FindGeneMarkers
#'
#' @description
#' Plot histogram of the ARIs.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param clustering.type object of class 'iloreg'
#' @param test object of class 'iloreg'
#' @param logfc.threshold object of class 'iloreg'
#' @param min.pct object of class 'iloreg'
#' @param min.diff.pct object of class 'iloreg'
#' @param min.cells.cluster object of class 'iloreg'
#' @param max.cells.per.cluster object of class 'iloreg'
#' @param random.seed object of class 'iloreg'
#' @param pseudocount.use object of class 'iloreg'
#' @param return.thresh object of class 'iloreg'
#' @param only.pos object of class 'iloreg'
#' @return iloreg Object
#' @keywords iterative logistic regression ILoReg consensus clustering

#' @export
#' @examples
#' a <- c(0,1,2)
setMethod("FindGeneMarkers", "iloreg", function(iloreg.object,
                                            clusters.1,
                                            clusters.2,
                                            clustering.type,
                                            test,
                                            logfc.threshold,
                                            min.pct,
                                            min.diff.pct,
                                            min.cells.cluster,
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
    stop("clustering.type must be 'optimal' or 'manual'")
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

  for (cluster in clusters)
  {
    cat(paste0("testing cluster ",cluster,"\n"))
    # Extract data

    data_cluster <- data[,clustering==cluster]
    data_other <- data[,clustering!=cluster]

    # Skip if the number of cells in the test or the reference set is lower than min.cells.group
    if (ncol(data_cluster) < min.cells.cluster | ncol(data_other) < min.cells.cluster)
    {
      next
    }

    # min.pct filter
    genes.pct_cluster <- apply(data_cluster,1,function(x) sum(x!=0))/ncol(data_cluster)
    genes.pct_other <- apply(data_other,1,function(x) sum(x!=0))/ncol(data_other)

    genes_to_include <- rownames(data_cluster)[genes.pct_cluster>=min.pct | genes.pct_other >= min.pct]

    data_cluster <- data_cluster[genes_to_include,]
    data_other <- data_other[genes_to_include,]

    # min.diff.pct filter
    if (!is.null(min.diff.pct))
    {
      genes.pct_cluster <- genes.pct_cluster[genes_to_include]
      genes.pct_other <- genes.pct_other[genes_to_include]

      genes_to_include <- rownames(data_cluster)[abs(genes.pct_cluster-genes.pct_other) >= min.diff.pct]

      data_cluster <- data_cluster[genes_to_include,]
      data_other <- data_other[genes_to_include,]
    }

    # logfc.threshold filter
    # Calculate log2 fold changes
    cluster_aves <- apply(data_cluster,1,mean)
    other_aves <- apply(data_other,1,mean)

    log2FC <- log2((cluster_aves+pseudocount.use)/(other_aves+pseudocount.use))

    genes_to_include <- rownames(data_cluster)[log2FC >= logfc.threshold]

    data_cluster <- data_cluster[genes_to_include,]
    data_other <- data_other[genes_to_include,]

    # Run DE test

    if (test=="wilcox")
    {
      wilcox.res <- lapply(rownames(data_cluster),function(x) wilcox.test(x=data_cluster[x,],y=data_other[x,]))
      p_values <- unlist(lapply(wilcox.res,function(x) x$p.value))
      names(p_values) <- rownames(data_cluster)

      # Adjust p-values
      adj_p_values <- p.adjust(p_values, method = "BH", n = length(p_values))

      res <- cbind(p_values,adj_p_values,log2FC[names(p_values)],genes.pct_cluster[names(p_values)],genes.pct_other[names(p_values)],abs(genes.pct_cluster-genes.pct_other)[names(p_values)])
      colnames(res) <- c("p.value","adj.p.value","log2FC","pct.1","pct.2","diff.pct")
      res <- as.data.frame(res)
      res$cluster <- cluster
      res$gene <- names(p_values)
    }

    results_list[[cluster]] <- res


  }

  results_df <- do.call(rbind,results_list)
  rownames(results_df) <- make.unique(unlist(lapply(results_list,rownames)))

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

#' @title Iterative logistic regression (ILoReg) consensus clustering
#'
#' @rdname VlnPlot
#' @name VlnPlot
#'
#' @description
#' Plot histogram of the ARIs.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param clustering.type object of class 'iloreg'
#' @param genes object of class 'iloreg'
#' @param return.plot object of class 'iloreg'
#' @return iloreg Object
#' @keywords iterative logistic regression ILoReg consensus clustering
#' @import ggplot2
#' @importFrom cowplot plot_grid

#' @export
#' @examples
#' a <- c(0,1,2)
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
    stop("clustering.type must be 'optimal' or 'manual'")
  }

  data <- iloreg.object@normalized.data

  df <- as.numeric(t(data[genes,]))
  df <- data.frame(matrix(df,ncol = 1,dimnames = list(1:length(df),"expression")))
  df$gene  <- unlist(lapply(genes,function(x) rep(x,ncol(data))))
  df$gene <- factor(df$gene)
  df$cluster <- rep(as.character(clustering),length(genes))
  df$cluster <- factor(df$cluster)

  plotlist <- lapply(genes,function(x) ggplot(df[df$gene==x,], aes(x=cluster, y=expression, fill=cluster))+geom_violin(trim=FALSE)+geom_jitter(height = 0, width = 0.1)+ggtitle(x))

  p <- plot_grid(plotlist = plotlist)
  print(p)

  if (return.plot)
  {
    return(p)
  }



})




setGeneric("GeneHeatmap", function(iloreg.object=NULL,
                               clustering.type="",
                               genes="",
                               return.plot=FALSE){
  standardGeneric("GeneHeatmap")
})

#' @title Iterative logistic regression (ILoReg) consensus clustering
#'
#' @rdname GeneHeatmap
#' @name GeneHeatmap
#'
#' @description
#' Plot histogram of the ARIs.
#'
#' @details
#' populates a Boolean matrix with the same dimension as nData.
#' The value is \code{TRUE} for an entry if it
#' is a dropout candidate; otherwise the value is \code{FALSE}.
#'
#' @param iloreg.object object of class 'iloreg'
#' @param clustering.type object of class 'iloreg'
#' @param genes object of class 'iloreg'
#' @param return.plot object of class 'iloreg'
#' @return iloreg Object
#' @keywords iterative logistic regression ILoReg consensus clustering
#' @import ggplot2
#' @importFrom grid textGrob gpar
#' @importFrom RColorBrewer brewer.pal

#' @export
#' @examples
#' a <- c(0,1,2)
setMethod("GeneHeatmap", "iloreg", function(iloreg.object,
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
    stop("clustering.type must be 'optimal' or 'manual'")
  }

  data <- iloreg.object@normalized.data
  data <- t(data[genes,])
  data <- scale(data,center = TRUE,scale = TRUE)
  data <- data[order(clustering),]

  # Clip values
  data[data < -2.5] <- -2.5
  data[data > 2.5] <- 2.5


  df <- as.numeric(data)
  df <- data.frame(matrix(df,ncol = 1,dimnames = list(1:length(df),"expression")))
  df$gene  <- unlist(lapply(genes,function(x) rep(x,nrow(data))))
  df$gene <- factor(df$gene)
  df$cluster <- rep(rownames(data),length(genes))
  df$cluster <- factor(df$cluster)


  text_high <- textGrob("Highest\nvalue", gp=gpar(fontsize=13, fontface="bold"))
  text_low <- textGrob("Lowest\nvalue", gp=gpar(fontsize=13, fontface="bold"))

  p <- ggplot(df, aes(cluster,gene)) +
    geom_tile(aes(fill = expression), colour = NA) +
    scale_fill_gradientn(colours=brewer.pal(11,"RdYlBu")) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_text(),
          axis.ticks.x = element_blank(),
          legend.position = "top",
          plot.margin = unit(c(1,1,2,1), "lines")) +
    annotation_custom(text_high,xmin=100,xmax=100,ymin=-0.07,ymax=-0.07) +
    annotation_custom(text_low,xmin=200,xmax=200,ymin=-0.07,ymax=-0.07) +
    coord_cartesian(clip = "off")

  print(p)

  if (return.plot)
  {
    return(p)
  }



})

