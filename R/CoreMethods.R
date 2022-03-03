#' @title Prepare \code{SingleCellExperiment} object for \code{ILoReg} analysis
#'
#' @description
#' This function prepares the \code{SingleCellExperiment} object for
#' \code{ILoReg} analysis. The only required input is an object of class
#' \code{SingleCellExperiment} with at least data in the \code{logcounts} slot.
#'
#' @param object an object of \code{SingleCellExperiment} class
#'
#' @name PrepareILoReg
#'
#' @return an object of \code{SingleCellExperiment} class
#'
#' @keywords prepare iloreg clean normalized data
#'
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom methods is
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#'
PrepareILoReg.SingleCellExperiment <- function(object) {

  # Check that there are data in `logcounts` slot
  if (!("logcounts" %in% assayNames(object))) {
    stop(paste("`Error: `logcounts` slot is missing from your ",
               "SingleCellExperiment object. This can be any kind of ",
               "normalized data matrix. Set it by executing ",
               "logcounts(object) <- norm_data",sep = ""))
    return(object)
  }

  # Remove duplicate features from the data in `logcounts` slot
  if (sum(duplicated(rownames(object))) != 0) {
    features_before <- length(rownames(object))
    object <- object[!duplicated(rownames(object)), ]
    features_after <- length(rownames(object))
    message(paste("data in SingleCellExperiment object contained duplicate ",
                  " features. ", features_before - features_after,
                  "/", features_before, " were filtered out."))
  }

  # Convert the data in `logcounts` slot into object of `dgCMatrix` class.
  if (is(logcounts(object), "matrix")) {
    logcounts(object) <- Matrix(logcounts(object),sparse = TRUE)
    message(paste("Converting object of `matrix` class into `dgCMatrix`.",
                  " Please note that ILoReg has been designed to work with ",
                  "sparse data, i.e. data with ",
                  "a high proportion of zero values! Dense data will likely " ,
                  "increase run time and memory usage drastically!",sep=""))
  }
  else if (is(logcounts(object), "data.frame")) {
    logcounts(object) <- Matrix(as.matrix(logcounts(object)),sparse = TRUE)
    message(paste("Converting object of `data.frame` class into `dgCMatrix`.",
                  " Please note that ILoReg has been designed to work with ",
                  "sparse data, i.e. data with ",
                  "a high proportion of zero values!",sep = ""))
  }
  else if (is(logcounts(object), "dgCMatrix")) {
    message("Data in `logcounts` slot already of `dgCMatrix` class...")
  }
  else {
    stop("Error: Data in `logcounts` slot is not of `matrix`, `data.frame` ",
         "or `dgCMatrix` class.")
    return(object)
  }

  # Filter genes that are not expressed in any of the cells
  genes_before_filtering <- nrow(object)
  non_expressing_genes <- rownames(object)[rowSums(logcounts(object)) != 0]
  object <- object[non_expressing_genes,]
  genes_after_filtering <- nrow(object)
  message(paste(genes_after_filtering,"/",genes_before_filtering,
                " genes remain after filtering genes with only zero values.",
                sep = ""))

  # Create a place into `metadata`` slot for the data from ILoReg
  metadata(object)$iloreg <- list()

  return(object)
}

#' @rdname PrepareILoReg
#' @aliases PrepareILoReg
setMethod("PrepareILoReg", signature(object = "SingleCellExperiment"),
          PrepareILoReg.SingleCellExperiment)

#' @title Run ICP runs parallerly
#'
#' @description
#' This functions runs in parallel \code{L} ICP runs, which is the computational
#' bottleneck of ILoReg. With ~ 3,000 cells this step should be completed
#' in ~ 2 h and ~ 1 h with 3 and 12 logical processors (threads), respectively.
#'
#' @param object An object of \code{SingleCellExperiment} class.
#' @param k A positive integer greater or equal to \code{2}, denoting
#' the number of clusters in Iterative Clustering Projection (ICP).
#' Decreasing \code{k} leads to smaller cell populations diversity
#' and vice versa. Default is \code{15}.
#' @param d A numeric greater than \code{0} and smaller than \code{1} that
#' determines how many cells \code{n} are down- or oversampled from each cluster
#' into the training data (\code{n=N/k*d}), where \code{N} is the total number
#' of cells, \code{k} is the number of clusters in ICP. Increasing above 0.3
#' leads greadually to smaller cell populations diversity.
#' Default is \code{0.3}.
#' @param L A positive integer greater than \code{1} denoting the number of
#' the ICP runs to run. Default is \code{200}. Increasing recommended with
#' a significantly larger sample size (tens of thousands of cells).
#' Default is \code{200}.
#' @param r A positive integer that denotes the number of reiterations
#' performed until the ICP algorithm stops.
#' Increasing recommended with a significantly larger sample size
#' (tens of thousands of cells). Default is \code{5}.
#' @param C A positive real number denoting the cost of constraints violation in
#' the L1-regularized logistic regression model from the LIBLINEAR library.
#' Decreasing leads to more stringent feature selection, i.e. less genes are
#' selected that are used to build the projection classifier. Decreasing to a
#' very low value (~ \code{0.01}) can lead to failure to identify central cell
#' populations. Default \code{0.3}.
#' @param reg.type "L1" or "L2". L2-regularization was not
#' investigated in the manuscript, but it leads to a more conventional
#' outcome (less subpopulations). Default is "L1".
#' @param max.iter A positive integer that denotes
#' the maximum number of iterations performed until ICP stops. This parameter
#' is only useful in situations where ICP converges extremely slowly, preventing
#' the algorithm to run too long. In most cases, reaching
#' the number of reiterations (\code{r=5}) terminates the algorithm.
#' Default is \code{200}.
#' @param threads A positive integer that specifies how many logical processors
#' (threads) to use in parallel computation.
#' Set \code{1} to disable parallelism altogether or \code{0} to use all
#' available threas except one. Default is \code{0}.
#' @param icp.batch.size A positive integer that specifies how many cells 
#' to randomly select for each ICP run from the complete data set. 
#' This is a new feature intended to speed up the process
#' with larger data sets. Default is \code{Inf}, which means using all cells.
#'
#' @name RunParallelICP
#'
#' @return an object of \code{SingleCellExperiment} class
#'
#' @keywords iterative clustering projection ICP logistic regression LIBLINEAR
#'
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doRNG %dorng%
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom doSNOW registerDoSNOW
#' @import Matrix
#' @import aricode
#' @import LiblineaR
#' @import SparseM
#' @importFrom SingleCellExperiment logcounts
#' @importFrom methods is
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,r=1,k=5)
#'
RunParallelICP.SingleCellExperiment <- function(object, k, d, L, r, C,
                                                reg.type, max.iter,
                                                threads,icp.batch.size){

  if (!is(object,"SingleCellExperiment")) {
    stop("object must of 'sce' class")
    return(object)
  }

  if (!is.numeric(k) | k < 2 | k%%1 != 0)
  {
    stop("k must be a positive integer and greater than 1")
  } else {
    metadata(object)$iloreg$k <- k
  }

  if (!is.numeric(d) | d >= 1 | d <= 0)
  {
    stop("d must be a numeric and in the range of (0,1)")
  } else {
    metadata(object)$iloreg$d <- d
  }

  if (!is.numeric(L) | L <= 0 | L%%1!=0)
  {
    stop("L must be a positive integer and greater than 0")
  } else {
    metadata(object)$iloreg$L <- L
  }

  if (!is.numeric(r) | r <= 0 | r%%1!=0)
  {
    stop("r must be a positive integer and greater than 0")
  } else {
    metadata(object)$iloreg$r <- r
  }

  if (!is.numeric(C) | C <= 0)
  {
    stop("C must be a numeric and greater than 0")
  } else {
    metadata(object)$iloreg$C <- C
  }

  if (!is.character(reg.type) | (reg.type != "L1" & reg.type != "L2"))
  {
    stop("reg.type parameter must be either 'L1' or 'L2'")
  } else {
    metadata(object)$iloreg$reg.type <- reg.type
  }

  if (!is.numeric(max.iter) | max.iter <= 0 | max.iter%%1 != 0)
  {
    stop("max.iter must be a positive integer and greater than 0")
  } else {
    metadata(object)$iloreg$max.iter <- max.iter
  }

  if (!is.numeric(threads) | threads < 0 | threads%%1 != 0)
  {
    stop("threads must be a positive integer or 0 (0 = use all available - 1)")
  } else {
    metadata(object)$iloreg$threads <- threads
  }
  
  

  if (!is.infinite(icp.batch.size))
  {
    if (!is.numeric(icp.batch.size) | icp.batch.size <= 2 | icp.batch.size%%1 != 0)
    {
      stop("icp.batch.size must be a positive integer > 2 or Inf (0 = use all cells in ICP)")
    } else {
      metadata(object)$iloreg$icp.batch.size <- icp.batch.size
    }
    
  }
  
  

  parallelism <- TRUE

  if (threads == 0) {
    cl <- makeCluster(detectCores(logical=TRUE)-1)
    # registerDoParallel(cl)
    registerDoSNOW(cl)
  } else if(threads == 1) {
    message("Parallelism disabled, because threads = 1")
    parallelism <- FALSE
  } else {
    cl<-makeCluster(threads)
    # registerDoParallel(cl)
    registerDoSNOW(cl)
  }

  dataset <- logcounts(object)

  if (parallelism) {
    pb <- txtProgressBar(min = 1, max = L, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    out <- foreach(task = seq_len(L),
                   .verbose = FALSE,
                   .combine = list,
                   .maxcombine = 1000,
                   .inorder = FALSE,
                   .multicombine = TRUE,
                   .options.snow = opts)  %dorng% {
                     try({
                       RunICP(normalized.data = dataset, k = k, d = d, r = r,
                              C = C, reg.type = reg.type, max.iter = max.iter,
                              icp.batch.size=icp.batch.size)
                     })
                   }
    close(pb)
    # stop local cluster
    stopCluster(cl)

  } else {
    out <- list()
    for (l in seq_len(L)) {
      try({
        message(paste0("ICP run: ",l))
        res <- RunICP(normalized.data = dataset, k = k, d = d, r = r,
                      C = C, reg.type = reg.type, max.iter = max.iter,
                      icp.batch.size=icp.batch.size)
        out[[l]] <- res
      })
    }
  }
  metadata(object)$iloreg$joint.probability <-
    lapply(out,function(x) x$probabilities)
  
  sds <- unlist(lapply(metadata(object)$iloreg$joint.probability,sd))
  
  metadata(object)$iloreg$joint.probability <-
    metadata(object)$iloreg$joint.probability[order(sds)]
  
  metadata(object)$iloreg$metrics <-
    lapply(out,function(x) x$metrics)

  return(object)
}
#' @rdname RunParallelICP
#' @aliases RunParallelICP
setMethod("RunParallelICP", signature(object = "SingleCellExperiment"),
          RunParallelICP.SingleCellExperiment)


#' @title PCA transformation of the joint probability matrix
#'
#' @description
#' Perform the PCA transformation of the joint probability matrix,
#' which reduces the dimensionality from k*L to p
#'
#' @param object object of \code{SingleCellExperiment} class
#' @param p a positive integer denoting the number of principal
#' components to calculate and select. Default is \code{50}.
#' @param scale a logical specifying whether the probabilities should be
#' standardized to unit-variance before running PCA. Default is \code{FALSE}.
#' @param threshold a thresfold for filtering out ICP runs before PCA with
#' the lower terminal projection accuracy below the threshold.
#' Default is \code{0}.
#'
#' @name RunPCA
#'
#' @return object of \code{SingleCellExperiment} class
#'
#' @keywords PCA eigendecomposition
#'
#' @importFrom RSpectra eigs_sym
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom S4Vectors metadata metadata<-
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#'
#'
RunPCA.SingleCellExperiment <- function(object, p, scale, threshold) {

  if (p > metadata(object)$iloreg$L*metadata(object)$iloreg$k) {
    stop(paste0("p larger than number of joint probabilities. Decrease p"))
  }

  metadata(object)$iloreg$p <- p
  metadata(object)$iloreg$scale.pca <- scale

  if (threshold == 0)
  {
    X <- do.call(cbind,metadata(object)$iloreg$joint.probability)
  } else {
    icp_runs_logical <- unlist(lapply(metadata(object)$iloreg$metrics,
                                      function(x) x["ARI",])) >= threshold
    X <- do.call(cbind,
                 metadata(object)$iloreg$joint.probability[icp_runs_logical])
  }
  X <- scale(X, scale = scale, center = TRUE)

  # X^T %*% X
  A = crossprod(X)

  # Perform eigendecomposition
  eigs_sym_out <- eigs_sym(A, p, which = "LM")

  rotated <- X %*% eigs_sym_out$vectors
  colnames(rotated) <- paste0("PC", seq_len(ncol(rotated)))

  reducedDim(object,type = "PCA") <- rotated

  return(object)

}

#' @rdname RunPCA
#' @aliases RunPCA
setMethod("RunPCA", signature(object = "SingleCellExperiment"),
          RunPCA.SingleCellExperiment)



#' @title Elbow plot of the standard deviations of the principal components
#'
#' @description
#' Draw an elbow plot of the standard deviations of the principal components
#' to deduce an appropriate value for p.
#'
#' @param object object of class 'iloreg'
#' @param return.plot logical indicating if the ggplot2 object
#' should be returned (default FALSE)
#'
#' @name PCAElbowPlot
#'
#' @return ggplot2 object if return.plot=TRUE
#'
#' @keywords PCA elbow plot
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom stats sd
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' PCAElbowPlot(sce)
#'
PCAElbowPlot.SingleCellExperiment <- function(object, return.plot) {

  df <- matrix(apply(reducedDim(object,"PCA"),2,sd),
               nrow = metadata(object)$iloreg$p,
               ncol = 1,
               dimnames =
                 list(seq_len(metadata(object)$iloreg$p),"SD"))
  df <- melt(df)

  p <- ggplot(df, aes_string(x = 'Var1', y = 'value')) +
    geom_line(color = "blue") +
    geom_point(color = "black") +
    theme_bw() +
    ylab("Standard Deviation") +
    xlab("PC")

  if (return.plot) {
    return(p)
  } else {
    print(p)
  }
}

#' @rdname PCAElbowPlot
#' @aliases PCAElbowPlot
setMethod("PCAElbowPlot", signature(object = "SingleCellExperiment"),
          PCAElbowPlot.SingleCellExperiment)


#' @title Uniform Manifold Approximation and Projection (UMAP)
#'
#' @description
#' Run nonlinear dimensionality reduction using UMAP with
#' the PCA-transformed consensus matrix as input.
#'
#' @param object of \code{SingleCellExperiment} class
#'
#' @name RunUMAP
#'
#' @return object of \code{SingleCellExperiment} class
#'
#' @keywords Uniform Manifold Approximation and Projection UMAP
#'
#' @importFrom umap umap
#' @importFrom SingleCellExperiment reducedDim reducedDim<-
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' sce <- RunUMAP(sce)
#'
RunUMAP.SingleCellExperiment <- function(object) {

  umap_out <- umap(reducedDim(object,"PCA"))

  reducedDim(object,"UMAP") <- umap_out$layout

  return(object)
}

#' @rdname RunUMAP
#' @aliases RunUMAP
setMethod("RunUMAP", signature(object = "SingleCellExperiment"),
          RunUMAP.SingleCellExperiment)


#' @title Barnes-Hut implementation of t-Distributed Stochastic
#' Neighbor Embedding (t-SNE)
#'
#' @description
#' Run nonlinear dimensionality reduction using t-SNE with the
#' PCA-transformed consensus matrix as input.
#'
#' @param object of \code{SingleCellExperiment} class
#' @param perplexity perplexity of t-SNE
#'
#' @name RunTSNE
#'
#' @return object of \code{SingleCellExperiment} class
#'
#' @keywords Barnes-Hut implementation of t-Distributed
#' Stochastic Neighbor Embedding t-SNE
#'
#' @importFrom Rtsne Rtsne
#' @importFrom SingleCellExperiment reducedDim reducedDim<-
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' sce <- RunTSNE(sce)
#'
RunTSNE.SingleCellExperiment <- function(object, perplexity) {

  rtsne_out <- Rtsne(reducedDim(object,"PCA"),
                     is_distance=FALSE,
                     perplexity=perplexity,
                     pca=FALSE)

  reducedDim(object,"TSNE") <- rtsne_out$Y

  return(object)
}

#' @rdname RunTSNE
#' @aliases RunTSNE
setMethod("RunTSNE", signature(object = "SingleCellExperiment"),
          RunTSNE.SingleCellExperiment)



#' @title Hierarchical clustering using the Ward's method
#'
#' @description
#' Perform Hierarchical clustering using the Ward's method.
#'
#' @param object of \code{SingleCellExperiment} class
#'
#' @name HierarchicalClustering
#'
#' @return object of \code{SingleCellExperiment} class
#'
#' @keywords ward hierarchical clustering
#'
#' @importFrom fastcluster hclust.vector
#' @importFrom S4Vectors metadata<-
#' @importFrom SingleCellExperiment reducedDim
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' sce <- HierarchicalClustering(sce)
#'
HierarchicalClustering.SingleCellExperiment <- function(object) {

  hc <- hclust.vector(reducedDim(object,"PCA"), method = "ward")

  metadata(object)$iloreg$hc <- hc

  return(object)
}

#' @rdname HierarchicalClustering
#' @aliases HierarchicalClustering
setMethod("HierarchicalClustering", signature(object = "SingleCellExperiment"),
          HierarchicalClustering.SingleCellExperiment)


#' @title Estimating optimal K using silhouette
#'
#' @description
#' The function estimates the optimal number of clusters K from the dendrogram
#' of the hierarhical clustering using the silhouette method.
#'
#' @param object of \code{SingleCellExperiment} class
#' @param K.start a numeric for the smallest
#' K value to be tested. Default is \code{2}.
#' @param K.end a numeric for the largest
#' K value to be tested. Default is \code{50}.
#'
#' @name CalcSilhInfo
#'
#' @return object of \code{SingleCellExperiment} class
#'
#' @keywords ward hierarchical clustering
#'
#' @importFrom S4Vectors metadata<- metadata
#' @importFrom parallelDist parDist
#' @importFrom cluster silhouette
#' @importFrom dendextend cutree
#' @importFrom stats as.dendrogram
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' sce <- HierarchicalClustering(sce)
#' sce <- CalcSilhInfo(sce)
#'
CalcSilhInfo.SingleCellExperiment <-
  function(object, K.start, K.end) {

    distance_matrix <- parDist(reducedDim(object,"PCA"),
                               method = "euclidean", threads = 1)
    distance_matrix <- as.matrix(distance_matrix)
    sis <- c()
    for (k in seq(K.start,K.end,1))
    {
      clustering <- cutree(metadata(object)$iloreg$hc,k=k)

      si <- silhouette(clustering,dmatrix = distance_matrix)
      avgsi <- summary(si)$avg.width
      sis <- c(sis,avgsi)
    }
    # Select optimal K and cluster the data
    k_optimal <- which.max(sis)+1

    message(paste0("optimal K: ",
                   k_optimal,
                   ", average silhouette score: ",
                   sis[which.max(sis)]))

    clustering <- factor(cutree(as.dendrogram(metadata(object)$iloreg$hc),
                                k = k_optimal))
    names(clustering) <- colnames(object)

    metadata(object)$iloreg$clustering.optimal <- clustering
    metadata(object)$iloreg$K.optimal <- k_optimal

    names(sis) <- seq(K.start,K.end,1)
    metadata(object)$iloreg$silhouette.information <- sis

    return(object)
  }

#' @rdname CalcSilhInfo
#' @aliases CalcSilhInfo
setMethod("CalcSilhInfo", signature(object = "SingleCellExperiment"),
          CalcSilhInfo.SingleCellExperiment)

#' @title Silhouette curve
#'
#' @description
#' Draw the silhouette curve: the average silhouette value across
#' the cells for a range of different K values.
#'
#' @param object of \code{SingleCellExperiment} class
#' @param return.plot a logical denoting whether the ggplot2 object
#' should be returned. Default is \code{FALSE}.
#'
#' @name SilhouetteCurve
#'
#' @return ggplot2 object if return.plot=TRUE
#'
#' @keywords ward hierarchical clustering
#'
#' @importFrom S4Vectors metadata
#' @import ggplot2
#' @importFrom DescTools AUC
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' sce <- HierarchicalClustering(sce)
#' sce <- CalcSilhInfo(sce)
#' SilhouetteCurve(sce)
#'
SilhouetteCurve.SingleCellExperiment <- function(object, return.plot) {

  sis <- metadata(object)$iloreg$silhouette.information
  df <- data.frame(cbind(names(sis),sis),
                   stringsAsFactors = FALSE)
  colnames(df) <- c("K","AvgSilhouette")
  df$AvgSilhouette <- as.numeric(df$AvgSilhouette)
  df$K <- as.numeric(df$K)

  auc <- round(AUC(df$K,df$AvgSilhouette),3)

  p<-ggplot(df, aes_string(x='K', y='AvgSilhouette')) +
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
}

#' @rdname SilhouetteCurve
#' @aliases SilhouetteCurve
setMethod("SilhouetteCurve", signature(object = "SingleCellExperiment"),
          SilhouetteCurve.SingleCellExperiment)

#' @title Selecting K clusters from hierarchical clustering
#'
#' @description
#' Selects K clusters from the dendrogram.
#'
#' @param object of \code{SingleCellExperiment} class
#' @param K a positive integer denoting how many clusters to select
#'
#' @name SelectKClusters
#'
#' @return object of \code{SingleCellExperiment} class
#'
#' @keywords select clusters
#'
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom dendextend cutree
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' sce <- HierarchicalClustering(sce)
#' sce <- SelectKClusters(sce,K=5)
#'
SelectKClusters.SingleCellExperiment <- function(object, K) {

  clustering <- factor(cutree(as.dendrogram(metadata(object)$iloreg$hc),k=K))
  names(clustering) <- colnames(object)

  metadata(object)$iloreg$clustering.manual <- clustering
  metadata(object)$iloreg$K.manual <- K

  return(object)

}

#' @rdname SelectKClusters
#' @aliases SelectKClusters
setMethod("SelectKClusters", signature(object = "SingleCellExperiment"),
          SelectKClusters.SingleCellExperiment)

#' @title Merge clusters
#'
#' @description
#' MergeClusters function enables merging clusters and naming the newly
#' formed cluster.
#'
#' @param object of \code{SingleCellExperiment} class
#' @param clusters.to.merge a character or numeric vector for the names of
#' the clusters to merge
#' @param new.name a character for the new name of the merged cluster.
#' If left empty, the new cluster name is formed by separating
#' the cluster names by "_".
#'
#' @name MergeClusters
#'
#' @return object of \code{SingleCellExperiment} class
#'
#' @keywords merge clusters
#'
#' @importFrom S4Vectors metadata metadata<-
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' sce <- HierarchicalClustering(sce)
#' sce <- SelectKClusters(sce,K=5)
#' sce <- MergeClusters(sce,clusters.to.merge=c(1,2),new.name="merged1")
#'
MergeClusters.SingleCellExperiment <- function(object,
                                               clusters.to.merge,
                                               new.name) {

  clusters.to.merge <- as.character(clusters.to.merge)

  clustering_old <- metadata(object)$iloreg$clustering.manual
  clusters_old <- levels(clustering_old)

  if (sum(clusters.to.merge %in% clusters_old)!=length(clusters.to.merge))
  {
    stop("invalid `clusters.to.merge argument`")
    return(object)
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

  metadata(object)$iloreg$clustering.manual <- clustering_new
  metadata(object)$iloreg$K.manual <- length(levels(clustering_new))

  return(object)

}

#' @rdname MergeClusters
#' @aliases MergeClusters
setMethod("MergeClusters", signature(object = "SingleCellExperiment"),
          MergeClusters.SingleCellExperiment)

#' @title Renaming all clusters at once
#'
#' @description
#' RenameAllClusters function enables renaming all cluster at once.
#'
#' @param object of \code{SingleCellExperiment} class
#' @param new.cluster.names object of class 'iloreg'
#'
#' @name RenameAllClusters
#'
#' @return object of \code{SingleCellExperiment} class
#'
#' @keywords rename all clusters
#'
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom plyr mapvalues
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' sce <- HierarchicalClustering(sce)
#' sce <- SelectKClusters(sce,K=5)
#' sce <- RenameAllClusters(sce,new.cluster.names=LETTERS[seq_len(5)])
#'
RenameAllClusters.SingleCellExperiment <- function(object, new.cluster.names) {

  new.cluster.names <- as.character(new.cluster.names)

  clustering_old <- metadata(object)$iloreg$clustering.manual
  clusters_old <- levels(clustering_old)

  if (length(clusters_old) != length(new.cluster.names))
  {
    stop(paste0("number of elements in clusters.to.merge is ",
                "unqual to the current number of clusters"))
    return(object)
  }

  clustering_new <- mapvalues(clustering_old,clusters_old,new.cluster.names)

  metadata(object)$iloreg$clustering.manual <- clustering_new

  return(object)

}

#' @rdname RenameAllClusters
#' @aliases RenameAllClusters
setMethod("RenameAllClusters", signature(object = "SingleCellExperiment"),
          RenameAllClusters.SingleCellExperiment)

#' @title Renaming one cluster
#'
#' @description
#' RenameCluster function enables renaming
#' a cluster in `clustering.manual` slot.
#'
#' @param object of \code{SingleCellExperiment} class
#' @param old.cluster.name a character variable denoting  the
#' old name of the cluster
#' @param new.cluster.name a character variable the
#' new name of the cluster
#'
#' @name RenameCluster
#'
#' @return object of \code{SingleCellExperiment} class
#'
#' @keywords rename one cluster
#'
#' @importFrom S4Vectors metadata metadata<-
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' sce <- HierarchicalClustering(sce)
#' sce <- SelectKClusters(sce,K=5)
#' sce <- RenameCluster(sce,1,"cluster1")
#'
RenameCluster.SingleCellExperiment <- function(object,
                                               old.cluster.name,
                                               new.cluster.name) {

  old.cluster.name <- as.character(old.cluster.name)
  new.cluster.name <- as.character(new.cluster.name)

  if (old.cluster.name=="" | new.cluster.name=="")
  {
    stop("'old.cluster.name' or 'new.cluster.name' empty\n")
  }

  clustering_old <- metadata(object)$iloreg$clustering.manual
  clusters_old <- levels(clustering_old)
  clustering_old <- as.character(clustering_old)
  names(clustering_old) <- colnames(object)

  if (!(old.cluster.name %in% clusters_old))
  {
    stop("'old.cluster.name' unvalid cluster name\n")
  }

  clustering_new <- clustering_old
  clustering_new[clustering_new==old.cluster.name] <- new.cluster.name

  clustering_new <- factor(clustering_new)
  names(clustering_new) <- names(clustering_old)

  metadata(object)$iloreg$clustering.manual <- clustering_new

  return(object)

}

#' @rdname RenameCluster
#' @aliases RenameCluster
setMethod("RenameCluster", signature(object = "SingleCellExperiment"),
          RenameCluster.SingleCellExperiment)

#' @title Visualize gene expression over nonlinear dimensionality reduction
#'
#' @description
#' GeneScatterPlot enables visualizing gene expression of a gene over
#' nonlinear dimensionality reduction with t-SNE or UMAP.
#'
#' @param object of \code{SingleCellExperiment} class
#' @param genes a character vector of the genes to be visualized
#' @param return.plot whether to return the ggplot2 object or just
#' draw it (default \code{FALSE})
#' @param dim.reduction.type "tsne" or "umap" (default "tsne")
#' @param point.size point size (default 0.7)
#' @param title text to write above the plot
#' @param plot.expressing.cells.last whether to plot the expressing genes
#' last to make the points more visible
#' @param nrow a positive integer that specifies the number of rows in
#' the plot grid. Default is \code{NULL}.
#' @param ncol a positive integer that specifies the number of columns
#' in the plot grid. Default is \code{NULL}.
#'
#' @name GeneScatterPlot
#'
#' @return ggplot2 object if return.plot=TRUE
#'
#' @keywords gene scatter plot visualization
#'
#' @importFrom SingleCellExperiment reducedDim logcounts
#' @importFrom S4Vectors metadata metadata<-
#' @import ggplot2
#' @importFrom scales muted
#' @importFrom cowplot plot_grid
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' sce <- RunTSNE(sce)
#' GeneScatterPlot(sce,"CD14",dim.reduction.type="tsne")
#' sce <- RunUMAP(sce)
#' GeneScatterPlot(sce,"CD14",dim.reduction.type="umap")
#'
GeneScatterPlot.SingleCellExperiment <- function(object,
                                                 genes,
                                                 return.plot,
                                                 dim.reduction.type,
                                                 point.size,
                                                 title,
                                                 plot.expressing.cells.last,
                                                 nrow,
                                                 ncol) {

  if (dim.reduction.type=="umap")
  {
    two.dim.data <- reducedDim(object,"UMAP")
    xlab <- "UMAP_1"
    ylab <- "UMAP_2"
  } else if (dim.reduction.type=="tsne"){
    two.dim.data <- reducedDim(object,"TSNE")
    xlab <- "tSNE_1"
    ylab <- "tSNE_2"
  } else {
    stop("dim.reduction.type must be either 'tsne' or 'umap'")
  }

  if (length(genes)==1)
  {

    df <- as.data.frame(two.dim.data)

    if (!(genes %in% rownames(object)))
    {
      stop("invalid gene name")
    }

    color.by <- logcounts(object)[genes,]
    df$group <- color.by
    colnames(df) <- c("dim1","dim2","group")

    if (title=="")
    {

      if (plot.expressing.cells.last)
      {
        df <- df[order(df$group,decreasing = FALSE),]
      }
      p<-ggplot(df, aes_string(x='dim1', y='dim2')) +
        geom_point(size=point.size,aes_string(color='group')) +
        scale_colour_gradient2(low = muted("red"), mid = "lightgrey",
                               high = "blue",name = genes) +
        xlab(xlab) +
        ylab(ylab) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"))


    } else {
      p<-ggplot(df, aes_string(x='dim1', y='dim2')) +
        geom_point(size=point.size,aes_string(color='group')) +
        scale_colour_gradient2(low = muted("red"), mid = "lightgrey",
                               high = "blue",name = genes) +
        xlab(xlab) +
        ylab(ylab) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black")) +
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

      if (!(gene %in% rownames(object)))
      {
        stop(paste0("invalid gene name: ",gene))
      }

      color.by <- logcounts(object)[gene,]
      df$group <- color.by
      colnames(df) <- c("dim1","dim2","group")

      if (plot.expressing.cells.last)
      {
        df <- df[order(df$group,decreasing = FALSE),]
      }

      if (title=="") {
        p<-ggplot(df, aes_string(x='dim1', y='dim2')) +
          geom_point(size=point.size,aes_string(color='group')) +
          scale_colour_gradient2(low = muted("red"), mid = "lightgrey",
                                 high = "blue",name = gene) +
          xlab(xlab) +
          ylab(ylab) +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"))
      } else {
        p<-ggplot(df, aes_string(x='dim1', y='dim2')) +
          geom_point(size=point.size,aes_string(color='group')) +
          scale_colour_gradient2(low = muted("red"), mid = "lightgrey",
                                 high = "blue",name = gene) +
          xlab(xlab) +
          ylab(ylab) +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black")) +
          ggtitle(title) +
          theme(plot.title = element_text(hjust = 0.5))

      }

      plot_list[[gene]] <- p

    }

    p <- plot_grid(plotlist = plot_list,align = "hv",nrow = nrow, ncol = ncol)

    if (return.plot) {
      return(p)
    } else {
      print(p)
    }
  }

}

#' @rdname GeneScatterPlot
#' @aliases GeneScatterPlot
setMethod("GeneScatterPlot", signature(object = "SingleCellExperiment"),
          GeneScatterPlot.SingleCellExperiment)

#' @title Visualize the clustering over nonliner dimensionality reduction
#'
#' @description
#' ClusteringScatterPlot function enables visualizing the clustering over
#' nonliner dimensionality reduction (t-SNE or UMAP).
#'
#' @param object of \code{SingleCellExperiment} class
#' @param clustering.type "manual" or "optimal". "manual" refers to the
#' clustering formed using the "SelectKClusters" function and "optimal" to
#' the clustering formed using the "CalcSilhInfo" function.
#' Default is "manual".
#' @param return.plot a logical denoting whether to return the ggplot2 object.
#' Default is \code{FALSE}.
#' @param dim.reduction.type "tsne" or "umap". Default is "tsne".
#' @param point.size point size. Default is Default is \code{0.7}.
#' @param title text to write above the plot
#' @param show.legend whether to show the legend on the right side of the plot.
#' Default is \code{TRUE}.
#'
#' @name ClusteringScatterPlot
#'
#' @return ggplot2 object if return.plot=TRUE
#'
#' @keywords clustering scatter plot nonlinear dimensionality reduction
#'
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom S4Vectors metadata metadata<-
#' @import ggplot2
#' @importFrom stats median
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' sce <- HierarchicalClustering(sce)
#' sce <- SelectKClusters(sce,K=5)
#' sce <- RunTSNE(sce)
#' ClusteringScatterPlot(sce,"manual",dim.reduction.type="tsne")
#' sce <- RunUMAP(sce)
#' ClusteringScatterPlot(sce,"manual",dim.reduction.type="umap")
#'
ClusteringScatterPlot.SingleCellExperiment <- function(object,
                                                       clustering.type,
                                                       return.plot,
                                                       dim.reduction.type,
                                                       point.size,
                                                       title,
                                                       show.legend) {

  if (dim.reduction.type=="umap")
  {
    two.dim.data <- reducedDim(object,"UMAP")
    xlab <- "UMAP_1"
    ylab <- "UMAP_2"
  } else if (dim.reduction.type=="tsne"){
    two.dim.data <- reducedDim(object,"TSNE")
    xlab <- "tSNE_1"
    ylab <- "tSNE_2"
  } else {
    stop("dim.reduction.type must be either 'tsne' or 'umap'")
  }

  if (clustering.type=="manual")
  {
    color.by <- metadata(object)$iloreg$clustering.manual
  } else if (clustering.type=="optimal")
  {
    color.by <- metadata(object)$iloreg$clustering.optimal
  } else {
    clustering <- metadata(object)$iloreg$clustering.manual
    message("clustering.type='manual'")
  }

  df <- as.data.frame(two.dim.data)

  df$cluster <- color.by
  colnames(df) <- c("dim1","dim2","cluster")

  two.dim.data_ <- two.dim.data
  rownames(two.dim.data_) <- names(color.by)
  cluster_centers <- lapply(levels(color.by),function(x) apply(two.dim.data_[names(color.by)[color.by==x],,drop=FALSE],2,median))
  cluster_centers <- do.call(rbind,cluster_centers)

  if (title == "")
  {
    p<-ggplot(df, aes_string(x='dim1', y='dim2')) +
      geom_point(size=point.size,aes_string(color='cluster')) +
      xlab(xlab) +
      ylab(ylab) +
      theme_classic() +
      annotate("text", x = cluster_centers[,1], y = cluster_centers[,2],
               label = levels(color.by))

  } else {

    p<-ggplot(df, aes_string(x='dim1', y='dim2')) +
      geom_point(size=point.size,aes_string(color='cluster')) +
      xlab(xlab) +
      ylab(ylab) +
      theme_classic() +
      annotate("text", x = cluster_centers[,1], y = cluster_centers[,2],
               label = levels(color.by)) +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))

  }

  if (!show.legend)
  {
    p <- p + theme(legend.position = "none")
  }

  if (return.plot) {
    return(p)
  } else {
    print(p)
  }

}

#' @rdname ClusteringScatterPlot
#' @aliases ClusteringScatterPlot
setMethod("ClusteringScatterPlot", signature(object = "SingleCellExperiment"),
          ClusteringScatterPlot.SingleCellExperiment)

#' @title identification of gene markers for all clusters
#'
#' @description
#' FindAllGeneMarkers enables identifying gene markers for all clusters at once.
#' This is done by differential expresission analysis where cells from one
#' cluster are compared against the cells from the rest of the clusters.
#' Gene and cell filters can be applied to accelerate
#' the analysis, but this might lead to missing weak signals.
#'
#' @param object of \code{SingleCellExperiment} class
#' @param clustering.type "manual" or "optimal". "manual" refers to the
#' clustering formed using the "SelectKClusters" function and "optimal"
#' to the clustering formed using the "CalcSilhInfo" function.
#' Default is "manual".
#' @param test Which test to use. Only "wilcoxon" (the Wilcoxon rank-sum test,
#' AKA Mann-Whitney U test) is supported at the moment.
#' @param log2fc.threshold Filters out genes that have log2 fold-change of the
#' averaged gene expression values below this threshold.
#' Default is \code{0.25}.
#' @param min.pct Filters out genes that have dropout rate (fraction of cells
#' expressing a gene) below this threshold in both comparison groups
#' Default is \code{0.1}.
#' @param min.diff.pct Filters out genes that do not have this minimum
#' difference in the dropout rates (fraction of cells expressing a gene)
#' between the two comparison groups. Default is \code{NULL}.
#' @param min.cells.group The minimum number of cells in the two comparison
#' groups to perform the DE analysis. If the number of cells is below the
#' threshold, then the DE analysis of this cluster is skipped.
#' Default is \code{3}.
#' @param max.cells.per.cluster The maximun number of cells per cluster if
#' downsampling is performed to speed up the DE analysis.
#' Default is \code{NULL}, i.e. no downsampling.
#' @param return.thresh If only.pos=TRUE, then return only genes that have the
#' adjusted p-value (adjusted by the Bonferroni method) below or equal to this
#' threshold. Default is \code{0.01}.
#' @param only.pos Whether to return only genes that have an adjusted p-value
#' (adjusted by the Bonferroni method) below or equal to the threshold.
#' Default is \code{FALSE}.
#'
#' @name FindAllGeneMarkers
#'
#' @return a data frame of the results if positive results were found, else NULL
#'
#' @keywords differential expression DE analysis gene markers
#'
#' @importFrom S4Vectors metadata
#' @importFrom SingleCellExperiment logcounts
#' @importFrom stats wilcox.test p.adjust
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' sce <- HierarchicalClustering(sce)
#' sce <- SelectKClusters(sce,K=5)
#' gene_markers <- FindAllGeneMarkers(sce)
#'
FindAllGeneMarkers.SingleCellExperiment <- function(object,
                                                    clustering.type,
                                                    test,
                                                    log2fc.threshold,
                                                    min.pct,
                                                    min.diff.pct,
                                                    min.cells.group,
                                                    max.cells.per.cluster,
                                                    return.thresh,
                                                    only.pos) {

  number.of.expressed.genes <- nrow(object)

  if (clustering.type=="manual")
  {
    clustering <- metadata(object)$iloreg$clustering.manual
  } else if (clustering.type=="optimal")
  {
    clustering <- metadata(object)$iloreg$clustering.optimal
  } else {
    clustering <- metadata(object)$iloreg$clustering.manual
    cat("clustering.type='manual'")
  }

  data <- logcounts(object)

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
        inds <- sample(seq_len(cells_in_cluster),
                       size = max.cells.per.cluster,
                       replace = FALSE)
        names_cluster <- names(clustering[clustering==cluster])
        cells_downsampled <- c(cells_downsampled,names_cluster[inds])
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

    # Skip if the number of cells in the test
    # or the reference set is lower than min.cells.group
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

    log2FC <- cluster_aves - other_aves
    
    genes_to_include <- rownames(data_cluster)[log2FC >= log2fc.threshold | log2FC <= -log2fc.threshold]

    data_cluster <- data_cluster[genes_to_include,,drop=FALSE]
    data_other <- data_other[genes_to_include,,drop=FALSE]


    cat(paste0(nrow(data_cluster)," genes left after log2fc.threshold filtering\n"))
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

}

#' @rdname FindAllGeneMarkers
#' @aliases FindAllGeneMarkers
setMethod("FindAllGeneMarkers", signature(object = "SingleCellExperiment"),
          FindAllGeneMarkers.SingleCellExperiment)

#' @title Identification of gene markers for a cluster or two arbitrary
#' combinations of clusters
#'
#' @description
#' FindGeneMarkers enables identifying gene markers for one cluster or
#' two arbitrary combinations of clusters, e.g. 1_2 vs. 3_4_5.
#' Gene and cell filters can be applied to accelerate
#' the analysis, but this might lead to missing weak signals.
#'
#' @param object of \code{SingleCellExperiment} class
#' @param clusters.1 a character or numeric vector denoting which clusters
#' to use in the first group (named group.1 in the results)
#' @param clusters.2 a character or numeric vector denoting which clusters
#' to use in the second group (named group.2 in the results)
#' @param clustering.type "manual" or "optimal". "manual" refers to the
#' clustering formed using the "SelectKClusters" function and "optimal" to
#' the clustering formed using the "CalcSilhInfo" function.
#' Default is "manual".
#' @param test Which test to use. Only "wilcoxon" (the Wilcoxon rank-sum test,
#' AKA Mann-Whitney U test) is supported at the moment.
#' @param logfc.threshold Filters out genes that have log2 fold-change of the
#' averaged gene expression values below this threshold.
#' Default is \code{0.25}.
#' @param min.pct Filters out genes that have dropout rate (fraction of cells
#' expressing a gene) below this threshold in both comparison groups
#' Default is \code{0.1}.
#' @param min.diff.pct Filters out genes that do not have this minimum
#' difference in the dropout rates (fraction of cells expressing a gene)
#' between the two comparison groups. Default is \code{NULL}.
#' @param min.cells.group The minimum number of cells in the two comparison
#' groups to perform the DE analysis. If the number of cells is below the
#' threshold, then the DE analysis is not performed.
#' Default is \code{3}.
#' @param max.cells.per.cluster The maximun number of cells per cluster
#' if downsampling is performed to speed up the DE analysis.
#' Default is \code{NULL}, i.e. no downsampling.
#' @param return.thresh If only.pos=TRUE, then return only genes that
#' have the adjusted p-value (adjusted by the Bonferroni method) below or
#' equal to this threshold.  Default is \code{0.01}.
#' @param only.pos Whether to return only genes that have an adjusted
#' p-value (adjusted by the Bonferroni method) below or equal to the
#' threshold. Default is \code{FALSE}.
#'
#' @name FindGeneMarkers
#'
#' @return a data frame of the results if positive results were found, else NULL
#'
#' @keywords differential expression DE analysis gene markers
#'
#' @importFrom S4Vectors metadata
#' @importFrom stats wilcox.test p.adjust
#' @importFrom SingleCellExperiment logcounts
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' sce <- HierarchicalClustering(sce)
#' sce <- SelectKClusters(sce,K=5)
#' gene_markes_1 <- FindGeneMarkers(sce,clusters.1=1)
#' gene_markes_1_vs_2 <- FindGeneMarkers(sce,clusters.1=1,clusters.2=2)
#'
FindGeneMarkers.SingleCellExperiment <- function(object,
                                                 clusters.1,
                                                 clusters.2,
                                                 clustering.type,
                                                 test,
                                                 logfc.threshold,
                                                 min.pct,
                                                 min.diff.pct,
                                                 min.cells.group,
                                                 max.cells.per.cluster,
                                                 return.thresh,
                                                 only.pos) {

  if (clustering.type=="manual")
  {
    clustering <- metadata(object)$iloreg$clustering.manual
  } else if (clustering.type=="optimal")
  {
    clustering <- metadata(object)$iloreg$clustering.optimal
  } else {
    clustering <- metadata(object)$iloreg$clustering.manual
    cat("clustering.type='manual'")
  }

  data <- logcounts(object)

  cells_to_include_1 <- names(clustering)[clustering %in% clusters.1]
  clustering_1 <- factor(rep("group.1",length(cells_to_include_1)))
  names(clustering_1) <- cells_to_include_1

  if (is.null(clusters.2)) {
    clusters.2 <- setdiff(levels(clustering),clusters.1)
  }
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
        inds <- sample(seq_len(cells_in_cluster),
                       size = max.cells.per.cluster,
                       replace = FALSE)
        cells_downsampled <- c(cells_downsampled,
                               names(clustering[clustering==cluster])[inds])
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

  # Skip if the number of cells in the test or the reference set
  # is lower than min.cells.group
  if (ncol(data_cluster) < min.cells.group | ncol(data_other) < min.cells.group)
  {
    cat("-----------------------------------\n")
    return(NULL)
  }

  # min.pct filter
  genes.pct_cluster <-
    apply(data_cluster,1,function(x) sum(x!=0))/ncol(data_cluster)
  genes.pct_other <-
    apply(data_other,1,function(x) sum(x!=0))/ncol(data_other)

  genes_to_include <-
    rownames(data_cluster)[genes.pct_cluster>=min.pct |
                             genes.pct_other >= min.pct]


  data_cluster <- data_cluster[genes_to_include,,drop=FALSE]
  data_other <- data_other[genes_to_include,,drop=FALSE]

  cat(paste0(nrow(data_cluster)," genes left after min.pct filtering\n"))
  if (nrow(data_cluster)==0)
  {
    cat("-----------------------------------\n")
    return(NULL)
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
    return(NULL)
  }

  # logfc.threshold filter
  # Calculate log2 fold changes
  cluster_aves <- apply(data_cluster,1,mean)
  other_aves <- apply(data_other,1,mean)

  log2FC <- cluster_aves - other_aves

  genes_to_include <- rownames(data_cluster)[log2FC >= logfc.threshold | log2FC <= -logfc.threshold]

  data_cluster <- data_cluster[genes_to_include,,drop=FALSE]
  data_other <- data_other[genes_to_include,,drop=FALSE]

  cat(paste0(nrow(data_cluster)," genes left after logfc.threshold filtering\n"))
  if (nrow(data_cluster)==0)
  {
    cat("-----------------------------------\n")
    return(NULL)
  }

  # Run DE test

  if (test=="wilcox")
  {
    wilcox.res <- lapply(rownames(data_cluster),function(x) wilcox.test(x=data_cluster[x,],y=data_other[x,]))
    p_values <- unlist(lapply(wilcox.res,function(x) x$p.value))
    names(p_values) <- rownames(data_cluster)

    # Adjust p-values
    adj_p_values <- p.adjust(p_values, method = "bonferroni", n = nrow(object))

    res <- cbind(p_values,
                 adj_p_values,
                 log2FC[names(p_values)],
                 genes.pct_cluster[names(p_values)],
                 genes.pct_other[names(p_values)],
                 abs(genes.pct_cluster-genes.pct_other)[names(p_values)])
    colnames(res) <- c("p.value","adj.p.value","log2FC","pct.1","pct.2","diff.pct")
    res <- as.data.frame(res)
    res$cluster <- cluster
    res$gene <- names(p_values)
  }

  results_list[[cluster]] <- res

  results_df <- do.call(rbind,results_list)
  rownames(results_df) <- make.unique(unlist(lapply(results_list,rownames)))

  results_df$cluster <- NULL

  if(only.pos) {
    results_df <- results_df[results_df$adj.p.value <= return.thresh,]
    return(results_df)
  }
  return(results_df)

}

#' @rdname FindGeneMarkers
#' @aliases FindGeneMarkers
setMethod("FindGeneMarkers", signature(object = "SingleCellExperiment"),
          FindGeneMarkers.SingleCellExperiment)

#' @title Gene expression visualization using violin plots
#'
#' @description
#' The VlnPlot function enables visualizing expression levels of a gene,
#' or multiple genes, across clusters using Violin plots.
#'
#' @param object of \code{SingleCellExperiment} class
#' @param clustering.type "manual" or "optimal". "manual"
#' refers to the clustering formed using the "SelectKClusters" function
#' and "optimal" to the clustering formed using the
#' "CalcSilhInfo" function. Default is "manual".
#' @param genes a character vector denoting the gene names that are visualized
#' @param return.plot return.plot whether to return the ggplot2 object
#' @param rotate.x.axis.labels a logical denoting whether the x-axis
#' labels should be rotated 90 degrees.
#' or just draw it. Default is \code{FALSE}.
#'
#' @name VlnPlot
#'
#' @return ggplot2 object if return.plot=TRUE
#'
#' @keywords violin plot
#'
#' @importFrom S4Vectors metadata
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom SingleCellExperiment logcounts
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' sce <- HierarchicalClustering(sce)
#' sce <- SelectKClusters(sce,K=5)
#' VlnPlot(sce,genes=c("CD3D","CD79A","CST3"))
#'
VlnPlot.SingleCellExperiment <- function(object,
                                         clustering.type,
                                         genes,
                                         return.plot,
                                         rotate.x.axis.labels) {


  if (clustering.type=="manual")
  {
    clustering <- metadata(object)$iloreg$clustering.manual
  } else if (clustering.type=="optimal")
  {
    clustering <- metadata(object)$iloreg$clustering.optimal
  } else {
    clustering <- metadata(object)$iloreg$clustering.manual
    message("clustering.type='manual'")
  }

  data <- logcounts(object)

  df <- as.numeric(t(data[genes,]))
  df <- data.frame(matrix(df,ncol = 1,dimnames = list(seq_len(length(df)),"Expression")))
  df$gene  <- unlist(lapply(genes,function(x) rep(x,ncol(data))))
  df$gene <- factor(df$gene)
  df$Cluster <- rep(as.character(clustering),length(genes))
  df$Cluster <- factor(df$Cluster)


  if (rotate.x.axis.labels)
  {
    plotlist <- lapply(genes,function(x) ggplot(df[df$gene==x,], aes_string(x='Cluster', y='Expression', fill='Cluster'))+geom_violin(trim=TRUE)+geom_jitter(height = 0, width = 0.1)+theme_classic()+ggtitle(x)+theme(plot.title = element_text(hjust = 0.5),legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1)))
  } else {
    plotlist <- lapply(genes,function(x) ggplot(df[df$gene==x,], aes_string(x='Cluster', y='Expression', fill='Cluster'))+geom_violin(trim=TRUE)+geom_jitter(height = 0, width = 0.1)+theme_classic()+ggtitle(x)+theme(plot.title = element_text(hjust = 0.5),legend.position = "none"))

  }


  p <- plot_grid(plotlist = plotlist)

  if (return.plot)
  {
    return(p)
  } else {
    print(p)
  }

}

#' @rdname VlnPlot
#' @aliases VlnPlot
setMethod("VlnPlot", signature(object = "SingleCellExperiment"),
          VlnPlot.SingleCellExperiment)

#' @title Heatmap visualization of the gene markers identified by FindAllGeneMarkers
#'
#' @description
#' The GeneHeatmap function enables drawing a heatmap of the gene markers
#' identified by FindAllGeneMarkers, where the cell are grouped
#' by the clustering.
#'
#' @param object of \code{SingleCellExperiment} class
#' @param clustering.type "manual" or "optimal". "manual" refers to the
#' clustering formed using the "SelectKClusters" function and "optimal"
#' to the clustering using the "CalcSilhInfo" function.
#' Default is "manual".
#' @param gene.markers a data frame of the gene markers generated
#' by FindAllGeneMarkers function. To accelerate the drawing, filtering
#' the dataframe by selecting e.g. top 10 genes is recommended.
#'
#' @name GeneHeatmap
#'
#' @return nothing
#'
#' @keywords gene heatmap grouped
#'
#' @importFrom S4Vectors metadata
#' @import pheatmap
#' @importFrom SingleCellExperiment logcounts
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,r=1,k=5) # Use L=200
#' sce <- RunPCA(sce,p=5)
#' sce <- HierarchicalClustering(sce)
#' sce <- SelectKClusters(sce,K=5)
#' gene_markers <- FindAllGeneMarkers(sce,log2fc.threshold = 0.5,min.pct = 0.5)
#' top10_log2FC <- SelectTopGenes(gene_markers,top.N=10,
#' criterion.type="log2FC",inverse=FALSE)
#' GeneHeatmap(sce,clustering.type = "manual",
#'  gene.markers = top10_log2FC)
#'
GeneHeatmap.SingleCellExperiment <- function(object,
                                             clustering.type,
                                             gene.markers) {


  if (clustering.type=="manual")
  {
    clustering <- metadata(object)$iloreg$clustering.manual
  } else if (clustering.type=="optimal")
  {
    clustering <- metadata(object)$iloreg$clustering.optimal
  } else {
    clustering <- metadata(object)$iloreg$clustering.manual
    cat("clustering.type='manual'")
  }

  data <- logcounts(object)
  data <- data[unique(gene.markers$gene),]
  # data <- scale(data,center = TRUE,scale = TRUE)
  data <- data[,order(clustering)]


  # Generate column annotations
  annotation = data.frame(cluster=sort(clustering))

  pheatmap(data,show_colnames = FALSE,
           gaps_col = cumsum(table(clustering[order(clustering)])),
           gaps_row = cumsum(table(gene.markers[!duplicated(gene.markers$gene),"cluster"])),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           annotation_col = annotation)

}

#' @rdname GeneHeatmap
#' @aliases GeneHeatmap
setMethod("GeneHeatmap", signature(object = "SingleCellExperiment"),
          GeneHeatmap.SingleCellExperiment)

#' @title Visualiation of a custom annotation over nonlinear
#' dimensionality reduction
#'
#' @description
#' The AnnotationScatterPlot enables visualizing arbitrary class labels
#' over the nonliner dimensionality reduction, e.g. t-SNE or UMAP.
#'
#' @param object of \code{SingleCellExperiment} class
#' @param annotation a character vector, factor or numeric for the class labels.
#' @param return.plot return.plot whether to return the ggplot2 object or
#' just draw it. Default is \code{FALSE}.
#' @param dim.reduction.type "tsne" or "umap". Default is \code{tsne}.
#' @param point.size point size. Default is \code{0.7}.
#' @param show.legend a logical denoting whether to show the legend on the right
#' side of the plot. Default is \code{TRUE}.
#'
#' @name AnnotationScatterPlot
#'
#' @return ggplot2 object if return.plot=TRUE
#'
#' @keywords annotation custom visualization t-sne umap nonlinear
#' dimensionality reduction
#'
#' @importFrom S4Vectors metadata
#' @import ggplot2
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom stats median
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = pbmc3k_500))
#' sce <- PrepareILoReg(sce)
#' ## These settings are just to accelerate the example, use the defaults.
#' sce <- RunParallelICP(sce,L=2,threads=1,C=0.1,k=5,r=1)
#' sce <- RunPCA(sce,p=5)
#' sce <- RunTSNE(sce)
#' sce <- HierarchicalClustering(sce)
#' sce <- SelectKClusters(sce,K=5)
#' ## Change the names to the first five alphabets and Visualize the annotation.
#' custom_annotation <- plyr::mapvalues(metadata(sce)$iloreg$clustering.manual,
#'                                      c(1,2,3,4,5),
#'                                      LETTERS[1:5])
#' AnnotationScatterPlot(sce,
#'                       annotation = custom_annotation,
#'                       return.plot = FALSE,
#'                       dim.reduction.type = "tsne",
#'                       show.legend = FALSE)
#'
#'
AnnotationScatterPlot.SingleCellExperiment <- function(object,
                                                       annotation,
                                                       return.plot,
                                                       dim.reduction.type,
                                                       point.size,
                                                       show.legend) {


  if (dim.reduction.type == "umap")
  {
    two.dim.data <- reducedDim(object,"UMAP")
    xlab <- "UMAP_1"
    ylab <- "UMAP_2"
  } else if (dim.reduction.type == "tsne"){
    two.dim.data <- reducedDim(object,"TSNE")
    xlab <- "tSNE_1"
    ylab <- "tSNE_2"
  } else {
    stop("dim.reduction.type must be either 'tsne' or 'umap'")
  }

  annotation <- factor(as.character(annotation))
  names(annotation) <- colnames(object)

  df <- as.data.frame(two.dim.data)
  df$group <- annotation
  colnames(df) <- c("dim1","dim2","group")

  two.dim.data_ <- two.dim.data
  rownames(two.dim.data_) <- names(annotation)
  cluster_centers <- lapply(levels(annotation),function(x) apply(two.dim.data_[names(annotation)[annotation==x],,drop=FALSE],2, median))
  cluster_centers <- do.call(rbind,cluster_centers)


  p<-ggplot(df, aes_string(x='dim1', y='dim2')) +
    geom_point(size=point.size,aes_string(color='group')) +
    xlab(xlab) +
    ylab(ylab) +
    theme_classic() +
    annotate("text",
             x = cluster_centers[,1],
             y = cluster_centers[,2],
             label = levels(annotation)) +
    guides(colour = guide_legend(override.aes = list(size=2)))


  if (!show.legend)
  {
    p <- p + theme(legend.position = "none")
  }

  if (return.plot)
  {
    return(p)
  } else {
    print(p)
  }

}

#' @rdname AnnotationScatterPlot
#' @aliases AnnotationScatterPlot
setMethod("AnnotationScatterPlot", signature(object = "SingleCellExperiment"),
          AnnotationScatterPlot.SingleCellExperiment)

