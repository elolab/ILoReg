#' 500 cells downsampled from the pbmc3k dataset.
#'
#' The preprocessing was done using Cell Ranger v2.2.0 and
#' the GRCh37.p13 human reference genome.
#'
#' @docType data
#'
#' @usage data(pbmc3k_500)
#'
#' @format raw_data and data, both dgCMatrix objects
#'
#' @keywords datasets
#'
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression}
#'
#' @examples
#' data(pbmc3k_500)
"pbmc3k_500"



#' @title Iterative Clustering Projection (ICP) clustering
#'
#' @description
#' The function implements Iterative Clustering Projection (ICP): a
#' supervised learning -based clustering, which maximizes clustering similarity
#' between the clustering and its projection by logistic regression.
#'
#' @param normalized.data A sparse matrix (dgCMatrix) containing
#' normalized gene expression data with genes in rows and cells in columns.
#' Default is \code{NULL}.
#' @param k A positive integer greater or equal to 2, denoting the number of
#' clusters in ICP. Default is \code{15}.
#' @param d A numeric that defines how many cells per cluster should be
#' down- and oversampled (d in N/k*d), when stratified.downsampling=FALSE,
#' or what fraction should be downsampled in the stratified approach
#' ,stratified.downsampling=TRUE. Default is \code{0.3}.
#' @param r A positive integer that denotes the number of reiterations
#' performed until the algorithm stops. Default is \code{5}.
#' @param C Cost of constraints violation (\code{C}) for L1-regulatization.
#' Default is \code{0.3}.
#' @param reg.type "L1" for LASSO and "L2" for Ridge. Default is "L1".
#' @param max.iter A positive integer that denotes the maximum number of
#' iterations performed until the algorithm ends. Default is \code{200}.
#'
#' @return A list comprising the probability matrix and the clustering
#' similarity measures: ARI, NMI, etc.
#'
#' @keywords iterative clustering projection ICP clustering
#'
#' @import Matrix
#' @importFrom aricode clustComp
#' @import LiblineaR
#' @import SparseM
#'
#' @export
#'
RunICP <- function(normalized.data = NULL,k = 15, d = 0.3, r = 5, C = 5,
                   reg.type = "L1", max.iter = 200) {

  first_round <- TRUE
  metrics <- NULL
  idents <- list()
  iterations <- 1
  probs <- NULL

  while (TRUE) {

    # Step 1: initialize clustering (ident_1) randomly, ARI=0 and r=0
    if (first_round) {
      ident_1 <- factor(sample(seq_len(k),ncol(normalized.data),replace = TRUE))
      names(ident_1) <- colnames(normalized.data)
      idents[[1]] <- ident_1
      ari <- 0
      reiterations <- 0
    }

    # Step 2: train logistic regression model
    res <- LogisticRegression(training.sparse.matrix = t(normalized.data),
                              training.ident = ident_1, C = C,
                              reg.type=reg.type,
                              test.sparse.matrix = t(normalized.data), d=d)

    names(res$predictions) <- colnames(normalized.data)
    rownames(res$probabilities) <- colnames(normalized.data)

    # Projected cluster probabilities
    probs <- res$probabilities

    # Projected clusters
    ident_2 <- res$predictions

    # Safety procedure: If k drops to 1, start from the beginning.
    # k should NOT decrease during the iteration when
    # the down- and oversampling approach is used for balancing training data.
    if (length(levels(factor(as.character(ident_2)))) < 2)
    {
      first_round <- TRUE
      metrics <- NULL
      idents <- list()
      iterations <- 1
      next
    }

    # Step 3: compare clustering similarity between clustering and projection
    comp_clust <- clustComp(c1 = ident_1, c2 = ident_2)

    if(first_round & comp_clust$ARI <= 0)
    {
      next
    }

    # Step 3.1: If ARI did not increase, reiterate
    if (comp_clust$ARI <= ari & !(first_round))
    {
      reiterations <- reiterations + 1
    }
    # Step 3.2: If ARI increased, proceed to next iteration round
    else {
      # Update clustering to the predicted clusters
      ident_1 <- ident_2
      first_round <- FALSE
      metrics <- cbind(metrics,comp_clust)
      iterations = iterations + 1
      idents[[iterations]] <- ident_2
      ari <- comp_clust$ARI
      reiterations <- 0

    }
    # Step 4: If the maximum number of reiterations or iterations was
    # reached, break the while loop
    if (reiterations == r | iterations == max.iter)
    {
      break
    }
  }
  # Step 5: Return result
  return(list(probabilities=probs, metrics=metrics))
}

#' @title Down- and oversample data
#'
#' @description
#' The function implements a script down- and oversamples data to
#' include n cells.
#'
#' @param x A character or numeric vector of data to down-and oversample.
#' @param n How many cells to include per cluster.
#'
#' @return a list containing the output of the LiblineaR prediction
#'
#' @keywords downsampling oversampling
#'
#' @export
#'
DownOverSampling <- function(x, n = 50) {
  if (length(x) < n) {
    res <- sample(x, size = n, replace = TRUE)
  } else {
    res <- sample(x, size = n, replace = FALSE)
  }
  return(res)
}


#' @title Clustering projection using logistic regression from
#' the LiblineaR R package
#'
#' @description
#' The function implements a script that downsamples data a dataset, trains
#' a logistic regression classifier model
#' and then projects its clustering onto itself using a trained
#' L1-regularized logistic regression model.
#'
#' @param training.sparse.matrix A sparse matrix (dgCMatrix) containing training
#' sample's gene expression data with genes in rows and cells in columns.
#' Default is \code{NULL}.
#' @param training.ident A named factor containing sample's cluster labels for
#' each cell in training.sparse.matrix. Default is \code{NULL}.
#' @param C Cost of constraints violation in L1-regularized logistic
#' regression (C). Default is \code{0.3}.
#' @param reg.type "L1" for LASSO and "L2" for Ridge. Default is "L1".
#' @param test.sparse.matrix A sparse matrix (dgCMatrix) containing test
#' sample's gene expression data with genes in rows and cells in columns.
#' Default is \code{NULL}.
#' @param d A numeric smaller than \code{1} and greater than \code{0}
#' that determines how many cells per cluster should be
#' down- and oversampled (d in N/k*d), where N is the total number of cells
#' and k the number of clusters. Default is \code{0.3}.
#'
#' @return a list containing the output of the LiblineaR prediction
#'
#' @keywords logistic regression LiblineaR projection downsampling oversampling
#'
#' @import Matrix
#' @import SparseM
#' @importFrom methods as
#' @importFrom LiblineaR LiblineaR
#' @importFrom stats predict
#'
#' @export
#'
LogisticRegression <- function(training.sparse.matrix = NULL,
                               training.ident = NULL,
                               C = 0.3,
                               reg.type = "L1",
                               test.sparse.matrix = NULL,
                               d = 0.3) {

  # Downsample training data
  if (!is.null(d))
  {
    cells_per_cluster <- ceiling((length(training.ident) / (length(levels(training.ident)))) * d)

    training_ident_subset <- as.character(unlist(lapply(split(names(training.ident),training.ident), function(x) DownOverSampling(x,cells_per_cluster))))

    training.ident <- training.ident[training_ident_subset]
    training.sparse.matrix <- training.sparse.matrix[training_ident_subset,]
  }

  # Transform training and test data from dgCMatrix to matrix.csr
  training.sparse.matrix <- as(training.sparse.matrix,"matrix.csr")
  test.sparse.matrix <- as(test.sparse.matrix,"matrix.csr")

  if (reg.type=="L2")
  {
    type <- 7
  } else if (reg.type=="L1")
  {
    type <- 6 #L1
  } else {
    stop("'reg.type' must be either 'L1' or 'L2'")
  }

  model <- LiblineaR(training.sparse.matrix, training.ident,
                     type = type, cost = C)
  prediction <- predict(model,proba = TRUE,test.sparse.matrix)

  return(prediction)
}

#' @title Select top or bottom N genes based on a selection criterion
#'
#' @description
#' The SelectTopGenes function enables selecting top or bottom N genes based
#' on a criterion (e.g. log2FC or adj.p.value).
#'
#' @param gene.markers A data frame of the gene markers found by
#' FindAllGeneMarkers function.
#' @param top.N How many top or bottom genes to select. Default is \code{10}.
#' @param criterion.type Which criterion to use for selecting the genes.
#' Default is "log2FC".
#' @param reverse Whether to select bottom instead of top N genes.
#' Default is \code{FALSE}.
#'
#' @return an object of `data.frame` class
#'
#' @keywords select top bottom N genes
#'
#' @importFrom dplyr group_by %>% top_n
#'
#' @export
#'
SelectTopGenes <- function(gene.markers = NULL, top.N = 10,
                           criterion.type = "log2FC", reverse=FALSE)
{

  if (reverse)
  {
    top.N <- -top.N
  }

  gene.markers %>%
    group_by(.data$cluster) %>% top_n(top.N, get(criterion.type)) -> top_N

  return(top_N)
}
