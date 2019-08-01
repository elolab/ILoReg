#'  Classification using logistic regression from LiblineaR R package
#'
#' Enables training of a logistic regression model using LiblineaR package.
#' @param training.sparse.matrix A sparse matrix (dgCMatrix) containing training sample's gene expression data with genes in rows and cells in columns. (default NULL)
#' @param training.ident A named factor containing sample's cluster labels for each cell in training.sparse.matrix. (default NULL)
#' @param print.info A logical denoting if information about the run (e.g. running time) should be printed. (default TRUE)
#' @param cost Cost of constraints (C). Recommendation is in range [0.2,0.5] for LASSO.  (default 0.3)
#' @param type "L1" for LASSO and "L2" for Ridge.  (default "L1")
#' @param test.sparse.matrix A sparse matrix (dgCMatrix) containing test sample's gene expression data with genes in rows and cells in columns. (default NULL)
#' @param downsampling.fraction A numeric that defines how many cells per cluster should be down- and oversampled (d in N/k*d), when stratified.downsampling=FALSE,  or what fraction should be downsampled in the stratified approach ,stratified.downsampling=TRUE. (default 0.3)
#' @param stratified.downsampling A logical denoting if stratified downsampling approach should be used in generation of the training data: (default FALSE).
#' @return predictions for test data
#' @keywords logistic regression stratified LiblineaR
#' @import Matrix
#' @import tictoc
#' @import LiblineaR
#' @import SparseM
#' @export
#' @examples
#' a <- c(0,1,2)


LogisticRegression <- function(training.sparse.matrix = NULL, training.ident = NULL,
                               print.info=TRUE,cost = 0.3,type="L1",
                               test.sparse.matrix = NULL,downsampling.fraction=0.3)

{

  # Downsample training data
  if (!is.null(downsampling.fraction))
  {
    cells_per_cluster <- ceiling((length(training.ident)/(length(levels(training.ident))))*downsampling.fraction)
    training.ident.subset <- as.character(unlist(lapply(split(names(training.ident),training.ident),
                                                        function(x) if (length(x) < cells_per_cluster) {sample(x,size = cells_per_cluster,replace = T)} else {sample(x,size = cells_per_cluster,replace = F)})))

    training.ident <- training.ident[training.ident.subset]
    training.sparse.matrix <- training.sparse.matrix[training.ident.subset,]
  }

  # Transform training data from dgCMatrix to matrix.csr
  training_size <- nrow(training.sparse.matrix)

  colnames_matrix <- colnames(training.sparse.matrix)

  tic(msg="Converting Matrix to SparseM")
  training.sparse.matrix <- new("matrix.csc", ra = training.sparse.matrix@x,
                                ja = training.sparse.matrix@i + 1L,
                                ia = training.sparse.matrix@p + 1L,
                                dimension = training.sparse.matrix@Dim)
  training.sparse.matrix <- as.matrix.csr(training.sparse.matrix)
  toc(log=TRUE)

  # Transform test data from dgCMatrix to matrix.csr
  tic(msg="Converting Matrix to SparseM")
  test.sparse.matrix <- new("matrix.csc", ra = test.sparse.matrix@x,
                            ja = test.sparse.matrix@i + 1L,
                            ia = test.sparse.matrix@p + 1L,
                            dimension = test.sparse.matrix@Dim)
  test.sparse.matrix <- as.matrix.csr(test.sparse.matrix)
  toc(log=TRUE)

  if (type=="L2")
  {
    type <- 7
  } else if (type=="L1")
  {
    type <- 6 #LASSO
  } else {
    stop("'type' must be either 'L1' (LASSO) or 'L2' (RIDGE)")
  }


  if (print.info) {

    tic(msg="Train LIBLINEAR model")
    model <- LiblineaR(training.sparse.matrix, training.ident, type = type, cost = cost,
                       svr_eps = NULL, bias = 1, wi = NULL, cross = 0, verbose = FALSE,
                       findC = FALSE, useInitC = TRUE)
    toc(log=TRUE)


  }

  # Predict test data using the model
  tic(msg="Predict using LIBLINEAR model")
  prediction <- predict(model,proba = TRUE,test.sparse.matrix)
  toc(log=TRUE)

  # Return result
  return(prediction)

}
