#'  Iterative logistic regression (ILoReg) clustering
#'
#' Enables iterative logistic regression (ILoReg) clustering.
#' @param normalized.data A sparse matrix (dgCMatrix) containing gene expression with genes in rows and cells in columns. (default NULL)
#' @param k A positive integer greater or equal to 2, which denotes the number of clusters in iterative logistic regression. If stratified.downsampling=TRUE, the final k is likely to be smaller. (default 15)
#' @param d A numeric that defines how many cells per cluster should be down- and oversampled (d in N/k*d), when stratified.downsampling=FALSE,  or what fraction should be downsampled in the stratified approach ,stratified.downsampling=TRUE. (default 0.3)
#' @param r A positive integer that denotes the maximum number of reiterations performed until the algorithm ends. (default 5)
#' @param C Cost of constraints (C). Recommended C is in range [0.2,0.5] for LASSO.  (default 0.3)
#' @param type "L1" for LASSO and "L2" for Ridge.  (default "L1")
#' @param max.number.of.iterations A positive integer that denotes the maximum number of iterations performed until the algorithm ends. (default 200)
#' @return return the probability matrix, and clutering comparison measures if return.metrics=TRUE or UMAP layouts from each epoch if perform.umap=TRUE
#' @keywords iterative logistic regression ILoReg clustering
#' @import Matrix
#' @import tictoc
#' @import aricode
#' @import LiblineaR
#' @import SparseM
#' @export
#' @examples
#' a <- c(0,1,2)


RunICP <- function(normalized.data = NULL,k=15,d=0.3,r=5,C=5,type="L1",max.number.of.iterations=100)
{

  first_round = TRUE
  metrics <- NULL
  idents <- list()
  iterations = 1 # iteration

  probs <- c()

  # Main loop that runs until the maximum number of reiterations (default reiterations=5) or maximum number of iterations (default iterations=200) is reached
  while (TRUE) {
    # Measure run time
    tic(msg="Running iterative logistic regression")

    # Step 1: initialize clustering randomly, ARI=0 and r=0
    cat(paste0("Iteration round ",iterations,"\n"))
    if (first_round) {
      ident_1 <- factor(sample(1:k,ncol(normalized.data),replace = T))
      names(ident_1) <- colnames(normalized.data)
      idents[[1]] <- ident_1
      ari <- 0
      reiterations <- 0
    }

    # Step 2: train logistic regression model
    res <- LogisticRegression(training.sparse.matrix = t(normalized.data),
                              training.ident = ident_1,
                              print.info = TRUE,
                              cost = C,
                              type=type,
                              test.sparse.matrix = t(normalized.data),
                              downsampling.fraction=d)

    names(res$predictions) <- colnames(normalized.data)
    rownames(res$probabilities) <- colnames(normalized.data)

    # Predicted cluster probabilities
    probs <- res$probabilities

    # Predicted clusters
    ident_2 <- res$predictions

    # If the number of clusters in prediction dropped below 2, start from the beginning
    if (length(levels(factor(as.character(ident_2)))) < 2)
    {
      first_round = TRUE
      metrics <- NULL
      idents <- list()
      iterations <- 1
      next
    }

    # Step 3: compare clustering of training data and prediction using ARI
    compclust <- clustComp(c1 = ident_1,c2 = ident_2)

    if(first_round & compclust$ARI <= 0)
    {
      next
    }

    # Step 3.1: If ARI did not increase, reiterate
    if (compclust$ARI <= ari & !(first_round))
    {
      cat(paste0("ARI decreased... re-iterating\n"))
      reiterations <- reiterations + 1
    }
    # Step 3.2: If ARI increased, proceed to next iteration round
    else {
      cat(paste0("ARI: ",compclust$ARI,"\n"))

      # Update clustering to the predicted clusters
      ident_1 <- ident_2
      first_round <- FALSE
      metrics <- cbind(metrics,compclust)
      iterations = iterations + 1
      idents[[iterations]] <- ident_2

      ari <- compclust$ARI
      reiterations <- 0

    }
    # Step 4: If the maximum number of reiterations or iterations was reached, break the while loop
    if (reiterations == r | iterations == max.number.of.iterations)
    {
      toc(log = TRUE)
      break
    }
  }
  # Step 5: Return result

  return(list(probabilities=probs,metrics=metrics))
}
