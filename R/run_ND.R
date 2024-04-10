#' Network diffusion
#'
#' Network diffusion: input scores (X0) are smoothed by network-diffusion, obtaining the corresponding network-constrained scores Xs.
#' @param X0 matrix or a list of matrices (if X0 was permuted (see 'perm_X0' function), where the first element of the list is the one obtained with real data); each column (layer) of the matrix X0 is a score vector over all vertices of G.
#' @param W  matrix; symmetrically normalized adjacency matrix W = D^-1 A D^-1, see 'normalize_adj_mat' function
#' @param alpha numeric; the smothing factor
#' @param nMax numeric; maximum number of iterations
#' @param eps numeric; the iteration will stop when the maximum difference between matrix Xs between two consecutive iteraction is smaller than \code{eps}
#' @param BPPARAM optional BiocParallelParam instance determining the parallel back-end to be used during evaluation. If NULL, parallel evaluation is disabled using SerialParam(). See ?bplapply.
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom NPATools ND
#' @return a matrix or a list of matrices (if data were permuted) with network diffusion scores.
#' @export

run_ND <- function(X0=NULL, W=NULL, alpha=0.7, nMax=1e4, eps=1e-6, BPPARAM = NULL){


  if (is.null(BPPARAM)) {
    BPPARAM <- SerialParam()
  }
  
  if(!is.list(X0)){
    Xs <- ND(X0 = X0, W = W, alpha = alpha, nMax = nMax, eps = eps)
  }
  if(is.list(X0)){
    Xs <-  bplapply(X0 , function(x) ND(X0 = x, W = W, alpha = alpha, nMax = nMax, eps = eps), BPPARAM = BPPARAM)
  }
  return(Xs)
}
