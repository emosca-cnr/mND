#' Network diffusion
#'
#' Network diffusion: input scores (X0) are smoothed by network-diffusion, obtaining the corresponding network-constrained scores Xs.
#' @param X0 matrix or a list of matrices (if X0 was permuted (see 'perm_X0' function), where the first element of the list is the one obtained with real data); each column (layer) of the matrix X0 is a score vector over all vertices of G.
#' @param W  matrix; symmetrically normalized adjacency matrix W = D^-1 A D^-1, see 'normalize_adj_mat' function
#' @param alpha numeric; the smothing factor
#' @param nMax numeric; maximum number of iterations
#' @param eps numeric; the iteration will stop when the maximum difference between matrix Xs between two consecutive iteraction is smaller than \code{eps}
#' @param cores numeric; number of cores to run in parallel in case of permutations
#' @import parallel
#' @return a matrix or a list of matrices (if data were permuted) with network diffusion scores.
#' @export

ND <- function(X0, W, alpha=0.7, nMax=1e4, eps=1e-6, cores=1){


  if(is.matrix(X0)){
    Xs <- ND_s(X0, W,  alpha, nMax, eps)$Xs
  }
  if(is.list(X0)){
    Xs <-  parallel::mclapply(X0 , function(x) ND_s(as.matrix(x), W, alpha, nMax, eps)$Xs, mc.cores = cores)
  }
  return(Xs)
}
