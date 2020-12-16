#' mND global function
#'
#' It's calculate permutations of X0 (if r parameter is used), applies network-diffusion on data, computes the mND score and the relative empirical p-value (if data are permuted).
#' Outputs of function can be used to classify genes in each layer with the classification function.
#' @param X0 matrix; each column (layer) of the matrix X0 is a score vector over all vertices of G.
#' @param W matrix; symmetrically normalized adjacency matrix W = D^-1 A D^-1, see 'normalize_adj_mat' function
#' @param alpha numeric; the smothing factor
#' @param nMax numeric; maximum number of iterations
#' @param eps numeric; the iteration will stop when the maximum difference between matrix Xs between two consecutive iteraction is smaller than \code{eps}
#' @param k numeric; number of top k first neighbors
#' @param r numeric; number of permutations
#' @param seed_n numeric; a given 'seed_n' will determine the same permutations of X0 and ensures the reproducibility of the analysis
#' @param cores numeric; number of cores to run in parallel
#' @import parallel
#' @return list with:
#' \itemize{
#' \item{\code{mND}}{: mND score (mND); if r parameter was used, the following columns will be added: corrisponding empirical p-value (p), product of mND score and -log10(p) (mNDp)}
#' \item{\code{t}}{: sum of top k first neighbours}
#' \item{\code{pl}}{: corrisponding empirical p-values of (t), only if r parameter was used}
#' \item{\code{tp}}{: product of the sum of top k neighbours (t) and -log10(pl), only if r parameter was used}
#' }
#' @export
#'
#'
framework_mND <- function(X0, W, alpha=0.7, nMax=1e4, eps=1e-6, k=3, r = NULL, seed_n = NULL, cores=1){


  #1) create permutation of X0
  if(!is.null(r)){
    X0 <- perm_X0(X0, r, W, seed_n)
  }

  #2) ND
  Xs <- ND(X0, W, alpha, nMax, eps, cores)

  #3) neighbour_index
  ind_adj <- neighbour_index(W)

  #4) mND
  mND_list <- mND(Xs, ind_adj, k, cores)

  if(!is.null(r)){
    mND_list <- signif_assess(mND_list)
  }

return(mND_list)
}
