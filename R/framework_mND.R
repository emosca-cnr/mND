#' mND global function
#'
#' It calculates permutations of X0 (if r parameter is used), applies network-diffusion on data, computes the mND score and the relative empirical p-value (if data are permuted).
#' Outputs of function can be used to classify genes in each layer with the classification function.
#' @param X0 matrix; each column (layer) of the matrix X0 is a score vector over all vertices of G.
#' @param W matrix; symmetrically normalized adjacency matrix W = D^-1 A D^-1, see 'normalize_adj_mat' function
#' @param alpha numeric; the smothing factor
#' @param nMax numeric; maximum number of iterations
#' @param eps numeric; the iteration will stop when the maximum difference between matrix Xs between two consecutive iteraction is smaller than \code{eps}
#' @param k numeric; number of top k first neighbors
#' @param r numeric; number of permutations
#' @param seed_n numeric; a given 'seed_n' will determine the same permutations of X0 and ensures the reproducibility of the analysis
#' @param BPPARAM optional BiocParallelParam instance determining the parallel back-end to be used during evaluation. If NULL, parallel evaluation is disabled using SerialParam(). See ?bplapply.
#' @param ... additional parameter to NPATools::perm_vertices()
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom NPATools perm_vertices perm_X0
#' @importFrom stats setNames
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
framework_mND <- function(X0=NULL, W=NULL, alpha=0.7, nMax=1e4, eps=1e-6, k=3, r = 0, seed_n = NULL, BPPARAM=NULL, ...){

  if (is.null(BPPARAM)) {
    BPPARAM <- SerialParam()
  }
  
  #1) create permutation of X0
  if(r>0){
    vert_deg <- setNames(rowSums(sign(W)), rownames(W))
    perms <- perm_vertices(vert_deg = vert_deg, seed_n = seed_n, k = r, ...)
    X0 <- perm_X0(X0 = X0, perms=perms)
  }

  #2) ND
  Xs <- run_ND(X0 = X0, W = W, alpha = alpha, nMax = nMax, eps = eps, BPPARAM = BPPARAM)

  #3) neighbour_index
  ind_adj <- neighbour_index(W)

  #4) mND
  mND_list <- mND(Xs = Xs, ind_adj = ind_adj, k = k, BPPARAM = BPPARAM)

  if(!is.null(r)){
    mND_list <- signif_assess(mND = mND_list, BPPARAM = BPPARAM)
  }

return(mND_list)
}
