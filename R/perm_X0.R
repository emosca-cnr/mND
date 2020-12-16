#' Permutation of X0
#'
#' Permutation of input matrix X0
#' @param X0 matrix; a score matrix X0, in which each column (layer) is a score vector over all vertices of G.
#' @param r numeric; number of permutations.
#' @param W matrix; symmetrically normalized adjacency matrix W = D^-1 A D^-1, see 'normalize_adj_mat' function
#' @param seed_n optional numeric; to specify the seed for (pseudo) random number generation; using the same seed
#'
#' @return list with permutations of X0 (where the first element of the list is the one obtained with real data)
#'
#' @export
#'

perm_X0 <- function(X0, r, W, seed_n = NULL){

  if(!is.null(seed_n)){
    set.seed(seed_n)
  }
  X0_list <- c(list(X0), lapply(1:r, function(x) matrix(as.numeric(X0), ncol = ncol(X0), dimnames = list(sample(rownames(X0), nrow(X0))))))
  X0_list <- lapply(X0_list, function(x) x[match(rownames(W), rownames(x)), , drop = F])

  return(X0_list)
}
