#'  Index neighbours
#'
#' Indices of neighbours for each gene
#' @param W matrix; symmetrically normalized adjacency matrix W = D^-1 A D^-1, see 'normalize_adj_mat' function
#' @return list with indices of neighbours for each gene
#'
#' @export
#'

neighbour_index <- function(W=NULL){

  ind_adj <- apply(W, 1, FUN=function(x) which(x>0, arr.ind = T))
  
  return(ind_adj)
}
