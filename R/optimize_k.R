#' Optimization of k value
#'
#' Optimization of the 'k' value by selecting a value that yields connected networks enriched in initial scores.
#'
#' @param Xs a matrix or list of matrices (if X0 was permuted, where the first element of the list is the one obtained with real data) with network diffusion scores, see 'ND' function
#' @param X0 matrix; initial score matrix X0, in which each column (layer) represents a score vector over all vertices of G
#' @param k  array; array with k values
#' @param ind_adj list; index of neighbours obtained from adjacency matrix, see 'neighbour_index' function
#' @param W matrix; symmetrically normalized adjacency matrix W = D^-1 A D^-1, see 'normalize_adj_mat' function
#' @param top numeric; number of genes to be considered at the top of the list of genes ranked by mND
#' @param BPPARAM optional BiocParallelParam instance determining the parallel back-end to be used during evaluation. If NULL, parallel evaluation is disabled using SerialParam(). See ?bplapply.
#' @return list with calculation of w at varying k and plot with results.
#' @importFrom BiocParallel SerialParam
#' @importFrom igraph graph_from_adjacency_matrix
#' @export

optimize_k <- function(Xs=NULL, X0=NULL, k=NULL, ind_adj=NULL, W=NULL, top = 200, BPPARAM = NULL){
  
  
  if (is.null(BPPARAM)) {
    BPPARAM <- SerialParam()
  }
  
  mND_list_k <- vector("list", length(k))
  
  #Calculate mND
  for(i in 1:length(mND_list_k)){
    mND_list_k[[i]] <- mND(Xs, ind_adj, k =  k[i], BPPARAM = BPPARAM)
  }
  names(mND_list_k) <- k
  
  #Significance assessment
  for(i in 1:length(mND_list_k)){
    mND_list_k[[i]] <- signif_assess(mND_list_k[[i]])
  }
  
  #Create the graph from adjacency matrix
  g <- graph_from_adjacency_matrix(W, mode = "undirected",  weighted = T)
  
  k_results <- eval_k(mND_list_k, X0, g, top = 200)
  
  return(k_results)
}
