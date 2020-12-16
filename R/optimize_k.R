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
#' @param cores numeric; number of cores to run in parallel in case of permutations
#' @return list with calculation of w at varying k and plot with results.
#' @import igraph
#' @importFrom "graphics" "legend" "lines"
#' @export

optimize_k <- function(Xs, X0, k, ind_adj, W, top = 200, cores = 4){

mND_list_k <- vector("list", length(k))

#Calculate mND
for(i in 1:length(mND_list_k)){
  mND_list_k[[i]] <- mND(Xs, ind_adj, k =  k[i], cores = 4)
}
names(mND_list_k) <- k

#Significance assessment
for(i in 1:length(mND_list_k)){
  mND_list_k[[i]] <- signif_assess(mND_list_k[[i]])
}

#Create the graph from adjacency matrix
g <- igraph::graph_from_adjacency_matrix(W, mode = "undirected",  weighted = T)

k_results <- eval_k(mND_list_k, X0, g, top = 200)

return(k_results)
}
