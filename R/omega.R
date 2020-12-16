#' Omega
#'
#' Calculation of omega function
#' @param G graph; gene x gene undirected interaction graph
#' @param u list; ranked names that have to be those in V(G)$name
#' @import igraph

omega <- function (G, u)
{
  Gi <- igraph::induced.subgraph(G, match(names(u), V(G)$name))
  Ai <- as.matrix(get.adjacency(Gi))
  idx.norm <- match(names(u), rownames(Ai))
  Ai <- Ai[idx.norm, idx.norm]
  U <- matrix(u, ncol = 1, dimnames = list(names(u)))
  omega_vect <- U %*% t(U)
  omega_vect <- omega_vect * Ai
  omega_vect[upper.tri(omega_vect)] <- 0
  omega_vect <- 2 * cumsum(rowSums(omega_vect))
  return(omega_vect)
}
