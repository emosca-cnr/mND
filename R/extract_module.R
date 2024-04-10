#' Extraction and plotting of connected subnetworks from a graph and a set of selected vertices
#' 
#' @description Extraction and plotting of connected subnetworks from a 
#' graph and a set of selected vertices
#' @param graph igraph object
#' @param selectedVertices subset of vertices
#' @param minSubnetSize minimum size of subnetwork
#' @param res_table result table
#' @return subnetwork 
#' @export
#' @importFrom igraph V<- V induced_subgraph set_vertex_attr clusters degree
#' @importFrom NPATools get_nconn_comp

extract_module <- function(graph=NULL, selectedVertices=NULL, res_table=NULL, minSubnetSize=2){
  
  #network between the selected gene
  dm <- induced_subgraph(graph, V(graph)$name[match(selectedVertices, V(graph)$name)])
  
  V(dm)$mNDp <- res_table$mNDp[match(V(dm)$name, rownames(res_table))]
  V(dm)$mND <- res_table$mND[match(V(dm)$name, rownames(res_table))]
  V(dm)$p <- res_table$p[match(V(dm)$name, rownames(res_table))]
  
  X0_cols <- 2:which(colnames(res_table)=="mND")
  for(i in X0_cols){
    dm <- set_vertex_attr(dm, colnames(res_table)[i], value = res_table[match(V(dm)$name, rownames(res_table)), i])
  }

  
  dm <- get_nconn_comp(dm, n = minSubnetSize)
  
  #get components and plot the graphs
  dmClusters <- clusters(dm)
  V(dm)$subnetId <- dmClusters$membership
  V(dm)$subnetSize <- dmClusters$csize[dmClusters$membership]
  V(dm)$d <- degree(dm)
  
   
  return(dm)
  
}
