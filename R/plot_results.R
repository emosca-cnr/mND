#' Plot results
#'
#' Visualization of results.
#' @param mND mND scores; if data were permuted use the output of 'signif_assess' function, otherwise output of 'mND' function.
#' @param classification data.frame; gene classification, output of 'classification' function
#' @param W matrix; symmetrically normalized adjacency matrix W = D^-1 A D^-1, see 'normalize_adj_mat' function
#' @param n numeric; number of top n ranking genes to consider in the plot of gene networks
#' @param directory string; output directory
#' @importFrom grDevices dev.off jpeg
#' @importFrom graphics abline axis barplot image layout layout.show par plot
#' @import igraph
#' @details
#' The possible plots are:
#' \itemize{
#' \item{\code{genes ranked by mND score and the corresponding p-value (for permutated data)}}{: mND argument is required (output of 'signif_assess' function);}
#' \item{\code{gene networks composed of the top n ranking genes}}{: mND, W and n arguments are required;}
#' \item{\code{gene classification for the top 100 ranking genes across layers}}{: classification argument is required. Legend: brown: isolated; orange: linker; purple: module; grey: not significant.}
#' }
#' @export
#'
plot_results <- function(mND = NULL, classification = NULL, W = NULL, n = NULL, directory = NULL){

  if(is.null(mND) & is.null(classification)){
    stop("You must defined at least \'mND\' or \'classification'\ arguments.")
  }
  if(is.null(directory)){
    stop("You have to define the output directory")
  }
  if (!dir.exists(directory)) {
    dir.create(directory)
  }
  #plot gene network
  if(!is.null(mND) & !is.null(n) & is.null(W)){
   cat("W argument is not defined, gene networks cannot be plotted.")
  }

  if(is.null(mND) & !is.null(n) & !is.null(W)){
    cat("mND argument is not defined, gene networks cannot be plotted.")
  }

  if(!is.null(mND) & is.null(n) & !is.null(W)){
    cat("n argument is not defined, gene networks cannot be plotted.")
  }
  if(!is.null(mND) & is.null(n) & is.null(W) & dim(mND[[1]])[2]<=2){
   stop("genes ranked by mND score and the corresponding p-value cannot be plotted: mND input is not in a correct format.")
  }

  #permutation
  if(!is.null(mND)){
    if(dim(mND[[1]])[2]>2){
      mND_score <- mND$mND
      jpeg(paste0(directory,"/mND_vs_pvalue.jpeg"), width = 180, height = 180, units = "mm", res=300)
      plot(mND_score$mND, -log10(mND_score$p),col="black", cex=0.7, pch=16,ylab="-log10(p)",xlab="mND")
      abline(h=-log10(0.05), lty=2, lwd=2, col="orange")
      dev.off()
      if(!is.null(n)){
        mND_score <- mND_score[order(mND_score$mNDp, decreasing = T),]
        mND_score <- mND_score[1:n,]
      }
    } else{
      if(!is.null(n)){
        mND_score <- mND$mND
        mND_score <- mND_score[order(mND_score$mND, decreasing = T),,drop=F]
        mND_score <- mND_score[1:n,,drop=F]
      }
    }
    if(!is.null(W) & !is.null(n)){
      g <- igraph::graph_from_adjacency_matrix(W, weighted = T, mode = "undirected")
      g_sub <- igraph::induced_subgraph(g, V(g)$name %in% rownames(mND_score))
      jpeg(paste0(directory,"/gene_network.jpeg"), width = 180, height = 180, units = "mm", res = 300)
      par(mar=c(1,1,1,1))
      igraph::plot.igraph(g_sub,
                  vertex.label=NA,
                  vertex.size=5,
                  edge.width=0.8,
                  edge.color='hotpink',
                  edge.lty=2,
                  edge.width = 0.7)
      dev.off()
    }
  }
  if(!is.null(classification)){
    class_plot(classification$gene_class, directory, 100)
  }
}
