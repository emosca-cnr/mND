#'  mND score
#'
#' Calculation of mND score
#' @param Xs a matrix or list of matrices (if X0 was permuted, where the first element of the list is the one obtained with real data) with network diffusion scores, see 'ND' function
#' @param ind_adj list; index of neighbours obtained from adjacency matrix, see 'neighbour_index' function
#' @param k numeric; number of top k first neighbors
#' @param min.norm.val normalized ND values below this value will be set to 0; this is useful to denoise the steady state of ND and focus only on the highest scores 
#' @param BPPARAM optional BiocParallelParam instance determining the parallel back-end to be used during evaluation. If NULL, parallel evaluation is disabled using SerialParam(). See ?bplapply.
#' @importFrom BiocParallel SerialParam bplapply
#' @return list with:
#' \itemize{
#' \item{\code{mND}}{: mND score}
#' \item{\code{t}}{: sum of top k neighbours}
#' }
#' If data were permuted (therefore Xs is a list), the output list described above will be reported for each permutation of Xs.
#' @export

mND <- function(Xs=NULL, ind_adj=NULL, k=3, BPPARAM = NULL, min.norm.val=0){
  
  if (is.null(BPPARAM)) {
    BPPARAM <- SerialParam()
  }
  
  if(is.matrix(Xs)){
    mND_list <-  mND_s(output = Xs, ind_adj = ind_adj, k = k, min.norm.val = min.norm.val)
  }
  
  if(is.list(Xs)){
    mND_list <-  bplapply(Xs, function(x) mND_s(output = x, ind_adj = ind_adj, k = k, min.norm.val = min.norm.val), BPPARAM = BPPARAM)
  }
  return(mND_list)
}

