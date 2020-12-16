#'  mND score
#'
#' Calculation of mND score
#' @param Xs a matrix or list of matrices (if X0 was permuted, where the first element of the list is the one obtained with real data) with network diffusion scores, see 'ND' function
#' @param ind_adj list; index of neighbours obtained from adjacency matrix, see 'neighbour_index' function
#' @param k numeric; number of top k first neighbors
#' @param cores numeric; number of cores to run in parallel in case of permutations
#' @import parallel
#' @return list with:
#' \itemize{
#' \item{\code{mND}}{: mND score}
#' \item{\code{t}}{: sum of top k neighbours}
#' }
#' If data were permuted (therefore Xs is a list), the output list described above will be reported for each permutation of Xs.
#' @export
mND <- function(Xs, ind_adj, k=3, cores = 1){

  if(is.matrix(Xs)){
  mND_list <-  mND_s(Xs, ind_adj, k)
  }
  if(is.list(Xs)){
    mND_list <-  parallel::mclapply(Xs, function(x) mND_s(x, ind_adj,k), mc.cores = cores)
   }
  return(mND_list)
}

