#' Integration with mND function
#' @param output Xs
#' @param ind_adj index of neighbours obtained from adj matrix, see 'neighbour_index' function
#' @param k number of top k first neighbors
#' @return a list with:
#' \itemize{
#' \item{\code{mND}}{: mND score}
#' \item{\code{t}}{: sum of top k neighbours}
#' }
#'


mND_s <- function(output, ind_adj, k){

  Xs<-  apply(output, 2, function(x) x/max(x))
  Xs <- as.data.frame(Xs)
  Xs <- Xs[match(names(ind_adj), rownames(Xs)),]

  if(identical(rownames(Xs), names(ind_adj))==FALSE){
    stop("Rownames of Xs and W are to be in the same order.")
  }

  layer_number <- dim(Xs)[2]
  colnames(Xs) <-paste0("Xs_l", seq(1:layer_number))
  sum_df <- matrix(0, nrow = dim(Xs)[1], ncol = layer_number)
  c_n_df2 <- paste0("t",seq(1:layer_number))
  rownames(sum_df) <- rownames(Xs)
  colnames(sum_df) <- c_n_df2
  #sum_df <- data.frame(gene = rownames(Xs), stringsAsFactors = F)
  #cycle through genes
  for(i in 1:length(ind_adj)){
    temp_Xs <- Xs[ind_adj[[i]], , drop=FALSE]
    temp_Xs <- as.matrix(temp_Xs)
    ki <-length(ind_adj[[i]])
    if(ki <= k){
      kii = ki
    }else{
      kii = k
    }

    temp_Xs_rank <- apply(temp_Xs, 2, function(x) rank(-x, ties.method = "random"))
    temp_Xs[temp_Xs_rank>kii] <- 0
    sum_neigh <- colSums(temp_Xs)/kii
    names(sum_neigh) <- c_n_df2
    sum_df[i,] <- sum_neigh
  }
  colnames(sum_df) <- c_n_df2
  sum_df <- as.data.frame(sum_df)
  mND_temp <- data.frame(Xs_sum = rowSums(Xs), sum_df = rowSums(sum_df), stringsAsFactors = F)
  mND_temp$score <- mND_temp$Xs_sum*mND_temp$sum_df
  rownames(mND_temp) <- rownames(Xs)
  mND_m <- mND_temp[,"score", drop = F]

  list_mND <- list(mND = mND_m, t =  sum_df)

  return(list_mND)
}

