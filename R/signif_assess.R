#' Significance assessment
#'
#' Significance assessment: empirical p-values are calculated by using the pool of r permuted versions of X0 (clearly, use this function only if 'perm_X0' function was used).
#' @param mND list; output of 'mND' function
#' @return list with:
#' \itemize{
#' \item{\code{mND}}{: mND score (mND), corrisponding empirical p-value (p), product of mND score and -log10(p) (mNDp)}
#' \item{\code{t}}{: sum of top k first neighbours}
#' \item{\code{pl}}{: corrisponding empirical p-values of (t)}
#' \item{\code{tp}}{: product of the sum of top k neighbours (t) and -log10(pl)}
#' }
#' @export

signif_assess <- function(mND){

  mND_perm <- lapply(mND, function(x) x[[1]])
  p_val <- calc_p(mND_perm)
  vect_tru <- mND[[1]][[1]]
  p_val <- p_val[match(rownames(vect_tru), rownames(p_val)),]
  mND_score <- data.frame(mND=vect_tru,p=p_val,mNDp=vect_tru*-log10(p_val),stringsAsFactors = F)
  rownames(mND_score) <- rownames(vect_tru)
  colnames(mND_score) <- c("mND","p", "mNDp")


  sum_perm <- lapply(mND, function(x) x[[2]])
  ###pval sum
  p_val_sum <- vector("list", dim(sum_perm[[1]])[2])
  for(i in 1:dim(sum_perm[[1]])[2]){
    sum_perm_li <- lapply(sum_perm, function(x) matrix(x[,i], nrow=nrow(x), ncol=1, dimnames = list(rownames(x), names(x)[i])))
    p_val_sum[[i]] <- calc_p(sum_perm_li)
  }

  temp <- do.call("cbind",p_val_sum)

  sum_true <- sum_perm[[1]]
  temp <- temp[,match(colnames(sum_true), colnames(temp))]
  temp <- temp[match(rownames(sum_true), rownames(temp)),]

  sum_pval <- sum_true*-log10(temp)
  colnames(sum_pval) <- paste0("tp", seq(1:dim(sum_pval)[2]))

  mND_list <- list(mND = mND_score,  t = mND[[1]]$t, pl=temp, tp = sum_pval)

  return(mND_list)
}
