#' Evaluation of the parameter k
#'
#' @param mND list; list of mND score obtained for different k values, list components should be named as the k values used
#' @param X0 matrix; initial score matrix X0, in which each column (layer) represents a score vector over all vertices of G
#' @param top numeric; number of genes to be considered at the top of the list of genes ranked by mND
#' @param G graph; gene x gene undirected interaction graph
#' @return list with calculation of w at varying k and plot with results.
#' @importFrom graphics legend lines
#' @importFrom dmfind omega

eval_k <- function (mND=NULL, X0=NULL, G=NULL, top = 200){

  n_k <- length(mND)
  if(dim(mND[[1]][[1]])[2]>2){
    mND <- lapply(mND, function(x) x$mND[,"mNDp",drop=F])
  } else{
    mND <- lapply(mND, function(x) x$mND[,"mND",drop=F])
  }

  mND <- lapply(mND, function(x) x[match(rownames(mND[[1]]), rownames(x)), , drop = F])
  X0 <- X0[match(rownames(mND[[1]]), rownames(X0)), ,drop=F]

  ans <- vector("list", length(mND))
  csX0 <- ans
  csmND <- ans
  #for each k
  for (i in 1:n_k) {
    #order by mND
    ans[[i]] <- mND[[i]]
    ref_order <- order(-ans[[i]])[1:top]
    csmND[[i]] <- omega(G, array(ans[[i]][, 1], dimnames = list(rownames(ans[[i]])))[ref_order])

    #omega X0 for each layer
    csX0[[i]] <- sapply(1:ncol(X0), function(j) omega(G, array(X0[, j], dimnames = list(rownames(X0)))[ref_order]))

  }

  #normalization
  #max of each layer across all k-
  oX0_max <- do.call(rbind, lapply(csX0, function(x) apply(x, 2, max)))
  oX0_max <- apply(oX0_max, 2, function(icol) icol[which.max(icol)])
  for (i in 1:n_k){
    csX0[[i]] <- t(apply(csX0[[i]], 1, function(irow) irow/oX0_max))
  }

  #max of mND across layers
  o_mND_max <- max(unlist(csmND))
  csmND <- lapply(csmND, function(x) x / o_mND_max)

  #optimal: sum of x0 from the two layers
  optimal <- csmND
  for(i in 1:n_k){
    optimal[[i]] <- rowSums(csX0[[i]])
  }
  names(optimal) <- names(mND)
  j<-1
  plot(optimal[[j]], type="l", ylim=c(0, max(unlist(optimal))), ylab=expression(omega), xlab="n")
  for(j in 2:length(mND)){ #k
    lines(optimal[[j]], col=j)
  }
  k <- names(mND)
  legend("bottomright", legend = k, col=1:length(k), lty= 1, cex=0.7)

  return(optimal)
}
