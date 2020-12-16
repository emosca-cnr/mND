#' Gene classification
#'
#' Gene classification: genes are classified in every layer
#' @param mND mND scores; if data were permuted use the output of 'signif_assess' function, otherwise output of 'mND' function
#' @param X0 matrix; initial score matrix X0, in which each column (layer) represents a score vector over all vertices of G
#' @param Hl list; defined by the user, composed of high scoring genes names in each layer of X0.
#' @param top array; number of genes with the highest neighbourhoods to define the gene set Nl (see below); 1 number for each layer.
#' @param alpha array; significance level on the empirical p-value to define the gene set Nl (see below); 1 number for each layer.
#' @details
#' The gene set Nl contains the genes with the highest neighbourhoods in layer l. If both 'top' and 'alpha' are provided, a maximum of top genes with p < alpha will be considered.
#' @return list with:
#' \itemize{
#' \item{\code{gene_class}}{: data.frame with classifications of each gene in every layers.
#' The possible classification are: M (Module); L (Linker); I (Isolated); NS (Not Significant);}
#' \item{\code{occ_labels}}{: occurrence of the different labels (M; L; I; NS) for each gene across layers.}
#' }
#' @export
#'
classification <- function(mND, X0, Hl, top = NULL, alpha = NULL){

  if(is.null(top) & is.null(alpha)){
    stop("You must defined \'top\' or \'alpha'\ parameter.")
  }
  n_layer <- dim(X0)[2]
  if(length(Hl)!=n_layer){
    stop("Length of Hl list must be the same of layers' number")
  }

 if(dim(mND[[1]])[2]>2){
   mND_score <- mND$mND
   mND_score <- mND_score[order(mND_score$mNDp, decreasing = T),]
   sum_l <- mND$tp
   pl <- mND$pl
   if(!is.null(alpha) & is.null(top)){


     if(length(alpha)!=n_layer){
       stop("Length of \'alpha'\ array must be the same of layers' number")
     }

     #set variable
     ind_t =1
     val_sum <- rep(0,n_layer)
     for(i in 1:n_layer){
     val_sum[i] <- length(pl[pl[,i]<alpha[i],i])
     }
     ind_zero <- which(val_sum==0)
     if(length(ind_zero)!=0){
       string_l <- paste0(" l",ind_zero)
       stop("No genes with significant p-value found in layer:",string_l)
     }
   } else if(!is.null(alpha) & !is.null(top)){
     if(length(top)!=n_layer){
       stop("Length of \'top'\ array must be the same of layers' number")
     }
     if(length(alpha)!=n_layer){
       stop("Length of \'alpha'\ array must be the same of layers' number")
     }
     ind_t =2
     val_sum <- rep(0,n_layer)
     for(i in 1:n_layer){
       val_sum[i] <- length(pl[pl[,i]<alpha[i],i])
     }
     ind_zero <- which(val_sum==0)
     if(length(ind_zero)!=0){
       string_l <- paste0(" l",ind_zero)
       stop("No genes with significant p-value found in layer:",string_l)
     }
     val_sum <- top
   } else if(is.null(alpha) & !is.null(top)){
     ind_t =3
     if(length(top)!=n_layer){
       stop("Length of \'top'\ array  must be the same of layers' number")
     }
     val_sum <- top
   }
 }else if(dim(mND[[1]])[2]==1){
   mND_score <- mND$mND
   mND_score <- mND_score[order(mND_score$mND, decreasing = T),,drop=F]
   sum_l <- mND$t
   if(!is.null(alpha)){
     stop("Only \'top\' parameter is allowed: data were not permuted.")
   }
   if(length(top)!=n_layer){
     stop("Length of \'top'\ array must be the same of layers' number")
   }
   ind_t =3
   val_sum <- top
 }


  X0_new <- matrix(0, nrow=dim(X0)[1], ncol=dim(X0)[2], dimnames = list(rownames(X0), names(X0)))
  min_sum <- rep(0,n_layer)


  for(i in 1:n_layer){
    X0_new[which(rownames(X0_new) %in% Hl[[i]]), i] <- 1
  }

  gene_sign <- vector("list",n_layer)
  if(ind_t==1){
    for(i in 1:n_layer){
      gene_sign[[i]] <- rownames(pl[pl[,i]<alpha[i],i, drop=F])
    }
  } else if (ind_t==2){
    for(i in 1:n_layer){
      sum_l_ord <- sum_l[order(sum_l[,i], decreasing = T),]
      pl_ord <- pl[match(rownames(sum_l_ord), rownames(pl)),]
      dim_sign <- length(pl_ord[pl_ord[,i]<alpha[i],i])
      gene_sign[[i]] <- rownames(pl_ord[pl_ord[,i]<alpha[i],i,  drop=F])[1:min(dim_sign,val_sum[i])]
    }
  } else if (ind_t==3){
    for(i in 1:n_layer){
      sum_l_ord <- sum_l
      sum_l_ord$r <- rank(-sum_l[,i], ties.method = "min")
      gene_sign[[i]] <- rownames(sum_l_ord[sum_l_ord$r <=val_sum[i],])
    }
  }

  s_bin <- matrix(0, nrow=dim(sum_l)[1], ncol=dim(sum_l)[2], dimnames = list(rownames(sum_l), names(sum_l)))
  for(i in 1: (dim(sum_l)[2])){
     s_bin[which(rownames(s_bin) %in% gene_sign[[i]]),i] <- 2
  }


  X0_new <- X0_new[match(rownames(mND_score), rownames(X0_new)),]
  s_bin <- s_bin [match(rownames(mND_score), rownames(s_bin)),]
  calc_L <- X0_new+s_bin
  temp <-ifelse(calc_L==3, "M", ifelse(calc_L==2, "L",ifelse(calc_L==1, "I","NS" )))
  classification_mat <- data.frame(temp, stringsAsFactors = T)
  colnames(classification_mat) <- colnames(X0)

  levs <- unique( unlist( lapply(classification_mat, levels ) ) )
  tempBp <- apply(classification_mat,1,function(x) table(factor(x, levels=levs)))
  tempBp <- t(tempBp)
  tempBp <- as.data.frame(tempBp)
  tempBPN <- apply(tempBp,1,function(x) x/sum(x))
  tempBPN <- t(tempBPN)
  classification_res <- list(gene_class = classification_mat, occ_labels = tempBPN)
  return(classification_res)
}
