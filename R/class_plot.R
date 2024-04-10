#' Classification plot
#'
#' Visulation of classification
#' @param classification_mat data.frame; see 'classification' function
#' @param directory string; output directory
#' @importFrom grDevices dev.off jpeg
#' @importFrom graphics abline axis barplot image layout layout.show par plot
#' @param n numeric; number of top n genes ranked by mND to plot
#'
class_plot <- function(classification_mat=NULL, directory=NULL, n=NULL){

  calc_L <- classification_mat
  levs <- unique( unlist( lapply(classification_mat, levels ) ) )
  
  #  Set these as the levels for each column
  tempBp <- apply(classification_mat,1,function(x) table(factor(x, levels=levs)))
  tempBp <- t(tempBp)
  tempBp <- as.data.frame(tempBp)
  tempBPN <- apply(tempBp,1,function(x) x/sum(x))
  tempBPN <- t(tempBPN)
  tempN <- tempBPN[1:n,c(1:3)]
  rownames(tempN) <- NULL
  rotate <- function(x) t(apply(x, 2, rev))
  tempN <- rotate(tempN)
  tempL <- tempBp[tempBp$L>0 & tempBp$M==0 & tempBp$I==0,]

  nL <- dim(classification_mat)[2]
  temp2 <-ifelse(classification_mat=="M", 3, ifelse(classification_mat=="L", 2,ifelse(classification_mat=="I", 1,0 )))
  temp2 <- temp2[1:n,]
  tt <- rotate(temp2)
  
  jpeg(file.path(directory, "classification.jpeg"), width = 25*nL, height = 400, units = "mm", res=300)
  if(nL>2){
    par(mar = c(2.9,6,2.9,0))
    layout.show(layout(matrix(c(1,2), nrow=1, byrow = T), widths =c(0.8, 0.2)))
    image(tt, col =c('0' = "gray75", '1' = "chocolate4", '2' = "chocolate1", '3' = "purple"),xaxt="none", yaxt="none")
    xx <- seq(0, 1, length.out = ncol(temp2))
    yy <- seq(1, 0, length.out = nrow(temp2))
    aa<-seq(0,n)
    bb <- seq(0,nL)
    cc <- aa*yy[n-1]+yy[n-1]/2
    vv <- xx+xx[2]/2
    abline(v=vv, col="black", xpd=F)
    abline(h=cc, col="black", xpd=F)
    axis(1, xx,colnames(temp2),  cex.axis=1)
    axis(2, yy,rownames(temp2),  cex.axis=1, las=2)
    par(mar = c(0,0.7,0,0.5))
    barplot((tempN), horiz = T, xaxt="n",space=0, col=c("chocolate4", "chocolate1","purple"))
  } else {
    par(mar = c(2.9,6,2.9,0.5))
    image(tt, col =c('0' = "gray75", '1' = "chocolate4", '2' = "chocolate1", '3' = "purple"),xaxt="none", yaxt="none")
    xx <- seq(0, 1, length.out = ncol(temp2))
    yy <- seq(1, 0, length.out = nrow(temp2))
    aa<-seq(0,n)
    bb <- seq(0,nL)
    cc <- aa*yy[n-1]+yy[n-1]/2
    vv <- xx+xx[2]/2
    abline(v=vv, col="black", xpd=F)
    abline(h=cc, col="black", xpd=F)
    axis(1, xx,colnames(temp2),  cex.axis=1)
    axis(2, yy,rownames(temp2),  cex.axis=1, las=2)
  }
  dev.off()
}
