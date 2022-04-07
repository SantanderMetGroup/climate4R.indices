#' @title Threshold application for Sinusoidal Modified Growing Degree Days calculation 
#' @description Threshold application for Sinusoidal Modified Growing Degree Days calculation 
#' (see \code{mgdd.th}).
#' @param tm Vector with HOURLY minimum temperature data
#' @param th1 Numeric.
#' @param th2 Numeric.
#' @param th3 Numeric.
#' @param th4 Numeric.
#' @param last Numeric.
#' @details Defaults for the th* parameters are: th1 = 12, th2 = 18, th3 = 28, th4 = 32, last = 35.
#' @author M. Iturbide
#' @export

mgdd.aux <- function(tm, th1 = 12, th2 = 18, th3 = 28, th4 = 32, last = 35){
  arr <- tm
  ind1 <- which(tm > th1)
  ind0 <- which(tm <= th1)
  tm1 <- tm[ind1]
  arr1 <-tm1
  arr1[which(tm1 < th2)] <- 2/3 * (tm1[which(tm1 < th2)] - th1)
  arr1[which(tm1 >= th2 & tm1 < th3)] <- 2/3 * (th2 - th1) + (tm1[which(tm1 >= th2 & tm1 < th3)] - th2)
  arr1[which(tm1 >= th3 & tm1 < th4)] <-  2/3 * (th2 - th1) + (th3 - th2) - 1.25*(tm1[which(tm1 >= th3 & tm1 < th4)] - th3)
  arr1[which(tm1 >= th4 & tm1 < last)] <- 2/3 * (th2 - th1) + (th3 - th2) - 1.25*(th4 - th3) - 3*(tm1[which(tm1 >= th4 & tm1 < last)] - th4)
  arr1[which(tm1 >= last)] <- 0
  arr[ind1] <- arr1
  arr[ind0] <- 0
  arr
}


#end
