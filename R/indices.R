#' @title Days with maximum temperature above a given threshold
#' @description Annual count of days with maximum temperature above a given threshold 
#' @param tx Vector with minimum temperature data
#' @param th Threshold value
#' @author M. Iturbide
#' @export
 
tx.th <- function(tx, th) {sum(tx > th, na.rm = TRUE)}

#end

#' @title Days with minimum temperature below a given threshold
#' @description Annual count of days with maximum temperature above a given threshold 
#' @param tn Vector with minimum temperature data
#' @param th Threshold value (Default is 0)
#' @author M. Iturbide
#' @export

tn.th <- function(tn, th = 0) {sum(tn < th, na.rm = TRUE)}

#end


#' @title Growing Degree Days 
#' @description Accumulated sum of the difference between daily mean temperature and the threshold (when higher than the threshold) 
#' over the primary growing season for mid-latitude agricultural areas in the northern Hemisphere (April-September).
#' @param tm Vector with mean temperature data
#' @param th Threshold value (Defalut is 5º)
#' @details 
#' @author M. Iturbide
#' @export

gdd.th <- function(tm, th = 5) {
  ind <- which(tm > th)
  sum(tm[ind] - th, na.rm = T)
}
#end


#' @title Heating Degree Days (HDD)
#' @description Symmetrical to the Cooling Degree Day index (CDD), the HDD index is used for illustrating 
#' energy demand for heating.
#' @param tn Matrix or vector with minimum temperature data
#' @param tx Matrix or vector with maximum temperature data
#' @param tm Matrix or vector with mean temperature data
#' @param th Threshold value (Defalut is 15.5º)
#' @author M. Iturbide
#' @references Spinoni, J., Vogt, J., and Barbosa, P. European degree-day climatologies and trends for the period 1951–2011. Int. J. Climatol. 35, 25–36. doi:10.1002/joc.3959.
#' 
#' Spinoni, J., Vogt, J. V, Barbosa, P., Dosio, A., McCormick, N., Bigano, A., et al. Changes of heating and cooling degree-days in Europe from 1981 to 2100. Int. J. Climatol. 38, e191–e208. doi:10.1002/joc.5362.
#' @export

hdd.th <- function(tn, tx, tm, th = 15.5) {
  arr <- tm
  arr[which(th >= tx)] <- th - tm[which(th >= tx)]
  arr[which(th < tx & tm <= th)] <- ((th - tn)/2) - ((tx - th)/4)
  arr[which(th < tm & tn <= th)] <- (th - tn)/4
  arr[which(th <= tn)] <- 0
  apply(arr, MARGIN = 2:3, FUN = "sum", na.rm = TRUE)
}

#end

#' @title Cooling Degree Day (CDD)
#' @description Symmetrical to the Heating Degree Day index (HDD), the CDD index is used for illustrating 
#' energy demand for cooling.
#' @param tn Matrix or vector with minimum temperature data
#' @param tx Matrix or vector with maximum temperature data
#' @param tm Matrix or vector with mean temperature data
#' @param th Threshold value (Defalut is 22º)
#' @author M. Iturbide
#' @references Spinoni, J., Vogt, J., and Barbosa, P. European degree-day climatologies and trends for the period 1951–2011. Int. J. Climatol. 35, 25–36. doi:10.1002/joc.3959.
#' 
#' Spinoni, J., Vogt, J. V, Barbosa, P., Dosio, A., McCormick, N., Bigano, A., et al. Changes of heating and cooling degree-days in Europe from 1981 to 2100. Int. J. Climatol. 38, e191–e208. doi:10.1002/joc.5362.
#' @export

cdd.th <- function(tn, tx, tm, th = 22) {
  arr <- tm
  th <- 22
  arr[which(th >= tx)] <- 0
  arr[which(th < tx & tm <= th)] <- (tx - th)/4
  arr[which(th < tm & tn <= th)] <- ((tx - th)/2) - ((th - tn)/4)
  arr[which(th <= tn)] <- tm - th
  apply(arr, MARGIN = 2:3, FUN = "sum", na.rm = TRUE)
}

#end

#' @title Percentile
#' @description Calculate the value of a given probability value or the other way around.
#' @param var Vector 
#' @param percent Number between 0 and 1 to calculate the value of a given probability.
#' @param value Number of the value for which the probability is going to be calculated.
#' @author M. Iturbide
#' @export
#' @importFrom stats ecdf quantile
#' @importFrom abind abind

percentile <- function(var, percent = NULL, value = NULL){
  if (!is.null(percent) & !is.null(value)) {
    value <- NULL
    warning("Values were given to both percentile and value... value will be ignored (set to NULL)")
  }
  argument <- percent
  if (!is.null(value)) argument <- value
  if (!is.matrix(argument)) {
    if (!is.null(percent)) {
      tryCatch({quantile(var, probs = argument/100, na.rm = TRUE)}, error = function(err){NA})
    } else if (!is.null(value)) {
      tryCatch({ecdf(var)(argument)*100}, error = function(err){NA})
    }
  } else {
  li <- lapply(1:dim(argument)[1], function(di) {
    lj <- lapply(1:dim(argument)[2], function(dj) {
      if (!is.null(percent)) {
        tryCatch({quantile(var[, di, dj], probs = argument[di, dj]/100, na.rm = TRUE)}, error = function(err){NA})
      } else if (!is.null(value)) {
        tryCatch({ecdf(var[, di, dj])(argument[di, dj])*100}, error = function(err){NA})
      }
    })
    do.call("abind", lj)
  })  
  do.call("abind", list(li, along = 0))
  }
}


#end


