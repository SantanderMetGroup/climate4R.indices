#' @title Photovoltaic potential index (pvpot)
#' 
#' @description Photovoltaic potential index describing cell potential production
#' with respect to the optimal potential at a global downward shortwave radiation
#' of 1000 W m^-2.
#'
#' @param tm Numeric vector, array or matrix with mean air temperature in °C.
#' @param wss Numeric vector, array or matrix with mean wind speed in m s^-1.
#' @param rad Numeric vector, array or matrix with mean surface solar radiation in
#'   W m^-2. If radiation is provided as accumulated energy in J m^-2, it must be
#'   converted beforehand by dividing by the accumulation period in seconds
#'   (e.g. 86400 for daily data, 3600 for hourly data).
#' @param clamp Logical. If TRUE, clamp output to the [0, 1] interval.
#'
#' @return Numeric vector, array or matrix with the photovoltaic potential index.
#'
#' @export

pvpot <- function(tm, wss, rad, clamp = FALSE) {
  
  if (!identical(dim(tm), dim(wss)) || !identical(dim(tm), dim(rad))) {
    stop("tm, wss and rad must have the same dimensions.")
  }
  
  c1 <- 0.943
  c2 <- 0.028
  c3 <- 1.528
  c4 <- 4.3
  
  beta <- 0.005
  t_ref <- 25
  
  tcell <- c1 * tm + c2 * rad - c3 * wss + c4
  pr <- 1 - beta * (tcell - t_ref)
  
  out <- pr * rad / 1000
  
  if (clamp) {
    out <- pmax(0, pmin(1, out))
  }
  
  out
}