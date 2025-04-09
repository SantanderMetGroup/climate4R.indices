#' @title Photovoltaic potential index (pvpot)
#' 
#' @description Photovoltaic potential index describing cell potential production with respect to the optimal potential of a global downward shortwave radiation of 1000 W/m2
#' 
#' @param tm Matrix or vector with mean temperature data (in Celsius ºC)
#' @param wss Matrix or vector with mean wind speed data (in meters per second m/s)
#' @param rad Matrix or vector with the mean radiation data (in Watios per squared meters W/m²)
#' 
#' @return A numeric vector or matrix with the photovoltaic potential index, between 0-1.
#' 
#' @author O.Mirones
#' 
#' @references Chenni, R., M. Makhlouf, T. Kerbache, and A. Bouzid. 2007. A Detailed Modeling Method for Photovoltaic Cells. Energy 32 (9): 1724–30. https://doi.org/10.1016/j.energy.2006.12.006.
#' 
#' Crook, Julia A., Laura A. Jones, Piers M. Forster, and Rolf Crook. 2011. “Climate Change Impacts on Future Photovoltaic and Concentrated Solar Power Energy Output.” Energy & Environmental Science 4 (9): 3101. https://doi.org/10.1039/c1ee01495a.
#' 
#' Jerez, Sonia, Isabelle Tobin, Robert Vautard, Juan Pedro Montávez, Jose María López-Romero, Françoise Thais, Blanka Bartok, et al. 2015. “The Impact of Climate Change on Photovoltaic Power Generation in Europe.” Nature Communications 6 (1): 10014. https://doi.org/10.1038/ncomms10014.

#' @export

pvpot <- function(tm, wss, rad){
  
  #c1, c2, c3 and c4 are internal parameters of the function. 
  c1 <- 0.943
  c2 <- 0.028
  c3 <- 1.528
  c4 <- 4.3
  
  tcell <- tm*c1 + rad*c2 - wss*c3 + c4
  
  #beta is also an internal parameters of the function. By default, Ta is 25ºC
  beta <- 0.005
  Ta <- 25
  
  pr <- 1 - beta*(tcell - Ta)
  
  pvpot <- pr*rad/1000
  
  return(pvpot)
}