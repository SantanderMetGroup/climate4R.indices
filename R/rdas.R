#' @title Exemplary data for \code{circIndexGrid.R}
#'
#' @description This exemplary grid contains monthly mean data of geopotential height at 500hPa from NCEP renalysis, for the period 2001-2010.
#' This R data object was obtained from the UDG server (\url{http://www.meteo.unican.es/udg-tap}, 
#' log-in is requiered, 
#' see \code{\link[loadeR]{loginUDG}}) by means of function \code{\link[loadeR]{loadGridData}} 
#' (package \href{https://github.com/SantanderMetGroup/loadeR}{\code{loadeR}}) in the following manner:
#' 
#' \code{loginUDG("username", "pasword")}\cr 
#'  
#' \code{NCEP_hgt500_2001_2010 <- loadGridData(dataset = "http://meteo.unican.es/tds5/dodsC/ncepReanalysis1/ncepReanalysis1_4xDaily.ncml", var = "hgt@500", latLim = c(20,90), lonLim=c(-180,180), years = 2001:2010, time="DD", aggr.d="mean", aggr.m="mean")}\cr
#' 
#' @format A grid object.
#' @usage data(NCEP_hgt500_2001_2010)
#' @name NCEP_hgt500_2001_2010
#' @docType data
#' 
#' @source \url{https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html} and \url{http://www.meteo.unican.es/udg-tap}.
#' @seealso \code{\link[loadeR]{loadGridData}}
#' @examples 
#' data("NCEP_hgt500_2001_2010")
#' # Plot in longlat projection
#' visualizeR::spatialPlot(transformeR::climatology(NCEP_hgt500_2001_2010), 
#' backdrop.theme = "coastline")
NULL


#' @title Exemplary data for \code{circIndexGrid.R}
#'
#' @description This exemplary grid contains monthly mean data of sea surface temperature from ERA-Interim renalysis, for the period 1981-2010.
#' This R data object was obtained from the UDG server (\url{http://www.meteo.unican.es/udg-tap}, 
#' log-in is requiered, 
#' see \code{\link[loadeR]{loginUDG}}) by means of function \code{\link[loadeR]{loadGridData}} 
#' (package \href{https://github.com/SantanderMetGroup/loadeR}{\code{loadeR}}) in the following manner:
#' 
#' \code{loginUDG("username", "pasword")}\cr 
#'  
#' \code{ERAInterim_sst_1981_2010 <- loadGridData(dataset = "http://meteo.unican.es/tds5/dodsC/interim/daily/interim20_daily.ncml", var = "sst", latLim = c(-6,6), lonLim=c(-170,-120), years = 1981:2010, time="DD", aggr.d="mean", aggr.m="mean")}\cr
#' 
#' @format A grid object.
#' @usage data(ERAInterim_sst_1981_2010)
#' @name ERAInterim_sst_1981_2010
#' @docType data
#' 
#' @source \url{https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era-interim} and \url{http://www.meteo.unican.es/udg-tap}.
#' @seealso \code{\link[loadeR]{loadGridData}}
#' @examples 
#' data("ERAInterim_sst_1981_2010")
#' # Plot 
#' visualizeR::spatialPlot(transformeR::climatology(ERAInterim_sst_1981_2010), 
#' backdrop.theme = "coastline")
NULL


#' @title Exemplary data for \code{circIndexGrid.R}
#'
#' @description This exemplary grid contains daily mean data of sea level pressure from NCEP renalysis, for the period 2001-2010, for a domain in the North-Atlantic.
#' This R data object was obtained from the UDG server (\url{http://www.meteo.unican.es/udg-tap}, 
#' log-in is requiered, 
#' see \code{\link[loadeR]{loginUDG}}) by means of function \code{\link[loadeR]{loadGridData}} 
#' (package \href{https://github.com/SantanderMetGroup/loadeR}{\code{loadeR}}) in the following manner:
#' 
#' \code{loginUDG("username", "pasword")}\cr 
#'  
#' \code{NCEP_slp_2001_2010 <- loadGridData(dataset = "http://meteo.unican.es/tds5/dodsC/ncepReanalysis1/ncepReanalysis1_4xDaily.ncml", var = "slp", latLim = c(30,70), lonLim = c(-70,10), years = 2001:2010, time="DD", aggr.d="mean")}\cr
#' 
#' @format A grid object.
#' @usage data(NCEP_slp_2001_2010)
#' @name NCEP_slp_2001_2010
#' @docType data
#' 
#' @source \url{https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html} and \url{http://www.meteo.unican.es/udg-tap}.
#' @seealso \code{\link[loadeR]{loadGridData}}
#' @examples 
#' data("NCEP_slp_2001_2010")
#' # Plot in longlat projection
#' visualizeR::spatialPlot(transformeR::climatology(NCEP_slp_2001_2010), 
#' backdrop.theme = "coastline")
NULL
