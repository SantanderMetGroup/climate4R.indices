#     indicesENSO.R Calculation of the CPC clirculation indices from grid
#
#     Copyright (C) 2019 Santander Meteorology Group (http://www.meteo.unican.es)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title Calculation of ENSO indices from grid.
#' @description Calculate ENSO indices of grids or multimember grids. 
#' @param data A grid (gridded or station dataset), or multimember grid object of sea surface temperature.
#' @param index.code Circulation index (or vector of indices) to be computed. See \code{?indexCircShow()} for details.
#' @param base Baseline grid to be substracted for the calculation of anomalies. Default: NULL. See \code{?scaleGrid}.
#' @param ref Reference grid to be added for the calculation of anomalies. Default: NULL. See \code{?scaleGrid}.
#' @param season Selected month(s) for the calculation. Default: NULL (i.e. as input grid).
#' @param members Select number of members.

#' @return A list of circulation indices (and members, if applicable) with:
#' \itemize{
#' * index: vector with the time series of the teleconnection index.
#' * pattern: matrix with the spatial pattern of the teleconnection.
#' * dates and coordinates as list attributes.
#' }
#' 
#' @details 
#' The calculation of ENSO indices is based on \url{https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni}, consisting of SST anomalies, using a different size moving window.
#' @export
#' 
#' @examples \dontrun{ 
#' data("ERAInterim_sst_1981_2010")
#' nino <- circulationIndices(sst=ERAInterim_sst_1981_2010, index.code = "NINO3.4")
#' }


indicesENSO <- function(data, base, ref,
                       season, index.code, 
                       members=members){
  
  
  enso.index <- c("NINO3.4", "ONI")
  ind.tele <- which(enso.index %in% index.code)

  #  *** CALCULATE MONTHLY ANOMALIES *** 
  data.cen <- redim(scaleGrid(data, base, ref, time.frame = "monthly", type="center"), member = T)

  ls <- lapply(1:length(ind.tele), function(p){
 
    res <- vector("list", members)   
    count.tele <- ind.tele[p]
    if(enso.index[count.tele]=="NINO3.4") wd <- 5   else if(enso.index[count.tele]=="ONI") wd <- 3 else message("Unknown index")
    k <- as.integer(wd/2)

    # *** MEAN ANOMALIES WITH WD-MONTH MOVING WINDOW ***
    n.mon <- length(data.cen$Dates$start)
    ls <-lapply(1:n.mon, function(i){
      if(i < (k+1)) {window.mon <- c(1:(1+k))
      }else if(i > (n.mon-k)) {window.mon <- c((n.mon-k):n.mon)
      }else {window.mon <- seq((i-k),(i+k))}
      
      subMon <-subsetDimension(data.cen, dimension = "time", indices=window.mon)  
      subMon.clim <- suppressMessages(climatology(subMon)) 
      if(i != 1) {subMon.clim$Dates$start <- subMon$Dates$start[(k+1)]; subMon.clim$Dates$end <- subMon$Dates$end[(k+1)]} # assign dates of the center of the window
      if(i == 1) {subMon.clim$Dates$end <- subMon$Dates$end[1]}
      return(subMon.clim)
    }
    )
    anom <- redim(redim(bindGrid(ls, dimension = "time"), drop=TRUE), member=TRUE)

    # *** SUBSET SEASON ***
    sub <- subsetGrid(anom, season= season, drop=FALSE) # if season=NULL, it does nothing, keeps the full series
    anom.sub <- redim(suppressMessages(climatology(sub)), member=T)$Data
    attr(anom.sub, "climatology:fun") <- NULL
    
    # *** SPATIAL AVERAGE TO GET INDEX ***
    index <- aggregateGrid(anom, aggr.lat = list(FUN="mean", na.rm = TRUE), aggr.lon = list(FUN="mean", na.rm = TRUE), weight.by.lat = FALSE)
    index.sub <- subsetGrid(index, season= season, drop=FALSE)$Data
    
    # *** SAVE RESULT ***
    for(x in 1:members){
      res[[x]]$index <- index.sub[x,]
      attr(res[[1]]$index, "dimensions") <- attr(index.sub, "dimensions")[2]
      res[[x]]$pattern <- anom.sub[x,,,]
      attr(res[[1]]$pattern, "dimensions") <- attr(anom.sub, "dimensions")[3:4]
    }
    attr(res, "dates_start") <- getRefDates(sub)
    attr(res, "dates_end") <- getRefDates(sub, which = "end")
    attr(res, "season") <- season
    attr(res, "xCoords") <- sub$xyCoords$x
    attr(res, "yCoords") <- sub$xyCoords$y
    attr(res, "projection") <- attr(sub$xyCoords, "projection")
    if(members>1) names(res) <- data$Members
    message("NOTE: Calculated ", enso.index[count.tele])
    return(res)
    
  })
  names(ls) <- enso.index[ind.tele]
  return(ls)
  
}