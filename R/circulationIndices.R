#     circulationIndices.R Calculation of clirculation indices of grid data
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

#' @title Calculation of circulation indices of Grids
#' @description Calculate circulation indices of grids or multimember grids. 
#' @param hgt A grid or multimember grid object of geopotential height.
#' @param zg A grid or multimember grid object of geopotential.
#' @param sst A grid or multimember grid object of sea surface temperature.
#' @param slp A grid or multimember grid object of sea level pressure.
#' @inheritParams indicesCPC
#' @inheritParams indicesENSO
#' 
#' @return A list of circulation indices (and members, if applicable) with:
#' \itemize{
#' \item index: vector with the time series of the teleconnection index.
#' \item pattern: matrix with the spatial pattern of the teleconnection.
#' \item dates and coordinates as list attributes.
#' \item further arguments related to the CPC indices, such as the corresponding (r)EOF and (temporal or spatial, depending on \code{'match'}) correlation with the original index.
#' }
#' 
#' @details
#' \code{\link{indexCircShow()}} displays on screen the full list of available circulation indices and their codes.
#' Several indices can be calculated at the same time, as long as they depend on the same input variable(s) and spatial domain. All indices are calculated on a monthly basis. Therefore a temporal aggregation is performed if input data is daily.
#' Results for the desired months in \code{season} are provided, but it is recommended to have full series as input, since many indices use a moving window for the calculation.
#' 
#' \strong{CPC indices}
#' 
#' Either \code{zg} or \code{hgt} are valid input variables. CPC indices are obtained, by default, as the first 10 Varimax-rotated EOFs, as explained in \url{https://www.cpc.ncep.noaa.gov/data/teledoc/telecontents.shtml}. The core of this function is \code{stats::prcomp} including Varimax rotation.
#' The rotated EOFs are obtained from the monthly standardized anomalies of geopotential or geopotential height at 500hPa, with a 3-month moving window.
#' The argument \code{match} is used to assign each rEOF to a circulation index and pattern. Matching is based on 'temporal' or 'spatial' correction of the CPC original (NCEP Reanalysis-based) indices.
#' 
#' \strong{ENSO indices}
#' 
#' The calculation of ENSO indices is based on \url{https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni}, consisting of SST anomalies, using a moving window of different size for each index.
#' 
#' @export
#' @author A. Casanueva
#' @examples \dontrun{ 
#' data("NCEP_hgt500_1981_2010")
#' cpc <- circulationIndices(hgt=NCEP_hgt500_1981_2010, index.code = c("NAO", "EA","PNA"))
#' data("ERAInterim_sst_1981_2010")
#' nino <- circulationIndices(sst=ERAInterim_sst_1981_2010, index.code = "NINO3.4")
#' }
#' 
circulationIndices <- function(zg=NULL,
                               hgt=NULL, 
                               sst=NULL,
                               slp=NULL,
                               index.code, 
                               season=NULL, 
                               base=NULL, ref=NULL,
                               match="spatial", n.pcs=10, rot=TRUE, 
                               members=NULL){

  if(length(season)>1) stop("More than one season is not available yet", call. = FALSE)  
  
  # *** CHECK INDICES NAME *** # include here any implemented index
  index.code <- toupper(index.code)
  cpc.index <- c("NAO", "EA", "WP", "EP/NP", "PNA", "EA/WR", "SCA", "TNH", "POL", "PT")
  enso.index <- c("NINO3.4", "ONI")
  if(!any(match(c(cpc.index,enso.index) ,index.code, nomatch=FALSE))) stop("Check index name", call. = FALSE)  
  if(any(match(cpc.index,index.code, nomatch=FALSE)) & any(match(enso.index,index.code, nomatch=FALSE))) stop("Indices require different input variables and spatial domains, see ?indexCircShow()", call. = FALSE)  
  

  # *** CHECK VARIABLE NAME AND DOMAIN ***
  if(any(match(cpc.index,index.code, nomatch = FALSE))){

    # Variables
    assertthat::assert_that((!is.null(zg) | !is.null(hgt)), msg = "zg or hgt arguments needed")
    
    if(!is.null(zg)){ if(getVarNames(zg)=="z@500") grid <- gridArithmetics(zg, 9.8, operator="/") else stop("Check variable name", call. = FALSE)
    } else if(!is.null(hgt)){ if(getVarNames(hgt)=="hgt@500") grid <- hgt else stop("Check variable name", call. = FALSE)
    }

    # Domain 
    lat <- c(20,90)
    lon <- c(-180, 180)
  } 
  
  if(any(match(enso.index,index.code, nomatch = FALSE))){
    
    # Variables
    assertthat::assert_that(!is.null(sst), msg = "sst argument needed")
    if(!is.null(sst)){ if(getVarNames(sst)=="SST") grid <- sst else stop("Check variable name", call. = FALSE)}
    
    # Domain 
    lat <- c(-5,5)
    lon <- c(-170,-120)
   } 
  
  # *** SET DOMAIN ***
  assertthat::assert_that((min(grid$xyCoords$x) <= lon[1] | max(grid$xyCoords$x) >= lon[2]), msg = "Longitudes do not cover the domain")
  assertthat::assert_that((min(grid$xyCoords$y) <= lat[1] | max(grid$xyCoords$y) >= lat[2]), msg = "Latitudes do not cover the domain")
  grid <- subsetGrid(grid, latLim =  lat, lonLim=lon)

  # *** REQUIRE MONTHLY DATA (if not provided) ***
  if(getTimeResolution(grid)!="MM"){
    grid <- aggregateGrid(grid, aggr.m = list(FUN = "mean", na.rm = TRUE))
  }

  # *** CHECK MEMBER DIMENSION ***
  if(is.null(members)) members <- getShape(redim(grid, member=T), "member") else  members <- as.integer(member)
  if (!(members %in% 1:getShape(redim(grid, member=T), "member"))) {
    stop("'members' value must be between 1 and the total number of members (", getShape(redim(grid, member=T), "member"), " in this case)", call. = FALSE)
  }

  # *** INITIALIZE RESULT ***
  aux1 <- NULL
  aux2 <- NULL

  # *** REDIM GRID *** (in order to hace same dims even without members)
  data.redim <- redim(grid, member = T) 
  
  # *** GRIDBOXES WEIGTHING ***
  latr <- grid$xyCoords$y*pi/180 # latitude to radians
  coslat <- sqrt(cos(latr))
  
  if(any(getDim(data.redim) %in% "lat")){
    
    nlat <- length(coslat)
    ls <-lapply(1:nlat, function(i){
      lati <-subsetDimension(data.redim, dimension = "lat", indices=i)  
      data.weight <- gridArithmetics(lati, coslat[i], operator = "*")
    }
    )
    data.w <- redim(redim(bindGrid(ls, dimension = "lat"), drop=TRUE), member=TRUE)
  }else{
    stop("Variable lat not found", call. = FALSE)
  }

  # *** START INDICES CALCULATION ***
  if(any(match(cpc.index,index.code, nomatch=FALSE))){
  
    aux1 <- indicesCPC(data=data.w, season=season, 
                      index.code=index.code, base=base, ref=ref,
                      match=match, n.pcs=n.pcs, rot=rot, 
                      members= members)
      
  
  } 
  if(any(match(enso.index,index.code, nomatch = FALSE))){
    
    aux2 <- indicesENSO(data=data.w, season=season, 
                      index.code=index.code, 
                      base=base, ref=ref,
                      members= members)
    
    
  } 
  
  index <- c(aux1,aux2)
  attr(index, "class") <- "circulationIndex"
  return(index)
}
