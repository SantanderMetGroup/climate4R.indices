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
#' @param zg A grid or multimember grid object of geopotential height.
#' @param z A grid or multimember grid object of geopotential.
#' @param sst A grid or multimember grid object of sea surface temperature.
#' @param psl A grid or multimember grid object of sea level pressure.
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
#' \code{\link{circIndexShow}} displays on screen the full list of available circulation indices and their codes.
#' Several indices can be calculated at the same time, as long as they depend on the same input variable(s) and spatial domain. All indices are calculated on a monthly basis. Therefore a temporal aggregation is performed if input data is daily.
#' Results for the desired months in \code{season} are provided, but it is recommended to have full series as input, since many indices use a moving window for the calculation.
#' 
#' \strong{CPC indices}
#' 
#' Either \code{z} or \code{zg} are valid input variables. CPC indices are obtained, by default, as the first 10 Varimax-rotated EOFs, as explained in \url{https://www.cpc.ncep.noaa.gov/data/teledoc/telecontents.shtml}. The core of this function is \code{stats::prcomp} including Varimax rotation.
#' The rotated EOFs are obtained from the monthly standardized anomalies of geopotential or geopotential height at 500hPa, with a 3-month moving window.
#' The argument \code{match} is used to assign each rEOF to a circulation index and pattern. Matching is based on 'temporal' or 'spatial' correction of the CPC original (NCEP Reanalysis-based) indices.
#' 
#' \strong{ENSO indices}
#' 
#' The calculation of ENSO indices is based on \url{https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni}, consisting of SST anomalies, using a moving window of different size for each index.
#' 
#' @export
#' @importFrom utils data
#' 
#' @author A. Casanueva
#' @examples 
#' data(NCEP_hgt500_2001_2010)
#' cpc <- circIndexGrid(zg=NCEP_hgt500_2001_2010, index.code = c("NAO", "EA","PNA"), season=1)
#' data(ERAInterim_sst_1981_2010)
#' nino <- circIndexGrid(sst=ERAInterim_sst_1981_2010, index.code = "NINO3.4")
#' 
#' 
circIndexGrid <- function(zg=NULL,
                           z=NULL, 
                           sst=NULL,
                           psl=NULL,
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
  if(!any(match(c(cpc.index,enso.index) ,index.code, nomatch=FALSE))) stop("Non valid index selected: Use circIndexShow() to select an index.", call. = FALSE)  
  if(any(match(cpc.index,index.code, nomatch=FALSE)) & any(match(enso.index,index.code, nomatch=FALSE))) stop("Indices require different input variables and spatial domains, see circIndexShow()", call. = FALSE)  
  

  # *** CHECK VARIABLE NAME AND DOMAIN ***
  if(any(match(cpc.index,index.code, nomatch = FALSE))){
    # Variables, units and levels
    if(!is.null(zg)){grid <- zg
    } else if(!is.null(z)){ if(transformeR::getGridUnits(z)=="m2.s-2") grid <- gridArithmetics(z, 9.8, operator="/") else stop("Check units of z" , call. = FALSE)
    } else {stop("z or zg arguments needed", call. = FALSE)}
    if(transformeR::getGridVerticalLevels(grid)!=500) stop("Required vertical level is 500hPa", call. = FALSE)
    
    # Domain 
    lat <- c(20,90)
    lon <- c(-180, 180)
  } 
  
  if(any(match(enso.index,index.code, nomatch = FALSE))){
    
    # Variables
    if(is.null(sst)) stop("sst argument needed", call. = FALSE)
    if(!is.null(sst)){ if(getVarNames(sst)=="SST") grid <- sst else stop("Check variable name", call. = FALSE)}
    
    # Domain 
    lat <- c(-5,5)
    lon <- c(-170,-120)
   } 
  
  # *** SET DOMAIN ***
  if(min(grid$xyCoords$x) > lon[1] | max(grid$xyCoords$x) < lon[2]) stop("Longitudes do not cover the domain", call. = FALSE)
  if(min(grid$xyCoords$y) > lat[1] | max(grid$xyCoords$y) < lat[2]) stop("Latitudes do not cover the domain", call. = FALSE)
  grid <- subsetGrid(grid, latLim =  lat, lonLim=lon)

  # *** REQUIRE MONTHLY DATA ***
  if(getTimeResolution(grid)!="MM" | length(getSeason(grid))!=12) stop("Monthly data for full years are needed", call. = FALSE)
  
  # *** CHECK MEMBER DIMENSION ***
  if(is.null(members)) members <- getShape(redim(grid, member=T), "member") else  members <- as.integer(members)
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
  
    aux1 <- indicesCPC(grid=data.w, season=season, 
                      index.code=index.code, base=base, ref=ref,
                      match=match, n.pcs=n.pcs, rot=rot, 
                      members= members)
      
  
  } 
  if(any(match(enso.index,index.code, nomatch = FALSE))){
    
    aux2 <- indicesENSO(grid=data.w, season=season, 
                      index.code=index.code, 
                      base=base, ref=ref,
                      members= members)
    
    
  } 
  
  index <- c(aux1,aux2)
  attr(index, "class") <- "circulationIndex"
  return(index)
}


#' @title List all available circulation indices
#' @description Print a table with a summary of the available circulation indices
#' @return Print a table on the screen with the following columns:
#' \itemize{
#' \item \strong{code}: Code of the index. This is the character string used as input value
#' for the argument \code{index.code} in \code{\link{circIndexGrid}}
#' \item \strong{longname}: Long description of the circulation index.
#' \item \strong{fun}: The name of the function used to calculate it.
#' \item \strong{zg, hgt, sst, slp}: A logical value (0/1) indicating the input variables required for index calculation
#' \item \strong{reference}: Reference for the implemented calculation.
#' }
#' @author J. Bedia, M. Iturbide, A. Casanueva
#' @export

circIndexShow <- function() {
  read.masterCirc()
}

#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom utils read.table

read.masterCirc <- function() {
  system.file("master_circulation", package = "climate4R.indices") %>% read.table(header = TRUE,
                                                                                  sep = ";",
                                                                                  stringsAsFactors = FALSE,
                                                                                  na.strings = "")
}

#' @title Read CPC teleconnections
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom utils read.csv
#' @details Downloaded from wget ftp://ftp.cpc.ncep.noaa.gov/wd52dg/data/indices/tele_index.nh
#' These are monthly tabulated indices since 1950 to present.  Indices are standardized by the 1981-2010 climatology.

read.tele <- function() {
  names <- c("Year","Month", "NAO","EA","WP","EP/NP","PNA","EA/WR","SCA","TNH","POL","PT","Expl.Var.")
  system.file("tele_index.nh", package = "climate4R.indices") %>% read.csv(skip=17,strip.white=T,sep='', 
                                                                           na.strings ="-99.90", stringsAsFactors = FALSE,
                                                                           header=TRUE, col.names = names)
}