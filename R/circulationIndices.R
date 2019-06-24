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
#' @param hgt A grid (gridded or station dataset), or multimember grid object of geopotential height.
#' @param zg A grid (gridded or station dataset), or multimember grid object of geopotential height.
#' @param sst
#' @param slp
#' @param index.code Circulation index to be computed. Implemented indices: NAO, EA, WP, EP/NP, PNA, EA/WR, SCA, TNH, POL, PT. See https://www.cpc.ncep.noaa.gov/data/teledoc/telecontents.shtml for details.
#' @param season Month or season desired for the calculation.
#' @param n.pcs Integer vector. Number of EOFs to be retained. Default to 10 for the CPC indices. See details.
#' @param rot Default: TRUE.
#' @param members
#' @param match Options "spatial", "temporal". Default: "temporal".
#' @param v.exp Maximum fraction of explained variance, in the range (0,1]. Used to determine the number of EOFs 
#' to be retained, as an alternative to \code{n.eofs}. Default to \code{NULL}. See details.

#' @return A list (of length 1 for single grids and length Nmem for multimember grid) with:
#'  \itemize{
#'  \item \code{PCs}: A matrix of principal components, arranged in columns by decreasing importance order 
#'  \item \code{EOFs}: A matrix of EOFs, arranged in columns by decreasing importance order
#'  \item \code{orig}: Either the original variable in the form of a 2D-matrix (when \code{keep.orig = TRUE}),
#'  or \code{NA} when \code{keep.origin = FALSE} (the default). In both cases, the parameters used for input data standardization
#'  (mean and standard deviation) are returned as attributes of this component (see the examples).
#'  }

#' @note 
#' 
#' @details
#' CPC indices are obtained as the first 10 Varimax-rotated EOFs, as explained in \url{https://www.cpc.ncep.noaa.gov/data/teledoc/telecontents.shtml}. The core of this function is \code{stats::prcomp} including Varimax rotation.
#' The rotated EOFs are obtained from the monthly geopotential height standardize anomalies, with a 3-month moving window.
#' 
#' @export
#' @references
#' @author A. Casanueva
#' @examples 
#' 
#' 
circulationIndices <- function(zg=NULL,
                               hgt=NULL, 
                               sst=NULL,
                               slp=NULL,
                               season, 
                               index.code, match="spatial", n.pcs=10, rot=TRUE, 
                               members=NULL){

    
  # *** CHECK INDICES NAME *** # include here any implemented index
  cpc.index <- c("NAO", "EA", "WP", "EP/NP", "PNA", "EA/WR", "SCA", "TNH", "POL", "PT")
  if(!any(match(c(cpc.index,"ENSO") ,index.code, nomatch=FALSE))) stop(message("Check index name"))  
  
  # *** CHECK VARIABLE NAME AND DOMAIN ***
  if(any(match(cpc.index,index.code, nomatch = FALSE))){

    # Variables
    assertthat::assert_that((!is.null(zg) | !is.null(hgt)), msg = "zg or hgt arguments needed")
    if(!is.null(zg)){ if(getVarNames(zg)=="z@500") grid <- zg else stop("Error: Check variable name")
    } else if(!is.null(hgt)){ if(getVarNames(hgt)=="hgt@500") grid <- hgt else stop("Error: Check variable name")
    }
    
    # Set domain
    assertthat::assert_that((min(grid$xyCoords$x) < -175 | max(grid$xyCoords$x) > 175), msg = "Longitudes do not cover the domain")
    assertthat::assert_that((min(grid$xyCoords$y) < 20 | max(grid$xyCoords$y) > 88), msg = "Latitudes do not cover the domain")
    grid <- subsetGrid(grid, latLim =  c(20,90))
  } 
  
  if(any(match("ENSO",index.code, nomatch = FALSE))){
    
    # Variables
    assertthat::assert_that(!is.null(sst), msg = "sst argument needed")
    if(!is.null(sst)){ if(getVarNames(sst)=="SST") grid <- sst else stop("Error: Check variable name")}
    
    # Set domain
    #assertthat::assert_that((min(grid$xyCoords$x) < -175 | max(grid$xyCoords$x) > 175), msg = "Longitudes do not cover the domain")
    #assertthat::assert_that((min(grid$xyCoords$y) < 20 | max(grid$xyCoords$y) > 88), msg = "Latitudes do not cover the domain")
    #grid <- subsetGrid(grid, latLim =  c(20,90))
  } 
  
  # *** CHECK MEMBER DIMENSION ***
  if(is.null(members)) members <- getShape(redim(grid, member=T), "member") else  members <- as.integer(member)
  if (!(members %in% 1:getShape(redim(grid, member=T), "member"))) {
    stop("'members' value must be between 1 and the total number of members (", getShape(redim(grid, member=T), "member"), " in this case)", call. = FALSE)
  }

  # *** REDIM GRID *** (in order to hace same dims even without members)
  data.redim <- redim(data, member = T)
  
  # *** GRIDBOXES WEIGTHING ***
  latr <- data$xyCoords$y*pi/180 # latitude to radians
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
    stop(message("Error: variable lat not found"))
  }
  
  # *** START INDICES CALCULATION ***
  if(any(match(cpc.index,index.code, nomatch=FALSE))){
  
    ind.tele <- which(cpc.index %in% index.code)
    years <- unique(getYearsAsINDEX(grid))
    
    # *** READ CPC TELECONNECTION INDICES *** 
    # Downloaded from wget ftp://ftp.cpc.ncep.noaa.gov/wd52dg/data/indices/tele_index.nh
    # These are monthly tabulated indices since 1950 to present.  Indices are standardized by the 1981-2010 climatology)
    names <- c("Year","Month", "NAO","EA","WP","EP/NP","PNA","EA/WR","SCA","TNH","POL","PT","Expl.Var.")
    tele <- read.csv(file=paste0(workdir, "/DATA/tele_index.nh"), skip=17,strip.white=T,sep='', na.strings ="-99.90",
                     stringsAsFactors = FALSE, header=TRUE, col.names = names)
    ntele <- ncol(tele)-3
    
    
    # *** SUBSET MONTHLY WITH 3-MONTH MOVING WINDOW ***
    # 3-month moving window per month. without the moving window results do not agree with CPC.
    #for(mon in 1:12){
      mon <- 1
      if(mon==1) window.mon <- c(12,1,2) else window.mon <- seq((mon-1),(mon+1))
      data.mon <- redim(subsetGrid(data.w, season = window.mon), member = T)
      n.window <- length(data.mon$Dates$start)
      ind.mon <- seq(mon,n.window ,length(window.mon))
      
      # Filter month and years from CPC
      ind.mon.t <- which(tele$Month==mon & tele$Year>=years[1] & tele$Year <=years[length(years)])
      assertthat::assert_that(length(ind.mon)==length(ind.mon.t), msg = "Dates for CPC teleconnections and grid do not match")
    
      #  *** CALCULATE MONTHLY STANDARDIZE ANOMALIES *** 
      data.cen <- redim(scaleGrid(data.mon, time.frame = "monthly", type="standardize"), member = T) # from this step, it is equivalent to use hgt or z
      
      # *** PERFORM rEOFs ***
      pca <- prinComp(data.cen, n.eofs = n.pcs, keep.orig = TRUE, rot=rot)
     
      # *** MATCH CPC AND CALCULATED PATTERNS WITH DIFFERENT CRITERIA ***
      res <- vector("list",length(ind.tele))
      if(match=="spatial"){
        ls <- lapply(1:members, function(x){
          
        # *** CALCULATE CPC PATTERNS (projecting the indices in the geopotential height fields) ***
        ngrids <- dim(pca[[1]][[x]]$orig)[2] 
        pattern.cpc <- matrix(NA, nrow=ngrids, ncol=n.pcs)
        for(p in 1:ntele){
          pattern.cpc[,p] <- cor(pca[[1]][[x]]$orig[ind.mon,], tele[[p+2]][ind.mon.t]) 
        }
        
        # *** SPATIAL CORRELATION BETWEEN CALCULATED EOFs AND CPC TELECONNECTIONS ***
        for(p in 1:length(ind.tele)){
    
          count.tele <- ind.tele[p]
          cpc1 <- pattern.cpc[,count.tele]
          corr.pattern <- cor(pca[[1]][[x]]$EOFs, cpc1)
 
          # Select the EOF that maximizes spatial corr (and take sign)
          idx <- which.max(abs(corr.pattern))
          res[[p]]$index <- pca[[1]][[x]]$PCs[ind.mon, idx] * sign(corr.pattern[idx])
          res[[p]]$pattern <- pca[[1]][[x]]$EOFs[,idx] * sign(corr.pattern[idx])
          res[[p]]$ind.eof <- idx
          res[[p]]$sign.eof <- sign(corr.pattern[idx])
          res[[p]]$corr.eof <- corr.pattern[idx]
          message("NOTE: Calculated ", cpc.index[count.tele])
        }
        names(res) <- cpc.index[ind.tele]
        return(res)
        })
      
     
      } else if(match=="temporal"){
    
        ls <- lapply(1:members, function(x){
        # *** TEMPORAL CORRELATION BETWEEN CALCULATED PCs AND CPC TELECONNECTIONS ***
        for(p in 1:length(ind.tele)){
          
          count.tele <- ind.tele[p]
          corr.index <- cor(pca[[1]][[x]]$PCs[ind.mon,], tele[[(count.tele+2)]][ind.mon.t])

          # Select the PC that maximizes temporal corr (and take sign)
          idx <- which.max(abs(corr.index))
          res[[p]]$index <- pca[[1]][[x]]$PCs[ind.mon, idx] * sign(corr.index[idx])
          res[[p]]$pattern <- pca[[1]][[x]]$EOFs[,idx] * sign(corr.index[idx])
          res[[p]]$ind.eof <- idx
          res[[p]]$sign.eof <- sign(corr.index[idx])
          res[[p]]$corr.eof <- corr.index[idx]
          message("NOTE: Calculated ", cpc.index[count.tele])
        } 
        names(res) <- cpc.index[ind.tele]
        return(res)
        })
       
      } else{stop("Error: Unknown match criterion for rEOFs and CPC matching")}
      if(members>1) names(ls) <- grid$Members
    
       # resize pattern. create grid object? 
      # plotting index and eof in stereographic, include in transformeR.
    
    #} # end loop over months
    
  
  } else if(match("ENSO",index.code, nomatch=FALSE)){ 
    browser()
      # *** SUBSET MONTHLY WITH 3-MONTH MOVING WINDOW ***
      # 5-month moving window per month.
      mon <- 1
      if(mon==1) window.mon <- c(11,12,1,2,3) else window.mon <- seq((mon-2),(mon+2)) 
      data.mon <- redim(subsetGrid(data.w, season = window.mon), member = T)
      n.window <- length(data.mon$Dates$start)
      ind.mon <- seq(mon,n.window ,length(window.mon))
    
    #  *** CALCULATE MONTHLY STANDARDIZE ANOMALIES *** 
    data.cen <- redim(scaleGrid(data.mon, time.frame = "monthly", type="standardize"), member = T) # from this step, it is equivalent to use hgt or z
    
    #pca <- prinComp(data.mon, n.eofs = n.pcs, keep.orig = TRUE) # da errores por tener puntos siempre con cts o cero
    # añadir class al resultado.
    
  } else {stop(message("index not implemented yet"))}
  
  return(ls)
}
