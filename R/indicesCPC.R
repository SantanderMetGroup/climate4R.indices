#     indicesCPC.R Calculation of the CPC clirculation indices from grid
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

#' @title Calculation of CPC circulation indices from grid.
#' @description Calculate circulation indices of grids or multimember grids. 
#' @param grid A grid (gridded or station dataset), or multimember grid object of geopotential height or geopotential.
#' @param index.code Circulation index (or vector of indices with same input variables) to be computed. See \code{circIndexShow()} for details.
#' @param base Baseline grid to be substracted for the calculation of anomalies. Default: NULL. See \code{?scaleGrid}.
#' @param ref Reference grid to be added for the calculation of anomalies. Default: NULL. See \code{?scaleGrid}.
#' @param season Selected month(s) for the calculation. Default: NULL (i.e. as input grid).
#' @param n.pcs Integer vector with the number of EOFs to be retained for the CPC indices. Default: 10. See details.
#' @param rot Logical. Should VARIMAX-Rotation be performed? Default: TRUE. This argument is only relevant for CPC indices. See details.
#' @param members Select number of members.
#' @param match Character string with the criterion to be used for the detection of CPC indices. Options "spatial" or "temporal". Default: "spatial". See details.

#' @return A list of circulation indices (and members, if applicable) with:
#' \itemize{
#' * index: vector with the time series of the teleconnection index.
#' * pattern: matrix with the spatial pattern of the teleconnection.
#' * dates and coordinates as list attributes.
#' * further arguments related to the CPC indices, such as the corresponding (r)EOF and (temporal or spatial, depending on \code{'match'}) correlation with the original index.
#' }
#' 
#' @details 
#' CPC indices are obtained, by default, as the first 10 Varimax-rotated EOFs, as explained in \url{https://www.cpc.ncep.noaa.gov/data/teledoc/telecontents.shtml}. The core of this function is \code{stats::prcomp} including Varimax rotation.
#' The rotated EOFs are obtained from the monthly standardized anomalies of geopotential or geopotential height at 500hPa, with a 3-month moving window.
#' The argument \code{match} is used to assign each rEOF to a circulation index and pattern. Matching is based on 'temporal' or 'spatial' correction of the CPC original (NCEP Reanalysis-based) indices.

#' @export
#' @importFrom stats cor 
#' @importFrom utils data
#' @examples \dontrun{ 
#' data(NCEP_hgt500_2001_2010)
#' cpc <- indicesCPC(grid=NCEP_hgt500_2001_2010, index.code = c("NAO", "EA","PNA"), season=1)
#' }


indicesCPC <- function(grid, base, ref,
                       season, index.code, 
                       match=match, n.pcs=n.pcs, rot=rot, 
                       members=members){

    cpc <- NULL 
    cpc.index <- c("NAO", "EA", "WP", "EP/NP", "PNA", "EA/WR", "SCA", "TNH", "POL", "PT")
    ind.tele <- which(cpc.index %in% index.code)
    years <- unique(getYearsAsINDEX(grid))

    # *** READ CPC TELECONNECTION INDICES *** 
    # Downloaded from wget ftp://ftp.cpc.ncep.noaa.gov/wd52dg/data/indices/tele_index.nh
    tele <- read.tele()
    ntele <- ncol(tele)-3
    
    #  *** CALCULATE MONTHLY STANDARDIZE ANOMALIES *** 
    data.cen <- redim(scaleGrid(grid, base, ref, time.frame = "monthly", type="standardize"), member = T) 
    
    # *** SUBSET MONTHLY WITH 3-MONTH MOVING WINDOW FOR PCA CALCULATION ***
    pca <- lapply(getSeason(data.cen), function(mon){
      if(mon==1){
        window.mon <- c(12,1,2)
      } else if(mon==12){
        window.mon <- c(11,12,1)
      } else {window.mon <- seq((mon-1),(mon+1))}
      data.mon <- redim(subsetGrid(data.cen, season = window.mon), member = T)
      
      # *** PERFORM rEOFs ***
      pca1 <- prinComp(data.mon, n.eofs = n.pcs, keep.orig = TRUE, rot=rot)
    })

    # *** MATCH CPC AND CALCULATED PATTERNS WITH DIFFERENT CRITERIA ***
    if(match=="spatial"){
      
      # *** LOAD CPC PATTERNS ***
      load(system.file("tele_pattern.Rdata", package = "climate4R.indices")) # loaded cpc
      
      # *** SPATIAL CORRELATION BETWEEN CALCULATED EOFs AND CPC PATTERNS ***
      ls <- lapply(1:length(ind.tele), function(p){

        count.tele <- ind.tele[p]
        cpc.interp <- suppressMessages(interpGrid(cpc[[count.tele]], new.coordinates =  list(x=data.cen$xyCoords$x, y=data.cen$xyCoords$y), method = "bilinear"))
        browser()
        res <- vector("list", members)
        for(x in 1:members){
          for(mon in season){
            
            ind.mon <- seq(mon,dim(pca[[mon]][[1]][[x]]$PCs)[1],3)
            cpc.redim <- array3Dto2Dmat(cpc.interp$Data)[mon,] 
            corr.pattern <- cor(pca[[mon]][[1]][[x]]$EOFs, cpc.redim)
          
            # Select the EOF that maximizes spatial corr (and take sign)
            idx <- which.max(abs(corr.pattern))
            res[[x]]$index <- pca[[mon]][[1]][[x]]$PCs[ind.mon, idx] * sign(corr.pattern[idx])
            attr(res[[x]]$index, "dimensions") <- "time"
            patt <- pca[[mon]][[1]][[x]]$EOFs[,idx] * sign(corr.pattern[idx])
            res[[x]]$pattern <- matrix(patt, nrow = length(data.cen$xyCoords$y), ncol = length(data.cen$xyCoords$x), byrow = FALSE) #redim
            attr(res[[x]]$pattern, "dimensions") <- c("lat","lon")
            attr(res[[x]], "ind.eof") <- idx
            attr(res[[x]], "sign.eof") <- sign(corr.pattern[idx])
            attr(res[[x]], "corr.eof") <- corr.pattern[idx]
            attr(res[[x]], "match") <- match
          }
        }
        attr(res, "season") <- season
        attr(res, "dates_start") <- attr(pca[[mon]], "dates_start")[ind.mon]
        attr(res, "dates_end") <- attr(pca[[mon]], "dates_end")[ind.mon]
        attr(res, "xCoords") <- data.cen$xyCoords$x
        attr(res, "yCoords") <- data.cen$xyCoords$y
        attr(res, "projection") <- attr(data.cen$xyCoords, "projection")
        if(members>1) names(res) <- data$Members
        message("NOTE: Calculated ", cpc.index[count.tele])
        return(res)
      })
      names(ls) <- cpc.index[ind.tele]
      
    } else if(match=="temporal"){

        ls <- lapply(1:length(ind.tele), function(p){
          
          count.tele <- ind.tele[p]
          res <- vector("list", members)
          for(x in 1:members){
            for(mon in season){
            
              # *** TEMPORAL CORRELATION BETWEEN CALCULATED PCs AND CPC TELECONNECTIONS ***
              ind.mon <- seq(mon,dim(pca[[mon]][[1]][[x]]$PCs)[1],3)
              ind.mon.t <- which(tele$Month==mon & tele$Year>=years[1] & tele$Year <=years[length(years)])
              corr.index <- cor(pca[[mon]][[1]][[x]]$PCs[ind.mon,], tele[[(count.tele+2)]][ind.mon.t])
          
              # Select the PC that maximizes temporal corr (and take sign)
              idx <- which.max(abs(corr.index))
              res[[x]]$index <- pca[[mon]][[1]][[x]]$PCs[ind.mon, idx] * sign(corr.index[idx])
              attr(res[[x]]$index, "dimensions") <- "time"
              patt <- pca[[mon]][[1]][[x]]$EOFs[,idx] * sign(corr.index[idx])
              res[[x]]$pattern <- matrix(patt, nrow = length(data.cen$xyCoords$y), ncol = length(data.cen$xyCoords$x), byrow = FALSE) 
              attr(res[[x]]$pattern, "dimensions") <- c("lat","lon")
              attr(res[[x]], "ind.eof") <- idx
              attr(res[[x]], "sign.eof") <- sign(corr.index[idx])
              attr(res[[x]], "corr.eof") <- corr.index[idx]
              attr(res[[x]], "match") <- match
            } 
          }
          attr(res, "season") <- season
          attr(res, "dates_start") <- attr(pca[[mon]], "dates_start")[ind.mon]
          attr(res, "dates_end") <- attr(pca[[mon]], "dates_end")[ind.mon]
          attr(res, "xCoords") <- data.cen$xyCoords$x
          attr(res, "yCoords") <- data.cen$xyCoords$y
          attr(res, "projection") <- attr(data.cen$xyCoords, "projection")
          if(members>1) names(res) <- data$Members
          message("NOTE: Calculated ", cpc.index[count.tele])
          return(res)
      })
      names(ls) <- cpc.index[ind.tele]
      
    } else{stop("Error: Unknown match criterion for rEOFs and CPC matching")}

  
    return(ls)
}