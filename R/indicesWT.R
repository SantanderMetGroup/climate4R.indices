#     indicesWT.R Calculation of the Weather types (WT) circulation indices from grid
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

#' @title Calculation of Weather types (WT) circulation indices from grid.
#' @description Calculate circulation indices of grids or multimember grids. 
#' @param grid A grid (gridded or station dataset), or multimember grid object of geopotential height or geopotential.
#' @param season Selected month(s) for the calculation. Default: NULL (i.e. as input grid).
#' @param cluster.type Weather typing method. See details.
#' @param centers Integer value indicating the number of clusters, \strong{k}, or center points. See details.
#' @param base Baseline grid to be substracted for the calculation of anomalies. Default: NULL. See \code{?scaleGrid}.
#' @param ref Reference grid to be added for the calculation of anomalies. Default: NULL. See \code{?scaleGrid}.
#' @details The clustering parameters for weather typing are internally passed to \code{\link[transformeR]{clusterGrid}}.
#' The function calculates the weather types from the season especified as a whole. 
#' @return The WT circulation indices (and members, if applicable) with:
#' \itemize{
#' \item index: vector with the corresponding cluster from each point of the series.
#' \item pattern: Array with the spatial pattern of the clusters obtained.
#' \item dates and coordinates.
#' \item further arguments related to the WT indices, such as the corresponding cluster type, number of clusters, etc.
#' }
#' @export
#' @examples 
#' data(NCEP_slp_2001_2010)
#' wt <- indicesWT(grid=NCEP_slp_2001_2010, cluster.type="kmeans", centers = 3)



indicesWT <- function(grid, season = getSeason(grid), cluster.type, centers = NULL, base = NULL, ref = NULL) {
  
  cluster.type <- match.arg(cluster.type, choices = c("kmeans", "som", "hierarchical"))
  # *** SUBSET FROM ARGUMENT "SEASON" ***
  data.mon <- subsetGrid(grid = grid, season = season)
  message("Calculating weather types from seasons: ", paste0(season,", "))
  
  #  *** CALCULATE SEASON CENTER ANOMALIES *** 
  if (is.null(base) & is.null(ref)){
    data.cen<-data.mon
  }else {
    data.cen <- scaleGrid(grid = data.mon, base = base, ref = ref, type = "center")
  }
 
  #  *** CALCULATE CLUSTERS BY MONTH *** 
  clusters <- clusterGrid(grid = data.cen, type = cluster.type, centers = centers)
  
  #  *** PREPARE OUTPUT GRID *** 
  wt <- vector("list", 1)
  names(wt)<-cluster.type
  
  suppressMessages(members <- getShape(clusters, dimension = "member"))
  if (is.na(members)) {
    clusters<-redim(clusters)
    members <- getShape(clusters, dimension = "member")
  }
  
  wt[[1]] <- vector("list", members)
  if(members>1) names(wt[[1]]) <- paste0("Member_", 1:members)
  
  for (x in 1:members){
    memb <- vector("list", 1)
    memb[[1]]$index <- attr(clusters, "cluster")[x, ]
    memb[[1]]$pattern <- clusters$Data[x, , , ]
    attr(memb[[1]], "season") <- attr(clusters$Dates, "season")
    attr(memb[[1]], "dates_start") <- clusters$Dates$start
    attr(memb[[1]], "dates_end") <- clusters$Dates$end
    attr(memb[[1]], "centers") <-  attr(clusters, "centers")
    if (cluster.type == "kmeans") {
      attr(memb[[1]], "withinss") <- attr(clusters, "withinss")[x, ]
      attr(memb[[1]], "betweenss") <-  attr(clusters, "betweenss")[x]
    } else if (cluster.type == "hierarchical") {
      attr(memb[[1]], "height") <- attr(clusters, "height")[x, ]
      attr(memb[[1]], "cutree.at.height") <- attr(clusters, "cutree.at.height")[x]
      attr(memb[[1]], "diff.height.threshold") <- attr(clusters, "diff.height.threshold")
    } 
    wt[[1]][[x]] <- memb
  }
  
  attr(wt, "xCoords") <- clusters$xyCoords$x
  attr(wt, "yCoords") <- clusters$xyCoords$y
  attr(wt, "projection") <- attr(clusters$xyCoords, "projection")
  
  return(wt)
}
