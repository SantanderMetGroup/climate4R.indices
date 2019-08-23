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
#' @details The clustering parameters for weather typing are internally passed to \code{\link[transformeR]{clusterGrid}}. 
#' @return The WT circulation indices (and members, if applicable) with:
#' \itemize{
#' \item index: vector with the corresponding cluster from each point of the series.
#' \item pattern: Array with the spatial pattern of the clusters obtained.
#' \item dates and coordinates.
#' \item further arguments related to the WT indices, such as the corresponding cluster type, number of clusters, etc.
#' }
#' @export
#' @examples 
#' data(NCEP_hgt500_2001_2010)
#' wt <- indicesWT(grid=NCEP_hgt500_2001_2010, season=1,
#'                 cluster.type="kmeans", centers = 3)



indicesWT <- function(grid, season, cluster.type, centers = NULL
                      #rot=rot, members=members
                      ) {
  cluster.type <- match.arg(cluster.type, choices = c("kmeans", "som", "hierarchical", "lamb"))
  # *** SUBSET WITH MONTHS FROM ARGUMENT "SEASON" ***
  data.mon <- subsetGrid(grid = grid, season = season)
  #  *** CALCULATE CLUSTERS BY MONTH*** 
  # depending on the index.code call to a different clustering type.
  # argument k needs to be provided to circIndexGrid()
  #wt <- transformeR::clusterGrid(data.mon, type = cluster.type, centers=centers)
  wt <- clusterGrid(data.mon, type = cluster.type, centers = centers)
  wt$index <- attr(wt, "clusters")
  attr(wt, "clusters") <- NULL
  wt$pattern <- wt$Data 
  wt$Data <- NULL
  return(wt)
}
