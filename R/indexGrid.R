#     indexGrid.R Climate Indices in Climate4R
#
#     Copyright (C) 2018 Santander Meteorology Group (http://www.meteo.unican.es)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program. If not, see <http://www.gnu.org/licenses/>.

#' @title Climate Indices in Climate4R
#' @description Calculation of indices.
#' @param tn A climate4R dataset of daily minimum temperature (degrees C)
#' @param tx A climate4R dataset of daily maximum temperature (degrees C)
#' @param tm A climate4R dataset of daily maximum temperature (degrees C)
#' @param pr A climate4R dataset of daily precipitation (mm)
#' @param baseline Optional climate4R dataset. Only used if \code{index.code = "P"}, for calculating the relevant percentiles.
#' @param index.code Character string, indicating the specific code of the index (see Details).
#' @param time.resolution Output time resolution. Choices are "month", "year" (default) and "climatology".
#' @param ... Optional. A list of arguments internally passed to the functions displayed by \code{\link{indexShow}}.
#' @template templateParallelParams
#' @import transformeR
#' @importFrom parallel stopCluster
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom utils head
#' @details \code{\link{indexShow}} will display on screen the full list of available indices and their codes.
#' The names of the internal functions calculating each index is also displayed, whose help files can aid in
#' the definition of index-specific arguments.
#' 
#' @template templateParallel
#' @author M. Iturbide
#' @export
#' @examples 
#' require(transformeR)
#' fd <- indexGrid(tn = EOBS_Iberia_tas,
#'                 time.resolution = "year",
#'                 index.code = "FD")
#' per1 <- indexGrid(tn = EOBS_Iberia_tas,
#'                  time.resolution = "year",
#'                  index.code = "P",
#'                  percent = 90)
#' per2 <- indexGrid(tn = CFS_Iberia_tas,
#'                 time.resolution = "year",
#'                 index.code = "P",
#'                 baseline = CFS_Iberia_tas,
#'                 percent = 90)
#' hdd <- indexGrid(tn = CFS_Iberia_tas,
#'                  tx = CFS_Iberia_tas,
#'                  tm = CFS_Iberia_tas,
#'                  time.resolution = "year",
#'                  index.code = "HDD")



indexGrid <- function(tn = NULL,
                      tx = NULL,
                      tm = NULL,
                      pr = NULL,
                      baseline = NULL,
                      index.code,
                      time.resolution = "year",
                      ...,
                      parallel = FALSE,
                      max.ncores = 16,
                      ncores = NULL) {
  index.arg.list <- list(...)
  choices <- c("FD", "TNth", "TXth", "GDD", "CDD", "HDD", "P", "dt_st_rnagsn", "nm_flst_rnagsn", 
                       "dt_fnst_rnagsn", "dt_ed_rnagsn", "dl_agsn", "dc_agsn", "rn_agsn", 
                       "avrn_agsn", "dc_rnlg_agsn", "tm_agsn", "dc_txh_agsn", "dc_tnh_agsn")
  if (!index.code %in% choices) stop("Non valid index selected: Use indexShow() to select an index.")
  if (index.code == "FD") {
    index.arg.list[["th"]] <- 0
    message("[", Sys.time(), "] th = 0 for index FD. Use index.code = 'TNth' to set a different threshold")
  }
  time.resolution <- match.arg(time.resolution,
                               choices = c("month", "year", "climatology"))
  if (!is.null(tn)) {
    if (getTimeResolution(tn) != "DD") stop("Daily data is required as input", call. = FALSE)
  }
  if (!is.null(tx)) {
    if (getTimeResolution(tx) != "DD") stop("Daily data is required as input", call. = FALSE)
  }
  if (!is.null(tm)) {
    if (getTimeResolution(tm) != "DD") stop("Daily data is required as input", call. = FALSE)
  }
  if (!is.null(pr)) {
    if (getTimeResolution(pr) != "DD") stop("Daily data is required as input", call. = FALSE)
  }
  if (!is.null(baseline)) {
    if (!index.code %in% c("P")) {
      warning("Index.code is not 'P', baseline ignored")
      baseline <- NULL
    } else {
      if (getTimeResolution(baseline) != "DD") stop("Daily data is required as input", call. = FALSE)
    }
  }
  aux <- read.master()
  metadata <- aux[grep(paste0("^", index.code, "$"), aux$code, fixed = FALSE), ]
  a <- c(!is.null(tn), !is.null(tx), !is.null(tm), !is.null(pr)) %>% as.numeric()
  if (!index.code %in% c("P")) {
    b <- metadata[ , 4:7] %>% as.numeric()
    if (any(b - a > 0)) {
      stop("The required input variable(s) for ", index.code,
           " index calculation are missing\nType \'?",
           metadata$indexfun, "\' for help", call. = FALSE)
    }
  } else {
    b <- a
    if (sum(b) > 1) stop(index.code, " is applied to single variable.")
  }
  grid.list <- list("tn" = tn, "tx" = tx, "tm" = tm, "pr" = pr)[which(as.logical(b))]
  namesgridlist <- names(grid.list)
  # Operations for the consistency of the grids
  locs <- lapply(grid.list, isRegular)
  if (!any(sum(unlist(locs)) != 0,  sum(unlist(locs)) != length(grid.list))) stop("Regular and Irregular grids can not be combined. See function interpGrid")
  if (length(grid.list) > 1) {
    grid.list <- intersectGrid(grid.list, type = "temporal", which.return = 1:length(grid.list))
    names(grid.list) <- namesgridlist
    grid.list <- suppressMessages(lapply(grid.list, function(i) interpGrid(i, getGrid(grid.list[[1]]))))
  }
  grid.list <- lapply(grid.list, function(r) redim(r, drop = TRUE))
  grid.list <- lapply(grid.list, function(r) redim(r, loc = !unique(unlist(locs))))
  if (!unique(unlist(locs))) stop("The implementation of indexGrid for irregular grids is under development.")
  # Member loop preparation
  ns.mem <- lapply(grid.list, function(r) getShape(r)[["member"]])
  if (sum(unlist(ns.mem) - rep(ns.mem[[1]], length(ns.mem))) != 0) stop("Number of members is different")
  n.mem <- unique(unlist(ns.mem))
  if (n.mem > 1) {
    parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
    apply_fun <- selectPar.pplyFun(parallel.pars, .pplyFUN = "lapply")
    if (parallel.pars$hasparallel) on.exit(parallel::stopCluster(parallel.pars$cl))
  } else {
    if (isTRUE(parallel)) message("NOTE: Parallel processing was skipped (unable to parallelize one single member)")
    apply_fun <- lapply
  }
  # Member loop
  message("[", Sys.time(), "] Calculating ", index.code, " ...")
  out.m <- apply_fun(1:n.mem, function(m){
    if (sum(b) == 1 & is.null(baseline)) {
      # Indices from a single variable
      aggr.arg <- switch(time.resolution,
                         "month" = "aggr.m",
                         "year" = "aggr.y",
                         "climatology" = "clim.fun")
      fun.call <- switch(time.resolution,
                         "month" = "aggregateGrid",
                         "year" = "aggregateGrid",
                         "climatology" = "climatology")
      input.arg.list <- list()
      input.arg.list[["grid"]] <- subsetGrid(grid.list[[1]], members = m)
      input.arg.list[[aggr.arg]] <- c(list("FUN" = metadata$indexfun), index.arg.list)
      suppressMessages(do.call(fun.call, input.arg.list))
      
      
    } else {
      # Indices from multiple variables or for baseline methods
      grid.list.aux <- lapply(grid.list, function(x) subsetGrid(x, members = m))
      months <- switch(time.resolution,
                       "month" = as.list(getSeason(grid.list.aux[[1]])),
                       "year" = list(getSeason(grid.list.aux[[1]])),
                       "climatology" = list(getSeason(grid.list.aux[[1]])))
      years <- switch(time.resolution,
                      "month" = as.list(unique(getYearsAsINDEX(grid.list.aux[[1]]))),
                      "year" = as.list(unique(getYearsAsINDEX(grid.list.aux[[1]]))),
                      "climatology" = list(unique(getYearsAsINDEX(grid.list.aux[[1]]))))
      if (!is.null(baseline)) {
        baseline.sub <- suppressWarnings(subsetGrid(baseline, members = m))
        if (is.null(index.arg.list[["percent"]]) & is.null(index.arg.list[["value"]])) stop("Baseline provided but percent or value not specified.")
        if (!is.null(index.arg.list[["percent"]]) & !is.null(index.arg.list[["value"]])) {
          warning("Values were given to both percentile and value... value will be ignored (set to NULL)")
          index.arg.list[["value"]] <- NULL
        }
        if (!is.null(index.arg.list[["percent"]])) {
          index.arg.list[["value"]] <- suppressMessages(
            redim(climatology(baseline.sub, 
                              clim.fun = list(FUN = percentile, percent = index.arg.list[["percent"]])), drop = TRUE)[["Data"]])
          index.arg.list[["percent"]] <- NULL
        } else if (!is.null(index.arg.list[["value"]])) {
          index.arg.list[["percent"]] <- suppressMessages(
            redim(climatology(baseline.sub,
                              clim.fun = list(FUN = percentile, value = index.arg.list[["value"]])), drop = TRUE)[["Data"]])
          index.arg.list[["value"]] <- NULL
        }
      }
      # EXCEPTION for FAO INDICES (require lat, dates, and NO temporal subsetting)
      if (metadata$indexfun == "agroindexFAO") {
        if (time.resolution != "year") message(index.code, " is calculated yaear by year by definition. argument time.resolution ignored.")
        out.aux <- suppressMessages(aggregateGrid(grid.list.aux[[1]], aggr.y = list(FUN = "mean", na.rm = TRUE)))
        input.arg.list <- lapply(grid.list.aux, function(d) d[["Data"]])
        datess <- as.Date(grid.list.aux[[1]][["Dates"]][["start"]])
        datess <- cbind(as.numeric(format(datess, "%Y")), as.numeric(format(datess, "%m")), as.numeric(format(datess, "%d")))
        lats <- grid.list.aux[[1]][["xyCoords"]][["y"]]
        latloop <- lapply(1:length(lats), function(l) {
                              lonloop <- lapply(1:getShape(grid.list.aux[[1]])["lon"], function(lo) {
                                     do.call(metadata$indexfun, c(lapply(input.arg.list, function(z) z[, l, lo]), "lat" = list(lats[l]), "dates" = list(datess), "index.code" = list(index.code)))
                              })
                              do.call("abind", list(lonloop, along = 0))
                  })
        out.aux[["Data"]] <- unname(aperm(do.call("abind", list(latloop, along = 0)), c(3, 1, 2)))
        attr(out.aux[["Data"]], "dimensions") <- c("time", "lat", "lon")
        out.aux
      } else {
       yg <- lapply(years, function(yi){
          mg <- lapply(months, function(mi) {
            out.aux <- suppressMessages(climatology(grid.list.aux[[1]]))
            grid.list.sub <- lapply(grid.list.aux, function(x) subsetGrid(x, years = yi, season = mi))
            input.arg.list <- lapply(grid.list.sub, function(d) d[["Data"]])
            if (index.code == "P") names(input.arg.list) <- "var"
            input.arg.list <- c(input.arg.list, index.arg.list)
            out.aux[["Data"]] <- unname(do.call(metadata$indexfun, input.arg.list))
            attr(out.aux[["Data"]], "dimensions") <- c("lat", "lon")
            out.aux
          })
          tryCatch({bindGrid(mg, dimension = "time")}, error = function(err){unlist(mg, recursive = FALSE)})
        })
        tryCatch({bindGrid(yg, dimension = "time")}, error = function(err){unlist(yg, recursive = FALSE)})
      }
    }
  })
  out <- suppressMessages(suppressWarnings(bindGrid(out.m, dimension = "member")))
  out[["Variable"]] <- list("varName" = index.code, 
                            "level" = out[["Variable"]][["level"]])
  attr(out[["Variable"]], "description") <- metadata$description
  attr(out[["Variable"]], "units") <- metadata$units
  attr(out[["Variable"]], "longname") <- metadata$longname
  message("[", Sys.time(), "] Done")
  return(out)
}


#' @title List all available Indices
#' @description Print a table with a summary of the available indices
#' @return Print a table on the screen with the following columns:
#' \itemize{
#' \item \strong{code}: Code of the index. This is the character string used as input value
#' for the argument \code{index.code} in \code{\link{indexGrid}}
#' \item \strong{longname}: Long description of the index
#' \item \strong{index.fun}: The name of the internal function used to calculate it
#' \item \strong{tn, tx, tm, pr}: A logical value (0/1) indicating the input variables required for index calculation
#' \item \strong{units}: The units of the index (when different from those of the input variable)
#' }
#' @author J. Bedia, M. Iturbide
#' @export

indexShow <- function() {
  read.master()
}



#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom utils read.table

read.master <- function() {
  system.file("master", package = "climate4R.indices") %>% read.table(header = TRUE,
                                                                      sep = ";",
                                                                      stringsAsFactors = FALSE,
                                                                      na.strings = "")
}
