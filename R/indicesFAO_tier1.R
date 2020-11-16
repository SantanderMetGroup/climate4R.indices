#' @title Wrapper to call FAO_tier1 index calculation functions.
#' @param index.code To call to the atomic function of the same name
#' @param ... Other parameters that can be passed to the selected index function (see Details).
#' @details Index calculation is done by different functions that receive the same name 
#' as the corresponding index code. Therefore, type \code{?index.code} to know the arguments
#' that are used in each case, e.g. \code{?gsl}.
#' @section Index codes:
#'    \strong{gsl} 
#'    
#'    \strong{avg}  
#'    
#'    \strong{nd_thre}
#'  
#'    \strong{nhw}
#'    
#'    \strong{dr}
#'    
#'    \strong{prcptot}
#'    
#'    \strong{nrd}
#'    
#'    \strong{lds}
#'    
#'    \strong{sdii}
#'    
#'    \strong{prcptot_thre}
#'    
#'    \strong{ns}
#'    
#' @author M. Iturbide
#' @export

agroindexFAO_tier1 <- function(index.code, ...) {
  choices <- c("gsl", "avg", "nd_thre", "nhw", "dr", "prcptot", "nrd", "lds", "sdii", "prcptot_thre", "ns")
  if (!index.code %in% choices) stop("Non valid index selected: Use indexShow() to select an index.")
  do.call(index.code, list(...))
}

#########################
## auxiliary functions ##
#########################

##############
## binSpell ##
##############
#' @title Function to compute number and length of binary (e.g. 0=dry, 1=wet) spells
#' @return Sequence of spells, including the duration of each spell. The search is done within the "data" vector
#' @param data Vector with data (e.g. daily precipitation)
#' @author R. Manzanas
#' @export

binSpell <- function(data) {
  
  ix <- c(which(data[-length(data)] != data[-1]), length(data))  
  
  # output list
  out <- list()
  out$len <- diff(c(0, ix))
  out$val <- data[ix]
  return(out)
}

##################
## yearStartEnd ##
##################
#' @title Function to find the position marking the start and the end of a given year (or a user-defined portion of the year)
#' @return Indices marking the start and the end of a given year (or a user-defined portion of the year). The search is done within the "dates" matrix
#' @param dates Matrix containing the full range of dates (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param year Year of interest (e.g. 1995)
#' @param year.start User-defined start of the year [in "YYYY-MM-DD" format]
#' @param year.end User-defined end of the year [in "YYYY-MM-DD" format]
#' @author R. Manzanas
#' @export

yearStartEnd <- function(dates, year, year.start = NULL, year.end = NULL) {
  
  if (!is.null(year.start) & !is.null(year.end)) {
    ind.start = which(dates[, 1] == as.numeric(substr(year.start, 1, 4)) &   # start of the year (as defined by the user)
                        dates[, 2] == as.numeric(substr(year.start, 6, 7)) & 
                        dates[, 3] == as.numeric(substr(year.start, 9, 10)))
    ind.end = which(dates[, 1] == as.numeric(substr(year.end, 1, 4)) &   # end of the year (as defined by the user)
                      dates[, 2] == as.numeric(substr(year.end, 6, 7)) & 
                      dates[, 3] == as.numeric(substr(year.end, 9, 10)))
    
  } else {
    ind.start = which(dates[, 1] == year & dates[, 2] == 1 & dates[, 3] == 1)  # year's definition by default (1-Jan to 31-Dec)
    ind.end = which(dates[, 1] == year & dates[, 2] == 12 & dates[, 3] == 31)  # year's definition by default (1-Jan to 31-Dec)
  }
  out = list()
  out$start = ind.start
  out$end = ind.end
  return(out)
}

#################################################################################################################################
## tier1 indices (https://docs.google.com/document/d/1iCodxrYPhlvSwEK24zMGmjwshD0UkE2NFcBcx8lYo0k/edit#heading=h.omp2879zs40p) ##
#################################################################################################################################

#########
## gsl ##
#########
#' @title Function to compute the growing season length
#' @description See http://etccdi.pacificclimate.org/list_27_indices.shtml for details
#' @return Number of days (per year)
#' @param tm Vector with mean temperature data
#' @param lat Latitude (in degrees)
#' @param dates Matrix containing the full range of dates corresponding to "tm" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param pnan Maximum percentage of missing data allowed in one year to compute the indices
#' @author R. Manzanas
#' @export

gsl <- function(tm, dates, lat, pnan = 25) {
  
  year = unique(dates[, 1])  # years of analysis
  
  ## initializing output vectors
  init = rep(NA, length(year))  # first day of the season
  end = rep(NA, length(year))  # last day of the season
  GSL = rep(NA, length(year))  # growing season length 
  
  for (iyear in year) {
    
    ## definition of the year
    if (lat >= 0) {  # northern hemisphere (year: Jan-Dec)
      year.init = which(dates[, 1] == iyear & dates[, 2] == 1 & dates[, 3] == 1)
      year.end = which(dates[, 1] == iyear & dates[, 2] == 12 & dates[, 3] == 31)
    } else {  # southern hemisphere (year: Jul-Jun)
      year.init = which(dates[, 1] == iyear & dates[, 2] == 7 & dates[, 3] == 1)
      year.end = which(dates[, 1] == iyear + 1 & dates[, 2] == 6 & dates[, 3] == 30)
    }
    
    if (length(year.init) > 0 & length(year.end) > 0) {  # checking for complete year
      ind.year = year.init:year.end;  nday = length(ind.year)
      dates.year = dates[ind.year, ]
      
      ## tm in year
      tm.year = tm[ind.year]
      rm(ind.year)
      
      if (lat >= 0) {  # northern hemisphere (index determining 1-Jul)
        ind.mid.year = which(dates.year[, 1] == iyear & dates.year[, 2] == 7 & dates.year[, 3] == 1)
      } else {  # southern hemisphere (index determining 1-Jan)
        ind.mid.year  = which(dates.year[, 1] == iyear & dates.year[, 2] == 1 & dates.year[, 3] == 1)
      }
      
      if (100*sum(is.na(tm.year))/nday < pnan) {  # asking for a minimum of non-missing data for computing the index 
        
        ## first 6-day spell with tm>5degC (within the first half of the year)
        binwarm = binSpell(tm.year > 5)
        indblockwarm = which(binwarm$val & (binwarm$len >= 6))[1]
        if ((length(indblockwarm) != 0) & (!is.na(indblockwarm))) {
          aux.ind.dinit = sum(binwarm$len[1:(indblockwarm-1)]) + 6  # end of first 6-day spell with tm>5degC
          if (aux.ind.dinit < ind.mid.year) {
            dinit = aux.ind.dinit
          } else {
            dinit = NA
          }
        } else {
          dinit = NA
        }
        
        ## first 6-day spell with tm<5degC (within the second half of the year)
        bincold = binSpell(tm.year < 5)
        indblockcold = which(bincold$val & (bincold$len >= 6))
        if (length(indblockcold) == 0) {
          indblockcold = NA
        }
        if ((length(indblockcold) != 0) & (!is.na(indblockcold))) {
          aux.ind.dend = c()
          for (j in indblockcold){
            aux.ind.dend[indblockcold == j] = sum(binwarm$len[1:(j-1)]) + 6  # end of 6-day spells with tm<5degC
          }
          ind.dend = which(aux.ind.dend > ind.mid.year)[1]
          if (length(ind.dend) != 0 & (!is.na(ind.dend))) {
            dend = aux.ind.dend[ind.dend]  
          } else {
            dend = NA
          }
        } else {
          dend = NA
        }
        
        init[which(year == iyear)] = sprintf("%04d-%02d-%02d", dates.year[dinit, 1], dates.year[dinit, 2], dates.year[dinit, 3])  # first day of the season
        end[which(year == iyear)] = sprintf("%04d-%02d-%02d", dates.year[dend, 1], dates.year[dend, 2], dates.year[dend, 3])  # last day of the season
        GSL[which(year == iyear)] = dend - dinit  # growing season length 
      }
    } else {
      message(sprintf("... year %d not complete ...", iyear))
    }
  }
  
  # output list
  index = list()
  index$year = year
  index$init = init
  index$end = end
  index$GSL = GSL
  return(index)
}
#index = gsl(tm, dates, lat)  # call to the function


#########
## avg ##
#########
#' @title Function to compute the average value of a daily time-series 
#' @return Average value (per year)
#' @param tm Vector with data (e.g. daily mean temperature)
#' @param dates Matrix containing the full range of dates corresponding to "tm" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param year Vector with years of interest (e.g. 1990:1995)
#' @param year.start Vector of dates [in "YYYY-MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season)
#' @param year.end Vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored
#' @param lat Latitude (NULL) to indicate that latitude information is not used.
#' @author R. Manzanas
#' @export

avg <- function(tm, dates, year = NULL, year.start = NULL, year.end = NULL, pnan = 25, lat = NULL) {
  if(!is.null(lat)) warning("This index doesn't use latitude information.")
  if (is.null(year)) {
    year = unique(dates[, 1])  # years of analysis
  }
  
  # initializing output
  index = rep(NA, 1, length(year))  
  
  for (iyear in year) {
    
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year = yearStartEnd(dates, iyear, year.start = year.start[year == iyear], year.end = year.end[year == iyear])  # bounding dates defining the portion of year of interest
    } else {
      ind.year = yearStartEnd(dates, iyear, year.start = NULL, year.end = NULL)  # bounding dates defining the year of interest
    }
    
    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {
        tm.year = tm[ind.year$start:ind.year$end]
        if (sum(is.na(tm.year)) < 0.01*pnan*length(tm.year)) {  # asking for a minimum of pnan (%) of non-missing days 
          index[year == iyear] = mean(tm.year, na.rm = T)  
        }
      }
    }
  }
  return(index)
}

# index = avg(tm, dates)  # call to the function
# index = avg(tm, dates, year = 1994:2018,  # call to the function
#                year.start = AS.dates[[1]]$sdate$days[, 15], 
#                year.end = AS.dates[[1]]$edate$days[, 15])  

#############
## nd_thre ##
#############
#' @title Function to compute the number of days exceeding (either below or above) a given threshold 
#' @return Number of days (per year)
#' @param any Vector with data of ANY variable (e.g. daily mean temperature)
#' @param dates Matrix containing the full range of dates corresponding to "data" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param threshold Threshold considered. Must be in the same units of "data"
#' @param direction "geq" (greater or equal to) or "leq" (lower or equal to)
#' @param year Vector with years of interest (e.g. 1990:1995)
#' @param year.start Vector of dates [in "YYYY-MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season)
#' @param year.end Vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored
#' @param lat Latitude (NULL) to indicate that latitude information is not used.
#' @author R. Manzanas
#' @export

nd_thre <- function(any, dates, threshold, direction = "geq", year = NULL, year.start = NULL, year.end = NULL, pnan = 25, lat = NULL) {
  data <- any
  if(!is.null(lat)) warning("This index doesn't use latitude information.")
  
  if (is.null(year)) {
    year = unique(dates[, 1])  # years of analysis
  }
  
  # initializing output
  index = rep(NA, 1, length(year))  
  
  for (iyear in year) {
    
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year = yearStartEnd(dates, iyear, year.start = year.start[year == iyear], year.end = year.end[year == iyear])  # bounding dates defining the portion of year of interest
    } else {
      ind.year = yearStartEnd(dates, iyear, year.start = NULL, year.end = NULL)  # bounding dates defining the year of interest
    }
    
    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$end) & !is.na(ind.year$end)) {
        data.year = data[ind.year$start:ind.year$end]
        if (sum(is.na(data.year)) < 0.01*pnan*length(data.year)) {  # asking for a minimum of pnan (%) of non-missing days 
          if (direction == "geq") {
            index[year == iyear] = sum(data.year >= threshold, na.rm = T)
          } else if (direction == "leq") {
            index[year == iyear] = sum(data.year <= threshold, na.rm = T)
          }
        }
      }
    }
  }
  return(index)
}
# index = nd_thre(tm, dates, 26, direction = "geq")  # call to the function
# index = nd_thre(pr, dates, 50, direction = "geq")  
# index = nd_thre(tm, dates, 26, direction = "geq",  # call to the function
#                year = 1994:2018, 
#                year.start = AS.dates[[1]]$sdate$days[, 15], 
#                year.end = AS.dates[[1]]$edate$days[, 15])  
# index = nd_thre(pr, dates, 50, direction = "geq",  # call to the function
#                year = 1994:2018, 
#                year.start = AS.dates[[1]]$sdate$days[, 15], 
#                year.end = AS.dates[[1]]$edate$days[, 15])  

#########
## nhw ##
#########
#' @title Function to compute the number of heatwaves (as defined by threshold and duration)
#' @return Number of heatwaves (per year)
#' @param tx Vector with daily maximum temperature
#' @param dates Matrix containing the full range of dates corresponding to "tx" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param threshold Threshold considered. Must be in the same units of "tx"
#' @param duration Duration (in days) considered to define the heatwave
#' @param year Vector with years of interest (e.g. 1990:1995)
#' @param year.start Vector of dates [in "YYYY-MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season)
#' @param year.end Vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored
#' @param lat Latitude (NULL) to indicate that latitude information is not used.
#' @author R. Manzanas
#' @export

nhw <- function(tx, dates, threshold, duration, year = NULL, year.start = NULL, year.end = NULL, pnan = 25, lat = NULL) {
  if(!is.null(lat)) warning("This index doesn't use latitude information.")
  
  if (is.null(year)) {
    year = unique(dates[, 1])  # years of analysis
  }
  
  # initializing output
  index = rep(NA, 1, length(year))  
  
  for (iyear in year) {
    
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year = yearStartEnd(dates, iyear, year.start = year.start[year == iyear], year.end = year.end[year == iyear])  # bounding dates defining the portion of year of interest
    } else {
      ind.year = yearStartEnd(dates, iyear, year.start = NULL, year.end = NULL)  # bounding dates defining the year of interest
    }
    
    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {
        tx.year = tx[ind.year$start:ind.year$end]
        if (sum(is.na(tx.year)) < 0.01*pnan*length(tx.year)) {  # asking for a minimum of pnan (%) of non-missing days 
          
          bin = binSpell(tx.year > threshold)
          index[year == iyear] = sum(bin$len[which(bin$val)] >= duration, na.rm = T)
        }
      }
    }
  }
  return(index)
}
# index = nhw(tm, dates, 26, 4)  # call to the function
# index = nhw(tm, dates, 26, 4,  # call to the function
#             year = 1994:2018, 
#             year.start = AS.dates[[1]]$sdate$days[, 15], 
#             year.end = AS.dates[[1]]$edate$days[, 15])  

########
## dr ##
########
#' @title Function to compute the diurnal temperature range
#' @return Mean diurnal temperature range (per year)
#' @param tx Vector with daily maximum temperature
#' @param tn Vector with daily minimum temperature
#' @param dates Matrix containing the full range of dates corresponding to "tx" and "tn" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param year Vector with years of interest (e.g. 1990:1995)
#' @param year.start Vector of dates [in "YYYY-MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season)
#' @param year.end Vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored
#' @param lat Latitude (NULL) to indicate that latitude information is not used.
#' @author R. Manzanas
#' @export

dr <- function(tx, tn, dates, year = NULL, year.start = NULL, year.end = NULL, pnan = 25, lat = NULL) {
  if(!is.null(lat)) warning("This index doesn't use latitude information.")
  
  if (is.null(year)) {
    year = unique(dates[, 1])  # years of analysis
  }
  
  # initializing output
  index = rep(NA, 1, length(year))  
  
  for (iyear in year) {
    
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year = yearStartEnd(dates, iyear, year.start = year.start[year == iyear], year.end = year.end[year == iyear])  # bounding dates defining the portion of year of interest
    } else {
      ind.year = yearStartEnd(dates, iyear, year.start = NULL, year.end = NULL)  # bounding dates defining the year of interest
    }
    
    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {
        tx.year = tx[ind.year$start:ind.year$end]
        tn.year = tn[ind.year$start:ind.year$end]
        if (sum(is.na(tx.year)) < 0.01*pnan*length(tx.year) & sum(is.na(tn.year)) < 0.01*pnan*length(tn.year)) {  # asking for a minimum of pnan (%) of non-missing days 
          
          index[year == iyear] = mean(tx.year - tn.year, na.rm = T)
        }
      }
    }
  }
  return(index)
}
# index = dr(tx, tn, dates)  # call to the function
# index = dr(tx, tn, dates,   # call to the function
#             year = 1994:2018, 
#             year.start = AS.dates[[1]]$sdate$days[, 15], 
#             year.end = AS.dates[[1]]$edate$days[, 15])  

#############
## prcptot ##
#############
#' @title Function to compute total precipitation in wet days
#' @return Total precipitation in wet days (per year)
#' @param pr Vector with daily precipitation
#' @param dates Matrix containing the full range of dates corresponding to "pr" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param wet.threshold Threshold considered to define wet days. Must be in the same units of "pr"
#' @param year Vector with years of interest (e.g. 1990:1995)
#' @param year.start Vector of dates [in "YYYY-MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season)
#' @param year.end Vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored
#' @param lat Latitude (NULL) to indicate that latitude information is not used.
#' @author R. Manzanas
#' @export

prcptot <- function(pr, dates, wet.threshold = 1, year = NULL, year.start = NULL, year.end = NULL, pnan = 25, lat = NULL) {
  if(!is.null(lat)) warning("This index doesn't use latitude information.")
  
  if (is.null(year)) {
    year = unique(dates[, 1])  # years of analysis
  }
  
  # initializing output
  index = rep(NA, 1, length(year))  
  
  for (iyear in year) {
    
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year = yearStartEnd(dates, iyear, year.start = year.start[year == iyear], year.end = year.end[year == iyear])  # bounding dates defining the portion of year of interest
    } else {
      ind.year = yearStartEnd(dates, iyear, year.start = NULL, year.end = NULL)  # bounding dates defining the year of interest
    }
    
    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {
        pr.year = pr[ind.year$start:ind.year$end]
        if (sum(is.na(pr.year)) < 0.01*pnan*length(pr.year)) {  # asking for a minimum of pnan (%) of non-missing days 
          
          index[year == iyear] = sum(pr.year[pr.year >= wet.threshold], na.rm = T)
        }
      }
    }
  }
  return(index)
}
# index = prcptot(pr, dates, wet.threshold = 2.5)  # call to the function
# index = prcptot(pr, dates, wet.threshold = 2.5,  # call to the function
#            year = 1994:2018, 
#            year.start = AS.dates[[1]]$sdate$days[, 15], 
#            year.end = AS.dates[[1]]$edate$days[, 15])  

#########
## nrd ##
#########
#' @title Function to compute number of rainy days
#' @return Number of rainy days (per year)
#' @param pr Vector with daily precipitation
#' @param dates Matrix containing the full range of dates corresponding to "pr" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param wet.threshold Threshold considered to define wet days. Must be in the same units of "pr"
#' @param year Vector with years of interest (e.g. 1990:1995)
#' @param year.start Vector of dates [in "YYYY-MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season)
#' @param year.end Vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored
#' @param lat Latitude (NULL) to indicate that latitude information is not used.
#' @author R. Manzanas
#' @export

nrd <- function(pr, dates, wet.threshold = 1, year = NULL, year.start = NULL, year.end = NULL, pnan = 25, lat = NULL) {
  if(!is.null(lat)) warning("This index doesn't use latitude information.")
  
  if (is.null(year)) {
    year = unique(dates[, 1])  # years of analysis
  }
  
  # initializing output
  index = rep(NA, 1, length(year))  
  
  for (iyear in year) {
    
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year = yearStartEnd(dates, iyear, year.start = year.start[year == iyear], year.end = year.end[year == iyear])  # bounding dates defining the portion of year of interest
    } else {
      ind.year = yearStartEnd(dates, iyear, year.start = NULL, year.end = NULL)  # bounding dates defining the year of interest
    }
    
    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {
        pr.year = pr[ind.year$start:ind.year$end]
        if (sum(is.na(pr.year)) < 0.01*pnan*length(pr.year)) {  # asking for a minimum of pnan (%) of non-missing days 
          
          index[year == iyear] = sum(pr.year >= wet.threshold, na.rm = T)
        }
      }
    }
  }
  return(index)
}
# index = nrd(pr, dates, wet.threshold = 2.5)  # call to the function
# index = nrd(pr, dates, wet.threshold = 2.5,  # call to the function
#                 year = 1994:2018, 
#                 year.start = AS.dates[[1]]$sdate$days[, 15], 
#                 year.end = AS.dates[[1]]$edate$days[, 15])  

#########
## lds ##
#########
#' @title Function to compute the average/maximum length of dry (as defined by wet.threshold) spells
#' @return Average/maximum length of dry spells (per year)
#' @param pr Vector with daily precipitation
#' @param dates Matrix containing the full range of dates corresponding to "pr" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param length.spell Either "mean" or "maximum" length spell
#' @param wet.threshold Threshold considered to define wet days. Must be in the same units of "pr"
#' @param year Vector with years of interest (e.g. 1990:1995)
#' @param year.start Vector of dates [in "YYYY-MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season)
#' @param year.end Vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored
#' @param lat Latitude (NULL) to indicate that latitude information is not used.
#' @author R. Manzanas
#' @export

lds <- function(pr, dates, length.spell = "mean", wet.threshold = 1, year = NULL, year.start = NULL, year.end = NULL, pnan = 25, lat = NULL) {
  if(!is.null(lat)) warning("This index doesn't use latitude information.")
  
  if (is.null(year)) {
    year = unique(dates[, 1])  # years of analysis
  }
  
  # initializing output
  index = rep(NA, 1, length(year))  
  
  for (iyear in year) {
    
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year = yearStartEnd(dates, iyear, year.start = year.start[year == iyear], year.end = year.end[year == iyear])  # bounding dates defining the portion of year of interest
    } else {
      ind.year = yearStartEnd(dates, iyear, year.start = NULL, year.end = NULL)  # bounding dates defining the year of interest
    }
    
    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {
        pr.year = pr[ind.year$start:ind.year$end]
        if (sum(is.na(pr.year)) < 0.01*pnan*length(pr.year)) {  # asking for a minimum of pnan (%) of non-missing days 
          
          bin = binSpell(pr.year < wet.threshold)
          if (length.spell == "mean") {
            index[year == iyear] = mean(bin$len[which(bin$val)], na.rm = T)
          } else if (length.spell == "maximum") {
            index[year == iyear] = max(bin$len[which(bin$val)], na.rm = T)
          }
        }
      }
    }
  }
  return(index)
}
# index = lds(pr, dates, length.spell = "mean", wet.threshold = 2.5)  # call to the function
# index = lds(pr, dates, length.spell = "maximum", wet.threshold = 2.5)  
# index = lds(pr, dates, length.spell = "mean", wet.threshold = 2.5,  # call to the function
#             year = 1994:2018, 
#             year.start = AS.dates[[1]]$sdate$days[, 15], 
#             year.end = AS.dates[[1]]$edate$days[, 15])  
# index = lds(pr, dates, length.spell = "maximum", wet.threshold = 2.5,
#             year = 1994:2018, 
#             year.start = AS.dates[[1]]$sdate$days[, 15], 
#             year.end = AS.dates[[1]]$edate$days[, 15]) 

##########
## sdii ##
##########
#' @title Function to compute precipitation intensity in wet (as defined by wet.threshold) days
#' @return Precipitation intensity in wet days (per year)
#' @param pr Vector with daily precipitation
#' @param dates Matrix containing the full range of dates corresponding to "pr" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param wet.threshold Threshold considered to define wet days. Must be in the same units of "pr"
#' @param year Vector with years of interest (e.g. 1990:1995)
#' @param year.start Vector of dates [in "YYYY-MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season)
#' @param year.end Vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored
#' @param lat Latitude (NULL) to indicate that latitude information is not used.
#' @author R. Manzanas
#' @export

sdii <- function(pr, dates, wet.threshold = 1, year = NULL, year.start = NULL, year.end = NULL, pnan = 25, lat = NULL) {
  if(!is.null(lat)) warning("This index doesn't use latitude information.")
  
  
  if (is.null(year)) {
    year = unique(dates[, 1])  # years of analysis
  }
  
  # initializing output
  index = rep(NA, 1, length(year))  
  
  for (iyear in year) {
    
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year = yearStartEnd(dates, iyear, year.start = year.start[year == iyear], year.end = year.end[year == iyear])  # bounding dates defining the portion of year of interest
    } else {
      ind.year = yearStartEnd(dates, iyear, year.start = NULL, year.end = NULL)  # bounding dates defining the year of interest
    }
    
    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {
        pr.year = pr[ind.year$start:ind.year$end]
        if (sum(is.na(pr.year)) < 0.01*pnan*length(pr.year)) {  # asking for a minimum of pnan (%) of non-missing days
          index[year == iyear] = mean(pr.year[pr.year >= wet.threshold], na.rm = T)
        }
      }
    }
  }
  return(index)
}
# index = sdii(pr, dates, wet.threshold = 2.5)  # call to the function
# index = sdii(pr, dates, wet.threshold = 2.5,  # call to the function
#             year = 1994:2018, 
#             year.start = AS.dates[[1]]$sdate$days[, 15], 
#             year.end = AS.dates[[1]]$edate$days[, 15])  

##################
## prcptot_thre ##
##################
#' @title Function to compute total amount of precipitation fallen in strong (as defined by threshold) rainy days
#' @return Total amount of precipitation fallen in strong rainy days (per year)
#' @param pr Vector with daily precipitation
#' @param dates Matrix containing the full range of dates corresponding to "pr" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param threshold Threshold considered to define strong rain. Must be in the same units of "pr"
#' @param year Vector with years of interest (e.g. 1990:1995)
#' @param year.start Vector of dates [in "YYYY-MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season)
#' @param year.end Vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored
#' @param lat Latitude (NULL) to indicate that latitude information is not used.
#' @author R. Manzanas
#' @export

prcptot_thre <- function(pr, dates, threshold = 50, year = NULL, year.start = NULL, year.end = NULL, pnan = 25, lat = NULL) {
  if(!is.null(lat)) warning("This index doesn't use latitude information.")
  
  
  if (is.null(year)) {
    year = unique(dates[, 1])  # years of analysis
  }
  
  # initializing output
  index = rep(NA, 1, length(year))  
  
  for (iyear in year) {
    
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year = yearStartEnd(dates, iyear, year.start = year.start[year == iyear], year.end = year.end[year == iyear])  # bounding dates defining the portion of year of interest
    } else {
      ind.year = yearStartEnd(dates, iyear, year.start = NULL, year.end = NULL)  # bounding dates defining the year of interest
    }
    
    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {
        pr.year = pr[ind.year$start:ind.year$end]
        if (sum(is.na(pr.year)) < 0.01*pnan*length(pr.year)) {  # asking for a minimum of pnan (%) of non-missing days
          index[year == iyear] = sum(pr.year[pr.year >= threshold], na.rm = T)
        }
      }
    }
  }
  return(index)
}
# index = prcptot_thre(pr, dates, threshold = 50)  # call to the function
# index = prcptot_thre(pr, dates, threshold = 2.5,  # call to the function
#              year = 1994:2018, 
#              year.start = AS.dates[[1]]$sdate$days[, 15], 
#              year.end = AS.dates[[1]]$edate$days[, 15])  

########
## ns ##
########
#' @title Function to compute the number of spells (either dry or wet), as defined by threshold and duration
#' @return Number of spells (either dry or wet) (per year)
#' @param pr Vector with daily precipitation
#' @param dates Matrix containing the full range of dates corresponding to "pr" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param wet.threshold Threshold considered to define wet days. Must be in the same units of "pr"
#' @param duration Duration (in days) of spells
#' @param type.spell Either "dry" or "wet"
#' @param year Vector with years of interest (e.g. 1990:1995)
#' @param year.start Vector of dates [in "YYYY-MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season)
#' @param year.end Vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored
#' @param lat Latitude (NULL) to indicate that latitude information is not used.
#' @author R. Manzanas
#' @export

ns <- function(pr, dates, wet.threshold, duration, type.spell = "dry", year = NULL, year.start = NULL, year.end = NULL, pnan = 25, lat = NULL) {
  if(!is.null(lat)) warning("This index doesn't use latitude information.")
  
  
  if (is.null(year)) {
    year = unique(dates[, 1])  # years of analysis
  }
  
  # initializing output
  index = rep(NA, 1, length(year))  
  
  for (iyear in year) {
    
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year = yearStartEnd(dates, iyear, year.start = year.start[year == iyear], year.end = year.end[year == iyear])  # bounding dates defining the portion of year of interest
    } else {
      ind.year = yearStartEnd(dates, iyear, year.start = NULL, year.end = NULL)  # bounding dates defining the year of interest
    }
    
    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {
        pr.year = pr[ind.year$start:ind.year$end]
        if (sum(is.na(pr.year)) < 0.01*pnan*length(pr.year)) {  # asking for a minimum of pnan (%) of non-missing days 
          
          bin = binSpell(pr.year >= wet.threshold)
          if (type.spell == "wet") {
            index[year == iyear] = sum(bin$len[which(bin$val)] >= duration, na.rm = T)
          } else if (type.spell == "dry") {
            index[year == iyear] = sum(bin$len[which(!bin$val)] >= duration, na.rm = T) 
          }
        }
      }
    }
  }
  return(index)
}
# index = ns(pr, dates, wet.threshold = 2.5, duration = 5, type.spell = "dry")  # call to the function
# index = ns(pr, dates, wet.threshold = 2.5, duration = 5, type.spell = "wet")  
# index = ns(tm, dates, wet.threshold = 2.5, duration = 5, type.spell = "dry",   # call to the function
#             year = 1994:2018, 
#             year.start = AS.dates[[1]]$sdate$days[, 15], 
#             year.end = AS.dates[[1]]$edate$days[, 15]) 
# index = ns(tm, dates, wet.threshold = 2.5, duration = 5, type.spell = "wet",   # call to the function
#            year = 1994:2018, 
#            year.start = AS.dates[[1]]$sdate$days[, 15], 
#            year.end = AS.dates[[1]]$edate$days[, 15]) 
