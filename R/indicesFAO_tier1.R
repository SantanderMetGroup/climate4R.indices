#' @title Wrapper to call FAO_tier1 index calculation functions.
#' @param index.code To call to the atomic function of the same name
#' @param ... Other parameters that can be passed to the selected index function.
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
binSpell <- function(v) {
  # Function to compute length of binary (e.g. 0=dry, 1=wet) spells
  # Returns a 0/1 sequence representing spells, as well as the duration of each spell
  
  # v: vector of data (e.g. daily precipitation)
  ix <- c(which(v[-length(v)] != v[-1]), length(v))  
  
  # output list
  out <- list()
  out$len <- diff(c(0, ix))
  out$val <- v[ix]
  return(out)
}

##################
## yearStartEnd ##
##################
yearStartEnd <- function(dates, year, year.start = NULL, year.end = NULL) {
  # Function to find the position marking the start and the end of a given year (or a user-defined portion of the year). The search is done within the "dates" matrix
  
  # dates: matrix containing the full range of dates (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
  # year: year of interest
  # year.start: user-defined start of the year [in "YYYY-MM-DD" format]
  # year.end: user-defined end of the year [in "YYYY-MM-DD" format]
  
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
#' @description see http://etccdi.pacificclimate.org/list_27_indices.shtml
#' @value Number of days/year
#' @param lat Latitude (Integer)
#' @param dates Matrix containing the full range of dates (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param tm Vector of mean temperature data
#' @param pnan Maximum percentage of missing data (tm) allowed in one year to compute the indices
#' @author R. Manzanas
#' @export
#' 
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
  # return(index)
  return(GSL)
}
#index = gsl(tm, dates, lat)  # call to the function


#########
## avg ##
#########
avg <- function(tm, dates, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # Function to compute the average value of a time-series
  
  # tm: vector of data (e.g. daily mean temperature)
  # dates: matrix containing the full range of dates (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
  # year: vector with years of interest
  # year.start: vector of dates [in "YYYY-MM-DD" format] defining the beginnig of a portion of interest within each year (e.g., the agronomic season)
  # year.end: vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
  # pnan: any year with a percentage of NA data above pnan will be ignored
  
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
nd_thre <- function(data, dates, threshold, direction = "geq", year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # Function to compute the number of days exceeding (either below or above) a given threshold 
  
  # data: vector of data (e.g. daily mean temperature)
  # dates: matrix containing the full range of dates (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
  # direction: "geq" (greater or equal to) or "leq" (lower or equal to)
  # year: vector with years of interest
  # year.start: vector of dates [in "YYYY-MM-DD" format] defining the beginnig of a portion of interest within each year (e.g., the agronomic season)
  # year.end: vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
  # pnan: any year with a percentage of NA data above pnan will be ignored
  
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
nhw <- function(tx, dates, threshold, duration, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # Function to compute the number of heatwaves (as defined by threshold and duration)
  
  # tx: vector of data (e.g. daily maximum temperature)
  # dates: matrix containing the full range of dates (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
  # threshold: threshold used to define the heatwave
  # duration: duration (in days) used to define te heatwave
  # year: vector with years of interest
  # year.start: vector of dates [in "YYYY-MM-DD" format] defining the beginnig of a portion of interest within each year (e.g., the agronomic season)
  # year.end: vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
  # pnan: any year with a percentage of NA data above pnan will be ignored
  
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
dr <- function(tx, tn, dates, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # Function to compute the diurnal temperature range
  
  # tx: vector with daily maximum temperature
  # tn: vector with daily minimum temperature
  # dates: matrix containing the full range of dates (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
  # year: vector with years of interest
  # year.start: vector of dates [in "YYYY-MM-DD" format] defining the beginnig of a portion of interest within each year (e.g., the agronomic season)
  # year.end: vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
  # pnan: any year with a percentage of NA data above pnan will be ignored
  
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
prcptot <- function(pr, dates, wet.threshold = 1, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # Function to compute total precipitation in wet days
  
  # pr: vector with daily precipitation
  # dates: matrix containing the full range of dates (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
  # wet.threshold: threshold used to define wet days
  # year: vector with years of interest
  # year.start: vector of dates [in "YYYY-MM-DD" format] defining the beginnig of a portion of interest within each year (e.g., the agronomic season)
  # year.end: vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
  # pnan: any year with a percentage of NA data above pnan will be ignored
  
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
nrd <- function(pr, dates, wet.threshold = 1, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # Function to compute number of rainy days
  
  # pr: vector with daily precipitation
  # dates: matrix containing the full range of dates (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
  # wet.threshold: threshold used to define wet days
  # year: vector with years of interest
  # year.start: vector of dates [in "YYYY-MM-DD" format] defining the beginnig of a portion of interest within each year (e.g., the agronomic season)
  # year.end: vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
  # pnan: any year with a percentage of NA data above pnan will be ignored
  
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
lds <- function(pr, dates, length.spell = "mean", wet.threshold = 1, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # Function to compute the average/maximum length of dry (as defined by wet.threshold) spells
  
  # pr: vector with daily precipitation
  # dates: matrix containing the full range of dates (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
  # length.spell: either "mean" or "maximum" length spell
  # wet.threshold: threshold used to define wet days
  # year: vector with years of interest
  # year.start: vector of dates [in "YYYY-MM-DD" format] defining the beginnig of a portion of interest within each year (e.g., the agronomic season)
  # year.end: vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
  # pnan: any year with a percentage of NA data above pnan will be ignored
  
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
sdii <- function(pr, dates, wet.threshold = 1, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # Function to compute precipitation intensity in wet (as defined by wet.threshold) days
  
  # pr: vector with daily precipitation
  # dates: matrix containing the full range of dates (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
  # wet.threshold: threshold used to define wet days
  # year: vector with years of interest
  # year.start: vector of dates [in "YYYY-MM-DD" format] defining the beginnig of a portion of interest within each year (e.g., the agronomic season)
  # year.end: vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
  # pnan: any year with a percentage of NA data above pnan will be ignored
  
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
prcptot_thre <- function(pr, dates, threshold = 50, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # Function to compute total amount of precipitation fallen in strong (as defined by threshold) rainy days
  
  # pr: vector with daily precipitation
  # dates: matrix containing the full range of dates (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
  # threshold: threshold used to define strong rain
  # year: vector with years of interest
  # year.start: vector of dates [in "YYYY-MM-DD" format] defining the beginnig of a portion of interest within each year (e.g., the agronomic season)
  # year.end: vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
  # pnan: any year with a percentage of NA data above pnan will be ignored
    
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
ns <- function(pr, dates, wet.threshold, duration, type.spell = "dry", year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # Function to compute the number of spells (either dry or wet), as defined by threshold and duration
  
  # pr: vector with daily precipitation
  # dates: matrix containing the full range of dates (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
  # wet.threshold: threshold used to define wet days
  # duration: duration (in days) of spells
  # type.spell: either "dry" or "wet" 
  # year: vector with years of interest
  # year.start: vector of dates [in "YYYY-MM-DD" format] defining the beginnig of a portion of interest within each year (e.g., the agronomic season)
  # year.end: vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
  # pnan: any year with a percentage of NA data above pnan will be ignored
   
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
            index[year == iyear] = sum(bin$len[!which(bin$val)] >= duration, na.rm = T) 
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
