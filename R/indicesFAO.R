

#' @title FAO Agronomic indices
#' @description Compute the selected FAO Agronomic index
#' @param index.code Character of the index code to be computed. Options are:
#'   "dt_st_rnagsn", "nm_flst_rnagsn", "dt_fnst_rnagsn", "dt_ed_rnagsn", "dl_agsn = rep", "dc_agsn, "rn_agsn, 
#'   "avrn_agsn", "dc_rnlg_agsn", "tm_agsn", "dc_txh_agsn", "dc_tnh_agsn" (See Details).
#' @details 
#' dt_st_rnagsn = first start date of the Agronomic Season (AS)
#' nm_flst_rnagsn =  number of false starts of the AS
#' dt_fnst_rnagsn =  final start date of the AS
#' dt_ed_rnagsn =  end date of the AS
#' dl_agsn =  length of the AS
#' dc_agsn =  number of rainy days in AS
#' rn_agsn =  total amount of rain in AS
#' avrn_agsn =  average amount of rain in rainy days in AS
#' dc_rnlg_agsn =  number of days with "high" rain in AS
#' tm_agsn =  average daily tas in AS
#' dc_txh_agsn = number of "hot" tmax days in AS
#' dc_tnh_agsn = number of "hot" tmin days in AS
#' @param lat Latitude (Integer)
#' @param dates Vector of dates 
#' @param pr Vector of precipitation data
#' @param tx Vector of maximum temperature data
#' @param tn Vector of minimum temperature data
#' @param pnan maximum percentage of missing data (pr, tx, tn) allowed in one year to compute the indices
#' @param shc maximum soil holding capacity
#' @param rndy amount of rain in a day to be considered as “rainy” day (in mm/day)
#' @param rnlg “large” rainfall event (mm/day) 
#' @param txh “hot” tmax (degC)
#' @param tnh “hot” tmin (degC)
#' @import transformeR
#' @importFrom stats na.omit
#' @author R. Manzanas
#' @export

agroindexFAO = function(lat, dates, index.code, pr = NULL, tx = NULL, tn = NULL, 
                        pnan = 25, shc = 100, rndy = 2.5, 
                        rnlg = 50, txh = 30, tnh = 18) {
  ## FUNCTION BOD
  tm = (tx + tn) / 2
  year = unique(dates[, 1])  # years of analysis
  
  ## initializing output vectors
  dt_st_rnagsn = rep(NA, length(year)) # first start date of the Agronomic Season (AS)
  nm_flst_rnagsn = rep(NA, length(year)) # number of false starts of the AS
  dt_fnst_rnagsn = rep(NA, length(year)) # final start date of the AS
  dt_ed_rnagsn = rep(NA, length(year)) # end date of the AS
  dl_agsn = rep(NA, length(year)) # length of the AS
  
  
  out = rep(NA, length(year)) # number of rainy days in AS
  
  for (iyear in year) {
    ## definition of the agronomic year
    if (lat >= 0) {  # northern hemisphere (agronomic year: Jul-Jun)
      year.init = which(dates[, 1] == iyear & dates[, 2] == 1 & dates[, 3] == 1)
      year.end = which(dates[, 1] == iyear & dates[, 2] == 12 & dates[, 3] == 31)
    } else {  # southern hemisphere (agronomic year: Jan-Dec)
      year.init = which(dates[, 1] == iyear & dates[, 2] == 7 & dates[, 3] == 1)
      year.end = which(dates[, 1] == iyear + 1 & dates[, 2] == 6 & dates[, 3] == 30)
    }
    
    if (length(year.init) > 0 & length(year.end) > 0) {  # checking for complete agronomic year
      ind.year = year.init:year.end;  nday = length(ind.year)
      dates.year = dates[ind.year, ]
      
      ## pr, tmax, tmin and tm in agronomic year
      pr.year = pr[ind.year]
      tx.year = tx[ind.year]
      tn.year = tn[ind.year]
      tm.year = tm[ind.year]
      rm(ind.year)
      
      if (max(c(100*sum(is.na(pr.year))/nday,
                100*sum(is.na(tx.year))/nday,
                100*sum(is.na(tn.year))/nday)) < pnan) { # asking for a minimum of non-missing data for computing indices 
        
        #########################################
        ## ET0 (according to Hargreaves equation)
        #########################################
        ET0 = computeET0(lat, dates.year, tx.year, tn.year) 
        
        ET0[is.na(ET0)] = mean(ET0, na.rm = T)  # filling NA values with the average ET0 (to not affect greatly the computation of ET0)
        ET0[ET0 == 0] = mean(ET0, na.rm = T)  # filling 0 values with the average ET0 (to not affect greatly the computation of ET0)
        
        #######################################
        ## Water Soil Content (WSC), day by day
        #######################################
        pr.year.WSC = pr.year
        pr.year.WSC[pr.year.WSC > rnlg] = rnlg  # upper limit of pr that can contribute to WSC
        pr.year.WSC[is.na(pr.year.WSC)] = 0  # filling NA values with 0 (to not affect greatly the computation of WSC)
        ## first day of agronomic year
        WSC = rep(NA, nday)
        aux.WSC = pr.year.WSC[1] - ET0[1]
        if (!is.na(aux.WSC) & aux.WSC < 0) {
          aux.WSC = 0
        } else if (!is.na(aux.WSC) & aux.WSC > shc) {
          aux.WSC = shc  # anything above shc goes to runoff
        }
        WSC[1] = aux.WSC;  rm(aux.WSC)
        ## rest of days in agronomic year
        for (iday in 2:nday) {
          aux.WSC =  WSC[iday-1] + pr.year.WSC[iday] - ET0[iday]  # daily WSC is computed by adding any remaining moisture from the previous day to current day rainfall after removing PET (rainfall below large storm definition rnlg minus PET, and up to shc)
          if (!is.na(aux.WSC) & aux.WSC < 0) {
            aux.WSC = 0
          } else if (!is.na(aux.WSC) & aux.WSC > shc) {
            aux.WSC = shc  # anything above shc goes to runoff
          }
          WSC[iday] = aux.WSC;  rm(aux.WSC)
        }
        
        #############################################
        ## DETERMINATION OF THE AGRONOMIC SEASON (AS)
        #############################################
        ## first start date of the AS: dt_st_rnagsn
        ## number of false starts of the AS: nm_flst_rnagsn
        ## final start date of the AS: dt_fnst_rnagsn
        ## end date of the AS: dt_ed_rnagsn
        #############################################
        potSAG = which(WSC > 0.25*shc & WSC > 20)  # potential candidate dates for being the start of the Agronomic Season (AS)
        # potSAG = potSAG[potSAG > 124] # imposing 1-Nov cutoff
        
        if (length(potSAG) > 0) {  # there is at least one candidate date for being the start of the AS
          
          
          ## first start date of the AS
          dt_st_rnagsn[which(year == iyear)] = potSAG[1]
          
          ## number of false starts
          fSAG = potSAG[1]
          indW = seq(fSAG+1, fSAG+30) 
          bindWSC = WSC[indW]
          bindWSC[bindWSC > 0] = 1
          lv.spell = binSpell(bindWSC)
          nFS = 0
          
          while ((sum(lv.spell$len[na.omit(lv.spell$val == 0)] >= 5) >= 1)) {  # false start: this occurs when within the next 30 days, the WSC falls to zero and remains there for 5 or more consecutive days
            nFS = nFS + 1;
            aux = which(lv.spell$val == 0 & lv.spell$len >= 5)[1]
            
            if ((fSAG+sum(lv.spell$len[1:aux])+1) < max(potSAG, na.rm = T)) {
              fSAG = fSAG + sum(lv.spell$len[1:aux])+1
              
              fSAG = potSAG[potSAG >= fSAG]
              fSAG = fSAG[1]
              
              indW = seq(fSAG+1, fSAG+30)
              if (max(indW, na.rm = T) <= nday) {
                bindWSC = WSC[indW]
                bindWSC[bindWSC > 0] = 1
                
                lv.spell = binSpell(bindWSC)
              }
            } else {
              break
            }
          }
          nm_flst_rnagsn[which(year == iyear)] = nFS  # total number of false starts
          
          
          ## final start of the AS
          indW = seq(fSAG+1, fSAG+30)
          tryCatch({
            bindWSC = WSC[indW]
            bindWSC[bindWSC > 0] = 1
            lv.spell = binSpell(bindWSC)
            
            if (sum(lv.spell$len[lv.spell$val == 0] >= 1) >= 1) {  # asking for 1 or more consecutive days with WSC=0 to establish the end of the AS
              dt_fnst_rnagsn[which(year == iyear)] = NA  # the end of the AS does not occur
              dt_ed_rnagsn[which(year == iyear)] = NA  # the end of the AS does not occur
            } else {
              dt_fnst_rnagsn[which(year == iyear)] = fSAG
              aux = which(lv.spell$val == 1)
              
              auxbindWSC = WSC[seq(fSAG+sum(lv.spell$len[1:aux[length(aux)]]), length(WSC))]
              auxbindWSC[auxbindWSC > 0] = 1
              
              aux.lv.spell = binSpell(auxbindWSC)
              
              indEAS = which(aux.lv.spell$val == 0 & aux.lv.spell$len >= 5)  # asking for 1 or more consecutive days with WSC=0 to establish the end of the AS
              if (length(indEAS) > 0) {
                dt_ed_rnagsn[which(year == iyear)] = fSAG+sum(lv.spell$len[1:aux[length(aux)]])-1 + sum(aux.lv.spell$len[seq(1, indEAS[1]-1)])+1
              } else {
                dt_ed_rnagsn[which(year == iyear)] = nday
              }
            }
          }, error = function(x) {
            suppressWarnings(warning())
          })
          sAS = dt_fnst_rnagsn[which(year == iyear)]  # start of AS (simplified notation)
          eAS = dt_ed_rnagsn[which(year == iyear)]  # end of AS (simplified notation)
          
          ## length of AS
          lAS = eAS - sAS  
          dl_agsn[which(year == iyear)] = lAS 
          
          if (index.code %in% c("dt_st_rnagsn", "nm_flst_rnagsn", "dt_fnst_rnagsn", "dt_ed_rnagsn", "dl_agsn")) {
            out <- switch(index.code,
                          "dt_st_rnagsn" = dt_st_rnagsn,
                          "nm_flst_rnagsn" = nm_flst_rnagsn,
                          "dt_fnst_rnagsn" = dt_fnst_rnagsn,
                          "dt_ed_rnagsn" = dt_ed_rnagsn,
                          "dl_agsn" = dl_agsn)
          } else {
            ################
            ## indices in AS
            ################
            if (!is.na(lAS)) {
              
              pr.AS = pr.year[seq(sAS, eAS-1)]  # pr in AS
              tm.AS = tm.year[seq(sAS, eAS-1)]  # tm in AS
              tx.AS = tx.year[seq(sAS, eAS-1)]  # tx in AS
              tn.AS = tn.year[seq(sAS, eAS-1)]  # tn in AS
              
              out[which(year == iyear)] <- switch(index.code,
                            "dc_agsn" =  sum(pr.AS >= rndy, na.rm = T),  # number of rainy days in AS
                            "rn_agsn" = sum(pr.AS[pr.AS >= rndy], na.rm = T),  # total amount of rain in AS
                            "avrn_agsn" = mean(pr.AS[pr.AS >= rndy], na.rm = T),  # average amount of rain in rainy days in AS
                            "dc_rnlg_agsn" = sum(pr.AS >= rnlg, na.rm = T),  # number of days with "high" rain in AS
                            
                            ##############################
                            ## temperature-related indices
                            ##############################
                            "tm_agsn" = mean(tm.AS, na.rm = T),  # average daily tm in AS
                            "dc_txh_agsn" = sum(tx.AS >= txh, na.rm = T),  # number of "hot" tmax days in AS
                            "dc_tnh_agsn" = sum(tn.AS >= tnh, na.rm = T)  # number of "hot" tmin days in AS
              )
            }
          }
        } 
      }
    } else {
      message(sprintf("... year %d not complete ...", iyear))
    }
  }
  return(out)
}



#' #' @title Occurrence and length of binary spells
#' #' @description Computes ccurrence and length of binary spells)
#' #' @param v Vector
#' #' @author R. Manzanas
#' 
#' binSpell <-  function(v) {
#'   ix <- c(which(v[-length(v)] != v[-1]), length(v));
#'   out <- list()
#'   out$len <- diff(c(0, ix))
#'   out$val <- v[ix]
#'   return(out)
#' }

#' @title ET0
#' @description computes ET0, accroding to Hargreaves equation
#' @param lat Latitude (Integer)
#' @param dates Vector of dates 
#' @param tx Vector of maximum temperature data
#' @param tn Vector of minimum temperature data
#' @author R. Manzanas

computeET0 = function(lat, dates, tx, tn) {
  phi = lat*pi/180;  # latitude in radians; positive (negative) in the northern (sothern) hemisphere
  ref.J = cbind(dates[, 1], rep(1, dim(dates)[1]), rep(1, dim(dates)[1]))  # origin for Julian days
  J = as.numeric(as.Date(sprintf("%s-%s-%s", dates[,1], dates[,2], dates[,3])) -
                   as.Date(sprintf("%s-%s-%s", ref.J[,1], ref.J[,2], ref.J[,3])) + 1)  # Julian days (index with respect to origin) in agronomic year
  
  dr = 1 + (0.033*cos((2*pi/365)*J))
  delta = 0.409*sin((2*pi*J/365) - 1.39)
  ws = acos(-tan(phi)*tan(delta))
  auxRa = (24*60/pi)*0.0820*dr*((ws*sin(phi)*sin(delta)) + (cos(phi)*cos(delta)*sin(ws))) # radiation (M*J*m^-2*day^-1)
  ET0 = 0.0023*(((tx + tn) / 2) + 17.8)*((tx - tn)^0.5)*auxRa*0.408 # evapotranspiration (mm/day)
  return(ET0)
}
#######################
## call to the function
#######################
# load("/media/maialen/work/WORK/LOCAL/FAO_INDICES/data_example_mai.rda", verbose = T)
# AI = agroindexFAO(index.code = "dc_tnh_agsn", lat = args$lat, dates = args$dates, pr = args$pr, tx = args$tx, tn = args$tn)
# AI
