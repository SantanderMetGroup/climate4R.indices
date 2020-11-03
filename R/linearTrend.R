#     linearTrend. Computes linear trend slopes
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

#' @title Computes linear trend slopes
#' @description Computes linear trend slopes, their confidence 
#'  intervals and some other related statistics for timeseries using 
#'  Santer et al. (2008) method modified for dealing with missing data in
#'  time series.  
#' @param grid A climate4R predictor dataset
#' @param p A numeric value. Confidence level for the uncertainty interval (0<p<1).
#' @details ##  Uses a model y = a + b*x + e, where e is assumed to be an AR(1) process. 
#'  Applies ordinary least squares (OLS) regression. Estimates the trend slope's
#'  standard error, confidence interval, and p-value using an AR(1)-based 
#'  reduction in the number of degrees of freedom (DOF), as described by 
#'  Santer et al. (2008) paper in Int. J. of Climatology (S2008 hereinafter). 
#'  S2008 algorithm is modified here for the use with time series where some 
#'  data points are missing and also is augmented here by provisions for two 
#'  special cases: when the lag-1 autoregression coefficient is negative 
#'  (0 is used) and when reduced DOF number falls below 3 or cannot be 
#'  estimated at all ("irregularity" code % is raised and Inf or NaN are 
#'  returned for some parameters). 
#'  Algorithm modification to accomodate the presence of missing data is as 
#'  follows: (a) the slope and intercept of the regression line are computed 
#'  by applying the OLS only to the subset of points where data values are not 
#'  missing; (b) while S2008 compute the lag-1 autocorrelation coefficient 
#'  of data residuals e as a sample correlation coefficient between e(1:N-1) 
#'  and e(2:N), here it is computed for e(ind) and e(ind+1), where ind is the 
#'  subset of indices i in 1:N such that both e(i) and e(i+1) are available. 
#' 
#'  Note that the same formula for reduced number of degrees of freedom
#'  (DOFr) is used as in S2008; in the presence of missing data this
#'  formula underestimates DOFr and thus results in wider (more
#'  conservative) confidence intervals and larger p-values (compared to
#'  the cases without missing data).
#'  
#'  ##  Please note: 
#'  (1) the time grid is expected to be uniform, and 
#'  (2) missing data values contain NaN
#'  (3) if this function is called with 2 input parameters only, p is assumed 
#'      to be NaN, and NaN is returned for cinthw; other output parameters 
#'      are unaffected
#'  (4) this function uses an external function mklr.m for the OLS regression  
#' 
#'  Input: 
#'  x - vector of time values (a uniform grid) 
#'  y - vector of data values (NaN in place of missing values)
#'  p - confidence level for the uncertainty interval (0<p<1)
#' 
#'  Output:
#'  b      - estimated slope of the linear trend 
#'  cinthw - halfwidth of the confidence interval
#'  sig    - estimated standard error in b
#'  DOFr   - reduced number of DOF (a.k.a. effective sample size)
#'  rho    - lag-1 autocorrelation coefficient of data residuals 
#'              (w.r.t. trend line)
#'  pval   - p-value of the estimated b in the two-sided Student's t test for 
#'              the null hypothesis of no trend
#'  irrc   - "irregularity" code:
#'              irrc = 0, regular application of the algorithm
#' 	      irrc = 1, rho < 0 (rho=0 value is used in sig and cinthw 
#'                  calculations, but unmodified rho value is returned)
#'              irrc = 10, DOFr < 3; in this case results of the calculation 
#'                  are not recommended for use; when DOFr-->2-0, values of 
#'                  sig and cinthw tend to infinity and pval-->1 (Inf are 
#'                  returned for sig and cinthw and 1 for pval when DOFr <= 2)
#'              irrc = 100, rho could not be estimated (NaN are returned for 
#'                  rho, Inf for sig and cinthw, and 1 for pval)
#'              irrc = 1000, no calculation is done b/c number of available 
#'                  data points Na < 3
#'  N      - timeseries length 
#'  a      - intercept of the trend line
#'  Na     - number of available data points
#'  Nc     - sample size for calculating rho
#' 
#'  Detailed algorithm:
#' 
#'  0) output parameters are assigned missing values, corresponding to the 
#'       case of no calculation; set irrc=0; N and Na are determined;
#'       if Na<3, set irrc = 1000, and no calculation is done; otherwise, 
#'       vector of available data ya (length Na) and corresponding vector 
#'       xa are formed;
#'  1) Do the OLS regression of timeseries y on time x for the subset of 
#'       available values of y (ya) and corresponding subset of x (xa): 
#'       ya=a+b*xa; 
#'       returns estimated a, b, sb (= estimated standard error for b)
#'  2) compute data residuals e=y-a-b*x;
#'  3) compute rho (= lag-1 autocorrelation coefficient for residuals) as 
#'       the correlation coefficient between e(ind) and e(ind+1), where 
#'       ind is the subset of indices i in 1:N such that both e(i) and e(i+1) 
#'       are available; Nc = length(ind); if rho cannot be determined 
#'       (e.g., Nc < 2), set irrc = 100 and no further calculation is done. 
#'  4) provision for the negative rho: rhop = max(rho,0); if rho<0,  irrc = 1; 
#'  5) compute reduced DOF in the sample (a.k.a. effective sample size): 
#'            DOFr=Na*(1-rhop)/(1+rhop);
#'  6) provision for DOFr < 3 : irrc = 10;
#'  7) if DOFr>2, compute the following:
#' 
#'     (a) final standard error estimate for b 
#'            sig=sb*sqrt((N-2)/(DOFr-2))
#' 
#'     (b) p-value for b in the two-sided Student's t test:
#'            pval = 2*(1-tcdf(abs(b)/sig,DOFr-2))
#' 									     
#'     (c) half-width of the (p*100)% confidence interval for b  
#'            cinthw=sig*tinv(0.5+p/2,DOFr-2), 
#'         where tinv(alpha,n) is the (alpha*100)%-quantile of Student's 
#'         t distribution with n degrees of freedom.
#' @return Return a multigrid C4R object with the 'p', 'pval' and 'irrc' parameters stacked as variables.
#' 
#' @author J. Baño-Medina, S. Herrera
#' @export
#' @importFrom magrittr %<>% %>%
#' @importFrom transformeR redim
#' @importFrom stats pt qt
#' @examples
#' library(transformeR)
#' library(visualizeR)
#' require(climate4R.datasets)
#' # Regular Grid
#' x <- get(data("CFS_Iberia_tas")) 
#' r <- linearTrend(x,p = 0.9) 
#' b <- subsetGrid(r,var = "b")
#' spatialPlot(b, backdrop.theme = "coastline", main = "Estimated slope of the linear trend (CFS)")
#' 
#' # Irregular Grid
#' x <- get(data("VALUE_Iberia_tas")) # irregular Grid
#' r <- linearTrend(x,p = 0.9)
#' b <- subsetGrid(r,var = "b")
#' spatialPlot(b, backdrop.theme = "coastline", main = "Estimated slope of the linear trend (VALUE)")

linearTrend <- function(grid, p = 0.90){
  message("The slope of the linear trend is estimated based on the temporal resolution of the data. Please consider using the function aggregateGrid prior to the call of linearTrend to adequate the data to your resolution of interest.")
  grid <- if (isRegular(grid)) {
    redim(grid, var = TRUE)
  } else {
    redim(grid, var = TRUE, loc = TRUE)
  }

  if (getShape(grid,"var") > 1) stop("Multigrid objects are not accepted. Please consider using `subsetGrid` function to build single-variable C4R objects.")
  x <- 1:getShape(grid,"time")
  lapply(1:getShape(grid,"member"), FUN = function(z) {
    grid %<>% subsetGrid(members = z)
    b <- pval <- irrc <- climatology(grid)
    y <- if (isRegular(grid)){
      array3Dto2Dmat(grid$Data)
    } else {
      grid$Data
    }
    aux <- sapply(1:ncol(y), FUN = function(zz) { 
      y <- y[,zz]
      
      # If called with 2 input parameters only, p is assumed to be NaN: 
      # if nargin < 3, p=NaN; end
      # (in this case, p = NULL, NA is returned for cinthw)
      
      ## Step 0
      irrc <- 0
      b <- NA
      sig <- Inf
      pval <- 1.0
      cinthw <- Inf
      DOFr <- 0
      rho <- NA
      Nc <- NA
      a <- NA
      
      N <- length(x)
      ia <- which(!is.na(y))
      Na <- length(ia)
      if (Na < 3) {
        irrc <- irrc +  1000
      } else {
        ya <- y[ia]
        xa <- x[ia]
        
        ## Step 1
        ## calling the regular OLS (mklr() is an external function) 
        ## for available data
        ols.model <- mklr(xa,ya)
        a <- ols.model$r0
        sa <- ols.model$sig0
        b <- ols.model$r1
        sb <- ols.model$sig1
        
        ## Step 2
        ## computing data residuals w.r.t. the OLS trend line:
        ea <- ya - a - b*xa
        
        e <- array(data = NA, dim = N)
        e[ia] <- ea
        
        ## Step 3
        ## lag-1 autocorrelation coefficient
        EE <- matrix(data = NA, nrow = N-1, ncol = 2)
        EE[,1] <- e[c(1:N-1)]
        EE[,2] <- e[c(2:N)]
        
        ic <- which(!is.na(apply(EE, FUN = sum, MARGIN = 1)))
        Nc <- length(ic)
        if (Nc < 2) {
          irrc <- irrc +  100
        } else {
          Ec <- EE[ic,]
          rho <- cor(Ec[,1],Ec[,2])
          if (is.na(rho)){
            irrc <- irrc + 100
          } else {
            ## Step 4
            rhop=max(c(rho,0));
            if (rho < 0) {
              irrc=irrc+1;
            }
            ## Step 5
            ## Computing reduced number of DOF
            DOFr <- Na*(1-rhop)/(1+rhop)
            ## Step 6 
            if (DOFr < 3) {
              irrc <- irrc+10
            }
            ## Step 7
            if (DOFr > 2) {
              ## Adjust estimated standard error in the trend slope 
              sig <- sb*sqrt((Na-2)/(DOFr-2))
              pval <- 2*(1-pt(abs(b)/sig,DOFr-2))
              ## half-width of the (p*100)% confidence interval for b:  
              cinthw <- sig*qt(0.5+p/2,DOFr-2)
            }
          }
        }
      } 
      c(b,pval,irrc)
    })
    if (isRegular(b)) {
      b$Data <- mat2Dto3Darray(aux[1,,drop = FALSE],x = b$xyCoords$x, y = b$xyCoords$y)
      pval$Data <- mat2Dto3Darray(aux[2,,drop = FALSE],x = pval$xyCoords$x, y = pval$xyCoords$y)
      irrc$Data <- mat2Dto3Darray(aux[3,,drop = FALSE],x = irrc$xyCoords$x, y = irrc$xyCoords$y)
    } else {
      b$Data <- aux[1,,drop = FALSE]; attr(b$Data,"dimensions") <- c("time","loc")
      pval$Data <- aux[1,,drop = FALSE]; attr(pval$Data,"dimensions") <- c("time","loc")
      irrc$Data <- aux[1,,drop = FALSE]; attr(irrc$Data,"dimensions") <- c("time","loc")
    }
    out <- makeMultiGrid(b,pval,irrc)
    out$Variable$varName <- c("b","pval","irrc")
    return(out)
  }) %>% bindGrid(dimension = "member") 
}

#' @title Error estimates of the regression coefficients
#' @description Univariate linear regression with error estimates for regression coefficients.
#' @param x A vector of predictors. 
#' @param y A vector of observations.
#' @details The mathematics behind the calculation of the error estimates of the 
#' regression coefficients are based on the Dover's book by A.A.Sveshnikov 
#' "Problems in Probability Theory,...", p.326, --AK, Feb 27  2004.
#' Performs univariate linear regression of y on x:
#'         y = r0 + r1*x
#' @return Returns a list with the estimates (coefficients: r0, r1, 
#' and their standard errors: sig0 and sig1).
#' @keywords internal
#' @author J. Baño-Medina, S. Herrera
mklr <- function(x,y){
  n <- length(x)
  s0 <- n
  s1 <- sum(x, na.rm = TRUE)
  s2 <- sum(x^2, na.rm = TRUE)
  v0 <- sum(y, na.rm = TRUE)
  v1 <- sum(y*x, na.rm = TRUE)
  v2 <- sum(y*x^2, na.rm = TRUE)
  
  r0 <- (s2*v0-s1*v1)/(s2*s0-s1^2)
  r1 <- (-s1*v0+s0*v1)/(s2*s0-s1^2)
  
  e <- y-(r0+r1*x)
  Smin <- sum(e^2, na.rm = TRUE)
  
  sig0 <- sqrt(s2/(s2*s0-s1^2)*Smin/(n-2))
  sig1 <- sqrt(s0/(s2*s0-s1^2)*Smin/(n-2))
  
  return(list(r0 = r0, sig0 = sig0, r1 = r1, sig1 = sig1))
}

