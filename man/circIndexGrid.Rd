% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/circIndexGrid.R
\name{circIndexGrid}
\alias{circIndexGrid}
\title{Calculation of circulation indices of grids}
\usage{
circIndexGrid(
  zg = NULL,
  z = NULL,
  sst = NULL,
  psl = NULL,
  index.code,
  season = NULL,
  base = NULL,
  ref = NULL,
  match = "spatial",
  n.pcs = 10,
  rot = TRUE,
  centers = NULL,
  center.point = c(-5, 55),
  members = NULL
)
}
\arguments{
\item{zg}{A grid or multimember grid object of geopotential height.}

\item{z}{A grid or multimember grid object of geopotential.}

\item{sst}{A grid or multimember grid object of sea surface temperature.}

\item{psl}{A grid or multimember grid object of sea level pressure.}

\item{index.code}{Circulation index (or vector of indices with same input variables) to be computed. See \code{circIndexShow()} for details.}

\item{season}{Selected month(s) for the calculation. Default: NULL (i.e. as input grid).}

\item{base}{Baseline grid to be substracted for the calculation of anomalies. Default: NULL. See \code{?scaleGrid}.}

\item{ref}{Reference grid to be added for the calculation of anomalies. Default: NULL. See \code{?scaleGrid}.}

\item{match}{Character string with the criterion to be used for the detection of CPC indices. Options "spatial" or "temporal". Default: "spatial". See details.}

\item{n.pcs}{Integer vector with the number of EOFs to be retained for the CPC indices. Default: 10. See details.}

\item{rot}{Logical. Should VARIMAX-Rotation be performed? Default: TRUE. This argument is only relevant for CPC indices. See details.}

\item{centers}{Integer value indicating the number of clusters, \strong{k}, or center points. See details.}

\item{center.point}{A two value vector that must include lon and lat from a location that will work as center point for the Lamb WT.
See details.}

\item{members}{Select number of members. Default: NULL (all members of the grid).}
}
\value{
A list of circulation indices (and members, if applicable) with:
\itemize{
\item index: vector with the time series of the teleconnection/circulation index.
\item pattern: matrix with the spatial pattern of the teleconnection.
\item dates and coordinates as list attributes.
\item further arguments related to the CPC indices, such as the corresponding (r)EOF and (temporal or spatial, depending on \code{'match'}) correlation with the original index.
}
}
\description{
Calculate circulation indices of grids or multimember grids.
}
\details{
\code{\link{circIndexShow}} displays on screen the full list of available circulation indices and their codes.
Several indices can be calculated at the same time, as long as they depend on the same input variable(s) and spatial domain. CPC and ENSO indices are calculated on a monthly basis. Therefore a temporal aggregation is performed if input data is daily.
Results for the desired months in \code{season} are provided, but it is recommended to have full series as input, since CPC and ENSO indices use a moving window for the calculation.

\strong{CPC indices}

Either \code{z} or \code{zg} are valid input variables. CPC indices are obtained, by default, as the first 10 Varimax-rotated EOFs, as explained in \url{https://www.cpc.ncep.noaa.gov/data/teledoc/telecontents.shtml}. The core of this function is \code{stats::prcomp} including Varimax rotation.
The rotated EOFs are obtained from the monthly standardized anomalies of geopotential or geopotential height at 500hPa, with a 3-month moving window.
The argument \code{match} is used to assign each rEOF to a circulation index and pattern. Matching is based on 'temporal' or 'spatial' correction of the CPC original (NCEP Reanalysis-based) indices.

\strong{ENSO indices}

The calculation of ENSO indices is based on \url{https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni}, consisting of SST anomalies, using a moving window of different size for each index.

\strong{WT indices}

(to be added)
}
\examples{
\dontrun{
data(NCEP_hgt500_2001_2010)
cpc <- circIndexGrid(zg=NCEP_hgt500_2001_2010, index.code = c("NAO", "EA","PNA"), season=1)
data(ERAInterim_sst_1981_2010)
nino <- circIndexGrid(sst=ERAInterim_sst_1981_2010, index.code = "NINO3.4")
data(NCEP_slp_2001_2010)
wt <- circIndexGrid(psl = NCEP_slp_2001_2010, index.code = "WT.LAMB")
}
}
\author{
A. Casanueva
}
