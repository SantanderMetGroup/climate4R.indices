% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indicesFAO_tier1.R
\name{avg}
\alias{avg}
\title{Function to compute the average value of a daily time-series}
\usage{
avg(
  tm,
  dates,
  year = NULL,
  year.start = NULL,
  year.end = NULL,
  pnan = 25,
  lat = NULL
)
}
\arguments{
\item{tm}{Vector with data (e.g. daily mean temperature)}

\item{dates}{Matrix containing the full range of dates corresponding to "tm" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)}

\item{year}{Vector with years of interest (e.g. 1990:1995)}

\item{year.start}{Vector of dates [in "YYYY-MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season)}

\item{year.end}{Vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)}

\item{pnan}{Any year with a percentage of NA data above "pnan" will be ignored}

\item{lat}{Latitude (NULL) to indicate that latitude information is not used.}
}
\value{
Average value (per year)
}
\description{
Function to compute the average value of a daily time-series
}
\author{
R. Manzanas
}
