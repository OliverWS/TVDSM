% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpersonal.R
\name{as.group}
\alias{as.group}
\title{Combine participant data into a group dyad}
\usage{
as.group(
  pn,
  cols = c("EDA"),
  norm = F,
  na.interpolate = T,
  interpolation.method = na.spline
)
}
\arguments{
\item{pn}{A list of participant data frames.}

\item{cols}{A vector of column names to include (default "EDA").}

\item{norm}{Logical, whether to normalize data.}

\item{na.interpolate}{Logical, if TRUE interpolates missing values.}

\item{interpolation.method}{Function used for interpolation (default na.spline).}
}
\value{
A combined data frame of group time series data.
}
\description{
Merges overlapping time series data from multiple participants into a single data frame.
}
