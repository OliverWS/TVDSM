% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpersonal.R
\name{saveDataForGraphing}
\alias{saveDataForGraphing}
\title{Save data for graphing from TVDSM analysis}
\usage{
saveDataForGraphing(mdls, xfname = "Dyad_X.csv", yfname = "Dyad_X.csv")
}
\arguments{
\item{mdls}{A list containing TVDSM model outputs and summary data.}

\item{xfname}{Filename for saving participant X's R-squared data.}

\item{yfname}{Filename for saving participant Y's R-squared data.}
}
\value{
None.
}
\description{
Exports the R-squared summary metrics of the TVDSM analysis to CSV files.
}
