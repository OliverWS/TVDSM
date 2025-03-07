% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpersonal.R
\name{analyzeGroup}
\alias{analyzeGroup}
\title{Analyze a group of participant data files}
\usage{
analyzeGroup(
  files,
  window_size = 60,
  window_step = 10,
  readFunction = read.eda,
  scanDirs = T,
  ...
)
}
\arguments{
\item{files}{A vector of file paths or a directory path containing CSV files.}

\item{window_size}{Window size in seconds for TVDSM analysis.}

\item{window_step}{Window step size in seconds.}

\item{readFunction}{Function to read each data file (default read.eda).}

\item{scanDirs}{Logical, if TRUE, treats files as directory and scans for CSV files.}

\item{...}{Additional arguments passed to the analysis functions.}
}
\value{
A list containing analysis results for each participant pair.
}
\description{
Performs pairwise TVDSM analysis on a group of participant data files.
}
