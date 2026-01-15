# Analyze a group of participant data files

Performs pairwise TVDSM analysis on a group of participant data files.

## Usage

``` r
analyzeGroup(
  files,
  window_size = 60,
  window_step = 10,
  readFunction = read.eda,
  scanDirs = T,
  ...
)
```

## Arguments

- files:

  A vector of file paths or a directory path containing CSV files.

- window_size:

  Window size in seconds for TVDSM analysis.

- window_step:

  Window step size in seconds.

- readFunction:

  Function to read each data file (default read.eda).

- scanDirs:

  Logical, if TRUE, treats files as directory and scans for CSV files.

- ...:

  Additional arguments passed to the analysis functions.

## Value

A list containing analysis results for each participant pair.
