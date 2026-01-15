# Combine participant data into a group dyad

Merges overlapping time series data from multiple participants into a
single data frame.

## Usage

``` r
as.group(
  pn,
  cols = c("EDA"),
  norm = F,
  na.interpolate = T,
  interpolation.method = na.spline
)
```

## Arguments

- pn:

  A list of participant data frames.

- cols:

  A vector of column names to include (default "EDA").

- norm:

  Logical, whether to normalize data.

- na.interpolate:

  Logical, if TRUE interpolates missing values.

- interpolation.method:

  Function used for interpolation (default na.spline).

## Value

A combined data frame of group time series data.
