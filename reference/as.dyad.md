# Convert Two Participant Data into a Dyad

Aligns and merges two participant data frames into a dyadic data frame
based on overlapping timestamps.

## Usage

``` r
as.dyad(
  p1,
  p2,
  cols = c("EDA"),
  norm = F,
  verbose = F,
  na.interpolate = T,
  interpolation.method = na.spline
)
```

## Arguments

- p1:

  Data frame for participant 1.

- p2:

  Data frame for participant 2.

- cols:

  Columns to include (default "EDA").

- norm:

  Logical, if TRUE normalize the data (default FALSE).

- verbose:

  Logical for verbose output (default FALSE).

- na.interpolate:

  Logical, if TRUE interpolate NA values (default TRUE).

- interpolation.method:

  Function to interpolate NAs (default na.spline).

## Value

A dyadic data frame with aligned timestamps.
