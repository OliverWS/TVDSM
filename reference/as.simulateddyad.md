# Create Simulated Dyad Data

Simulates a dyad using two participant data frames.

## Usage

``` r
as.simulateddyad(p1, p2, cols = c("EDA"), norm = F, verbose = F)
```

## Arguments

- p1:

  Data frame for participant 1.

- p2:

  Data frame for participant 2.

- cols:

  Columns to include (default "EDA").

- norm:

  Logical indicating whether to normalize data (default FALSE).

- verbose:

  Logical for verbose output (default FALSE).

## Value

A data frame representing the simulated dyad.
