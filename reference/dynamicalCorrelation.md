# Compute dynamical correlation using bootstrap methods

Calculates the dynamical correlation between participant data using a
bootstrap approach.

## Usage

``` r
dynamicalCorrelation(
  a_files,
  b_files,
  read_func = read.eda,
  cols = c("EDA"),
  sr = 32,
  randPairs = F
)
```

## Arguments

- a_files:

  A vector of file paths for participant A.

- b_files:

  A vector of file paths for participant B.

- read_func:

  Function to read data files (default read.eda).

- cols:

  A vector of measurement columns (default "EDA").

- sr:

  Sampling rate (default 32).

- randPairs:

  Logical, if TRUE, uses random pairing.

## Value

Bootstrapped dynamical correlation test results.
