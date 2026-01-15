# Create Dyad Data Frame

Constructs a dyadic data frame from EDA data with generated timestamps.

## Usage

``` r
makeDyad(
  edaData,
  sample_rate = 0.1,
  start = strptime("2015-01-01 0:0:0", "%Y-%m-%d %H:%M:%S"),
  norm = T
)
```

## Arguments

- edaData:

  Matrix or data frame with two columns representing signals.

- sample_rate:

  Sampling rate (default: 0.1).

- start:

  Start time as POSIXlt (default: 2015-01-01 0:0:0).

- norm:

  Logical indicating whether to normalize data (default: TRUE).

## Value

A data frame with columns: Timestamp, Person A, and Person B.
