# Fit dynamic system model for dyad data

Downsamples dyad data and fits a dynamic system model for both
participants.

## Usage

``` r
fitDSModelForDyad(
  dyad,
  type = 4,
  downsample = 1,
  lag = 0,
  x_mu = NULL,
  y_mu = NULL,
  verbose = F
)
```

## Arguments

- dyad:

  A data frame with dyad time series data.

- type:

  Model type (default 4).

- downsample:

  Downsampling factor (default 1).

- lag:

  Lag value (default 0).

- x_mu:

  Baseline for participant 1.

- y_mu:

  Baseline for participant 2.

- verbose:

  Logical, whether to output verbose messages.

## Value

A list containing model outputs and timing information.
