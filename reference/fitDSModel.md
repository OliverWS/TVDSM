# Fit a dynamic system model

Fits a dynamic system model to time series data using systemfit and SUR.

## Usage

``` r
fitDSModel(
  x,
  y,
  x_mu = mean(x, na.rm = T),
  y_mu = mean(y, na.rm = T),
  type = 2,
  step = 0.25,
  p.value = 0.01,
  verbose = F,
  lag = 0
)
```

## Arguments

- x:

  Numeric vector representing x values.

- y:

  Numeric vector representing y values.

- x_mu:

  Baseline value for x; default is mean(x, na.rm=TRUE).

- y_mu:

  Baseline value for y; default is mean(y, na.rm=TRUE).

- type:

  Model type (default 2).

- step:

  Step size parameter (default 0.25).

- p.value:

  Significance level for the likelihood ratio test.

- verbose:

  Logical, whether to print verbose output.

- lag:

  Lag value to apply.

## Value

A list containing the fitted models, R-squared values, and model
coefficients.
