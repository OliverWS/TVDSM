# Simulate a Dyadic Time Series

This function simulates a dyadic time series for two interacting persons
using self- and co-regulation parameters.

## Usage

``` r
simulateDyad(
  duration = 600,
  fs = 32,
  selfReg.coef = 0.5,
  coReg.coef = 1,
  interaction.coef = 0,
  lag = 0,
  mu = 1,
  sd = 2,
  trend = 0,
  sr.ratio = 0.95
)
```

## Arguments

- duration:

  Duration of simulation in seconds (default: 600).

- fs:

  Sampling frequency (default: 32).

- selfReg.coef:

  Self-regulation coefficient (default: 0.5).

- coReg.coef:

  Co-regulation coefficient (default: 1.0).

- interaction.coef:

  Interaction coefficient (default: 0).

- lag:

  Lag to apply (default: 0).

- mu:

  Mean value for baseline signal (default: 1).

- sd:

  Standard deviation for noise (default: 2).

- trend:

  Linear trend added to the signal (default: 0).

- sr.ratio:

  Ratio of signal to noise (default: 0.95).

## Value

A list with elements: plt (ggplot object), data (simulated data frame),
and eq (model equation).
