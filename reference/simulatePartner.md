# Simulate a Partner Time Series

This function simulates a partner's time series based on an input signal
and regulation parameters.

## Usage

``` r
simulatePartner(
  x.signal,
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

- x.signal:

  Numeric vector representing the input signal.

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

  Mean value (default: 1).

- sd:

  Standard deviation (default: 2).

- trend:

  Linear trend to add (default: 0).

- sr.ratio:

  Ratio of signal to noise (default: 0.95).

## Value

A list with elements: data (simulated partner data frame), plt (ggplot
object), and eq (model equation).
