# Compare Simulated Signals

Compares simulated signals using partner simulation and dyadic analysis,
then plots the results.

## Usage

``` r
compareSim(
  x.signal = sampleSignal(),
  fs = 0.1,
  coReg.coef = 0.1,
  selfReg.coef = 0.05
)
```

## Arguments

- x.signal:

  Signal for person A (default: sampleSignal()).

- fs:

  Sampling frequency (default: 0.1).

- coReg.coef:

  Co-regulation coefficient (default: 0.1).

- selfReg.coef:

  Self-regulation coefficient (default: 0.05).

## Value

A plot grid comparing simulation outputs.
