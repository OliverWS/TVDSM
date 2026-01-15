# Run False Pairings Analysis

Runs TVDSM analysis on true and false pairings from lists of participant
signals.

## Usage

``` r
runFalsePairings(
  participantAList,
  participantBList,
  window_size = 60,
  window_step = window_size,
  downsample = 1,
  lag = 0,
  type = 3
)
```

## Arguments

- participantAList:

  List of signals for person A.

- participantBList:

  List of signals for person B.

- window_size:

  Window size for analysis (default: 60).

- window_step:

  Step size between windows (default: window_size).

- downsample:

  Downsample factor (default: 1).

- lag:

  Lag to be applied (default: 0).

- type:

  Analysis type specifier (default: 3).

## Value

A list containing results for true pairs and false pairs.
