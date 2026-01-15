# Compare Simulated Dyads

Compares true dyad pairings with simulated (false) dyad pairings using
statistical tests and visualization.

## Usage

``` r
compareSimulatedDyads(pairs, title = "True vs. Simulated Dyads")
```

## Arguments

- pairs:

  A list with elements 'truePairs' and 'falsePairs'.

- title:

  Title for the comparison plot (default: "True vs. Simulated Dyads").

## Value

A list containing the compiled data frame, plot object, and t-test
results.
