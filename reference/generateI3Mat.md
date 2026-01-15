# Generate an I3 matrix for interdependence analysis

Constructs a 3D matrix of I3 values across participants and time windows
and optionally saves it as JSON.

## Usage

``` r
generateI3Mat(data, outputFile = NULL)
```

## Arguments

- data:

  A nested list of dyad analysis results.

- outputFile:

  Optional filename to save the JSON output.

## Value

A list containing the legend, I3 matrix, and timestamps.
