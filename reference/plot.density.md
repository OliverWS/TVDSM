# Plot density contour of beta coefficients

Creates a contour plot from a 2D kernel density estimate of beta
coefficients.

## Usage

``` r
# S3 method for class 'density'
plot(
  data,
  group = "Dyad",
  xname = "Particpant 1 Beta",
  yname = "Participant 2 Beta"
)
```

## Arguments

- data:

  A data frame containing beta coefficient values.

- group:

  Group label for the plot.

- xname:

  Label for the x-axis (default: "Particpant 1 Beta").

- yname:

  Label for the y-axis (default: "Participant 2 Beta").

## Value

A plot of the density contours.
