# Plot dyadic data using ggplot2

This function reshapes data for two variables (dyad) and plots them as
lines with distinct colors.

## Usage

``` r
plotDyad(
  d,
  title = "",
  ylabel = "EDA",
  legend_label = "Participant",
  grayscale = FALSE,
  ...
)
```

## Arguments

- d:

  A data frame with at least a 'Timestamp' column and two other
  variables to compare.

- title:

  Plot title.

- ylabel:

  Label for the y-axis. Default is 'EDA'.

- legend_label:

  Label for the legend. Default is 'Participant'.

- grayscale:

  Logical. If TRUE, the plot is rendered in grayscale.

- ...:

  Additional parameters passed to ggplot.

## Value

A ggplot object.

## Examples

``` r
plotDyad(mydata, title = 'Dyad Plot')
#> Error: object 'mydata' not found
```
