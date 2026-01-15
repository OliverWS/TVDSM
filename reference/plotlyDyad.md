# Plot dyadic data using Plotly

This function reshapes data for two variables and creates an interactive
Plotly line plot.

## Usage

``` r
plotlyDyad(
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

  A data frame with at least a 'Timestamp' column and two variables.

- title:

  Plot title.

- ylabel:

  Label for the y-axis. Default is 'EDA'.

- legend_label:

  Label for the legend. Default is 'Participant'.

- grayscale:

  Logical. If TRUE, applies a grayscale color scheme.

- ...:

  Additional parameters.

## Value

A Plotly plot object.

## Examples

``` r
plotlyDyad(mydata, title = 'Dyad Interactive Plot')
#> Error in melt(d, id.vars = c("Timestamp"), variable_name = legend_label): could not find function "melt"
```
