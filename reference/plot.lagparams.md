# Plot lag parameters over time

Plots the lags and R-squared values over time using ggplot2.

## Usage

``` r
# S3 method for class 'lagparams'
plot(
  data,
  xname = "x",
  yname = "y",
  title = paste(xname, "&", yname),
  save = F,
  f = "plot.pdf",
  use.delta.rsquared = T,
  by.condition = T,
  autoscale = T
)
```

## Arguments

- data:

  A data frame containing timestamps and parameter values.

- xname:

  Label for the x variable.

- yname:

  Label for the y variable.

- title:

  The plot title.

- save:

  Logical, whether to save the plot.

- f:

  Filename to save the plot (default "plot.pdf").

- use.delta.rsquared:

  Logical, whether to use delta R-squared values.

- by.condition:

  Logical, whether to group by condition.

- autoscale:

  Logical, whether to auto-scale the plot.

## Value

A combined ggplot object.
