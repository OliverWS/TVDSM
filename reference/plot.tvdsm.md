# Plot TVDSM model results

Generates a combined plot displaying TVDSM model parameters and
R-squared values over time.

## Usage

``` r
# S3 method for class 'tvdsm'
plot(
  data,
  xname = "x",
  yname = "y",
  title = paste(xname, "&", yname),
  save = F,
  f = "plot.pdf",
  use.delta.rsquared = T,
  by.condition = T,
  autoscale = F,
  plotParams = T,
  grayscale = F,
  rawData = NULL,
  returnSeparatePlots = F,
  fontFamily = "Helvetica"
)
```

## Arguments

- data:

  A list of TVDSM model objects.

- xname:

  Friendly name for participant 1.

- yname:

  Friendly name for participant 2.

- title:

  Plot title.

- save:

  Logical, whether to save the plot.

- f:

  Filename to save the plot.

- use.delta.rsquared:

  Logical, whether to use delta R-squared values.

- by.condition:

  Logical, whether to group by condition.

- autoscale:

  Logical, whether to auto-scale the plot.

- plotParams:

  Logical, whether to include model parameter plots.

- grayscale:

  Logical, whether to use a grayscale color scheme.

- rawData:

  Optional raw data for plotting.

- returnSeparatePlots:

  Logical, if TRUE returns separate plots.

- fontFamily:

  Font family for text.

## Value

A ggplot object or a list of ggplot objects.
