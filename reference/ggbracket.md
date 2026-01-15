# Add bracket annotations for statistical comparisons on a ggplot

This function adds bracket annotations and p-value labels to a ggplot,
based on t-tests between groups.

## Usage

``` r
ggbracket(
  plt,
  pad = 0.1,
  offset = 0,
  sig.cutoff = 0.05,
  ci.func = mean_cl_normal
)
```

## Arguments

- plt:

  A ggplot object with a grouping variable mapped to the x-axis.

- pad:

  Numeric, padding value for the annotations.

- offset:

  Numeric, unused in this implementation.

- sig.cutoff:

  Numeric, significance level threshold. Comparisons with p-values below
  this are annotated.

- ci.func:

  Function to compute confidence intervals, default is mean_cl_normal.

## Value

A ggplot object with added annotation segments and text.

## Examples

``` r
ggbracket(my_plot)
#> Error in ggplot_build(plt): could not find function "ggplot_build"
```
