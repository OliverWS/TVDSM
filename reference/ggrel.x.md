# Calculate a relative x-coordinate in a ggplot

Given a ggplot object, this function returns an x-coordinate at a
relative position between the plot limits.

## Usage

``` r
ggrel.x(plt = last_plot(), x = 0.5)
```

## Arguments

- plt:

  A ggplot object. Defaults to the last plot.

- x:

  A numeric value between 0 and 1 representing the relative position.
  Default is 0.5.

## Value

A numeric value representing the x-coordinate.

## Examples

``` r
ggrel.x()
#> Error in ggplot_build(plt): could not find function "ggplot_build"
```
