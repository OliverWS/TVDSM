# Calculate a relative y-coordinate in a ggplot

Given a ggplot object, this function returns a y-coordinate at a
relative position between the plot limits.

## Usage

``` r
ggrel.y(plt = last_plot(), y = 0.5)
```

## Arguments

- plt:

  A ggplot object. Defaults to the last plot.

- y:

  A numeric value between 0 and 1 representing the relative position.
  Default is 0.5.

## Value

A numeric value representing the y-coordinate.

## Examples

``` r
ggrel.y()
#> Error in ggplot_build(plt): could not find function "ggplot_build"
```
