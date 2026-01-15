# Create a histogram with density and normality annotation

This function generates a density plot for a numeric vector. If groups
are provided, it overlays densities by group.

## Usage

``` r
gghist(x, groups = NULL)
```

## Arguments

- x:

  Numeric vector to plot.

- groups:

  Optional grouping variable. If provided, densities are colored by
  group.

## Value

A ggplot object with the density plot and normality test annotation (if
applicable).

## Examples

``` r
gghist(mydata$variable)
#> Error: object 'mydata' not found
```
