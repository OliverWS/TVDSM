# Annotate a ggplot with relative coordinates

This function adds an annotation (default text) at a relative position
on the ggplot.

## Usage

``` r
annotate.relative(plt = last_plot(), geom = "text", x = 0.5, y = 0.1, ...)
```

## Arguments

- plt:

  A ggplot object. Defaults to the last plot.

- geom:

  Character string specifying the type of geom to add; defaults to
  'text'.

- x:

  Relative x-position (0 to 1). Default is 0.5.

- y:

  Relative y-position (0 to 1). Default is 0.1.

- ...:

  Additional arguments passed to the annotate function.

## Value

The annotation layer added to the plot.

## Examples

``` r
annotate.relative(plt, label = 'Hello')
#> Error in ggplot_build(plt): could not find function "ggplot_build"
```
