# Arrange multiple ggplot objects in a grid

This function accepts multiple ggplot objects (either via ... or a
provided list) and displays them in a grid layout.

## Usage

``` r
multiplot(
  ...,
  plotlist = NULL,
  file,
  cols = 1,
  layout = NULL,
  heights = c(2, 1)
)
```

## Arguments

- ...:

  ggplot objects passed as individual arguments.

- plotlist:

  An optional list of ggplot objects.

- file:

  (unused) parameter reserved for file output.

- cols:

  Number of columns in the layout.

- layout:

  A matrix specifying the layout. If provided, cols is ignored.

- heights:

  A numeric vector specifying the relative heights of rows.

## Value

None. The function prints the plots directly.

## Examples

``` r
multiplot(p1, p2, cols = 2)
#> Loading required package: grid
#> Error: object 'p1' not found
```
