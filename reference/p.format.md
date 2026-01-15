# Format p-value for display

Formats a p-value into a string with appropriate formatting. If the
p-value is less than 0.001, it omits the equal sign.

## Usage

``` r
p.format(p)
```

## Arguments

- p:

  Numeric, the p-value to format.

## Value

A character string representing the formatted p-value.

## Examples

``` r
p.format(0.0005)
#> [1] "p <0.001"
```
