# ggCaterpillar: Plot random effects as QQ-plots or caterpillar dotplots.

This function creates either a normal QQ-plot or a caterpillar-style
dotplot for random effects.

## Usage

``` r
ggCaterpillar(re, QQ = TRUE, likeDotplot = TRUE)
```

## Arguments

- re:

  A list or matrix of random effects.

- QQ:

  Logical flag; if TRUE a QQ-plot is generated, otherwise a caterpillar
  plot is used (default: TRUE).

- likeDotplot:

  Logical flag; if TRUE the plot is styled like a dotplot (default:
  TRUE).

## Value

A list of ggplot objects.
