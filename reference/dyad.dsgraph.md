# Plot Dyadic Graph

Plot Dyadic Graph

## Usage

``` r
dyad.dsgraph(
  p1,
  p2,
  norm = T,
  epoch = 1,
  step = 0.25,
  title = "State Space Plot",
  type = 2,
  measure = "EDA"
)
```

## Arguments

- p1:

  Data for participant 1.

- p2:

  Data for participant 2.

- norm:

  Logical indicating whether to normalize the data (default TRUE).

- epoch:

  Epoch size (default 1).

- step:

  Step size (default 0.25).

- title:

  Plot title.

- type:

  Model type (default 2).

- measure:

  Column name for measurement (default "EDA").

## Value

A plot of the model params
