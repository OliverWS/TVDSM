# Plot State Space Graph

Plot State Space Graph

## Usage

``` r
statespacegraph(
  xt,
  yt,
  norm = T,
  step = 0.25,
  title = "State Space Plot",
  xlabel = "X Data",
  ylabel = "Y Data",
  type = 3
)
```

## Arguments

- xt:

  Numeric vector for x values.

- yt:

  Numeric vector for y values.

- norm:

  Logical to normalize data (default TRUE).

- step:

  Step size (default 0.25).

- title:

  Plot title.

- xlabel:

  Label for x-axis.

- ylabel:

  Label for y-axis.

- type:

  Model type (default 3).

## Value

A ggplot object representing the state space graph.
