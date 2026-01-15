# Plot a metric with grouping and error bars

This function generates a summary plot for a given metric partitioned by
a group, including mean points and error bars.

## Usage

``` r
plotMetric(data, metric, group, sig.cutoff = 0.05)
```

## Arguments

- data:

  A data frame containing the data to plot.

- metric:

  Character string naming the metric column to plot.

- group:

  Character string naming the grouping variable.

- sig.cutoff:

  Numeric significance cutoff for brackets. Default is 0.05.

## Value

A ggplot object with the metric plot and statistical brackets.

## Examples

``` r
plotMetric(mydata, 'response', 'group')
#> Error in ggplot(data = data, aes_string(x = group, y = metric, col = group)): could not find function "ggplot"
```
