# Create a scatter plot with a linear fit and correlation annotation

This function generates a scatter plot for two variables, adds a linear
regression line, and annotates the plot with the correlation coefficient
and p-value.

## Usage

``` r
scatter_plot(x, y, data)
```

## Arguments

- x:

  Character string representing the name of the x variable in the data
  frame.

- y:

  Character string representing the name of the y variable in the data
  frame.

- data:

  A data frame containing the variables.

## Value

A ggplot object with the scatter plot and annotations.

## Examples

``` r
scatter_plot('height', 'weight', mydata)
#> Error in ggplot(data, aes_string(x = x, y = y)): could not find function "ggplot"
```
