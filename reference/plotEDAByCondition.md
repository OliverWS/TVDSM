# Plot EDA by Condition

This function reads EDA data and overlays condition codes as shaded
rectangles on the plot.

## Usage

``` r
plotEDAByCondition(eda, codes, title = "")
```

## Arguments

- eda:

  A source for EDA data; expected to be processed by read.eda.

- codes:

  A data frame containing condition codes with columns 'Start Time',
  'End Time', and 'Condition'.

- title:

  Plot title.

## Value

A ggplot object representing the EDA plot with condition overlays.

## Examples

``` r
plotEDAByCondition(eda_file, condition_codes, title = "EDA by Condition")
#> Error: object 'eda_file' not found
```
