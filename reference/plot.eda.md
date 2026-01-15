# Plot EDA data with optional accelerometer data

This function creates a line plot for EDA data and, optionally, an
additional plot for accelerometer data.

## Usage

``` r
# S3 method for class 'eda'
plot(data, title = "EDA Data", show_acc = FALSE, show_raw_acc = FALSE)
```

## Arguments

- data:

  A data frame containing at least columns 'Timestamp' and 'EDA'. Also
  expected are 'X', 'Y', 'Z' if accelerometer is desired.

- title:

  Character string for the plot title. Default is 'EDA Data'.

- show_acc:

  Logical. If TRUE, an additional plot for acceleration is shown.

- show_raw_acc:

  Logical. If TRUE, raw accelerometer data (X, Y, Z) is plotted;
  otherwise, the magnitude (Movement) is computed and plotted.

## Value

A ggplot object or a cowplot object if accelerometer data is added.

## Examples

``` r
plot.eda(mydata, title = 'Subject 1 EDA', show_acc = TRUE)
#> Error in ggplot(data = data, aes(x = Timestamp, y = EDA)): could not find function "ggplot"
```
