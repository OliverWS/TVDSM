# Read iButton Data

Reads data from an iButton CSV file skipping header lines.

## Usage

``` r
read.ibutton(path, header = T)
```

## Arguments

- path:

  File path to the CSV file.

- header:

  Logical indicating if file has header (default TRUE).

## Value

A data frame with Timestamp and Temperature columns.
