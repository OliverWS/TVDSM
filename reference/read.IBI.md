# Read IBI Data from CSV

This function reads interbeat interval (IBI) data from a CSV file.

## Usage

``` r
read.IBI(
  file,
  startTime = "1970-01-01 0:0:0",
  start = strptime(startTime, "%Y-%m-%d %H:%M:%S", tz = "")
)
```

## Arguments

- file:

  File path to the CSV file.

- startTime:

  Starting time as a string (default "1970-01-01 0:0:0").

- start:

  Starting time as a POSIXlt object (default parsed from startTime).

## Value

A data frame with Timestamp and IBI columns.
