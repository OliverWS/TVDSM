# Read IBI Time Series Data

Reads IBI data and constructs a time series with computed heart rate
(HR).

## Usage

``` r
read.IBI.TS(
  file,
  startTime = "1970-01-01 0:0:0",
  start = strptime(startTime, "%Y-%m-%d %H:%M:%S", tz = "")
)
```

## Arguments

- file:

  File path to the CSV file.

- startTime:

  Starting time string (default "1970-01-01 0:0:0").

- start:

  Starting time as POSIXlt (default parsed from startTime).

## Value

A data frame with Timestamp, IBI, and HR columns.
