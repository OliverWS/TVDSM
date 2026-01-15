# Read RTCodes

## Usage

``` r
read.RTCodes(
  path = "ConditionTimes.csv",
  startCol = "Start.Time",
  endCol = "End.Time",
  fmt = "%y-%m-%d %H:%M:%S"
)
```

## Arguments

- path:

  File path (default "ConditionTimes.csv").

- startCol:

  Column name for start time.

- endCol:

  Column name for end time.

- fmt:

  Time format (default "

A data frame with RTCodes. Reads RTCodes from a CSV or Excel file.
