# Read Condition Codes

## Usage

``` r
read.Codes(
  path = "ConditionTimes.csv",
  startCol = "Start.Time",
  endCol = "End.Time",
  conditionCol = "Condition",
  fmt = "%H:%M:%S",
  ...
)
```

## Arguments

- path:

  File path (default "ConditionTimes.csv").

- startCol:

  Column name for start time.

- endCol:

  Column name for end time.

- conditionCol:

  Column name for condition.

- fmt:

  Format string for time (default "

  ...Additional parameters.

A data frame with condition codes. Reads condition codes from a CSV or
Excel file.
