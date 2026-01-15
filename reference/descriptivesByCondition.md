# Descriptives By Condition

Descriptives By Condition

## Usage

``` r
descriptivesByCondition(
  filename,
  codes,
  col = "EDA",
  title = filename,
  codeFunc = read.RTCodes,
  startCol = "Start.Time",
  endCol = "End.Time",
  conditionCol = "Condition",
  timeformat = "%Y-%m-%d %H:%M:%S"
)
```

## Arguments

- filename:

  Path to the EDA data file.

- codes:

  Codes data (either file path or dataframe).

- col:

  Column for analysis (default "EDA").

- title:

  Title for the plot.

- codeFunc:

  Function to read codes, default is read.RTCodes.

- startCol:

  Start time column name, default "Start.Time".

- endCol:

  End time column name, default "End.Time".

- conditionCol:

  Condition column name, default "Condition".

- timeformat:

  Time format string.

## Value

Data frame of descriptive statistics.
