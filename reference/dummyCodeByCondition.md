# Dummy Code By Condition

Dummy Code By Condition

## Usage

``` r
dummyCodeByCondition(
  filename,
  codes,
  outputFile = paste(filename, "_DummyCoded.csv", sep = ""),
  saveFile = F,
  codeFunc = read.RTCodes,
  startCol = "Start.Time",
  endCol = "End.Time",
  conditionCol = "Condition",
  timeformat = "%Y-%m-%d %H:%M:%S"
)
```

## Arguments

- filename:

  Filename for the EDA data.

- codes:

  Either a path to a file with codes or a dataframe with codes.

- outputFile:

  Output file name (default is filename appended with
  "\_DummyCoded.csv").

- saveFile:

  Logical indicating whether to save the dummy coded file.

- codeFunc:

  Function to read codes, default is read.RTCodes.

- startCol:

  Column name for start time, default "Start.Time".

- endCol:

  Column name for end time, default "End.Time".

- conditionCol:

  Column name for condition, default "Condition".

- timeformat:

  Time format for conversion.

## Value

Data frame with dummy coded data.
