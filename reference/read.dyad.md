# Read Dyad Data from File

Reads dyad data from a file and constructs a data frame with Timestamps
assigned based on the sample rate.

## Usage

``` r
read.dyad(
  f,
  sample_rate = 32,
  p1.name = NULL,
  p2.name = NULL,
  start = strptime("2016-01-01 0:0:0", "%Y-%m-%d %H:%M:%S"),
  readFunction = read.csv
)
```

## Arguments

- f:

  File path to the dyad data file.

- sample_rate:

  Sampling rate (default 32.0).

- p1.name:

  Optional name for participant 1 column (if NULL, first column name is
  used).

- p2.name:

  Optional name for participant 2 column (if NULL, second column name is
  used).

- start:

  Starting time as POSIXlt (default "2016-01-01 0:0:0").

- readFunction:

  Function to read the file (default read.csv).

## Value

A data frame with Timestamp and dyad data.
