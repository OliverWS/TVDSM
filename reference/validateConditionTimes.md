# Validate condition times within data range

## Usage

``` r
validateConditionTimes(d, codes, timeformat = "%Y-%m-%d %H:%M:%S")
```

## Arguments

- d:

  A data frame with a Timestamp column.

- codes:

  A data frame with condition timing (with Start.Time, End.Time, and
  Condition).

- timeformat:

  A string specifying the timestamp format (default "

Stops with an error if a condition time is out-of-range. Checks each
condition's start and end times to confirm they fall within the dyad's
timestamp range.
