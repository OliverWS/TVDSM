# Save interpersonal analysis data to a CSV file

Exports the summary data (and optionally condition-based data) from a
dyad analysis to a CSV file.

## Usage

``` r
saveInterpersonalData(
  dyadData,
  outputFilename = NULL,
  I3.only = F,
  aname = NULL,
  bname = NULL,
  by.condition = F
)
```

## Arguments

- dyadData:

  A list containing TVDSM dyad analysis output.

- outputFilename:

  Optional filename for the output CSV.

- I3.only:

  Logical, if TRUE, only saves the I3 metric.

- aname:

  Friendly name for participant 1.

- bname:

  Friendly name for participant 2.

- by.condition:

  Logical, if TRUE, includes condition information.

## Value

The data frame that was saved.
