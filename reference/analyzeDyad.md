# Analyze dyad data for TVDSM analysis

Processes two participant datasets or a dyad object, performs TVDSM
analysis, generates plots, and returns model summaries.

## Usage

``` r
analyzeDyad(
  f1 = "",
  f2 = "",
  dyad = c(),
  xname = f1,
  yname = f2,
  norm = F,
  window_size = 60 * 5,
  window_step = window_size,
  start = "",
  end = "",
  func = fitDSModelForDyad,
  na.rm = T,
  simulate = F,
  dname = paste(xname, yname, sep = "+"),
  lag = 0,
  noPlots = F,
  plotParams = T,
  pltTitle = paste(dname, "(", "Lag", "=", lag, ")"),
  measure = "EDA",
  type = 4,
  downsample = 1,
  x.baseline = NA,
  y.baseline = NA,
  verbose = F,
  codes = "",
  useRealTime = F,
  incTSD = T,
  ...
)
```

## Arguments

- f1:

  File path for participant 1.

- f2:

  File path for participant 2.

- dyad:

  Optional dyad object.

- xname:

  Friendly name for participant 1.

- yname:

  Friendly name for participant 2.

- norm:

  Logical, whether to standardize data.

- window_size:

  Window size in seconds.

- window_step:

  Step size in seconds.

- start:

  Start time for analysis.

- end:

  End time for analysis.

- func:

  Function for model fitting (default fitDSModelForDyad).

- na.rm:

  Logical.

- simulate:

  Logical, whether to simulate a dyad.

- dname:

  Friendly name for the dyad.

- lag:

  Lag value for analysis.

- noPlots:

  Logical, if TRUE, no plots are produced.

- plotParams:

  Logical, whether to include model parameter plots.

- pltTitle:

  Title for the plots.

- measure:

  Measurement column name.

- type:

  Model type.

- downsample:

  Downsampling factor.

- x.baseline:

  Baseline for participant 1.

- y.baseline:

  Baseline for participant 2.

- verbose:

  Logical, whether to be verbose.

- codes:

  File path for condition codes.

- useRealTime:

  Logical, whether to interpret condition codes in real time.

- incTSD:

  Logical, whether to include time series descriptives.

- ...:

  Additional arguments.

## Value

A list containing TVDSM model outputs, summary, and plots.
