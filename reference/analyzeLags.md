# Analyze model performance over varying lags

Iterates over a range of lag values to fit models and determine optimal
lag based on R-squared.

## Usage

``` r
analyzeLags(
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
  measure = "EDA",
  downsample = 1,
  minLag = 0,
  maxLag = 5,
  noPlots = F,
  relativeToLag = -1,
  type = 4,
  plotParams = T,
  pltTitle = NULL
)
```

## Arguments

- f1:

  File path for participant 1.

- f2:

  File path for participant 2.

- dyad:

  Optional dyad object (if provided, f1 and f2 are ignored).

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

  Start time.

- end:

  End time.

- func:

  Function to use for fitting the model.

- na.rm:

  Logical, whether to remove NA values.

- simulate:

  Logical, whether to simulate dyad data.

- dname:

  Friendly name for the dyad.

- measure:

  Measurement column (e.g., "EDA").

- downsample:

  Downsampling factor.

- minLag:

  Minimum lag value.

- maxLag:

  Maximum lag value.

- noPlots:

  Logical, if TRUE no plots are generated.

- relativeToLag:

  Reference lag index (default -1).

- type:

  Model type.

- plotParams:

  Logical, whether to plot model parameters.

- pltTitle:

  Title for the plots.

## Value

A list with optimal lag summary and model outputs.
