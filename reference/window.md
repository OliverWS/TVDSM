# Apply Moving Window Operation

Apply Moving Window Operation

## Usage

``` r
window(x, window_size, window_overlap = 0, FUN, na.rm = T)
```

## Arguments

- x:

  Numeric vector.

- window_size:

  Size of each window.

- window_overlap:

  Overlap between windows.

- FUN:

  Function to apply on each window.

- na.rm:

  Logical indicating whether to remove NAs (default TRUE).

## Value

A vector of computed values.
