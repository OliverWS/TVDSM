# Extract Coefficient with Significance Check

Extract Coefficient with Significance Check

## Usage

``` r
o.coef(mdl, c = 2, cutoff = 0.05, useD = F)
```

## Arguments

- mdl:

  Linear model object.

- c:

  Parameter index (default 2).

- cutoff:

  Significance cutoff (default 0.05).

- useD:

  Logical, if TRUE use standardized coefficient (default FALSE).

## Value

Coefficient value or 0 if not significant.
