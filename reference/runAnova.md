# Run ANOVA on model set

Runs ANOVA comparisons on delta R-squared values from a list of models
using a specified condition.

## Usage

``` r
runAnova(mdls, key = "Condition")
```

## Arguments

- mdls:

  A list of model objects.

- key:

  The key to extract the condition variable (default "Condition").

## Value

A fitted lmer model object.
