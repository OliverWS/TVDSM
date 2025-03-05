# tvdsm: Tools for Dyadic Data Analysis

`tvdsm` is an R package designed to help you analyze dyadic data efficiently. This README serves as the main page for the pkgdown website and provides support documentation including quick examples and guidance on using the key function `analyzeDyad`.

More information about the package can be found at [https://oliverws.github.io/TVDSM/](https://oliverws.github.io/TVDSM/).

To read in-depth information about the TVDSM approach, please see [Saunders Wilder, O (2017)](https://doi.org/10.17760/D20316235).

## Quick Start

Below is a simple example of using `analyzeDyad` to analyze dyadic data:

```r
# Load the tvdsm package
library(tvdsm)


# Run the analysis using analyzeDyad
PersonA <- read.csv("PersonA.csv")
PersonB <- read.csv("PersonB.csv")
dyad <- as.dyad(PersonA, PersonB, cols = c("EDA"))

result <- analyzeDyad(dyad=dyad, window_size=60, lag=0, measure="EDA")


```

## What does analyzeDyad do?

The `analyzeDyad` function performs the following steps:

- **Data Processing:** Cleans and pre-processes the dyadic data for analysis.
- **Statistical Analysis:** Applies the chosen method (e.g., "default") to uncover patterns and relationships in the data.
- **Output Generation:** Produces a summary of the analysis including statistical measures, model fits, and diagnostic plots.

For example, the output typically includes:

- A **summary table** of fit indices and parameter estimates.
- **Diagnostic plots** for assessing model assumptions and fit.
- Additional **supporting metrics** to help interpret the analysis results.

## Learn More and Get Help

For more detailed documentation and support, please visit our [Documentation Site](https://oliverws.github.io/TVDSM/)

*tvdsm* is maintained by Oliver Saunders Wilder, PhD. Contributions and issues are welcome! 