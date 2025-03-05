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
PersonA <- read.csvdata("PersonA.csv")
PersonB <- read.csvdata("PersonB.csv")
dyad <- as.dyad(PersonA, PersonB, cols = c("EDA"))

result <- analyzeDyad(dyad=dyad, window_size=60, lag=0,xname = "PersonA",yname="PersonB",plotParams = "raw")

head(result$summary)
#             timestamp            name               Start                 End dx.r.squared dy.r.squared   x.selfreg   x.coreg x.interaction   y.selfreg   y.coreg y.interaction
# 1 1970-01-01 00:00:00 PersonA+PersonB 1970-01-01 00:00:00 1970-01-01 00:01:00    0.8609977    0.7678614  0.00000000 0.7595278            NA  0.00000000 1.3113947            NA
# 2 1970-01-01 00:01:00 PersonA+PersonB 1970-01-01 00:01:00 1970-01-01 00:02:00    0.8972776    0.6761051  0.02346823 0.7968545            NA -0.05435490 1.2197551            NA
# 3 1970-01-01 00:02:00 PersonA+PersonB 1970-01-01 00:02:00 1970-01-01 00:03:00    0.8029374    0.7413561  0.00000000 0.5092651            NA  0.00000000 1.9555250            NA
# 4 1970-01-01 00:03:00 PersonA+PersonB 1970-01-01 00:03:00 1970-01-01 00:04:00    0.8384866    0.7174741  0.00000000 0.3875892            NA  0.00000000 2.5704267            NA
# 5 1970-01-01 00:04:00 PersonA+PersonB 1970-01-01 00:04:00 1970-01-01 00:05:00    0.7615926    0.5771764  0.00000000 0.2225588            NA -0.05550322 4.3719645            NA
# 6 1970-01-01 00:05:00 PersonA+PersonB 1970-01-01 00:05:00 1970-01-01 00:06:00    0.6726747    0.4688025 -0.06390888 1.7151281            NA -0.14787452 0.4763664            NA

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