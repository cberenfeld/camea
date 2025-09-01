## Description

A package for causal meta-analysis estimation of a binary response and binary treatment. The output is the aggregated estimate and its standard deviation, the package allows for simple and fast computation of different measures such as risk difference, risk ratio, odds ratio and survival ratio.

The package consists of a the function "causalmeta". This function can take the measures Risk Difference, Risk Ratio and log-Risk Ratio, Odds Ratio and log-Odds Ratio, Survival Ratio and log-Survival Ratio.

The data can be input in the following formats: vectors, data.frame, matrix or sparse matrix.

The output of the function is a list of two elements: the first element is a dataframe that has, for each study, its related estimate and standard error; the second element of the list is the aggregated causal-meta analysis estimate and its standard error.

If "plot = TRUE" the function "causalmeta" prints a forestplot that displays the estimate of interest for each study, the aggregated causal-meta estimate and the random-effects model estimate. All estimates are shown with their associated 95%-level CI.

## How to install

The latest release of the package can be installed through CRAN (soon):

```{r, eval = FALSE}
install.packages("causalmeta")
```

The development version can be installed from github

```{r, eval = FALSE}
devtools::install_github("cberenfeld/causal_meta")
```

Another installation possibility is to clone the repo, and then within the r-package folder run

```{r, eval = FALSE}
Rscript build_package.R
```

## Usage Example

```{r}
require(causalmeta)

# Example 1: With vectors of total numbers and events in treated and control

# Generate data
n <- 10

treated_total <- sample.int(1000, n)
control_total <- sample.int(1000, n)
treated_events <- rbinom(n,treated_total,0.5)
control_events <- rbinom(n,control_total,0.5)

# Choose measure: RD
measure <- "RD"

# Apply function
result <- causalmeta(measure = measure, ai = treated_events, n1i = treated_total,
                    ci = control_events, n2i = control_total)

# Example 2: With contingency tables entries in data.frame format

# Generate data
n <- 10

treated_total <- sample.int(1000, n)
control_total <- sample.int(1000, n)
treated_events <- rbinom(n,treated_total,0.5)
control_events <- rbinom(n,control_total,0.5)
treated_negatives <- treated_total - treated_events
control_negatives <- control_total - control_events

dat <- data.frame(treated_events,control_events,treated_negatives,control_negatives)

# Choose measure: log-OR
measure <- "OR"

# Apply function
result <- causalmeta(measure = measure, ai = treated_events, bi = treated_negatives,
                    ci = control_events, di = control_negatives, data = dat, log.scale = TRUE)


```

## Plot

If "plot = TRUE" it is possible to print the forestplot of the meta-analysis.

```{r my-forest-plot, fig.width=9, fig.height=6, out.width='100%'}
# apply the same function but print the plot
result <- causalmeta(measure = measure, ai = treated_events, n1i = treated_total,
                    ci = control_events, n2i = control_total, plot = TRUE)
```

## Use of metafor package

If the measure used is Risk Difference, Risk Ratio or Odds Ratio, the forest plot, along with the causal meta-analysis estimate, will also have the estimate made with the random-effects model. The latter is estimated via the package "metafor" by Viechtbauer W. . The package "metafor" provides only estimates of the Risk Ratio and Odds Ratio in (natural) logarithmic scale, for the natural scale our package takes the exponential of the estimates and of the error term (to compute the confidence intervals).

## References

Berenfeld, C., Boughdiri, A., Colnet, B., van Amsterdam, W. A. C., Bellet, A., Khellaf, R., Scornet, E., & Josse, J. (2025). Causal Meta-Analysis: Rethinking the Foundations of Evidence-Based Medicine. arXiv:2505.20168.

Viechtbauer, W. (2010). Conducting Meta-Analyses in R with the metafor Package. Journal of Statistical Software, 36(3), 1–48. https://doi.org/10.18637/jss.v036.i03
