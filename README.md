## Description

A package for causal meta-analysis estimation of a binary response and binary treatment using aggregated data, based on Berenfeld et al. (https://arxiv.org/abs/2505.20168). The output is the aggregated estimate and its standard deviation, the package allows for simple and fast computation of different measures such as Risk Difference (RD), Risk Ratio (RR), Odds Ratio (OR) and Survival Ratio (SR).

The package consists of the function `camea`. This function can take the measures Risk Difference, Risk Ratio and log-Risk Ratio, Odds Ratio and log-Odds Ratio, Survival Ratio and log-Survival Ratio.

The data can be input in the following formats: vectors, data.frame, matrix or sparse matrix.

The output of the function is a list of two elements: the first element is a dataframe that has, for each study, its related estimate and standard error; the second element of the list is the aggregated causal-meta analysis estimate and its standard error.

If `plot = TRUE` the function `camea` prints a forestplot that displays the estimate of interest for each study and the aggregated estimate(s). All estimates are shown with their associated 95%-level CI.

If `random.effects = TRUE`, the output and plot will include both the causal meta-analysis estimate and the random-effects meta-analysis estimate.

If `log.scale = TRUE`, the meta-analysis will be performed on the logarithmic transformation of the chosen measure.

## How to install

The latest release of the package can be installed through CRAN (soon):

```{r, eval = FALSE}
install.packages("CaMeA")
```

The development version can be installed from github

```{r, eval = FALSE}
devtools::install_github("cberenfeld/camea", build_vignettes = TRUE)
```

Another installation possibility is to clone the repo, and then within the r-package folder run

```{r, eval = FALSE}
Rscript build_package.R
```

## Usage Example

```{r}
require(camea)

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
result <- camea(measure = measure, ai = treated_events, n1i = treated_total,
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
result <- camea(measure = measure, ai = treated_events, bi = treated_negatives,
                    ci = control_events, di = control_negatives, data = dat, log.scale = TRUE)


```

## Plot

If `plot = TRUE` it is possible to print the forestplot of the meta-analysis.

```{r my-forest-plot, fig.width=9, fig.height=6, out.width='100%'}
# apply the same function but print the plot
result <- camea(measure = measure, ai = treated_events, n1i = treated_total,
                    ci = control_events, n2i = control_total, plot = TRUE)
```

## Use of metafor package

If the selected effect measure is Risk Difference (RD), Risk Ratio (RR) or Odds Ratio (OR) and the option `random.effects = TRUE` is set, the output and the forest plot will display two summary estimates: the primary causal meta-analysis estimate and the conventional random-effects estimate.

The random-effects estimate is calculated via the package `metafor` by W. Viechtbauer (2010), which is a well known tool for conducting meta-analysis on R. 

It is important to notice that the package `metafor` only provides estimates of the Risk Ratio and Odds Ratio on the logarithmic scale. Our package automatically converts these log-scale estimates back to the natural scale by applying the exponential function to both the point estimate and the bounds of the confidence interval.

## References

Berenfeld, C., Boughdiri, A., Colnet, B., van Amsterdam, W. A. C., Bellet, A., Khellaf, R., Scornet, E., & Josse, J. (2025). Causal Meta-Analysis: Rethinking the Foundations of Evidence-Based Medicine. arXiv:2505.20168.

Viechtbauer, W. (2010). Conducting Meta-Analyses in R with the metafor Package. Journal of Statistical Software, 36(3), 1–48. https://doi.org/10.18637/jss.v036.i03
