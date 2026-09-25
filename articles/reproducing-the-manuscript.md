# Reproducing the Manuscript's Analyses with POEM

This vignette walks through **each analysis in the manuscript** “Power
Enhancement in High-Dimensional Heterogeneous Mediation Analysis” (in
press at the *Journal of the American Statistical Association*) and
reproduces it with POEM. It follows the manuscript’s own structure: the
Monte Carlo simulation study (its Section 5) and the empirical data
analysis (its Section 6). The companion vignette *Power Enhancement for
High-Dimensional Mediation: A Guided Tour* introduces the method itself;
here the focus is faithful replication.

**Two things to keep in mind while reading.**

1.  *Scope of the comparison.* The manuscript benchmarks the
    power-enhanced tests (PE-HDMM for linear outcomes, PE-HDGMM for
    generalized outcomes) against several competitors: HDMM, HILMA, and
    GlobalTest. POEM implements the proposed PE tests and the **HDMM**
    benchmark (the total-indirect-effect Wald test that the PE component
    augments). Every POEM result therefore reports the central contrast
    of the manuscript’s figures and tables, **PE versus HDMM**; HILMA
    and GlobalTest are separate methods outside this package.

2.  *Run size.* The manuscript’s simulations use `n = 300`, `p = 500`,
    and **1,000 replications** per design point. Reproducing that here
    would take far too long for a vignette, so the live code uses small
    `n`, `p`, and `n_rep`; the curves have the right *shape* but are
    noisier. Each simulation section also shows the exact manuscript
    settings in a `eval = FALSE` block you can run yourself.

``` r

set.seed(113)
```

## Part I. Monte Carlo simulation (manuscript Section 5)

The simulation study asks two questions of every test. Under the global
null of no active mediator (`c1 = 0`) the rejection rate estimates the
**Type I error rate**, which should sit near the nominal `alpha = 0.05`.
Away from the null (`c1` $`\neq 0`$) the rejection rate estimates
**power**.
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md)
sweeps `c1` and returns both the benchmark (`rejection_hdmm`) and
power-enhanced (`rejection_pe`) rejection rates at each grid point.

A small helper draws a manuscript-style rejection-rate figure for one
design:

``` r

plot_rejection <- function(pc, main) {
  plot(pc$c1, pc$rejection_pe, type = "b", pch = 19, ylim = c(0, 1),
       xlab = expression(c[1] * "  (signal strength)"),
       ylab = "empirical rejection rate", main = main)
  lines(pc$c1, pc$rejection_hdmm, type = "b", pch = 1, lty = 2)
  abline(h = 0.05, col = "grey60", lty = 3)
  legend("right", c("PE (proposed)", "HDMM (benchmark)"),
         pch = c(19, 1), lty = c(1, 2), bty = "n")
}
```

### Section 5.1 Linear mediation models

The manuscript studies a univariate exposure with `n = 300` continuous
outcomes and `p = 500` candidate mediators whose noise has an
autoregressive `0.5^|i-j|` covariance. The mediator-exposure coefficient
is `Gamma_x = c1 * tau` with `tau = (0.1, 0.2, 0.3, 0.4, 0.5, 0, ...)`,
the direct effect is `c2 = 0.5`, and two outcome-mediator patterns are
contrasted:

- **Setting (i), homogeneous:**
  `alpha_m = (1, 0.8, 0.6, 0.4, 0.2, 0, ...)`, so all active mediators
  push the same way and the total indirect effect is `beta = 0.7 * c1`
  (nonzero). This is the regime where ordinary tests already work; the
  question is whether PE *adds* power.
- **Setting (ii), contrasting:**
  `alpha_m = (1, -0.5, 0.4, -0.3, 0, ...)`, four active mediators whose
  effects cancel so the total indirect effect is `beta = 0` for every
  `c1`. This is the regime that defeats a total-indirect-effect test,
  and where PE is designed to rescue power.

POEM maps these to `pattern = "homogeneous"` and
`pattern = "contrasting"` in
[`simulate_mediation_data()`](https://yelleknek.github.io/POEM/reference/simulate_mediation_data.md)
and
[`pe_power_curve()`](https://yelleknek.github.io/POEM/reference/pe_power_curve.md).
The reduced live run (the manuscript’s Figure for linear models):

``` r

set.seed(113)
lm_homo <- pe_power_curve(n = 150, p = 50, outcome = "continuous",
                          pattern = "homogeneous",
                          c1_grid = c(0, 0.25, 0.5, 0.75, 1),
                          c2 = 0.5, n_rep = 25)
lm_cont <- pe_power_curve(n = 150, p = 50, outcome = "continuous",
                          pattern = "contrasting",
                          c1_grid = c(0, 0.25, 0.5, 0.75, 1),
                          c2 = 0.5, n_rep = 25)

op <- par(mfrow = c(1, 2))
plot_rejection(lm_homo, "(i) homogeneous")
plot_rejection(lm_cont, "(ii) contrasting")
```

![Reduced reproduction of the manuscript's linear-model rejection-rate
figure for the homogeneous and contrasting
settings](reproducing-the-manuscript_files/figure-html/lm-curves-1.png)

``` r

par(op)
```

Both panels reproduce the manuscript’s qualitative findings. At `c1 = 0`
both tests sit near `0.05` (correct size). Under the **homogeneous**
setting both tests gain power as the signal grows, with PE a step ahead
of HDMM. Under the **contrasting** setting HDMM stays flat near `0.05`
forever, it is blind to canceling effects, while PE climbs toward one.
That gap is the whole point of the power enhancement.

To reproduce the manuscript’s figure at its actual resolution (including
the negative half of the `c1` grid, which the symmetric design fills
in), run:

``` r

c1_grid <- seq(-1, 1, by = 0.1)
lm_homo_full <- pe_power_curve(n = 300, p = 500, outcome = "continuous",
                               pattern = "homogeneous", c1_grid = c1_grid,
                               c2 = 0.5, n_rep = 1000, cores = 4)
lm_cont_full <- pe_power_curve(n = 300, p = 500, outcome = "continuous",
                               pattern = "contrasting", c1_grid = c1_grid,
                               c2 = 0.5, n_rep = 1000, cores = 4)
```

[`pe_simulation_study()`](https://yelleknek.github.io/POEM/reference/pe_simulation_study.md)
is a convenience wrapper that runs several patterns at once and stacks
the results, which is how the manuscript’s figure panels are produced:

``` r

study <- pe_simulation_study(n = 300, p = 500, outcome = "continuous",
                             patterns = c("homogeneous", "contrasting"),
                             c1_grid = seq(-1, 1, by = 0.1), n_rep = 1000)
```

#### The real-data-motivated heterogeneous setting (supplement)

The supplement adds a third, harder pattern: a setting calibrated to the
DNA-methylation case study of Guo et al. (2022), with **eleven active
mediators of mixed sign** that neither all agree (homogeneous) nor
exactly cancel (contrasting).
[`simulate_guo_mediation()`](https://yelleknek.github.io/POEM/reference/simulate_guo_mediation.md)
generates it from the calibration constants shipped as
\[`guo_calibration`\]; the total indirect effect is
`guo_calibration$beta_per_c1 * c1`, about `-1.597 * c1`. (The manuscript
reports `-1.5977`, computed from the full-precision Guo coefficients;
the package ships those coefficients rounded to the three decimals
printed in the supplement, which give `-1.597`.)

``` r

d <- simulate_guo_mediation(c1 = 1, seed = 113)
dim(d$M)                              # 85 x 1008, the case-study dimensions
#> [1]   85 1008
d$active_mediators                    # the 11 active loci
#>  [1]  1  2  3  4  5  6  7  8  9 10 11
round(d$beta, 4)                      # beta_per_c1 * c1, i.e. about -1.597
#> [1] -1.597
```

Because `p = 1008` greatly exceeds `n`, the penalized fit needs a wider
tuning grid than the small-example default. With it, the power-enhanced
screen recovers the mixed-sign mediators and improves on the benchmark
Wald test:

``` r

dn <- simulate_guo_mediation(c1 = 1, n = 150, seed = 113)
fit <- pe_mediation(dn$X, dn$Y, dn$M, outcome = "continuous",
                    lambda_grid = seq(0.1, 10, length.out = 50))
fit[fit$term %in% c("pval_hdmm", "pval_pe", "n_active_mediators"), ]
#>  term               value   
#>  pval_hdmm          < 0.0001
#>  pval_pe            < 0.0001
#>  n_active_mediators 6       
#> 
#> Outcome model: continuous (linear)
#> Active mediators identified (6): 1, 4, 5, 6, 9, 10
```

The manuscript’s figure sweeps `c1` over `c(0, +/- 0.1, ..., +/- 1)` at
the case-study size `n = 85`, with and without the nine calibrated
confounders, averaged over 1,000 replications:

``` r

rej <- function(c1, confounders) {
  ps <- replicate(1000, {
    d <- simulate_guo_mediation(c1 = c1, confounders = confounders)
    f <- pe_mediation(d$X, d$Y, d$M, Z = d$Z, outcome = "continuous",
                      lambda_grid = seq(0.1, 10, length.out = 50))
    c(f$value[f$term == "pval_hdmm"], f$value[f$term == "pval_pe"])
  })
  rowMeans(ps <= 0.05)
}
sapply(seq(-1, 1, by = 0.1), rej, confounders = FALSE)
```

### Section 5.2 Mediation models with generalized outcomes

The same construction is repeated for non-continuous outcomes through a
link function: a **logistic** model for binary outcomes (PE-HDGMM versus
HDGMM) and a **Poisson** model for counts. Again `n = 300`, `p = 500`,
with a homogeneous and a contrasting pattern apiece; the manuscript
fixes the direct effect at `c2 = 1` for the logistic models and
`c2 = 0.4` for the Poisson models.

Here is the dramatic **contrasting** case for each outcome type (the
panel where PE most clearly rescues power), at reduced size:

``` r

set.seed(113)
logit_cont <- pe_power_curve(n = 180, p = 50, outcome = "binary",
                             pattern = "contrasting",
                             c1_grid = c(0, 0.25, 0.5, 0.75, 1),
                             c2 = 1, n_rep = 20)
pois_cont <- pe_power_curve(n = 180, p = 50, outcome = "count",
                            pattern = "contrasting",
                            c1_grid = c(0, 0.25, 0.5, 0.75, 1),
                            c2 = 0.4, n_rep = 20)

op <- par(mfrow = c(1, 2))
plot_rejection(logit_cont, "Logistic, contrasting")
plot_rejection(pois_cont, "Poisson, contrasting")
```

![Reduced reproduction of the manuscript's logistic and Poisson
rejection-rate figures under the contrasting
setting](reproducing-the-manuscript_files/figure-html/glm-curves-1.png)

``` r

par(op)
```

As in the manuscript, the benchmark HDGMM has little power against
contrasting effects while PE-HDGMM gains substantial power, all while
holding the size at `c1 = 0`. The full manuscript designs, including the
homogeneous panels:

``` r

c1_grid <- seq(-1, 1, by = 0.1)
# Logistic: c2 = 1
for (pat in c("homogeneous", "contrasting"))
  pe_power_curve(n = 300, p = 500, outcome = "binary", pattern = pat,
                 c1_grid = c1_grid, c2 = 1, n_rep = 1000, cores = 4)
# Poisson: c2 = 0.4
for (pat in c("homogeneous", "contrasting"))
  pe_power_curve(n = 300, p = 500, outcome = "count", pattern = pat,
                 c1_grid = c1_grid, c2 = 0.4, n_rep = 1000, cores = 4)
```

### Section 5.3 Identifying individual mediators (FWER / FDR)

Beyond the global test, the manuscript’s supplement studies how well the
PE machinery **identifies which** individual mediators are active, and
whether it controls the family-wise error rate (FWER) or false discovery
rate (FDR). POEM exposes this through the `method` argument
(`"Bonferroni"` for FWER, `"BH"` or `"BY"` for FDR) and the
`report_all_methods` switch, which reports all three side by side as the
supplement does.

For a *single* fit,
[`pe_selection()`](https://yelleknek.github.io/POEM/reference/pe_selection.md)
reports the active set under each method side by side, together with the
power-enhanced statistic and p-value each one produces (the screened set
feeds `J_m`, so the methods can give different global p-values, exactly
the PE / PE_BH / PE_BY columns of the manuscript’s extended tables):

``` r

set.seed(113)
d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
                             pattern = "contrasting", c1 = 1)
d$active_mediators                                  # ground truth
#> [1] 1 2 3 4

fit_all <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                        report_all_methods = TRUE)
pe_selection(fit_all)                               # per-method sets + p-values
#>       method n_active      j_pe   stat_pe       pval_pe active_mediators
#> 1 Bonferroni        3  864.1774  864.4981 5.124789e-190          2, 3, 4
#> 2         BH        4 1105.1162 1105.4370 2.173658e-242      2, 3, 4, 33
#> 3         BY        3  864.1774  864.4981 5.124789e-190          2, 3, 4
```

To reproduce the supplement’s identification *study* (its empirical
FWER, FDR, precision, and recall over many replications),
[`pe_identification_study()`](https://yelleknek.github.io/POEM/reference/pe_identification_study.md)
scores each method’s selections against the known active set across a
grid of `c1`. At `c1 = 0` there are no active mediators, so the `fwer`
column is the probability of selecting anything (the empirical FWER
under the null) and `recall` is `NA`:

``` r

set.seed(113)
pe_identification_study(n = 150, p = 60, outcome = "continuous",
                        pattern = "contrasting", c1_grid = c(0, 0.5, 1),
                        n_rep = 30)
#>  c1  method     fwer    fdr     precision recall n_valid
#>  0   Bonferroni 0.03333 0.03333 0         <NA>   30     
#>  0   BH         0.03333 0.03333 0         <NA>   30     
#>  0   BY         0.03333 0.03333 0         <NA>   30     
#>  0.5 Bonferroni 0       0       1         0.15   30     
#>  0.5 BH         0       0       1         0.1667 30     
#>  0.5 BY         0       0       1         0.1333 30     
#>  1   Bonferroni 0.03333 0.01111 0.9885    0.4917 30     
#>  1   BH         0.2333  0.08611 0.9109    0.575  30     
#>  1   BY         0.06667 0.01944 0.9792    0.4917 30     
#> 
#> Outcome model: continuous
```

Consistent with the manuscript’s report of “close-to-zero empirical FWER
and FDR, high precision, but medium-to-low recall,” the conservative
Bonferroni rule recovers the strongest active mediators with very few
false positives, and the FDR rules typically recover somewhat more. The
manuscript-scale study raises `n`, `p`, and `n_rep` to 300 / 500 / 1000:

``` r

pe_identification_study(n = 300, p = 500, outcome = "continuous",
                        pattern = "contrasting", c1_grid = seq(0, 1, by = 0.1),
                        n_rep = 1000, cores = 4)
```

## Part II. Empirical data analysis (manuscript Section 6)

The manuscript asks whether health-care expenditure mediates the
relationship between economic growth and population health across the
world’s countries. POEM ships that exact data set,
`WHO_health_mediation`: a country-by-year panel (2000–2021) with

- **Exposure** `X`: annual growth in GDP per capita (World Bank),
- **Mediators** `M`: 57 health-expenditure indicators (WHO Global Health
  Expenditure Database),
- **Outcomes** `Y`: five health outcomes, infant mortality (`imr`),
  under-five mortality (`u5mr`), life expectancy (`leb`), low
  birthweight (`lbw`), and undernourishment (`pou`),

with WHO region and World Bank income group available as the grouping
and confounding variables.

``` r

data(WHO_health_mediation)
dim(WHO_health_mediation)
#> [1] 2002   68
table(WHO_health_mediation$region)
#> 
#>  AFR  AMR  EMR  EUR SEAR  WPR 
#>  704  616  176  198  110  198
```

### Section 6, Table for infant mortality (IMR)

[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_analysis.md)
reproduces the manuscript’s per-outcome table: it fits the global model
(region, income, and year as confounders), then a separate model within
each WHO region (income and year as confounders) and each income group
(region and year as confounders), and reports the HDMM and PE p-values
for the total indirect effect together with the active mediators the PE
component identifies.

The function uses two tuning grids that match the manuscript: a narrower
`lambda_grid_global = seq(0.1, 5, length.out = 100)` for the global
all-country model and a wider
`lambda_grid = seq(0.1, 10, length.out = 100)` for every subgroup model.
Both grids are the defaults, so no explicit arguments are needed to
reproduce the manuscript exactly.

``` r

imr <- WHO_mediation_analysis("imr", groupings = c("global", "region", "income"))
imr[, c("group", "n_countries", "pval_hdmm", "pval_pe", "active_mediators")]
#>  group        n_countries pval_hdmm pval_pe    active_mediators 
#>  ALL          91          0.04094   7.633e-29  gge_gdp          
#>  AFR          32          0.4672    0.4672     none             
#>  AMR          28          0.4129    0.4129     none             
#>  EMR          8           0.03045   0.03045    none             
#>  EUR          9           0.003479  1.579e-44  gge_gdp          
#>  SEAR         5           0.3087    4.128e-41  ext_usd2021_pc   
#>  WPR          9           0.8721    3.766e-83  chi_che, pvtd_gdp
#>  Low          16          0.3157    2.101e-62  pvtd_usd2021     
#>  Lower-middle 27          0.9546    0.9546     none             
#>  Upper-middle 26          0.4945    0.4945     none             
#>  High         22          0.02378   3.609e-103 oops_che, shi_che
```

Read this against the manuscript’s IMR table. The headline
reproductions:

- **Global (ALL):** general government health expenditure as a share of
  GDP (`gge_gdp`) is the active mediator; HDMM is borderline
  (`pval_hdmm` $`\approx 0.041`$) while PE is essentially zero.
- **SEAR (South-East Asia):** HDMM sees nothing (`pval_hdmm` $`\approx
  0.31`$) yet PE rejects decisively, naming external health expenditure
  per capita (`ext_usd2021_pc`). This is a textbook
  heterogeneous-mediation rescue.
- **WPR (Western Pacific):** HDMM is near 1 (`pval_hdmm`
  $`\approx 0.87`$) while PE detects compulsory health insurance and
  private expenditure (`chi_che`, `pvtd_gdp`).
- **EUR (Europe):** both tests reject; PE identifies general government
  expenditure as a share of GDP (`gge_gdp`).
- **AFR, AMR:** no active mediation, by either test.

These are precisely the cases the manuscript highlights, where
contrasting or heterogeneous mediation hides the signal from the
conventional test and the PE component recovers it.

The income-group rows are returned in the same call; the strongest
signals appear in the low- and high-income groups, matching the
manuscript.

[`summary()`](https://rdrr.io/r/base/summary.html) digests the whole
table, highlighting the groups where the power-enhanced test detects
mediation the benchmark misses (the rows the manuscript prints in bold):

``` r

summary(imr)
#> POEM comparison: IMR
#>   groups: 11 (11 with data), alpha = 0.05
#>   detected by benchmark (HDMM): 4   by power-enhanced (PE): 7
#>   PE detects mediation in 3 groups the benchmark misses:
#>  group n_countries pval_hdmm  pval_pe  active_mediators
#>   SEAR           5     0.309 4.13e-41    ext_usd2021_pc
#>    WPR           9     0.872 3.77e-83 chi_che, pvtd_gdp
#>    Low          16     0.316 2.10e-62      pvtd_usd2021
#>   most-flagged mediators:
#>     gge_gdp            2
#>     chi_che            1
#>     ext_usd2021_pc     1
#>     oops_che           1
#>     pvtd_gdp           1
#>     pvtd_usd2021       1
#>     shi_che            1
```

The supplement also gives an **extended** version of each table that
breaks the power-enhanced result out by multiplicity method. Pass
`full_table = TRUE` to get that layout, with the PE p-value and active
mediators reported separately under Bonferroni (FWER), BH, and BY (FDR):

``` r

WHO_mediation_analysis("imr", groupings = "global", full_table = TRUE)
#>  grouping group n_countries pval_hdmm pval_pe_bonferroni active_bonferroni
#>  global   ALL   91          0.04094   7.633e-29          gge_gdp          
#>  pval_pe_bh active_bh        pval_pe_by active_by
#>  3.136e-57  chi_che, gge_gdp 7.633e-29  gge_gdp  
#> 
#> Outcome model: imr
```

### The other outcomes, and the region-by-income cells

The remaining manuscript outcomes (under-five mortality, life
expectancy, low birthweight, and undernourishment) run identically; only
the outcome code changes:

``` r

for (y in c("u5mr", "leb", "lbw", "pou"))
  print(WHO_mediation_analysis(y, groupings = c("global", "region", "income")))
```

The manuscript also drills into the **region-by-income** cells (for
example, AFR-Low or WPR-High). Build any one cell with
[`WHO_mediation_design()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_design.md),
which assembles the exposure, outcome, mediator, and confounder blocks
for a chosen subset, then call
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md).
A small cell can contain indicators that never vary inside it; passing
`drop_constant = TRUE` drops those zero-variance columns with a warning
and keeps the reported active set numbered against the original
mediators
([`WHO_mediation_analysis()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_analysis.md)
does the same internally):

``` r

cell <- WHO_mediation_design("imr", region = "AFR", income = "Upper-middle")
c(n_countries= cell$n_countries, mediators = ncol(cell$M))
#> n_countries   mediators 
#>           4          57

keep <- apply(cell$M, 2L, sd) > 0
fit <- tryCatch(
  pe_mediation(scale(cell$X), cell$Y - mean(cell$Y),
               scale(cell$M[, keep, drop = FALSE]),
               Z = cell$Z, outcome = "continuous", scale = FALSE,
               lambda_grid = seq(0.1, 5, length.out = 100)),
  error = function(e) NULL)
if (is.null(fit)) "singular covariance (cell too small)" else
  cell$indicators[keep][attr(fit, "active_mediators")]
#> [1] "pvtd_ncu2021_pc"
```

### Why a mediator? Opening up the screening step

For any fit,
[`pe_mediators()`](https://yelleknek.github.io/POEM/reference/pe_mediators.md)
shows the evidence the PE component acted on, the two path statistics,
their p-values, the combined screening p-value, and the active flag for
every candidate the penalized fit selected. This is the per-mediator
detail behind the last column of the manuscript’s tables.

``` r

des <- WHO_mediation_design("imr")                  # global IMR design
gfit <- pe_mediation(scale(des$X), des$Y - mean(des$Y), scale(des$M),
                     Z = des$Z, outcome = "continuous", scale = FALSE,
                     lambda_grid = seq(0.1, 5, length.out = 100))
pe_mediators(gfit)
#>  mediator t_outcome t_exposure screen_p  selected
#>  2        8.234     -1.686     0.09183   FALSE   
#>  4        2.874     0.6476     0.5172    FALSE   
#>  9        1.188     0.537      0.5913    FALSE   
#>  10       2.003     0.9471     0.3436    FALSE   
#>  14       -0.7627   1.231      0.4457    FALSE   
#>  16       2.302     -0.8475    0.3967    FALSE   
#>  17       -8.491    -1.315     0.1886    FALSE   
#>  19       -2.944    0.6753     0.4995    FALSE   
#>  21       -5.535    3.111      0.001862  FALSE   
#>  25       -2.673    -2.494     0.01264   FALSE   
#>  26       -4.26     -3.732     0.0001901  TRUE   
#>  31       2.809     3.451      0.004962  FALSE   
#>  35       -4.063    -1.373     0.1697    FALSE   
#>  38       -0.2667   -0.2293    0.8187    FALSE   
#>  39       -1.438    2.224      0.1505    FALSE   
#>  41       -3.921    0.325      0.7452    FALSE   
#>  45       -5.506    1.731      0.08343   FALSE   
#>  47       -4.697    0.3339     0.7384    FALSE   
#>  50       -0.7842   0.18       0.8571    FALSE   
#>  52       0.2591    0.1229     0.9022    FALSE   
#>  53       6.643     -1.532     0.1255    FALSE   
#>  57       7.506     1.361      0.1735    FALSE
```

## Reproducibility notes

- The global, regional, and income analyses above are deterministic.
  [`WHO_mediation_analysis()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_analysis.md)
  uses a tuning grid whose defaults are set to be
  `lambda_grid_global = seq(0.1, 5, length.out = 100)` for the
  all-country model. Calling the function with no explicit grid
  arguments reproduces the active mediators and p-values shown here,
  which the package’s test suite checks as a fixed benchmark.
- The Monte Carlo sections are stochastic. The reduced runs reproduce
  the *shape* of the manuscript’s curves; the `eval = FALSE` blocks give
  the exact `n = 300`, `p = 500`, 1,000-replication settings. Set
  `cores` above 1 to parallelize and `seed` for exact reproducibility.
- Exact real-data p-values depend on the preprocessing the manuscript
  fixes (standardizing the exposure and mediators, centering the
  outcome, leaving the confounder indicators unscaled, and using
  treatment contrasts with alphabetical reference levels for the region
  and income factors);
  [`WHO_mediation_analysis()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_analysis.md)
  applies all of this for you. See the guided-tour vignette for doing it
  by hand.

## Reference

Yu, X., and Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.
