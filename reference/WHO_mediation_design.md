# Build a mediation design from the WHO health-expenditure data

Assembles the exposure, outcome, mediator, and confounder matrices for
one of the analyses in the manuscript from
[WHO_health_mediation](https://yelleknek.github.io/POEM/reference/WHO_health_mediation.md),
ready to pass to
[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md).
The exposure is GDP per capita growth, the mediators are the 57
health-expenditure indicators, and the outcome is the chosen health
indicator. Optionally restrict to a WHO region and/or a World Bank
income group, reproducing the manuscript's region-specific,
income-specific, and region-by-income models.

## Usage

``` r
WHO_mediation_design(
  outcome = c("imr", "u5mr", "leb", "lbw", "pou"),
  region = NULL,
  income = NULL,
  data = WHO_health_mediation
)
```

## Arguments

- outcome:

  Which health outcome to use as \\Y\\: `"imr"` (infant mortality rate),
  `"u5mr"` (under-five mortality rate), `"leb"` (life expectancy at
  birth), `"lbw"` (prevalence of low birthweight), or `"pou"`
  (prevalence of undernourishment).

- region:

  Optional character vector of WHO regions to keep (`"AFR"`, `"AMR"`,
  `"EMR"`, `"EUR"`, `"SEAR"`, `"WPR"`). Default `NULL` keeps all
  regions.

- income:

  Optional character vector of income groups to keep (`"Low"`,
  `"Lower-middle"`, `"Upper-middle"`, `"High"`). Default `NULL` keeps
  all groups.

- data:

  The source data frame. Defaults to
  [WHO_health_mediation](https://yelleknek.github.io/POEM/reference/WHO_health_mediation.md).

## Value

A list with `X` (exposure, \\n \times 1\\), `Y` (outcome vector), `M`
(the \\n \times 57\\ mediator matrix), `Z` (the confounder matrix, or
`NULL` if no covariate varies), and the metadata `outcome`, `n`
(country-year observations), `n_countries` (unique countries),
`indicators`, `region`, and `income`.

## Details

Rows with a missing outcome or missing exposure are dropped (the article
handles each outcome's coverage separately, which is why the
low-birthweight and undernourishment outcomes have fewer observations).
The confounder matrix \\Z\\ contains the WHO region and income group (as
indicator variables) and the year (numeric), exactly as in the article,
but only for covariates that still vary after subsetting: a
region-specific model drops the (now constant) region, and unused factor
levels are dropped so no all-zero indicator columns remain.

To reproduce the manuscript's exact test statistics, fit with the
article's preprocessing, standardizing the exposure and mediators and
centering the outcome while leaving the confounder indicators unscaled:


      des <- WHO_mediation_design("imr")
      pe_mediation(scale(des$X), des$Y - mean(des$Y), scale(des$M),
                   Z = des$Z, outcome = "continuous", scale = FALSE,
                   lambda_grid = seq(0.1, 5, length.out = 100))

The simpler call
`pe_mediation(des$X, des$Y, des$M, des$Z, outcome = "continuous")`
(which standardizes every block, including \\Z\\) identifies the same
active mediators; only the benchmark Wald p-value shifts, because
scaling the indicator columns changes that test.

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEM/reference/pe_mediation.md),
[WHO_health_mediation](https://yelleknek.github.io/POEM/reference/WHO_health_mediation.md),
[WHO_indicator_codebook](https://yelleknek.github.io/POEM/reference/WHO_indicator_codebook.md).

Other WHO mediation data:
[`WHO_health_mediation`](https://yelleknek.github.io/POEM/reference/WHO_health_mediation.md),
[`WHO_indicator_codebook`](https://yelleknek.github.io/POEM/reference/WHO_indicator_codebook.md),
[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEM/reference/WHO_mediation_analysis.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
# Global infant-mortality model: GDP growth -> health spending -> IMR.
des <- WHO_mediation_design("imr")
c(n = des$n, mediators = ncol(des$M))
#>         n mediators 
#>      2002        57 
fit <- pe_mediation(scale(des$X), des$Y - mean(des$Y), scale(des$M),
  Z = des$Z, outcome = "continuous", scale = FALSE,
  lambda_grid = seq(0.1, 5, length.out = 100))
fit
#>  term                  value   
#>  stat_hdmm             4.178   
#>  pval_hdmm             0.0409  
#>  stat_pe               124.2   
#>  j_pe                  120     
#>  pval_pe               < 0.0001
#>  total_indirect_effect 0.4034  
#>  total_indirect_lower  0.01661 
#>  total_indirect_upper  0.7901  
#>  n_active_mediators    1       
#>  df                    1       
#>  n_candidate_mediators 57      
#>  n_observations        2002    
#> 
#> Outcome model: continuous (linear)
#> Active mediators identified (1): 26
# The active mediator's full name:
WHO_indicator_codebook$description[
  match(des$indicators[attr(fit, "active_mediators")],
        WHO_indicator_codebook$indicator)]
#> [1] "General Government Expenditure (GGE) as % of GDP"

# A region-specific model (Europe).
des_eur <- WHO_mediation_design("imr", region = "EUR")
```
