# Assemble the WHO / World Bank health-expenditure mediation data shipped with
# POEM (objects `WHO_health_mediation` and `WHO_indicator_codebook`). The raw
# source files live OUTSIDE the package and are not redistributed; the shipped
# objects (data/WHO_health_mediation.rda) are the reproducible output. This
# script records exactly how they were built. To re-run it, set `src` below to
# the directory holding the raw merged files.
#
# Provenance. The exposure (annual % growth of GDP per capita) and the five
# health outcomes come from World Bank Open Data; the 57 health-expenditure
# indicators (the candidate mediators) and the WHO region / World Bank income
# group labels come from the WHO Global Health Expenditure Database (GHED). The
# two sources were merged and entities / indicators with excessive missingness
# dropped, leaving a country-by-year panel of 91 WHO members, 2000--2021.

src <- Sys.getenv("POEM_WHO_RAW", unset = "path/to/raw-source-files")

# Melt a wide country-by-year table (Country.Code, X2000, ..., X20YY) into a
# long (code, year, value) frame.
melt_wide <- function(df, valname) {
  yrcols <- grep("^X[0-9]{4}$", names(df), value = TRUE)
  years  <- as.integer(sub("^X", "", yrcols))
  data.frame(
    code = rep(df$Country.Code, each = length(yrcols)),
    year = rep(years, times = nrow(df)),
    value = as.numeric(t(as.matrix(df[, yrcols]))),
    stringsAsFactors = FALSE) |>
    `names<-`(c("code", "year", valname))
}

# The GHED panel (mediators + metadata) and the GDP exposure are the same
# across the five outcome files; take them from the IMR file, which has the
# full 91-country coverage.
imr_env <- new.env(); load(file.path(src, "data_GDPR_GHED_IMR.RData"), imr_env)
ghed <- imr_env$data_ghed
indicators <- names(ghed)[6:ncol(ghed)]          # the 57 health-spending codes
stopifnot(length(indicators) == 57L)

master <- ghed[, c("code", "country", "region", "income", "year", indicators)]
master$code   <- as.character(master$code)
master$country<- as.character(master$country)
master$region <- as.character(master$region)
# Order the income groups from low to high (World Bank ordering).
master$income <- factor(as.character(master$income),
                        levels = c("Low", "Lower-middle", "Upper-middle", "High"))

# Exposure: GDP per capita growth (annual %).
gdp <- melt_wide(imr_env$data_gdpr, "gdp_growth")

# Each outcome comes from its own file (different country / year coverage).
outcome_files <- list(
  imr  = c("data_GDPR_GHED_IMR.RData",  "data_imr"),
  u5mr = c("data_GDPR_GHED_U5MR.RData", "data_u5mr"),
  leb  = c("data_GDPR_GHED_LEB.RData",  "data_leb"),
  lbw  = c("data_GDPR_GHED_LBW.RData",  "data_lbw"),
  pou  = c("data_GDPR_GHED_POU.RData",  "data_pou"))

WHO_health_mediation <- merge(master, gdp, by = c("code", "year"), all.x = TRUE)
for (nm in names(outcome_files)) {
  e <- new.env(); load(file.path(src, outcome_files[[nm]][1]), e)
  long <- melt_wide(get(outcome_files[[nm]][2], e), nm)
  WHO_health_mediation <- merge(WHO_health_mediation, long,
                                by = c("code", "year"), all.x = TRUE)
}

# Tidy column order: keys, covariates, exposure, the five outcomes, then the
# 57 mediators. Sort by country then year.
key_cols <- c("code", "country", "region", "income", "year")
out_cols <- c("gdp_growth", "imr", "u5mr", "leb", "lbw", "pou")
WHO_health_mediation <- WHO_health_mediation[
  order(WHO_health_mediation$code, WHO_health_mediation$year),
  c(key_cols, out_cols, indicators)]
rownames(WHO_health_mediation) <- NULL

# Codebook: the 57 indicator codes and their full names, transcribed from the
# article's supplementary codebook table.
WHO_indicator_codebook <- data.frame(
  indicator = indicators,
  description = c(
    "Current Health Expenditure (CHE) as % of Gross Domestic Product (GDP)",
    "Current Health Expenditure (CHE) per capita in US$",
    "Current Health Expenditure (CHE)",
    "Domestic General Government Health Expenditure (GGHE-D)",
    "Domestic Private Health Expenditure (PVT-D)",
    "External Health Expenditure (EXT)",
    "Domestic Health Expenditure (DOM) as % of Current Health Expenditure (CHE)",
    "Domestic General Government Health Expenditure (GGHE-D) as % of CHE",
    "Domestic Private Health Expenditure (PVT-D) as % of CHE",
    "Out-of-pocket Spending (OOPS) as % of CHE",
    "Voluntary Prepayments as % of CHE",
    "External Health Expenditure (EXT) as % of CHE",
    "Domestic General Government Health Expenditure (GGHE-D) as % of GDP",
    "Domestic General Government Health Expenditure (GGHE-D) as % of General Government Expenditure (GGE)",
    "Domestic General Government Health Expenditure (GGHE-D) per capita in US$",
    "Domestic Private Health Expenditure (PVT-D) per capita in US$",
    "Out-of-Pocket Expenditure (OOPS) per capita in US$",
    "External Health Expenditure (EXT) per capita in US$",
    "Compulsory Financing Arrangements (CFA) as % of CHE",
    "Government Financing Arrangements (GFA) as % of CHE",
    "Compulsory Health Insurance (CHI) as % of CHE",
    "Social Health Insurance (SHI) as % of CHE",
    "Compulsory Private Health Insurance (CHI-PVT) as % of CHE",
    "Voluntary Financing Arrangements (VFA) as % of CHE",
    "Voluntary Health Insurance (VHI) as % of CHE",
    "General Government Expenditure (GGE) as % of GDP",
    "Gross Domestic Product (GDP) per capita in US$",
    "Current Health Expenditure (CHE), in million current US$",
    "Domestic General Government Health Expenditure (GGHE-D), in million current US$",
    "Domestic Private Health Expenditure (PVT-D), in million current US$",
    "External Health Expenditure (EXT), in million current US$",
    "Current Health Expenditure (CHE), in current NCU per capita",
    "Domestic General Government Health Expenditure (GGHE-D), in current NCU per capita",
    "Domestic Private Health Expenditure (PVT-D), in current NCU per capita",
    "External Health Expenditure (EXT), in current NCU per capita",
    "Current Health Expenditure (CHE), in current PPP per capita",
    "Domestic General Government Health Expenditure (GGHE-D), in current PPP per capita",
    "Domestic Private Health Expenditure (PVT-D), in current PPP per capita",
    "External Health Expenditure (EXT), in current PPP per capita",
    "Domestic Private Health Expenditure (PVT-D), as % of GDP",
    "External Health Expenditure (EXT), as % of GDP",
    "Current Health Expenditure (CHE), in million constant (2021) NCU",
    "Domestic General Government Health Expenditure (GGHE-D), in million constant (2021) NCU",
    "Domestic Private Health Expenditure (PVT-D), in million constant (2021) NCU",
    "External Health Expenditure (EXT), in million constant (2021) NCU",
    "Current Health Expenditure (CHE), in million constant (2021) US$",
    "Domestic General Government Health Expenditure (GGHE-D), in million constant (2021) US$",
    "Domestic Private Health Expenditure (PVT-D), in million constant (2021) US$",
    "External Health Expenditure (EXT), in million constant (2021) US$",
    "Current Health Expenditure (CHE), in constant (2021) NCU per capita",
    "Domestic General Government Health Expenditure (GGHE-D), in constant (2021) NCU per capita",
    "Domestic Private Health Expenditure (PVT-D), in constant (2021) NCU per capita",
    "External Health Expenditure (EXT), in constant (2021) NCU per capita",
    "Current Health Expenditure (CHE), in constant (2021) US$ per capita",
    "Domestic General Government Health Expenditure (GGHE-D), in constant (2021) US$ per capita",
    "Domestic Private Health Expenditure (PVT-D), in constant (2021) US$ per capita",
    "External Health Expenditure (EXT), in constant (2021) US$ per capita"),
  stringsAsFactors = FALSE)
stopifnot(nrow(WHO_indicator_codebook) == 57L)

# Structural verification against the article before saving.
region_counts <- vapply(
  split(WHO_health_mediation$code, WHO_health_mediation$region),
  function(x) length(unique(x)), integer(1))
stopifnot(
  nrow(WHO_health_mediation) == 2002L,                 # 91 members x 22 years
  length(unique(WHO_health_mediation$code)) == 91L,
  all(range(WHO_health_mediation$year) == c(2000, 2021)),
  !anyNA(WHO_health_mediation[, indicators]),           # mediators complete
  # WHO region country counts match the paper's tallies.
  region_counts[["AFR"]] == 32L, region_counts[["AMR"]] == 28L,
  region_counts[["EMR"]] == 8L,  region_counts[["EUR"]] == 9L,
  region_counts[["SEAR"]] == 5L, region_counts[["WPR"]] == 9L)

if (requireNamespace("usethis", quietly = TRUE)) {
  usethis::use_data(WHO_health_mediation, overwrite = TRUE, compress = "xz")
  usethis::use_data(WHO_indicator_codebook, overwrite = TRUE, compress = "xz")
} else {
  save(WHO_health_mediation, file = "data/WHO_health_mediation.rda", compress = "xz")
  save(WHO_indicator_codebook, file = "data/WHO_indicator_codebook.rda", compress = "xz")
}
