# Pre-release quality control for POEM: the static invariants.
#
# Runs every invariant the project conventions state (the package's
# working conventions, and the portfolio QC master where the package does not
# override it), mechanically, and prints one PASS/FAIL line per check with
# the offending files when a check fails. Exit status is nonzero on any
# failure, so tools/release_gate.R can gate on it. Modeled on DMAR's
# tools/release_checks.R; the POEM-specific checks are the citation status
# of the in-press article, the anonymization residue from the double-blind
# review, the seed save-and-restore idiom, and the hex logo.
#
# Usage, from the package root:
#   Rscript tools/release_checks.R
#
# This file is not shipped: tools/ is excluded by .Rbuildignore. It never
# carries the literal strings it polices; those patterns are assembled from
# character codes.

root <- normalizePath(".")
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  stop("Run from the package root.", call. = FALSE)
}

failures <- 0L
check <- function(name, ok, detail = character()) {
  status <- if (ok) "PASS" else "FAIL"
  cat(sprintf("[%s] %s\n", status, name))
  if (!ok) {
    failures <<- failures + 1L
    for (d in detail) cat("       ", d, "\n", sep = "")
  }
  invisible(ok)
}

r_files    <- list.files("R", pattern = "[.]R$", full.names = TRUE)
rd_files   <- list.files("man", pattern = "[.]Rd$", full.names = TRUE)
vig_rmd    <- list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE)
prose      <- c(r_files, vig_rmd, "DESCRIPTION", "NEWS.md", "README.md", "inst/CITATION")

grep_files <- function(pattern, files, ...) {
  hits <- character()
  for (f in files) {
    lines <- readLines(f, warn = FALSE)
    m <- grep(pattern, lines, ...)
    if (length(m)) hits <- c(hits, sprintf("%s:%d", f, m))
  }
  hits
}
code_lines <- function(f) {
  lines <- readLines(f, warn = FALSE)
  lines[!grepl("^\\s*#", lines)]
}
grep_code <- function(pattern, files = r_files) {
  hits <- character()
  for (f in files) {
    m <- grep(pattern, code_lines(f))
    if (length(m)) hits <- c(hits, sprintf("%s (%d line%s)", f, length(m), if (length(m) == 1) "" else "s"))
  }
  hits
}

## ---- Patterns, and their positive controls -----------------------------
# A detector that cannot fire proves nothing, so each one is first run on a
# known-bad string it must flag (the portfolio's rule 8 for gates).
dash_pat    <- paste0("[", intToUtf8(c(0x2014, 0x2013)), "]")
hyphen_pat  <- "(^|[^-!<0-9])--($|[^->0-9a-zA-Z])"
residue_pat <- paste(c(intToUtf8(c(65,110,111,110,121,109,111,117,115)),        # the review-era author
                       intToUtf8(c(97,110,111,110,121,109,105,122,101,100)),    # "... for review"
                       intToUtf8(c(100,111,117,98,108,101,45,98,108,105,110,100)),
                       "\\(minor revision\\)", "\\(submitted\\)", "\\(under review\\)",
                       "\\(in preparation\\)", "conditionally accepted",
                       "JASA-T&M"), collapse = "|")
is_code <- function(txt) {
  parsed <- tryCatch(parse(text = txt, keep.source = FALSE), error = function(e) NULL)
  if (is.null(parsed) || !length(parsed)) return(FALSE)
  any(vapply(parsed, function(e) is.call(e) && !identical(e[[1L]], as.name("~")), logical(1)))
}
controls <- c(
  dash    = grepl(dash_pat, paste0("a ", intToUtf8(0x2014), " b")),
  hyphen  = grepl(hyphen_pat, "the benchmark -- the one") && !grepl(hyphen_pat, "2000--2021"),
  residue = grepl(residue_pat, "Yu and Kelley (minor revision)", ignore.case = TRUE),
  code    = is_code("fit <- pe_mediation(X, Y, M)") && !is_code("the mediators all push the same way"))
check("detector positive controls (each detector flags its known-bad case)",
      all(controls), names(controls)[!controls])

## ---- Examples ---------------------------------------------------------------

h <- grep_files("\\\\dontrun", c(r_files, rd_files))
check("no \\dontrun anywhere", length(h) == 0, h)

# No \donttest{} anywhere (the portfolio master, decided 2026-08-10 and
# adopted for POEM on 2026-09-25): one block makes --as-cran run the whole
# example corpus twice, so every example is kept fast enough to run.
h <- grep_files("\\\\donttest", c(r_files, rd_files))
check("no \\donttest anywhere", length(h) == 0, h)

# No commented-out code in any @examples block (CRAN: "Some code lines in
# examples are commented out. Please never do that."). Each contiguous
# window of comment lines is parsed with the markers stripped; prose does
# not parse and is never flagged.
commented_hits <- character()
for (f in r_files) {
  lines <- readLines(f, warn = FALSE)
  for (s in grep("^#' @examples", lines)) {
    i <- s + 1L; ex <- character(); idx <- integer()
    while (i <= length(lines) && grepl("^#'", lines[i]) && !grepl("^#' @[a-zA-Z]", lines[i])) {
      ex <- c(ex, sub("^#' ?", "", lines[i])); idx <- c(idx, i); i <- i + 1L
    }
    cm <- grep("^\\s*#", ex)
    if (!length(cm)) next
    for (run in split(cm, cumsum(c(1L, diff(cm) != 1L)))) {
      txt <- sub("^\\s*#+\\s?", "", ex[run]); n <- length(run); flagged <- integer()
      for (a in seq_len(n)) for (b in a:n) {
        w <- paste(txt[a:b], collapse = "\n")
        if (nzchar(trimws(w)) && is_code(w)) flagged <- c(flagged, run[a:b])
      }
      if (length(flagged)) commented_hits <- c(commented_hits, sprintf("%s:%d", f, idx[sort(unique(flagged))]))
    }
  }
}
check("no commented-out code in any @examples block", length(commented_hits) == 0, commented_hits)

ex_lines <- function(f) {
  lines <- readLines(f, warn = FALSE); out <- integer(); in_ex <- FALSE
  for (i in seq_along(lines)) {
    if (grepl("^#' @examples", lines[i])) { in_ex <- TRUE; next }
    if (in_ex && (!grepl("^#'", lines[i]) || grepl("^#' @[a-zA-Z]", lines[i]))) in_ex <- FALSE
    if (in_ex) out <- c(out, i)
  }
  out
}
guard_hits <- character(); seed_hits <- character()
for (f in r_files) {
  lines <- readLines(f, warn = FALSE)
  for (i in ex_lines(f)) {
    if (grepl("requireNamespace", lines[i]) && grepl("^#' [^#]", lines[i]))
      guard_hits <- c(guard_hits, sprintf("%s:%d", f, i))
    if (grepl("set\\.seed\\(", lines[i]) && !grepl("set\\.seed\\(113\\)", lines[i]))
      seed_hits <- c(seed_hits, sprintf("%s:%d", f, i))
  }
}
check("no requireNamespace guard executing in @examples", length(guard_hits) == 0, guard_hits)
check("every example seed is set.seed(113)", length(seed_hits) == 0, seed_hits)

## ---- Run-time policies (code lines only) -------------------------------------

# Any function using randomness takes seed = NULL and goes through
# .poem_local_seed() (withr::local_seed underneath), which sets the seed for
# the call and restores the caller's state on exit. Nothing in R/ calls
# set.seed() or reads or writes .Random.seed or the global environment: that
# idiom is what CRAN's review returned DMAR for on 2026-09-04.
h <- grep_code("set\\.seed\\(")
check("no set.seed() anywhere in R/ (seeds go through .poem_local_seed())", length(h) == 0, h)
h <- grep_code("\\.Random\\.seed|\\.GlobalEnv|globalenv\\(\\)")
check("no function reads or writes .Random.seed or the global environment", length(h) == 0, h)
seed_default <- grep_files("^[^#]*function\\(.*seed *= *[0-9]", r_files)
check("no baked-in default seed in any signature", length(seed_default) == 0, seed_default)
h <- grep_code("<<-")
check("no <<- anywhere in R/", length(h) == 0, h)
h <- grep_code("options\\(\\s*warn")
check("no options(warn = ...) anywhere in R/", length(h) == 0, h)
h <- grep_code("(filename|file|path)\\s*=\\s*[\"']")
check("no default file path in any function signature or call", length(h) == 0, h)

# Full precision: round(), signif(), and formatC() only in the display layer
# or in an integer-validation comparison.
display <- c("R/poem_tbl.R", "R/summary_poem_tbl.R")
rnd <- character()
for (f in setdiff(r_files, display)) {
  cl <- code_lines(f)
  m <- grep("\\b(round|signif|formatC)\\(", cl)
  m <- m[!grepl("(==|!=|-)\\s*round\\(", cl[m])]
  if (length(m)) rnd <- c(rnd, sprintf("%s: %s", f, trimws(cl[m])))
}
check("no round()/signif()/formatC() outside the display layer and validation",
      length(rnd) == 0, rnd)

## ---- DESCRIPTION, NEWS, citation -----------------------------------------------

desc <- read.dcf("DESCRIPTION")
ver <- unname(desc[, "Version"])
news_head <- grep("^# ", readLines("NEWS.md", warn = FALSE), value = TRUE)[1]
check(sprintf("NEWS.md leads with version %s", ver),
      identical(sub("^# [A-Za-z]+ ", "", news_head), ver), news_head)

auth <- eval(parse(text = desc[, "Authors@R"]))
roles <- lapply(auth, function(p) p$role)
fam <- vapply(auth, function(p) p$family, "")
auth_ok <- identical(fam, c("Yu", "Kelley")) &&
  identical(roles[[1]], "aut") && setequal(roles[[2]], c("aut", "cre")) &&
  identical(unname(auth[[2]]$email), "kkelley@nd.edu")
check("Authors@R: Yu (aut), Kelley (aut, cre, kkelley@nd.edu)", auth_ok,
      format(auth, include = c("given", "family", "role", "email")))
url_ok <- "URL" %in% colnames(desc) && grepl("github.com/yelleKneK/POEM", desc[, "URL"]) &&
  "BugReports" %in% colnames(desc) && grepl("github.com/yelleKneK/POEM/issues", desc[, "BugReports"])
check("URL and BugReports name github.com/yelleKneK/POEM", url_ok)

# The article's citation status is one fact stated in several places, and
# the places must agree: while inst/CITATION carries note = "In press", the
# DESCRIPTION and README cite it as in press; once it carries a volume, no
# "in press" may remain anywhere.
cit_txt <- paste(readLines("inst/CITATION", warn = FALSE), collapse = "\n")
published <- grepl("\\bvolume\\s*=", cit_txt)
in_press_hits <- grep_files("in press", c(prose, man = rd_files), ignore.case = TRUE)
if (published) {
  check("article published: no 'in press' left anywhere", length(in_press_hits) == 0, in_press_hits)
} else {
  need <- c(inst_CITATION = grepl('note\\s*=\\s*"In press"', cit_txt),
            DESCRIPTION = grepl("(in press)", desc[, "Description"], fixed = TRUE),
            README = any(grepl("(in press)", readLines("README.md", warn = FALSE), fixed = TRUE)),
            CITATION.cff = !file.exists("CITATION.cff") ||
              any(grepl("status: in-press", readLines("CITATION.cff", warn = FALSE), fixed = TRUE)))
  check("article in press: CITATION, DESCRIPTION, README, and CITATION.cff agree",
        all(need), names(need)[!need])
}
h <- grep_files(residue_pat, c(prose, rd_files), ignore.case = TRUE)
check("no review-era residue (anonymized authors, pre-acceptance citation status)",
      length(h) == 0, h)

cit <- tryCatch(utils::readCitationFile("inst/CITATION", meta = as.list(desc[1, ])),
                error = function(e) conditionMessage(e))
check("inst/CITATION readable from DESCRIPTION metadata alone (CRAN incoming state)",
      inherits(cit, "bibentry") && !grepl("packageVersion", cit_txt),
      if (!inherits(cit, "bibentry")) cit)

# The working-conventions file (untracked; its name is assembled from
# character codes so this tracked script never carries it) states the
# version in its status block and must agree.
conv_file <- intToUtf8(c(67, 76, 65, 85, 68, 69, 46, 109, 100))
lock <- c(CITATION.cff = !file.exists("CITATION.cff") ||
            any(grepl(sprintf("^version: %s$", ver), readLines("CITATION.cff", warn = FALSE))),
          conventions = !file.exists(conv_file) ||
            any(grepl(sprintf("Version: %s", ver), readLines(conv_file, warn = FALSE), fixed = TRUE)))
check(sprintf("version %s in lockstep (DESCRIPTION, NEWS, CITATION.cff, the conventions status line)", ver),
      all(lock), names(lock)[!lock])

## ---- House style ---------------------------------------------------------------

h <- grep_files(dash_pat, prose)
check("no em or en dashes", length(h) == 0, h)
h <- grep_files(hyphen_pat, setdiff(prose, "DESCRIPTION"))
h <- h[!grepl("R CMD|--as-cran|--run-donttest", vapply(h, function(x) {
  p <- strsplit(x, ":")[[1]]; readLines(p[1], warn = FALSE)[as.integer(p[2])] }, ""))]
check("no double hyphen used as a dash in prose (numeric ranges excepted)", length(h) == 0, h)
nonascii <- character()
for (f in c(r_files, "DESCRIPTION")) {
  m <- grep("[^\\x01-\\x7F]", readLines(f, warn = FALSE, encoding = "UTF-8"), perl = TRUE)
  if (length(m)) nonascii <- c(nonascii, sprintf("%s:%d", f, m))
}
check("no non-ASCII characters in R/ or DESCRIPTION", length(nonascii) == 0, nonascii)

## ---- Structure, data, vignettes, logo ----------------------------------------

data_objs <- sub("[.]rda$", "", list.files("data", pattern = "[.]rda$"))
rd_db <- tools::Rd_db(dir = ".")
aliases <- unlist(lapply(rd_db, function(rd) {
  unlist(lapply(rd[vapply(rd, function(x) identical(attr(x, "Rd_tag"), "\\alias"), logical(1))],
                function(x) paste(unlist(x), collapse = "")))
}))
has_fs <- vapply(rd_db, function(rd) {
  tags <- vapply(rd, function(x) attr(x, "Rd_tag"), "")
  all(c("\\format", "\\source") %in% tags)
}, logical(1))
fs_aliases <- unlist(lapply(names(rd_db)[has_fs], function(n) {
  rd <- rd_db[[n]]
  unlist(lapply(rd[vapply(rd, function(x) identical(attr(x, "Rd_tag"), "\\alias"), logical(1))],
                function(x) paste(unlist(x), collapse = "")))
}))
# example_* objects share one page (example_data.Rd).
undoc <- data_objs[!data_objs %in% fs_aliases]
check("every shipped data set documented with \\format and \\source", length(undoc) == 0, undoc)

vig_names <- sub("[.]Rmd$", "", basename(vig_rmd))
dangling <- character()
for (hit in grep_files('vignette\\("([^"]+)"', c(r_files, vig_rmd, "README.md"))) {
  parts <- strsplit(hit, ":")[[1]]
  line <- readLines(parts[1], warn = FALSE)[as.integer(parts[2])]
  for (call in regmatches(line, gregexpr('vignette\\("([^"]+)"', line))[[1]]) {
    nm <- sub('vignette\\("([^"]+)".*', "\\1", call)
    if (!nm %in% vig_names) dangling <- c(dangling, paste0(hit, " -> ", nm))
  }
}
check("no vignette() call names a vignette that does not ship", length(dangling) == 0, dangling)

# A git that actually runs. On this machine /usr/local/bin/git is a stale
# Intel binary that precedes /usr/bin in the PATH R's shell sees; it exits
# 126, and a git check fed its empty output would pass having checked
# nothing. So resolve a working git first and fail loudly without one.
git_bin <- local({
  cands <- unique(c(Sys.which("git"), "/usr/bin/git", "/opt/homebrew/bin/git"))
  cands <- cands[nzchar(cands) & file.exists(cands)]
  ok <- vapply(cands, function(g) identical(suppressWarnings(
    tryCatch(system2(g, "--version", stdout = FALSE, stderr = FALSE), error = function(e) 1L)), 0L), logical(1))
  if (any(ok)) cands[ok][1] else NA_character_
})
tracked <- if (is.na(git_bin)) character() else
  suppressWarnings(system2(git_bin, "ls-files", stdout = TRUE))
check("a working git is available (the git checks below are not vacuous)",
      !is.na(git_bin) && length(tracked) > 0 && is.null(attr(tracked, "status")),
      "no git binary on this machine ran; checked Sys.which, /usr/bin, /opt/homebrew/bin")
gen <- grep("^vignettes/.*[.](R|html)$|^doc/|^Meta/|^inst/doc/|[.]Rcheck/|[.]tar[.]gz$|Rplots[.]pdf$",
            tracked, value = TRUE)
check("no generated build products tracked in git", length(gen) == 0, gen)
private <- grep(paste0("(^|/)", conv_file, "$|[.]docx$|GAP-ANALYSIS|HANDOFF|QC_|_QC|RELEASE_QC|evaluation|audit"),
                tracked, value = TRUE, ignore.case = TRUE)
check("no private working documents tracked in git", length(private) == 0, private)

# The hex logo: shipped, referenced, and still the family master's art.
logo <- c(png = file.exists("man/figures/logo.png"), svg = file.exists("man/figures/logo.svg"),
          readme = any(grepl("man/figures/logo.png", readLines("README.md", warn = FALSE), fixed = TRUE)))
check("hex logo shipped (man/figures/logo.png, logo.svg) and shown in README", all(logo),
      names(logo)[!logo])
family_svg <- file.path(dirname(root), "KenKelleyHexFamily", "hexes", "hex_POEM.svg")
if (file.exists(family_svg) && logo[["svg"]]) {
  check("logo.svg identical to the KenKelleyHexFamily master (hexes/hex_POEM.svg)",
        identical(unname(tools::md5sum("man/figures/logo.svg")), unname(tools::md5sum(family_svg))))
}
if (file.exists("dev/make_POEM_hex.R") && logo[["svg"]]) {
  gen_dir <- tempfile("hex"); dir.create(file.path(gen_dir, "dev"), recursive = TRUE)
  file.copy("dev/make_POEM_hex.R", file.path(gen_dir, "dev"))
  old <- setwd(gen_dir)
  st <- system2("Rscript", "dev/make_POEM_hex.R", stdout = FALSE, stderr = FALSE)
  setwd(old)
  regen <- file.path(gen_dir, "man", "figures", "logo.svg")
  check("dev/make_POEM_hex.R regenerates the shipped logo.svg byte for byte",
        st == 0L && file.exists(regen) &&
          identical(unname(tools::md5sum(regen)), unname(tools::md5sum("man/figures/logo.svg"))))
}

# Tarball hygiene: build one (vignettes off, so the size is a floor) and
# confirm nothing ships that should not.
cat("       building the tarball to inspect it...\n")
bld <- tempfile("bld"); dir.create(bld)
old_wd <- setwd(bld)
ok_build <- system2("R", c("CMD", "build", "--no-manual", "--no-build-vignettes", shQuote(root)),
                    stdout = FALSE, stderr = FALSE) == 0
setwd(old_wd)
if (ok_build) {
  tb <- list.files(bld, pattern = "[.]tar[.]gz$", full.names = TRUE)[1]
  contents <- untar(tb, list = TRUE)
  bad <- c(grep("tools/|dev/|[.]github|_pkgdown|CITATION[.]cff|[.]Rcheck|cran-comments|logo-hires|Rplots",
                contents, value = TRUE),
           grep("^[^/]+/(?!NEWS[.]md$|README[.]md$|LICENSE[.]md$)[^/]+[.]md$",
                contents, value = TRUE, perl = TRUE))
  check("tarball contains no working files", length(bad) == 0, bad)
  sz <- file.info(tb)$size / 1048576
  check(sprintf("tarball floor %.2f Mb below the 5 Mb guideline", sz), sz < 5)
} else {
  check("R CMD build succeeds", FALSE, "build failed; run it by hand for the error")
}

cat(sprintf("\n%d failure(s).\n", failures))
if (failures > 0) quit(status = 1L)
