# Documentation and run-time hygiene that R CMD check only catches with the
# manuals built, or not at all. Ported from the DMAR package after two CRAN
# returns in 2026: Rd markup inside \eqn{} (the PDF manual stops on it),
# commented-out code in examples, and functions that write to the global
# environment. Each detector proves itself on a known-bad snippet first, so a
# detector that cannot fire never passes vacuously.

rd_sources <- function() {
  man <- file.path("..", "..", "man")
  if (dir.exists(man)) tools::Rd_db(dir = file.path("..", "..")) else tools::Rd_db("POEM")
}

math_text <- function(rd) {
  out <- character()
  walk <- function(x) {
    tag <- attr(x, "Rd_tag")
    if (!is.null(tag) && tag %in% c("\\eqn", "\\deqn"))
      out <<- c(out, paste(unlist(x), collapse = ""))
    if (is.list(x)) for (el in x) walk(el)
  }
  walk(rd)
  out
}

rd_tags <- function(rd) {
  tags <- character()
  walk <- function(x) {
    tag <- attr(x, "Rd_tag")
    if (!is.null(tag)) tags <<- c(tags, tag)
    if (is.list(x)) for (el in x) walk(el)
  }
  walk(rd)
  unique(tags)
}

examples_text <- function(rd) {
  ex <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\examples"), rd)
  if (!length(ex)) return("")
  paste(unlist(ex), collapse = "")
}

markup_in_math <- "\\\\(code|emph|link|strong|pkg|bold|verb|var|samp|file|env|dQuote|sQuote|href|url)\\{"

test_that("the math-markup detector flags the construct CRAN rejected", {
  bad <- tools::parse_Rd(textConnection(paste0(
    "\\name{x}\\alias{x}\\title{x}\\description{The region is ",
    "\\eqn{(-\\code{delta_lower}, \\code{delta_upper})}.}")))
  expect_true(any(grepl(markup_in_math, math_text(bad))))
  good <- tools::parse_Rd(textConnection(paste0(
    "\\name{x}\\alias{x}\\title{x}\\description{The region is ",
    "\\eqn{(-\\delta_L, +\\delta_U)}, with \\code{delta_lower} as ",
    "\\eqn{\\delta_L}.}")))
  expect_false(any(grepl(markup_in_math, math_text(good))))
})

test_that("no help page carries Rd markup inside \\eqn{} or \\deqn{}", {
  db <- rd_sources()
  skip_if(length(db) == 0L, "no Rd sources available")
  offenders <- names(Filter(function(rd) any(grepl(markup_in_math, math_text(rd))), db))
  expect_identical(offenders, character(0))
})

test_that("no help page uses \\donttest{} or \\dontrun{}", {
  db <- rd_sources()
  skip_if(length(db) == 0L, "no Rd sources available")
  offenders <- names(Filter(function(rd) any(c("\\donttest", "\\dontrun") %in% rd_tags(rd)), db))
  expect_identical(offenders, character(0))
})

test_that("no example block gates on requireNamespace()", {
  db <- rd_sources()
  skip_if(length(db) == 0L, "no Rd sources available")
  offenders <- names(Filter(function(rd) grepl("requireNamespace", examples_text(rd), fixed = TRUE), db))
  expect_identical(offenders, character(0))
})

# Comment lines inside \examples whose text is R code ("Some code lines in
# examples are commented out. Please never do that.", CRAN, 2026-09-04).
commented_code <- function(rd) {
  ex <- examples_text(rd)
  if (!nzchar(ex)) return(character())
  lines <- strsplit(ex, "\n", fixed = TRUE)[[1]]
  cm <- grep("^\\s*#", lines)
  if (!length(cm)) return(character())
  is_code <- function(txt) {
    parsed <- tryCatch(parse(text = txt, keep.source = FALSE), error = function(e) NULL)
    if (is.null(parsed) || !length(parsed)) return(FALSE)
    is_call <- vapply(parsed, is.call, logical(1))
    is_formula <- vapply(parsed, function(e)
      is.call(e) && identical(e[[1L]], as.name("~")), logical(1))
    any(is_call & !is_formula)
  }
  hits <- integer()
  for (run in split(cm, cumsum(c(1L, diff(cm) != 1L)))) {
    txt <- sub("^\\s*#+\\s?", "", lines[run])
    n <- length(run)
    for (a in seq_len(n)) for (b in a:n) {
      window <- paste(txt[a:b], collapse = "\n")
      if (nzchar(trimws(window)) && is_code(window)) hits <- c(hits, run[a:b])
    }
  }
  lines[sort(unique(hits))]
}

test_that("the commented-code detector separates code from prose", {
  rd <- tools::parse_Rd(textConnection(paste0(
    "\\name{x}\\alias{x}\\title{x}\\examples{\n",
    "# The fit is shown for comparison.\n",
    "# pe_mediation(X, Y, M, outcome = \"continuous\")\n",
    "# start <- Sys.time()\n",
    "# Sys.sleep(0.2)\n",
    "pe_lambda_grid(300, 500)\n}")))
  hits <- commented_code(rd)
  expect_length(hits, 3L)
  expect_false(any(grepl("shown for comparison", hits)))
})

test_that("no help page carries commented-out code in its examples", {
  db <- rd_sources()
  skip_if(length(db) == 0L, "no Rd sources available")
  offenders <- Filter(length, lapply(db, commented_code))
  expect_identical(names(offenders), character(0),
                   info = paste(unlist(lapply(names(offenders), function(f)
                     paste(f, offenders[[f]], sep = ": "))), collapse = "\n"))
})

# Run-time policies from the same CRAN review: no writing to the global
# environment (any assign() or rm() aimed at .GlobalEnv, `.Random.seed`, and
# any <<-), no options(warn = -1), and no default path in a function that
# writes a file. The namespace is deparsed and inspected.
namespace_functions <- function() {
  ns <- asNamespace("POEM")
  nms <- ls(ns, all.names = TRUE)
  fns <- Filter(is.function, mget(nms, envir = ns, inherits = FALSE))
  fns[!vapply(fns, is.primitive, logical(1))]
}

test_that("no function touches the global environment, uses <<-, or sets options(warn)", {
  fns <- namespace_functions()
  text <- vapply(fns, function(f) paste(deparse(f, width.cutoff = 500L), collapse = "\n"),
                 character(1))
  offenders <- names(text)[grepl("\\.GlobalEnv|globalenv\\(\\)|\\.Random\\.seed|<<-|options\\(warn", text)]
  expect_identical(offenders, character(0))
})

test_that("no function that writes a file has a default path", {
  fns <- namespace_functions()
  bad <- character()
  for (nm in names(fns)) {
    fm <- formals(fns[[nm]])
    for (arg in intersect(names(fm), c("file", "filename", "path", "dir", "directory"))) {
      if (!is.null(fm[[arg]]) && !identical(fm[[arg]], quote(expr = )))
        bad <- c(bad, sprintf("%s(%s = %s)", nm, arg, deparse(fm[[arg]])))
    }
  }
  expect_identical(bad, character(0))
})
