#!/usr/bin/env Rscript
###############################################################################
##  GS4PB — App → Package Function Sync (Phase 1)
##
##  Extracts updated function bodies from the monolithic app file and patches
##  the matching definitions in GS4PB package R/ files, preserving all
##  roxygen2 headers and @export tags.
##
##  Usage:
##    Rscript inst/sync/sync_app_to_pkg.R \
##      --app  "C:/.../GS4PBAppDockCyV/App/GS_Pipeline_Jan_2026_FnsApp.R" \
##      --pkg  "C:/.../GS4PB" \
##      --rebuild        # optional: run devtools::document() + build() after sync
##      --dry-run        # optional: show diffs without writing files
##
##  Or from within an R session:
##    source(system.file("sync", "sync_app_to_pkg.R", package = "GS4PB"))
##    sync_app_to_pkg(app_file = "...", pkg_root = "...", rebuild = TRUE)
###############################################################################

## ── Function → package file mapping ─────────────────────────────────────────
## Add new functions here as the package grows.

.FUNC_MAP <- list(
  "R/GS_Pipeline_FnsApp.r" = c(
    "VCFtoDF", "getTasObj", "getGenoTas_to_DF", "getGenoQCStats",
    "getFilteredSitesGenoData", "getFilteredTaxaGenoData",
    "getGenoQCStatsFilt1", "getGenoQCStatsFilt2",
    "getMergedData", "getPhenoMEData", "getMergedDataME",
    "getProcessedData", "getPredictionData",
    "getLibStats", "writeCVOutTable", "writeGPOutTable"
  ),
  "R/ST_Functions.R" = c(
    "getemCVR", "getRankedPredictedValues"
  ),
  "R/MT_Functions.R" = c(
    "getMTCVR", "getRankedPredictedValuesMT", "summarize_MT_Fits"
  ),
  "R/ME_Functions.R" = c(
    "getMEPred", "stdizePhVarNames", "fitMEModels_LOFO_Pred",
    "summarize_ME_CV_fits", "stdizeCVCovarNames"
  ),
  "R/Imputation_Functions.R" = c(
    "getImputedData_Num", "getImputedData_LDKNNI",
    "getGenoData_API", "getImpGenoData_API"
  ),
  "R/Opt_TS_Functions.r" = c(
    "getOptimalTS", "getRandomTS", "getTSComparisons", "getTSComparisonsMT"
  ),
  "R/Enviromics_Functions.R" = c(
    "getEnvData", "getEnvKernel", "syncEnvPhenoDat", "plotEnvRel"
  ),
  "R/Genomic_Kinship_Functions.R" = c(
    "getQCFilteredM", "getGmat", "getKinshipDiag",
    "removeDuplicatesFromG", "getTunedG", "getKinshipPCA"
  )
)

## ── Core helpers ─────────────────────────────────────────────────────────────

#' Find the start line of a function definition in a character vector of lines.
#' Matches `funcName <- function(` or `funcName=function(` patterns.
.find_fn_start <- function(lines, fn_name) {
  pat <- paste0("^\\s*", fn_name, "\\s*(<-|=)\\s*function\\s*\\(")
  which(grepl(pat, lines))
}

#' Locate a top-level assignment to `fn_name` and return its 1-based
#' c(start, end) line range, using R's own parser (srcrefs).
#'
#' This replaces a brace-counting scan. Counting `{` and `}` per line
#' mis-detects the closing brace whenever one appears inside a string or a
#' comment, and the failure is silent: the body is cut short and the truncated
#' text then overwrites correct package code. Observed in practice —
#' MT_Functions.R was truncated from 433 lines to 42 (getMTCVR ~184 lines ->
#' ~20), and summarize_MT_Fits was reported "not in app" although it is
#' defined in the monolith. Parsing is exact by construction and additionally
#' reports files that are not valid R instead of mangling them.
.fn_line_range <- function(lines, fn_name, what = "source") {
  exprs <- tryCatch(parse(text = lines, keep.source = TRUE),
                    error = function(e) {
                      warning("  [!] Could not parse ", what, ": ",
                              conditionMessage(e))
                      NULL
                    })
  if (is.null(exprs)) return(NULL)
  srcs <- attr(exprs, "srcref")
  for (i in seq_along(exprs)) {
    e <- exprs[[i]]
    if (is.call(e) && length(e) >= 3L &&
        as.character(e[[1L]]) %in% c("<-", "=", "<<-") &&
        is.symbol(e[[2L]]) &&
        identical(as.character(e[[2L]]), fn_name)) {
      sr <- srcs[[i]]
      return(c(start = as.integer(sr[[1L]]), end = as.integer(sr[[3L]])))
    }
  }
  NULL
}

#' Extract a function's signature + body as a character vector of lines
#' from the monolith text lines.
.extract_fn_from_lines <- function(lines, fn_name) {
  rng <- .fn_line_range(lines, fn_name, what = paste0("app monolith [", fn_name, "]"))
  if (is.null(rng)) return(NULL)
  lines[rng[["start"]]:rng[["end"]]]
}

#' Replace a function's signature + body in pkg_lines with new_body_lines.
#' Lines before the function start (roxygen2 block) are preserved.
.replace_fn_in_lines <- function(pkg_lines, fn_name, new_body_lines) {
  rng <- .fn_line_range(pkg_lines, fn_name, what = paste0("package file [", fn_name, "]"))
  if (is.null(rng)) {
    message("  [not found] ", fn_name, " — not present in package file, skipping")
    return(list(lines = pkg_lines, changed = FALSE))
  }
  start <- rng[["start"]]
  end   <- rng[["end"]]

  old_body <- pkg_lines[start:end]
  if (identical(old_body, new_body_lines)) {
    return(list(lines = pkg_lines, changed = FALSE))
  }

  ## Safety gate: never write a replacement that is not itself valid R.
  if (!isTRUE(tryCatch({ parse(text = new_body_lines); TRUE },
                       error = function(e) FALSE))) {
    warning("  [!] Refusing to replace ", fn_name,
            " — extracted text does not parse on its own")
    return(list(lines = pkg_lines, changed = FALSE))
  }

  updated <- c(
    if (start > 1L) pkg_lines[seq_len(start - 1L)] else character(0),
    new_body_lines,
    if (end < length(pkg_lines)) pkg_lines[seq(end + 1L, length(pkg_lines))] else character(0)
  )
  list(lines = updated, changed = TRUE)
}

## ── Main sync function ────────────────────────────────────────────────────────

#' Sync pipeline functions from app monolith → GS4PB package R/ files.
#'
#' @param app_file  Path to the monolithic app R file.
#' @param pkg_root  Path to the GS4PB package root directory.
#' @param rebuild   If TRUE, run devtools::document() + devtools::build() after sync.
#' @param dry_run   If TRUE, report diffs but do not write any files.
#' @param func_map  Named list mapping pkg R/ file paths → character vector of
#'                  function names to sync. Defaults to the built-in .FUNC_MAP.
#'
#' @return Invisibly, a data.frame with columns: pkg_file, fn_name, status.
#'
#' @export
sync_app_to_pkg <- function(
    app_file  = NULL,
    pkg_root  = NULL,
    rebuild   = FALSE,
    dry_run   = FALSE,
    func_map  = .FUNC_MAP
) {
  ## Defaults when called from CLI or interactively without args
  if (is.null(app_file)) stop("app_file must be specified.")
  if (is.null(pkg_root)) stop("pkg_root must be specified.")

  app_file <- normalizePath(app_file, mustWork = TRUE)
  pkg_root <- normalizePath(pkg_root, mustWork = TRUE)

  cat("\n── GS4PB App → Package Sync ───────────────────────────────────────────\n")
  cat("  App file : ", app_file, "\n")
  cat("  Pkg root : ", pkg_root, "\n")
  if (dry_run) cat("  [DRY RUN] No files will be written.\n")
  cat("\n")

  ## Load app monolith into a clean environment
  cat("  Sourcing app monolith...\n")
  app_env <- new.env(parent = emptyenv())
  # Suppress messages/plots from the app source
  suppressMessages(suppressWarnings(
    tryCatch(source(app_file, local = app_env, echo = FALSE),
             error = function(e) {
               message("  [!] source() failed: ", conditionMessage(e),
                       "\n      Falling back to text-only parsing.")
             })
  ))
  app_lines <- readLines(app_file, warn = FALSE)

  ## Process each pkg file
  results <- list()

  for (rel_path in names(func_map)) {
    pkg_file <- file.path(pkg_root, rel_path)
    if (!file.exists(pkg_file)) {
      message("  [missing] Package file not found: ", rel_path)
      next
    }
    pkg_lines <- readLines(pkg_file, warn = FALSE)
    fn_names  <- func_map[[rel_path]]
    file_changed <- FALSE

    cat(sprintf("  %-40s\n", rel_path))

    for (fn in fn_names) {
      new_body <- .extract_fn_from_lines(app_lines, fn)

      if (is.null(new_body)) {
        cat(sprintf("    %-40s [not in app]\n", fn))
        results[[length(results) + 1L]] <- data.frame(
          pkg_file = rel_path, fn_name = fn, status = "not_in_app",
          stringsAsFactors = FALSE)
        next
      }

      res <- .replace_fn_in_lines(pkg_lines, fn, new_body)

      if (!res$changed) {
        cat(sprintf("    %-40s [unchanged]\n", fn))
        results[[length(results) + 1L]] <- data.frame(
          pkg_file = rel_path, fn_name = fn, status = "unchanged",
          stringsAsFactors = FALSE)
      } else {
        cat(sprintf("    %-40s [CHANGED]\n", fn))
        pkg_lines    <- res$lines
        file_changed <- TRUE
        results[[length(results) + 1L]] <- data.frame(
          pkg_file = rel_path, fn_name = fn, status = "changed",
          stringsAsFactors = FALSE)
      }
    }

    if (file_changed && !dry_run) {
      ## Final gate: the whole rewritten file must parse before it replaces
      ## working code on disk. Per-function checks can still leave the file as
      ## a whole broken (e.g. an extraction that swallows a neighbouring
      ## definition), so validate the end result and refuse to write otherwise.
      parses <- isTRUE(tryCatch({ parse(text = pkg_lines); TRUE },
                                error = function(e) {
                                  warning("  [!] NOT WRITING ", rel_path,
                                          " — result would not parse: ",
                                          conditionMessage(e))
                                  FALSE
                                }))
      if (!parses) {
        cat(sprintf("    => SKIPPED %s (result would not parse)\n", rel_path))
      } else {
        writeLines(pkg_lines, pkg_file)
        cat(sprintf("    => wrote %s\n", rel_path))
      }
    }
    cat("\n")
  }

  ## Summary
  report <- do.call(rbind, results)
  n_changed   <- sum(report$status == "changed")
  n_unchanged <- sum(report$status == "unchanged")
  n_missing   <- sum(report$status == "not_in_app")

  cat(sprintf("  Summary: %d changed, %d unchanged, %d not in app\n",
              n_changed, n_unchanged, n_missing))

  if (n_changed > 0 && !dry_run && rebuild) {
    cat("\n  Running devtools::document() + devtools::build()...\n")
    old_wd <- setwd(pkg_root)
    on.exit(setwd(old_wd))
    devtools::document()
    devtools::build(vignettes = TRUE)
  } else if (n_changed > 0 && !dry_run) {
    cat("\n  To rebuild the package:\n")
    cat(sprintf('    setwd("%s"); devtools::document(); devtools::build()\n', pkg_root))
  }

  cat("── Done ────────────────────────────────────────────────────────────────\n\n")
  invisible(report)
}

## ── CLI entry point ──────────────────────────────────────────────────────────

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  .get_arg <- function(flag, default = NULL) {
    idx <- which(args == flag)
    if (length(idx) && idx < length(args)) args[[idx + 1L]] else default
  }
  .has_flag <- function(flag) any(args == flag)

  app_file <- .get_arg("--app")
  pkg_root <- .get_arg("--pkg",
                        default = normalizePath(
                          system.file(package = "GS4PB"), mustWork = FALSE))
  rebuild  <- .has_flag("--rebuild")
  dry_run  <- .has_flag("--dry-run")

  if (is.null(app_file)) {
    cat("Usage: Rscript sync_app_to_pkg.R --app <monolith.R> --pkg <pkg_root> [--rebuild] [--dry-run]\n")
    quit(status = 1)
  }

  sync_app_to_pkg(
    app_file = app_file,
    pkg_root = pkg_root,
    rebuild  = rebuild,
    dry_run  = dry_run
  )
}
