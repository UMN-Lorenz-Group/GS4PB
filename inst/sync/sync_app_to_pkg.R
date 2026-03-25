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
    "getMEPred", "stdizePhVarNames", "fitMEModels_LOFO_Pred"
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

#' Given the index of a `function(` line, find the matching closing `}` line
#' using brace counting.
.find_fn_end <- function(lines, start_idx) {
  depth <- 0L
  for (i in seq(start_idx, length(lines))) {
    depth <- depth +
      nchar(gsub("[^{]", "", lines[[i]])) -
      nchar(gsub("[^}]", "", lines[[i]]))
    if (depth == 0L && i > start_idx) return(i)
    # Handle single-line functions like  f <- function() NULL
    if (i == start_idx && !grepl("\\{", lines[[i]])) return(i)
  }
  NA_integer_
}

#' Extract a function's signature + body as a character vector of lines
#' from the monolith text lines.
.extract_fn_from_lines <- function(lines, fn_name) {
  starts <- .find_fn_start(lines, fn_name)
  if (length(starts) == 0L) return(NULL)
  start <- starts[[1L]]
  end   <- .find_fn_end(lines, start)
  if (is.na(end)) {
    warning("  [!] Could not find closing } for ", fn_name, " (start=", start, ")")
    return(NULL)
  }
  lines[start:end]
}

#' Replace a function's signature + body in pkg_lines with new_body_lines.
#' Lines before the function start (roxygen2 block) are preserved.
.replace_fn_in_lines <- function(pkg_lines, fn_name, new_body_lines) {
  starts <- .find_fn_start(pkg_lines, fn_name)
  if (length(starts) == 0L) {
    message("  [not found] ", fn_name, " — not present in package file, skipping")
    return(list(lines = pkg_lines, changed = FALSE))
  }
  start <- starts[[1L]]
  end   <- .find_fn_end(pkg_lines, start)
  if (is.na(end)) {
    warning("  [!] Could not find closing } for ", fn_name, " in pkg file")
    return(list(lines = pkg_lines, changed = FALSE))
  }

  old_body <- pkg_lines[start:end]
  if (identical(old_body, new_body_lines)) {
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
      writeLines(pkg_lines, pkg_file)
      cat(sprintf("    => wrote %s\n", rel_path))
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
