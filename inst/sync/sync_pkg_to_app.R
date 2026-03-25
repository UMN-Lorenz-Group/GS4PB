#!/usr/bin/env Rscript
###############################################################################
##  GS4PB — Package → App Sync (Phase 2 — stub)
##
##  This script is a stub for the FUTURE direction once both apps are stable
##  enough to consume the GS4PB package directly.
##
##  Phase 2 migration steps (do manually or extend this script):
##
##  Step 1 — Add GS4PB to each app's renv.lock
##    renv::install("UMN-Lorenz-Group/GS4PB")   # from GitHub
##    renv::snapshot()
##
##  Step 2 — Replace source() call in app.R
##    Replace:
##      FN <- paste(getwd(), "/GS_Pipeline_Jan_2026_FnsApp.R", sep="")
##      source(FN)
##    With:
##      library(GS4PB)
##      rTASSEL::initializeTASSEL()
##
##  Step 3 — Remove the monolith from the app directory
##    file.remove("GS_Pipeline_Jan_2026_FnsApp.R")
##
##  Step 4 — Update any direct function calls
##    All exported GS4PB functions are now available after library(GS4PB).
##    Non-exported internal functions will need to be prefixed with GS4PB:::
##    if still needed, or exported from the package.
##
##  Step 5 — Verify the app runs
##    shiny::runApp(".")
##
##  Once both apps are on library(GS4PB), there is nothing to sync —
##  the sync_app_to_pkg.R script becomes obsolete.
##
###############################################################################

message("This is a Phase 2 migration stub. See comments in this file for instructions.")
message("Current sync direction: App -> Package (sync_app_to_pkg.R)")
message("Run Phase 2 when both apps are ready to consume library(GS4PB).")
