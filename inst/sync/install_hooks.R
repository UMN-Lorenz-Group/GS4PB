###############################################################################
##  GS4PB — Install local git hooks
##
##  Run this once from the GS4PB package root to install the post-commit hook
##  that reminds you to sync changed pipeline functions to the app repos.
##
##  Usage (from GS4PB package root):
##    source("inst/sync/install_hooks.R")
##  OR after installing the package:
##    source(system.file("sync", "install_hooks.R", package = "GS4PB"))
###############################################################################

.install_gs4pb_hooks <- function() {

  ## Find the .git directory — works whether called from pkg root or installed
  pkg_root <- tryCatch(
    rprojroot::find_root(rprojroot::has_file("DESCRIPTION")),
    error = function(e) getwd()
  )

  git_dir  <- file.path(pkg_root, ".git")
  hooks_dir <- file.path(git_dir, "hooks")

  if (!dir.exists(git_dir)) {
    message("[!] No .git directory found at: ", pkg_root)
    message("    Run this script from inside the GS4PB git repository.")
    return(invisible(FALSE))
  }

  if (!dir.exists(hooks_dir)) dir.create(hooks_dir, recursive = TRUE)

  ## Locate hook source (works installed or from source)
  hook_src <- system.file("sync", "post-commit", package = "GS4PB")
  if (!nzchar(hook_src)) {
    hook_src <- file.path(pkg_root, "inst", "sync", "post-commit")
  }
  if (!file.exists(hook_src)) {
    message("[!] post-commit hook source not found at: ", hook_src)
    return(invisible(FALSE))
  }

  hook_dst <- file.path(hooks_dir, "post-commit")
  file.copy(hook_src, hook_dst, overwrite = TRUE)

  ## Make executable (no-op on Windows but harmless)
  tryCatch(Sys.chmod(hook_dst, "755"), error = function(e) NULL)

  message("[ok] post-commit hook installed at: ", hook_dst)
  message("     It will remind you to sync R/ changes to app repos after each commit.")
  invisible(TRUE)
}

.install_gs4pb_hooks()
