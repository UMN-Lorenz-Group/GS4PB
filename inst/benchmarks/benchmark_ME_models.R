## =============================================================================
## benchmark_ME_models.R
##
## Benchmark multi-environment genomic prediction models on the GS4PB ME trial
## dataset, under identical cross-validation folds:
##
##   BGGE          MM / MDs / MDe            (the app's current backend; GB and GK)
##   BGLR RKHS     the same kernels            (same models, different sampler)
##   BGLR marker   BRR / BayesA / BayesB / BayesC / BL, main-effect and marker x
##                 environment (Lopez-Cruz et al. 2015)
##   bWGR          mrr / mkr / MRR3-XFA / MEGA / SEM ("megavariate": one column
##                 per environment, Xavier & Habier 2022; Xavier et al. 2024)
##
## Data preparation follows vignettes/GS4PB_ME_Workflow.Rmd. The model fitters
## are written out here rather than routed through getME_CV(), so that every
## backend sees the same masked cells, the same kernels and the same fixed
## effects, and so that the benchmark can be read top to bottom.
##
## Run:        Rscript benchmark_ME_models.R
## Quick test: GS4PB_BENCH_QUICK=1 Rscript benchmark_ME_models.R
##             (3 folds, 1 rep, 400-iteration chains, 8 models: ~7 min on 4 cores)
##
## Sizing a full run. Measured on PYT23 (1531 obs, 996 lines, 7 envs, 923
## markers) at 400 iterations: BGGE MDs 128 s, BGLR RKHS MDs 117 s, BGGE MM
## 72 s, BGLR BayesB_GxE 6 s, BRR 3 s, bWGR < 7 s. Kernel MCMC scales linearly
## in iterations, so at the default 6000 one (scheme, rep, fold) round of all 25
## models is roughly 5 CPU-hours, and the default 2 schemes x 3 reps x 5 folds
## is ~150 CPU-hours: an HPC job, not a laptop one. Override without editing:
##   GS4PB_BENCH_WORKERS=20 GS4PB_BENCH_ITER=6000 GS4PB_BENCH_REPS=3 \
##   GS4PB_BENCH_FOLDS=5 GS4PB_BENCH_SCHEMES=CV1,CV2 GS4PB_BENCH_MODELS=BGGE:MDs_GB,bWGR:mkr
##
## Outputs (in cfg$out_dir):
##   prepared_data.rds            cached genotype/phenotype prep (delete to redo)
##   benchmark_predictions.csv    every test-cell prediction, all models
##   benchmark_metrics_by_env.csv cor / RMSE per rep x fold x model x env
##   benchmark_summary.csv        mean +/- sd across folds, per model x env
##   benchmark_overall.csv        one row per model: mean within-env cor, runtime
##   benchmark_runtime.csv        seconds per fit
##   benchmark_cor_by_model.png   boxplot of within-env correlation by model
##   benchmark_results.rds        everything above as one R object
## =============================================================================

suppressPackageStartupMessages({
  library(GS4PB)
  library(BGGE)
  library(BGLR)
  library(bWGR)
  library(foreach)
  library(doParallel)
})

QUICK <- identical(Sys.getenv("GS4PB_BENCH_QUICK"), "1")

## ---- 0. Configuration -------------------------------------------------------

cfg <- list(
  ## Paths. GS4PB_BENCH_DATA overrides the data directory (e.g. a Docker mount).
  data_dir = Sys.getenv("GS4PB_BENCH_DATA",
                        unset = "C:/Users/vramasub/Desktop/UMN_GIT/Genotyping_Data/MultiEnv_Trial_Data"),
  vcf_file   = "CB2019_Lorenz_PYT2023_mergedVCF_SNPs_multiAlleleAny.vcf",
  pheno_file = "PYT23_BLUEs_Locwise_4Traits_V3.csv",
  out_dir    = Sys.getenv("GS4PB_BENCH_OUT", unset = "benchmark_ME_out"),

  ## Phenotype columns (PYT23 layout: StrLocTst, Strain, Test, Loc, Year, traits)
  trait      = "Yield",
  strain_col = "Strain",
  env_col    = "Loc",
  loc_col    = "Loc",
  yr_col     = "Year",
  id_col     = "StrLocTst",

  ## Genotype QC (vignette section 3)
  MAF_thresh      = 0.05,
  site_call_rate  = 0.80,
  taxa_min_called = 0.60,
  ldknni_l = 30L, ldknni_k = 10L,

  ## Cross-validation, split by the package's getTstIndices_CV_Data():
  ##   CV1  novel genotypes         CV2  novel cells (each genotype stays observed)
  ##   CV0  novel environments      CV00 novel genotypes x novel environments
  ## CV0/CV00 use one fold per environment; k_folds applies to CV1/CV2.
  cv_schemes = strsplit(Sys.getenv("GS4PB_BENCH_SCHEMES",
                          unset = if (QUICK) "CV1" else "CV1,CV2"), ",")[[1]],
  k_folds    = as.integer(Sys.getenv("GS4PB_BENCH_FOLDS", unset = if (QUICK) "3" else "5")),
  n_reps     = as.integer(Sys.getenv("GS4PB_BENCH_REPS",  unset = if (QUICK) "1" else "3")),
  seed       = 2026L,

  ## MCMC settings shared by BGGE and BGLR so the samplers are compared fairly.
  ## burnIn is a quarter of nIter; thin fixed at 5.
  mcmc = local({
    it <- as.integer(Sys.getenv("GS4PB_BENCH_ITER", unset = if (QUICK) "400" else "6000"))
    list(nIter = it, burnIn = as.integer(it / 4), thin = 5L)
  }),

  ## Gaussian kernel (matches get_kinship_kernels defaults)
  gk_bandwidth = 1, gk_quantil = 0.5,

  ## Which models to run. Names are "<family>:<variant>"; see MODEL_SPECS.
  ## GS4PB_BENCH_MODELS (comma-separated) overrides; otherwise a representative
  ## subset in quick mode, everything otherwise.
  models = local({
    m <- Sys.getenv("GS4PB_BENCH_MODELS", unset = "")
    if (nzchar(m)) strsplit(m, ",")[[1]]
    else if (QUICK) c("BGGE:MM_GB", "BGGE:MDs_GB", "BGLR_RKHS:MDs_GB",
                      "BGLR:BRR_G", "BGLR:BayesB_GxE", "bWGR:mrr", "bWGR:mkr",
                      "bWGR:MRR3_XFA")
    else NULL   # NULL = every entry in MODEL_SPECS
  }),

  n_workers = max(1L, min(parallel::detectCores() - 1L,
                          as.integer(Sys.getenv("GS4PB_BENCH_WORKERS", unset = "8"))))
)

dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("== GS4PB ME benchmark %s==\n  out: %s | workers: %d | schemes: %s | folds: %d x reps: %d\n",
            if (QUICK) "(QUICK) " else "", normalizePath(cfg$out_dir), cfg$n_workers,
            paste(cfg$cv_schemes, collapse = ","), cfg$k_folds, cfg$n_reps))

## Every model the script knows. `family` selects the fitter; the remaining
## fields are its arguments.
MODEL_SPECS <- list(
  ## BGGE, the backend behind getME_CV(): main effect (MM), plus a single GxE
  ## kernel (MDs) or one per environment (MDe). GB = linear, GK = Gaussian.
  "BGGE:MM_GB"  = list(family = "BGGE", model = "MM",  kernel = "GB"),
  "BGGE:MDs_GB" = list(family = "BGGE", model = "MDs", kernel = "GB"),
  "BGGE:MDe_GB" = list(family = "BGGE", model = "MDe", kernel = "GB"),
  "BGGE:MM_GK"  = list(family = "BGGE", model = "MM",  kernel = "GK"),
  "BGGE:MDs_GK" = list(family = "BGGE", model = "MDs", kernel = "GK"),
  "BGGE:MDe_GK" = list(family = "BGGE", model = "MDe", kernel = "GK"),
  ## BGLR fitting the identical kernels as RKHS terms
  "BGLR_RKHS:MM_GB"  = list(family = "BGLR_RKHS", model = "MM",  kernel = "GB"),
  "BGLR_RKHS:MDs_GB" = list(family = "BGLR_RKHS", model = "MDs", kernel = "GB"),
  "BGLR_RKHS:MDe_GB" = list(family = "BGLR_RKHS", model = "MDe", kernel = "GB"),
  "BGLR_RKHS:MDs_GK" = list(family = "BGLR_RKHS", model = "MDs", kernel = "GK"),
  ## BGLR marker-based priors. _G = markers + env fixed; _GxE adds one
  ## environment-specific marker block per environment (Lopez-Cruz 2015).
  "BGLR:BRR_G"      = list(family = "BGLR_marker", prior = "BRR",    gxe = FALSE),
  "BGLR:BRR_GxE"    = list(family = "BGLR_marker", prior = "BRR",    gxe = TRUE),
  "BGLR:BayesA_G"   = list(family = "BGLR_marker", prior = "BayesA", gxe = FALSE),
  "BGLR:BayesA_GxE" = list(family = "BGLR_marker", prior = "BayesA", gxe = TRUE),
  "BGLR:BayesB_G"   = list(family = "BGLR_marker", prior = "BayesB", gxe = FALSE),
  "BGLR:BayesB_GxE" = list(family = "BGLR_marker", prior = "BayesB", gxe = TRUE),
  "BGLR:BayesC_G"   = list(family = "BGLR_marker", prior = "BayesC", gxe = FALSE),
  "BGLR:BayesC_GxE" = list(family = "BGLR_marker", prior = "BayesC", gxe = TRUE),
  "BGLR:BL_G"       = list(family = "BGLR_marker", prior = "BL",     gxe = FALSE),
  "BGLR:BL_GxE"     = list(family = "BGLR_marker", prior = "BL",     gxe = TRUE),
  ## bWGR multivariate solvers on the genotype x environment matrix. These
  ## treat each environment as a trait, so GxE is the between-column
  ## covariance; there is no separate main-effect term to switch on or off.
  "bWGR:mrr"      = list(family = "bWGR", solver = "mrr"),
  "bWGR:mkr"      = list(family = "bWGR", solver = "mkr"),
  "bWGR:MRR3_XFA" = list(family = "bWGR", solver = "MRR3_XFA"),
  "bWGR:MEGA"     = list(family = "bWGR", solver = "MEGA"),
  "bWGR:SEM"      = list(family = "bWGR", solver = "SEM")
)
if (is.null(cfg$models)) cfg$models <- names(MODEL_SPECS)
unknown <- setdiff(cfg$models, names(MODEL_SPECS))
if (length(unknown)) stop("Unknown model(s): ", paste(unknown, collapse = ", "))

## ---- 1. Data preparation (vignette sections 3-4), cached -------------------

prep_path <- file.path(cfg$out_dir, "prepared_data.rds")

prepare_data <- function(cfg) {
  vcf_path   <- file.path(cfg$data_dir, cfg$vcf_file)
  pheno_path <- file.path(cfg$data_dir, cfg$pheno_file)
  for (f in c(vcf_path, pheno_path)) if (!file.exists(f)) stop("Not found: ", f)

  cat("[prep] VCF -> rTASSEL\n")
  tas     <- getTasObj(vcf_path)
  geno_df <- getGenoTas_to_DF(tas)
  n_markers <- nrow(geno_df)
  cat("      ", getGenoQCStats(geno_df), "\n")

  cat("[prep] filter sites / taxa\n")
  tas_f1 <- getFilteredSitesGenoData(tas,
              siteMinCnt = round(cfg$site_call_rate * n_markers), MAF = cfg$MAF_thresh)
  tas_f2 <- getFilteredTaxaGenoData(tas_f1, MinNotMissing = cfg$taxa_min_called)

  ## LD-KNNi. Not getImputedData_Num(): rTASSEL's imputeNumeric() is a no-op
  ## on a discrete genotype table and leaves every missing call in place.
  cat("[prep] LD-KNNi imputation\n")
  geno_imp <- getGenoTas_to_DF(getImputedData_LDKNNI(tas_f2, l = cfg$ldknni_l, k = cfg$ldknni_k))
  cat("      ", getGenoQCStats(geno_imp), "\n")

  cat("[prep] phenotype\n")
  ph <- read.csv(pheno_path, stringsAsFactors = FALSE)
  id_cols <- setdiff(colnames(ph), cfg$trait)
  std <- stdizePhVarNames(ph, IDColsME = id_cols, StrainME = cfg$strain_col,
                          EnvVarME = cfg$env_col, LocMEColVar = cfg$loc_col, YrMEColVar = cfg$yr_col)
  ph_std <- std[[1]]; id_cols_mod <- std[[2]]
  ph_fmt <- getPhenoMEData(PhenoME = ph_std, TraitME = cfg$trait, nSelTraits = 1L,
                           IDColsME = id_cols_mod, StrainME = cfg$strain_col)
  merged <- getMergedDataME(phData = ph_fmt, genoImp = geno_imp, TargetIDs = NULL)
  DT_raw <- as.data.frame(merged[[1]][[1]])
  X_raw  <- as.matrix(merged[[2]][[1]])

  ## Long table: one row per genotype x environment. Ordered by env then gid so
  ## BGGE's block-diagonal kernels and its `ne` argument line up without any
  ## reordering downstream.
  DT <- data.frame(uniqID = as.character(DT_raw[[cfg$id_col]]),
                   gid    = as.character(DT_raw$Strain),
                   env    = as.character(DT_raw$Env),
                   y      = as.numeric(DT_raw[[cfg$trait]]),
                   stringsAsFactors = FALSE)
  DT <- DT[!is.na(DT$y), ]
  DT <- DT[order(DT$env, DT$gid), ]
  rownames(DT) <- NULL
  if (any(duplicated(DT[c("gid", "env")])))
    stop("Replicated (gid, env) cells; this benchmark expects one BLUE per cell.")

  ## Marker matrix: gid x marker, centred. Keep only genotypes present in DT.
  if (is.null(rownames(X_raw)) || !all(unique(DT$gid) %in% rownames(X_raw))) {
    miss <- setdiff(unique(DT$gid), rownames(X_raw))
    stop("Genotypes in phenotype table but not in marker matrix: ", paste(head(miss), collapse = ", "))
  }
  gids <- sort(unique(DT$gid))
  X <- X_raw[gids, , drop = FALSE]
  storage.mode(X) <- "numeric"
  X <- X[, apply(X, 2, function(v) length(unique(v)) > 1), drop = FALSE]  # drop fixed markers
  Xc <- scale(X, center = TRUE, scale = FALSE)
  attr(Xc, "scaled:center") <- NULL

  list(DT = DT, Xc = Xc, gids = gids, envs = sort(unique(DT$env)),
       n_markers_raw = n_markers, n_markers = ncol(Xc))
}

if (file.exists(prep_path)) {
  cat("[prep] using cached", prep_path, "\n")
  dat <- readRDS(prep_path)
} else {
  dat <- prepare_data(cfg)
  saveRDS(dat, prep_path)
}
cat(sprintf("[data] %d obs | %d genotypes | %d environments (%s) | %d markers\n",
            nrow(dat$DT), length(dat$gids), length(dat$envs),
            paste(dat$envs, collapse = ","), dat$n_markers))

## ---- 2. Design matrices and kernels (shared by every model) -----------------

build_designs <- function(dat, cfg) {
  DT <- dat$DT
  env_f <- factor(DT$env, levels = dat$envs)
  gid_f <- factor(DT$gid, levels = dat$gids)
  Ze <- model.matrix(~ 0 + env_f); colnames(Ze) <- dat$envs
  Zg <- model.matrix(~ 0 + gid_f); colnames(Zg) <- dat$gids
  ne <- as.vector(table(env_f))      # observations per environment, in env order

  ## Genotype relationship matrices (gid x gid)
  p <- ncol(dat$Xc)
  G_GB <- tcrossprod(dat$Xc) / p
  D    <- as.matrix(dist(dat$Xc))^2
  G_GK <- exp(-cfg$gk_bandwidth * D / quantile(D, probs = cfg$gk_quantil))
  dimnames(G_GB) <- dimnames(G_GK) <- list(dat$gids, dat$gids)

  ## Observation-level kernels, as get_kinship_kernels() builds them:
  ##   MM : K_G = Zg G Zg'
  ##   MDs: + K_G o (Ze Ze')           one GxE kernel, block-diagonal by env
  ##   MDe: + K_G o (Ze_i Ze_i')       one kernel per environment
  obs_kernels <- function(G) {
    KG  <- Zg %*% G %*% t(Zg)
    KGE <- KG * tcrossprod(Ze)
    KGE_env <- lapply(seq_along(dat$envs), function(i) KG * tcrossprod(Ze[, i, drop = FALSE]))
    names(KGE_env) <- dat$envs
    list(KG = KG, KGE = KGE, KGE_env = KGE_env)
  }
  list(Ze = Ze, Zg = Zg, ne = ne,
       GB = obs_kernels(G_GB), GK = obs_kernels(G_GK),
       G_GB = G_GB, G_GK = G_GK,
       Xn = Zg %*% dat$Xc)                     # n x p marker matrix in obs order
}

cat("[kernels] building designs and kernels ...\n")
des <- build_designs(dat, cfg)

## ---- 3. Cross-validation folds ---------------------------------------------

## The folds come from the package's own splitter, getTstIndices_CV_Data(),
## the one getME_CV() uses and that has been validated against the app, so a
## result here is comparable to a result from the app. It is internal (not
## exported), hence the `:::`. Its CV2 is not "random cells": each genotype's
## rows are spread across folds so every genotype stays observed somewhere,
## and single-observation genotypes are never masked.
##
## Returns, for each replicate, a list of integer vectors: the DT rows held out
## in each fold. Every model in a run sees exactly these rows as test cells.
make_folds <- function(DT, scheme, k, n_reps, seed) {
  stopifnot(scheme %in% c("CV1", "CV2", "CV0", "CV00"))
  ## CV0 folds over environments; one fold per environment makes it a genuine
  ## leave-one-environment-out whatever k_folds says. CV00 keeps k_folds: it
  ## masks a share of the genotypes that appear in BOTH the held-out and the
  ## training environments, and the splitter refuses a fold with too few such
  ## cells — on a thinly connected trial (each line in 1-2 environments) that
  ## can be every fold, in which case the scheme is skipped below.
  if (scheme == "CV0") k <- length(unique(DT$env))
  sp <- utils::capture.output(
    out <- GS4PB:::getTstIndices_CV_Data(DT, CVMet = scheme, nIter = n_reps, k = k,
                                         id_col = "gid", env_col = "env", y_col = "y",
                                         base_seed = seed))
  ## Drop empty folds (CV2 can produce them when few rows are maskable).
  lapply(out$test_idx, function(rep_folds) Filter(length, rep_folds))
}

## ---- 4. Model fitters --------------------------------------------------------
##
## Each takes the masked response (NA at test cells) and returns predictions for
## every row of DT, in DT order. Test predictions are then picked out by index,
## so the fitters never see which rows are being scored.

kernel_list <- function(des, kernel, model) {
  K <- des[[kernel]]
  out <- list(KG = list(Kernel = K$KG, Type = "D"))
  if (model == "MDs") out$KGE <- list(Kernel = K$KGE, Type = "BD")
  if (model == "MDe") for (e in names(K$KGE_env))
    out[[paste0("KGE_", e)]] <- list(Kernel = K$KGE_env[[e]], Type = "BD")
  out
}

fit_BGGE <- function(y_masked, spec, des, cfg) {
  K <- kernel_list(des, spec$kernel, spec$model)
  ## Full ~0+env design, as the package passes it: BGGE has no separate mu.
  fm <- BGGE::BGGE(y = y_masked, K = K, XF = des$Ze, ne = des$ne,
                   ite = cfg$mcmc$nIter, burn = cfg$mcmc$burnIn, thin = cfg$mcmc$thin,
                   verbose = FALSE, tol = 1e-10, R2 = 0.5)
  as.numeric(fm$yHat)
}

fit_BGLR_RKHS <- function(y_masked, spec, des, cfg) {
  K <- kernel_list(des, spec$kernel, spec$model)
  ETA <- c(list(env = list(X = des$Ze[, -1, drop = FALSE], model = "FIXED")),  # BGLR fits mu
           lapply(K, function(k) list(K = k$Kernel, model = "RKHS")))
  fm <- BGLR::BGLR(y = y_masked, ETA = ETA, nIter = cfg$mcmc$nIter,
                   burnIn = cfg$mcmc$burnIn, thin = cfg$mcmc$thin,
                   saveAt = tempfile("bglr_"), verbose = FALSE)
  as.numeric(fm$yHat)
}

fit_BGLR_marker <- function(y_masked, spec, des, cfg) {
  ETA <- list(env = list(X = des$Ze[, -1, drop = FALSE], model = "FIXED"),
              G   = list(X = des$Xn, model = spec$prior))
  if (isTRUE(spec$gxe)) {
    ## Marker x environment: a copy of the marker matrix zeroed outside each
    ## environment, so each environment gets its own deviation from the main
    ## marker effects, all under the same prior.
    for (i in seq_len(ncol(des$Ze)))
      ETA[[paste0("GxE_", colnames(des$Ze)[i])]] <-
        list(X = des$Xn * des$Ze[, i], model = spec$prior)
  }
  fm <- BGLR::BGLR(y = y_masked, ETA = ETA, nIter = cfg$mcmc$nIter,
                   burnIn = cfg$mcmc$burnIn, thin = cfg$mcmc$thin,
                   saveAt = tempfile("bglr_"), verbose = FALSE)
  as.numeric(fm$yHat)
}

fit_bWGR <- function(y_masked, spec, des, cfg, dat) {
  ## Wide genotype x environment matrix; NA marks both unobserved and test cells.
  Y <- matrix(NA_real_, length(dat$gids), length(dat$envs), dimnames = list(dat$gids, dat$envs))
  Y[cbind(dat$DT$gid, dat$DT$env)] <- y_masked
  fit <- switch(spec$solver,
    mrr      = bWGR::mrr(Y, dat$Xc),
    mkr      = bWGR::mkr(Y, des$G_GB),
    MRR3_XFA = bWGR::MRR3(Y, dat$Xc, XFA = TRUE, NumXFA = min(3L, ncol(Y) - 1L)),
    MEGA     = bWGR::MEGA(Y, dat$Xc),
    SEM      = bWGR::SEM(Y, dat$Xc),
    stop("Unknown bWGR solver: ", spec$solver))
  hat <- fit$hat
  dimnames(hat) <- dimnames(Y)
  as.numeric(hat[cbind(dat$DT$gid, dat$DT$env)])
}

fit_one <- function(y_masked, spec, des, cfg, dat) {
  switch(spec$family,
    BGGE        = fit_BGGE(y_masked, spec, des, cfg),
    BGLR_RKHS   = fit_BGLR_RKHS(y_masked, spec, des, cfg),
    BGLR_marker = fit_BGLR_marker(y_masked, spec, des, cfg),
    bWGR        = fit_bWGR(y_masked, spec, des, cfg, dat),
    stop("Unknown family: ", spec$family))
}

## ---- 5. Run every (scheme, rep, fold, model) ---------------------------------

fold_sets <- setNames(lapply(cfg$cv_schemes, function(sch)
  make_folds(dat$DT, sch, cfg$k_folds, cfg$n_reps, cfg$seed)), cfg$cv_schemes)
for (sch in cfg$cv_schemes) {
  n_test <- vapply(fold_sets[[sch]][[1]], length, 1L)
  if (!length(n_test)) {
    cat(sprintf("[folds] %-4s produced NO test cells — the splitter found nothing it could mask\n",
                sch), "        under this scheme on this trial design; skipping it.\n")
    next
  }
  cat(sprintf("[folds] %-4s rep 1: %d folds, test cells per fold %s (of %d)\n",
              sch, length(n_test), paste(n_test, collapse = "/"), nrow(dat$DT)))
}
empty <- vapply(fold_sets, function(f) all(lengths(f) == 0), TRUE)
fold_sets <- fold_sets[!empty]; cfg$cv_schemes <- names(fold_sets)
if (!length(cfg$cv_schemes)) stop("No CV scheme produced any test cells; nothing to run.")

tasks <- do.call(rbind, lapply(cfg$cv_schemes, function(sch) {
  folds <- fold_sets[[sch]]
  do.call(rbind, lapply(seq_along(folds), function(r)
    do.call(rbind, lapply(seq_along(folds[[r]]), function(f)
      data.frame(scheme = sch, rep = r, fold = f, model = cfg$models,
                 stringsAsFactors = FALSE)))))
}))
## CV0/CV00 hold out whole environments; a wide-format solver has nothing to
## estimate for an all-NA column, so those combinations are skipped, not failed.
skip <- tasks$scheme %in% c("CV0", "CV00") &
        vapply(tasks$model, function(m) MODEL_SPECS[[m]]$family == "bWGR", TRUE)
if (any(skip)) cat(sprintf("[run] skipping %d bWGR x CV0/CV00 tasks (whole-environment hold-out)\n", sum(skip)))
tasks <- tasks[!skip, ]

cat(sprintf("[run] %d tasks (%d models) on %d workers\n", nrow(tasks), length(cfg$models), cfg$n_workers))
registerDoParallel(cfg$n_workers)
t_all <- proc.time()

results <- foreach(i = seq_len(nrow(tasks)), .packages = c("BGGE", "BGLR", "bWGR"),
                   .errorhandling = "pass") %dopar% {
  tk   <- tasks[i, ]
  spec <- MODEL_SPECS[[tk$model]]
  test <- fold_sets[[tk$scheme]][[tk$rep]][[tk$fold]]
  y_masked <- dat$DT$y; y_masked[test] <- NA
  t0 <- proc.time()[["elapsed"]]
  pred <- tryCatch(fit_one(y_masked, spec, des, cfg, dat),
                   error = function(e) structure(conditionMessage(e), class = "fit_error"))
  secs <- proc.time()[["elapsed"]] - t0
  if (inherits(pred, "fit_error"))
    return(data.frame(tk, status = "error", error = unclass(pred), secs = secs,
                      uniqID = NA, gid = NA, env = NA, obs = NA, pred = NA))
  data.frame(tk, status = "ok", error = "", secs = secs,
             uniqID = dat$DT$uniqID[test], gid = dat$DT$gid[test], env = dat$DT$env[test],
             obs = dat$DT$y[test], pred = pred[test], row.names = NULL)
}
stopImplicitCluster()
cat(sprintf("[run] done in %.1f min\n", (proc.time() - t_all)[["elapsed"]] / 60))

## foreach with .errorhandling="pass" hands back the condition object for a
## task that failed outside the tryCatch; normalise those too.
results <- lapply(seq_along(results), function(i) {
  r <- results[[i]]
  if (inherits(r, "error")) data.frame(tasks[i, ], status = "error", error = conditionMessage(r),
                                       secs = NA, uniqID = NA, gid = NA, env = NA, obs = NA, pred = NA)
  else r
})
pred_df <- do.call(rbind, results)
failed  <- unique(pred_df[pred_df$status == "error", c("scheme", "rep", "fold", "model", "error")])
if (nrow(failed)) {
  cat("[run] FAILED fits:\n"); print(failed, row.names = FALSE)
}
pred_ok <- pred_df[pred_df$status == "ok", ]

## ---- 6. Metrics ---------------------------------------------------------------

## Within-environment predictive ability is what a breeder sees when ranking
## lines at a site, so correlations are computed per environment and then
## averaged, not pooled across environments (pooling rewards a model that
## merely gets the environment means right).
##
## Under CV0 read only the correlation. A held-out environment's fixed effect
## is unidentifiable from the training data, so every model's predictions for
## it carry an arbitrary level shift: rankings survive, RMSE does not (it was
## ~16 vs ~6 under CV1/CV2 on PYT23). The marker x environment models are hit
## harder still — their environment-specific blocks for the unseen environment
## are sampled from the prior and add noise, where the kernel models' block
## structure shrinks that term towards zero.
metric_rows <- function(d) {
  ok <- is.finite(d$obs) & is.finite(d$pred)
  n  <- sum(ok)
  data.frame(n = n,
             cor  = if (n > 2) cor(d$obs[ok], d$pred[ok]) else NA_real_,
             rmse = if (n > 0) sqrt(mean((d$obs[ok] - d$pred[ok])^2)) else NA_real_,
             secs = d$secs[1])
}
by_env <- do.call(rbind, lapply(split(pred_ok, pred_ok[c("scheme", "rep", "fold", "model", "env")], drop = TRUE),
  function(d) cbind(d[1, c("scheme", "rep", "fold", "model", "env")], metric_rows(d))))
rownames(by_env) <- NULL

agg <- function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
summary_df <- do.call(rbind, lapply(split(by_env, by_env[c("scheme", "model", "env")], drop = TRUE), function(d) {
  cbind(d[1, c("scheme", "model", "env")],
        n_fits = nrow(d),
        cor_mean = agg(d$cor)[["mean"]],  cor_sd = agg(d$cor)[["sd"]],
        rmse_mean = agg(d$rmse)[["mean"]], rmse_sd = agg(d$rmse)[["sd"]])
}))
rownames(summary_df) <- NULL

overall_df <- do.call(rbind, lapply(split(by_env, by_env[c("scheme", "model")], drop = TRUE), function(d) {
  ## one cor per fit = mean across its environments; then mean +/- sd over fits
  per_fit <- tapply(d$cor, interaction(d$rep, d$fold, drop = TRUE), mean, na.rm = TRUE)
  secs    <- tapply(d$secs, interaction(d$rep, d$fold, drop = TRUE), function(s) s[1])
  cbind(d[1, c("scheme", "model")],
        family = MODEL_SPECS[[d$model[1]]]$family,
        n_fits = length(per_fit),
        cor_mean = mean(per_fit), cor_sd = sd(per_fit),
        rmse_mean = mean(d$rmse, na.rm = TRUE),
        secs_mean = mean(secs), secs_sd = sd(secs))
}))
rownames(overall_df) <- NULL
overall_df <- overall_df[order(overall_df$scheme, -overall_df$cor_mean), ]

runtime_df <- unique(pred_ok[c("scheme", "rep", "fold", "model", "secs")])

## ---- 7. Write outputs -----------------------------------------------------------

w <- function(x, f) write.csv(x, file.path(cfg$out_dir, f), row.names = FALSE)
w(pred_df,    "benchmark_predictions.csv")
w(by_env,     "benchmark_metrics_by_env.csv")
w(summary_df, "benchmark_summary.csv")
w(overall_df, "benchmark_overall.csv")
w(runtime_df, "benchmark_runtime.csv")
saveRDS(list(cfg = cfg, model_specs = MODEL_SPECS, tasks = tasks, folds = fold_sets,
             predictions = pred_df, by_env = by_env, summary = summary_df,
             overall = overall_df, runtime = runtime_df, failed = failed,
             data_dims = dat[c("n_markers_raw", "n_markers", "envs")]),
        file.path(cfg$out_dir, "benchmark_results.rds"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  by_env$family <- vapply(by_env$model, function(m) MODEL_SPECS[[m]]$family, "")
  ord <- overall_df$model[order(overall_df$cor_mean)]
  by_env$model <- factor(by_env$model, levels = unique(ord))
  p <- ggplot(by_env, aes(x = model, y = cor, fill = family)) +
    geom_boxplot(outlier.size = 0.6) +
    facet_wrap(~ scheme, ncol = 1) +
    coord_flip() +
    labs(y = "Within-environment correlation (obs vs pred, test cells)", x = NULL,
         title = sprintf("GS4PB ME benchmark — %s, %d envs, %d genotypes, %d markers",
                         cfg$trait, length(dat$envs), length(dat$gids), dat$n_markers),
         subtitle = sprintf("%d folds x %d reps per scheme; MCMC %d/%d/%d",
                            cfg$k_folds, cfg$n_reps, cfg$mcmc$nIter, cfg$mcmc$burnIn, cfg$mcmc$thin)) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(cfg$out_dir, "benchmark_cor_by_model.png"), p,
         width = 9, height = 2.5 + 0.28 * length(cfg$models) * length(cfg$cv_schemes), dpi = 150)
}

cat("\n== Overall (mean within-env cor across fits; higher is better) ==\n")
print(format(overall_df[c("scheme", "model", "cor_mean", "cor_sd", "rmse_mean", "secs_mean")],
             digits = 3), row.names = FALSE)
cat("\nOutputs in", normalizePath(cfg$out_dir), "\n")
