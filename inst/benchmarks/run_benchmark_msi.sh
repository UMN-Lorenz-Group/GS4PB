#!/bin/bash
# =============================================================================
# run_benchmark_msi.sh — submit benchmark_ME_models.R on MSI (SLURM)
#
#   sbatch run_benchmark_msi.sh
#
# Sizing below comes from measured per-fit timings on PYT23 (1531 obs, 996
# genotypes, 7 environments, 923 markers) scaled to the configured iteration
# count, and from the memory the kernels actually occupy. See "SIZING" at the
# bottom before changing --cpus-per-task, --mem or --time.
#
# PREREQUISITES
#   1. GS4PB installed in your R library on MSI, from a source tree that has
#      BOTH the benchmark script and the importFrom(dplyr,"%>%") fix. At the
#      time of writing those are local commits, not on GitHub: install from a
#      copy of the working tree, not from a fresh clone of the remote.
#   2. GS4PB's NAMESPACE imports rJava, so *loading* the package starts a JVM
#      even when the data is already prepared. Java must work. It also pulls
#      EnvRtype, STPGA, golem, shiny, plotly, DT, reticulate, vcfR, NAM,
#      rrBLUP and more — that install is the real work here, not the benchmark.
#      rTASSEL is NOT in the import list: it is needed only to build
#      prepared_data.rds, so staging that file avoids it entirely.
#   3. prepared_data.rds staged in $OUTDIR (see "DATA" below), OR rTASSEL
#      working so the job can build it.
# The preflight below checks 1-3 and exits before burning allocation.
# =============================================================================
#SBATCH --job-name=gs4pb_me_bench
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=24g
#SBATCH --time=12:00:00
#SBATCH --partition=msismall
#SBATCH --mail-type=END,FAIL
#SBATCH --output=gs4pb_me_bench_%j.out
#SBATCH --error=gs4pb_me_bench_%j.err

set -euo pipefail

# ---- Paths -----------------------------------------------------------------
# Edit these three. Everything else is derived.
PROJECT=/home/lorenzgr/$USER/GS4PB            # where the benchmark script lives
DATADIR=/home/lorenzgr/shared/MultiEnv_Trial_Data   # VCF + phenotype CSV
OUTDIR=/scratch.global/$USER/gs4pb_me_bench   # results; scratch, not home

SCRIPT="$PROJECT/inst/benchmarks/benchmark_ME_models.R"
mkdir -p "$OUTDIR"

# ---- Modules ---------------------------------------------------------------
# Match the vignette's HPC section. Java is needed only when the job has to
# build prepared_data.rds itself (rTASSEL); harmless otherwise.
module load R/4.4.2-openblas-rocky8
module load udunits gdal proj netcdf geos
# rTASSEL starts a JVM; give it room and a writable cache.
export JAVA_TOOL_OPTIONS="-Xmx30g"
export XDG_CACHE_HOME="$OUTDIR/.cache"
mkdir -p "$XDG_CACHE_HOME"

# ---- Threading -------------------------------------------------------------
# One BLAS thread per worker. Without this, 20 foreach workers x 20 BLAS
# threads oversubscribe the allocation and the job runs slower than serial.
# The R script sets these too; belt and braces, since the module may export
# its own defaults.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# ---- Benchmark configuration ----------------------------------------------
export GS4PB_BENCH_DATA="$DATADIR"
export GS4PB_BENCH_OUT="$OUTDIR"
export GS4PB_BENCH_WORKERS="${SLURM_CPUS_PER_TASK:-20}"
export GS4PB_BENCH_SCHEMES="CV1,CV2"   # CV00 is infeasible on this trial design
export GS4PB_BENCH_FOLDS=5
export GS4PB_BENCH_REPS=3
export GS4PB_BENCH_ITER=6000
# export GS4PB_BENCH_MODELS="BGGE:MDs_GB,BGLR:BayesB_GxE,bWGR:mkr"  # subset

echo "host=$(hostname) job=${SLURM_JOB_ID:-none} cpus=${SLURM_CPUS_PER_TASK:-?} mem=${SLURM_MEM_PER_NODE:-?}"

# ---- Preflight -------------------------------------------------------------
# Fail in seconds rather than after an hour of queueing plus an hour of work.
echo "-- preflight --"
Rscript -e '
  ok <- TRUE
  need <- c("GS4PB","BGGE","BGLR","bWGR","foreach","doParallel","ggplot2")
  for (p in need) {
    v <- tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
    cat(sprintf("  %-11s %s\n", p, if (is.na(v)) "MISSING" else v))
    if (is.na(v)) ok <- FALSE
  }
  cat(sprintf("  %-11s %s\n", "rTASSEL",
      tryCatch(as.character(packageVersion("rTASSEL")), error = function(e) "absent (fine if prepared_data.rds is staged)")))
  # Loading GS4PB starts the JVM via rJava; prove it works now.
  jvm <- tryCatch({ suppressMessages(library(GS4PB)); "loaded" },
                  error = function(e) paste("FAILED -", conditionMessage(e)))
  cat(sprintf("  library(GS4PB): %s\n", jvm))
  if (!identical(jvm, "loaded")) ok <- FALSE
  # The fold splitter and the pipe fix are what make this run reproducible.
  cat(sprintf("  getTstIndices_CV_Data present: %s\n",
      exists("getTstIndices_CV_Data", envir = asNamespace("GS4PB"), inherits = FALSE)))
  pipe_ok <- "%>%" %in% getNamespaceImports("GS4PB")$dplyr
  cat(sprintf("  importFrom(dplyr,\"%%>%%\") present: %s%s\n", pipe_ok,
      if (isTRUE(pipe_ok)) "" else "   <- stale GS4PB: reinstall from the fixed source"))
  if (!isTRUE(pipe_ok)) ok <- FALSE
  prep <- file.path(Sys.getenv("GS4PB_BENCH_OUT"), "prepared_data.rds")
  cat(sprintf("  prepared_data.rds: %s\n",
      if (file.exists(prep)) "staged" else "absent - job will build it (needs rTASSEL)"))
  if (!file.exists(prep) && !requireNamespace("rTASSEL", quietly = TRUE)) {
    cat("  ERROR: no prepared_data.rds and no rTASSEL; nothing can build the data.\n"); ok <- FALSE
  }
  if (!ok) quit(status = 1)
' || { echo "preflight FAILED — not starting the benchmark"; exit 1; }

echo "start $(date)"
Rscript "$SCRIPT"
echo "end   $(date)"

# =============================================================================
# DATA — two ways to get prepared_data.rds into $OUTDIR
#
#   (a) Let the job build it. Needs rTASSEL + Java on MSI and takes a few
#       minutes (VCF read, filtering, LD-KNNi imputation). It is cached, so
#       only the first job pays.
#   (b) Build it once elsewhere (laptop, Docker) and copy it over:
#         scp prepared_data.rds $USER@agate.msi.umn.edu:$OUTDIR/
#       The job then never touches rTASSEL or Java. Prefer this if the
#       rTASSEL install on MSI is not already known-good.
#
# SIZING — why these numbers
#
#   Memory. Kernels are dense n x n doubles, n = 1531 observations, so one
#   kernel is 17.9 MB. The parent builds GB and GK up front (~359 MB total).
#   The costly fit is BGGE MDe: 8 kernels, each eigen-decomposed into U and tU,
#   ~429 MB resident. 20 such fits at once is ~8.6 GB, plus the parent, plus R
#   itself. --mem=24g leaves comfortable headroom; 12g would be tight and 8g
#   would OOM partway through, after hours of work. Workers are forked, so the
#   359 MB parent is shared copy-on-write, not multiplied.
#
#   Time. Measured at 400 iterations: BGGE MDs 128 s, BGLR RKHS MDs 117 s,
#   BGGE MM 72 s, marker models ~3-6 s, bWGR ~2 s. Kernel MCMC scales linearly
#   in iterations, so at 6000 iterations one (scheme, rep, fold) round of all
#   25 models is ~4.6 CPU-hours. CV1 and CV2 at 5 folds x 3 reps is 30 rounds
#   = ~138 CPU-hours, ~7 h wall on 20 cores with perfect packing. Packing is
#   not perfect — the long kernel fits dominate the tail — so --time=12:00:00.
#
#   Scaling up. Doubling iterations to 12000 roughly doubles both. If you add
#   CV0 (7 more folds per rep), add ~50% to CV1+CV2 time.
#
# RESUMING — there is no checkpointing. If the job is killed, the only thing
# preserved is prepared_data.rds. To make a long run restartable, submit it as
# several jobs over disjoint model subsets via GS4PB_BENCH_MODELS, each with
# its own GS4PB_BENCH_OUT, and combine the benchmark_*.csv files afterwards.
# =============================================================================
