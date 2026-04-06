###############################################################################
##  GS4PB — Reference Pipeline Script
##  Covers both SE (Single-Environment: ST and MT) and ME (Multi-Environment)
##  workflows from VCF loading through GEBV predictions.
##
##  Prerequisites:
##    install.packages("path/to/GS4PB_0.0.0.9000.tar.gz", repos = NULL,
##                     type = "source")
##    rTASSEL::initializeTASSEL()   # start the JVM before any rTASSEL calls
##
##  Replace every <...> placeholder with real paths / column names.
###############################################################################

library(GS4PB)

## ── User inputs ──────────────────────────────────────────────────────────────
vcf_path       <- "<path/to/genotype.vcf>"      # input VCF file
se_pheno_path  <- "<path/to/se_phenotype.csv>"  # SE phenotype CSV
me_pheno_path  <- "<path/to/me_phenotype.csv>"  # ME phenotype CSV (wide, one row per obs)
target_id_path <- "<path/to/target_ids.txt>"    # optional: one strain ID per line
out_dir        <- "<path/to/output_directory>"  # results directory

## SE phenotype
trait_se       <- "Yield"         # single trait column name
traits_mt      <- c("Yield", "Height", "Maturity")  # multi-trait vector

## ME phenotype column names (will be standardised internally)
strain_col     <- "Strain"
env_col        <- "Env"           # environment factor column
loc_col        <- "Location"      # location column (or NULL)
yr_col         <- "Year"          # year column (or NULL)
trait_me       <- c("Yield")      # ME traits

## ── VCF filter / imputation parameters ───────────────────────────────────────
MAF_thresh     <- 0.05   # minor allele frequency cutoff
site_min_cnt   <- 0.80   # minimum site call rate (fraction)
taxa_min_nr_miss <- 0.80 # minimum individual call rate

## ── CV parameters ────────────────────────────────────────────────────────────
k_folds  <- 5    # number of CV folds
n_iter   <- 5    # number of CV replications
cv_met   <- "CV1"  # "CV1" (by observation) or "CV2" (by factor level)


###############################################################################
## SECTION 1 — GENOTYPIC DATA PROCESSING  (shared by SE and ME)
###############################################################################

## 1a — Load VCF ---------------------------------------------------------------
# rTASSEL reads the VCF into JVM memory; getGenoTas_to_DF converts it to R.
tas        <- getTasObj(vcf_path)
geno_df    <- getGenoTas_to_DF(tas)        # tibble: SNPID, Chrom, Pos, REF, ALT, <samples>
cat(getGenoQCStats(geno_df))               # baseline QC message

## 1b — Filter markers (Filt1) ------------------------------------------------
# Remove SNPs with low call rate or below MAF threshold.
tas_filt1  <- getFilteredSitesGenoData(tas, siteMinCnt = site_min_cnt, MAF = MAF_thresh)
geno_filt1 <- getGenoTas_to_DF(tas_filt1)
cat(getGenoQCStatsFilt1(geno_df, geno_filt1,
                        missSitesTH = 1 - site_min_cnt, MAFTH = MAF_thresh))

## 1c — Filter samples (Filt2) ------------------------------------------------
# Remove individuals with high missingness.
tas_filt2  <- getFilteredTaxaGenoData(tas_filt1, MinNotMissing = taxa_min_nr_miss)
geno_filt2 <- getGenoTas_to_DF(tas_filt2)
cat(getGenoQCStatsFilt2(geno_filt1, geno_filt2, missIndTH = 1 - taxa_min_nr_miss))

## 1d — Imputation -------------------------------------------------------------
# Choose ONE method; comment out the others.

## Option A: Numeric (mean) imputation via rTASSEL
tas_imp    <- getImputedData_Num(tas_filt2, nN = 5, Dist = "Euclidean")
geno_imp   <- getGenoTas_to_DF(tas_imp)    # imputed genotype data frame

## Option B: LD-KNNi imputation via rTASSEL
# tas_imp  <- getImputedData_LDKNNI(tas_filt2, l = 50, k = 5)
# geno_imp <- getGenoTas_to_DF(tas_imp)

## Option C: AlphaPlantImpute (external tool required)
# getGenoData_API(tas_filt2)              # writes current_GenoTable.genotypes
# # ... run AlphaPlantImpute externally to produce imputed_out.genotypes ...
# geno_imp <- getImpGenoData_API(vcfIDTab)  # reads back imputed results

## 1e — (Optional) Genomic Kinship / G matrix ---------------------------------
# Useful when supplying a custom kernel (A_ext) to SE:MT or ME models.

qc_result  <- getQCFilteredM(geno_imp, maf = MAF_thresh)   # ASRgenomics QC filter
M_clean    <- qc_result$M.clean                            # cleaned numeric marker matrix

gmat_out   <- getGmat(M_clean, method = "VanRaden")        # compute G (VanRaden)
G          <- gmat_out$G

diag_out   <- getKinshipDiag(G, duplicate.thr = 0.98)      # detect duplicate individuals
cat(diag_out$msg)

# Remove duplicates (greedy: keep first in each pair, remove second)
if (!is.null(diag_out$list.duplicate) && nrow(diag_out$list.duplicate) > 0) {
  dedup_out <- removeDuplicatesFromG(G, diag_out$list.duplicate)
  G         <- dedup_out$G_clean
  cat(dedup_out$msg)
}

G_tuned    <- getTunedG(G, method = "blend", pblend = 0.05)$Gb  # tune for PD

pca_result <- getKinshipPCA(G_tuned, ncp = 20)                  # PCA for visualisation
# plot(pca_result$ind$coord[, 1], pca_result$ind$coord[, 2], ...)


###############################################################################
## SECTION 2 — SE WORKFLOW: SINGLE-ENVIRONMENT (ST and MT)
###############################################################################

## 2a — Load SE phenotype -------------------------------------------------------
# CSV must have one row per individual, with a strain ID column and trait columns.
pheno_se   <- read.csv(se_pheno_path, stringsAsFactors = FALSE)

## (Optional) Load target (prediction) strain IDs
# target_ids <- readLines(target_id_path)
target_ids <- NULL   # set to NULL to predict only on training strains

## 2b — Merge genotype + SE phenotype ------------------------------------------
# Aligns strain IDs across geno and pheno; returns numeric 0/1/2 coded genotypes.
merged     <- getMergedData(gt2d = geno_imp, Pheno = pheno_se, testIDs = target_ids)
# merged[[1]] = train data frame   (strain rows × marker + trait columns)
# merged[[2]] = test data frame    (target strains; NULL if target_ids = NULL)

## 2c — Process data (remove rows with NA trait values) ------------------------
proc_data  <- getProcessedData(merged, trait = trait_se)

## 2d — Sample candidate set for CV/GP ----------------------------------------
# noCandidates: how many training lines to include (set equal to nrow to use all)
n_cand     <- nrow(proc_data[[1]])   # use all training individuals
pred_data  <- getPredictionData(proc_data, noCandidates = n_cand)

## 2e — (Optional) Training Set Optimization via Genetic Algorithm -------------
# Uncomment to use STPGA to select an optimal training subset.

# ga_params <- list(
#   errorstat = "MSE", InitPop = NULL, npop = 50, nelite = 5,
#   mutprob = 0.01, mutintensity = 1, niterations = 100,
#   minitbefstop = 10, tabu = FALSE, tabumemsize = 5,
#   plotiters = FALSE, lambda = 0.5, mc.cores = 2
# )
# n_train_select <- 80   # number of lines to select
#
# opt_ts  <- getOptimalTS(pred_data, trait = trait_se, nTraits = 1,
#                         noCandidates = n_cand, nTrainToSelect = n_train_select,
#                         GAParameters = ga_params)
# rand_ts <- getRandomTS(pred_data, trait = trait_se, nTraits = 1,
#                        noCandidates = n_cand, nTrainToSelect = n_train_select)
# ts_cmp  <- getTSComparisons(pred_data, Train_STPGA = opt_ts, Train_Random = rand_ts,
#                              trait = trait_se, nTraits = 1,
#                              testIds = rownames(pred_data[[2]]))
# print(ts_cmp)
# opt_ids <- opt_ts[[2]]   # character vector of selected training strain IDs
opt_ids <- NULL            # set NULL to use all training data

## ── SE:ST (Single-Environment, Single-Trait) ──────────────────────────────────

## 2f — SE:ST Cross-Validation -------------------------------------------------
# k-fold CV using bWGR (emRR, emBB, emBL) — returns prediction ability per model.
cv_st  <- getemCVR(pred_data, trait = trait_se, nTraits = 1, k = k_folds, nIter = n_iter)
print(cv_st)   # named vector: emRR, emBB, emBL prediction abilities

## 2g — SE:ST GEBV Prediction --------------------------------------------------
# Trains on training set (or opt_ids subset) and ranks the prediction set.
gebv_st <- getRankedPredictedValues(
  pred_data, nTraits = 1, trait = trait_se,
  GPModel = "rrBLUP (rrBLUP)",   # options: "rrBLUP (rrBLUP)", "BayesB (bWGR)", etc.
  optTS   = opt_ids               # NULL = use all; opt_ids = GA-selected subset
)
# gebv_st[[1]] = ranked training predictions  (LineID, Predicted Value)
# gebv_st[[2]] = ranked target predictions    (LineID, Predicted Value, Reliability)
writeCVOutTable(as.data.frame(cv_st),   type = "ST", outDirPath = out_dir)
writeGPOutTable(gebv_st[[1]],           type = "ST", outDirPath = out_dir)

## ── SE:MT (Single-Environment, Multi-Trait) ───────────────────────────────────

## 2h — SE:MT Cross-Validation -------------------------------------------------
# Bayesian Multitrait MCMC via BGLR.  A_ext supplies the G matrix as kernel.
cv_mt  <- getMTCVR(
  pred_data, trait = traits_mt, nTraits = length(traits_mt),
  k = k_folds, nIter = n_iter, CVMet = cv_met,
  b_iter = 5000, b_burnin = 1000, b_thin = 10, b_R2 = 0.5, b_digits = 4,
  A_ext = G_tuned   # set to NULL to use A.mat() internally
)
print(summarize_MT_Fits(cv_mt, CVMet = cv_met, traits = traits_mt))

## 2i — SE:MT GEBV Prediction --------------------------------------------------
gebv_mt <- getRankedPredictedValuesMT(
  pred_data, nTraits = length(traits_mt), trait = traits_mt,
  GPModelMT = "BRR (BGLR)",   # options: "BRR (BGLR)", "RKHS (BGLR)", "Spike-Slab (BGLR)"
  optTS     = opt_ids,
  b_iter = 5000, b_burnin = 1000, b_thin = 10, b_R2 = 0.5, b_digits = 4,
  A_ext = G_tuned
)
# gebv_mt[[1]] = ranked training predictions  (LineID, Pred_Trait1, Pred_Trait2, ...)
# gebv_mt[[2]] = ranked target predictions
writeCVOutTable(cv_mt,       type = "MT", outDirPath = out_dir)
writeGPOutTable(gebv_mt[[1]], type = "MT", outDirPath = out_dir)


###############################################################################
## SECTION 3 — ME WORKFLOW: MULTI-ENVIRONMENT
###############################################################################

## 3a — Load ME phenotype & standardise column names ---------------------------
# CSV: one row per observation (strain × environment combination).
pheno_me_raw <- read.csv(me_pheno_path, stringsAsFactors = FALSE)
id_cols_me   <- colnames(pheno_me_raw)[!colnames(pheno_me_raw) %in% trait_me]

# Standardise column names so that downstream ME functions find "Strain", "Env", etc.
std_out      <- stdizePhVarNames(pheno_me_raw, IDColsME  = id_cols_me,
                                  StrainME    = strain_col,
                                  EnvVarME    = env_col,
                                  LocMEColVar = loc_col,
                                  YrMEColVar  = yr_col)
pheno_me_mod <- std_out[[1]]      # standardised phenotype data frame
id_cols_mod  <- std_out[[2]]      # updated ID column names
env_var_mod  <- std_out[[3]]      # standardised environment variable name

# Format ME phenotype into per-trait data frames required by getMergedDataME.
pheno_me_fmt <- getPhenoMEData(PhenoME   = pheno_me_mod,
                               TraitME   = trait_me,
                               nSelTraits = length(trait_me),
                               IDColsME  = id_cols_mod,
                               StrainME  = strain_col)

## 3b — Merge ME genotype + phenotype ------------------------------------------
# Returns DT_Filt_List (pheno per trait) and genoDat_List (geno matrix per trait).
merged_me  <- getMergedDataME(phData   = pheno_me_fmt,
                               genoImp  = geno_imp,
                               TargetIDs = NULL)
DT_list    <- merged_me[[1]]    # list of phenotype data frames, one per trait
geno_list  <- merged_me[[2]]    # list of genotype matrices, one per trait

## 3c — (Optional) Enviromics — build environmental relationship matrix (KE) ---
# Required for FitEnvModels = TRUE in ME CV / GP.  Skip and set KE = NULL otherwise.

# coords_df <- read.csv("<path/to/location_coordinates.csv>")
# # Expected columns: Location, Country, Lat, Long
# env_data  <- getEnvData(coords_df, startDate = "2020-04-01", endDate = "2020-10-31")
# env_ker   <- getEnvKernel(env_data, process = TRUE, Gaussian = FALSE)
# env_ker   <- syncEnvPhenoDat(env_ker, LocCoord = coords_df,
#                               OtherLoc = FALSE)   # TRUE if using alternate location names
# KE        <- env_ker$W                            # symmetric env × env kernel matrix
# plotEnvRel(env_ker, outDirPath = out_dir)         # save heatmap
KE <- NULL   # set to env_ker$W if enviromics data is available

## 3d — ME Cross-Validation ----------------------------------------------------
# Multi-environment CV using BGGE kernels.  KG = G_tuned uses the genomic kinship
# matrix computed in Section 1e; set KG = NULL to estimate internally.
cv_me <- getME_CV(
  DT_1_Filt_List  = DT_list,
  genoDat_List    = geno_list,
  traits          = trait_me,
  KG              = G_tuned,          # genomic kernel; NULL = estimated internally
  KE              = KE,               # env kernel; NULL = no env modeling
  CVMet           = cv_met,
  factVar         = env_var_mod,      # CV stratification factor (e.g., "Env")
  K               = k_folds,
  NIter           = n_iter,
  KMethod         = "Linear",         # "Linear" or "Gaussian"
  FitEnvModels    = !is.null(KE),
  fixedME         = NULL,
  envVar          = env_var_mod,
  IDColsME        = id_cols_mod,
  IDColME         = env_var_mod
)
# Summarise CV results (prediction ability, RMSE, variance components)
cv_me_summary <- summarize_ME_CV_fits(cv_me, cvMet = cv_met)
print(cv_me_summary)
writeCVOutTable(cv_me_summary, type = "ME", outDirPath = out_dir)

## 3e — ME GEBV Prediction -----------------------------------------------------
# Leave-one-factor-out (LOFO) prediction.  maskFactrLev specifies which
# environment(s) to mask; maskProp = 1 masks all genotypes within that level.
me_pred <- fitMEModels_LOFO_Pred(
  DT_1_Filt_List  = DT_list,
  genoDat_List    = geno_list,
  traits          = trait_me,
  KG              = G_tuned,
  KE              = KE,
  factrVar        = env_var_mod,
  maskFactrLev    = unique(pheno_me_mod[[env_var_mod]])[1],  # mask first env as example
  maskProp        = 1,
  nSampleReps     = 1,
  fixedME         = NULL,
  envVar          = env_var_mod,
  IDColsME        = id_cols_mod,
  IDColME         = env_var_mod,
  KMethod         = "Linear",
  FitEnvModels    = !is.null(KE)
)
# me_pred: list of per-trait prediction data frames (LineID, Env, Predicted_Value)

## 3f — Combine and write ME output --------------------------------------------
me_combined <- getCombinedTab(
  outputListME = me_pred,
  TraitME      = trait_me,
  IDColsME     = id_cols_mod,
  IDColME      = env_var_mod,
  envVar       = env_var_mod,
  fitEnvCov    = !is.null(KE),
  reaction     = FALSE
)
writeGPOutTable(me_combined, type = "ME", outDirPath = out_dir)
cat("Pipeline complete. Results written to:", out_dir, "\n")




