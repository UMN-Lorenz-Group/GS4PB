

getME_CV(DT_1_Filt_List,genoDat_List,traits,KG=NULL,KE=NULL,CVMet,factVar="Env",K,NIter,KMethod="Linear",FitEnvModels=FALSE,Reaction=FALSE,dimension_KE=NULL,
  bandwidth = 1,quantil = 0.5,fixedME,envVar,IDColsME,IDColME,LocME,YrME,b_iter = 1000,b_burnin = 200,b_thin = 10,b_tol = 1e-10,b_R2 = 0.5, b_digits = 4){


#### CV schemes 
## CV parameters
k_folds <- 5
n_iter  <- 5
cv_met  <- "CV1" 
							   
G_tuned <- NULL
KE <- NULL
kMet <- "Gaussian"			
fixME <- "Env"
id_col_me <- "StrLocTst" 
locME <- "All"
yrME <- "All"
bys_iter <- 10000 
bys_burnin <- 2000 
bys_thin <-  5
  
  
  DT_1_Filt_List = DT_list
  genoDat_List   = geno_list
  traits         = trait_me
  KG             = G_tuned        # NULL = estimated internally
  KE             = KE
  CVMet          = cv_met
  factVar        = env_var_mod
  K              = k_folds
  NIter          = n_iter
  KMethod        = kMet       # "Linear" or "Gaussian"
  FitEnvModels   = !is.null(KE)
  Reaction		 = FALSE
  dimension_KE	 = NULL
  bandwidth 	 = 1
  quantil 		 = 0.5
  fixedME        = fixME
  envVar         = env_var_mod
  IDColsME       = id_cols_mod
  IDColME        = id_col_me
  LocME 		 = locME
  YrME			 = yrME
  b_iter		 = bys_iter
  b_burnin		 = bys_burnin
  b_thin		 = bys_thin
  b_tol 		 = 1e-10
  b_R2 			 = 0.5 
  b_digits 		 = 4  
  

cv_me <- getME_CV(
  DT_1_Filt_List = DT_list,
  genoDat_List   = geno_list,
  traits         = trait_me,
  KG             = G_tuned,        # NULL = estimated internally
  KE             = KE,
  CVMet          = cv_met,
  factVar        = env_var_mod,
  K              = k_folds,
  NIter          = n_iter,
  KMethod        = kMet,       # "Linear" or "Gaussian"
  FitEnvModels   = !is.null(KE),
  Reaction		 = FALSE,
  dimension_KE	 = NULL,
  bandwidth 	 = 1,
  quantil 		 = 0.5,
  fixedME        = fixME,
  envVar         = env_var_mod,
  IDColsME       = id_cols_mod,
  IDColME        = id_col_me,
  LocME 		 = locME,
  YrME			 = yrME,
  b_iter		 = bys_iter,
  b_burnin		 = bys_burnin,
  b_thin		 = bys_thin,
  b_tol 		 = 1e-10,
  b_R2 			 = 0.5, 
  b_digits 		 = 4
)

 DT, genoDat, strainGeno,
	  KG = NULL, KE = NULL,Fixed,
	  CVMet, k, nIter,FactVar,method,


fitMEModels_CV_fits(
DT=DT_2_ord
genoDat
strainGeno
KG=KG
KE=ke
Fixed=fixedMEMod
CVMet
k=K
nIter=NIter
FactVar=factVarMod
method=KMethod
Reaction= Reaction
dimension_KE
bandwidth = bandwidth
quantil = quantil
trait = NULL
store_full_fits = FALSE
 b_iter = b_iter
b_burnin = b_burnin
b_thin = b_thin
b_tol = b_tol
b_R2 = b_R2
 b_digits = b_digits 
 
 
 Dat_Out <- GS4PB:::getTstIndices_CV_Data(DT_2B, CVMet, nIter, k,cvFactVar=FactVar)
  print(is.null(Dat_Out))
 #)
 
 
 # run all tasks in parallel
 task_results <- foreach(k = seq_len(nrow(idx)),
                        .packages=c("BGGE","BGLR"),
                        .export = c("fit_DTfold","fit_kernel_models","get_kinship_kernels")
                       
  ) %dopar% {
  
  irep <- idx$irep[k]
  fold <- idx$fold[k]

  # isolate data for this task
  DT_fold <- Dat_Out$DT_masked[[irep]][[fold]]
  tst_idx <- Dat_Out$test_idx[[irep]][[fold]]

  DT_fold$isTest <- 0
  DT_fold$isTest[tst_idx] <- 1

  # fit_DTfold <- function(DT_fold,genoDat,trait, KG, KE,fixedCovar,cvFactVar=NULL,Reaction, method,intcpt.rand=FALSE,dimension_KE=NULL,bandwidth = 1,quantil = 0.5,b_iter = 1000,b_burnin = 200,b_thin = 10,b_tol = 1e-10,b_R2 = 0.5, b_digits = 4) {

  # do NOT print from workers; capture if you need it
  # run safely
  res <- try({
    ofit <- GS4PB:::fit_DTfold(
      droplevels(DT_fold), genoDat, trait, KG, KE, Fixed,
      cvFactVar = FactVar, Reaction, method, intcpt.rand=FALSE,
      dimension_KE = dimension_KE, bandwidth = bandwidth, quantil = quantil,
      b_iter = b_iter, b_burnin = b_burnin, b_thin = b_thin,
      b_tol = b_tol, b_R2 = b_R2, b_digits = b_digits
    )
    list(
      status = "ok",
      models = ofit$models,
      raw    = ofit$raw,
      tst_idx = as.integer(tst_idx)
    )
  }, silent = TRUE)

  # normalize output
  if (inherits(res, "try-error")) {
    out <- list(
      status = "error",
      error  = as.character(res),
      models = NULL, raw = NULL, tst_idx = integer(0)
    )
  } else {
    out <- res
  }

  # always include coordinates so we can rebuild structure
  out$irep <- irep
  out$fold <- fold
  out
}