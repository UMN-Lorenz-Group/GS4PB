

 
#' Multi-Environment Genomic Prediction (getME_Pred)
#'
#' Performs  genomic prediction in multi-environment trials,
#' allowing flexible kernel estimation, environmental modeling, and user-defined
#' cross-validation schemes.
#'
#'
#' @param DT_1_Filt_List List of phenotype data frames after merging with genotype data.
#' @param genoDat_List List of genotype data frames after merging with phenotype data.
#' @param traits Character vector/ Character specifying the traits to analyze.
#' @param KG Matrix or NULL. Genotypic relationship matrix. Default = NULL, in which case
#'   the kernel is estimated internally. A user-defined matrix can be supplied.
#' @param KE Matrix or NULL. Environmental relationship matrix. Default = NULL.
#'   Required if `FitEnvModels = TRUE`.
#' @param KMethod Character. Kernel method for estimating genotypic and
#'   genotype × environment kernels. Default = "Linear". Option: "Gaussian".
#' @param FitEnvModels Logical. Whether to include environmental covariates in the model.
#'   Default = FALSE. If TRUE, a user-provided `KE` must be supplied.
#' @param fixedME Character. Variable specifying fixed effects in the model.
#' @param envVar Character. Variable specifying the environmental factor.
#' @param IDColsME Character vector of ID columns (e.g., "Strain", "Location", "Strain-Location").
#' @param LocME Character vector. Locations to include. Default = "All" to use all locations,
#'   or specify one or more locations.
#' @param YrME Character vector. Years to include. Default = "All" to use all years,
#'   or specify one or more years.
#'
#' @return A list containing results from prediction, including predicted values,
#'   model fit statistics, variance components and performance metrics..
#' @examples
#' \dontrun{
#' # Example usage of getMEPred
#' result <- getMEPred(...)
#' }
#' @export
getMEPred <- function(DT_1_Filt_List,genoDat_List,traits,KG=NULL,KE=NULL,KMethod="Linear",FitEnvModels=FALSE,fixedME=fixME,envVar=varEnv,IDColsME,LocME=LocationME,YrME=YearME){

 nTraits <- length(traits)
 GK_Pred_Trt <- list()
 
 for(nTrt in 1:nTraits){
 
  DT_1 <- DT_1_Filt_List[[nTrt]]
  genoDat <- genoDat_List[[nTrt]]
  trait <- traits[nTrt]
 
   if(LocME == "All"){
		DT_1A <- DT_1
  }else if(LocME != "All"){
	    locInd <- which(DT_1$Loc %in% LocME)
		DT_1A <- DT_1[locInd,]
  } 
 
	 if(YrME =="All"){ 
		DT_1B <- DT_1A
	 }else if(YrME != "All"){ 
		yrInd <- which(DT_1A$Year %in% YrME)
		DT_1B <- DT_1A[yrInd,] 
	 }

  
  DT_2 <- DT_1B
  dim(DT_2)
  
  DT_2$Loc <- as.factor(DT_2$Loc) 
  DT_2$Location <- DT_2$Loc
  Loc <- levels(factor(DT_2$env))
  
  nanInd <-   which(is.nan(DT_2[,trait]))
  if(length(nanInd)>0){DT_2 <- DT_2[-nanInd,]}
  
  DT_2 <- droplevels(DT_2)
  
 #EnvVar <- "Loc"
 
  EnvVar <- envVar
  LineID <- "Strain"
  Trait <- trait

### Set Env and Trait Value Columns   

  selColInd <- match(c(EnvVar,LineID,Trait),colnames(DT_2))
  selIDColInd <- match(c(EnvVar,LineID),IDColsME)
  
  colnames(DT_2)[selColInd] <- c("env","gid","value") 
  loc <- levels(factor(DT_2$Loc))
  
  IDColsMEMod <- IDColsME
  IDColsMEMod[selIDColInd] <- c("env","gid")
  
  Fixed= fixedME
  ke <- KE

 IDColsMEList <- list(IDColsME,IDColsMEMod)
 names(IDColsMEList) <- c("IDCols","IDColsMod")
 
 
  if(FitEnvModels==FALSE){
    GK_Pred <- fitMEModels_Predictions(DT_2,genoDat,strainGeno,KG=NULL,KE=NULL,method=KMethod,fitEnvModels=FitEnvModels,FixedTerm=Fixed,IDColsList=IDColsMEList)
  }else if(FitEnvModels==TRUE){ 
    GK_Pred <- fitMEModels_Predictions(DT_2,genoDat,strainGeno,KG=NULL,KE=ke,method=KMethod,fitEnvModels=FitEnvModels,FixedTerm=Fixed,IDColsList=IDColsMEList)
  } 
 
  GK_Pred_Trt[[nTrt]] <- GK_Pred 
   
  }
	return(GK_Pred_Trt)
	

} 







####
#' Standardize ME phenotype column names for downstream pipeline use
#' @param PhenoME Data frame. Raw ME phenotype table.
#' @param IDColsME Character vector. All ID columns in PhenoME.
#' @param StrainME Character. Column name for genotype/strain ID.
#' @param EnvVarME Character. Column name for the environment factor.
#' @param LocMEColVar Character or NULL. Location column name.
#' @param YrMEColVar Character or NULL. Year column name.
#' @return List: [[1]] modified phenotype data frame, [[2]] updated IDColsME, [[3]] standardized env variable name.
#' @export
stdizePhVarNames <- function(PhenoME,IDColsME,StrainME,EnvVarME,LocMEColVar,YrMEColVar){

    print(IDColsME)

	cNamesPhME <- colnames(PhenoME) ## FIX bug1: initialise cNamesPhME before any if-block so it is always defined

	if(!is.null(LocMEColVar) & !is.null(YrMEColVar)){

		locInd <- which(colnames(PhenoME)%in% LocMEColVar)
		PhenoME$LocNew <- PhenoME[,locInd]

		cNamesPhME <- colnames(PhenoME) ## FIX bug1: refresh cNamesPhME after adding LocNew column
		cNamesPhME[which(cNamesPhME %in% YrMEColVar)] <- "Year"
		##
		IDColsME <- c(IDColsME,"LocNew")
		IDColsME[which(IDColsME %in% YrMEColVar)] <- "Year"

		
	}

	## If the selected env var is not "Env" , change any Env into EnvSplName
	if(!is.null(EnvVarME)){
		if(EnvVarME != "Env"){
			preEnvInd <- which(cNamesPhME %in% "Env")
			if(length(preEnvInd) >0) {cNamesPhME[preEnvInd] <- "EnvSplName" }

			idPreEnvInd <- which(IDColsME %in% "Env")
			if(length(idPreEnvInd)>0){ IDColsME[idPreEnvInd] <- "EnvSplName"}
		}
	}

	### Reset EnvVar to a standardized "Env" and strainvarME to "Strain
	cNamesPhME[which(cNamesPhME %in% StrainME)] <- "Strain"
	cNamesPhME[which(cNamesPhME %in% EnvVarME)] <- "Env"

	##
	IDColsME[which(IDColsME %in% StrainME)] <- "Strain"
	IDColsME[which(IDColsME %in% EnvVarME)] <- "Env"

	if(length(grep("LocNew",cNamesPhME))>0){ ## FIX bug2: was cNamePhME (typo), now cNamesPhME
		cNamesPhME <- gsub("LocNew","Loc",cNamesPhME)
		IDColsME <- gsub("LocNew","Loc",IDColsME)
	}
	colnames(PhenoME) <- cNamesPhME ## FIX bug2: moved outside if-block so colnames are always written back

	EnvVarMEMod <- "Env" 
	outList <- list(PhenoME,IDColsME,EnvVarMEMod)
    print(IDColsME)

	return(outList)
}

####



### PC Analysis helper ###
# genoDF: rows = markers; cols 1:5 = SNPID/Chrom/Position/REF/ALT; cols 6+ = sample IDs
# nPCs: number of principal components to compute (default 2)
# Returns: list(pcDF = data.frame with Strain + PC columns, d = singular values)

getGenoPC <- function(genoDF, nPCs = 2){

	# Extract genotype matrix (rows=markers, cols=samples)
	genoMat <- genoDF[, -c(1:5)]

	# Convert to numeric
	genoNum <- apply(genoMat, 2, as.numeric)  # rows=markers, cols=samples

	# Transpose so rows=samples, cols=markers
	ZZ <- t(genoNum)
	strainIDs <- colnames(genoMat)

	# Center markers; replace any residual NAs with 0
	ZZ_centered <- scale(ZZ, center = TRUE, scale = FALSE)
	ZZ_centered[is.na(ZZ_centered)] <- 0

	# Cap nPCs to matrix dimensions
	nPCs <- min(nPCs, nrow(ZZ_centered) - 1, ncol(ZZ_centered) - 1)

	# SVD via RSpectra
	pc <- RSpectra::svds(ZZ_centered, k = nPCs)

	# Build PC data frame
	pcDF <- as.data.frame(pc$u)
	colnames(pcDF) <- paste0("PC", seq_len(nPCs))
	pcDF$Strain <- strainIDs
	pcDF <- pcDF[, c("Strain", paste0("PC", seq_len(nPCs)))]

	return(list(pcDF = pcDF, d = pc$d))
}




stdizeGPCovarNames <- function(IDColsME,StrainME,EnvVarME,varEnv_MEGP,fixMEGP,mskFactrGP){

	IDColsME[which(IDColsME %in% StrainME)] <- "Strain"

	if(!is.null(EnvVarME)){ ## FIX bug8b: wrap both env conditionals in null-guard

		if(varEnv_MEGP == EnvVarME){
			IDColsME[which(IDColsME %in% varEnv_MEGP)] <- "Env"
			varEnv_MEGP_Mod <- "Env"
		}else{varEnv_MEGP_Mod <- varEnv_MEGP}

		if(fixMEGP == EnvVarME){
			IDColsME[which(IDColsME %in% EnvVarME)] <- "Env" ## FIX bug7b: use EnvVarME as search key (fixMEGP==EnvVarME here, and prior block may have renamed it)
			fixMEGP_Mod <- "Env"
		}else{fixMEGP_Mod <- fixMEGP}

		if(mskFactrGP == EnvVarME){
			IDColsME[which(IDColsME %in% EnvVarME)] <- "Env" ## FIX bug7b: use EnvVarME as search key (mskFactrGP==EnvVarME here, and prior block may have renamed it)
			mskFactrGP_Mod <- "Env"
		}else{mskFactrGP_Mod <- mskFactrGP}

	} ## end if(!is.null(EnvVarME)) FIX bug8b

	return(list(varEnv_MEGP_Mod,fixMEGP_Mod,mskFactrGP_Mod,IDColsME)) ## FIX bug6b: return IDColsME as third element

}
		 



#' Leave-One-Factor-Out Prediction for Multi-Environment Models
#'
#' Fits ME genomic prediction models using a leave-one-factor-level-out CV
#' scheme. Supports MM, MDs, MDe (multi-env) and SM (single-env) BGGE models.
#'
#' @param DT_1_Filt_List List of filtered data tables, one per trait.
#' @param genoDat_List List of genotype matrices, one per trait.
#' @param traits Character vector of trait names.
#' @param KG Genomic kinship matrix (optional).
#' @param KE Environmental kernel matrix (optional).
#' @param factrVar Factor variable name for leave-one-factor-out CV.
#' @param maskFactrLev Factor level to mask (leave out).
#' @param maskProp Proportion of factor level to mask (default 1).
#' @param nSampleReps Number of sampling replicates (default 1).
#' @param fixedME Fixed effect variable names.
#' @param envVar Environment variable name.
#' @param IDColsME Vector of ID column names.
#' @param IDColME Unique ID column name.
#' @param LocME Location filter ("All" or specific location).
#' @param YrME Year filter ("All" or specific year).
#' @param KMethod Kernel method for building kinship kernels.
#' @param FitEnvModels Logical; include environmental kernel (default FALSE).
#' @param Reaction Logical; use reaction norm model (default FALSE).
#' @param dimension_KE Dimension for KE (optional).
#' @param bandwidth Bandwidth parameter (default 1).
#' @param quantil Quantile parameter (default 0.5).
#' @param b_iter BGGE/BGLR iterations (default 1000).
#' @param b_burnin BGGE/BGLR burn-in (default 200).
#' @param b_thin BGGE/BGLR thinning (default 10).
#' @param b_tol BGGE/BGLR tolerance (default 1e-10).
#' @param b_R2 BGGE/BGLR R2 prior (default 0.5).
#' @param b_digits Rounding digits for Maturity.rating (default 4).
#'
#' @return A list with elements \code{Cor}, \code{vcomp}, and \code{preds},
#'   each a list of length \code{nTraits}.
#'
#' @export
fitMEModels_LOFO_Pred <- function(DT_1_Filt_List,genoDat_List,traits,KG = NULL,KE = NULL,factrVar, maskFactrLev, maskProp=1,nSampleReps=1,fixedME,envVar,IDColsME,IDColME,LocME,YrME,KMethod,FitEnvModels=FALSE,Reaction = FALSE,dimension_KE=NULL,bandwidth = 1,quantil = 0.5,b_iter = 1000,b_burnin = 200,b_thin = 10,b_tol = 1e-10,b_R2 = 0.5, b_digits = 4){


 print("ME_GP_In")
 nTraits <- length(traits)
 UniqID <- IDColME
 
 ME_Out_CV_Trt <- list()
 DT_2_List <- list() 
 IDColsMEList_Out <- list() 
  
 if(FitEnvModels){
    ke <- KE
 }else{
    ke <- NULL 
 }
  
 cor_LOT_CV_List_Trt <- list()
 var_LOT_CV_List_Trt <- list()
 fit_Out_LOT_CV_List_Trt <- list()
  
 
 for(nTrt in 1:nTraits){
 
  DT_1 <- DT_1_Filt_List[[nTrt]]
  genoDat <- genoDat_List[[nTrt]]
  trait <- traits[nTrt]
 
  if(LocME == "All"){
		DT_1A <- DT_1
  }else if(LocME != "All"){
	    locInd <- which(DT_1$Loc %in% LocME)
		DT_1A <- DT_1[locInd,]
  } 
 
  if(YrME =="All"){ 
		DT_1B <- DT_1A
  }else if(YrME != "All"){ 
		yrInd <- which(DT_1A$Year %in% YrME)
		DT_1B <- DT_1A[yrInd,] 
  }

  
  nanInd <-   which(is.nan(DT_1B[,trait]))
  if(length(nanInd)>0){DT_1C <- DT_1B[-nanInd,]}else{DT_1C <- DT_1B}
  
  DT_2 <- DT_1C 
  dim(DT_2)
 
  # DT_2$Loc <- as.factor(DT_2$Loc) 
  # DT_2$Location <- DT_2$Loc
  # Loc <- levels(factor(DT_2$Loc))
  
  DT_2 <- droplevels(DT_2)
  EnvVar <- envVar
  LineID <- "Strain"
  Trait <- trait
 
## 
  selColInd <- match(c(EnvVar,LineID,Trait,UniqID),colnames(DT_2))
  selIDColInd <- match(c(EnvVar,LineID,UniqID),IDColsME)
  
  colnames(DT_2)[selColInd] <- c("env","gid","value","uniqID") 
  loc <- levels(factor(DT_2$env))
  
  IDColsMEMod <- IDColsME
  IDColsMEMod[selIDColInd] <- c("env","gid","uniqID")
  
##  
## Set Env and Trait Value Columns   
## Modify all modified vars
  
  factInd <- which(IDColsME %in% factrVar)
  factVarMod <- factrVar
  factVarMod <- IDColsMEMod[factInd]

##
  fixInd <- which(IDColsME %in% fixedME)
  fixedMEMod <-  IDColsMEMod[fixInd]

  strainGeno <- rownames(genoDat)

##
  IDColsMEList <- list(IDColsME,IDColsMEMod)
  names(IDColsMEList) <- c("IDCols","IDColsMod")
  
### 

  fact_col <- which(colnames(DT_2) %in% factVarMod) 
  id_col = "gid"
  env_col = "env"
  y_col = "value"
  	 
  if(length(fact_col) ==0){
	stop("Factor for leave-one-factor-level-out CV should be present in data table.")
  }

  DT_2$Obs <- DT_2$value
  N <- nrow(DT_2)
  unq_env <- sort(unique(DT_2[[env_col]]))
  n_env   <- length(unq_env)
	  
  unq_factrs <- levels(factor(DT_2[[fact_col]]))
  n_factrs   <- length(unq_factrs)
  
###
# --- Prepare sampling replicates  & masks ---
  
  if(maskProp ==1){
	   nSampleReps <- 1
  }
 
#### 

 
  DT_masked_rep <- vector("list", nSampleReps)
  train_idx_rep <- vector("list", nSampleReps)
  test_idx_rep  <- vector("list", nSampleReps)

  base_seed <- 125

  ord <- order(DT_2$env,DT_2$gid)
  DT_2_ord <- DT_2[ord,,drop=FALSE]
  
  
  factr_mask <- DT_2_ord[[fact_col]] == maskFactrLev
  factrIndices <- which(factr_mask) # same rows in 2B; we only NA-out the value col 
  DT_2_Factr <- DT_2_ord[factrIndices,]

  k <- round(1/maskProp, digits=0)

 # Within a factor level perform k-fold CV-1 like split, where  
 # CV1: leave-genotypes-out across ALL environments
 # Test = all rows of the held-out genotype set
 # ============================================================
  
  unq_ids <- unique(DT_2_Factr[[id_col]])
  n_ids   <- length(unq_ids) 
	
  if (k > n_ids) stop("'k' cannot exceed the number of unique genotypes. Adjust mask proportion. 
                       To decrease 'k' increase mask proportion and vice versa")

  .split_kfold <- function(n, k) {
    perm <- sample.int(n)
    base <- n %/% k; r <- n %% k
    sizes <- rep(base, k); if (r > 0) sizes[seq_len(r)] <- sizes[seq_len(r)] + 1L
    split(perm, rep(seq_len(k), sizes))
  }

  for (rep in seq_len(nSampleReps)){
     
	  set.seed(base_seed + rep)

      # disjoint folds over genotypes
      idx_in_ids <- .split_kfold(n_ids, k)
      
      # per-fold results
      DT_list   <- vector("list", k)
      tr_list   <- vector("list", k)
      tst_list  <- vector("list", k)

      for(f in seq_len(k)){

	    DT_2_ord$isTest <- 0
	   
        test_ids_idx <- sort(idx_in_ids[[f]])
        test_ids     <- unq_ids[test_ids_idx]

        DT_2_Factr <- DT_2_ord[DT_2_ord[[fact_col]]== maskFactrLev,]
		DT_2_NonFactr <- DT_2_ord[DT_2_ord[[fact_col]]!= maskFactrLev,]

###
        test_rows_factr <- which(DT_2_Factr[[id_col]] %in% test_ids)
        
        DT_2_Factr$isTest[test_rows_factr] <- 1	
		DT_2_Factr[[y_col]][test_rows_factr] <- NA 
		
		DT_f <- rbind.data.frame(DT_2_Factr,DT_2_NonFactr)
			
		ord <- order(DT_f$env, DT_f$gid)
        DT_f_ord <- DT_f[ord, , drop = FALSE]
		
		N <- nrow(DT_f_ord) 
		test_rows <- which(DT_f_ord$isTest == 1)
		train_rows <- setdiff(seq_len(N), test_rows)

        DT_list[[f]]  <- DT_f_ord
        tr_list[[f]]  <- train_rows
        tst_list[[f]] <- test_rows
      }

      DT_masked_rep[[rep]] <- DT_list
      train_idx_rep[[rep]] <- tr_list
      test_idx_rep [[rep]] <- tst_list
    }

    MaskedData <- list(DT_masked = DT_masked_rep,
                       train_idx = train_idx_rep,
					   test_idx  = test_idx_rep,
                       factor_levels = unq_factrs)
 
   # after building MaskedData {...}

	uniqid <- "uniqID"
	addColName <- factVarMod  # single column name string (character scalar)

	fit_MM_df_Out  <- NULL
	fit_MDs_df_Out <- NULL
	fit_MDe_df_Out <- NULL
	fit_SM_df_Out  <- NULL

	# VarComp collectors — accumulate per fold, average after loops
	varMM_folds  <- list(); varMDs_folds <- list(); varMDe_folds <- list()
	varSM_folds  <- list()
	corSM_LOT_CV <- NULL

	for (sampRepNo in seq_len(nSampleReps)) {

	  DT_lists  <- MaskedData$DT_masked[[sampRepNo]]

	  for (f in seq_along(DT_lists)) {

		DT_f <- DT_lists[[f]]

		# Response and keys
		y <- "value"; gid <- "gid"; env <- "env"

		# Align X to DT_f$gid
		idx <- match(DT_f[[gid]], rownames(genoDat))
		if (anyNA(idx)) {
		  missing_g <- unique(DT_f[[gid]][is.na(idx)])
		  stop("These gids are missing in genoDat: ",
			   paste(head(missing_g, 10), collapse = ", "),
			   if (length(missing_g) > 10) " ...")
		}
		X <- genoDat[idx, , drop = FALSE]

		# Optional kernels method if KG missing
		KMethod2 <- if (!is.null(KG)) NULL else KMethod

		ne <- length(unique(DT_f[[env]]))

		# Build fixed-effects matrix
		fixedMEMod_terms <- paste(fixedMEMod, collapse = " + ")
		fixed <- model.matrix(as.formula(paste("~ 0 +", fixedMEMod_terms)), data = DT_f)

		if (ne > 1) {
		  MM  <- get_kinship_kernels(DT_f[, c(env, gid, y)], X,
				   K_G = KG, K_E = ke, intercept.random = FALSE,
				   model = "MM",  method = KMethod2,
				   dimension_KE = dimension_KE, bandwidth = bandwidth,
				   quantil = quantil, reaction = Reaction)

		  MDs <- get_kinship_kernels(DT_f[, c(env, gid, y)], X,
				   K_G = KG, K_E = ke, intercept.random = FALSE,
				   model = "MDs", method = KMethod2,
				   dimension_KE = dimension_KE, bandwidth = bandwidth,
				   quantil = quantil, reaction = Reaction)

		  MDe <- get_kinship_kernels(DT_f[, c(env, gid, y)], X,
				   K_G = KG, K_E = ke, intercept.random = FALSE,
				   model = "MDe", method = KMethod2,
				   dimension_KE = dimension_KE, bandwidth = bandwidth,
				   quantil = quantil, reaction = Reaction)

		  fit_MM  <- fit_kernel_models(y = y, env = env, gid = gid, uniqID = uniqid,
						Data = DT_f, random = MM,  fixed = fixed, engine = "BGGE",
						addCol = addColName, verbose = TRUE,
						iterations = b_iter, burnin = b_burnin,
						thining = b_thin, tol = b_tol, R2 = b_R2, digits = b_digits)

		  fit_MDs <- fit_kernel_models(y = y, env = env, gid = gid, uniqID = uniqid,
						Data = DT_f, random = MDs, fixed = fixed, engine = "BGGE",
						addCol = addColName, verbose = TRUE,
						iterations = b_iter, burnin = b_burnin,
						thining = b_thin, tol = b_tol, R2 = b_R2, digits = b_digits)

		  fit_MDe <- fit_kernel_models(y = y, env = env, gid = gid, uniqID = uniqid,
						Data = DT_f, random = MDe, fixed = fixed, engine = "BGGE",
						addCol = addColName, verbose = TRUE,
						iterations = b_iter, burnin = b_burnin,
						thining = b_thin, tol = b_tol, R2 = b_R2, digits = b_digits)

		  if (identical(trait, "Maturity.rating")) {
			fit_MM$yHat  <- round(fit_MM$yHat,  0)
			fit_MDs$yHat <- round(fit_MDs$yHat, 0)
			fit_MDe$yHat <- round(fit_MDe$yHat, 0)
		  }

		  DT_Out  <- fit_MM$dat
		  Test_Idx <- which(DT_Out$isTest == 1)

		  obs_factr      <- DT_Out$Obs[Test_Idx]
		  pred_MM_factr  <- fit_MM$yHat [Test_Idx]
		  pred_MDs_factr <- fit_MDs$yHat[Test_Idx]
		  pred_MDe_factr <- fit_MDe$yHat[Test_Idx]

		  varMM_folds[[length(varMM_folds)+1]]   <- fit_MM$VarComp
		  varMDs_folds[[length(varMDs_folds)+1]] <- fit_MDs$VarComp
		  varMDe_folds[[length(varMDe_folds)+1]] <- fit_MDe$VarComp

		  mm_df <- data.frame(
			UniqID = DT_Out[[uniqid]][Test_Idx],
			Strain = DT_Out[[gid]][Test_Idx],
			Env    = DT_Out[[env]][Test_Idx],
			Rep    = sampRepNo,
			Obs    = obs_factr,
			Pred   = pred_MM_factr
		  )
		  mm_df[[factrVar]] <- DT_Out[[addColName]][Test_Idx]

		  mds_df <- data.frame(
			UniqID = DT_Out[[uniqid]][Test_Idx],
			Strain = DT_Out[[gid]][Test_Idx],
			Env    = DT_Out[[env]][Test_Idx],
			Rep    = sampRepNo,
			Obs    = obs_factr,
			Pred   = pred_MDs_factr
		  )
		  mds_df[[factrVar]] <- DT_Out[[addColName]][Test_Idx]

		  mde_df <- data.frame(
			UniqID = DT_Out[[uniqid]][Test_Idx],
			Strain = DT_Out[[gid]][Test_Idx],
			Env    = DT_Out[[env]][Test_Idx],
			Rep    = sampRepNo,
			Obs    = obs_factr,
			Pred   = pred_MDe_factr
		  )
		  mde_df[[factrVar]] <- DT_Out[[addColName]][Test_Idx]

		  fit_MM_df_Out  <- rbind(fit_MM_df_Out,  mm_df)
		  fit_MDs_df_Out <- rbind(fit_MDs_df_Out, mds_df)
		  fit_MDe_df_Out <- rbind(fit_MDe_df_Out, mde_df)

		} else {
		  # Single-env
		  SM <- get_kinship_kernels(DT_f[, c(env, gid, y)], X,
				   K_G = KG, K_E = NULL, intercept.random = FALSE,
				   model = "SM", method = KMethod2,
				   dimension_KE = dimension_KE, bandwidth = bandwidth,
				   quantil = quantil, reaction = Reaction)

		  fixed1 <- model.matrix(~ 1, data = DT_f)

		  fit_SM <- fit_kernel_models(y = y, env = env, gid = gid, uniqID = uniqid,
					   Data = DT_f, random = SM, fixed = fixed1, engine = "BGGE",
					   addCol = addColName, verbose = TRUE,
					   iterations = b_iter, burnin = b_burnin,
					   thining = b_thin, tol = b_tol, R2 = b_R2, digits = b_digits)

		  if (identical(trait, "Maturity.rating")) fit_SM$yHat <- round(fit_SM$yHat, 0)

		  DT_Out  <- fit_SM$dat
		  # test rows already marked via masking (isTest == 1)
		  Test_Idx <- which(DT_Out$isTest == 1)

		  obs_factr   <- DT_Out$Obs[Test_Idx]
		  pred_SM     <- fit_SM$yHat[Test_Idx]

		  varSM_folds[[length(varSM_folds)+1]] <- fit_SM$VarComp

		  sm_df <- data.frame(
			UniqID = DT_Out[[uniqid]][Test_Idx],
			Strain = DT_Out[[gid]][Test_Idx],
			Env    = DT_Out[[env]][Test_Idx],
			Rep    = sampRepNo,
			Obs    = obs_factr,
			Pred   = pred_SM
		  )
		  sm_df[[factrVar]] <- DT_Out[[addColName]][Test_Idx]

		  fit_SM_df_Out <- rbind(fit_SM_df_Out, sm_df)
		}
	  } # folds
	} # reps

	# Average VarComp across folds (each element is a data frame with numeric cols)
	.avg_vcomp <- function(lst) {
	  if (length(lst) == 1L) return(lst[[1]])
	  num_cols <- sapply(lst[[1]], is.numeric)
	  out <- lst[[1]]
	  for (i in seq_along(lst)[-1]) out[, num_cols] <- out[, num_cols] + lst[[i]][, num_cols]
	  out[, num_cols] <- out[, num_cols] / length(lst)
	  out
	}

	# Aggregate correlations and pack return
	if (!is.null(fit_MM_df_Out)) {
	  corMM_LOT_CV  <- cor(fit_MM_df_Out$Pred,  fit_MM_df_Out$Obs,  use = "pairwise.complete.obs")
	  corMDs_LOT_CV <- cor(fit_MDs_df_Out$Pred, fit_MDs_df_Out$Obs, use = "pairwise.complete.obs")
	  corMDe_LOT_CV <- cor(fit_MDe_df_Out$Pred, fit_MDe_df_Out$Obs, use = "pairwise.complete.obs")

	  cor_LOT_CV_List  <- list(MM  = corMM_LOT_CV, MDs = corMDs_LOT_CV, MDe = corMDe_LOT_CV)
	  var_LOT_CV_List  <- list(MM  = .avg_vcomp(varMM_folds),
	                           MDs = .avg_vcomp(varMDs_folds),
	                           MDe = .avg_vcomp(varMDe_folds))
	  fit_Out_LOT_CV_List <- list(MM = fit_MM_df_Out, MDs = fit_MDs_df_Out, MDe = fit_MDe_df_Out)
	} else {
	  corSM_LOT_CV <- cor(fit_SM_df_Out$Pred, fit_SM_df_Out$Obs, use = "pairwise.complete.obs")
	  cor_LOT_CV_List  <- list(SM = corSM_LOT_CV)
	  var_LOT_CV_List  <- list(SM = .avg_vcomp(varSM_folds))
	  fit_Out_LOT_CV_List <- list(SM = fit_SM_df_Out)
	}

   
   cor_LOT_CV_List_Trt[[nTrt]] <- cor_LOT_CV_List
   var_LOT_CV_List_Trt[[nTrt]] <- var_LOT_CV_List
   fit_Out_LOT_CV_List_Trt[[nTrt]] <- fit_Out_LOT_CV_List
  
  }
    
  # Named return for clarity
  return(list(
    Cor = cor_LOT_CV_List_Trt,
    vcomp = var_LOT_CV_List_Trt,
    preds = fit_Out_LOT_CV_List_Trt
  ))
}
	


#### 
#' Compute Kernel Matrices for Genomic Models
#'
#' The `getK` function computes kernel matrices for genomic models, supporting
#' different types of kernels and models. It is designed for applications in
#' genomic selection and prediction.
#'
#' @param Y A data frame containing the phenotypic data. The first two columns 
#'   should represent the environment and genotype identifiers.
#' @param X A numeric matrix of marker data, with genotypes as rows and markers 
#'   as columns.
#' @param kernel A character string specifying the kernel type. Options are 
#'   `"GK"` (Gaussian kernel) or `"GB"` (Genomic relationship matrix). Default 
#'   is `"GK"`.
#' @param setKernel A list of precomputed kernel matrices. If provided, this 
#'   overrides the `kernel` parameter. Default is `NULL`.
#' @param bandwidth A numeric vector specifying the bandwidth(s) for the 
#'   Gaussian kernel. Default is `1`.
#' @param model A character string specifying the model type. Options are:
#'   - `"SM"`: Single model
#'   - `"MM"`: Multiple model
#'   - `"MDs"`: Model with genotype-by-environment interaction (shared)
#'   - `"MDe"`: Model with genotype-by-environment interaction (exclusive)
#'   Default is `"SM"`.
#' @param quantil A numeric value representing the quantile used for scaling in 
#'   the Gaussian kernel. Default is `0.5`.
#' @param intercept.random A logical value indicating whether to include a 
#'   random intercept in the model. Default is `FALSE`.
#'
#' @return A list of kernel matrices. Each element of the list contains:
#'   - `Kernel`: The kernel matrix.
#'   - `Type`: A character string indicating the type of kernel, either `"D"` 
#'     (diagonal) or `"BD"` (block-diagonal for interaction models).
#'
#' @details
#' - For the `"GK"` kernel, the function calculates a Gaussian kernel matrix 
#'   based on pairwise distances between genotypes.
#' - For the `"GB"` kernel, it computes a genomic relationship matrix using 
#'   the marker data.
#' - The `setKernel` parameter allows the user to provide custom kernel matrices 
#'   for advanced use cases.
#'
#' @examples
#' \dontrun{
#' # Example phenotypic data
#' Y <- data.frame(
#'   Environment = rep(c("Env1", "Env2"), each = 5),
#'   Genotype = factor(rep(letters[1:5], 2)),
#'   Value = rnorm(10)
#' )
#'
#' # Example marker data
#' X <- matrix(runif(100), nrow = 5)
#' rownames(X) <- letters[1:5]
#'
#' # Compute Gaussian kernels
#' kernels <- getK(Y = Y, X = X, kernel = "GK", model = "SM")
#'
#' # Compute kernels with genotype-by-environment interaction
#' kernels_interaction <- getK(Y = Y, X = X, kernel = "GK", model = "MDs")
#' }
#'
#' @export

getK <- 
function (Y, X, kernel = c("GK", "GB"), setKernel = NULL, bandwidth = 1, 
          model = c("SM", "MM", "MDs", "MDe"), quantil = 0.5, intercept.random = FALSE) 
{
  Y <- as.data.frame(Y)
  Y[colnames(Y)[1:2]] <- lapply(Y[colnames(Y)[1:2]], factor)
  subjects <- levels(Y[, 2])
  env <- levels(Y[, 1])
  nEnv <- length(env)
  if (any(table(Y[, c(1:2)]) > 1)) 
    warning("There are repeated genotypes in some environment. They were kept")
    switch(model, SM = {
     if (nEnv > 1) stop("Single model choosen, but more than one environment is in the phenotype file")
      Zg <- model.matrix(~factor(Y[, 2L]) - 1)
     }, Cov = {
      Zg <- model.matrix(~factor(subjects) - 1)
     }, {
      Ze <- model.matrix(~factor(Y[, 1L]) - 1)
      Zg <- model.matrix(~factor(Y[, 2L]) - 1)
    })
  if (is.null(setKernel)){
    if (is.null(rownames(X))) 
      stop("Genotype names are missing")
    if (!all(subjects %in% rownames(X))) 
      stop("Not all genotypes presents in the phenotypic file are in marker matrix")
    X <- X[subjects, ]
    switch(kernel, GB = {
      ker.tmp <- tcrossprod(X)/ncol(X)
      G <- list(list(Kernel = Zg %*% tcrossprod(ker.tmp, 
                                                Zg), Type = "D"))
    }, GK = {
      D <- (as.matrix(dist(X)))^2
      G <- list()
      for (i in 1:length(bandwidth)) {
        ker.tmp <- exp(-bandwidth[i] * D/quantile(D, 
                                                  quantil))
        G[[i]] <- list(Kernel = Zg %*% tcrossprod(ker.tmp, 
                                                  Zg), Type = "D")
      }
    }, {
      stop("kernel selected is not available. Please choose one method available or make available other kernel through argument K")
    })
    names(G) <- seq(length(G))
  }
  else {
    nullNames <- sapply(setKernel, function(x) any(sapply(dimnames(x), 
                                                          is.null)))
    if (any(nullNames)) 
      stop("Genotype names are missing in some kernel")
    equalNames <- sapply(setKernel, function(x) mapply(function(z, 
                                                                y) all(z %in% y), z = list(subjects), y = dimnames(x)))
    if (!all(equalNames)) 
      stop("Not all genotypes presents in phenotypic file are in the kernel matrix.\n             Please check dimnames")
    K <- lapply(setKernel, function(x) x[subjects, subjects])
    ker.tmp <- K
    G <- lapply(ker.tmp, function(x) list(Kernel = Zg %*% 
                                            tcrossprod(x, Zg), Type = "D"))
    if (is.null(names(K))) {
      names(G) <- seq(length(G))
    }
    else {
      names(G) <- names(setKernel)
    }
  }
  # tmp.names <- names(G)
  # names(G) <- if (length(G) > 1) 
  #   paste("G", tmp.names, sep = "_")
  # else "G"
  # 
  if (length(G) == 0) {
    names(G) <- character(0)
  } else {
    tmp.names <- names(G)
    if (is.null(tmp.names)) tmp.names <- seq_along(G)
    names(G) <- if (length(G) > 1) 
      paste("G", tmp.names, sep = "_")
    else "G"
  }
  switch(model, SM = {
    out <- G
  }, MM = {
    out <- G
  }, MDs = {
    E <- tcrossprod(Ze)
    GE <- lapply(G, function(x) list(Kernel = x$Kernel * 
                                       E, Type = "BD"))
    names(GE) <- if (length(G) > 1) paste("GE", tmp.names, 
                                          sep = "_") else "GE"
    out <- c(G, GE)
  }, MDe = {
    ZEE <- matrix(data = 0, nrow = nrow(Ze), ncol = ncol(Ze))
    out.tmp <- list()
    for (j in 1:length(G)) {
      out.tmp <- c(out.tmp, lapply(1:nEnv, function(i) {
        ZEE[, i] <- Ze[, i]
        ZEEZ <- ZEE %*% t(Ze)
        K3 <- list(Kernel = G[[j]]$Kernel * ZEEZ, Type = "BD")
        return(K3)
      }))
    }
    if (length(G) > 1) {
      names(out.tmp) <- paste(rep(env, length(G)), rep(tmp.names, 
                                                       each = nEnv), sep = "_")
    } else {
      names(out.tmp) <- env
    }
    out <- c(G, out.tmp)
  }, {
    stop("Model selected is not available ")
  })
  if (intercept.random) {
    Gi <- list(Kernel = Zg %*% tcrossprod(diag(length(subjects)), 
                                          Zg), Type = "D")
    out <- c(out, list(Gi = Gi))
  }
  return(out)
}






#######

#' Multi-Environment Cross-Validation (getME_CV)
#'
#' Performs cross-validation for genomic prediction in multi-environment trials,
#' allowing flexible kernel estimation, environmental modeling, and user-defined
#' cross-validation schemes.
#'
#' @param DT_1_Filt_List List of phenotype data frames after merging with genotype data.
#' @param genoDat_List List of genotype data frames after merging with phenotype data.
#' @param traits Character vector/ Character specifying the traits to analyze.
#' @param KG Matrix or NULL. Genotypic relationship matrix. Default = NULL, in which case
#'   the kernel is estimated internally. A user-defined matrix can be supplied.
#' @param KE Matrix or NULL. Environmental relationship matrix. Default = NULL.
#'   Required if `FitEnvModels = TRUE`.
#' @param CVMet Character. Cross-validation scheme. Default = "CV1".
#'   Options: "CV1", "CV2", "CV0", "CV00","CV-LOFO".
#' @param factVar Character. Factor variable to be used in the model for leave-one-factor-out CV.
#' @param K Integer. Number of folds in k-fold cross-validation.
#' @param NIter Integer. Number of replicates for cross-validation.
#' @param KMethod Character. Kernel method for estimating genotypic and
#'   genotype × environment kernels. Default = "Linear". Option: "Gaussian".
#' @param FitEnvModels Logical. Whether to include environmental covariates in the model.
#'   Default = FALSE. If TRUE, a user-provided `KE` must be supplied.
#' @param fixedME Character. Variable specifying fixed effects in the model.
#' @param envVar Character. Variable specifying the environmental factor.
#' @param IDColsME Character vector of ID columns (e.g., "Strain", "Location", "Strain-Location").
#' @param IDColME Character. Unique identifier column.
#' @param LocME Character vector. Locations to include. Default = "All" to use all locations,
#'   or specify one or more locations.
#' @param YrME Character vector. Years to include. Default = "All" to use all years,
#'   or specify one or more years.
#'
#' @return A list containing cross-validation results, including predicted values,
#'   model fit statistics, and performance metrics.
#'
#' @export
#' @importFrom foreach foreach
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' result <- getME_CV(
#'   DT_1_Filt_List, genoDat_List, traits = c("Yield"),
#'   CVMet = "CV1", K = 5, NIter = 10
#' )
#' }


getME_CV <- function(DT_1_Filt_List,genoDat_List,traits,KG=NULL,KE=NULL,CVMet,factVar="Env",K,NIter,KMethod="Linear",FitEnvModels=FALSE,Reaction=FALSE,dimension_KE=NULL,
  bandwidth = 1,quantil = 0.5,fixedME,envVar,IDColsME,IDColME,LocME,YrME,b_iter = 1000,b_burnin = 200,b_thin = 10,b_tol = 1e-10,b_R2 = 0.5, b_digits = 4){

 print("ME_CV_In")
 nTraits <- length(traits)
 UniqID <- IDColME
 EnvVar <- envVar
 LineID <- "Strain"
  
 ### 
 
 ME_Out_CV_Trt <- list()
 DT_2_List <- list() 
 IDColsMEList_Out <- list() 
 
 
 if(FitEnvModels){
    ke <- KE
 }else{
    ke <- NULL 
 }
 
 for(nTrt in 1:nTraits){
 
  DT_1 <- DT_1_Filt_List[[nTrt]]
  genoDat <- genoDat_List[[nTrt]]
  trait <- traits[nTrt]
 
  if(LocME == "All"){
		DT_1A <- DT_1
  }else if(LocME != "All"){
	    locInd <- which(DT_1$Loc %in% LocME)
		DT_1A <- DT_1[locInd,]
  } 
 
  if(YrME =="All"){ 
		DT_1B <- DT_1A
  }else if(YrME != "All"){ 
		yrInd <- which(DT_1A$Year %in% YrME)
		DT_1B <- DT_1A[yrInd,] 
  }

#### 
  
  DT_2 <- DT_1B
  #dim(DT_2)
  
  DT_2$Loc <- as.factor(DT_2$Loc) 
  Loc <- levels(factor(DT_2$Loc))
  
  #nanInd <-   which(is.nan(DT_2[,trait]))
  nanInd <- which(is.nan(suppressWarnings(as.numeric(DT_2[[trait]])))) 
  if(length(nanInd)>0){DT_2 <- DT_2[-nanInd,]}
  
  DT_2 <- droplevels(DT_2)
  Trait <- trait
  
#### 
### Set Env and Trait Value Columns   

  selColInd <- match(c(EnvVar,LineID,Trait,UniqID),colnames(DT_2))
  selIDColInd0 <- match(c(EnvVar,LineID,UniqID),IDColsME)
  if(anyNA(selIDColInd0)) selIDColInd <- selIDColInd0[!is.na(selIDColInd0)] else selIDColInd <- selIDColInd0
  
  
  colnames(DT_2)[selColInd] <- c("env","gid","value","uniqID") 
  print(selIDColInd0)
  
  
  IDColsMEMod <- IDColsME
  IDColsMEMod[selIDColInd] <- c("env","gid","uniqID")
  print(IDColsMEMod)
  print(IDColsMEMod[selIDColInd])
  
  
## Modify all modified vars
  
  factInd <- which(IDColsME %in% factVar)
  factVarMod <- factVar
  factVarMod <- IDColsMEMod[factInd]

##
  fixInd <- which(IDColsME %in% fixedME)
  fixedMEMod <-  IDColsMEMod[fixInd]

  strainGeno <- rownames(genoDat)

##
  IDColsMEList <- list(IDColsME,IDColsMEMod)
  names(IDColsMEList) <- c("IDCols","IDColsMod")
  
### 
  DT_2$Obs <- DT_2$value
  ord <- order(DT_2$env, DT_2$gid)
  DT_2_ord <- DT_2[ord, , drop = FALSE]
  
### 
    
  if(CVMet=="CV_LOFO"){
	   print(paste("Running Leave One ",factVarMod," Out CV",sep=""))
	   factorsVar <- levels(factor(DT_2[,factVarMod]))
  }else if(CVMet!= "CV_LOFO"){
	   print(paste("Running ",CVMet,"cross validation",sep=""))
  }
 
###

 
	ME_CV <- fitMEModels_CV_fits(DT_2_ord,genoDat,strainGeno,KG=KG,KE=ke,Fixed=fixedMEMod,CVMet,K,NIter,FactVar=factVarMod,method=KMethod,Reaction= Reaction,dimension_KE,bandwidth = bandwidth,quantil = quantil,trait = NULL,store_full_fits = FALSE, b_iter = b_iter,b_burnin = b_burnin,b_thin = b_thin,b_tol = b_tol,b_R2 = b_R2, b_digits = b_digits)
	
	ME_Out_CV_Trt[[nTrt]] <- ME_CV
   
    DT_2_List[[nTrt]] <- DT_2_ord
	IDColsMEList_Out[[nTrt]] <- IDColsMEList
  }
  
  outList <- list(ME_Out_CV_Trt,DT_2_List,IDColsMEList_Out)
  names(outList) <- c("results","DT","IDColsME")
  return(outList)
}




#' Obtain test/train data for the various cross validation schemes
#'
#' A method to split input data to obtain test and training data for various cross validation schemes
#'
#' @param DT_2B  Data table which is input to fitMEModels .
#' @param CVMet Character. Cross-validation scheme. Default = "CV1".
#'   Options: "CV1", "CV2", "CV0", "CV00","CV-LOFO".
#' @param K Integer. Number of folds in k-fold cross-validation.
#' @param NIter Integer. Number of replicates for cross-validation.
#' @param id_col Character. Variable specifying the colname with strain id with default value of "gid".
#' @param env_col Character. Variable specifying the colname with environments with default value of "env"..
#' @param y_col Character. Variable specifying the colname with trait value 
#' @param base_seed seed value for set.seed for sampling 
#' @param cv00_partial_prop proportion for masking in soft environments with default value set to 0.5.
#'
#' @return a list of list containing the data table, train and test indices for each of the fold for each of the reps. 
#' In addition, the CV method and unique environments in the data are printed.
#' @examples
#' \dontrun{
#' # Example usage of getTstIndices_CV_Data
#' result <- getTstIndices_CV_Data(...)
#' }

  getTstIndices_CV_Data <- function(DT_2B, CVMet, nIter = 5, k = 5,cvFactVar = NULL,
                                   id_col = "gid", env_col = "env", y_col = "value",
                                   base_seed = 125L,
                                   cv00_partial_prop = 0.2){

 
   
## CV1, where novel genotypes in tested environments are predicted.
## CV2, where tested genotypes in tested environments are predicted.
## CV0, where tested genotypes in untested novel environments are predicted.
## CV00, where novel genotypes in novel environments are predicted.
## CV LOFO (Leave One Factor Out), eg: Leave One Test Out/ Leave One Line Out cross validation.
   
  print("TstInd_In")

  stopifnot(is.data.frame(DT_2B))
  stopifnot(all(c(id_col, env_col, y_col) %in% names(DT_2B)))
  if (k < 2 && CVMet %in% c("CV1","CV2","CV0","CV00")) {
    stop("k must be >= 2 for the requested CV method.")
  }

  N <- nrow(DT_2B)

  # ---- helpers ----
  # split a set of integer indices (1..n) into k disjoint, near-equal folds
  .split_kfold <- function(n, k) {
    perm <- sample.int(n)
    base <- n %/% k; r <- n %% k
    sizes <- rep(base, k); if (r > 0) sizes[seq_len(r)] <- sizes[seq_len(r)] + 1L
    split(perm, rep(seq_len(k), sizes))
  }
    

  split_kfold_indices <- function(items, k){
	  m <- length(items)
	  out <- vector("list", k)
	  if (m == 0L) return(out)

	  perm  <- sample(items, m)
	  base  <- m %/% k
	  r     <- m %% k

	  sizes <- rep.int(base, k)
	  if (r > 0L) {
		# Randomly choose which folds receive the r extra rows (spreads load across folds)
		extra <- sample.int(k, r)
		sizes[extra] <- sizes[extra] + 1L
	  }

	  idx <- rep.int(seq_len(k), sizes)
	  tmp <- split(perm, idx)

	  for (f in seq_len(k)) out[[f]] <- if (!is.null(tmp[[as.character(f)]])) as.integer(tmp[[as.character(f)]]) else integer(0L)
	  out
   }

 

  # Make containers (rep x fold)
  DT_masked_rep <- vector("list", nIter)
  train_idx_rep <- vector("list", nIter)
  test_idx_rep  <- vector("list", nIter)

  # ============================================================
  # CV1: leave-genotypes-out across ALL environments
  # Test = all rows of the held-out genotype set
  # ============================================================
  if (CVMet == "CV1") {
    unq_ids <- unique(DT_2B[[id_col]])
    n_ids   <- length(unq_ids)
    if (k > n_ids) stop("k cannot exceed the number of unique genotypes (CV1).")

    for (rep in seq_len(nIter)) {
      set.seed(base_seed + rep)

      # disjoint folds over genotypes
      idx_in_ids <- .split_kfold(n_ids, k)

      # per-fold results
      DT_list   <- vector("list", k)
      tr_list   <- vector("list", k)
      tst_list  <- vector("list", k)

      for (f in seq_len(k)) {
        test_ids_idx <- sort(idx_in_ids[[f]])
        test_ids     <- unq_ids[test_ids_idx]

        test_rows <- which(DT_2B[[id_col]] %in% test_ids)
        train_rows <- setdiff(seq_len(N), test_rows)

        DT_f <- DT_2B
        DT_f[[y_col]][test_rows] <- NA

        DT_list[[f]]  <- DT_f
        tr_list[[f]]  <- train_rows
        tst_list[[f]] <- test_rows
      }

      DT_masked_rep[[rep]] <- DT_list
      train_idx_rep[[rep]] <- tr_list
      test_idx_rep [[rep]] <- tst_list
    }

    return(list(DT_masked = DT_masked_rep,
                train_idx = train_idx_rep,
                test_idx  = test_idx_rep,
                cv = "CV1",
                id_lookup = unq_ids))
  }

  # ============================================================
  # CV2: predict missing G×E cells; each genotype remains observed somewhere
  # Test = a subset of rows per genotype; small m (<=1) never masked
  # Works for unbalanced data (different m per genotype)
  # ============================================================
  if (CVMet == "CV2") {
    # Pre-split rows by genotype
    rows_by_gid <- split(seq_len(N), DT_2B[[id_col]])

    for (rep in seq_len(nIter)) {
      set.seed(base_seed + rep)

      # fold assignment per row (0 = always train; 1..k = test fold)
      fold_assign <- integer(N)

      for (g in names(rows_by_gid)) {
        rows_g <- rows_by_gid[[g]]
        m <- length(rows_g)
        if (m <= 1L) {
          # not enough rows to mask -> keep in training for all folds
          fold_assign[rows_g] <- 0L
        } else {
          
		  # distribute rows of this genotype over folds (disjoint)
          # perm <- sample(rows_g, m)
          # base <- m %/% k; r <- m %% k
          # sizes <- rep(base, k); if (r > 0) sizes[seq_len(r)] <- sizes[seq_len(r)] + 1L
          # splits <- split(perm, rep(seq_len(k), sizes))
		  
		  splits <- split_kfold_indices(rows_g, k)
 	  
		  
          for (f in seq_len(k)) {
            if (length(splits[[f]]) > 0L) fold_assign[splits[[f]]] <- f
          }
        }
      }

      DT_list   <- vector("list", k)
      tr_list   <- vector("list", k)
      tst_list  <- vector("list", k)

      for (f in seq_len(k)) {
        test_rows <- which(fold_assign == f)
        train_rows <- setdiff(seq_len(N), test_rows)

        DT_f <- DT_2B
        if (length(test_rows)) DT_f[[y_col]][test_rows] <- NA

        DT_list[[f]]  <- DT_f
        tr_list[[f]]  <- train_rows
        tst_list[[f]] <- test_rows
      }

      DT_masked_rep[[rep]] <- DT_list
      train_idx_rep[[rep]] <- tr_list
      test_idx_rep [[rep]] <- tst_list
    }

    return(list(DT_masked = DT_masked_rep,
                train_idx = train_idx_rep,
                test_idx  = test_idx_rep,
                cv = "CV2"))
  }

  # ============================================================
  # CV0: leave-environment-out
  # Test = ALL rows from a held-out environment fold
  # Handles unbalanced presence/absence of genotypes by construction
  # ============================================================
  if (CVMet == "CV0") {
    unq_env <- unique(DT_2B[[env_col]])
    n_env   <- length(unq_env)
    if (k > n_env) stop("k cannot exceed the number of unique environments (CV0).")

    for (rep in seq_len(nIter)) {
      set.seed(base_seed + rep)

      # disjoint folds over environments
      idx_in_env <- .split_kfold(n_env, k)

      DT_list   <- vector("list", k)
      tr_list   <- vector("list", k)
      tst_list  <- vector("list", k)

      for (f in seq_len(k)) {
        test_env_idx <- sort(idx_in_env[[f]])
        test_envs    <- unq_env[test_env_idx]

        test_rows <- which(DT_2B[[env_col]] %in% test_envs)
        train_rows <- setdiff(seq_len(N), test_rows)

        DT_f <- DT_2B
        if (length(test_rows)) DT_f[[y_col]][test_rows] <- NA

        DT_list[[f]]  <- DT_f
        tr_list[[f]]  <- train_rows
        tst_list[[f]] <- test_rows
      }

      DT_masked_rep[[rep]] <- DT_list
      train_idx_rep[[rep]] <- tr_list
      test_idx_rep [[rep]] <- tst_list
    }

    return(list(DT_masked = DT_masked_rep,
                train_idx = train_idx_rep,
                test_idx  = test_idx_rep,
                cv = "CV0",
                env_lookup = unq_env))
  }

  # ============================================================
  # CV00: “hard + partial” environment masking (variant you coded)
  # For each env fold f:
  #   - mask ALL rows from the envs in fold f (hard removal)
  #   - PLUS randomly mask a proportion (cv00_partial_prop) of genotypes 
  #     from the remaining envs (partial removal) to create novel geno x novel environment
  # Works for unbalanced data.
  # ============================================================
  if (CVMet == "CV00") {
	  if (!(cv00_partial_prop > 0 && cv00_partial_prop < 1)) {
		stop("cv00_partial_prop must be in (0,1).")
	  }
	  N <- nrow(DT_2B)
	  unq_env <- sort(unique(DT_2B[[env_col]]))
	  n_env   <- length(unq_env)
	  if (k > n_env) stop("k cannot exceed the number of unique environments (CV00).")

	  DT_masked_rep <- vector("list", nIter)
	  train_idx_rep <- vector("list", nIter)
	  test_idx_rep  <- vector("list", nIter)
	  novel_gid_rep <- vector("list", nIter)

	  for (rep in seq_len(nIter)) {
		set.seed(base_seed + rep)
		
		env_folds <- .split_kfold(n_env, k)

		DT_list <- vector("list", k)
		tr_list <- vector("list", k)
		tst_list <- vector("list", k)
		gid_list <- vector("list", k)

		for (f in seq_len(k)) {
		  hard_envs <- unq_env[ sort(env_folds[[f]]) ]
		  soft_envs <- setdiff(unq_env, hard_envs)

		  hard_rows <- which(DT_2B[[env_col]] %in% hard_envs)
		  soft_rows <- which(DT_2B[[env_col]] %in% soft_envs)
		  
		  ensure_overlap <- TRUE
		
		  # Candidate genotypes to potentially make novel
		  gids_soft <- unique(DT_2B[[id_col]][soft_rows])
		  if (ensure_overlap) {
			# keep only those that also occur in the hard envs so we can compute metrics
			gids_hard <- unique(DT_2B[[id_col]][hard_rows])
			gids_cand <- intersect(gids_soft, gids_hard)
		  } else {
			gids_cand <- gids_soft
		  }

          n_gids_hard <- length(gids_hard)

		  n_cand <- length(gids_cand)
		  n_mask <- max(0L, floor(n_cand * cv00_partial_prop))
		  masked_gids <- if (n_mask > 0) sample(gids_cand, n_mask, replace = FALSE) else character(0)

		  # Rows of masked genotypes across ALL envs → remove from training (keep them novel)
		  rm_rows_all <- if (length(masked_gids)) which(DT_2B[[id_col]] %in% masked_gids) else integer(0)

		  # Evaluation set = masked genotypes × HARD (left-out) envs
		  test_rows <- if (length(masked_gids)) {
			which(DT_2B[[env_col]] %in% hard_envs & DT_2B[[id_col]] %in% masked_gids)
		  } else integer(0)
		  
		  min_test_rows <- length(unique(DT_2B[[id_col]]))/(2*k)
		  

		  if (length(test_rows) < min_test_rows) {
			# not enough eval rows—skip masking in this fold
			masked_gids <- character(0)
			rm_rows_all <- integer(0)
			test_rows   <- integer(0)
		  }

		  # Mask values (training is the complement; models see NAs)
		  DT_f <- DT_2B
		  if (length(rm_rows_all)) DT_f[[y_col]][rm_rows_all] <- NA

		  train_rows <- setdiff(seq_len(N), rm_rows_all)

		  DT_list[[f]] <- DT_f
		  tr_list[[f]] <- train_rows
		  tst_list[[f]] <- sort(unique(test_rows))
		  gid_list[[f]] <- masked_gids
		}

		DT_masked_rep[[rep]] <- DT_list
		train_idx_rep[[rep]] <- tr_list
		test_idx_rep [[rep]] <- tst_list
		novel_gid_rep[[rep]] <- gid_list
	  }

	  return(list(
		DT_masked   = DT_masked_rep,
		train_idx   = train_idx_rep,
		test_idx    = test_idx_rep,
		novel_gids  = novel_gid_rep,
		cv          = "CV00",
		env_levels  = unq_env,
		partial_prop = cv00_partial_prop
	  ))
    }

####	
  if (CVMet == "CV_LOFO"){
	  fact_col <- which(colnames(DT_2B) %in% cvFactVar) 
	 
	  if (length(fact_col) ==0) {
		stop("Factor for leave-one-factor-level-out CV should be present in data table.")
	  }
	  
	  N <- nrow(DT_2B)
	  unq_env <- sort(unique(DT_2B[[env_col]]))
	  n_env   <- length(unq_env)
	  
	  unq_factrs <- levels(factor(DT_2B[[fact_col]]))
	  n_factrs   <- length(unq_factrs)
	  
	  nIter <- 1

	  DT_masked_rep <- vector("list", nIter)
	  train_idx_rep <- vector("list", nIter)
	  test_idx_rep  <- vector("list", nIter)
	  novel_gid_rep <- vector("list", nIter)
	  
	  for(rep in seq_len(nIter)){
	   set.seed(base_seed + rep)
     
       k <- n_factrs
	   fact_lev <- unq_factrs

       # per-fold results
       DT_list   <- vector("list", k)
       tr_list   <- vector("list", k)
       tst_list  <- vector("list", k)

       for (f in seq_len(k)) {
        
		test_rows <- which(DT_2B[[fact_col]] %in% fact_lev[f])
        train_rows <- setdiff(seq_len(N), test_rows)

        DT_f <- DT_2B
        DT_f[[y_col]][test_rows] <- NA

        DT_list[[f]]  <- DT_f
        tr_list[[f]]  <- train_rows
        tst_list[[f]] <- test_rows
      }

      DT_masked_rep[[rep]] <- DT_list
      train_idx_rep[[rep]] <- tr_list
      test_idx_rep [[rep]] <- tst_list
    }

    return(list(DT_masked = DT_masked_rep,
                train_idx = train_idx_rep,
                test_idx  = test_idx_rep,
                cv = "CV_LOFO",
                factor_levels = unq_factrs))

	}

  stop("Unknown CVMet. Use one of: 'CV1', 'CV2', 'CV0', 'CV00'.'CV_LOFO'")
} 



# split_kfold_indices (commented out — orphaned roxygen removed)
# split_kfold_indices <- function(items, k) {
	  # # items: vector of indices to split (e.g., row ids for a genotype)
	  # m <- length(items)
	  # if (m == 0L) return(replicate(k, integer(0L), simplify = FALSE))

	  # perm  <- sample(items, m, replace = FALSE)
	  # base  <- m %/% k
	  # r     <- m %% k

	  # sizes <- rep.int(base, k)
	  # if (r > 0L) sizes[seq_len(r)] <- sizes[seq_len(r)] + 1L

	  # idx <- rep.int(seq_len(k), sizes)        # length m (may omit some folds if size 0)
	  # tmp <- if (length(idx)) split(perm, idx) else list()

	  # # ensure length k with integer(0) for empty folds
	  # out <- vector("list", k)
	  # for (f in seq_len(k)) {
		# key <- as.character(f)
		# out[[f]] <- if (!is.null(tmp[[key]])) as.integer(tmp[[key]]) else integer(0L)
	  # }
	  # out
   # }
  
  
  
  

  
####### 

### Function to get kernel matrix for EMDe 

#' Title of get_kernel_MDe
#'
#' Description of what get_kernel_MDe does.
#'
#' @param Y Description of Y.
#' @param KG Description of KG.
#' @param KE Description of KE.
#' @param intercept.random=FALSE Description of intercept.random=FALSE.
#' @param dimension_KE=NULL Description of dimension_KE=NULL.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of get_kernel_MDe
#' result <- get_kernel_MDe(...)
#' }
get_kernel_MDe <- function(Y,KG,KE,intercept.random=FALSE,dimension_KE=NULL){
  
  ne = length(unique(Y$env))
  Zg <- stats::model.matrix(~0 + gid, Y)
  colnames(Zg) = gsub(colnames(Zg), pattern = "gid", replacement = "")
  ng <- length(unique(Y$gid))
  
  K_G = KG 
  K_E = KE
  
  model_b <- "MDe"
  
 # Kde = getK(Y = Y, setKernel = K_G, model = model_b, intercept.random = intercept.random) 
  
  Kde = getK(Y = Y, X= X,kernel = "GK", model = model_b, intercept.random = intercept.random) 
  
  
  names(Kde)[grep(names(Kde),pattern="^G$")] <- paste0("KG_", names(Kde)[grep("G",names(Kde))])
  names(Kde)[grep("KG_G$",names(Kde),invert = T)] <- paste0("KGE_", names(Kde)[grep("KG_G$",names(Kde),invert = T)])
  
  
  obs_GxE = paste0(Y$env, ":", Y$gid) 
  
  K <- Kde
  for (j in 1:length(K)) colnames(K[[j]]$Kernel) = rownames(K[[j]]$Kernel) = obs_GxE
  
  if(!is.null(KE)){
    if (is.null(dimension_KE))
      dimension_KE <- "q"
    
    if (dimension_KE == "q") {
      K_Em = list()
      Jg = matrix(1, ncol = ng, nrow = ng)
      colnames(Jg) = rownames(Jg) = levels(Y$gid)
      for (q in 1:length(K_E)) K_Em[[q]] = kronecker(K_E[[q]], 
                                                     Jg, make.dimnames = T)
      for (q in 1:length(K_Em)) K_Em[[q]] = K_Em[[q]][rownames(K_Em[[q]]) %in% 
                                                        obs_GxE, colnames(K_Em[[q]]) %in% obs_GxE]
    }
    
    if (dimension_KE == "n") {
      K_Em <- K_E
      for (j in 1:length(K_Em)) colnames(K_Em[[j]]) = rownames(K_Em[[j]]) = obs_GxE
      for (q in 1:length(K_Em)) K_Em[[q]] = K_Em[[q]][rownames(K_Em[[q]]) %in% 
                                                        obs_GxE, colnames(K_Em[[q]]) %in% obs_GxE]
    }
    
    
    h <- length(K_E)
    n <- length(K)
    K_e <- c()
    for (q in 1:h) K_e[[q]] = list(Kernel = K_Em[[q]], Type = "D")
    names(K_e) <- paste0("KE_", names(K_E))
    
    K_f <- Map(c, c(K, K_e)) 
    
  }else {K_f <- K } 
  
  return(K_f)
}




#### Function to extract variance components and estimate CI..

#' Title of Vcomp.BGGE
#'
#' Description of what Vcomp.BGGE does.
#'
#' @param model Description of model.
#' @param env Description of env.
#' @param gid Description of gid.
#' @param digits = digits Description of digits = digits.
#' @param alfa = 0.1 Description of alfa = 0.1.
#'
#' @return .
#' @examples
#' \dontrun{
#' # Example usage of Vcomp.BGGE
#' result <- Vcomp.BGGE(...)
#' }

Vcomp.BGGE <- function(model, env, gid, digits = digits,alfa = 0.1){
  t = length(unique(gid))
  e = length(unique(env))
  n = t * e
  K = model$K
  size = length(K)
  comps = data.frame(matrix(NA, ncol = 3, nrow = size))
  VarE = data.frame(matrix(NA, ncol = 3, nrow = 1))
  names(comps) = names(VarE) = c("K", "Var", "SD.var")
  for (k in 1:size) {
    comps[k, 1] = names(K)[k]
    comps[k, 2] = round(K[[k]]$varu, digits)
    comps[k, 3] = round(K[[k]]$varu.sd, digits)
  }
  VarE[1, 1] = "Residual"
  VarE[1, 2] = round(model$varE, digits)
  VarE[1, 3] = round(model$varE.sd, digits)
  comps = rbind(comps, VarE)
  comps$Type <- comps$K
  comps$Type[grep(comps$K, pattern = "KGE_")] = "GxE"
  comps$Type[grep(comps$K, pattern = "KG_")] = "Genotype (G)"
  # Override: KG_GE* kernels are GxE interactions, not main G
  # KG_GE -> "GxE", KG_GEbw1 -> "GxE_bw1", KG_GE_Env1 -> "GxE_Env1"
  gxe_idx <- grep(comps$K, pattern = "KG_GE")
  if (length(gxe_idx)) {
    sfx <- sub("KG_GE_?", "", comps$K[gxe_idx])
    comps$Type[gxe_idx] <- ifelse(nchar(sfx) == 0, "GxE", paste0("GxE_", sfx))
  }
  comps$Type[grep(comps$K, pattern = "^E$")] = "GxE"
  comps$Type[grep(comps$K, pattern = "KE_")] = "Environment (E)"
  comps$CI_upper = NA
  comps$CI_lower = NA
  ENV = which(comps$Type %in% "Environment (E)")
  GID = which(comps$Type %in% "Genotype (G)")
  GE = which(comps$Type %in% "GxE")
  R = which(comps$Type %in% "Residual")

  #### Upper CI
  comps$CI_upper[ENV] = (n - e) * comps$Var[ENV]/qchisq((alfa/2),
                                                        n - e)
  comps$CI_upper[GID] = (n - t) * comps$Var[GID]/qchisq((alfa/2),
                                                        n - t)
  comps$CI_upper[GE] = (n - t - e) * comps$Var[GE]/qchisq((alfa/2),
                                                          n - t - e)
  comps$CI_upper[R] = (n - t - e) * comps$Var[R]/qchisq((alfa/2),
                                                        n - t - e)

  #### Lower CI
  comps$CI_lower[ENV] = (n - e) * comps$Var[ENV]/qchisq((1 -
                                                           alfa/2), n - e)
  comps$CI_lower[GID] = (n - t) * comps$Var[GID]/qchisq((1 -
                                                           alfa/2), n - t)
  comps$CI_lower[GE] = (n - t - e) * comps$Var[GE]/qchisq((1 -
                                                             alfa/2), n - t - e)
  comps$CI_lower[R] = (n - t - e) * comps$Var[R]/qchisq((1 -
                                                           alfa/2), n - t - e)


  comps$CI_upper = round(comps$CI_upper, digits)
  comps$CI_lower = round(comps$CI_lower, digits)
  comps <- comps[, c(4, 1:2, 6, 5, 3)]

  p = c("KG_", "KE_", "KGE_")
  for (i in 1:3) comps$K = gsub(x = comps$K, pattern = p[i],
                                replacement = "")
  return(comps)
}


 
#' summarize_ME_CV_fits 
#'
#' Summarize output from fitMEModels_CV_fits
#' @param fit_obj output from fitMEModels_CV()
#' @param DT data table used for fitting the models in fitMEModels_CV
#' @param IDColsList a vector containing the ID cols in ME data table
#'
#' @return A list containing various metrics: predictive ability, RMSE,
#'   variance components, fit values, and summary of these metrics.
#' @export
#' @examples
#' \dontrun{
#' result <- summarize_ME_CV_fits(fit_obj, cvMet = "CV1")
#' }
 
 
summarize_ME_CV_fits <- function(
  fit_obj,          # output from fitMEModels_CV_fits()
  cvMet = CVMet
){
  # --- helpers ---
  safe_metrics <- function(pred, obs) {
    ok <- is.finite(pred) & is.finite(obs); n <- sum(ok)
    if (!n) return(list(cor = NA_real_, rmse = NA_real_, n_test = 0L))
	
    list(cor = if (n > 1) stats::cor(pred[ok], obs[ok],use="pairwise.complete.obs") else NA_real_,
         rmse = sqrt(mean((pred[ok] - obs[ok])^2)),
         n_test = n)
  }
  vc_df <- function(vc) {
    if (is.null(vc)) return(data.frame(component=character(0), value=numeric(0)))
    if (is.numeric(vc) && is.null(dim(vc))) {
      comps <- names(vc); if (is.null(comps)) comps <- paste0("VC", seq_along(vc))
      return(data.frame(component = comps, value = as.numeric(vc), check.names = FALSE))
    }
    if (is.data.frame(vc)) {
      nm <- tolower(names(vc))
      comp_col <- if ("component" %in% nm) names(vc)[match("component", nm)] else names(vc)[1]
      num_cols <- names(vc)[vapply(vc, is.numeric, logical(1))]; if (!length(num_cols)) return(data.frame(component=character(0), value=numeric(0)))
      return(data.frame(component = as.character(vc[[comp_col]]),
                        value = as.numeric(vc[[num_cols[1]]]), check.names = FALSE))
    }
    if (is.list(vc)) {
      v <- unlist(vc, use.names = TRUE); comps <- names(v); if (is.null(comps)) comps <- paste0("VC", seq_along(v))
      return(data.frame(component = comps, value = as.numeric(v), check.names = FALSE))
    }
    data.frame(component=character(0), value=numeric(0))
  }
 
  # --- accumulators ---
  metrics_all <- list(); var_all <- list(); fits_all <- list(); env_metrics_all <- list()

  nReps <- length(fit_obj$results)
  N     <- fit_obj$nrow
  IDCols <- c("gid","env")

  for (irep in seq_len(nReps)){
    folds <- fit_obj$results[[irep]]$folds
    nFolds <- length(folds)
	met_dat_list_nF <- list()
	
    for (fold in seq_len(nFolds)) {
      item <- folds[[fold]]
      if (!length(item) || is.null(item$status) || item$status != "ok") next
      
	  tst_idx <- as.integer(item$tst_idx)
      
	  # clamp to valid DT rows
      tst_idx <- tst_idx[is.finite(tst_idx) & tst_idx >= 1L & tst_idx <= N]
      if (!length(tst_idx)) next
	  
	  met_dat_list <- list()

	    for (m in names(item$models)){

		 # VC
		 vc_long <- vc_df(item$models[[m]]$VarComp)
		  
		 if (nrow(vc_long)) {
			vc_long$replicate <- irep; vc_long$fold <- fold; vc_long$model <- m
			vc_long <- vc_long[, c("replicate","fold","model","component","value")]
			var_all[[length(var_all)+1L]] <- vc_long
		 }
		
		 ##
		
		  yHat <- item$models[[m]]$yHat
          DT <- item$models[[m]]$dat
		  
		  # Coerce to numeric vector and guard length
		  yHat <- as.numeric(yHat)
		  if (!length(yHat)) next
		  
		  ok_idx <- tst_idx[tst_idx <= length(yHat)]
		  if (!length(ok_idx)) next
		  ok_idx <- ok_idx[is.finite(yHat[ok_idx])]
		  if (!length(ok_idx)) next

		  # Predictions & obs aligned
		  pred <- yHat[ok_idx]
		  obs  <- DT$Obs[ok_idx]
		  
		  if (!length(pred) || !length(obs)) next
	
          met_dat <- cbind.data.frame(pred, obs, env = as.character(DT$env[ok_idx]))
		  met_dat_list[[m]] <- met_dat
		}
		met_dat_list_nF[[fold]] <- met_dat_list

	}


  if(cvMet != "CV_LOFO"){
    met_data_list <- list()
    for(m in names(item$models)){

		   met_data <- do.call(rbind.data.frame,lapply(met_dat_list_nF, function(x) x[[m]]))
		   met_data_list[[m]] <- met_data

		   # Per-env PA across all folds combined
		   envs <- unique(met_data$env)
		   env_met_list <- lapply(envs, function(e) {
		     idx <- met_data$env == e
		     safe_metrics(met_data$pred[idx], met_data$obs[idx])
		   })
		   names(env_met_list) <- envs

		   mean_cor  <- mean(vapply(env_met_list, function(x) x$cor,    numeric(1)), na.rm = TRUE)
		   mean_rmse <- mean(vapply(env_met_list, function(x) x$rmse,   numeric(1)), na.rm = TRUE)
		   n_total   <- sum(vapply(env_met_list,  function(x) x$n_test, integer(1)))

		   metrics_all[[length(metrics_all)+1L]] <-
			data.frame(replicate = irep, model = m,
					   cor = mean_cor, rmse = mean_rmse,
					   n_test = n_total, n_envs = length(envs),
					   stringsAsFactors = FALSE)

		   # Per-env breakdown rows
		   for (e in envs) {
		     em <- env_met_list[[e]]
		     env_metrics_all[[length(env_metrics_all)+1L]] <- data.frame(
		       replicate = irep, model = m, env = e,
		       cor = em$cor, rmse = em$rmse, n_test = em$n_test,
		       stringsAsFactors = FALSE
		     )
		   }

		   # Obs/Pred table — all columns must have the same #rows
		   tab_core <- DT[ok_idx, IDCols, drop = FALSE]
		   if (nrow(tab_core)) {
			 tab <- cbind.data.frame(
			  tab_core,
			  Obs = obs,
			  Pred = pred,
			  model = rep(m, length(ok_idx)),
			  replicate = rep(irep, length(ok_idx)),
			  stringsAsFactors = FALSE
			 )
			 #if (length(IDCols)) names(tab)[seq_along(IDCols)] <- IDColsMod
			 fits_all[[length(fits_all)+1L]] <- tab
		   }
	}
	
   }else if(cvMet == "CV_LOFO"){ 
   
    met_data_list <- list()
	
    for(m in names(item$models)){ 
		 
		 met_data <- lapply(met_dat_list_nF, function(x) x[[m]])
		 met_data_list[[m]] <- met_data
		 		   
	     met_list <- list()
		 for(nLev in seq_along(met_data)){	
		 
           pred_dat <- met_data[[nLev]]$pred
		   obs_dat <- met_data[[nLev]]$obs
		   met  <- safe_metrics(pred_dat, obs_dat)
		   
		   met_list[[nLev]] <- met
		  
		  if(nLev ==1){
		   met_comb <- data.frame(factorLevel = nLev, model = m,
					   cor = met$cor, rmse = met$rmse, n_test = met$n_test,
					   stringsAsFactors = FALSE)
		  }else if(nLev >1){
					   
		   met_comb <- rbind.data.frame(met_comb,data.frame(factorLevel = nLev, model = m,
					   cor = met$cor, rmse = met$rmse, n_test = met$n_test,
					   stringsAsFactors = FALSE))
		  }
		  
		 }
		 metrics_all[[length(metrics_all)+1L]] <- met_comb

         	
 
		# Obs/Pred table — all columns must have the same #rows
		   tab_core <- DT[ok_idx, IDCols, drop = FALSE]
		   if (nrow(tab_core)) {
			 tab <- cbind.data.frame(
			  tab_core,
			  Obs = obs,
			  Pred = pred,
			  model = rep(m, length(ok_idx)),
			  replicate = rep(irep, length(ok_idx)),
			  stringsAsFactors = FALSE
			 )
			 #if (length(IDCols)) names(tab)[seq_along(IDCols)] <- IDColsMod
			 fits_all[[length(fits_all)+1L]] <- tab
		   }
       }
    }

}	
	metrics_long     <- if (length(metrics_all))     do.call(rbind, metrics_all)     else data.frame()
	var_long         <- if (length(var_all))         do.call(rbind, var_all)         else data.frame()
	fits_long        <- if (length(fits_all))        do.call(rbind, fits_all)        else data.frame()
	env_metrics_long <- if (length(env_metrics_all)) do.call(rbind, env_metrics_all) else data.frame()

    metrics_summary <- if (nrow(metrics_long)) {
     aggregate(cbind(cor, rmse) ~ model, data = metrics_long, function(x)
      c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), n = sum(is.finite(x))))
    }else data.frame()


 level_mapping_mod <- setNames(c("G+E", "G+E+GxE", "G+E+GxEi", "G+E+W", "G+E+GxE+W", "G+E+GxEi+W","G+E+W+GxW","G+E+GxE+W+GxW","G+E+GxEi+W+GxW"),
                              c("MM", "MDs", "MDe", "E-MM", "E-MDs", "E-MDe","RN-MM","RN-MDs","RN-MDe"))

 metrics_long$model <- factor(level_mapping_mod[metrics_long$model])
 var_long$model <- factor(level_mapping_mod[var_long$model])
 fits_long$model <- factor(level_mapping_mod[fits_long$model])
 if (nrow(env_metrics_long)) {
   env_metrics_long$model <- factor(level_mapping_mod[env_metrics_long$model])
 }
 
 
# safe builder for metrics_summary from aggregate(cbind(cor, rmse) ~ model, ...)
  build_metrics_summary <- function(ms) {
   if (!nrow(ms)) return(data.frame())

   to_mat <- function(x) {
    if (is.matrix(x)) return(x)
    if (is.list(x))  return(do.call(rbind, x))
    # single named vector
    m <- matrix(as.numeric(x), nrow = 1)
    colnames(m) <- names(x)
    m
   }

  cor_mat  <- to_mat(ms$cor)
  rmse_mat <- to_mat(ms$rmse)

  # ensure expected colnames
  if (is.null(colnames(cor_mat)))  colnames(cor_mat)  <- c("mean","sd","n")[seq_len(ncol(cor_mat))]
  if (is.null(colnames(rmse_mat))) colnames(rmse_mat) <- c("mean","sd","n")[seq_len(ncol(rmse_mat))]

  data.frame(
    model      = ms$model,
    cor_mean   = cor_mat[,  "mean", drop = TRUE],
    cor_sd     = cor_mat[,  "sd",   drop = TRUE],
    cor_n      = cor_mat[,  "n",    drop = TRUE],
    rmse_mean  = rmse_mat[, "mean", drop = TRUE],
    rmse_sd    = rmse_mat[, "sd",   drop = TRUE],
    rmse_n     = rmse_mat[, "n",    drop = TRUE],
    row.names  = NULL
  )
}

# Build metrics summary
metrics_summary <- if (nrow(metrics_long)) {
  agg <- aggregate(cbind(cor, rmse) ~ model, data = metrics_long, function(x)
    c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), n = sum(is.finite(x))))
  build_metrics_summary(agg)
} else data.frame()


metrics_summary$model <- factor(level_mapping_mod[metrics_summary$model]) 

 print(metrics_summary)

  list(
    metrics_long     = metrics_long,
    env_metrics_long = env_metrics_long,
    var_long         = var_long,
    fits_long        = fits_long,
    metrics_summary  = metrics_summary,
    cv               = fit_obj$cv
  )

}




#' Fit multi-environment genomic prediction models for cross validation schemes 
#'
#' @param DT Data frame with phenotypic data.
#' @param genoDat data frame with genotypic data.
#' @param strainGeno Character vector containing strain IDs of genotypes.
#' @param KG Matrix or NULL. Genotypic relationship matrix. Default = NULL, in which case
#'   the kernel is estimated internally. A user-defined matrix can be supplied.
#' @param KE Matrix or NULL. Environmental relationship matrix. Default = NULL.
#'   Required if `FitEnvModels = TRUE`.
#' @param CVMet Character. Cross-validation scheme. Default = "CV1".
#'   Options: "CV1", "CV2", "CV0", "CV00","CV-LOFO".
#' @param K Integer. Number of folds in k-fold cross-validation.
#' @param NIter Integer. Number of replicates for cross-validation.
#' @param method Character. Kernel method for estimating genotypic and
#'   genotype × environment kernels. Default = "Linear". Option: "Gaussian".
#' @param fitEnvModels Logical. Whether to include environmental covariates in the model.
#'   Default = FALSE. If TRUE, a user-provided `KE` must be supplied.
#' @param trait Character. Variable specifying the trait name.
#' @param store_full_fits logical value to store full fits in results. Default value is set to FALSE.
#' @importFrom BGGE BGGE
#' @importFrom EnvRtype get_kernel kernel_model
#' @export
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of fitMEModels_CV_fits
#' result <- fitMEModels_CV_fits(...)
#' }

  fitMEModels_CV_fits <- function(
	  DT, genoDat, strainGeno,
	  KG = NULL, KE = NULL,Fixed,
	  CVMet, k, nIter,FactVar,method,
	  Reaction = FALSE, dimension_KE=NULL,
	  bandwidth = 1,quantil = 0.5,  
	  trait = NULL,                # e.g., "Maturity.rating" to round predictions
	  store_full_fits = FALSE ,     # set TRUE only if you really need raw fit objects (memory-heavy)
      b_iter = 1000,b_burnin = 200,b_thin = 10,b_tol = 1e-10,b_R2 = 0.5, b_digits =4
  ){
  
  requireNamespace("foreach", quietly = TRUE)
  print("MECV Fit - in")

  DT_2A <- DT
  DT_2B <- droplevels(DT)
  # CV splits
  
  Dat_Out <- getTstIndices_CV_Data(DT_2B, CVMet, nIter, k,cvFactVar=FactVar)
  print(is.null(Dat_Out))
  
  nReps   <- length(Dat_Out$DT_masked)
  N       <- nrow(DT_2B)
  
 #If KG is provided, ignore GB/GK method
  
  if (!is.null(KG)) method <- NULL

# backend — explicit cluster avoids the doSNOW/.doSnowGlobals code path
cl_par <- makeCluster(min(max(detectCores() - 1, 1), 10))
registerDoParallel(cl_par)
on.exit(stopCluster(cl_par), add = TRUE)

# figure out folds per replication (robust even if they vary)
nFolds_by_rep <- vapply(Dat_Out$DT_masked, length, integer(1))
stopifnot(length(nFolds_by_rep) == nReps)
maxFolds <- max(nFolds_by_rep)

# index all (irep, fold) pairs that exist
idx <- do.call(rbind, lapply(seq_len(nReps), function(r) {
  cbind(irep = r, fold = seq_len(nFolds_by_rep[r]))
}))
idx <- as.data.frame(idx)

# run all tasks in parallel
 task_results <- foreach(k = seq_len(nrow(idx)),
                        .packages=c("BGGE","BGLR"),
                        .export = c("fit_kernel_models","get_kinship_kernels")
                       
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
    ofit <- fit_DTfold(
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

# rebuild: rep_results[[irep]] -> list(folds = list( ... per fold ... ))
  rep_results <- vector("list", nReps)

  for (irep in seq_len(nReps)){
	  nFolds <- nFolds_by_rep[irep]
	  fold_list <- vector("list", nFolds)

	  # pick results for this irep
	  items <- Filter(function(x) x$irep==irep,task_results)

	  # place by fold index
	  for (it in items){
		fold_list[[it$fold]] <- list(
		  status = it$status,
		  models = it$models,
		  raw    = it$raw,
		  tst_idx = it$tst_idx
		)
	  }

	  rep_results[[irep]] <- list(folds = fold_list)
  }

  #rep_results_c <- lapply(rep_results,function(x) do.call(c,lapply(x,function(y) y)))
  
  # Shape:
  # results[[rep]]$folds[[fold]]$models[["MM"|"MDs"|"MDe"| "SM"]] = list(yHat=..., VarComp=...)
  # results[[rep]]$folds[[fold]]$tst_idx = integer vector (original row indices)
  # results[[rep]]$folds[[fold]]$status / $error
  return(list(
    cv      = Dat_Out$cv,
    results = rep_results,
    nrow    = N
  ))
}




###

# helper to fit one fold and return predictions + VC (not raw fits unless asked)
 
  
  fit_DTfold <- function(DT_fold,genoDat,trait, KG, KE,fixedCovar,cvFactVar=NULL,Reaction, method,intcpt.rand=FALSE,dimension_KE=NULL,bandwidth = 1,quantil = 0.5,b_iter = 1000,b_burnin = 200,b_thin = 10,b_tol = 1e-10,b_R2 = 0.5, b_digits = 4) {
    
	print("fit_fold_in")
	Y   <- DT_fold[, c("env","gid","value")]
    y   <- "value"; gid <- "gid"; env <- "env"
    
	ne  <- length(unique(Y$env))
    X   <- if (!is.null(rownames(genoDat))) genoDat[which(rownames(genoDat) %in% Y$gid), , drop = FALSE] else genoDat
    KMethod <- method 
	fitModels <- list(); raw_store <- list()
	addCol <- cvFactVar
	
    if(ne>1){

        # --- Build kernels ---
	
		MM <- get_kinship_kernels(Y,X,K_G = KG, K_E = KE,intercept.random=intcpt.rand,model="MM",method=KMethod,dimension_KE=dimension_KE,bandwidth = bandwidth,quantil = quantil,reaction=Reaction)
		MDs <- get_kinship_kernels(Y,X,K_G = KG,K_E = KE,intercept.random=intcpt.rand,model="MDs",method=KMethod,dimension_KE=dimension_KE,bandwidth = bandwidth,quantil = quantil,reaction=Reaction)
		MDe <- get_kinship_kernels(Y,X,K_G = KG,K_E= KE,intercept.random=intcpt.rand,model="MDe",method=KMethod,dimension_KE=dimension_KE,bandwidth = bandwidth,quantil = quantil,reaction=Reaction)
		
       	### Set fixed effect matrix         
        lev_env <- levels(DT_fold$env)
        #fixed = model.matrix(as.formula(paste(" ~ 0+", fixedCovar,sep="")),DT_fold)
		
		fixed = model.matrix(as.formula(paste("~ 0 +", paste(fixedCovar, collapse = " + "))), DT_fold)

        NE_tab  <- table(DT_fold$env)
        NE <- as.integer(NE_tab[lev_env])  # aligned to design
	    uniqid <- "uniqID"
		
	   # --- Fit kernel models ---
       fit_MM  <- fit_kernel_models(y = y, env = env, gid = gid, uniqID=uniqid, Data = DT_fold, random = MM,  fixed = fixed, engine="BGGE",addCol= addCol,verbose=TRUE,iterations = b_iter, burnin = b_burnin, thining = b_thin, tol = b_tol, R2 = b_R2,digits= b_digits)
       fit_MDs <- fit_kernel_models(y = y, env = env, gid = gid, uniqID=uniqid, Data = DT_fold, random = MDs, fixed = fixed, engine="BGGE",addCol= addCol,verbose=TRUE,iterations = b_iter, burnin = b_burnin, thining = b_thin, tol = b_tol, R2 = b_R2,digits= b_digits)
       fit_MDe <- fit_kernel_models(y = y, env = env, gid = gid, uniqID=uniqid, Data = DT_fold, random = MDe, fixed = fixed, engine="BGGE",addCol= addCol,verbose=TRUE,iterations = b_iter, burnin = b_burnin, thining = b_thin, tol = b_tol, R2 = b_R2,digits= b_digits)
       
	 
       fitModels$MM  <- list(yHat = fit_MM$yHat,dat=fit_MM$dat, VarComp = fit_MM$VarComp)
	   fitModels$MDs <- list(yHat = fit_MDs$yHat,dat=fit_MDs$dat, VarComp = fit_MDs$VarComp)
	   fitModels$MDe <- list(yHat = fit_MDe$yHat,dat=fit_MDe$dat, VarComp = fit_MDe$VarComp)

	  if (store_full_fits) raw_store <- list(MM = fit_MM, MDs = fit_MDs, MDe = fit_MDe)
 	   
	}else if(ne==1){
        SM <- get_kinship_kernels(Y,X,KG,KE,intercept.random=intcpt.rand,model="SM",method=KMethod,dimension_KE=dimension_KE,bandwidth = bandwidth,quantil = quantil,reaction=Reaction)
        fixed <- model.matrix(~1, DT_fold)
		uniqid <- "uniqID"
		fit_SM <- fit_kernel_models(y = y, env = env, gid = gid,uniqID = uniqid,Data = DT_fold, random = SM, fixed = fixed,engine="BGGE",addCol= addCol,verbose=TRUE,iterations = b_iter, burnin = b_burnin, thining = b_thin, tol = b_tol, R2 = b_R2,digits= b_digits)
        fitModels$SM <- list(yHat = fit_SM$yHat,dat=fit_SM$dat, VarComp = fit_SM$VarComp)
        if (store_full_fits) raw_store <- list(SM = fit_SM)
    }
	
	# Optional rounding for ordinal trait
    if (identical(trait, "Maturity.rating")) {
      for (m in names(fitModels)) fitModels[[m]]$yHat <- round(fitModels[[m]]$yHat)
	  
    }
  
   list(models = fitModels, raw = if (store_full_fits) raw_store else NULL)
  }



# ================================================================
# get_kinship_kernel
# Standardized kernel construction for multi-environment GP
# ================================================================
# Y: data.frame with columns "env", "gid" (and optional response)
# X: optional marker matrix (rownames = gid); used when K_G is NULL
# K_G: optional named list of genotype kernels (gid x gid)
# K_E: optional named list of environment kernels:
#      - if dimension_KE = "q": env x env
#      - if dimension_KE = "n": obs x obs (rownames/colnames = paste0(env, ":", gid))
# model: one of c("SM","MM","MDs","MDe","EMM","EMDs","EMDe","RNMM","RNMDs","RNMDe")
# method: "GB" or "GK" (only used when K_G is NULL and X is provided)
# dimension_KE: "q" or "n" (default "q" if K_E given)
# bandwidth: numeric or vector (for GK)
# quantil: numeric in (0,1) (for GK distance scale)
# reaction: optional logical override; if NULL it is inferred from 'model'
# intercept.random: add KG_Gi random intercept (Type="D")
#
# Returns: a flat named list of kernel entries:
#   Each entry: list(Kernel = <n x n matrix>, Type = "D" or "BD")
#

get_kinship_kernels <- function(
  Y,
  X = NULL,
  K_G = NULL,
  K_E = NULL,
  intercept.random = FALSE,
  model = "MM",
  method = c("GB","GK"),
  dimension_KE = NULL,
  bandwidth = 1,
  quantil = 0.5,
  reaction = NULL
) {

  # ------------------------- checks & setup -------------------------
  stopifnot(all(c("env","gid") %in% names(Y)))
  Y <- as.data.frame(Y)
  Y$env <- factor(Y$env)
  Y$gid <- factor(Y$gid)

  n  <- nrow(Y)
  env_lv <- levels(Y$env)
  gid_lv <- levels(Y$gid)
  nEnv   <- length(env_lv)
  nGid   <- length(gid_lv)

  # Design matrices (observation order!)
  Ze <- stats::model.matrix(~ 0 + Y$env)  ; colnames(Ze) <- env_lv   # n x nEnv
  Zg <- stats::model.matrix(~ 0 + Y$gid)  ; colnames(Zg) <- gid_lv   # n x nGid
  obs <- paste0(Y$env, ":", Y$gid)

  # Warn if repeated (gid, env) cells
  if (any(table(Y$env, Y$gid) > 1)) {
    warning("There are replicated (env, gid) cells; kernels are built on the row-level observation index.")
  }

  # Parse model -> base and reaction
  .parse_model <- function(m, rxn) {
    m <- match.arg(m, c("SM","MM","MDs","MDe","EMM","EMDs","EMDe","RNMM","RNMDs","RNMDe"))
    if (grepl("^RN", m)) {
      base <- sub("^RN", "", m); rn <- TRUE
    } else if (grepl("^EM", m)) {
      base <- sub("^E", "", m); rn <- FALSE
    } else {
      base <- m; rn <- isTRUE(rxn)  # allow manual override only for bare models
    }
    list(base = base, rn = rn)
  }
  
  pm <- .parse_model(model, reaction)
  base_model <- pm$base   # "SM"|"MM"|"MDs"|"MDe"
  do_rn      <- pm$rn

  method <- match.arg(method)

  # Helper: align matrix to obs order using dimnames (for KE "n")
  .align_to_obs <- function(M, obs_names) {
    if (is.null(rownames(M)) || is.null(colnames(M))) {
      stop("Alignment requires row/col names on the matrix.")
    }
    ix <- match(obs_names, rownames(M))
    iy <- match(obs_names, colnames(M))
    if (anyNA(ix) || anyNA(iy)) stop("Obs names of KE (dimension_KE='n') do not match Y order.")
    M[ix, iy, drop = FALSE]
  }

  # ------------------------- build genetic side (G) -------------------------
  G_list <- list()

  if (!is.null(K_G)) {
    if (is.matrix(K_G)) K_G <- list(G = K_G)
    if (is.null(names(K_G))) names(K_G) <- paste0("G", seq_along(K_G))

    for (g in seq_along(K_G)) {
      Kg <- K_G[[g]]
      if (is.null(rownames(Kg)) || is.null(colnames(Kg))) {
        stop("Each K_G matrix must have gid row/col names.")
      }
      
	  Kg <- Kg[gid_lv, gid_lv, drop = FALSE]                       # reorder by subjects
      ## Kernel for user provided G
	  Gobs <- Zg %*% Kg %*% t(Zg)                                  # n x n
      rownames(Gobs) <- colnames(Gobs) <- obs
      nm <- if (length(K_G) == 1) "KG_G" else paste0("KG_G_", names(K_G)[g])
      G_list[[nm]] <- list(Kernel = Gobs, Type = "D")
    }

  } else if (is.null(K_G) & !is.null(X)) {
    if (is.null(rownames(X))) stop("X must have rownames equal to gid.")
    if (!all(gid_lv %in% rownames(X))) {
      miss <- setdiff(gid_lv, rownames(X)); stop("Genotypes missing in X: ", paste(miss, collapse=", "))
    }
    X <- X[gid_lv, , drop = FALSE]

    if (method == "GB") {
	## Kernel for GB method
      Kg <- tcrossprod(X) / ncol(X)                                # gid x gid
      Gobs <- Zg %*% Kg %*% t(Zg)
      rownames(Gobs) <- colnames(Gobs) <- obs
      G_list[["KG_G"]] <- list(Kernel = Gobs, Type = "D")
    } else {  
	 ## Gaussian Kernel 
      D <- as.matrix(stats::dist(X))^2
      if (length(bandwidth) < 1) bandwidth <- 1
      for (i in seq_along(bandwidth)) {
        Kg <- exp(-bandwidth[i] * D / stats::quantile(D, probs = quantil))
        Gobs <- Zg %*% Kg %*% t(Zg)
        rownames(Gobs) <- colnames(Gobs) <- obs
        nm <- if (length(bandwidth) == 1) "KG_G" else paste0("KG_G_bw", i)
        G_list[[nm]] <- list(Kernel = Gobs, Type = "D")
      }
    }

  } else {
    # default Ig
    Gobs <- Zg %*% diag(nGid) %*% t(Zg)                            # same-gid across envs
    rownames(Gobs) <- colnames(Gobs) <- obs
    G_list[["KG_G"]] <- list(Kernel = Gobs, Type = "D")
    message("K_G and X were NULL; using identity in gid space (KG_G).")
  }

  # ------------------------- base MDs / MDe GE -------------------------
  OUT <- G_list

  if (base_model == "MDs") {
  ## 
    E_all <- Ze %*% t(Ze)  # == tcrossprod(Ze)
    for (nm in names(G_list)) {
      OUT[[paste0("KG_GE", sub("^KG_G_?", "", nm))]] <-
        list(Kernel = G_list[[nm]]$Kernel * E_all, Type = "BD")
    }
  } else if (base_model == "MDe") {
    for (nm in names(G_list)) {
      Gobs <- G_list[[nm]]$Kernel
	  ## Block-wise (BD) kernel for each of the envs
      for (i in seq_len(nEnv)) {
        Bi <- Ze[, i, drop = FALSE] %*% t(Ze[, i, drop = FALSE])   # per-env block
        OUT[[paste0("KG_GE_", env_lv[i], sub("^KG_G_?", "", nm))]] <-
          list(Kernel = Gobs * Bi, Type = "BD")
      }
    }
  } else if (base_model %in% c("MM","SM")) {
    # nothing extra; keep G only
  } else {
    stop("Unknown base model: ", base_model)
  }

  # ------------------------- environment side (KE) -------------------------
  E_obs <- list()
  if (!is.null(K_E)) {
    if (is.null(dimension_KE)) dimension_KE <- "q"
    if (is.matrix(K_E)) K_E <- list(E = K_E)
    if (is.null(names(K_E))) names(K_E) <- paste0("E", seq_along(K_E))

    if (dimension_KE == "q") {
      for (e in seq_along(K_E)) {
	  
        KEe <- K_E[[e]]
        if (is.null(rownames(KEe)) || is.null(colnames(KEe))) {
          stop("K_E[['", names(K_E)[e], "']] must have env row/col names when dimension_KE='q'.")
        }
		### Env Kernel for enviromics enriched models 
        KEe <- KEe[env_lv, env_lv, drop = FALSE]
        Eo  <- Ze %*% KEe %*% t(Ze)
        rownames(Eo) <- colnames(Eo) <- obs
        E_obs[[names(K_E)[e]]] <- Eo
        OUT[[paste0("KE_", names(K_E)[e])]] <- list(Kernel = Eo, Type = "D")
      }
    } else if (dimension_KE == "n") {
      for (e in seq_along(K_E)) {
        Eo <- .align_to_obs(K_E[[e]], obs)
        rownames(Eo) <- colnames(Eo) <- obs
        E_obs[[names(K_E)[e]]] <- Eo
        OUT[[paste0("KE_", names(K_E)[e])]] <- list(Kernel = Eo, Type = "D")
      }
    } else {
      stop("dimension_KE must be 'q' or 'n'.")
    }
  }

  # ------------------------- reaction-norm terms (RN*) -------------------------
  if (do_rn && length(E_obs) > 0) {
    # RNMM/RNMDs: dense across envs (Type = "D")
    if (base_model %in% c("MM","MDs")) {
      for (g in names(G_list)) {
        Gobs <- G_list[[g]]$Kernel
        gkey <- sub("^KG_G_?", "", g)
		
        for (e in names(E_obs)) {
          GE <- Gobs * E_obs[[e]]
		  
          OUT[[paste0("KGE", gkey, "_G", e)]] <- list(Kernel = GE, Type = "D")
        }
      }
    }
    # RNMDe: per-environment blocks (Type = "BD")
    if (base_model == "MDe") {
      for (g in names(G_list)) {
        Gobs <- G_list[[g]]$Kernel
        gkey <- sub("^KG_G_?", "", g)
        for (e in names(E_obs)) {
		### Reaction norm kernel 
          GE_full <- Gobs * E_obs[[e]]
		## BD kernels for RN 
          for (i in seq_len(nEnv)) {
            Bi <- Ze[, i, drop = FALSE] %*% t(Ze[, i, drop = FALSE])
            OUT[[paste0("KGE", gkey, "_G", e, "_", env_lv[i])]] <-
              list(Kernel = GE_full * Bi, Type = "BD")
          }
        }
      }
    }
  } else if (do_rn && length(E_obs) == 0) {
    message("Reaction-norm requested but K_E is NULL; no RN terms were added.")
  }

  # ------------------------- random intercept (optional) -------------------------
  if (isTRUE(intercept.random)) {
    Gi <- Zg %*% diag(nGid) %*% t(Zg)
    rownames(Gi) <- colnames(Gi) <- obs
    OUT[["KG_Gi"]] <- list(Kernel = Gi, Type = "D")
  }

  OUT
}
  
  
  
  

  
# ================================================================
# fit_kernel_models
#   Unified wrapper to fit multi-kernel, multi-environment models
#   using either BGGE (internal Gibbs sampler) or BGLR (RKHS).
#
# Inputs
#   y            : response column name (character) OR numeric vector
#   data         : data.frame with at least columns env, gid, and y (if y is a name)
#   random       : list of kernels from get_geno_kernel_v3()
#                  each element: list(Kernel=<n x n>, Type="D"|"BD")
#   fixed        : optional design matrix (n x p) of fixed effects (e.g. ~0+env)
#   env, gid     : column names in 'data' for environment and genotype
#   engine       : "BGGE" or "BGLR"
#   iterations   : MCMC total iterations
#   burnin       : burn-in iterations
#   thining      : thinning interval
#   tol, R2      : BGGE internals (spectral threshold; prior R2 split)
#   digits       : rounding for variance components table
#   verbose      : logical or integer for print cadence
#
# Returns
#   list(yHat, varE, random, fit, VarComp, engine)
#     - 'random' mirrors BGGE style:
#         $<kernel>$u     : fitted random effect per observation (n)
#         $<kernel>$varu  : variance component estimate (if available)
#         $<kernel>$u.sd  : (BGGE) posterior sd of u; (BGLR) NA if not available
#         $<kernel>$varu.sd : (BGGE) posterior sd; (BGLR) NA if not available
# ================================================================

fit_kernel_models <- function(
  y,
  env,
  gid,
  uniqID,
  Data = NULL,
  random = NULL,
  fixed = NULL,
  engine = c("BGGE","BGLR"),
  addCol = NULL,
  verbose = FALSE,
  iterations = 1000,
  burnin = 200,
  thining = 10,
  tol = 1e-10,
  R2 = 0.5,
  digits = 4
){

  engine <- match.arg(engine)
  if (is.null(random) || !length(random)) stop("Missing 'random' kernels (from get_geno_kernel_v3).")

  Obs <- "Obs"
  isTest <- "isTest"
  
  print(addCol)


  # ----- Build Y (env, gid, y) in observation order -----
  if (is.character(y)) {
    if (is.null(Data)) stop("When 'y' is a column name, 'Data' must be provided.")
	
    if (!all(c(env, gid, y) %in% names(Data))) stop("env/gid/y not found in 'Data'.")
	
	if(!is.null(addCol) & any(colnames(Data) %in% addCol)){
			Y <- data.frame(env = Data[[env]], gid = Data[[gid]], y = Data[[y]], Obs=Data[[Obs]], uniqID=Data[[uniqID]], isTest=Data[[isTest]], addCol = Data[[addCol]], stringsAsFactors = FALSE)
	}else if(is.null(addCol)){ 
			Y <- data.frame(env = Data[[env]], gid = Data[[gid]], y = Data[[y]], Obs=Data[[Obs]], uniqID=Data[[uniqID]],isTest=Data[[isTest]], stringsAsFactors = FALSE)
	}
  }else {
    if (is.null(Data) || !all(c(env, gid) %in% names(Data))) {
      stop("When 'y' is a numeric vector, provide 'Data' with env and gid columns.")
    }
    if (length(y) != nrow(Data)) stop("Length of 'y' must match nrow(Data).")
    if(!is.null(addCol) & any(colnames(Data) %in% addCol)){   
      Y <- data.frame(env = Data[[env]], gid = Data[[gid]], y = as.numeric(y), Obs=Data[[Obs]], uniqID=Data[[uniqID]], isTest=Data[[isTest]],addCol= Data[[addCol]],stringsAsFactors = FALSE)
	}else{
	  Y <- data.frame(env = Data[[env]], gid = Data[[gid]], y = as.numeric(y), Obs=Data[[Obs]], uniqID=Data[[uniqID]], isTest=Data[[isTest]], stringsAsFactors = FALSE)
	}
  }
  
  
  Y$env <- factor(Y$env); Y$gid <- factor(Y$gid)
  n  <- nrow(Y)
  obs <- paste0(Y$env, ":", Y$gid)

  # ----- Align every kernel to current obs order (safeguard) -----
  align_to_obs <- function(M, obs_names) {
    if (is.null(rownames(M)) || is.null(colnames(M)))
      stop("All kernels must have row/col names matching 'paste(env,\":\",gid)'.")
    ix <- match(obs_names, rownames(M)); iy <- match(obs_names, colnames(M))
    if (anyNA(ix) || anyNA(iy)) stop("Kernel names do not match current obs order.")
    M[ix, iy, drop = FALSE]
  }
  
### Align kernel to observations 
  
  K_aligned <- lapply(random, function(k) {
    K <- align_to_obs(k$Kernel, obs)
    list(Kernel = K, Type = k$Type)
  })
  names(K_aligned) <- names(random)

# ----- Fixed effects matrix (optional) -----

  XF <- NULL
  if (!is.null(fixed)) {
    if (!is.matrix(fixed)) fixed <- as.matrix(fixed)
    if (nrow(fixed) != n) stop("fixed must have nrow equal to number of observations.")
    XF <- fixed
  }

  # ======================================================================
  # BGGE BACKEND (internal Gibbs; groups rows by environment for BD)
  # ======================================================================
  if (engine == "BGGE"){

    # ---- Reorder to make BD blocks contiguous by environment ----
    ord <- order(Y$env, Y$gid)
    Y_ord <- Y[ord, , drop = FALSE]
    y_ord <- Y_ord$y
    
	if (!is.null(XF)) XF_ord <- XF[ord, , drop = FALSE] else XF_ord <- NULL
    ne_vec <- as.vector(table(Y_ord$env))

    K_ord <- lapply(K_aligned, function(k) {
      Ko <- k$Kernel[ord, ord, drop = FALSE]
      list(Kernel = Ko, Type = k$Type)
    })
	
 ### BGGE Function (from the original package)

    BGGE <- function(y, K, XF = NULL, ne, ite = 1000, burn = 200, 
        thin = 3, verbose = FALSE, tol = 1e-10, R2 = 0.5) {
        
		dcondb <- function(n, media, vari) {
            sd <- sqrt(vari)
            return(rnorm(n, media, sd))
        }
        dcondsigb <- function(b, deltav, n, nu, Sc) {
            z <- sum(b^2 * deltav)
            return(1/rgamma(1, (n + nu)/2, (z + nu * Sc)/2))
        }
		
        dcondsigsq <- function(Aux, n, nu, Sce) {
            return(1/rgamma(1, (n + nu)/2, crossprod(Aux)/2 + 
                Sce/2))
        }
		
        rmvnor <- function(n, media, sigma) {
            z <- rnorm(n)
            return(media + crossprod(chol(sigma), z))
        }
		
        eig <- function(K, tol) {
            ei <- eigen(K)
            fil <- which(ei$values > tol)
            return(list(ei$values[fil], ei$vectors[, fil]))
        }
		
        setDEC <- function(K, tol, ne) {
            sK <- vector("list", length = length(K))
            typeM <- sapply(K, function(x) x$Type)
            if (!all(typeM %in% c("BD", "D"))) 
                stop("Matrix should be of types BD or D")
            if (missing(ne)) {
                if (any(typeM == "BD")) 
                  stop("For type BD, number of subjects in each sub matrices should be provided")
            }
            else {
                if (length(ne) <= 1 & any(typeM == "BD")) 
                  stop("ne invalid. For type BD, number of subjects in each sub matrices should be provided")
                nsubK <- length(ne)
                if (nsubK > 1) {
                  posf <- cumsum(ne)
                  posi <- cumsum(c(1, ne[-length(ne)]))
                }
            }
            for (i in 1:length(K)){
                if (K[[i]]$Type == "D") {
                  tmp <- list()
                  ei <- eig(K[[i]]$Kernel, tol)
                  tmp$s <- ei[[1]]
                  tmp$U <- ei[[2]]
                  tmp$tU <- t(ei[[2]])
                  tmp$nr <- length(ei[[1]])
                  tmp$type <- "D"
                  tmp$pos <- NA
                  sK[[i]] <- list(tmp)
                }
		        if (K[[i]]$Type == "BD") {
                  cont <- 0
                  tmp <- list()
                  for (j in 1:nsubK) {
                    Ktemp <- K[[i]]$Kernel[(posi[j]:posf[j]), 
                      (posi[j]:posf[j])]
                    ei <- eig(Ktemp, tol)
                    if (length(ei[[1]]) != 0) {
                      cont <- cont + 1
                      tmp[[cont]] <- list()
                      tmp[[cont]]$s <- ei[[1]]
                      tmp[[cont]]$U <- ei[[2]]
                      tmp[[cont]]$tU <- t(ei[[2]])
                      tmp[[cont]]$nr <- length(ei[[1]])
                      tmp[[cont]]$type <- "BD"
                      tmp[[cont]]$pos <- c(posi[j], posf[j])
                    }
                  }
                  if (length(tmp) > 1) {
                    sK[[i]] <- tmp
                  }
                  else {
                    sK[[i]] <- list(tmp[[1]])
                  }
                }
            }
            return(sK)
        }
        if (as.numeric(verbose) != 0) {
            cat("Setting parameters...", "\n", "\n")
        }
		
        y <- as.numeric(y)
        yNA <- is.na(y)
        whichNa <- which(yNA)
        nNa <- length(whichNa)
        mu <- mean(y, na.rm = TRUE)
        n <- length(y)
        if (nNa > 0) {
            y[whichNa] <- mu
        }
		
        if (is.null(names(K))) {
            names(K) <- paste("K", seq(length(K)), sep = "")
        }
		
        if (!is.null(XF)){
            Bet <- solve(crossprod(XF), crossprod(XF, y))
            nBet <- length(Bet)
            tXX <- solve(crossprod(XF))
        }
		
        nk <- length(K)
        nr <- numeric(length(K))
        typeM <- sapply(K, function(x) x$Type)
        if (!all(typeM %in% c("BD", "D"))) 
            stop("Matrix should be of types BD or D")
        if (length(ne) == 1 & any(typeM == "BD")) 
            stop("Type BD should be used only for block diagonal matrices")
			
        Ei <- setDEC(K = K, ne = ne, tol = tol)
        nCum <- sum(seq(1, ite)%%thin == 0)
        chain <- vector("list", length = 3)
        names(chain) <- c("mu", "varE", "K")
        chain$varE <- numeric(nCum)
        chain$mu <- numeric(nCum)
        chain$K <- vector("list", length = nk)
        names(chain$K) <- names(K)
        chain$K[seq(nk)] <- list(list(varU = numeric(nCum)))
        cpred <- vector("list", length = nk)
        names(cpred) <- names(K)
        cpred[seq(nk)] <- list(U = matrix(NA_integer_, nrow = nCum, 
            ncol = n))
			
	    
        nu <- 3
        Sce <- (nu + 2) * (1 - R2) * var(y, na.rm = TRUE)
        Sc <- numeric(length(K))
        for (i in 1:length(K)) {
            Sc[i] <- (nu + 2) * R2 * var(y, na.rm = T)/mean(diag(K[[i]]$Kernel))
        }
        tau <- 0.01
        u <- list()
        for (j in 1:nk) {
            u[[j]] <- rnorm(n, 0, 1/(2 * n))
        }
        sigsq <- var(y)
        sigb <- rep(0.2, nk)
        temp <- y - mu
        if (!is.null(XF)) {
            B.mcmc <- matrix(0, nrow = nCum, ncol = nBet)
            temp <- temp - XF %*% Bet
        }
        temp <- temp - Reduce("+", u)
        nSel <- 0
        i <- 1
		
        while (i <= ite) {
            time.init <- proc.time()[3]
            temp <- temp + mu
            mu <- rnorm(1, mean(temp), sqrt(sigsq/n))
            temp <- temp - mu
            if (!is.null(XF)) {
                temp <- temp + XF %*% Bet
                vari <- sigsq * tXX
                media <- tXX %*% crossprod(XF, temp)
                Bet <- rmvnor(nBet, media, vari)
                temp <- temp - XF %*% Bet
            }
            for (j in 1:nk) {
                if (typeM[j] == "D") {
                  temp <- temp + u[[j]]
                  d <- crossprod(Ei[[j]][[1]]$U, temp)
                  s <- Ei[[j]][[1]]$s
                  deltav <- 1/s
                  lambda <- sigb[j]
                  vari <- s * lambda/(1 + s * lambda * tau)
                  media <- tau * vari * d
                  nr <- Ei[[j]][[1]]$nr
                  b <- dcondb(nr, media, vari)
                  u[[j]] <- crossprod(Ei[[j]][[1]]$tU, b)
                  temp <- temp - u[[j]]
                }
                if (typeM[j] == "BD") {
                  nsk <- length(Ei[[j]])
                  if (length(nsk > 1)) {
                    temp <- temp + u[[j]]
                    d <- NULL
                    s <- NULL
                    neiv <- numeric(nsk)
                    pos <- matrix(NA, ncol = 2, nrow = nsk)
                    for (k in 1:nsk) {
                      pos[k, ] <- Ei[[j]][[k]]$pos
                      d <- c(d, crossprod(Ei[[j]][[k]]$U, temp[pos[k, 
                        1]:pos[k, 2]]))
                      neiv[k] <- length(Ei[[j]][[k]]$s)
                      s <- c(s, Ei[[j]][[k]]$s)
                    }
                    deltav <- 1/s
                    lambda <- sigb[j]
                    vari <- s * lambda/(1 + s * lambda * tau)
                    media <- tau * vari * d
                    nr <- length(s)
                    b <- dcondb(nr, media, vari)
                    posf <- cumsum(neiv)
                    posi <- cumsum(c(1, neiv[-length(neiv)]))
                    utmp <- numeric(n)
                    for (k in 1:nsk) {
                      utmp[pos[k, 1]:pos[k, 2]] <- crossprod(Ei[[j]][[k]]$tU,b[posi[k]:posf[k]])
                    }
                    u[[j]] <- utmp
                    temp <- temp - u[[j]]
                  }
                  else {
                    temp <- temp + u[[j]]
                    pos <- Ei[[j]]$pos
                    d <- crossprod(Ei[[j]][[1]]$U, temp[pos[1]:pos[2]])
                    s <- Ei[[j]][[1]]$s
                    deltav <- 1/s
                    lambda <- sigb[j]
                    vari <- s * lambda/(1 + s * lambda * tau)
                    media <- tau * vari * d
                    nr <- Ei[[j]][[1]]$nr
                    b <- dcondb(nr, media, vari)
                    utmp <- numeric(n)
                    utmp[pos[1]:pos[2]] <- crossprod(Ei[[j]][[1]]$tU, 
                      b)
                    u[[j]] <- utmp
                    temp <- temp - u[[j]]
                  }
                }
                sigb[j] <- dcondsigb(b, deltav, nr, nu, Sc[j])
            }
            res <- temp
            sigsq <- dcondsigsq(res, n, nu, Sce)
            tau <- 1/sigsq
            if (nNa > 0){
                uhat <- Reduce("+", u)
                if (!is.null(XF)) {
                  aux <- XF[yNA, ] %*% Bet
                }
                else {
                  aux <- 0
                }
                y[whichNa] <- aux + mu + uhat[whichNa] + rnorm(n = nNa, 
                  sd = sqrt(sigsq))
                temp[whichNa] <- y[whichNa] - uhat[whichNa] - 
                  aux - mu
            }
            if (i%%thin == 0) {
                nSel <- nSel + 1
                chain$varE[nSel] <- sigsq
                chain$mu[nSel] <- mu
                if (!is.null(XF)) {
                  B.mcmc[nSel, ] <- Bet
                }
                for (j in 1:nk) {
                  cpred[[j]][nSel, ] <- u[[j]]
                  chain$K[[j]]$varU[nSel] <- sigb[j]
                }
            }
            if (as.numeric(verbose) != 0 & i%%as.numeric(verbose) == 
                0) {
                time.end <- proc.time()[3]
                cat("Iter: ", i, "time: ", round(time.end - time.init, 
                  3), "\n")
            }
            i <- i + 1
        }
        
		draw <- seq(ite)[seq(ite)%%thin == 0] > burn
        mu.est <- mean(chain$mu[draw])
        yHat <- mu.est
        
		if (!is.null(XF)) {
            B <- colMeans(B.mcmc[draw, ])
            yHat <- yHat + XF %*% B
        }
        u.est <- sapply(cpred, FUN = function(x) colMeans(x[draw, 
            ]))
        yHat <- yHat + rowSums(u.est)
        out <- list()
        out$yHat <- yHat
        out$varE <- mean(chain$varE[draw])
        out$varE.sd <- sd(chain$varE[draw])
        out$K <- vector("list", length = nk)
        names(out$K) <- names(cpred)
        for (i in 1:nk) {
            out$K[[i]]$u <- colMeans(cpred[[i]][draw, ])
            out$K[[i]]$u.sd <- apply(cpred[[i]][draw, ], MARGIN = 2, 
                sd)
            out$K[[i]]$varu <- mean(chain$K[[i]]$varU[draw])
            out$K[[i]]$varu.sd <- sd(chain$K[[i]]$varU[draw])
        }
        out$chain <- chain
        out$ite <- ite
        out$burn <- burn
        out$thin <- thin
        out$model <- K$model
        out$kernel <- K$kernel
        out$y <- y
        class(out) <- "BGGE"
        return(out)
    }
   
    
   Vcomp.BGGE <- function(model, env, gid, digits = digits,alfa = 0.1){
        t = length(unique(gid))
        e = length(unique(env))
        n = t * e
        K = model$K
        size = length(K)
        comps = data.frame(matrix(NA, ncol = 3, nrow = size))
        VarE = data.frame(matrix(NA, ncol = 3, nrow = 1))
        names(comps) = names(VarE) = c("K", "Var", "SD.var")
        for (k in 1:size) {
            comps[k, 1] = names(K)[k]
            comps[k, 2] = round(K[[k]]$varu, digits)
            comps[k, 3] = round(K[[k]]$varu.sd, digits)
        }
        VarE[1, 1] = "Residual"
        VarE[1, 2] = round(model$varE, digits)
        VarE[1, 3] = round(model$varE.sd, digits)
        comps = rbind(comps, VarE)
        comps$Type <- comps$K
        comps$Type[grep(comps$K, pattern = "KGE_")] = "GxE"
        comps$Type[grep(comps$K, pattern = "KG_")] = "Genotype (G)"
        # Override: KG_GE* kernels are GxE interactions, not main G
        # KG_GE -> "GxE", KG_GEbw1 -> "GxE_bw1", KG_GE_Env1 -> "GxE_Env1"
        gxe_idx <- grep(comps$K, pattern = "KG_GE")
        if (length(gxe_idx)) {
          sfx <- sub("KG_GE_?", "", comps$K[gxe_idx])
          comps$Type[gxe_idx] <- ifelse(nchar(sfx) == 0, "GxE", paste0("GxE_", sfx))
        }
        comps$Type[grep(comps$K, pattern = "^E$")] = "GxE"
        comps$Type[grep(comps$K, pattern = "KE_")] = "Environment (E)"
        comps$CI_upper = NA
        comps$CI_lower = NA
        ENV = which(comps$Type %in% "Environment (E)")
        GID = which(comps$Type %in% "Genotype (G)")
        GE = which(comps$Type %in% "GxE")
        R = which(comps$Type %in% "Residual")
		
		#### Upper CI
        comps$CI_upper[ENV] = (n - e) * comps$Var[ENV]/qchisq((alfa/2), 
            n - e)
        comps$CI_upper[GID] = (n - t) * comps$Var[GID]/qchisq((alfa/2), 
            n - t)
        comps$CI_upper[GE] = (n - t - e) * comps$Var[GE]/qchisq((alfa/2), 
            n - t - e)
        comps$CI_upper[R] = (n - t - e) * comps$Var[R]/qchisq((alfa/2), 
            n - t - e)
		
		#### Lower CI
        comps$CI_lower[ENV] = (n - e) * comps$Var[ENV]/qchisq((1 - 
            alfa/2), n - e)
        comps$CI_lower[GID] = (n - t) * comps$Var[GID]/qchisq((1 - 
            alfa/2), n - t)
        comps$CI_lower[GE] = (n - t - e) * comps$Var[GE]/qchisq((1 - 
            alfa/2), n - t - e)
        comps$CI_lower[R] = (n - t - e) * comps$Var[R]/qchisq((1 - 
            alfa/2), n - t - e)
			
		
        comps$CI_upper = round(comps$CI_upper, digits)
        comps$CI_lower = round(comps$CI_lower, digits)
        comps <- comps[, c(4, 1:2, 6, 5, 3)]
        
		p = c("KG_", "KE_", "KGE_")
        for (i in 1:3) comps$K = gsub(x = comps$K, pattern = p[i], 
            replacement = "")
        return(comps)
    }
  
  
    # ---- Fit BGGE on env-ordered data ----
    start <- Sys.time()
    fit_bgge <- BGGE(y = y_ord, K = K_ord, XF = XF_ord, ne = ne_vec,
                     ite = iterations, burn = burnin, thin = thining,
                     verbose = verbose, tol = tol, R2 = R2)
    end <- Sys.time()
    cat("--------------------------------------------------------\n")
    cat("Running BGGE (Bayesian Genotype + GxE)\n")
    cat("More Detail in Granato et al. (2018) G3\n")
    cat("--------------------------------------------------------\n")
    cat(paste0("Start at: ", start, " | Ended at: ", end, "\n"))

    # Re-map predictions to original order
    #yHat <- numeric(n); #yHat[ord] <- fit_bgge$yHat
	
	vComp <- Vcomp.BGGE(model = fit_bgge, env = Y$env, gid = Y$gid, digits = digits)

    ret <- list(
      yHat   = fit_bgge$yHat,
      varE   = fit_bgge$varE,
      random = fit_bgge$K,      # already contains u, u.sd, varu, varu.sd
      fit    = fit_bgge,
      VarComp= vComp, 
	  dat = Y_ord,	  
      engine = "BGGE"
    )
    return(ret)
  }

  # ======================================================================
  # BGLR BACKEND (RKHS using provided kernels)
  # ======================================================================
  if (engine == "BGLR") {
    if (!requireNamespace("BGLR", quietly = TRUE)) {
      stop("Package 'BGLR' is required for engine='BGLR'. Please install.packages('BGLR').")
    }

    # Build ETA: one RKHS term per kernel; add fixed as model='FIXED' if provided
    ETA <- list()
    if (!is.null(XF)) {
      ETA[["FIXED"]] <- list(X = XF, model = "FIXED")
    }
    for (nm in names(K_aligned)) {
      ETA[[nm]] <- list(K = K_aligned[[nm]]$Kernel, model = "RKHS")
    }

    start <- Sys.time()
    fit_bglr <- BGLR::BGLR(
      y       = Y$y,
      ETA     = ETA,
      nIter   = iterations,
      burnIn  = burnin,
      thin    = thining,
      verbose = verbose
    )
    end <- Sys.time()
    cat("--------------------------------------------------------\n")
    cat("Running BGLR (multi-kernel RKHS)\n")
    cat("--------------------------------------------------------\n")
    cat(paste0("Start at: ", start, " | Ended at: ", end, "\n"))

    # Compose per-kernel outputs similar to BGGE
    random_out <- list()
    for (nm in names(K_aligned)) {
      u_vec <- fit_bglr$ETA[[nm]]$u
      # varU may be a scalar or vector depending on BGLR settings; guard
      varu  <- tryCatch(fit_bglr$ETA[[nm]]$varU, error = function(e) NA_real_)
      random_out[[nm]] <- list(
        u        = as.numeric(u_vec),
        u.sd     = rep(NA_real_, length(u_vec)),
        varu     = if (length(varu)) as.numeric(varu[length(varu)]) else NA_real_,
        varu.sd  = NA_real_
      )
    }

    # yHat: if BGLR supplies fitted values use that; else reconstruct
    yHat <- fit_bglr$yHat
    if (is.null(yHat) || anyNA(yHat)) {
      mu  <- if (!is.null(fit_bglr$mu)) mean(fit_bglr$mu) else 0
      yHat <- rep(mu, n)
      if (!is.null(XF)) {
        # BGLR does not always store beta for FIXED; if absent, add the fitted contribution from X %*% b
        if (!is.null(fit_bglr$ETA$FIXED$beta)) {
          yHat <- yHat + as.numeric(XF %*% fit_bglr$ETA$FIXED$beta)
        }
      }
      for (nm in names(K_aligned)) yHat <- yHat + as.numeric(fit_bglr$ETA[[nm]]$u)
    }

    # Residual variance (if available)
    varE <- tryCatch({
      ve <- fit_bglr$varE
      if (length(ve)) as.numeric(ve[length(ve)]) else NA_real_
    }, error = function(e) NA_real_)

    # Simple VarComp table (no CIs from BGLR without extra bookkeeping)
    VarComp <- data.frame(
      Type    = ifelse(grepl("^KE_", names(K_aligned)), "Environment (E)",
                       ifelse(grepl("^KG_", names(K_aligned)), "Genotype (G)", "GxE")),
      K       = names(K_aligned),
      Var     = sapply(names(K_aligned), function(nm) {
        v <- tryCatch(fit_bglr$ETA[[nm]]$varU, error=function(e) NA_real_)
        if (length(v)) as.numeric(v[length(v)]) else NA_real_
      }),
      CI_upper= NA_real_,
      CI_lower= NA_real_,
      SD.var  = NA_real_,
      row.names = NULL
    )
    VarComp <- rbind(VarComp, data.frame(
      Type="Residual", K="Residual", Var=varE, CI_upper=NA, CI_lower=NA, SD.var=NA
    ))

    ret <- list(
      yHat   = as.numeric(yHat),
      varE   = varE,
      random = random_out,
      fit    = fit_bglr,
      VarComp= VarComp,
      engine = "BGLR"
    )
    return(ret)
  }

}



###
### fitMEModels_CV
### Legacy CV fitting function (original EnvRtype-based approach)
###

fitMEModels_CV <- function(DT,genoDat,strainGeno,KG=NULL,KE=NULL,CVMet,k,nIter,method,fitEnvModels=FALSE,FixedTerm=Fixed,IDColsList=IDColsMEList){

### Prepare Data for CV


  DT_2A <- DT
  DT_2B <- DT

  dim(DT_2B)
  DT_2B <- droplevels(DT_2B)


## CV1, where novel genotypes in tested environments are predicted.
## CV2, where tested genotypes in tested environments are predicted.
## CV0, where tested genotypes in untested novel environments are predicted.
## CV00, where novel genotypes in novel environments are predicted.
## CV LOFO (Leave One Factor Out), eg: Leave One Test Out/ Leave One Line Out cross validation.



 Dat_Out_List <- getTstIndices_CV_Data(DT_2B,CVMet,nIter,k)

##### Here IDCols is the vector with 'gid', 'env' and 'value'

  IDCols <- IDColsList$`IDColsMod`

 ### IDColsMod is the vector with the original IDs  for output
  IDColsMod <- IDColsList$`IDCols`

 ###

  if(!is.null(KG)){
    method <- NULL
  }

  Fixed <- FixedTerm


   cor_CV_List_Reps <- list()
   var_CV_List_Reps <- list()
   fit_Out_CV_List_Reps <- list()

  nReps <- length(Dat_Out_List[[1]])

  for(nrep in 1:nReps){


	cor_CV_List_nF <- list()
	var_CV_List_nF <- list()
	fit_Out_CV_List_nF <- list()

	k <- length(Dat_Out_List[[1]][[nrep]])

   for(nF in 1:k){

    DT_2B <- Dat_Out_List[[1]][[nrep]][[nF]]
	tstIndices2 <-  Dat_Out_List[[3]][[nrep]][[nF]]
    DT_tstIndices2 <- tstIndices2

    dim(DT_2B)
    DT_2B <- droplevels(DT_2B)


    Y <- DT_2B[,c("env","gid","value")]

    y <- "value"
    gid <- "gid"
    env <- "env"

    X <- genoDat


	if(fitEnvModels==FALSE & is.null(KE) & !is.null(KG) & is.null(method)){

		## Creating kernel models with EnvRtype::get_kernel

		system.time({
		  MM = EnvRtype::get_kernel(K_G = KG, y=y, gid = gid, env = env, data = DT_2B, model = "MM")
		})

		system.time({
		  MDs = EnvRtype::get_kernel(K_G = KG, y=y, gid = gid, env = env, data = DT_2B, model = "MDs")
		})

		### Fit heterogeneous variance using BGGE

		MDe <- get_kernel_MDe(Y,KG,KE,intercept.random=FALSE,dimension_KE=NULL)

		fixed= model.matrix(~0+env,DT_2B)

		system.time({
		  fit_MM= EnvRtype::kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=MM, fixed=fixed)
		})


		system.time({
		  fit_MDs= EnvRtype::kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=MDs, fixed=fixed)
		})


		y_MDe <- DT_2B[,"value"]
		NE <- table(DT_2B$env)

		system.time({
		  fit_MDe <- BGGE(y_MDe,K=MDe,XF=fixed,ne=NE)
		})

		digits <- 3
		VarComp = Vcomp.BGGE(model = fit_MDe, env = Y$env,gid = Y$gid, digits = digits)

		fit_MDe_Out = list(yHat = fit_MDe$yHat, varE = fit_MDe$varE, random = fit_MDe$K,
						   BGGE = fit_MDe, VarComp = VarComp)

		corMM_CV <- cor(fit_MM$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")
		corMDs_CV <- cor(fit_MDs$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")
		corMDe_CV <- cor(fit_MDe_Out$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")

		varMM_CV <- fit_MM$VarComp
		varMDs_CV <- fit_MDs$VarComp
		varMDe_CV <- fit_MDe_Out$VarComp

	  }else if(fitEnvModels ==FALSE & is.null(KE) & is.null(KG) & !is.null(method)){
		if(method=="GB"){
		  ## Creating kernel models with EnvRtype::get_kernel

		  system.time({
			MM_GB = get_geno_kernel(Y,X,KG,KE,intercept.random=FALSE,model="MM",method="GB",dimension_KE=NULL)
		  })

		  system.time({
			MDs_GB =get_geno_kernel(Y,X,KG,KE,intercept.random=FALSE,model="MDs",method="GB",dimension_KE=NULL)
		  })


		  system.time({
			MDe_GB <- get_geno_kernel(Y,X,KG,KE,intercept.random=FALSE,model="MDe",method="GB",dimension_KE=NULL)

		  })

		  fixed=model.matrix(~0+env,DT_2B)

		  system.time({
			fit_MM_GB=EnvRtype::kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=MM_GB, fixed=fixed)
		  })


		  system.time({
			fit_MDs_GB=EnvRtype::kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=MDs_GB, fixed=fixed)
		  })


		  ###Fit heterogeneous variance using BGGE

		  y_MDe <- DT_2B[,"value"]
		  NE <- table(DT_2B$env)
		  system.time({
			fit_MDe_GB <- BGGE::BGGE(y_MDe,K=MDe_GB,XF=fixed,ne=NE)
		  })

		  digits <- 3
		  VarComp = Vcomp.BGGE(model = fit_MDe_GB, env = Y$env,
							   gid = Y$gid, digits = digits)


		  fit_MDe_GB_Out = list(yHat = fit_MDe_GB$yHat, varE = fit_MDe_GB$varE, random = fit_MDe_GB$K,
								BGGE = fit_MDe_GB, VarComp = VarComp)


		  corMM_CV <- cor(fit_MM_GB$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")
		  corMDs_CV <- cor(fit_MDs_GB$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")
		  corMDe_CV <- cor(fit_MDe_GB_Out$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")

		  varMM_CV <- fit_MM_GB$VarComp
		  varMDs_CV <- fit_MDs_GB$VarComp
		  varMDe_CV <- fit_MDe_GB_Out$VarComp

		  fit_MM_Out <- cbind.data.frame(DT_2A[DT_tstIndices2,IDCols],DT_2A[DT_tstIndices2,"value"],fit_MM_GB$yHat[tstIndices2])
		  fit_MDs_Out <- cbind.data.frame(DT_2A[DT_tstIndices2,IDCols],DT_2A[DT_tstIndices2,"value"],fit_MDs_GB$yHat[tstIndices2])
		  fit_MDe_Out <- cbind.data.frame(DT_2A[DT_tstIndices2,IDCols],DT_2A[DT_tstIndices2,"value"],fit_MDe_GB_Out$yHat[tstIndices2])

		  colnames(fit_MM_Out) <- c(IDColsMod,"Obs","Pred")
		  colnames(fit_MDs_Out) <- c(IDColsMod,"Obs","Pred")
		  colnames(fit_MDe_Out) <- c(IDColsMod,"Obs","Pred")

    }
	if(method=="GK"){
      ## Creating kernel models with EnvRtype::get_kernel

		  system.time({
			MM_GK = get_geno_kernel(Y,X,KG,KE,intercept.random=FALSE,model="MM",method="GK",dimension_KE=NULL)
		  })

		  system.time({
			MDs_GK = get_geno_kernel(Y,X,KG,KE,intercept.random=FALSE,model="MDs",method="GK",dimension_KE=NULL)
		  })

		  ### Fit heterogeneous variance using BGGE
		  system.time({
			MDe_GK <- get_geno_kernel(Y,X,KG,KE,intercept.random=FALSE,model="MDe",method="GK",dimension_KE=NULL)

		  })

####

		  fixed=model.matrix(~0+env,DT_2B)

		  system.time({
			fit_MM_GK=EnvRtype::kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=MM_GK, fixed=fixed)
		  })


		  system.time({
			fit_MDs_GK=EnvRtype::kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=MDs_GK, fixed=fixed)
		  })

		  ### MDe

		  y_MDe <- DT_2B[,"value"]
		  NE <- table(DT_2B$env)
		  system.time({
			fit_MDe_GK <- BGGE::BGGE(y_MDe,K=MDe_GK,XF=fixed,ne=NE)
		  })

		  digits <- 3
		  VarComp = Vcomp.BGGE(model = fit_MDe_GK, env = Y$env,
							   gid = Y$gid, digits = digits)


		  fit_MDe_GK_Out = list(yHat = fit_MDe_GK$yHat, varE = fit_MDe_GK$varE, random = fit_MDe_GK$K,
								BGGE = fit_MDe_GK, VarComp = VarComp)

		  ###

		  corMM_CV <- cor(fit_MM_GK$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")
		  corMDs_CV <- cor(fit_MDs_GK$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")
		  corMDe_CV <- cor(fit_MDe_GK_Out$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")

		  varMM_CV <- fit_MM_GK$VarComp
		  varMDs_CV <- fit_MDs_GK$VarComp
		  varMDe_CV <- fit_MDe_GK_Out$VarComp



		  fit_MM_Out <- cbind.data.frame(DT_2A[DT_tstIndices2,IDCols],DT_2A[DT_tstIndices2,"value"],fit_MM_GK$yHat[tstIndices2])
		  fit_MDs_Out <- cbind.data.frame(DT_2A[DT_tstIndices2,IDCols],DT_2A[DT_tstIndices2,"value"],fit_MDs_GK$yHat[tstIndices2])
		  fit_MDe_Out <- cbind.data.frame(DT_2A[DT_tstIndices2,IDCols],DT_2A[DT_tstIndices2,"value"],fit_MDe_GK_Out$yHat[tstIndices2])

		  colnames(fit_MM_Out) <- c(IDColsMod,"Obs","Pred")
		  colnames(fit_MDs_Out) <- c(IDColsMod,"Obs","Pred")
		  colnames(fit_MDe_Out) <- c(IDColsMod,"Obs","Pred")

    }
  }else if(fitEnvModels ==TRUE & !is.null(KE) & !is.null(KG) & is.null(method)){

    EMM = EnvRtype::get_kernel(K_G = KG, K_E = KE, y= y, gid = gid, env = env, data = DT_2B,model = "EMM",ne=nrow(KE))

    EMDs = EnvRtype::get_kernel(K_G = KG, K_E = KE, y=y, gid = gid, env = env, data = DT_2B,model = "EMDs")

    EMDe <- get_kernel_MDe(Y,KG,KE,intercept.random=FALSE,dimension_KE=NULL)

    fixed=model.matrix(~0+env,DT_2B)

    system.time({
      fit_EMM =EnvRtype::kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=EMM, fixed=fixed)
    })

    system.time({
      fit_EMDs=EnvRtype::kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=EMDs, fixed=fixed)
    })

    y <- DT_2B[,"value"]
    NE <- table(DT_2B$env)

    system.time({
      fit_EMDe <- BGGE::BGGE(y,K=EMDe,XF=fixed,ne=NE)
    })

    digits <-3
    VarComp = Vcomp.BGGE(model = fit_EMDe, env = Y$env,
                         gid = Y$gid, digits = digits)

    fit_EMDe_Out = list(yHat = fit_EMDe$yHat, varE = fit_EMDe$varE, random = fit_EMDe$K,
                        BGGE = fit_EMDe, VarComp = VarComp)


    corMM_CV <- cor(fit_EMM$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")
    corMDs_CV <- cor(fit_EMDs$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")
    corMDe_CV <- cor(fit_EMDe$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")

	fit_MM_Out <- fit_EMM
	fit_MDs_Out <- fit_EMDs
	fit_MDe_Out <- fit_EMDe_Out

    varMM_CV <- fit_EMM$VarComp
    varMDs_CV <- fit_EMDs$VarComp
    varMDe_CV <- fit_EMDe_Out$VarComp

	fit_MM_Out <- cbind.data.frame(DT_2A[DT_tstIndices2,IDCols],DT_2A[DT_tstIndices2,"value"],fit_EMM$yHat[tstIndices2])
    fit_MDs_Out <- cbind.data.frame(DT_2A[DT_tstIndices2,IDCols],DT_2A[DT_tstIndices2,"value"],fit_EMDs$yHat[tstIndices2])
    fit_MDe_Out <- cbind.data.frame(DT_2A[DT_tstIndices2,IDCols],DT_2A[DT_tstIndices2,"value"],fit_EMDe_Out$yHat[tstIndices2])

    colnames(fit_MM_Out) <- c(IDColsMod,"Obs","Pred")
    colnames(fit_MDs_Out) <- c(IDColsMod,"Obs","Pred")
    colnames(fit_MDe_Out) <- c(IDColsMod,"Obs","Pred")

  }else if(fitEnvModels ==TRUE & !is.null(KE) & is.null(KG) & !is.null(method)){

    if(method=="GK"){
      system.time({
        EMM_GK = get_geno_kernel(Y,X,KG,KE,intercept.random=FALSE,model="MM",method="GK",dimension_KE=NULL)
      })
      system.time({
        EMDs_GK = get_geno_kernel(Y,X,KG,KE,intercept.random=FALSE,model="MDs",method="GK",dimension_KE=NULL)
      })
      system.time({
        EMDe_GK <- get_geno_kernel(Y,X,KG,KE,intercept.random=FALSE,model="MDe",method="GK",dimension_KE=NULL)
      })

      fixed= model.matrix(~0+env,DT_2B)

      system.time({
        fit_EMM_GK =EnvRtype::kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=EMM_GK, fixed=fixed)
      })

      system.time({
        fit_EMDs_GK =EnvRtype::kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=EMDs_GK, fixed=fixed)
      })

      y <- DT_2B[,"value"]
      NE <- table(DT_2B$env)

      system.time({
        fit_EMDe_GK <- BGGE::BGGE(y,K=EMDe_GK,XF=fixed,ne=NE)
      })

      digits <-3
      VarComp = Vcomp.BGGE(model = fit_EMDe_GK, env = Y$env,
                           gid = Y$gid, digits = digits)

      fit_EMDe_GK_Out = list(yHat = fit_EMDe_GK$yHat, varE = fit_EMDe_GK$varE, random = fit_EMDe_GK$K,
                             BGGE = fit_EMDe_GK, VarComp = VarComp)


      corMM_CV <- cor(fit_EMM$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")
      corMDs_CV <- cor(fit_EMDs$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")
      corMDe_CV <- cor(fit_EMDe$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")

	  fit_MM_Out <- fit_EMM
	  fit_MDs_Out <- fit_EMDs
	  fit_MDe_Out <- fit_EMDe_Out

      varMM_CV <- fit_EMM$VarComp
      varMDs_CV <- fit_EMDs$VarComp
      varMDe_CV <- fit_EMDe_Out$VarComp

	  fit_MM_Out <- cbind.data.frame(DT_2A[DT_tstIndices2,IDCols],DT_2A[DT_tstIndices2,"value"],fit_EMM$yHat[tstIndices2])
      fit_MDs_Out <- cbind.data.frame(DT_2A[DT_tstIndices2,IDCols],DT_2A[DT_tstIndices2,"value"],fit_EMDs$yHat[tstIndices2])
      fit_MDe_Out <- cbind.data.frame(DT_2A[DT_tstIndices2,IDCols],DT_2A[DT_tstIndices2,"value"],fit_EMDe_Out$yHat[tstIndices2])


      colnames(fit_MM_Out) <- c(IDColsMod,"Obs","Pred")
      colnames(fit_MDs_Out) <- c(IDColsMod,"Obs","Pred")
      colnames(fit_MDe_Out) <- c(IDColsMod,"Obs","Pred")


    }
	if(method=="GB"){
      system.time({
        EMM_GB = get_geno_kernel(Y,X,KG,KE,intercept.random=FALSE,model="MM",method="GB",dimension_KE=NULL)
      })

      system.time({
        EMDs_GB =get_geno_kernel(Y,X,KG,KE,intercept.random=FALSE,model="MDs",method="GB",dimension_KE=NULL)
      })

      ### Fit heterogeneous variance using BGGE
      system.time({
        EMDe_GB <- get_geno_kernel(Y,X,KG,KE,intercept.random=FALSE,model="MDe",method="GB",dimension_KE=NULL)

      })


      fixed=model.matrix(~0+env,DT_2B)

      system.time({
        fit_EMM_GB =EnvRtype::kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=EMM_GB, fixed=fixed)
      })

      system.time({
        fit_EMDs_GB =EnvRtype::kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=EMDs_GB, fixed=fixed)
      })

      y <- DT_2B[,"value"]
      NE <- table(DT_2B$env)

      system.time({
        fit_EMDe_GB <- BGGE::BGGE(y,K=EMDe_GB,XF=fixed,ne=NE)
      })

      digits <-3
      VarComp = Vcomp.BGGE(model = fit_EMDe_GB, env = Y$env,
                           gid = Y$gid, digits = digits)

      fit_EMDe_GB_Out = list(yHat = fit_EMDe_GB$yHat, varE = fit_EMDe_GB$varE, random = fit_EMDe_GB$K,
                             BGGE = fit_EMDe_GB, VarComp = VarComp)



      corMM_CV <- cor(fit_EMM$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")
      corMDs_CV <- cor(fit_EMDs$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")
      corMDe_CV <- cor(fit_EMDe$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"],use="pairwise.complete.obs")

	  fit_MM_Out <- fit_EMM
	  fit_MDs_Out <- fit_EMDs
	  fit_MDe_Out <- fit_EMDe_Out

      varMM_CV <- fit_EMM$VarComp
      varMDs_CV <- fit_EMDs$VarComp
      varMDe_CV <- fit_EMDe_Out$VarComp

	  fit_MM_Out <- cbind.data.frame(DT_2A[DT_tstIndices2,IDCols],DT_2A[DT_tstIndices2,"value"],fit_EMM$yHat[tstIndices2])
      fit_MDs_Out <- cbind.data.frame(DT_2A[DT_tstIndices2,IDCols],DT_2A[DT_tstIndices2,"value"],fit_EMDs$yHat[tstIndices2])
      fit_MDe_Out <- cbind.data.frame(DT_2A[DT_tstIndices2,IDCols],DT_2A[DT_tstIndices2,"value"],fit_EMDe_Out$yHat[tstIndices2])



      colnames(fit_MM_Out) <- c(IDColsMod,"Obs","Pred")
      colnames(fit_MDs_Out) <- c(IDColsMod,"Obs","Pred")
      colnames(fit_MDe_Out) <- c(IDColsMod,"Obs","Pred")

   }
  }

	  cor_CV_List <- list(corMM_CV,corMDs_CV,corMDe_CV)
	  var_CV_List <- list(varMM_CV,varMDs_CV,varMDe_CV)
	  fit_Out_CV_List <- list(fit_MM_Out,fit_MDs_Out,fit_MDe_Out)


	  cor_CV_List_nF[[nF]] <- cor_CV_List
	  var_CV_List_nF[[nF]] <- var_CV_List
	  fit_Out_CV_List_nF[[nF]] <- fit_Out_CV_List
  }

   cor_CV_List_Reps[[nrep]] <- cor_CV_List_nF
   var_CV_List_Reps[[nrep]] <- var_CV_List_nF
   fit_Out_CV_List_Reps[[nrep]] <- fit_Out_CV_List_nF

  }


  return(list(cor_CV_List_Reps,var_CV_List_Reps,fit_Out_CV_List_Reps))

}
