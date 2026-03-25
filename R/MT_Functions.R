
#' Perform cross validation for multi-trait genomic prediction in CV1 and CV2 schemes #'
#'
#'
#' @param Data_Table_Num_Filt_List List. A list containing the data tables.
#' @param trait Character. Character specifying the traits to be included in the model.
#' @param nTraits Ineger. Integer specifying the number of traits.
#' @param k Integer. Integer specifying the number of folds in k-fold CV
#' @param nIter Integer. Integer value specifying the number of iterations.
#' @param CVMet Character. Character specifying the cross validation scheme with the default value set to "CV1".
#' Options include "CV1" and "CV2"
#'
#' @return a data frame that contains output stats such as prediction ability for the mlti-trait models .
#' @examples
#' \dontrun{
#' # Example usage of getMTCVR
#' result <- getMTCVR(...)
#' }
#' @export
getMTCVR <- function(Data_Table_Num_Filt_List, trait, nTraits, k, nIter, CVMet,
                     b_iter=5000, b_burnin=1000, b_thin=10, b_R2=0.5, b_digits=4,
                     A_ext = NULL) {

  Train <- Data_Table_Num_Filt_List[[1]]

  # --- data extraction ---
  PH <- as.matrix( apply(Train[, trait, drop = FALSE], 2, function(x) as.numeric(as.character(x))) )
  if (is.null(rownames(Train))) stop("Row names (genotype IDs) are required.")
  ids <- rownames(Train)

  # select marker cols
  genoInd <- grep("^ss", colnames(Train))
  if (length(genoInd) == 0) stop("No marker columns found with pattern '^ss\\d+$'.")
  Xraw <- as.matrix( apply(Train[, genoInd, drop = FALSE], 2, as.numeric) )
  Xraw <- Xraw - 1  # map  0/1/2 -> -1,0,1 (assumes your coding)

  # remove any rows with NA in any target trait
  keep <- stats::complete.cases(PH)
  PH   <- PH[keep, , drop = FALSE]
  Xraw <- Xraw[keep, , drop = FALSE]
  ids  <- ids[keep]

  n <- nrow(Xraw)
  if (n < k) stop("k-fold CV requires n >= k.")

  # global imputation (NOTE: for strict CV, move this inside each fold using only train stats)
  if (anyNA(Xraw)) {
    Ximp <- snpQC(Xraw, impute = TRUE, remove = FALSE, MAF = 0.02)
  } else {
    Ximp <- Xraw
  }
  rownames(Ximp) <- ids
  rownames(PH)   <- ids

  # precompute kernels where applicable
  if (!is.null(A_ext)) {
    common <- intersect(ids, rownames(A_ext))
    A_tot  <- A_ext[common, common]
  } else {
    A_tot <- rrBLUP::A.mat(Ximp)
  }

  # CV fold generator (balanced; uses all samples)
  make_folds <- function(n, k, seed) {
    set.seed(seed)
    idx <- sample.int(n)
    split(idx, rep_len(1:k, n))
  }

  # prediction extractor
  get_yhat <- function(fm) {
    if (!is.null(fm$yHat)) return(fm$yHat)
    if (!is.null(fm$YHat)) return(fm$YHat)
    if (!is.null(fm$ETAHat)) return(fm$ETAHat)
    # last resort: sum ETA + mu
    if (!is.null(fm$ETA) && !is.null(fm$mu)) {
      # try to assemble fitted values
      parts <- lapply(fm$ETA, function(e) if (!is.null(e$u)) e$u else if (!is.null(e$b)) Ximp %*% e$b else 0)
      fit <- Reduce(`+`, parts)
      fit <- sweep(fit, 2, fm$mu, `+`)
      return(fit)
    }
    stop("Could not find predictions (yHat/YHat/ETAHat) in BGLR Multitrait fit.")
  }

  # correlation helper that works for 1 or many traits
  colPA <- function(pred, obs) {
    if (is.null(dim(pred))) return(stats::cor(pred, obs, use = "pairwise.complete.obs"))
    diag(stats::cor(pred, obs, use = "pairwise.complete.obs"))
  }

  GPModelMT_List <- c("BRR (BGLR)", "RKHS (BGLR)", "Spike-Slab (BGLR)")

  # Pre-generate all fold sets so workers receive deterministic splits
  fold_sets <- lapply(seq_len(nIter), function(r) make_folds(n, k, seed = 125 + r))

  # Build flat task index: (model, rep, trait, fold)
  if (CVMet == "CV1") {
    idx <- do.call(rbind, lapply(seq_along(GPModelMT_List), function(m)
      do.call(rbind, lapply(seq_len(nIter), function(r)
        data.frame(m = m, r = r, t = NA_integer_, fold = seq_len(k),
                   stringsAsFactors = FALSE)))))
  } else if (CVMet == "CV2") {
    idx <- do.call(rbind, lapply(seq_along(GPModelMT_List), function(m)
      do.call(rbind, lapply(seq_len(nIter), function(r)
        do.call(rbind, lapply(seq_len(nTraits), function(t)
          data.frame(m = m, r = r, t = t, fold = seq_len(k),
                     stringsAsFactors = FALSE)))))))
  } else {
    stop("CVMet must be 'CV1' or 'CV2'.")
  }

  # Parallel backend — FORK on Linux (numeric registerDoParallel works inside callr);
  # PSOCK inside callr on Windows causes .doSnowGlobals errors, so fall back to sequential.
  library(doParallel); library(foreach)
  n_cores <- min(max(detectCores() - 1, 1), 10)
  if (.Platform$OS.type == "unix") {
    registerDoParallel(n_cores)
    on.exit(stopImplicitCluster(), add = TRUE)
  } else {
    registerDoSEQ()   # sequential on Windows — PSOCK inside callr hits .doSnowGlobals
  }

  task_results <- foreach(
    j        = seq_len(nrow(idx)),
    .packages = c("BGLR", "rrBLUP"),
    .export   = c("get_yhat", "colPA")
  ) %dopar% {

    m  <- idx$m[j];  r <- idx$r[j];  t <- idx$t[j];  fi <- idx$fold[j]
    tst        <- fold_sets[[r]][[fi]]
    model_name <- GPModelMT_List[m]

    yNA <- PH
    if (CVMet == "CV1") {
      yNA[tst, ] <- NA
    } else {
      yNA[tst, t] <- NA
    }

    ETA <-
      if (model_name == "BRR (BGLR)")        list(list(X = Ximp, model = "BRR"))      else
      if (model_name == "RKHS (BGLR)")       list(list(K = A_tot, model = "RKHS"))    else
      if (model_name == "Spike-Slab (BGLR)") list(list(X = Ximp, model = "SpikeSlab")) else NULL

    # Unique saveAt prefix per task prevents BGLR temp file conflicts across workers
    save_prefix <- file.path(tempdir(), paste0("mt_m", m, "_r", r, "_t",
                                               ifelse(is.na(t), 0L, t), "_f", fi, "_"))
    # First fold of each rep runs verbose so at least one MCMC trace is visible
    verb <- (fi == 1L)

    fm <- tryCatch(
      BGLR::Multitrait(y = yNA, ETA = ETA, nIter = b_iter, burnIn = b_burnin,
                       thin = b_thin, R2 = b_R2, saveAt = save_prefix, verbose = verb),
      error = function(e) NULL
    )

    if (is.null(fm)) return(list(m=m, r=r, t=t, fold=fi, pred=NULL, obs=NULL))

    yhat <- get_yhat(fm)
    if (CVMet == "CV1") {
      list(m=m, r=r, t=NA_integer_, fold=fi,
           pred = yhat[tst, , drop=FALSE],
           obs  =   PH[tst, , drop=FALSE])
    } else {
      list(m=m, r=r, t=t, fold=fi,
           pred = yhat[tst, t],
           obs  =   PH[tst, t])
    }
  }

  # Reassemble flat task results into model × trait accuracy matrix
  out_models <- vector("list", length(GPModelMT_List))
  for (m in seq_along(GPModelMT_List)) {
    rep_list <- vector("list", nIter)
    for (r in seq_len(nIter)) {
      if (CVMet == "CV1") {
        tasks_mr <- Filter(function(x) x$m == m && x$r == r, task_results)
        tasks_mr <- tasks_mr[order(sapply(tasks_mr, `[[`, "fold"))]
        pred_tab <- do.call(rbind, lapply(tasks_mr, `[[`, "pred"))
        obs_tab  <- do.call(rbind, lapply(tasks_mr, `[[`, "obs"))
        rep_list[[r]] <- colPA(pred_tab, obs_tab)
      } else {
        pa_vec <- numeric(nTraits)
        for (t in seq_len(nTraits)) {
          tasks_mrt <- Filter(function(x) x$m==m && x$r==r && !is.na(x$t) && x$t==t,
                              task_results)
          tasks_mrt <- tasks_mrt[order(sapply(tasks_mrt, `[[`, "fold"))]
          pred_vec  <- unlist(lapply(tasks_mrt, `[[`, "pred"), use.names=FALSE)
          obs_vec   <- unlist(lapply(tasks_mrt, `[[`, "obs"),  use.names=FALSE)
          pa_vec[t] <- stats::cor(pred_vec, obs_vec, use = "pairwise.complete.obs")
        }
        rep_list[[r]] <- pa_vec
      }
    }
    M <- do.call(rbind, rep_list)
    out_models[[m]] <- colMeans(M, na.rm = TRUE)
  }

  OUT <- do.call(rbind, out_models)
  colnames(OUT) <- colnames(PH)
  OUT <- round(OUT, 2)
  data.frame(GPModel = GPModelMT_List, OUT, row.names = NULL, check.names = FALSE)
}

####


####
#' Perform multi-trait genomic prediction (getRankedPredictedValuesMT)
#'
#' Description of what getRankedPredictedValuesMT does.
#'
#' @param Data_Table_Num_Filt_List List. A list containing the data tables.
#' @param nTraits Ineger. Integer specifying the number of traits.
#' @param trait Character. Character specifying the traits to be included in the model.
#' @param GPModelMT Character. Character specifying the MT model for fitting the genomic prediction model. Default value is set to "BRR".
#' Other options include "RKHS" and "Spike-Slab"
#' @param optTS Logical.  Logical value specifying whether the optimal train set is used. The defaul NULL value allows the use of the complete training set for training the model.
#'
#' @return List. A list containing the output data frames containing the predicted and observed values for the training and target populations
#' @examples
#' \dontrun{
#' # Example usage of getRankedPredictedValuesMT
#' result <- getRankedPredictedValuesMT(...)
#' }
#' @export
getRankedPredictedValuesMT <- function(Data_Table_Num_Filt_List,nTraits,trait,GPModelMT,optTS=NULL,
                                       b_iter=5000, b_burnin=1000, b_thin=10, b_R2=0.5, b_digits=4, A_ext=NULL){

	 TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	 TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]

### Train data processing

	 set.seed(125)
	 strainID <- as.character(TrainData_Table_Num_Filt[,1])

	 trainSetID <- as.character(as.vector(strainID))
	 optTS <- cbind(trainSetID,c(1:length(trainSetID)))

	 initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	 finalCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[length(grep("ss",colnames(TrainData_Table_Num_Filt)))]


	 if(is.null(optTS)){

		 trainPheno0 <- TrainData_Table_Num_Filt[,trait]
		 geno_012 <- apply(TrainData_Table_Num_Filt[,c(initCol:finalCol)],2,as.numeric)
		 trainGeno <- apply(geno_012,2,function(x) x-1)
		 Geno <- as.character(TrainData_Table_Num_Filt[,1])
	 }else if(!is.null(optTS)){

	     strainID <- as.character(TrainData_Table_Num_Filt[,1])
		 trainSetID <- as.character(as.vector(optTS[,1]))
	     trainIndices <- which(strainID %in% trainSetID)
		 trainPheno0 <- TrainData_Table_Num_Filt[trainIndices,trait]
		 geno_012 <- apply(TrainData_Table_Num_Filt[trainIndices,c(initCol:finalCol)],2,as.numeric)
		 trainGeno <- apply(geno_012,2,function(x) x-1)
		 Geno <- as.character(TrainData_Table_Num_Filt[trainIndices,1])
	 }


### Check if the markers are missing in > 99% of the lines

	 allNAMarkers <- apply(trainGeno,2,function(x) if(length(which(is.na(x)))==nrow(trainGeno)){1}else{0})
     allNAMarkerIndices <- which(allNAMarkers!=0)

	 if(length(allNAMarkerIndices)>=1){
		trainGeno0 <- trainGeno[,-allNAMarkerIndices]
	 }else{trainGeno0 <- trainGeno}


	 if(anyNA(trainGeno0)){
		trainGeno_Imp0 <- snpQC(trainGeno0,impute=TRUE,remove=FALSE)
	 }else{trainGeno_Imp0 <- trainGeno0}


 ### Remove lines with missing pheno

	 if(anyNA(trainPheno0)){
	   ph_NA_Indices_Tab <- apply(trainPheno0,2,function(x) which(is.na(x)))
	   ph_NA_Indices <- unique(c(unlist(ph_NA_Indices_Tab)))
	   trainPheno <- trainPheno0[-ph_NA_Indices,]
	   trainGeno_Imp <- trainGeno_Imp0[-ph_NA_Indices,]
	   Geno <- as.character(TrainData_Table_Num_Filt[-ph_NA_Indices,1])
	 }else if(!anyNA(trainPheno0)){

	   trainPheno <- trainPheno0
	   trainGeno_Imp <- trainGeno_Imp0
	   Geno <- as.character(TrainData_Table_Num_Filt[,1])
	 }

####### Set the same set of markers for both train and test sets in the same order
##Check with sorted markers and check with test set
	  trainssIDs <- colnames(trainGeno_Imp)


#### Test geno table
	if(!is.null(TestData_Table_Num_Filt)){
	 TestGenoTable <- apply(TestData_Table_Num_Filt[,-1],2,as.numeric)
	 rownames(TestGenoTable) <- TestData_Table_Num_Filt[,1]

	 if(anyNA(TestGenoTable)){
	    testGeno_Imp0 <- snpQC(TestGenoTable,impute=TRUE,remove=FALSE)
		testGeno_Imp <- apply(testGeno_Imp0,2,function(x) as.numeric(x-1))

		testNA <- apply(testGeno_Imp,2,function(x) length(which(is.na(x))))
		testNAIndices <-(which(unlist(testNA) !=0))
		if(length(testNAIndices)<1){
		   testNAIndices <- NULL
		}
	 }else{
	    TestGenoTable <- apply(TestData_Table_Num_Filt[,-1],2,as.numeric)
	    rownames(TestGenoTable) <- TestData_Table_Num_Filt[,1]
	    testGeno_Imp <- apply(TestGenoTable,2,function(x) as.numeric(x-1))
	    testNAIndices <- NULL
	 }

### Make sure trainGenoImp and testGenoImp have the same markers in matching order and ensure no NA in testGeno_Imp

  	  testssIDs <- colnames(testGeno_Imp)
	  mtchInd <- match(trainssIDs,testssIDs)
	  matchInd <- mtchInd[which(!is.na(mtchInd))]
	  testGeno_Imp_Ord <- testGeno_Imp[,matchInd]
	  testGeno_Imp <- testGeno_Imp_Ord

    }else{testGeno_Imp <- NULL}

## Kinship Matrix

	 A <- if (!is.null(A_ext) && all(Geno %in% rownames(A_ext))) {
	   A_ext[Geno, Geno]
	 } else {
	   A.mat(trainGeno_Imp)
	 }
	 colnames(A) <- Geno
	 rownames(A) <- Geno

	 Y <- trainPheno

    if(!is.null(testGeno_Imp)){


	### Make sure trainGenoImp and testGenoImp have the same markers in matching order and ensure no NA in X
	 testssIDs <- colnames(testGeno_Imp)
	 trainssIDs <- colnames(trainGeno_Imp)
	 mtchInd2 <- match(testssIDs,trainssIDs)
	 matchInd2 <- mtchInd2[which(!is.na(mtchInd2))]
	 trainGeno_Imp_Ord <- trainGeno_Imp[,matchInd]
	 trainGeno_Imp <- trainGeno_Imp_Ord

	 X <- rbind(trainGeno_Imp,testGeno_Imp)
	 #rownames(X) <- c(rownames(trainGeno_Imp),rownames(TestGenoTable))
	 n <- nrow(X)
	 p <- ncol(X)

	 nTst <- nrow(testGeno_Imp)
	 Ytst <- matrix(rep(NA,nTst*ncol(Y)),nrow=nTst,ncol=ncol(Y))
	 colnames(Ytst) <- colnames(Y)
	 yNA <- as.matrix(rbind(Y,Ytst))
	 yNA_DF <- as.data.frame(yNA)
	 #yNA_DF$id <- as.factor(paste(c(Geno,rownames(testGeno_Imp)),"_",yNA_DF[,1],sep=""))

	}else if(is.null(testGeno_Imp)){
	 X <- trainGeno_Imp
	 n <- nrow(X)
	 p <- ncol(X)
	 yNA <- as.matrix(Y)
	 yNA_DF <- as.data.frame(yNA)
	 #yNA_DF$id <- as.factor(paste(Geno,"_",yNA_DF[,1],sep=""))
	}

	if(GPModelMT == "BRR (BGLR)"){

		ETA <- list(list(X=X,model="BRR"))
		fm3 <- BGLR::Multitrait(y=yNA,ETA=ETA,nIter=b_iter,burnIn=b_burnin,thin=b_thin,R2=b_R2,verbose=TRUE)
	}else if(GPModelMT == "RKHS (BGLR)"){
	        A.Tot <- if (!is.null(A_ext) && all(rownames(X) %in% rownames(A_ext))) {
	          A_ext[rownames(X), rownames(X)]
	        } else { A.mat(X) }
		ETA <- list(list(K=A.Tot,model="RKHS"))
		fm3 <- BGLR::Multitrait(y=yNA,ETA=ETA,nIter=b_iter,burnIn=b_burnin,thin=b_thin,R2=b_R2,verbose=TRUE)
	}else if(GPModelMT == "Spike-Slab (BGLR)"){

		ETA <- list(list(X=X,model="SpikeSlab"))
		fm3 <- BGLR::Multitrait(y=yNA,ETA=ETA,nIter=b_iter,burnIn=b_burnin,thin=b_thin,R2=b_R2,verbose=TRUE)
	 }else if(GPModelMT == "GBLUP (SOMMER)"){

		A.Tot <- if (!is.null(A_ext) && all(rownames(X) %in% rownames(A_ext))) {
	  A_ext[rownames(X), rownames(X)]
	} else { A.mat(X) }
		    rownames(A.Tot) <- yNA_DF$id
		colnames(A.Tot) <- yNA_DF$id
		fm3 <- mmer(as.formula(paste("cbind(",paste(trait,collapse=","),")~1",sep="")),
            random=~vs(id,Gu=A.Tot),
            rcov=~units,data=yNA_DF,verbose = TRUE)
	 }


### Upper bound of Reliability

 cleanREPV2 <- function(y,gen,fam=NULL,thr=0.95){

  if(is.vector(y)) y=matrix(y,ncol=1)
  if(is.null(fam)) fam = rep(1,nrow(y))

  ibs = function(gen){
  f1 = function(x,gen) apply(gen,1,function(y,x) mean(y==x,na.rm = T),x=x)
  f2 = apply(gen,1,f1,gen=gen)
  return(f2)}

  GG = function(gen, r = 1) {
    a1 = (gen - 1)
    a1[a1 == -1] = 0
    A1 = (tcrossprod(a1))
    a2 = -(gen - 1)
    a2[a2 == -1] = 0
    A2 = (tcrossprod(a2))
    d = round(exp(-abs(gen - 1)))
    D = tcrossprod(d)
    G = A1 + A2 + D
    G = (G/ncol(gen))^r
    return(G)
  }
  cat("solving identity matrix\n")
  G = GG(gen)

  rownames(G) = 1:nrow(G)

  lt = G*lower.tri(G) # lower triang
  r = 1* lt>thr # logical matrix: repetitions

  # starting point of new data
  #rownames(gen) = 1:nrow(gen)
  Ny=y;  Nfam=fam;  Ngen=gen
  NGen <- gen
  # summary
  cs = colSums(r) # how many times id will be repeated
  while(any(cs>0)){

    i = which(cs>0)[1]
    cat("indiviual",rownames(gen)[i],"had",cs[i],'duplicate(s)\n')
    w = which(r[,i])
    if(ncol(y)>1){y[i,] = colMeans(y[c(i,w),],na.rm=T)
    }else{y[i] = mean(y[c(i,w)],na.rm=T)}
    if(ncol(y)>1){Ny=Ny[-w,] ;rownames(Ny) <- rownames(Ny)[-w]}else{Ny=Ny[-w];names(Ny) <- names(Ny)[-w]}
    Nfam=Nfam[-w]
    Ngen=NGen[-w,]
	rownames(Ngen) <- rownames(NGen)[-w]
	NGen <- Ngen
    r = r[-w,]
    cs = colSums(r)
  }
  return(list(y=Ny,gen=Ngen,fam=Nfam))
}


    cleanData <- cleanREPV2(trainPheno[,1],apply(trainGeno_Imp,2,function(x) x+1))
    M <-  cleanData[[2]]
    M.Pdt <- t(M)%*% solve(M %*% t(M) + diag(nrow(M))) %*% M

    getU <- function(M.Pdt,v){
	   v.hat <- M.Pdt %*% v
	   U <- (t(v.hat)%*% v.hat)/ (t(v) %*% v)
	   return(U)
    }


	trnInd <- c(1:length(trainSetID))
	Train_PredictedValues <- fm3$ETAHat[trnInd,]
    Train_SortedPredictedValues <- sort.int(Train_PredictedValues[,1],decreasing=TRUE,index.return=TRUE)
	Train_StrainID <- rownames(yNA_DF)

   # U <- apply(trainGeno_Imp,1,function(x) getU(M.Pdt,x))

### Sort Output

	Train_SortedIndices <- Train_SortedPredictedValues[[2]]
	Sorted_Train_StrainID <- Train_StrainID[Train_SortedIndices]

	#U_Sorted <- U[(Train_SortedIndices)]

	Train_SortedPredValuesTab <- apply(Train_PredictedValues[Train_SortedIndices,],2,function(x) round(x,digits=2))

## Output DF

	TrainOutputDF <- cbind.data.frame(Sorted_Train_StrainID,Train_SortedPredValuesTab)
	colnames(TrainOutputDF) <- c("LineID",colnames(Train_SortedPredValuesTab))

	if(!is.null(testGeno_Imp)){

	  ### grep tstIDs and trainIDs
		tst <- c((length(trainSetID)+1):nrow(X))

		if(GPModelMT != "GBLUP (SOMMER)"){
		  Test_PredictedValues <- fm3$ETAHat[tst,]
		}else if(GPModelMT == "GBLUP (SOMMER)"){
		  UMat <- do.call(cbind,lapply(fm3$U$`u:id`,function(x) x))
		  Test_PredictedValues <- c()
		  for(nT in 1:length(trait)){
			Test_PredictedValues <- cbind.data.frame(Test_PredictedValues,(UMat[,nT]+fm3$Beta[nT,3])[tst])
		  }
		}

	   Test_SortedPredictedValues <- sort.int(Test_PredictedValues[,1],decreasing=TRUE,index.return=TRUE)
	   Test_StrainID <- rownames(TestData_Table_Num_Filt)


	   U <- apply(testGeno_Imp,1,function(x) getU(M.Pdt,x))

### Sort Output

		Test_SortedIndices <- Test_SortedPredictedValues[[2]]
		Sorted_Test_StrainID <- Test_StrainID[Test_SortedIndices]
		U_Sorted <- U[(Test_SortedIndices)]
		Test_SortedPredValuesTab <- apply(Test_PredictedValues[Test_SortedIndices,],2,function(x) round(x,digits=2))

		### Output DF

		TargetOutputDF <- cbind.data.frame(Sorted_Test_StrainID ,Test_SortedPredValuesTab,round(U_Sorted,digits=5))
		colnames(TargetOutputDF) <- c("LineID",paste("Predicted Values for ",colnames(Test_SortedPredValuesTab),sep=""),"Upper Bound of Reliability")

        outputDF_List <- list(TrainOutputDF,TargetOutputDF)

       }else{

	     TargetOutputDF <- "Target set not defined. Select training set to view from side panel."
	     outputDF_List <- list(TrainOutputDF,TargetOutputDF)

	   }

	 return(outputDF_List)

 }

###
###

####
#' Summarize multi-trait CV fit results
#' @param PA_df Data frame. CV results from getMTCVR.
#' @param CVMet Character. "CV1" or "CV2" (default "CV1").
#' @param traits Character vector. Trait names (default NULL).
#' @return Data frame with mean/SD/min/max PA and RMSE per model and trait.
#' @export
summarize_MT_Fits <- function(PA_df, CVMet = "CV1", traits = NULL) {

  # traits defaults to all non-GPModel columns
  if (is.null(traits)) {
    traits <- setdiff(colnames(PA_df), "GPModel")
  }

  # --- build PA_long via manual loop (no tidyverse dependency) ---
  long_rows <- vector("list", nrow(PA_df) * length(traits))
  idx <- 1L
  for (i in seq_len(nrow(PA_df))) {
    for (tr in traits) {
      long_rows[[idx]] <- data.frame(
        GPModel = PA_df$GPModel[i],
        Trait   = tr,
        PA      = PA_df[[tr]][i],
        CVMet   = CVMet,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  PA_long <- do.call(rbind, long_rows)

  # --- best model per trait ---
  best_model <- vapply(traits, function(tr) {
    vals <- PA_df[[tr]]
    best_row <- which.max(vals)
    if (length(best_row)) PA_df$GPModel[best_row] else NA_character_
  }, character(1L))

  list(
    PA_wide    = PA_df,
    PA_long    = PA_long,
    best_model = best_model,
    CVMet      = CVMet
  )

}
