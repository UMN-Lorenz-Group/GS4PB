
#' Perform cross validation for multi-trait genomic prediction in CV1 and CV2 schemes #'
#' 
#'
#' @param Data_Table_Num_Filt_List List. A list containing the data tables.
#' @param trait Character. Character specifying the traits to be included in the model.
#' @param nTraits Ineger. Integer specifying the number of traits.
#' @param k Integer. Integer specifying the number of folds in k-fold CV
#' @param nIter Integer. Integer value specifying the number of iterations.
#' @Param CVMet Character. Character specifying the cross validation scheme with the default value set to "CV1". 
#' Options include "CV1" and "CV2"
#'
#' @return a data frame that contains output stats such as prediction ability for the mlti-trait models .
#' @examples
#' \dontrun{
#' # Example usage of getMTCVR
#' result <- getMTCVR(...)
#' }

getMTCVR <- function(Data_Table_Num_Filt_List, trait, nTraits, k, nIter, CVMet) {

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
  A_tot <- rrBLUP::A.mat(Ximp)

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
  out_models <- vector("list", length(GPModelMT_List))

  for (m in seq_along(GPModelMT_List)) {
    model_name <- GPModelMT_List[m]
    rep_list <- vector("list", nIter)

    for (r in seq_len(nIter)) {
      folds <- make_folds(n, k, seed = 125 + r)

      # storage per fold
      pred_all <- list()
      obs_all  <- list()

      if (CVMet == "CV1") {
        for (i in seq_along(folds)) {
          tst <- folds[[i]]
          yNA <- PH
          yNA[tst, ] <- NA

          ETA <-
            if (model_name == "BRR (BGLR)")      list(list(X = Ximp, model = "BRR")) else
            if (model_name == "RKHS (BGLR)")     list(list(K = A_tot, model = "RKHS")) else
            if (model_name == "Spike-Slab (BGLR)") list(list(X = Ximp, model = "SpikeSlab")) else NULL

          fm <- BGLR::Multitrait(y = yNA, ETA = ETA, nIter = 1000, burnIn = 500, verbose = FALSE)
          yhat <- get_yhat(fm)

          pred_all[[i]] <- yhat[tst, , drop = FALSE]
          obs_all[[i]]  <- PH[tst, ,  drop = FALSE]
        }

        pred_tab <- do.call(rbind, pred_all)
        obs_tab  <- do.call(rbind,  obs_all)
        rep_list[[r]] <- colPA(pred_tab, obs_tab)

      } else if (CVMet == "CV2") {
        # evaluate each trait by considering that trait as focal trait 
		# and masking ONLY that trait in the test set
        pa_vec <- numeric(nTraits)
        for (t in seq_len(nTraits)) {
          pred_t <- list()
          obs_t  <- list()
          for (i in seq_along(folds)) {
            tst <- folds[[i]]
            yNA <- PH
            yNA[tst, t] <- NA
            ETA <-
              if (model_name == "BRR (BGLR)")      list(list(X = Ximp, model = "BRR")) else
              if (model_name == "RKHS (BGLR)")     list(list(K = A_tot, model = "RKHS")) else
              if (model_name == "Spike-Slab (BGLR)") list(list(X = Ximp, model = "SpikeSlab")) else NULL
ls()

            fm <- BGLR::Multitrait(y = yNA, ETA = ETA, nIter = 1000, burnIn = 500, verbose = FALSE)
            yhat <- get_yhat(fm)
            pred_t[[i]] <- yhat[tst, t]
            obs_t[[i]]  <- PH[tst, t]
          }
          pred_vec <- unlist(pred_t, use.names = FALSE)
          obs_vec  <- unlist(obs_t,  use.names = FALSE)
          pa_vec[t] <- stats::cor(pred_vec, obs_vec, use = "pairwise.complete.obs")
        }
        rep_list[[r]] <- pa_vec

      } else {
        stop("CVMet must be 'CV1' or 'CV2'.")
      }
    } # end replicates

    # average across replicates
    M <- do.call(rbind, rep_list)
    out_models[[m]] <- colMeans(M, na.rm = TRUE)
  }

  OUT <- do.call(rbind, out_models)
  colnames(OUT) <- colnames(PH)
  OUT <- round(OUT, 2)
  data.frame(GPModel = GPModelMT_List, OUT, row.names = NULL, check.names = FALSE)
}



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
getRankedPredictedValuesMT <- function(Data_Table_Num_Filt_List,nTraits,trait,GPModelMT,optTS=NULL){ 
   
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
    	 
	 A <- A.mat(trainGeno_Imp)	
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
			fm3 <- Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
	}else if(GPModelMT == "RKHS (BGLR)"){ 
	        A.Tot <- A.mat(X)
			ETA <- list(list(K=A.Tot,model="RKHS"))
			fm3 <- Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
	}else if(GPModelMT == "Spike-Slab (BGLR)"){ 
				  
			ETA <- list(list(X=X,model="SpikeSlab"))
			fm3 <- Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
	 }
	 

### Upper bound of Reliability
 
    cleanData <- cleanREPV2(trainPheno[,1],apply(trainGeno_Imp,2,function(x) x))
    M <-  cleanData[[2]]
    M.Pdt <- t(M)%*% solve(M %*% t(M) + diag(nrow(M))) %*% M
      
   
	
	trnInd <- c(1:length(trainSetID))
	Train_PredictedValues <- fm3$ETAHat[trnInd,]
    Train_SortedPredictedValues <- sort.int(Train_PredictedValues[,1],decreasing=TRUE,index.return=TRUE) 
	Train_StrainID <- rownames(yNA_DF)	 

   

### Sort Output
	
	Train_SortedIndices <- Train_SortedPredictedValues[[2]]
	Sorted_Train_StrainID <- Train_StrainID[Train_SortedIndices]
	
	Train_SortedPredValuesTab <- apply(Train_PredictedValues[Train_SortedIndices,],2,function(x) round(x,digits=2))
		
## Output DF
		
	TrainOutputDF <- cbind.data.frame(Sorted_Train_StrainID,Train_SortedPredValuesTab)
	colnames(TrainOutputDF) <- c("LineID",colnames(Train_SortedPredValuesTab))
	
	if(!is.null(testGeno_Imp)){
	
	  ### grep tstIDs and trainIDs
		tst <- c((length(trainSetID)+1):nrow(X))
		
		if(GPModelMT != "GBLUP (SOMMER)"){ 
		  Test_PredictedValues <- fm3$ETAHat[tst,]
		}
		# else if(GPModelMT == "GBLUP (SOMMER)"){  
		  # UMat <- do.call(cbind,lapply(fm3$U$`u:id`,function(x) x))
		  # Test_PredictedValues <- c()
		  # for(nT in 1:length(trait)){
			# Test_PredictedValues <- cbind.data.frame(Test_PredictedValues,(UMat[,nT]+fm3$Beta[nT,3])[tst])
		  # }
		# }
	 
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
 