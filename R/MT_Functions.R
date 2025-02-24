
#' Title of getMTCVR
#'
#' Description of what getMTCVR does.
#'
#' @param Data_Table_Num_Filt_List Description of Data_Table_Num_Filt_List.
#' @param trait Description of trait.
#' @param nTraits Description of nTraits.
#' @param k Description of k.
#' @param nIter Description of nIter.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getMTCVR
#' result <- getMTCVR(...)
#' }
getMTCVR <- function(Data_Table_Num_Filt_List,trait,nTraits,k,nIter){
  
     TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	
	 pheno_wNA <- apply(TrainData_Table_Num_Filt[,trait],2,function(x)  as.numeric(as.character(x)))
	 genoInd <- grep("ss",colnames(TrainData_Table_Num_Filt))
	 geno_012 <- apply(TrainData_Table_Num_Filt[,genoInd],2,as.numeric)
	 geno_wNA <- apply(geno_012,2,function(x) x-1)
	 Geno_wNA <- rownames(TrainData_Table_Num_Filt)
	 
	 NAIndices <- unlist(apply(pheno_wNA,2,function(x)which(is.na(x))))
     if(length(NAIndices) >0){
		 pheno <- pheno_wNA[-NAIndices,]
		 geno <- geno_wNA[-NAIndices,]
		 Geno <- Geno_wNA[-NAIndices]
	 }else if(length(NAIndices)==0){ 
	     pheno <- pheno_wNA
		 geno <- geno_wNA
		 Geno <- Geno_wNA
	 }
	 if(anyNA(geno)){	 
		geno_Imp <- snpQC(geno,impute=TRUE,remove=FALSE,MAF=0.02)
	 }else{geno_Imp <- geno}

     n <- nrow(geno_Imp)
##############################################################
	 
	 # nTr <- round(n*((k-1)/k),digits=0) 
	 # trainIndices <- sample(c(1:n),nTr)
	 # testIndices <- setdiff(c(1:n),trainIndices)
	  
	 # trPheno <- pheno[trainIndices,]
	 # trGeno <- geno_Imp[trainIndices,]
	 # testPheno <- pheno[testIndices,]
	 # testGeno <- geno_Imp[testIndices,]
	 
	 # trGenoNames <- Geno[trainIndices]
	 # tstGenoNames <- Geno[testIndices]
	 
############################################################

  GPModelMT_List <- list("BRR (BGLR)","RKHS (BGLR)","Spike-Slab (BGLR)") #,"GBLUP (SOMMER)")
  PA_Out_GP <- list()
  
  

  for(nGP in 1:length(GPModelMT_List)){ 
    
	GPModelMT <- GPModelMT_List[[nGP]]
	
	PA_List <- list()
	
	for(nrep in 1:nIter){
     nK <- floor(n/k)
	 k_List <- list()
	 set.seed(125+nrep) 
	 tot <- c(1:n)
	 Test_Pred <- list()
	 Y_Tst <- list() 
     for(i in 1:k){
	   k_List[[i]] <- sample(tot,nK)
	   tot <- setdiff(tot,k_List[[i]]) 
	 }
	 trIndices_List <- list()
	 tstIndices_List <- list()
     for(i in 1:k){ 
	   trIndices_List[[i]] <- unlist(k_List[-i])
	   tstIndices_List[[i]] <- k_List[[i]]
	 }	  
	 PA <- list()
	 
	for(i in 1:k){
	 trainIndices <-  trIndices_List[[i]]
	 testIndices <- tstIndices_List[[i]]
	#########################################################
	 Y <- pheno
	 yNA <- pheno
	 nTst <- length(testIndices)
	 yNA[testIndices,] <- matrix(rep(NA,nTst*ncol(yNA)),nrow=nTst,ncol=ncol(yNA))
	 
## Kinship Matrix 
     rownames(geno_Imp) <- Geno
	 rownames(Y) <- Geno
	 rownames(yNA) <- Geno
	
	 X <- geno_Imp
	 rownames(X) <- rownames(geno_Imp)
	 n <- nrow(X)
	 p <- ncol(X)
	 
	 yNA_DF <- as.data.frame(yNA)
	 yNA_DF$id <- as.factor(rownames(X))
	 
	 if(GPModelMT == "BRR (BGLR)"){
				  
			ETA <- list(list(X=X,model="BRR"))
			fm3 <- BGLR::Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
	 }else if(GPModelMT == "RKHS (BGLR)"){ 
	        A.Tot <- rrBLUP::A.mat(X)
			ETA <- list(list(K=A.Tot,model="RKHS"))
			fm3 <- BGLR::Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
	 }else if(GPModelMT == "Spike-Slab (BGLR)"){ 
				  
			ETA <- list(list(X=X,model="SpikeSlab"))
			fm3 <- BGLR::Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
	 }	
	 # if(GPModelMT == "GBLUP (SOMMER)"){ 
	         
			# A.Tot <- A.mat(X)
			# rownames(A.Tot) <- rownames(X) 
			# colnames(A.Tot) <- rownames(X)
			# fm3 <- sommer::mmer(as.formula(paste("cbind(",paste(trait,collapse=","),")~1",sep="")),
            # random=~vs(id,Gu=A.Tot),
            # rcov=~units,
            # data=yNA_DF,verbose = TRUE)
	 # }
	 
	 tst <- testIndices
	
	 if(GPModelMT != "GBLUP (SOMMER)"){ 
	  Test_PredictedValues <- fm3$ETAHat[tst,]
	 }
	 
	 # if(GPModelMT == "GBLUP (SOMMER)"){  
	  # UMat <- do.call(cbind,lapply(fm3$U$`u:id`,function(x) x))
	  # Test_PredictedValues <- c()
	  # for(nT in 1:length(trait)){
	    # Test_PredictedValues <- cbind(Test_PredictedValues,(UMat[,nT]+fm3$Beta[nT,3])[tst])
		
	  # }
	 # }
	 	 
	 # PA[[i]] <- diag(cor(Test_PredictedValues,Y[tst,]))
	 
	  Test_Pred[[i]] <- Test_PredictedValues
	  Y_Tst[[i]] <- Y[tst,]
	 
	}
	
	Test_Pred_Tab <- do.call(rbind,lapply(Test_Pred,function(x) x))
	Y_Tst_Tab <- do.call(rbind,lapply(Y_Tst,function(x) x))
	
	PA <- diag(cor(Test_Pred_Tab,Y_Tst_Tab))
	  
	#PA_List[[nrep]] <- do.call(rbind,lapply(PA,function(x) x))
	 
	 
	 PA_List[[nrep]] <- PA
	  
	}
	
	PA_Out_GP[[nGP]] <- apply(do.call(rbind,lapply(PA_List,function(x)x)),2,mean)
	
  }
	
	PA_Out1 <-  apply(do.call(rbind,lapply(PA_Out_GP,function(x) x)),2,function(x) round(x,digits=2))
	PA_Out2 <- cbind(unlist(GPModelMT_List),PA_Out1)
	
	#PA_Out <- rbind(c("GPModel",colnames(PA_Out1)),PA_Out2)
	PA_Out <- PA_Out2
	colnames(PA_Out) <- c("GPModel",colnames(PA_Out1))
	
	return(PA_Out)
   
 }


####
#' Title of getRankedPredictedValuesMT
#'
#' Description of what getRankedPredictedValuesMT does.
#'
#' @param Data_Table_Num_Filt_List Description of Data_Table_Num_Filt_List.
#' @param nTraits Description of nTraits.
#' @param trait Description of trait.
#' @param GPModelMT Description of GPModelMT.
#' @param optTS=NULL Description of optTS=NULL.
#'
#' @return Description of what is returned.
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
 