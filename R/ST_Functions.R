
	
#' Title of getemCVR
#'
#' Description of what getemCVR does.
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
#' # Example usage of getemCVR
#' result <- getemCVR(...)
#' }
getemCVR <- function(Data_Table_Num_Filt_List,trait,nTraits,k,nIter){
  
     TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	 initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	 finalCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[length(grep("ss",colnames(TrainData_Table_Num_Filt)))]
	 
	 pheno0 <- as.numeric(as.character(TrainData_Table_Num_Filt[,trait]))
	 
	  
	 geno_012 <- apply(TrainData_Table_Num_Filt[,c(initCol:finalCol)],2,as.numeric)
	 geno <- apply(geno_012,2,function(x) x-1)
	 Geno <- rownames(TrainData_Table_Num_Filt)
  
     if(anyNA(geno)){	 
		geno_Imp0 <- snpQC(geno,impute=TRUE,remove=FALSE)
	 }else{geno_Imp0 <- geno}
	 
	 if(anyNA(pheno0)){ 
	  ph_NA_Indices  <- which(is.na(pheno0))
	  pheno <- pheno0[-ph_NA_Indices] 
	  geno_Imp <- geno_Imp0[-ph_NA_Indices,] 
     } 
     if(!anyNA(pheno0)){ 
	   pheno <- pheno0
	   geno_Imp <- geno_Imp0
     }
	 
     emCVR <- emCV(pheno,geno_Imp,k,nIter)
  
     return(emCVR)
  }
  
################################# 

	
### 


 #' Title of getRankedPredictedValues
#'
#' Description of what getRankedPredictedValues does.
#'
#' @param Data_Table_Num_Filt_List Description of Data_Table_Num_Filt_List.
#' @param nTraits Description of nTraits.
#' @param trait Description of trait.
#' @param GPModel Description of GPModel.
#' @param fixedX=NULL Description of fixedX=NULL.
#' @param fixedData=NULL Description of fixedData=NULL.
#' @param optTS=NULL Description of optTS=NULL.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getRankedPredictedValues
#' result <- getRankedPredictedValues(...)
#' }
getRankedPredictedValues <- function(Data_Table_Num_Filt_List,nTraits,trait,GPModel,fixedX=NULL,fixedData=NULL,optTS=NULL){ 
    
     if(!is.null(fixedX) & fixedX != "NULL" & fixedData !="NULL"){ 
	   GPModel <- "rrBLUP (rrBLUP)"
	   Fixed.X <- fixedData[[1]]
	   Test.X <- fixedData[[2]]
	 }
	   
	 TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	 TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
	 
	
## TrainData Table Processing 
	 
	  initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	  finalCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[length(grep("ss",colnames(TrainData_Table_Num_Filt)))]
	  
	
	  if(is.null(optTS)){
	    
		 trainPheno0 <- as.numeric(as.character(TrainData_Table_Num_Filt[,trait]))
		 geno_012 <- apply(TrainData_Table_Num_Filt[,c(initCol:finalCol)],2,as.numeric)
		 trainGeno <- apply(geno_012,2,function(x) x-1)
		 Geno <- as.character(TrainData_Table_Num_Filt[,1])
		 rownames(trainGeno) <- Geno
	  }
	 
	  if(!is.null(optTS)){
	     strainID <- as.character(TrainData_Table_Num_Filt[,1])
		 StrainID <- as.character(TrainData_Table_Num_Filt[,ncol(TrainData_Table_Num_Filt)])
		 trainSetID <- as.character(as.vector(optTS))
	     trainIndices <- which(StrainID %in% trainSetID)
		 trainPheno0 <- as.numeric(as.character(TrainData_Table_Num_Filt[trainIndices,trait]))
		 geno_012 <- apply(TrainData_Table_Num_Filt[trainIndices,c(initCol:finalCol)],2,as.numeric)
		 trainGeno <- apply(geno_012,2,function(x) x-1)
		 Geno <- as.character(TrainData_Table_Num_Filt[trainIndices,1])
		 rownames(trainGeno) <- Geno
	  } 
	 
### Check if the markers are missing in > 99% of the lines 
	 
	  allNAMarkers <- apply(trainGeno,2,function(x) if(length(which(is.na(x)))==nrow(trainGeno)){1}else{0})
     
	  allNAMarkerIndices <- which(allNAMarkers!=0)
	  
	  if(length(allNAMarkerIndices)>=1){
		 trainGeno0 <- trainGeno[,-allNAMarkerIndices]
	  }else{ trainGeno0 <- trainGeno}
	 
	 
	 
	   if(anyNA(trainGeno0)){	 
		 trainGeno_Imp0 <- snpQC(trainGeno0,impute=TRUE,remove=FALSE)
	   }else{trainGeno_Imp0 <- trainGeno0}
	 
	 
	
	 
	

 ### Remove lines with missing pheno	 
	  
	if(anyNA(trainPheno0)){
	   ph_NA_Indices <- which(is.na(trainPheno0))
	   trainPheno <- trainPheno0[-ph_NA_Indices] 
	   trainGeno_Imp <- trainGeno_Imp0[-ph_NA_Indices,]
	  # Geno <- as.character(TrainData_Table_Num_Filt[-ph_NA_Indices,1])
	   Geno <- as.character(rownames(trainGeno_Imp))
	   if(!is.null(fixedX)  & fixedX !="NULL" & fixedData != "NULL" ){
	     Fixed.X.Mod <- Fixed.X[-ph_NA_Indices,]
	   }
	 }
	 if(!anyNA(trainPheno0)){
	  
	   trainPheno <- trainPheno0
	   trainGeno_Imp <- trainGeno_Imp0
	  # Geno <- as.character(TrainData_Table_Num_Filt[,1])
	   Geno <- as.character(rownames(trainGeno_Imp))
	   if(!is.null(fixedX)  & fixedX !="NULL" & fixedData != "NULL"){
	     Fixed.X.Mod <- Fixed.X
	   }
	 }
	 
	 trainssIDs <- colnames(trainGeno_Imp)
	 
####### Set the same set of markers for both train and test sets with the same order 

#### Test Geno Table
	
   	 if(!is.null(TestData_Table_Num_Filt)){
	   TestGenoTable <- apply(TestData_Table_Num_Filt[,-1],2,as.numeric)
	 }else{TestGenoTable <- NULL}
	
	 if(!is.null(TestGenoTable)){
	 
		 if(anyNA(TestGenoTable)){ 
			testGeno_Imp0 <- snpQC(TestGenoTable,impute=TRUE,remove=FALSE)
			testGeno_Imp <- apply(testGeno_Imp0,2,function(x) as.numeric(x-1))
			
			testNA <- apply(testGeno_Imp,2,function(x) length(which(is.na(x))))
			testNAIndices <-(which(unlist(testNA) !=0))
			if(length(testNAIndices)<1){
			   testNAIndices <- NULL
			}
		 }else{testGeno_Imp <- apply(TestGenoTable,2,function(x) as.numeric(x-1)) 
			testNAIndices <- NULL
		 }
	 
	  ##Check with sorted markers and check with test set 
		  
	  testssIDs <- colnames(testGeno_Imp) 
	  testGeno_Imp_Ord <- testGeno_Imp[, match(trainssIDs,testssIDs)]
	  testGeno_Imp <- testGeno_Imp_Ord
	 
	 }else if(is.null(TestGenoTable)){testGeno_Imp <- NULL}
      


	  
## Kinship Matrix
	 # A <- A.mat(trainGeno_Imp)	
	 # colnames(A) <- Geno
	 # rownames(A) <- Geno

## Prepare Data Table for GP 
	 # Data <- cbind.data.frame(Geno,trainPheno)
	 # colnames(Data) <- c("Geno","Pheno")
	 # Geno <- "Geno"
	 # Pheno <- "Pheno"
	
#### Impute trainGeno and train using mixed.solve
		 
	if(GPModel == "rrBLUP (rrBLUP)"){ 
	  	
	 if(!is.null(fixedX) & fixedX !="NULL" & fixedData != "NULL"){
	   pred <- mixed.solve(trainPheno,Z=trainGeno_Imp,X=Fixed.X.Mod,SE=FALSE,return.Hinv =FALSE) 
	   Mean <- Fixed.X.Mod %*% pred$beta
	 }
	 if(is.null(fixedX) | fixedX =="NULL" | fixedData == "NULL"){
	     pred <- mixed.solve(trainPheno,Z=trainGeno_Imp,SE=FALSE,return.Hinv =FALSE) 
	     Mean <- as.numeric(pred$beta)
	 }
	 
 	 Effects <- pred$u

	
 	 PredictedValues <- Mean + (trainGeno_Imp %*% Effects)
 	 SortedPredictedValues <- sort.int(PredictedValues,decreasing=TRUE,index.return=TRUE)
 	  
   
	 if(!is.null(fixedX) & fixedX !="NULL" & fixedData != "NULL" ){
	        Mean.tst <- Test.X %*% pred$beta
	 }
	 if(is.null(fixedX) | fixedX =="NULL" | fixedData == "NULL" ){
		    Mean.tst <- as.numeric(pred$beta)
	  }
     
	 if(!is.null(testGeno_Imp)){
	   
	   while(anyNA(testGeno_Imp)){
		  testGeno_Imp_Mod <- testGeno_Imp
		  trainGeno_Imp_Mod <- trainGeno_Imp
		  Effects_Mod <- Effects
		  testNA <- apply(testGeno_Imp_Mod,2,function(x) length(which(is.na(x))))
		  testNAIndices <-(which(unlist(testNA) !=0))
		  testGeno_Imp <-testGeno_Imp_Mod[,-testNAIndices]
		  Effects <- Effects_Mod[-testNAIndices]
		  trainGeno_Imp <- trainGeno_Imp_Mod[,-testNAIndices]
       }
   		 		  
	   Test_PredictedValues <-  Mean.tst + (testGeno_Imp %*% Effects)
	   Test_SortedPredictedValues <- sort.int(Test_PredictedValues,decreasing=TRUE,index.return=TRUE) 
	   Test_StrainID <- rownames(TestData_Table_Num_Filt)
    
	  }
	 }
	 
	 
	if(GPModel == "rrBLUP (bWGR)"){
	 
		 pred <- emRR(trainPheno,trainGeno_Imp) 
		 Mean <- as.numeric(pred$mu)
		 Effects <- pred$b
			 
		 PredictedValues <- Mean + (trainGeno_Imp %*% Effects)
		 SortedPredictedValues <- sort.int(PredictedValues,decreasing=TRUE,index.return=TRUE)
		
		if(!is.null(testGeno_Imp)){
			 while(anyNA(testGeno_Imp)){
			   testGeno_Imp_Mod <- testGeno_Imp
			   trainGeno_Imp_Mod <- trainGeno_Imp
			   Effects_Mod <- Effects
			   testNA <- apply(testGeno_Imp_Mod,2,function(x) length(which(is.na(x))))
			   testNAIndices <-(which(unlist(testNA) !=0))
			   testGeno_Imp <-testGeno_Imp_Mod[,-testNAIndices]
			   Effects <- Effects_Mod[-testNAIndices]
			   trainGeno_Imp <- trainGeno_Imp_Mod[,-testNAIndices]
			 }

			 Test_PredictedValues <-  Mean + (testGeno_Imp %*% Effects)
			 Test_SortedPredictedValues <- sort.int(Test_PredictedValues,decreasing=TRUE,index.return=TRUE) 
			 Test_StrainID <- rownames(TestData_Table_Num_Filt)
		
		}
	  }
	  
    if(GPModel == "BayesB (bWGR)"){ 
	 
		 pred <- emBB(trainPheno,trainGeno_Imp) 
		 Mean <- as.numeric(pred$mu)
		 Effects <- pred$b
		 
		 PredictedValues <- Mean + (trainGeno_Imp %*% Effects)
		 SortedPredictedValues <- sort.int(PredictedValues,decreasing=TRUE,index.return=TRUE)
		 
		 
		if(!is.null(testGeno_Imp)){
		 
			 while(anyNA(testGeno_Imp)){
			   testGeno_Imp_Mod <- testGeno_Imp
			   trainGeno_Imp_Mod <- trainGeno_Imp
			   Effects_Mod <- Effects
			   testNA <- apply(testGeno_Imp_Mod,2,function(x) length(which(is.na(x))))
			   testNAIndices <-(which(unlist(testNA) !=0))
			   testGeno_Imp <-testGeno_Imp_Mod[,-testNAIndices]
			   Effects <- Effects_Mod[-testNAIndices]
			   trainGeno_Imp <- trainGeno_Imp_Mod[,-testNAIndices]
			 }
			 
			 Test_PredictedValues <-  Mean + (testGeno_Imp %*% Effects)
			 Test_SortedPredictedValues <- sort.int(Test_PredictedValues,decreasing=TRUE,index.return=TRUE) 
			 Test_StrainID <- rownames(TestData_Table_Num_Filt)
		 
		}
     }
	  
    if(GPModel == "BayesLASSO (bWGR)"){ 
	 
		 pred <- emBL(trainPheno,trainGeno_Imp) 
		 Mean <- as.numeric(pred$mu)
		 Effects <- pred$b
		 PredictedValues <- Mean + (trainGeno_Imp %*% Effects)
		 SortedPredictedValues <- sort.int(PredictedValues,decreasing=TRUE,index.return=TRUE)
		 
		 
		 if(!is.null(testGeno_Imp)){
		 
			 while(anyNA(testGeno_Imp)){
			   testGeno_Imp_Mod <- testGeno_Imp
			   trainGeno_Imp_Mod <- trainGeno_Imp
			   Effects_Mod <- Effects
			   testNA <- apply(testGeno_Imp_Mod,2,function(x) length(which(is.na(x))))
			   testNAIndices <- (which(unlist(testNA) !=0))
			   testGeno_Imp <- testGeno_Imp_Mod[,-testNAIndices]
			   Effects <- Effects_Mod[-testNAIndices]
			   trainGeno_Imp <- trainGeno_Imp_Mod[,-testNAIndices]
			 }
			 
			 Test_PredictedValues <-  Mean + (testGeno_Imp %*% Effects)
			 Test_SortedPredictedValues <- sort.int(Test_PredictedValues,decreasing=TRUE,index.return=TRUE) 
			 Test_StrainID <- rownames(TestData_Table_Num_Filt)
		 }
    }

### Train Output Table
	
	Train_StrainID <- rownames(trainGeno_Imp)
	Train_SortedStrainID <- Train_StrainID[SortedPredictedValues[[2]]]
	Train_SortedPredValues <- SortedPredictedValues[[1]]
	TrainOutputDF <- cbind.data.frame(Train_SortedStrainID,Train_SortedPredValues)
	colnames(TrainOutputDF) <- c("LineID",paste("Predicted Value for ",trait,sep=""))


### Upper bound of Reliability
   
	# trainGeno_Imp2 <- apply(trainGeno_Imp,2,function(x) as.numeric(x)+1)
	# rownames(trainGeno_Imp2) <- rownames(trainGeno_Imp)
	# print(dim(trainGeno_Imp2))
	# print(length(trainPheno))
   	
# function to remove repeated genotypes


	cleanData <- cleanREPV2(trainPheno,trainGeno_Imp)
	
	M <-  cleanData[[2]]
    M.Pdt <- t(M)%*% solve(M %*% t(M) + diag(nrow(M))) %*% M
      
   

### 
	
	if(!is.null(testGeno_Imp)){
		  
		U <- apply(testGeno_Imp,1,function(x) getU(M.Pdt,x)) 
			
		Test_SortedIndices <- Test_SortedPredictedValues[[2]]
		Sorted_Test_StrainID <- Test_StrainID[Test_SortedIndices]
		U_Sorted <- U[(Test_SortedIndices)]
			
		TrgtoutputDF <- cbind.data.frame(Sorted_Test_StrainID ,round(Test_SortedPredictedValues[[1]],digits=2),round(U_Sorted,digits=5))
		colnames(TrgtoutputDF) <- c("LineID",paste("Predicted Value for ",trait,sep=""),"Upper Bound of Reliability")
		outputDF_List <- list(TrainOutputDF,TrgtoutputDF)
	}else if(is.null(testGeno_Imp)){
	
	    TrgtoutputDF <- "Target set has not been defined. Select training set to view from the side panel." 
		outputDF_List <- list(TrainOutputDF,TrgtoutputDF)
	}
   
      
    return(outputDF_List)
 }
 