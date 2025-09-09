#' Generate an Optimal Training Subset Using Genetic Algorithm
#'
#' The `getOptimalTS` function identifies an optimal subset of training samples
#' for genomic prediction using a genetic algorithm (GA). It reduces or keeps candidate
#' data, imputes genotypic data, computes principal components (PCs), and selects
#' the optimal subset based on specified parameters.
#'
#' @param Data_Table_Num_Filt_List A list containing two data frames:
#'   - The first data frame (`TrainData_Table_Num_Filt`) represents the training data.
#'   - The second data frame (`TestData_Table_Num_Filt`) represents the test data.
#' @param trait A character vector specifying the trait(s) to analyze. Can be a single trait or multiple traits.
#' @param nTraits An integer specifying the number of traits in the dataset.
#' @param noCandidates An integer specifying the number of candidates to consider for training.
#' @param nTrainToSelect An integer specifying the number of training samples to select optimally.
#' @param GAParameters A list containing parameters for the genetic algorithm:
#'   - `errorstat`: A single or multiple error statistics to optimize.
#'   - `InitPop`: Initial population for the GA.
#'   - `npop`: Number of individuals in the population.
#'   - `nelite`: Number of elite individuals to retain.
#'   - `mutprob`: Mutation probability.
#'   - `mutintensity`: Intensity of mutation.
#'   - `niterations`: Number of iterations for the GA.
#'   - `minitbefstop`: Minimum iterations before stopping.
#'   - `tabu`: Whether to use tabu search.
#'   - `tabumemsize`: Size of the tabu memory.
#'   - `plotiters`: Whether to plot iterations.
#'   - `lambda`: Regularization parameter.
#'   - `mc.cores`: Number of cores for parallel execution.
#'
#' @return A list containing:
#'   - `Train_STPGA`: Results of the genetic algorithm, including evaluated solutions.
#'   - `Train_STPGA_Ids`: Character vector of IDs for the selected optimal training samples.
#'   - `Train_STPGA_Indices`: Integer vector of indices corresponding to the selected samples.
#'
#' @details
#' - The function processes training and test data to compute genotypic data.
#' - Genotypic data is imputed using the `snpQC` function.
#' - Principal components are calculated from the genotypic data.
#' - A genetic algorithm is applied to select an optimal training subset based on user-specified parameters.
#' - Supports both single-objective and multi-objective optimization.
#'
#' @examples
#' \dontrun{
#' # Example input data
#' TrainData <- data.frame(
#'   ss1 = c(1, 0, 1),
#'   ss2 = c(0, 1, 0),
#'   trait1 = c(5.2, 6.1, 4.8)
#' )
#' TestData <- data.frame(
#'   ss1 = c(1, 1),
#'   ss2 = c(0, 0),
#'   trait1 = c(5.9, 4.6)
#' )
#'
#' Data_Table_Num_Filt_List <- list(TrainData, TestData)
#'
#' GAParams <- list(
#'   errorstat = "MSE",
#'   InitPop = NULL,
#'   npop = 50,
#'   nelite = 5,
#'   mutprob = 0.01,
#'   mutintensity = 1,
#'   niterations = 100,
#'   minitbefstop = 10,
#'   tabu = FALSE,
#'   tabumemsize = 5,
#'   plotiters = TRUE,
#'   lambda = 0.5,
#'   mc.cores = 2
#' )
#'
#' # Generate an optimal training subset
#' result <- getOptimalTS(
#'   Data_Table_Num_Filt_List,
#'   trait = "trait1",
#'   nTraits = 1,
#'   noCandidates = 3,
#'   nTrainToSelect = 2,
#'   GAParameters = GAParams
#' )
#'
#' print(result$Train_STPGA_Ids)
#' }
#'
#' @importFrom NAM snpQC
#' @import STPGA
#' @export

getOptimalTS <- function(Data_Table_Num_Filt_List,trait,nTraits,noCandidates,nTrainToSelect,GAParameters){
  
    TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
	
	initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	finalCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[length(grep("ss",colnames(TrainData_Table_Num_Filt)))]
	 
	initColTst <- grep("ss",colnames(TestData_Table_Num_Filt))[1]
	finalColTst <- grep("ss",colnames(TestData_Table_Num_Filt))[length(grep("ss",colnames(TestData_Table_Num_Filt)))]
	 
	
	trainPheno <- TrainData_Table_Num_Filt[,trait]
	trainGeno_012 <- apply(TrainData_Table_Num_Filt[,c(initCol:finalCol)],2,as.numeric)
	testGeno_012 <- apply(TestData_Table_Num_Filt[,c(initColTst:finalColTst)],2,as.numeric)
	
	
	
## Complete Genotypes table
	
	geno_012 <- rbind(trainGeno_012,testGeno_012)
    totGeno <-  apply(geno_012,2,function(x) x-1)
	rownames(totGeno) <- c(rownames(TrainData_Table_Num_Filt),rownames(TestData_Table_Num_Filt))
	if(length(trait)==1){
	 names(trainPheno) <- rownames(TrainData_Table_Num_Filt)
    }
	if(length(trait)>1){
	 rownames(trainPheno) <- rownames(TrainData_Table_Num_Filt)
    }
	
    totCandidates <-  rownames(totGeno)
	
    #Test <- testIds
	
	Test <- rownames(TestData_Table_Num_Filt)
	testIndices <- which(totCandidates %in% Test)
 
### Set1 with 99 PCs, 100 candidates , 100 test and 50 nTrain 
    set.seed(125)
	
## Reduce or keep candidate genotypic data 
    
    nTrain <- nrow(TrainData_Table_Num_Filt)
	
	 if(noCandidates!= nTrain){
	   trainSetIndices <-sample(c(1:nTrain),noCandidates)
	 }
	 if(noCandidates== nTrain){
	   trainSetIndices <- c(1:nTrain)
	 }
	
	G <- rbind(totGeno[trainSetIndices,],totGeno[testIndices,])
	rownames(G) <- c(rownames(totGeno)[trainSetIndices],rownames(totGeno)[testIndices])
	Candidates <- rownames(totGeno)[trainSetIndices]
	
	G_Imp <- snpQC(G,remove=FALSE,impute=TRUE)
	rownames(G_Imp) <- rownames(G)
    GenoSVD <- svd(G_Imp,nu=99,nv=99)
    PC99 <- G_Imp%*%GenoSVD$v
    rownames(PC99)<-rownames(G_Imp)
    Train_STPGA <- c()
	
    system.time({
	
	  if(length(GAParameters$errorstat)==1){
         Train_STPGA <- GenAlgForSubsetSelection(P=PC99,Candidates,Test,ntoselect=nTrainToSelect,npop=GAParameters$npop, 
		 nelite=GAParameters$nelite, mutprob=GAParameters$mutprob, mutintensity = GAParameters$mutintensity,
         niterations=GAParameters$niterations,minitbefstop=GAParameters$minitbefstop, tabu=GAParameters$tabu,
         tabumemsize = GAParameters$tabumemsize,plotiters=GAParameters$plotiters,errorstat=GAParameters$errorstat,
		 lambda=GAParameters$lambda,InitPop=GAParameters$InitPop, mc.cores=GAParameters$mc.cores)
		 
	   }else if(length(GAParameters$errorstat)>1){
         Train_STPGA <- GenAlgForSubsetSelectionMO(Pcs=PC99,Candidates,Test,ntoselect=nTrainToSelect,InitPop=GAParameters$InitPop,
         npopGA=GAParameters$npop, nelite=GAParameters$nelite, mutprob=GAParameters$mutprob, mutintensity = GAParameters$mutintensity,
         niterations=GAParameters$niterations,minitbefstop=GAParameters$minitbefstop, tabu=GAParameters$tabu,
         tabumemsize = GAParameters$tabumemsize,plotiters=GAParameters$plotiters,errorstat=GAParameters$errorstat,
		 lambda=GAParameters$lambda, mc.cores=GAParameters$mc.cores)
		}
	  
    })
	
	Train_STPGA_Ids <- as.character(Train_STPGA$`Solution with rank 1`)
	Train_STPGA_Indices <-  which(Candidates %in% Train_STPGA_Ids)
    return(list(Train_STPGA,Train_STPGA_Ids,Train_STPGA_Indices))
 }
 
 

#' Generate a Random Training Subset
#'
#' The `getRandomTS` function selects a random subset of training samples from a 
#' given dataset for use in training genomic models. It supports filtering and 
#' imputing genotypic data for selected candidates.
#'
#' @param Data_Table_Num_Filt_List A list containing two data frames:
#'   - The first data frame (`TrainData_Table_Num_Filt`) represents the training data.
#'   - The second data frame (`TestData_Table_Num_Filt`) represents the test data.
#' @param trait A character vector specifying the trait(s) to analyze. Can be a single trait or multiple traits.
#' @param nTraits An integer specifying the number of traits in the dataset.
#' @param noCandidates An integer specifying the number of candidates to consider for training.
#' @param nTrainToSelect An integer specifying the number of training samples to randomly select.
#'
#' @return A list with the following components:
#'   - `Train_Random`: A character vector containing the names of the randomly selected training samples.
#'   - `trainRandomIndices`: An integer vector containing the indices of the randomly selected training samples.
#'
#' @details
#' - The function processes the input data to extract genotypic and phenotypic information.
#' - Genotypic data is imputed using the `snpQC` function.
#' - A random subset of training samples is selected from the candidate set.
#' - The number of candidates and training samples can be controlled using `noCandidates` and `nTrainToSelect`.
#'
#' @examples
#' \dontrun{
#' # Example input data
#' TrainData <- data.frame(
#'   ss1 = c(1, 0, 1),
#'   ss2 = c(0, 1, 0),
#'   trait1 = c(5.2, 6.1, 4.8)
#' )
#' TestData <- data.frame(
#'   ss1 = c(1, 1),
#'   ss2 = c(0, 0),
#'   trait1 = c(5.9, 4.6)
#' )
#'
#' Data_Table_Num_Filt_List <- list(TrainData, TestData)
#' 
#' # Generate a random training subset
#' result <- getRandomTS(
#'   Data_Table_Num_Filt_List,
#'   trait = "trait1",
#'   nTraits = 1,
#'   noCandidates = 3,
#'   nTrainToSelect = 2
#' )
#'
#' print(result$Train_Random)
#' print(result$trainRandomIndices)
#' }
#'
#' @importFrom NAM snpQC
#' @import STPGA
#' @export
 
 getRandomTS <-  function(Data_Table_Num_Filt_List,trait,nTraits,noCandidates,nTrainToSelect){ 
 
    TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
	
	trainPheno <- TrainData_Table_Num_Filt[,trait]
	if(length(trait)==1){
	 names(trainPheno) <- rownames(TrainData_Table_Num_Filt)
    }
	if(length(trait)>1){
	 rownames(trainPheno) <- rownames(TrainData_Table_Num_Filt)
    }
	
	initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	finalCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[length(grep("ss",colnames(TrainData_Table_Num_Filt)))]
	 
	initColTst <- grep("ss",colnames(TestData_Table_Num_Filt))[1]
	finalColTst <- grep("ss",colnames(TestData_Table_Num_Filt))[length(grep("ss",colnames(TestData_Table_Num_Filt)))]
	
	trainGeno_012 <- apply(TrainData_Table_Num_Filt[,c(initCol:finalCol)],2,as.numeric)
	testGeno_012 <- apply(TestData_Table_Num_Filt[,c(initColTst:finalColTst)],2,as.numeric)
	
## Complete genotypic table of all candidates

    trainGeno <-  apply(trainGeno_012,2,function(x) x-1)
	rownames(trainGeno) <- rownames(TrainData_Table_Num_Filt)
	
	Test <- rownames(TestData_Table_Num_Filt)
	
### Set1 with 99 PCs, 100 candidates , 100 test and 50 nTrain 
   
	
## Reduce or keep Candidate Genotypic data 
    nTrain <- nrow(TrainData_Table_Num_Filt)
	set.seed(125)
	# trainSetIndices <- c(1:nTrain)
	
	if(noCandidates!= nTrain){
	  trainSetIndices <-sample(c(1:nTrain),noCandidates)
	}
	if(noCandidates== nTrain){
	 trainSetIndices <- c(1:nTrain)
	}
	
### Define candidate training set
	
    G_Train <- trainGeno[trainSetIndices,]
	Candidates<- rownames(trainGeno)[trainSetIndices]
		
## Impute G_Train 

	G_Train_Imp <- snpQC(G_Train,remove=FALSE,impute=TRUE)
	rownames(G_Train_Imp) <- rownames(G_Train)
	
	trainRandomIndices <- sample(c(1:nrow(G_Train_Imp)),nTrainToSelect) 
	Train_Random <- rownames(G_Train_Imp)[trainRandomIndices]
		
   return(list(Train_Random,trainRandomIndices))

  }
  
  
   
#' Title of getTSComparisons
#'
#' Description of what getTSComparisons does.
#'
#' @param Data_Table_Num_Filt_List Description of Data_Table_Num_Filt_List.
#' @param Train_STPGA Description of Train_STPGA.
#' @param Train_Random Description of Train_Random.
#' @param trait Description of trait.
#' @param nTraits Description of nTraits.
#' @param testIds Description of testIds.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getTSComparisons
#' result <- getTSComparisons(...)
#' }

getTSComparisons <- function(Data_Table_Num_Filt_List,Train_STPGA,Train_Random,trait,nTraits,testIds){ 
  
    TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
	
	trainPheno <- TrainData_Table_Num_Filt[,trait]
	names(trainPheno) <- rownames(TrainData_Table_Num_Filt)
	
	initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	finalCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[length(grep("ss",colnames(TrainData_Table_Num_Filt)))]
	 
	initColTst <- grep("ss",colnames(TestData_Table_Num_Filt))[1]
	finalColTst <- grep("ss",colnames(TestData_Table_Num_Filt))[length(grep("ss",colnames(TestData_Table_Num_Filt)))]
	 
   	
	trainGeno_012 <- apply(TrainData_Table_Num_Filt[,c(initCol:finalCol)],2,as.numeric)
	testGeno_012 <- apply(TestData_Table_Num_Filt[,c(initColTst:finalColTst)],2,as.numeric)
	
		
	trainGeno <-  apply(trainGeno_012,2,function(x) x-1)
	testGeno <- apply(testGeno_012,2,function(x) x-1)
	rownames(trainGeno) <- rownames(TrainData_Table_Num_Filt)

    ## TotGeno
	
	geno_012 <- rbind(trainGeno_012,testGeno_012)
    totGeno <-  apply(geno_012,2,function(x) x-1)
	rownames(totGeno) <- c(rownames(TrainData_Table_Num_Filt),rownames(TestData_Table_Num_Filt))
	
    Test <- testIds
    Candidates <- as.character(rownames(trainGeno))
 
   
   ## Train Set for STPGA and Random   
 
	
	  
	 train_STPGA_Ind <- which(Candidates %in% as.character(Train_STPGA[[1]]$`Solution with rank 1`))
	
	
	
     train_Random_Ind <- which(Candidates %in% as.character(Train_Random[[1]]))
     trainGeno_STPGA <- trainGeno[train_STPGA_Ind,]
	 trainGeno_Random <- trainGeno[train_Random_Ind,]
  
	 
	## Impute missing values 
	
	 trainGeno_STPGA_Imp <- snpQC(trainGeno_STPGA,remove=FALSE,impute=TRUE)
	 trainGeno_Random_Imp <- snpQC(trainGeno_Random,remove=FALSE,impute=TRUE)

     trainPheno_STPGA <- trainPheno[train_STPGA_Ind]
     trainPheno_Random <- trainPheno[train_Random_Ind]
             
     pred_STPGA <- emCV(trainPheno_STPGA,trainGeno_STPGA_Imp,5,5)
	 pred_Random <- emCV(trainPheno_Random,trainGeno_Random_Imp,5,5)
		
	 
	 PA_Table <- rbind(pred_STPGA[c("emRR","emBB","emBL")],pred_Random[c("emRR","emBB","emBL")])
 	 colnames(PA_Table) <- c("emRR","emBB","emBL")
	 rownames(PA_Table) <-  c("STPGA Training Set","Random Training Set")
	 
	 return(PA_Table)
  }
  

   
#' Title of getTSComparisonsMT
#'
#' Description of what getTSComparisonsMT does.
#'
#' @param Data_Table_Num_Filt_List Description of Data_Table_Num_Filt_List.
#' @param Train_STPGA Description of Train_STPGA.
#' @param Train_Random Description of Train_Random.
#' @param trait Description of trait.
#' @param nTraits Description of nTraits.
#' @param testIds Description of testIds.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getTSComparisonsMT
#' result <- getTSComparisonsMT(...)
#' }
getTSComparisonsMT <- function(Data_Table_Num_Filt_List,Train_STPGA,Train_Random,trait,nTraits,testIds){ 
  
    TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
		
	
	initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	finalCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[length(grep("ss",colnames(TrainData_Table_Num_Filt)))]
	 
	initColTst <- grep("ss",colnames(TestData_Table_Num_Filt))[1]
	finalColTst <- grep("ss",colnames(TestData_Table_Num_Filt))[length(grep("ss",colnames(TestData_Table_Num_Filt)))]
	 
	
	trainGeno_012 <- apply(TrainData_Table_Num_Filt[,c(initCol:finalCol)],2,as.numeric)
	testGeno_012 <- apply(TestData_Table_Num_Filt[,c(initColTst:finalColTst)],2,as.numeric)
		
	trainGeno <-  apply(trainGeno_012,2,function(x) x-1)
	testGeno <- apply(testGeno_012,2,function(x) x-1)
	rownames(trainGeno) <- rownames(TrainData_Table_Num_Filt)

   ## TotGeno
	
	geno_012 <- rbind(trainGeno_012,testGeno_012)
    totGeno <-  apply(geno_012,2,function(x) x-1)
	rownames(totGeno) <- c(rownames(TrainData_Table_Num_Filt),rownames(TestData_Table_Num_Filt))
	
    Test <- testIds
    Candidates <- as.character(rownames(trainGeno))
 
    ## Train Set for STPGA and Random   
 

	train_STPGA_Ind <- which(Candidates %in% as.character(Train_STPGA[[1]]$`Solution with rank 1`))
	
    train_Random_Ind <- which(Candidates %in% as.character(Train_Random[[1]]))
    trainGeno_STPGA <- trainGeno[train_STPGA_Ind,]
	trainGeno_Random <- trainGeno[train_Random_Ind,]
  
	 
	## Impute missing values 
	
	 trainGeno_STPGA_Imp <- snpQC(trainGeno_STPGA,remove=FALSE,impute=TRUE)
	 trainGeno_Random_Imp <- snpQC(trainGeno_Random,remove=FALSE,impute=TRUE)
	
	### ST Models

  STModels<- FALSE
  if(STModels==TRUE){	
	 
	 PA_TableComb <- c()
	 for(nT in 1:length(trait)){
	 
	  trainPheno <- TrainData_Table_Num_Filt[,trait[nT]]
	  names(trainPheno) <- rownames(TrainData_Table_Num_Filt)
   	
 	  trainPheno_STPGA <- trainPheno[train_STPGA_Ind]
      trainPheno_Random <- trainPheno[train_Random_Ind]
      NAIndices_STPGA <- which(is.na(trainPheno_STPGA))   
      NAIndices_Random<- which(is.na(trainPheno_Random))
	  if(length(NAIndices_STPGA) >= 1){
        pred_STPGA <- emCV(trainPheno_STPGA[-NAIndices_STPGA],trainGeno_STPGA_Imp[-NAIndices_STPGA,],5,5)
	  }
	  if(length(NAIndices_Random) >= 1){
	    pred_Random <- emCV(trainPheno_Random[-NAIndices_Random],trainGeno_Random_Imp[-NAIndices_Random,],5,5)
	  }
	  if(length(NAIndices_STPGA) <1){
        pred_STPGA <- emCV(trainPheno_STPGA,trainGeno_STPGA_Imp,5,5)
	  }
	  if(length(NAIndices_Random) < 1){
	   pred_Random <- emCV(trainPheno_Random,trainGeno_Random_Imp,5,5)
	  }
	  PA_Table <- rbind(rep(trait[nT],3),c("RR","BB","BL"),pred_STPGA[c("emRR","emBB","emBL")],pred_Random[c("emRR","emBB","emBL")])
 	  #colnames(PA_Table) <- c("emRR","emBB","emBL")
	  if(nT==1){
	   PA_Table2 <-  cbind(c("Trait","GPModel","STPGA","Random"),PA_Table)
	  }
	  if(nT>1){
	  
	    PA_Table2 <- PA_Table
	  }
	  
	  PA_TableComb <- cbind(PA_TableComb,PA_Table2)
	 } 
	  
	  PA_Out_Table <- PA_TableComb
    }  

  
	 
###MT Models 

  MTModels <- TRUE

  if(MTModels==TRUE){
	
	 TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	 trainPheno <- TrainData_Table_Num_Filt[,trait]
     rownames(trainPheno) <- rownames(TrainData_Table_Num_Filt)
   	
 	 	 
	 NAIndices <- unlist(apply(TrainData_Table_Num_Filt[,trait],2,function(x)which(is.na(x))))
   
    
	 trainPheno_STPGA <- trainPheno[train_STPGA_Ind,]
     trainPheno_Random <- trainPheno[train_Random_Ind,]
	 trainSet_Pheno <- list(trainPheno_STPGA,trainPheno_Random)
	 trainSet_Geno <-  list(trainGeno_STPGA_Imp,trainGeno_Random_Imp)
	 
	 PA_Out_List <- list()
	 
	 nIter <- 2
	 k<- 2
	 
     for(nTS in 1:2){ 	 
	
	  pheno_wNA<- trainSet_Pheno[[nTS]]
      geno_Imp_wNA <- trainSet_Geno[[nTS]]
	
	  
	  pheno <- pheno_wNA[-NAIndices,]
	  geno_Imp <- geno_Imp_wNA[-NAIndices,]
	    
	  Geno <- rownames(geno_Imp)
	  
	  rownames(geno_Imp) <- Geno
	  rownames(pheno) <- Geno
	  
	  
	  Y <- pheno
	  yNA <- as.matrix(pheno)
	  X <- geno_Imp
	  rownames(X) <- rownames(geno_Imp)
	  n <- nrow(X)
	  p <- ncol(X)
	  	 
		 
	  yNA_DF <- as.data.frame(yNA)
	  yNA_DF$id <- as.factor(rownames(X))
	  
	  PA_TableComb <- c()
	
##############################################################
	 
      GPModelMT_List <- list("BRR (BGLR)","RKHS (BGLR)","Spike-Slab (BGLR)") #,"GBLUP (SOMMER)")
      PA_Out_GP <- list()
  
      n <- nrow(trainPheno_STPGA)

      for(nGP in 1:length(GPModelMT_List)){ 
    
		GPModelMT <- GPModelMT_List[[nGP]]
		
		PA_List <- list()
		
		for(nrep in 2:nIter){
	 
		 nK <- floor(n/k)
		 k_List <- list()
		 set.seed(125+nrep) 
		 tot <- c(1:n)
		  for(i in 1:k){
			set.seed(125+nrep+k+5) 
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
		
		 nTst <- length(testIndices)
		 
		 pheno[testIndices,] <- matrix(rep(NA,nTst*ncol(pheno)),nrow=nTst,ncol=ncol(pheno))
			 
	## Kinship Matrix 
	     
			
		 
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
		 # if(GPModelMT == "GBLUP (SOMMER)"){ 
				 
				# A.Tot <- A.mat(X)
				# rownames(A.Tot) <- rownames(X) 
				# colnames(A.Tot) <- rownames(X)
				# fm3 <- mmer(as.formula(paste("cbind(",paste(trait,collapse=","),")~1",sep="")),
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
	 	 
	    PA[[i]] <- diag(cor(Test_PredictedValues,Y[tst,]))
	 
	}
	  
	   PA_List[[nrep]] <- do.call(rbind,lapply(PA,function(x) x))
	  
	}
	PA_Out_GP[[nGP]] <- apply(do.call(rbind,lapply(PA_List,function(x)x)),2,function(x) mean(x,na.rm=TRUE))
	
  }
  
	PA_Out1 <-  apply(do.call(rbind,lapply(PA_Out_GP,function(x) x)),2,function(x) round(x,digits=2))
	
	
	PA_Out2 <- cbind(unlist(GPModelMT_List),PA_Out1)
	PA_Out <- rbind(c("GPModel",colnames(PA_Out1)),PA_Out2)
    PA_Out_List[[nTS]] <- PA_Out

  }
  
  STPGAOut <- cbind(rep("STPGA",nrow(PA_Out_List[[1]])),PA_Out_List[[1]])
  RandomOut <- cbind(rep("Random",nrow(PA_Out_List[[2]])),PA_Out_List[[2]])
  
  PA_Out_Table1 <- sapply(c(1:nrow(STPGAOut)),function(x) rbind(STPGAOut[,x],RandomOut[,x]))
  PA_Out_Table <- PA_Out_Table1[-1,]
  PA_Out_Table[1,1] <- "TrainSet"
  } 
#########################	 
	 
	return(PA_Out_Table)
 }
 
 
 
  