

 
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





#' Multi-Environment Genomic Prediction for Leave-half-the-Test-Out mode(fitMEModels_KFoldCV_V5)
#'
#' Performs  genomic prediction in multi-environment trials,
#' allowing flexible kernel estimation, environmental modeling, and user-defined
#' cross-validation schemes.
#'
#'
#' @param DT_1_Filt_List List of phenotype data frames after merging with genotype data.
#' @param genoDat_List List of genotype data frames after merging with phenotype data.
#' @param strainGeno Character vector/ Character specifying the strain names.
#' @param KG Matrix or NULL. Genotypic relationship matrix. Default = NULL, in which case
#'   the kernel is estimated internally. A user-defined matrix can be supplied.
#' @param KE Matrix or NULL. Environmental relationship matrix. Default = NULL.
#'   Required if `FitEnvModels = TRUE`.
#' @param tstVar Character. Character specifying the factor that is dropped in leave-one-out-factor 
#' @param tst Integer/Character. Integer/ Character specifying the test/factor level that is removed from the training set. 
#' @param K Integer. Number of folds in k-fold CV
#' @parma NRep Integer. Number of iterations in k-fold CV 
#' @param KMethod Character. Kernel method for estimating genotypic and
#'   genotype × environment kernels. Default = "Linear". Option: "Gaussian".
#' @param Reaction Logical. Whether to fit reaction norm models when KE is not null.
#'   Default = FALSE. If TRUE, a user-provided `KE` must be supplied.
#' @param trait Character. Variable specifying trait value to be fitted.
#'
#'
#' @return A list containing results from prediction, including predicted values,
#'   model fit statistics, variance components and performance metrics..
#' @examples
#' \dontrun{
#' # Example usage of getMEPred
#' result <- getMEPred(...)
#' }


fitMEModels_KFoldCV_V5 <- function(DT,genoDat,strainGeno,KG=NULL,KE=NULL,tstVar,tst,K,NRep,KMethod,Reaction=FALSE,trait){ 
  
  ###
  cor_LOT_CV_List <- list() 
  var_LOT_CV_List <- list() 
  fit_Out_LOT_CV_List <- list() 
  
  if(is.null(KE)){
   Reaction <- FALSE
  }
  
  ### Prepare Data for CVs
  
  DT_2A <- DT
  DT_2B <- DT
    
  k <- K
  nrep <- NRep
  
  DT_tstIndices2 <- which(DT_2A[,tstVar] %in% tst)
  tstIndices2 <- which(DT_2B[,tstVar] %in% tst)
    
  DT_Tst <- DT_2B[tstIndices2,] 
  
  DT_2A_Tst <- DT_Tst
  
  
## Select GIDs for K-fold CV and Pred  in CV2 type of cross validation 
  
  nLns <- round(length(unique(DT_Tst$gid))/k,digits=0)
  set.seed(1250+nrep)
  rmStrains <- sample(unique(DT_Tst$gid),nLns)
  
  NAInd <- which(DT_Tst$gid %in% rmStrains)
  nNAInd <- setdiff(c(1:length(DT_Tst$gid)),NAInd)
  RmInd <- list(NAInd,nNAInd)
   
  if(!is.null(KG)){
    KMethod <- NULL
  }
    
 for(nfold in 1:k){
  
  rmInd <- RmInd[[nfold]]
  DT_2B_Tst <- DT_Tst
  DT_2B_Tst[rmInd,"value"] <- NA
  dim(DT_2B_Tst)
  
  DT_2B_Tst <- droplevels(DT_2B_Tst)
  DT_2B[tstIndices2,] <- DT_2B_Tst
   
  Y <- DT_2B[,c("env","gid","value")]
  y <- "value"
  gid <- "gid"
  env <- "env"
  ne = length(unique(Y$env))
 
  X <- genoDat[which(rownames(genoDat) %in% Y$gid),]
  
  
      
  ### Creating kernel models with get_kernel
  if(ne>1){
      system.time({
        MM = get_kinship_kernels_v3(Y,X,K_G= KG, K_E = KE,intercept.random=FALSE,model="MM",method=KMethod,dimension_KE=NULL,bandwidth = 1,quantil = 0.5,reaction=Reaction)
      })
      
      system.time({
        MDs = get_kinship_kernels_v3(Y,X,K_G = KG,K_E = KE,intercept.random=FALSE,model="MDs",method=KMethod,dimension_KE=NULL,bandwidth = 1,quantil = 0.5,reaction=Reaction)
      })
      
      ### Fit heterogeneous variance using BGGE    
      system.time({
       MDe <- get_kinship_kernels_v3(Y,X,K_G = KG,K_E= KE,intercept.random=FALSE,model="MDe",method=KMethod,dimension_KE=NULL,bandwidth = 1,quantil = 0.5,reaction=Reaction)
        
      })
      
      fixed=model.matrix(~0+env,DT_2B)
      
      system.time({
        fit_MM=kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=MM, fixed=fixed)
      })
      
      
      system.time({
        fit_MDs=kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=MDs, fixed=fixed)
      })
      
      ###
      
      y_MDe <- DT_2B[,"value"]
      NE <- table(DT_2B$env)
      system.time({
        fit_MDe <- BGGE(y_MDe,K=MDe,XF=fixed,ne=NE)
      })
      
      digits <- 3
      VarComp = Vcomp.BGGE(model = fit_MDe, env = Y$env, 
                           gid = Y$gid, digits = digits)
      
      
      fit_MDe_Out = list(yHat = fit_MDe$yHat, varE = fit_MDe$varE, random = fit_MDe$K, 
                            BGGE = fit_MDe, VarComp = VarComp)
							
	
       
	  if(trait=="Maturity.rating"){
	     fit_MM$yHat <- round(fit_MM$yHat,digits=0)
		 fit_MDs$yHat <- round(fit_MDs$yHat,digits=0)
		 fit_MDe_Out$yHat <- round(fit_MDe_Out$yHat,digits=0)
      }
	
 ####   
     
      
      varMM_LOT_CV <- fit_MM$VarComp
      varMDs_LOT_CV <- fit_MDs$VarComp
      varMDe_LOT_CV <- fit_MDe_Out$VarComp
      
      pred_MM_Tst <- fit_MM$yHat[tstIndices2]
	  pred_MDs_Tst <- fit_MDs$yHat[tstIndices2]
	  pred_MDe_Tst <- fit_MDe$yHat[tstIndices2]
	  
	  corMM_LOT_CV <- cor(pred_MM_Tst[rmInd],DT_2A[rmInd,"value"])
      corMDs_LOT_CV <- cor(pred_MDs_Tst[rmInd],DT_2A[rmInd,"value"])
      corMDe_LOT_CV <- cor(pred_MDe_Tst[rmInd],DT_2A[rmInd,"value"])
	  
	
      fit_MM_Out <- cbind.data.frame(DT_2A_Tst[rmInd,"gid"],DT_2A_Tst[rmInd,"env"],DT_2A_Tst[rmInd,"Test"],DT_2A_Tst[rmInd,"LocTest"],DT_2A_Tst[rmInd,"value"],pred_MM_Tst[rmInd])
      fit_MDs_Out <- cbind.data.frame(DT_2A_Tst[rmInd,"gid"],DT_2A_Tst[rmInd,"env"],DT_2A_Tst[rmInd,"Test"],DT_2A_Tst[rmInd,"LocTest"],DT_2A_Tst[rmInd,"value"],pred_MDs_Tst[rmInd])
      fit_MDe_Out <- cbind.data.frame(DT_2A_Tst[rmInd,"gid"],DT_2A_Tst[rmInd,"env"],DT_2A_Tst[rmInd,"Test"],DT_2A_Tst[rmInd,"LocTest"],DT_2A_Tst[rmInd,"value"],pred_MDe_Tst[rmInd])
   
      colnames(fit_MM_Out) <- c("Strain","Loc","Test","LocTest","Obs","Pred")
      colnames(fit_MDs_Out) <- c("Strain","Loc","Test","LocTest","Obs","Pred")
      colnames(fit_MDe_Out) <- c("Strain","Loc","Test","LocTest","Obs","Pred")
	  
	 }else if(ne==1){ 
	    K_E <- NULL
		Reaction <- NULL
	   
	   system.time({
        SM = get_geno_kernel_v2(Y,X,KG,KE,intercept.random=FALSE,model="SM",method=KMethod,dimension_KE=NULL,reaction=Reaction)
       })
	  
	   fixed=model.matrix(~1,DT_2B)
	
	   system.time({
			fit_SM=kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=SM, fixed=fixed)
	   })
	   
	   if(trait=="Maturity.rating"){
	     fit_SM$yHat <- round(fit_SM$yHat,digits=0)
       }
	   
	   
	   varSM_LOT_CV <- fit_SM$VarComp
	   pred_SM_Tst <- fit_SM$yHat[tstIndices2]
	   corSM_LOT_CV <- cor(pred_SM_Tst,DT_2A[rmInd,"value"])
	   
	   fit_SM_Out <- cbind.data.frame(DT_2A_Tst[rmInd,"gid"],DT_2A_Tst[rmInd,"env"],DT_2A_Tst[rmInd,"Test"],DT_2A_Tst[rmInd,"Location"],DT_2A_Tst[rmInd,"value"],pred_SM_Tst[rmInd])        
	   colnames(fit_SM_Out) <- c("Strain","Loc","Test","Location","Obs","Pred")
  
	}
   
  if(ne>1){
		cor_LOT_CV_List[[nfold]] <- list(corMM_LOT_CV,corMDs_LOT_CV,corMDe_LOT_CV)
		var_LOT_CV_List[[nfold]] <- list(varMM_LOT_CV,varMDs_LOT_CV,varMDe_LOT_CV)
		fit_Out_LOT_CV_List[[nfold]] <- list(fit_MM_Out,fit_MDs_Out,fit_MDe_Out)
  }else if(ne==1){
		cor_LOT_CV_List[[nfold]] <- corSM_LOT_CV
		var_LOT_CV_List[[nfold]] <- varSM_LOT_CV
		fit_Out_LOT_CV_List[[nfold]] <- fit_SM_Out
  }
  
 }
  
  return(list(cor_LOT_CV_List,var_LOT_CV_List,fit_Out_LOT_CV_List,RmInd))
  
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





#' Fit multi-environment genomic prediction models to obtain predicted values  
#'
#' @param DT Data frame with phenotypic data.
#' @param genoDat data frame with genotypic data.
#' @param strainGeno Character vector containing strain IDs of genotypes.
#' @param KG Matrix or NULL. Genotypic relationship matrix. Default = NULL, in which case
#'   the kernel is estimated internally. A user-defined matrix can be supplied.
#' @param KE Matrix or NULL. Environmental relationship matrix. Default = NULL.
#'   Required if `FitEnvModels = TRUE`.
#' @param method Character. Kernel method for estimating genotypic and
#'   genotype × environment kernels. Default = "Linear". Option: "Gaussian".
#' @param fitEnvModels Logical. Whether to include environmental covariates in the model.
#'   Default = FALSE. If TRUE, a user-provided `KE` must be supplied.
#' @param FixedTerm= Character. Variable specifying fixed effects in the model.
#' @param IDColsList List with old and modified character vector of ID columns.
#'
#' @return a list containing correlation between training and target sets, variance components and fit values
#' @examples
#' \dontrun{
#' # Example usage of fitMEModels_Predictions
#' result <- fitMEModels_Predictions(...)
#' }

fitMEModels_Predictions <- function(DT,genoDat,strainGeno,KG=NULL,KE=NULL,method,fitEnvModels=FALSE,FixedTerm=Fixed,IDColsList=IDColsMEList){ 
  
### Prepare Data for CV
  
  DT_2A <- DT
  DT_2B <- DT
     
  dim(DT_2B)
  DT_2B <- droplevels(DT_2B)
  
  Y <- DT_2B[,c("env","gid","value")]
  
  y <- "value"
  gid <- "gid"
  env <- "env"
  
  X <- genoDat
  
 ### Here IDCols is the vector with 'gid', 'env' and 'value'
 
  IDCols <- IDColsList$`IDColsMod`
  
 ### IDColsMod is the vector with the original IDs  for output
  IDColsMod <- IDColsList$`IDCols`
  
 ###
  
  if(!is.null(KG)){
    method <- NULL
  }
  
  Fixed <- FixedTerm
  
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
    VarComp = Vcomp.BGGE(model = fit_MDe, env = Y$env, 
                         gid = Y$gid, digits = digits)
    
    fit_MDe_Out = list(yHat = fit_MDe$yHat, varE = fit_MDe$varE, random = fit_MDe$K, 
                       BGGE = fit_MDe, VarComp = VarComp)
    
    corMM_LOT_CV <- cor(fit_MM$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"])
    corMDs_LOT_CV <- cor(fit_MDs$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"])
    corMDe_LOT_CV <- cor(fit_MDe_Out$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"])
    
    varMM_LOT_CV <- fit_MM$VarComp
    varMDs_LOT_CV <- fit_MDs$VarComp
    varMDe_LOT_CV <- fit_MDe_Out$VarComp
   
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
      
      
      corMM_LOT_CV <- cor(fit_MM_GB$yHat,DT_2A[,"value"])
      corMDs_LOT_CV <- cor(fit_MDs_GB$yHat,DT_2A[,"value"])
      corMDe_LOT_CV <- cor(fit_MDe_GB_Out$yHat,DT_2A[,"value"])
      
      varMM_LOT_CV <- fit_MM_GB$VarComp
      varMDs_LOT_CV <- fit_MDs_GB$VarComp
      varMDe_LOT_CV <- fit_MDe_GB_Out$VarComp
          
   
      fit_MM_Out <- cbind.data.frame(DT_2A[,IDCols],DT_2A[,"value"],fit_MM_GB$yHat)
      fit_MDs_Out <- cbind.data.frame(DT_2A[,IDCols],DT_2A[,"value"],fit_MDs_GB$yHat)
      fit_MDe_Out <- cbind.data.frame(DT_2A[,IDCols],DT_2A[,"value"],fit_MDe_GB_Out$yHat)
   
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
      
      corMM_LOT_CV <- cor(fit_MM_GK$yHat,DT_2A[,"value"])
      corMDs_LOT_CV <- cor(fit_MDs_GK$yHat,DT_2A[,"value"])
      corMDe_LOT_CV <- cor(fit_MDe_GK_Out$yHat,DT_2A[,"value"])
      
      varMM_LOT_CV <- fit_MM_GK$VarComp
      varMDs_LOT_CV <- fit_MDs_GK$VarComp
      varMDe_LOT_CV <- fit_MDe_GK_Out$VarComp
      
    	  
	  fit_MM_Out <- cbind.data.frame(DT_2A[,IDCols],DT_2A[,"value"],fit_MM_GK$yHat)
      fit_MDs_Out <- cbind.data.frame(DT_2A[,IDCols],DT_2A[,"value"],fit_MDs_GK$yHat)
      fit_MDe_Out <- cbind.data.frame(DT_2A[,IDCols],DT_2A[,"value"],fit_MDe_GK_Out$yHat)
   
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
    
    
    corMM_LOT_CV <- cor(fit_EMM$yHat,DT_2A[,"value"])
    corMDs_LOT_CV <- cor(fit_EMDs$yHat,DT_2A[,"value"])
    corMDe_LOT_CV <- cor(fit_EMDe$yHat,DT_2A[,"value"])
    
	fit_MM_Out <- fit_EMM
	fit_MDs_Out <- fit_EMDs
	fit_MDe_Out <- fit_EMDe_Out
	
    varMM_LOT_CV <- fit_EMM$VarComp
    varMDs_LOT_CV <- fit_EMDs$VarComp
    varMDe_LOT_CV <- fit_EMDe_Out$VarComp
	
	fit_MM_Out <- cbind.data.frame(DT_2A[,IDCols],DT_2A[,"value"],fit_EMM$yHat)
    fit_MDs_Out <- cbind.data.frame(DT_2A[,IDCols],DT_2A[,"value"],fit_EMDs$yHat)
    fit_MDe_Out <- cbind.data.frame(DT_2A[,IDCols],DT_2A[,"value"],fit_EMDe_Out$yHat)
   
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
      
      
      corMM_LOT_CV <- cor(fit_EMM_GK$yHat,DT_2A[,"value"])
      corMDs_LOT_CV <- cor(fit_EMDs_GK$yHat,DT_2A[,"value"])
      corMDe_LOT_CV <- cor(fit_EMDe_GK_Out$yHat,DT_2A[,"value"])
      
      varMM_LOT_CV <- fit_EMM_GK$VarComp
      varMDs_LOT_CV <- fit_EMDs_GK$VarComp
      varMDe_LOT_CV <- fit_EMDe_GK_Out$VarComp
	  
	  # fit_MM_Out <- fit_EMM_GK
	  # fit_MDs_Out <- fit_EMDs_GK
	  # fit_MDe_Out <- fit_EMDe_GK_Out
	  
	  fit_MM_Out <- cbind.data.frame(DT_2A[,IDCols],DT_2A[,"value"],fit_EMM_GK$yHat)
      fit_MDs_Out <- cbind.data.frame(DT_2A[,IDCols],DT_2A[,"value"],fit_EMDs_GK$yHat)
      fit_MDe_Out <- cbind.data.frame(DT_2A[,IDCols],DT_2A[,"value"],fit_EMDe_GK_Out$yHat)
   
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
      
      
      corMM_LOT_CV <- cor(fit_EMM_GB$yHat,DT_2A[,"value"])
      corMDs_LOT_CV <- cor(fit_EMDs_GB$yHat,DT_2A[,"value"])
      corMDe_LOT_CV <- cor(fit_EMDe_GB_Out$yHat,DT_2A[,"value"])
      
      varMM_LOT_CV <- fit_EMM_GB$VarComp
      varMDs_LOT_CV <- fit_EMDs_GB$VarComp
      varMDe_LOT_CV <- fit_EMDe_GB_Out$VarComp
	  
	  
	  fit_MM_Out <- cbind.data.frame(DT_2A[,IDCols],DT_2A[,"value"],fit_EMM_GB$yHat)
      fit_MDs_Out <- cbind.data.frame(DT_2A[,IDCols],DT_2A[,"value"],fit_EMDs_GB$yHat)
      fit_MDe_Out <- cbind.data.frame(DT_2A[,IDCols],DT_2A[,"value"],fit_EMDe_GB_Out$yHat)
   
      colnames(fit_MM_Out) <- c(IDColsMod,"Obs","Pred")
      colnames(fit_MDs_Out) <- c(IDColsMod,"Obs","Pred")
      colnames(fit_MDe_Out) <- c(IDColsMod,"Obs","Pred")
    
   }
  }  
  
  cor_LOT_CV_List <- list(corMM_LOT_CV,corMDs_LOT_CV,corMDe_LOT_CV)
  var_LOT_CV_List <- list(varMM_LOT_CV,varMDs_LOT_CV,varMDe_LOT_CV)
  fit_Out_LOT_CV_List <- list(fit_MM_Out,fit_MDs_Out,fit_MDe_Out)
  
  return(list(cor_LOT_CV_List,var_LOT_CV_List,fit_Out_LOT_CV_List))
  
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


getME_CV <- function(DT_1_Filt_List,genoDat_List,traits,KG=NULL,KE=NULL,CVMet,factVar,K,NIter,KMethod="Linear",FitEnvModels=FALSE,fixedME=fixME,envVar=varEnv,IDColsME,IDColME,LocME=LocationME,YrME=YearME){

 print("ME_CV_In")
 nTraits <- length(traits)
 UniqID <- IDColME
 
 ME_Out_CV_Trt <- list()
 DT_2_List <- list() 
 IDColsMEList_Out <- list()
 
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
  Loc <- levels(factor(DT_2$Loc))
  
  nanInd <-   which(is.nan(DT_2[,trait]))
  if(length(nanInd)>0){DT_2 <- DT_2[-nanInd,]}
  
  DT_2 <- droplevels(DT_2)
  EnvVar <- envVar
  LineID <- "Strain"
  Trait <- trait
  
#### 
 
### Set Env and Trait Value Columns   

  selColInd <- match(c(EnvVar,LineID,Trait,UniqID),colnames(DT_2))
  selIDColInd <- match(c(EnvVar,LineID,UniqID),IDColsME)
  
  colnames(DT_2)[selColInd] <- c("env","gid","value","uniqID") 
  loc <- levels(factor(DT_2$Loc))
  
  IDColsMEMod <- IDColsME
  IDColsMEMod[selIDColInd] <- c("env","gid","uniqID")
  
  Fixed= fixedME
  ke <- KE

  IDColsMEList <- list(IDColsME,IDColsMEMod)
  names(IDColsMEList) <- c("IDCols","IDColsMod")
	 
	 if(CVMet=="CV_LOFO"){
	 
	    print(paste("Running Leave One ",factVar," Out CV",sep=""))
		 factorsVar <- levels(factor(DT_2[,factVar]))
	 
		 if(FitEnvModels==FALSE){
			ME_LOFCV <-   foreach(nFact=1:length(factorsVar),.export=c("fitMEModels_LOF_CV","get_kernel_MDe","get_geno_kernel","Vcomp.BGGE"),.packages=c("BGGE","EnvRtype")) %dopar%
			   (fitMEModels_LOF_CV(DT_2,genoDat,strainGeno,KG=NULL,KE=NULL,factVar,factorsVar[nFact],method=KMethod,fitEnvModels=FitEnvModels,FixedTerm=Fixed,IDColsList=IDColsMEList))
		 }else if(FitEnvModels==TRUE){ 
			ME_LOFCV <- foreach(nFact=1:length(factorsVar),.export=c("fitMEModels_LOF_CV","get_kernel_MDe","get_geno_kernel","Vcomp.BGGE"),.packages=c("BGGE","EnvRtype")) %dopar% (fitMEModels_LOF_CV(DT_2,genoDat,strainGeno,KG=NULL,KE=ke,factVar,factorsVar[nFact],method=KMethod,fitEnvModels=FitEnvModels,FixedTerm=Fixed,IDColsList=IDColsMEList))
		 }
		
		 ME_Out_CV_Trt[[nTrt]] <- ME_LOFCV
	  
	 }else if(CVMet!= "CV_LOFO"){
		 if(FitEnvModels==FALSE){
			ME_CV <- fitMEModels_CV_fits(DT_2,genoDat,strainGeno,KG=NULL,KE=NULL,CVMet,K,NIter,method=KMethod,fitEnvModels=FitEnvModels,trait = NULL,store_full_fits = FALSE)
		 }else if(FitEnvModels==TRUE){ 
			ME_CV <- fitMEModels_CV_fits(DT_2,genoDat,strainGeno,KG=NULL,KE=ke,CVMet,K,NIter,method=KMethod,fitEnvModels=FitEnvModels,trait = NULL,store_full_fits = FALSE)
		 }
		 ME_Out_CV_Trt[[nTrt]] <- ME_CV
	  }
	  
	  DT_2_List[[nTrt]] <- DT_2
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


getTstIndices_CV_Data <- function(DT_2B, CVMet, nIter = 5, k = 5,
                                  id_col = "gid", env_col = "env", y_col = "value",
                                  base_seed = 125L,
                                  cv00_partial_prop = 0.2) {

 
   
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
    

  split_kfold_indices <- function(items, k) {
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

  stop("Unknown CVMet. Use one of: 'CV1', 'CV2', 'CV0', 'CV00'.")
}



#' Title of split_kfold_indices
#'
#' A helper function to split a vector of items into k-fold for cross validation
#'
#' @param items a vector of genotypes or environments .
#' @param k number of folds.
#' 
#' @return a list of items for each fold. 
#' \dontrun{
#' # Example usage of split_kfold_indices
#' result <- split_kfold_indices(...)
#' }

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
  
  
  
  

  
#' Title of fitMEModels_LOF_CV
#'
#' Description of what fitMEModels_LOF_CV does.
#'
#' @param DT Description of DT.
#' @param genoDat Description of genoDat.
#' @param strainGeno Description of strainGeno.
#' @param KG=NULL Description of KG=NULL.
#' @param KE=NULL Description of KE=NULL.
#' @param factVar Description of factVar.
#' @param factr Description of factr.
#' @param method Description of method.
#' @param fitEnvModels=FALSE Description of fitEnvModels=FALSE.
#' @param FixedTerm=Fixed Description of FixedTerm=Fixed.
#' @param IDColsList=IDColsMEList Description of IDColsList=IDColsMEList.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of fitMEModels_LOF_CV
#' result <- fitMEModels_LOF_CV(...)
#' }

fitMEModels_LOF_CV <- function(DT,genoDat,strainGeno,KG=NULL,KE=NULL,factVar,factr,method,fitEnvModels=FALSE,FixedTerm=Fixed,IDColsList=IDColsMEList){ 
  
### Prepare Data for CV
  
  DT_2A <- DT
  DT_2B <- DT
     
  dim(DT_2B)
  DT_2B <- droplevels(DT_2B)
  
  DT_factIndices2 <- which(DT_2A[,factVar] %in% factr)
  factIndices2 <- which(DT_2B[,factVar] %in% factr)
  
  
  DT_2B[factIndices2 ,"value"] <- NA
  dim(DT_2B)
  DT_2B <- droplevels(DT_2B)
  
  
  Y <- DT_2B[,c("env","gid","value")]
  
  y <- "value"
  gid <- "gid"
  env <- "env"
  
  X <- genoDat
  
 ### Here IDCols is the vector with 'gid', 'env' and 'value'
 
  IDCols <- IDColsList$`IDColsMod`
  
 ### IDColsMod is the vector with the original IDs  for output
  IDColsMod <- IDColsList$`IDCols`
  
 ###
  
  if(!is.null(KG)){
    method <- NULL
  }
  
  Fixed <- FixedTerm
  
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
    VarComp = Vcomp.BGGE(model = fit_MDe, env = Y$env, 
                         gid = Y$gid, digits = digits)
    
    fit_MDe_Out = list(yHat = fit_MDe$yHat, varE = fit_MDe$varE, random = fit_MDe$K, 
                       BGGE = fit_MDe, VarComp = VarComp)
    
    corMM_LOT_CV <- cor(fit_MM$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"])
    corMDs_LOT_CV <- cor(fit_MDs$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"])
    corMDe_LOT_CV <- cor(fit_MDe_Out$yHat[tstIndices2],DT_2A[DT_tstIndices2,"value"])
    
    varMM_LOT_CV <- fit_MM$VarComp
    varMDs_LOT_CV <- fit_MDs$VarComp
    varMDe_LOT_CV <- fit_MDe_Out$VarComp
   
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
      
      
      corMM_LOT_CV <- cor(fit_MM_GB$yHat[factIndices2],DT_2A[DT_factIndices2,"value"],use="pairwise.complete.obs")
      corMDs_LOT_CV <- cor(fit_MDs_GB$yHat[factIndices2],DT_2A[DT_factIndices2,"value"],use="pairwise.complete.obs")
      corMDe_LOT_CV <- cor(fit_MDe_GB_Out$yHat[factIndices2],DT_2A[DT_factIndices2,"value"],use="pairwise.complete.obs")
      
      varMM_LOT_CV <- fit_MM_GB$VarComp
      varMDs_LOT_CV <- fit_MDs_GB$VarComp
      varMDe_LOT_CV <- fit_MDe_GB_Out$VarComp
          
   
   
      fit_MM_Out <- cbind.data.frame(DT_2A[DT_factIndices2,IDCols],DT_2A[DT_factIndices2,"value"],fit_MM_GB$yHat[factIndices2])
      fit_MDs_Out <- cbind.data.frame(DT_2A[DT_factIndices2,IDCols],DT_2A[DT_factIndices2,"value"],fit_MDs_GB$yHat[factIndices2])
      fit_MDe_Out <- cbind.data.frame(DT_2A[DT_factIndices2,IDCols],DT_2A[DT_factIndices2,"value"],fit_MDe_GB_Out$yHat[factIndices2])
   
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
      
      corMM_LOT_CV <- cor(fit_MM_GK$yHat[factIndices2],DT_2A[DT_factIndices2,"value"],use="pairwise.complete.obs")
      corMDs_LOT_CV <- cor(fit_MDs_GK$yHat[factIndices2],DT_2A[DT_factIndices2,"value"],use="pairwise.complete.obs")
      corMDe_LOT_CV <- cor(fit_MDe_GK_Out$yHat[factIndices2],DT_2A[DT_factIndices2,"value"],use="pairwise.complete.obs")
      
      varMM_LOT_CV <- fit_MM_GK$VarComp
      varMDs_LOT_CV <- fit_MDs_GK$VarComp
      varMDe_LOT_CV <- fit_MDe_GK_Out$VarComp
      
	     	  
	  fit_MM_Out <- cbind.data.frame(DT_2A[DT_factIndices2,IDCols],DT_2A[DT_factIndices2,"value"],fit_MM_GK$yHat[factIndices2])
      fit_MDs_Out <- cbind.data.frame(DT_2A[DT_factIndices2,IDCols],DT_2A[DT_factIndices2,"value"],fit_MDs_GK$yHat[factIndices2])
      fit_MDe_Out <- cbind.data.frame(DT_2A[DT_factIndices2,IDCols],DT_2A[DT_factIndices2,"value"],fit_MDe_GK_Out$yHat[factIndices2])
   
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
    
    
    corMM_LOT_CV <- cor(fit_EMM$yHat[factIndices2],DT_2A[DT_factIndices2,"value"])
    corMDs_LOT_CV <- cor(fit_EMDs$yHat[factIndices2],DT_2A[DT_factIndices2,"value"])
    corMDe_LOT_CV <- cor(fit_EMDe$yHat[factIndices2],DT_2A[DT_factIndices2,"value"])
    
	fit_MM_Out <- fit_EMM
	fit_MDs_Out <- fit_EMDs
	fit_MDe_Out <- fit_EMDe_Out
	
    varMM_LOT_CV <- fit_EMM$VarComp
    varMDs_LOT_CV <- fit_EMDs$VarComp
    varMDe_LOT_CV <- fit_EMDe_Out$VarComp
	
	fit_MM_Out <- cbind.data.frame(DT_2A[DT_factIndices2,IDCols],DT_2A[DT_factIndices2,"value"],fit_EMM$yHat[factIndices2])
    fit_MDs_Out <- cbind.data.frame(DT_2A[DT_factIndices2,IDCols],DT_2A[DT_factIndices2,"value"],fit_EMDs$yHat[factIndices2])
    fit_MDe_Out <- cbind.data.frame(DT_2A[DT_factIndices2,IDCols],DT_2A[DT_factIndices2,"value"],fit_EMDe_Out$yHat[factIndices2])
   
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
      
      
      corMM_LOT_CV <- cor(fit_EMM$yHat[factIndices2],DT_2A[DT_factIndices2,"value"],use="pairwise.complete.obs")
      corMDs_LOT_CV <- cor(fit_EMDs$yHat[factIndices2],DT_2A[DT_factIndices2,"value"],use="pairwise.complete.obs")
      corMDe_LOT_CV <- cor(fit_EMDe$yHat[factIndices2],DT_2A[DT_factIndices2,"value"],use="pairwise.complete.obs")
    
	  fit_MM_Out <- fit_EMM
	  fit_MDs_Out <- fit_EMDs
	  fit_MDe_Out <- fit_EMDe_Out
	
      varMM_LOT_CV <- fit_EMM$VarComp
      varMDs_LOT_CV <- fit_EMDs$VarComp
      varMDe_LOT_CV <- fit_EMDe_Out$VarComp
	
	  fit_MM_Out <- cbind.data.frame(DT_2A[DT_factIndices2,IDCols],DT_2A[DT_factIndices2,"value"],fit_EMM$yHat[factIndices2])
      fit_MDs_Out <- cbind.data.frame(DT_2A[DT_factIndices2,IDCols],DT_2A[DT_factIndices2,"value"],fit_EMDs$yHat[factIndices2])
      fit_MDe_Out <- cbind.data.frame(DT_2A[DT_factIndices2,IDCols],DT_2A[DT_factIndices2,"value"],fit_EMDe_Out$yHat[factIndices2])
   
   
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
      
      
       
      corMM_LOT_CV <- cor(fit_EMM$yHat[factIndices2],DT_2A[DT_factIndices2,"value"],use="pairwise.complete.obs")
      corMDs_LOT_CV <- cor(fit_EMDs$yHat[factIndices2],DT_2A[DT_factIndices2,"value"],use="pairwise.complete.obs")
      corMDe_LOT_CV <- cor(fit_EMDe$yHat[factIndices2],DT_2A[DT_factIndices2,"value"],use="pairwise.complete.obs")
    
	  fit_MM_Out <- fit_EMM
	  fit_MDs_Out <- fit_EMDs
	  fit_MDe_Out <- fit_EMDe_Out
	
      varMM_LOT_CV <- fit_EMM$VarComp
      varMDs_LOT_CV <- fit_EMDs$VarComp
      varMDe_LOT_CV <- fit_EMDe_Out$VarComp
	
	  fit_MM_Out <- cbind.data.frame(DT_2A[DT_factIndices2,IDCols],DT_2A[DT_factIndices2,"value"],fit_EMM$yHat[factIndices2])
      fit_MDs_Out <- cbind.data.frame(DT_2A[DT_factIndices2,IDCols],DT_2A[DT_factIndices2,"value"],fit_EMDs$yHat[factIndices2])
      fit_MDe_Out <- cbind.data.frame(DT_2A[DT_factIndices2,IDCols],DT_2A[DT_factIndices2,"value"],fit_EMDe_Out$yHat[factIndices2])
   
      
   
      colnames(fit_MM_Out) <- c(IDColsMod,"Obs","Pred")
      colnames(fit_MDs_Out) <- c(IDColsMod,"Obs","Pred")
      colnames(fit_MDe_Out) <- c(IDColsMod,"Obs","Pred")
    
   }
  }  
  
  cor_LOT_CV_List <- list(corMM_LOT_CV,corMDs_LOT_CV,corMDe_LOT_CV)
  var_LOT_CV_List <- list(varMM_LOT_CV,varMDs_LOT_CV,varMDe_LOT_CV)
  fit_Out_LOT_CV_List <- list(fit_MM_Out,fit_MDs_Out,fit_MDe_Out)
  
  return(list(cor_LOT_CV_List,var_LOT_CV_List,fit_Out_LOT_CV_List))
  
}


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


####

#MM_GB = get_geno_kernel(Y,X,KG,KE,intercept.random=FALSE,model="MM",method="GB",dimension_KE=NULL)

#' Title of get_geno_kernel
#'
#' Description of what get_geno_kernel does.
#'
#' @param Y Description of Y.
#' @param X Description of X.
#' @param KG Description of KG.
#' @param KE Description of KE.
#' @param intercept.random=FALSE Description of intercept.random=FALSE.
#' @param model Description of model.
#' @param method Description of method.
#' @param dimension_KE=NULL Description of dimension_KE=NULL.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of get_geno_kernel
#' result <- get_geno_kernel(...)
#' }

get_geno_kernel <- function(Y,X,KG,KE,intercept.random=FALSE,model,method,dimension_KE=NULL){
  
  ne = length(unique(Y$env))
  Zg <- stats::model.matrix(~0 + gid, Y)
  colnames(Zg) = gsub(colnames(Zg), pattern = "gid", replacement = "")
  ng <- length(unique(Y$gid))
  
  K_G = KG 
  K_E = KE
  
  model_b <- model
  
  if(is.null(KG)){
    if(method=="GB"){
      K= getK(Y = Y,X=X, kernel = "GB", model = model_b, intercept.random = intercept.random) 
    }else if(method=="GK"){
      K = getK(Y = Y, X= X,kernel = "GK", model = model_b, intercept.random = intercept.random) 
    }
  }else if(!is.null(KG)){
    
    K= getK(Y = Y, setKernel = K_G, model = model_b, intercept.random = intercept.random) 
  }
  
  names(K)[grep(names(K),pattern="^G$")] <- paste0("KG_", names(K)[grep("^G$",names(K))])
  names(K)[grep("KG_G$",names(K),invert = T)] <- paste0("KGE_", names(K)[grep("KG_G$",names(K),invert = T)])
  
  obs_GxE = paste0(Y$env, ":", Y$gid) 
  
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
#' @return 
#' @export
#' @examples a list containing various metrics such as predictive ability, RMSE, variance components, fit values and summary of these metrics
#' \dontrun{ 
#' # Examples of summarize_ME_CV_fits
#' result <- 
#' }
 
 summarize_ME_CV_fits <- function(
  fit_obj,          # output from fitMEModels_CV_fits()
  DT,               # original DT used in fitting (for Obs)
  IDColsList = IDColsMEList
){
  # --- helpers ---
  safe_metrics <- function(pred, obs) {
    ok <- is.finite(pred) & is.finite(obs); n <- sum(ok)
    if (!n) return(list(cor = NA_real_, rmse = NA_real_, n_test = 0L))
    list(cor = if (n > 1) stats::cor(pred[ok], obs[ok]) else NA_real_,
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
  # ID columns validate/align
  resolve_idcols <- function(DT, in_cols, out_cols) {
    if (is.null(in_cols) || !length(in_cols)) { in_cols <- c("gid","env"); out_cols <- c("gid","env") }
    if (is.null(out_cols)) out_cols <- in_cols
    if (length(out_cols) != length(in_cols)) {
      n <- min(length(in_cols), length(out_cols))
      in_cols  <- in_cols[seq_len(n)]
      out_cols <- out_cols[seq_len(n)]
    }
    df <- data.frame(col_in=in_cols, col_out=out_cols, stringsAsFactors=FALSE)
    keep <- df$col_in %in% names(DT)
    if (!all(keep)) message("[summarize_ME_CV_fits] Dropping missing ID columns: ",
                            paste(df$col_in[!keep], collapse=", "))
    df <- df[keep, , drop=FALSE]
    if (!nrow(df)) {
      fb <- c("gid","env"); fb <- fb[fb %in% names(DT)]
      df <- data.frame(col_in=fb, col_out=fb, stringsAsFactors=FALSE)
    }
    list(IDCols=df$col_in, IDColsMod=df$col_out)
  }

  ID_in  <- IDColsList$`IDCols`
  ID_out <- IDColsList$`IDColsMod`
  ids <- resolve_idcols(DT, ID_in, ID_out)
  IDCols <- ids$IDCols; IDColsMod <- ids$IDColsMod

  # --- accumulators ---
  metrics_all <- list(); var_all <- list(); fits_all <- list()

  nReps <- length(fit_obj$results)
  N     <- nrow(DT)
 
  for (irep in seq_len(nReps)) {
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

		  # Coerce to numeric vector and guard length
		  yHat <- as.numeric(yHat)
		  if (!length(yHat)) next

		  # Build safe index:
		  #  - within DT
		  #  - within yHat length
		  #  - finite predictions
		  
		  ok_idx <- tst_idx[tst_idx <= length(yHat)]
		  if (!length(ok_idx)) next
		  ok_idx <- ok_idx[is.finite(yHat[ok_idx])]
		  if (!length(ok_idx)) next

		  # Predictions & obs aligned
		  pred <- yHat[ok_idx]
		  obs  <- DT$value[ok_idx]
		  if (!length(pred) || !length(obs)) next
	
          met_dat <- cbind.data.frame(pred,obs)
		  
		  met_dat_list[[m]] <- met_dat  
		}
		  met_dat_list_nF[[fold]] <- met_dat_list
		  
	}

    met_data_list <- list()
    for(m in names(item$models)){ 
		 
		   met_data <- do.call(rbind.data.frame,lapply(met_dat_list_nF, function(x) x[[m]]))
		   met_data_list[[m]] <- met_data
		 		   
           pred_dat <- met_data$pred
		   obs_dat <- met_data$obs
		   met  <- safe_metrics(pred_dat, obs_dat)

		   metrics_all[[length(metrics_all)+1L]] <-
			data.frame(replicate = irep, model = m,
					   cor = met$cor, rmse = met$rmse, n_test = met$n_test,
					   stringsAsFactors = FALSE)
 
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
			 if (length(IDCols)) names(tab)[seq_along(IDCols)] <- IDColsMod
			 fits_all[[length(fits_all)+1L]] <- tab
		    }
	}
	
   }

  metrics_long <- if (length(metrics_all)) do.call(rbind, metrics_all) else data.frame()
  var_long     <- if (length(var_all))     do.call(rbind, var_all)     else data.frame()
  fits_long    <- if (length(fits_all))    do.call(rbind, fits_all)    else data.frame()

  metrics_summary <- if (nrow(metrics_long)) {
    aggregate(cbind(cor, rmse) ~ model, data = metrics_long, function(x)
      c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), n = sum(is.finite(x))))
  } else data.frame()

 
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

# --- use it like this ---
metrics_summary <- if (nrow(metrics_long)) {
  agg <- aggregate(cbind(cor, rmse) ~ model, data = metrics_long, function(x)
    c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), n = sum(is.finite(x))))
  build_metrics_summary(agg)
} else data.frame()


  list(
    metrics_long    = metrics_long,
    var_long        = var_long,
    fits_long       = fits_long,
    metrics_summary = metrics_summary,
    cv              = fit_obj$cv
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
  KG = NULL, KE = NULL,
  CVMet, k, nIter,
  method,
  fitEnvModels = FALSE,
  trait = NULL,                # e.g., "Maturity.rating" to round predictions
  store_full_fits = FALSE      # set TRUE only if you really need raw fit objects (memory-heavy)
){
  requireNamespace("foreach", quietly = TRUE)
  print("MECV Fit - in")

  DT_2A <- DT
  DT_2B <- droplevels(DT)

  
  # CV splits
  Dat_Out <- getTstIndices_CV_Data(DT_2B, CVMet, nIter, k)
  nReps   <- length(Dat_Out$DT_masked)
  N       <- nrow(DT_2A)

  # If KG is provided, ignore GB/GK method
  if (!is.null(KG)) method <- NULL

  # helper to fit one fold and return predictions + VC (not raw fits unless asked)
  .fit_fold <- function(DT_fold, KG, KE, method, fitEnvModels, genoDat, trait) {
    Y   <- DT_fold[, c("env","gid","value")]
    y   <- "value"; gid <- "gid"; env <- "env"
    ne  <- length(unique(Y$env))
    X   <- if (!is.null(rownames(genoDat))) genoDat[rownames(genoDat) %in% Y$gid, , drop = FALSE] else genoDat

    fitted <- list(); raw_store <- list()

    if (!fitEnvModels && is.null(KE) && !is.null(KG) && is.null(method)) {
      if (ne > 1) {
        MM  <- get_kernel(K_G = KG, y = y, gid = gid, env = env, data = DT_fold, model = "MM")
        MDs <- get_kernel(K_G = KG, y = y, gid = gid, env = env, data = DT_fold, model = "MDs")
        fixed <- model.matrix(~0 + env, DT_fold)
        fit_MM  <- kernel_model(y = y, env = env, gid = gid, data = DT_fold, random = MM,  fixed = fixed)
        fit_MDs <- kernel_model(y = y, env = env, gid = gid, data = DT_fold, random = MDs, fixed = fixed)
        y_MDe <- DT_fold[["value"]]; NE <- table(DT_fold$env)
        MDe  <- get_kernel_MDe(Y, KG, KE, intercept.random = FALSE, dimension_KE = NULL)
        fit_MDe <- BGGE::BGGE(y_MDe, K = MDe, XF = fixed, ne = NE)
        VC_MDe  <- Vcomp.BGGE(model = fit_MDe, env = Y$env, gid = Y$gid, digits = 3)

        fitted$MM  <- list(yHat = fit_MM$yHat,   VarComp = fit_MM$VarComp)
        fitted$MDs <- list(yHat = fit_MDs$yHat,  VarComp = fit_MDs$VarComp)
        fitted$MDe <- list(yHat = fit_MDe$yHat,  VarComp = VC_MDe)

        if (store_full_fits) raw_store <- list(MM = fit_MM, MDs = fit_MDs, MDe = fit_MDe)

      } else {
        SM    <- get_kernel(K_G = KG, y = y, gid = gid, env = env, data = DT_fold, model = "SM")
        fixed <- model.matrix(~1, DT_fold)
        fit_SM <- kernel_model(y = y, env = env, gid = gid, data = DT_fold, random = SM, fixed = fixed)
        fitted$SM <- list(yHat = fit_SM$yHat, VarComp = fit_SM$VarComp)
        if (store_full_fits) raw_store <- list(SM = fit_SM)
      }

    } else if (!fitEnvModels && is.null(KE) && is.null(KG) && !is.null(method)) {
      stopifnot(method %in% c("GB","GK"))
      if (ne > 1) {
        MMg  <- get_geno_kernel(Y, X, KG, KE, intercept.random = FALSE, model = "MM",  method = method, dimension_KE = NULL)
        MDsg <- get_geno_kernel(Y, X, KG, KE, intercept.random = FALSE, model = "MDs", method = method, dimension_KE = NULL)
        MDeg <- get_geno_kernel(Y, X, KG, KE, intercept.random = FALSE, model = "MDe", method = method, dimension_KE = NULL)
        fixed <- model.matrix(~0 + env, DT_fold)
        fit_MM  <- kernel_model(y = y, env = env, gid = gid, data = DT_fold, random = MMg,  fixed = fixed)
        fit_MDs <- kernel_model(y = y, env = env, gid = gid, data = DT_fold, random = MDsg, fixed = fixed)
        y_MDe <- DT_fold[["value"]]; NE <- table(DT_fold$env)
        fit_MDe <- BGGE::BGGE(y_MDe, K = MDeg, XF = fixed, ne = NE)
        VC_MDe  <- Vcomp.BGGE(model = fit_MDe, env = Y$env, gid = Y$gid, digits = 3)

        fitted$MM  <- list(yHat = fit_MM$yHat,   VarComp = fit_MM$VarComp)
        fitted$MDs <- list(yHat = fit_MDs$yHat,  VarComp = fit_MDs$VarComp)
        fitted$MDe <- list(yHat = fit_MDe$yHat,  VarComp = VC_MDe)

        if (store_full_fits) raw_store <- list(MM = fit_MM, MDs = fit_MDs, MDe = fit_MDe)

      } else {
        SMg   <- get_geno_kernel(Y, X, KG, KE, intercept.random = FALSE, model = "SM", method = method, dimension_KE = NULL)
        fixed <- model.matrix(~1, DT_fold)
        fit_SM <- kernel_model(y = y, env = env, gid = gid, data = DT_fold, random = SMg, fixed = fixed)
        fitted$SM <- list(yHat = fit_SM$yHat, VarComp = fit_SM$VarComp)
        if (store_full_fits) raw_store <- list(SM = fit_SM)
      }

    } else if (fitEnvModels && !is.null(KE) && !is.null(KG) && is.null(method)) {
      EMM  <- get_kernel(K_G = KG, K_E = KE, y = y, gid = gid, env = env, data = DT_fold, model = "EMM",  ne = nrow(KE))
      EMDs <- get_kernel(K_G = KG, K_E = KE, y = y, gid = gid, env = env, data = DT_fold, model = "EMDs")
      EMDe <- get_kernel_MDe(Y, KG, KE, intercept.random = FALSE, dimension_KE = NULL)
      fixed <- if (ne > 1) model.matrix(~0 + env, DT_fold) else model.matrix(~1, DT_fold)
      fit_EMM  <- kernel_model(y = y, env = env, gid = gid, data = DT_fold, random = EMM,  fixed = fixed)
      fit_EMDs <- kernel_model(y = y, env = env, gid = gid, data = DT_fold, random = EMDs, fixed = fixed)
      yM <- DT_fold[["value"]]; NE <- table(DT_fold$env)
      fit_EMDe <- BGGE::BGGE(yM, K = EMDe, XF = fixed, ne = NE)
      VC_EMDe  <- Vcomp.BGGE(model = fit_EMDe, env = Y$env, gid = Y$gid, digits = 3)

      fitted$MM  <- list(yHat = fit_EMM$yHat,   VarComp = fit_EMM$VarComp)
      fitted$MDs <- list(yHat = fit_EMDs$yHat,  VarComp = fit_EMDs$VarComp)
      fitted$MDe <- list(yHat = fit_EMDe$yHat,  VarComp = VC_EMDe)

      if (store_full_fits) raw_store <- list(MM = fit_EMM, MDs = fit_EMDs, MDe = fit_EMDe)
    }

    # Optional rounding for ordinal trait
    if (identical(trait, "Maturity.rating")) {
      for (m in names(fitted)) fitted[[m]]$yHat <- round(fitted[[m]]$yHat)
    }

    list(models = fitted, raw = if (store_full_fits) raw_store else NULL)
  }

  # foreach over reps (folds inside)
  # rep_results <- foreach::foreach(
    # irep = seq_len(nReps),
    # .combine = 'c',
    # .multicombine = TRUE,
    # .errorhandling = "pass",
    # .packages = c("BGGE")
  # ) %dopar% {
  
  rep_results <- list() 
  
  for(irep in seq_len(nReps)){
    nFolds <- length(Dat_Out$DT_masked[[irep]])
    fold_list <- vector("list", nFolds)

    for (fold in seq_len(nFolds)) {
      res <- try({
        DT_fold <- Dat_Out$DT_masked[[irep]][[fold]]
        tst_idx <- Dat_Out$test_idx [[irep]][[fold]]
        # keep original test indices (for later scoring)
        out_fit <- .fit_fold(droplevels(DT_fold), KG, KE, method, fitEnvModels, genoDat, trait)
        list(status = "ok", models = out_fit$models, raw = out_fit$raw, tst_idx = as.integer(tst_idx))
      }, silent = TRUE)

      if (inherits(res, "try-error")) {
        fold_list[[fold]] <- list(status = "error", error = as.character(res), models = NULL, raw = NULL, tst_idx = integer(0))
      } else {
        fold_list[[fold]] <- res
      }
    }

    
	rep_results[[irep]] <- list(list(folds = fold_list))
  }
 
  rep_results_c <- lapply(rep_results,function(x) do.call(c,lapply(x,function(y) y)))
  # Shape:
  # results[[rep]]$folds[[fold]]$models[["MM"|"MDs"|"MDe"| "SM"]] = list(yHat=..., VarComp=...)
  # results[[rep]]$folds[[fold]]$tst_idx = integer vector (original row indices)
  # results[[rep]]$folds[[fold]]$status / $error
  list(
    cv      = Dat_Out$cv,
    results = rep_results_c,
    nrow    = N
  )
}





###

#' Title of fitMEModels_CV
#'
#' Description of what fitMEModels_CV does.
#'
#' @param DT Description of DT.
#' @param genoDat Description of genoDat.
#' @param strainGeno Description of strainGeno.
#' @param KG=NULL Description of KG=NULL.
#' @param KE=NULL Description of KE=NULL.
#' @param CVMet Description of CVMet.
#' @param k Description of k.
#' @param nIter Description of nIter.
#' @param method Description of method.
#' @param fitEnvModels=FALSE Description of fitEnvModels=FALSE.
#' @param FixedTerm=Fixed Description of FixedTerm=Fixed.
#' @param IDColsList=IDColsMEList Description of IDColsList=IDColsMEList.
#' @importFrom BGGE BGGE
#' @importFrom EnvRtype get_kernel kernel_model
#' @return .
#' @examples
#' \dontrun{
#' # Example usage of fitMEModels_CV
#' result <- fitMEModels_CV(...)
#' }

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
    
		## Creating kernel models with get_kernel
		
		system.time({
		  MM = get_kernel(K_G = KG, y=y, gid = gid, env = env, data = DT_2B, model = "MM")
		})
		
		system.time({
		  MDs = get_kernel(K_G = KG, y=y, gid = gid, env = env, data = DT_2B, model = "MDs")
		})
		
		### Fit heterogeneous variance using BGGE    
		
		MDe <- get_kernel_MDe(Y,KG,KE,intercept.random=FALSE,dimension_KE=NULL)

		fixed= model.matrix(~0+env,DT_2B)
		
		system.time({
		  fit_MM= kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=MM, fixed=fixed)
		})
		
		
		system.time({
		  fit_MDs= kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=MDs, fixed=fixed)
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
		  ## Creating kernel models with get_kernel
		  
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
			fit_MM_GB=kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=MM_GB, fixed=fixed)
		  })
		  
		  
		  system.time({
			fit_MDs_GB=kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=MDs_GB, fixed=fixed)
		  })
		  
		  
		  ###Fit heterogeneous variance using BGGE    
		  
		  y_MDe <- DT_2B[,"value"]
		  NE <- table(DT_2B$env)
		  system.time({
			fit_MDe_GB <- BGGE(y_MDe,K=MDe_GB,XF=fixed,ne=NE)
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
      ## Creating kernel models with get_kernel
      
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
			fit_MM_GK=kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=MM_GK, fixed=fixed)
		  })
		  
		  
		  system.time({
			fit_MDs_GK=kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=MDs_GK, fixed=fixed)
		  })
		  
		  ### MDe
		  
		  y_MDe <- DT_2B[,"value"]
		  NE <- table(DT_2B$env)
		  system.time({
			fit_MDe_GK <- BGGE(y_MDe,K=MDe_GK,XF=fixed,ne=NE)
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
    
    EMM = get_kernel(K_G = KG, K_E = KE, y= y, gid = gid, env = env, data = DT_2B,model = "EMM",ne=nrow(KE))
    
    EMDs = get_kernel(K_G = KG, K_E = KE, y=y, gid = gid, env = env, data = DT_2B,model = "EMDs")
    
    EMDe <- get_kernel_MDe(Y,KG,KE,intercept.random=FALSE,dimension_KE=NULL)
    
    fixed=model.matrix(~0+env,DT_2B)
    
    system.time({
      fit_EMM =kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=EMM, fixed=fixed)
    })
    
    system.time({
      fit_EMDs=kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=EMDs, fixed=fixed)
    })
    
    y <- DT_2B[,"value"]
    NE <- table(DT_2B$env)
    
    system.time({
      fit_EMDe <- BGGE(y,K=EMDe,XF=fixed,ne=NE)
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
        fit_EMM_GK =kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=EMM_GK, fixed=fixed)
      })
      
      system.time({
        fit_EMDs_GK =kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=EMDs_GK, fixed=fixed)
      })
      
      y <- DT_2B[,"value"]
      NE <- table(DT_2B$env)
      
      system.time({
        fit_EMDe_GK <- BGGE(y,K=EMDe_GK,XF=fixed,ne=NE)
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
        fit_EMM_GB =kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=EMM_GB, fixed=fixed)
      })
      
      system.time({
        fit_EMDs_GB =kernel_model(y=y, env=env, gid=gid, data=DT_2B, random=EMDs_GB, fixed=fixed)
      })
      
      y <- DT_2B[,"value"]
      NE <- table(DT_2B$env)
      
      system.time({
        fit_EMDe_GB <- BGGE(y,K=EMDe_GB,XF=fixed,ne=NE)
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






# ================================================================
# get_kinship_kernel_v3
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
get_kinship_kernels_v3 <- function(
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
      base <- sub("^EM", "", m); rn <- FALSE
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
      Kg <- tcrossprod(X) / ncol(X)                                # gid x gid
      Gobs <- Zg %*% Kg %*% t(Zg)
      rownames(Gobs) <- colnames(Gobs) <- obs
      G_list[["KG_G"]] <- list(Kernel = Gobs, Type = "D")
    } else {  # GK
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
    E_all <- Ze %*% t(Ze)  # == tcrossprod(Ze)
    for (nm in names(G_list)) {
      OUT[[paste0("KG_GE_", sub("^KG_G_?", "", nm))]] <-
        list(Kernel = G_list[[nm]]$Kernel * E_all, Type = "BD")
    }
  } else if (base_model == "MDe") {
    for (nm in names(G_list)) {
      Gobs <- G_list[[nm]]$Kernel
      for (i in seq_len(nEnv)) {
        Bi <- Ze[, i, drop = FALSE] %*% t(Ze[, i, drop = FALSE])   # per-env block
        OUT[[paste0("KG_GE_", env_lv[i], "_", sub("^KG_G_?", "", nm))]] <-
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
          OUT[[paste0("KGE_", gkey, "_", e)]] <- list(Kernel = GE, Type = "D")
        }
      }
    }
    # RNMDe: per-environment blocks (Type = "BD")
    if (base_model == "MDe") {
      for (g in names(G_list)) {
        Gobs <- G_list[[g]]$Kernel
        gkey <- sub("^KG_G_?", "", g)
        for (e in names(E_obs)) {
          GE_full <- Gobs * E_obs[[e]]
          for (i in seq_len(nEnv)) {
            Bi <- Ze[, i, drop = FALSE] %*% t(Ze[, i, drop = FALSE])
            OUT[[paste0("KGE_", gkey, "_", e, "_", env_lv[i])]] <-
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
