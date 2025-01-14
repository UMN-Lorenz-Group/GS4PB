
#' Title of getMEPred
#'
#' Description of what getMEPred does.
#'
#' @param DT_1_Filt_List Description of DT_1_Filt_List.
#' @param genoDat_List Description of genoDat_List.
#' @param traits Description of traits.
#' @param KG=NULL Description of KG=NULL.
#' @param KE=NULL Description of KE=NULL.
#' @param KMethod="Linear" Description of KMethod="Linear".
#' @param FitEnvModels=FALSE Description of FitEnvModels=FALSE.
#' @param fixedME=fixME Description of fixedME=fixME.
#' @param envVar=varEnv Description of envVar=varEnv.
#' @param IDColsME Description of IDColsME.
#' @param LocME=LocationME Description of LocME=LocationME.
#' @param YrME=YearME Description of YrME=YearME.
#'
#' @return Description of what is returned.
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


 
 #' Title of fitMEModels_Predictions
#'
#' Description of what fitMEModels_Predictions does.
#'
#' @param DT Description of DT.
#' @param genoDat Description of genoDat.
#' @param strainGeno Description of strainGeno.
#' @param KG=NULL Description of KG=NULL.
#' @param KE=NULL Description of KE=NULL.
#' @param method Description of method.
#' @param fitEnvModels=FALSE Description of fitEnvModels=FALSE.
#' @param FixedTerm=Fixed Description of FixedTerm=Fixed.
#' @param IDColsList=IDColsMEList Description of IDColsList=IDColsMEList.
#'
#' @return Description of what is returned.
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

######

#######

#' Title of getME_CV
#'
#' Description of what getME_CV does.
#'
#' @param DT_1_Filt_List Description of DT_1_Filt_List.
#' @param genoDat_List Description of genoDat_List.
#' @param traits Description of traits.
#' @param KG=NULL Description of KG=NULL.
#' @param KE=NULL Description of KE=NULL.
#' @param CVMet Description of CVMet.
#' @param factVar Description of factVar.
#' @param K Description of K.
#' @param NIter Description of NIter.
#' @param KMethod="Linear" Description of KMethod="Linear".
#' @param FitEnvModels=FALSE Description of FitEnvModels=FALSE.
#' @param fixedME=fixME Description of fixedME=fixME.
#' @param envVar=varEnv Description of envVar=varEnv.
#' @param IDColsME Description of IDColsME.
#' @param IDColME Description of IDColME.
#' @param LocME=LocationME Description of LocME=LocationME.
#' @param YrME=YearME Description of YrME=YearME.
#'
#' @return Description of what is returned.
#'
#' @importFrom foreach foreach
#' @examples
#' \dontrun{
#' # Example usage of getME_CV
#' result <- getME_CV(...)
#' }
getME_CV <- function(DT_1_Filt_List,genoDat_List,traits,KG=NULL,KE=NULL,CVMet,factVar,K,NIter,KMethod="Linear",FitEnvModels=FALSE,fixedME=fixME,envVar=varEnv,IDColsME,IDColME,LocME=LocationME,YrME=YearME){

 print("ME_CV_In")
 nTraits <- length(traits)
 ME_Out_CV_Trt <- list()
 
 UniqID <- IDColME
 
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
  EnvVar <- envVar
  LineID <- "Strain"
  Trait <- trait

  
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
			ME_CV <- fitMEModels_CV(DT_2,genoDat,strainGeno,KG=NULL,KE=NULL,CVMet,K,NIter,method=KMethod,fitEnvModels=FitEnvModels,FixedTerm=Fixed,IDColsList=IDColsMEList)
		 }else if(FitEnvModels==TRUE){ 
			ME_CV <- fitMEModels_CV(DT_2,genoDat,strainGeno,KG=NULL,KE=ke,CVMet,K,NIter,method=KMethod,fitEnvModels=FitEnvModels,FixedTerm=Fixed,IDColsList=IDColsMEList)
		 }
		 ME_Out_CV_Trt[[nTrt]] <- ME_CV
	  }
  }
  
  return(ME_Out_CV_Trt)
}

### 


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
#' @return Description of what is returned.
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




#' Title of getTstIndices_CV_Data
#'
#' Description of what getTstIndices_CV_Data does.
#'
#' @param DT_2B Description of DT_2B.
#' @param CVMet Description of CVMet.
#' @param nIter Description of nIter.
#' @param k Description of k.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getTstIndices_CV_Data
#' result <- getTstIndices_CV_Data(...)
#' }
getTstIndices_CV_Data <- function(DT_2B,CVMet,nIter,k){
  
# CV1
 
	if(CVMet =="CV1"){
	
   ## A fraction of the strains are removed from all locations, so n corresponds the total no of unique strains 
	
		  n <- length(unique(DT_2B[,"gid"]))
		  
		  trIndices_List_Rep <- list()
		  tstIndices_List_Rep <- list()
		  DT_2B_List_Rep <- list()
		 
		  unqStrains <- unique(DT_2B[,"gid"])
		   
		   for(nrep in 1:nIter){
			 
			
			 k_List <- list()
			 trIndices_List <- list()
			 tstIndices_List <- list()
			 DT_2B_List <- list()
			 
			 
			 nK <- floor(n/k)
			 set.seed(125+nrep) 
			 tot <- c(1:n)
			
			 Test_Pred <- list()
			 Y_Tst <- list() 
			 
			 trn_List <- list()
			 
			 
			 
			 for(nF in 1:k){
			   k_List[[nF]] <- sample(tot,nK)
			   trn_List[[nF]] <- setdiff(tot,k_List[[nF]])
			 }
			 
			 for(nF in 1:k){ 
			   trIndices_List[[nF]] <- unlist(trn_List[[nF]])
			   tstIndices_List[[nF]] <- k_List[[nF]]
			   
			   testIndices <- tstIndices_List[[nF]]
			   testStrains <- unqStrains[testIndices]
			   DT_2B_tst <- DT_2B
			   tstInd_in_DT <-  which(DT_2B_tst$Strain %in% testStrains)
			   DT_2B_tst[tstInd_in_DT,"value"] <- NA
			   DT_2B_List[[nF]] <- DT_2B_tst
			   
			 }	  
		   
		   
			 trIndices_List_Rep[[nrep]] <- trIndices_List
			 tstIndices_List_Rep[[nrep]] <- tstIndices_List
			 DT_2B_List_Rep[[nrep]] <- DT_2B_List
			 
		   }
	   
	    outList <- list(DT_2B_List_Rep,trIndices_List_Rep,tstIndices_List_Rep)
		
		 
	  }else if(CVMet=="CV2"){
	   
		   n <- length(unique(DT_2B[,"uniqID"]))
		   trIndices_List_Rep <- list()
		   tstIndices_List_Rep <- list()
		   DT_2B_List_Rep <- list()
			  
		   for(nrep in 1:nIter){
			 
			 nK <- floor(n/k)
			 k_List <- list()
			 tot <- c(1:n)
			 set.seed(125+nrep) 
			 Test_Pred <- list()
			 Y_Tst <- list() 
			 
			 DT_2B_List <- list()
			 trIndices_List <- list()
			 tstIndices_List <- list()
			 
			 trn_List <- list()
			 
			 for(nF in 1:k){
			   k_List[[nF]] <- sample(tot,nK)
			   trn_List[[nF]] <- setdiff(tot,k_List[[nF]])
			 }
			
			 
			 for(nF in 1:k){ 
			   trIndices_List[[nF]] <- unlist(trn_List[[nF]])
			   tstIndices_List[[nF]] <- k_List[[nF]]
			   
			   tstIndices <- tstIndices_List[[nF]]
			  
			   DT_2B_tst <- DT_2B
			  
			   DT_2B_tst[tstIndices,"value"] <- NA
			   DT_2B_List[[nF]] <- DT_2B_tst
			   
			 }  
		   
		   
			 trIndices_List_Rep[[nrep]] <- trIndices_List
			 tstIndices_List_Rep[[nrep]] <- tstIndices_List
			 DT_2B_List_Rep[[nrep]] <- DT_2B_List
			 
		   }
		   
			outList <- list(DT_2B_List_Rep,trIndices_List_Rep,tstIndices_List_Rep)
		
	  }else if(CVMet=="CV0"){
	  
	  #GID_Env_Tab <-  table(DT_1_Filt_List[[1]][,"env"],DT_1_Filt_List[[1]][,"gid"])
	   
		   GID_Env_Tab <-  table(DT_2B[,"env"],DT_2B[,"gid"])
		   Tot_Env_GId <- apply(GID_Env_Tab,2,sum)
		   
		   fractEnv <- table(Tot_Env_GId)/sum(table(Tot_Env_GId))
		   fractEnvSel <- fractEnv[which(fractEnv>=0.2)]
		
		   nTotEnv <- as.numeric(names(fractEnvSel))
		   nTotEnvSel <- nTotEnv[nTotEnv>=2]
		  
		   selStrainIndices <- which(Tot_Env_GId %in% nTotEnvSel)
		   selGID_Tab <- GID_Env_Tab[,selStrainIndices]
		 
		   selEnvList <- lapply(c(1:ncol(selGID_Tab)),function(x) which(selGID_Tab[,x] ==1))
		   selEnvTab <- do.call(rbind,lapply(selEnvList,function(x)names(x)))
		  
		   LocCombns <- apply(selEnvTab,1,function(x) paste(x,sep="-",collapse="-"))
		   LocCombsTab <- table(LocCombns)
		   LocCombList <- lapply(names(LocCombsTab),function(x) strsplit(x,"-")[[1]])
		   
		   nLoc_in_Comb <- lapply(LocCombList,length)
		   names(nLoc_in_Comb) <- names(LocCombsTab)
		   
	### max no.of env levels where genotypes are present 
		  
		   k <-  max(unlist(nLoc_in_Comb))
			  
		  
		   trIndices_List_Rep <- list()
		   tstIndices_List_Rep <- list()
		   DT_2B_List_Rep <- list()

           nrep <- 1			  
			  
		   trIndices_List <- list()
		   tstIndices_List <- list()
		   DT_2B_List <- list()	   

		   for(nF in 1:k){
			 
				locSubS <- unlist(lapply(LocCombList,function(x) unlist(x)[nF]))
				tstIndices <- which(DT_2B[,"env"] %in% locSubS)
				trnIndices <- which(!DT_2B[,"env"] %in% locSubS)
					
				tstIndices_List[[nF]] <- tstIndices
				trIndices_List[[nF]] <- trnIndices
				
				DT_2B_tst <- DT_2B
			    DT_2B_tst[tstIndices,"value"] <- NA
				DT_2B_List[[nF]] <- DT_2B_tst

			}
			  
			trIndices_List_Rep[[nrep]] <- trIndices_List
			tstIndices_List_Rep[[nrep]] <- tstIndices_List
			DT_2B_List_Rep[[nrep]] <- DT_2B_List
			
			outList <- list(DT_2B_List_Rep,trIndices_List_Rep,tstIndices_List_Rep)
		  	  
	  }else if(CVMet=="CV00"){
	  
	       GID_Env_Tab <-  table(DT_2B[,"env"],DT_2B[,"gid"])
		   Tot_Env_GId <- apply(GID_Env_Tab,2,sum)
		   
		   fractEnv <- table(Tot_Env_GId)/sum(table(Tot_Env_GId))
		   fractEnvSel <- fractEnv[which(fractEnv>=0.2)]
		
		   nTotEnv <- as.numeric(names(fractEnvSel))
		   nTotEnvSel <- nTotEnv[nTotEnv>=2]
		  
		   selStrainIndices <- which(Tot_Env_GId %in% nTotEnvSel)
		   selGID_Tab <- GID_Env_Tab[,selStrainIndices]
		 
		   selEnvList <- lapply(c(1:ncol(selGID_Tab)),function(x) which(selGID_Tab[,x] ==1))
		   selEnvTab <- do.call(rbind,lapply(selEnvList,function(x)names(x)))
		  
		   LocCombns <- apply(selEnvTab,1,function(x) paste(x,sep="-",collapse="-"))
		   LocCombsTab <- table(LocCombns)
		   LocCombList <- lapply(names(LocCombsTab),function(x) strsplit(x,"-"))
		   
		   nLoc_in_Comb <- lapply(LocCombList,length)
		   names(nLoc_in_Comb) <- names(LocCombsTab)
		   
	### max no.of env levels where genotypes are present 
		  
		   k <-  max(unlist(nLoc_in_Comb))  
			  
		  
		   trIndices_List_Rep <- list()
		   tstIndices_List_Rep <- list()
		   DT_2B_List_Rep <- list()
        


           for(nrep in 1:nIter){			 
			   
			   set.seed(125+nrep)
			   
			   trIndices_List <- list()
			   tstIndices_List <- list()
			   DT_2B_List <- list()	   

			   for(nF in 1:k){
				 
				   
					totRmLoc <- unlist(lapply(LocCombList,function(x) unlist(x)[nF]))
					partRmLoc <- unlist(lapply(LocCombList,function(x) unlist(x)[-nF]))
										
					testIndicesTot <- which(DT_2B[,"env"] %in% totRmLoc)
					testIndicesPart0 <- which(DT_2B[,"env"] %in% partRmLoc)
					testIndicesPart <- sample(testIndicesPart0,round(length(testIndicesPart0)/2,digits=0))
					
					testIndices <- c(testIndicesTot,testIndicesPart)
					trnIndices <- setdiff(c(1:nrow(DT_2B)),testIndices)
					
					tstIndices_List[[nF]] <- testIndices
					trIndices_List[[nF]] <- trnIndices
												
					DT_2B_tst <- DT_2B
					  
					DT_2B_tst[testIndices,"value"] <- NA
					DT_2B_List[[nF]] <- DT_2B_tst

				}
				  
				trIndices_List_Rep[[nrep]] <- trIndices_List
				tstIndices_List_Rep[[nrep]] <- tstIndices_List
				DT_2B_List_Rep[[nrep]] <- DT_2B_List
			}
			outList <- list(DT_2B_List_Rep,trIndices_List_Rep,tstIndices_List_Rep)
	  }
	  
   return(outList)
 
  }
  

#### #foreach(nTst=1:nTsts,.packages=c("BGGE","EnvRtype")) %dopar%

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
#' @return Description of what is returned.
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


