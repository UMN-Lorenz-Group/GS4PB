
#' Title of getImputedData_LDKNNI
#'
#' Description of what getImputedData_LDKNNI does.
#'
#' @param FiltGeno Description of FiltGeno.
#' @param l Description of l.
#' @param k Description of k.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getImputedData_LDKNNI
#' result <- getImputedData_LDKNNI(...)
#' }
getImputedData_LDKNNI <- function(FiltGeno,l,k){ 

   tasGenoImp <- rTASSEL::imputeLDKNNi(FiltGeno, highLDSSites = l, knnTaxa = k, maxDistance = 1e+07)
  
  return(tasGenoImp)
} 

#' Title of getImputedData_Num
#'
#' Description of what getImputedData_Num does.
#'
#' @param FiltGeno Description of FiltGeno.
#' @param nN Description of nN.
#' @param Dist Description of Dist.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getImputedData_Num
#' result <- getImputedData_Num(...)
#' }
getImputedData_Num <- function(FiltGeno,nN,Dist){
  
  
   tasGenoImp <- rTASSEL::imputeNumeric(FiltGeno,byMean = TRUE,nearestNeighbors = nN, distance = Dist)
  
  return(tasGenoImp)
}


#' Title of getGenoData_API
#'
#' Description of what getGenoData_API does.
#'
#' @param FiltGeno Description of FiltGeno.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getGenoData_API
#' result <- getGenoData_API(...)
#' }
getGenoData_API <- function(FiltGeno){
  
    filtGenoMat <- as.matrix(FiltGeno)
    filtGenoMat[is.na(filtGenoMat)] <- 9
    currentGenoMat <- cbind.data.frame(rownames(filtGenoMat),filtGenoMat)

    tableReport <- rJava::new(
    rJava::J("net.maizegenetics.dna.map.PositionListTableReport"),
    FiltGeno %>% rTASSEL:::getPositionList()) %>% 
    rTASSEL:::tableReportToDF() %>% as.data.frame()
      	
	varSplit <- strsplit(tableReport[,"VARIANT"],"/")
	varSplitTab <- cbind.data.frame(unlist(lapply(varSplit,function(x) x[1])),unlist(lapply(varSplit,function(x) x[2])))

    vcfIDTab <- cbind.data.frame(tableReport[,c("Name","Chromosome","Position")],varSplitTab)
	colnames(vcfIDTab) <- c("SNPID","Chrom","Position","REF","ALT")
	

    # print(dim(vcfIDTab))
	# print(is.data.frame(vcfIDTab))
    write.table(vcfIDTab,"currentVCFIDTab.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
	
    write.table(currentGenoMat,"current_GenoTable.genotypes",quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)
 
   return(vcfIDTab)
}

#####

#' Title of getImpGenoData_API
#'
#' Description of what getImpGenoData_API does.
#'
#' @param vcfIDTab Description of vcfIDTab.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getImpGenoData_API
#' result <- getImpGenoData_API(...)
#' }
getImpGenoData_API <- function(vcfIDTab){
    
	# vcfIDTab_DF <- as.data.frame(vcfIDTab)
    genoImp <-  read.table("imputed_out.genotypes",sep=" ",header=FALSE)
	vcfIDTab_DF <- as.data.frame(read.table("currentVCFIDTab.txt",sep="\t",header=TRUE))
	
	## Add check to verify that the dimensions match and send a message to the app.   
	genoImpDF <- as.data.frame(genoImp[,-1])
    rownames(genoImpDF) <- genoImp[,1]
	colnames(genoImpDF) <- as.vector(unlist(vcfIDTab_DF[,"SNPID"]))
    genoImpDFT <- as.data.frame(t(genoImpDF))
	genoImpDFT$SNPID <- as.vector(unlist(vcfIDTab_DF[,"SNPID"]))
	# print(dim(genoImpDFT))
    genoImpDFOut <- merge(vcfIDTab_DF,genoImpDFT,by="SNPID")
	genoImpDFOut_Sort <- genoImpDFOut[order(genoImpDFOut$Chrom,genoImpDFOut$Position),]
	
	
	# print(dim(genoImpDFOut))
    return(as_tibble(genoImpDFOut_Sort))
}

######
