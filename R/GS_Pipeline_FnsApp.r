  
###########################################################################################  
#' Title of VCFtoDF
#'
#' Description of what VCFtoDF does.
#'
#' @param infile Description of infile.
#'
#' @importFrom vcfR read.vcfR extract.gt getFIX
#' @importFrom dplyr select rename left_join as_tibble mutate
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of VCFtoDF
#' result <- VCFtoDF(...) 
#' }

VCFtoDF <- function(infile){
  
  vcf <- read.vcfR(infile, verbose = FALSE)
  gt <- extract.gt(vcf, element = "GT", as.numeric = F)
  fix_T <- as_tibble(getFIX(vcf))
  
  gt2 <- matrix(0, ncol = ncol(gt), nrow = nrow(gt))
  colnames(gt2) <- colnames(gt)
  
  gt2a <- apply(gt,2, function(x) gsub("1/1","BB",x))
  gt2b <- gsub("0[/|]0","AA",gt2a)
  gt2c <- gsub("[10][/|][10]","AB",gt2b)
  gt2d <- gsub("\\.[/|]\\.","NA",gt2c)
  
    
  gt2 <- as_tibble(gt2d) %>%
    mutate(SNPID = fix_T$ID)
  
  gt.simple <- fix_T %>% select("ID", "CHROM", POS, REF, ALT) %>%
    rename(SNPID=ID, Chrom=CHROM, Position=POS) %>%
    left_join(gt2, by = 'SNPID')
  
  return(gt.simple)
}
 

######

#' Title of getProcessedData_NUST_withFilters
#'
#' Description of what getProcessedData_NUST_withFilters does.
#'
#' @param gt2d Description of gt2d.
#' @param NUST_BLUEs Description of NUST_BLUEs.
#' @param NUST_Meta_Table Description of NUST_Meta_Table.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getProcessedData_NUST_withFilters
#' result <- getProcessedData_NUST_withFilters(...)
#' }

getProcessedData_NUST_withFilters <- function(gt2d,NUST_BLUEs,NUST_Meta_Table){

	NAIndices <- apply(gt2d[,6:ncol(gt2d)],2,function(x) which(is.na(x)))
	b <- unlist(lapply(NAIndices,length))
	Non_ZeroIndices <- which(b!=0)
	markerID <- as.vector(unlist(gt2d[,1]))



	NUST_Genotypes_VCF_gt2d <- gt2d[,-c(1:5)]
	NUST_Genotypes_VCF <- NUST_Genotypes_VCF_gt2d[,-Non_ZeroIndices]
	NUST_Meta_Table_StrainID <- NUST_Meta_Table[,1]
	NUST_Genotypes_VCF_ID <- colnames(NUST_Genotypes_VCF)


### Remove special characters from strain ID in meta table and genotype table

	NUST_Meta_Table_StrainID_Mod <- gsub("-","",NUST_Meta_Table_StrainID) 
	NUST_Meta_Table_StrainID_Mod1 <- gsub("\\.","",NUST_Meta_Table_StrainID_Mod) 
	NUST_Meta_Table_StrainID_Mod2 <- gsub("\\(","",NUST_Meta_Table_StrainID_Mod1) 
	NUST_Meta_Table_StrainID_Mod3 <- gsub("\\)","",NUST_Meta_Table_StrainID_Mod2) 
	NUST_Meta_Table_StrainID_Mod4 <- gsub("\\_","",NUST_Meta_Table_StrainID_Mod3) 


	NUST_Genotypes_VCF_ID <- gsub("-","",NUST_Genotypes_VCF_ID) 
	NUST_Genotypes_VCF_ID1 <- gsub("\\.","",NUST_Genotypes_VCF_ID) 
	NUST_Genotypes_VCF_ID2 <- gsub("\\(","",NUST_Genotypes_VCF_ID1) 
	NUST_Genotypes_VCF_ID3 <- gsub("\\)","",NUST_Genotypes_VCF_ID2) 
	NUST_Genotypes_VCF_ID4 <- gsub("\\_","",NUST_Genotypes_VCF_ID3)

### Checks 1

	length(which(NUST_Meta_Table_StrainID_Mod4 %in%  NUST_Genotypes_VCF_ID4))
	length(NUST_Genotypes_VCF_ID4)
	length(unique(NUST_Meta_Table_StrainID_Mod4))
	length(unique(NUST_Genotypes_VCF_ID4))

### Remove duplicated IDs from NUST genotype table

	NUST_Meta_Table_Mod <- NUST_Meta_Table
	NUST_Meta_Table_Mod[,1] <- NUST_Meta_Table_StrainID_Mod4
	NUST_Genotypes_Table_Mod <- cbind(NUST_Genotypes_VCF_ID4,t(NUST_Genotypes_VCF)) 

	colnames(NUST_Genotypes_Table_Mod)[1] <- "Strain" 
	colnames(NUST_Meta_Table_Mod)[1] <- "Strain"

	NUST_Meta_Table_Mod2 <- NUST_Meta_Table_Mod[which(!duplicated(NUST_Meta_Table_StrainID_Mod4)),]

### Checks 2 (Equal length)

	length(NUST_Meta_Table_Mod2[,1]) 
	length(unique(NUST_Meta_Table_Mod2[,1]))

## NUST merged table comprising meta data and genotypes table  


	NUST_Merged_Table <- merge(NUST_Meta_Table_Mod2,NUST_Genotypes_Table_Mod,by="Strain")
	init <- ncol(NUST_Meta_Table_Mod2)+1
	final <- ncol(NUST_Merged_Table)
	colnames(NUST_Merged_Table)[init:final] <- markerID

	########## IDs that are persent in genotypes table and not in the merged table 
	diffIDs <- setdiff(NUST_Genotypes_Table_Mod[,1],NUST_Merged_Table[,1])
	diffIndices <- which(NUST_Genotypes_Table_Mod[,1] %in% diffIDs)

	### Numeric coded genotype table from merged data table 


	NUST_Genotypes_Table_Mod_Merged <- NUST_Merged_Table[,c(1,init:final)]
	NUST_Genotypes_Table_Mod_Num1 <- apply(NUST_Genotypes_Table_Mod_Merged[,-1],2,function(x) gsub("BB","-1",x)) 
	NUST_Genotypes_Table_Mod_Num2 <- apply(NUST_Genotypes_Table_Mod_Num1,2,function(x) gsub("AB","0",x)) 
	NUST_Genotypes_Table_Mod_Num3 <- apply(NUST_Genotypes_Table_Mod_Num2,2,function(x) gsub("AA","1",x)) 
	NUST_Genotypes_Table_Mod_Num <- apply(NUST_Genotypes_Table_Mod_Num3,2,function(x) as.numeric(x)+1)

	### Merged Numeric Table

	NUST_Genotypes_Table_Mod_Num_Comb <- cbind(NUST_Genotypes_Table_Mod_Merged[,1],NUST_Genotypes_Table_Mod_Num)
	colnames(NUST_Genotypes_Table_Mod_Num_Comb)[1] <- "Strain" 

	NUST_Merged_Table_Num <- merge(NUST_Meta_Table_Mod2,NUST_Genotypes_Table_Mod_Num_Comb,by="Strain")
	dim(NUST_Merged_Table_Num)
	colnames(NUST_Merged_Table_Num)
	
	
# Filter Genotype Table based on year and Prog_Breeder_Factor

	dupIndices <- which(duplicated(as.character(NUST_Merged_Table[,1])))
	dupStrain <- as.character(NUST_Merged_Table[dupIndices,1])

	dupIndices_in_Table <- which(as.character(NUST_Merged_Table[,1]) %in% dupStrain)

	NUST_Merged_Table_noDup <- NUST_Merged_Table[-dupIndices,]

	rownames(NUST_Merged_Table_noDup) <- NUST_Merged_Table_noDup[,1]
	year_filt_indices <- (which(as.numeric(as.character(NUST_Merged_Table_noDup[,7])) >= 2010))
	NUST_Merged_Table_Filt1 <- NUST_Merged_Table_noDup[year_filt_indices,] 

	dim(NUST_Merged_Table_noDup)
	dim(NUST_Merged_Table_Filt1)


	dupIndices_Num <- which(duplicated(as.character(NUST_Merged_Table_Num[,1])))
	dupStrain_Num <- as.character(NUST_Merged_Table_Num[dupIndices,1])
	dupIndices_in_Table_Num <- which(as.character(NUST_Merged_Table_Num[,1]) %in% dupStrain_Num)

	NUST_Merged_Table_Num_noDup <- NUST_Merged_Table_Num[-dupIndices_Num,]
	rownames(NUST_Merged_Table_Num_noDup) <- NUST_Merged_Table_Num_noDup[,1]
	year_filt_indices <- (which(as.numeric(as.character(NUST_Merged_Table_Num_noDup[,7])) >= 2010))
	NUST_Merged_Table_Num_Filt1 <- NUST_Merged_Table_Num_noDup[year_filt_indices,] 


### Apply Filter(min program size 50 to reduce the number of programs to 12) 
####

	minProgramSize <- 50
	BrProg_col <- 9
	NUST_Merged_Table_noDup[,BrProg_col] <- gsub("_","-",NUST_Merged_Table_noDup[,BrProg_col])
	BrProg_factor <- factor(NUST_Merged_Table[,BrProg_col]) #,labels=c(1:length(levels(factor(NUST_Merged_Table_Filt1[,9])))))
	BrProg_factor_rm <- which(table(BrProg_factor)< minProgramSize)
	BrProg_filt_name <- names(table(BrProg_factor))[-BrProg_factor_rm]
	BrProg_factor_filt_indices <- which(as.character(BrProg_factor) %in% BrProg_filt_name)

	BrProg_factor_filt <- as.character(BrProg_factor)[BrProg_factor_filt_indices]
	BrProg_factor_filt_table <- table(as.character(BrProg_factor)[BrProg_factor_filt_indices])
	length(BrProg_factor_filt_table)
	BrProg_filt <- BrProg_factor[BrProg_factor_filt_indices]



	MG_col <- 10 
	MG_factor <- factor(NUST_Merged_Table_noDup[,MG_col])
	MG_filt <- MG_factor[BrProg_factor_filt_indices]

### Filtered Genotype Table 

	init <- ncol(NUST_Meta_Table_Mod2)+1
	final <- ncol(NUST_Merged_Table_noDup)
	NUST_Genotypes_Table_Mod_Filt1 <- NUST_Merged_Table_noDup[BrProg_factor_filt_indices,init:final]


    NUST_Genotypes_Table_Mod_Num_Filt1 <- NUST_Merged_Table_Num[BrProg_factor_filt_indices,init:final]

    rownames(NUST_Genotypes_Table_Mod_Filt1) <- rownames(NUST_Merged_Table_noDup[BrProg_factor_filt_indices,init:final])
    rownames(NUST_Genotypes_Table_Mod_Num_Filt1) <- rownames(NUST_Merged_Table_Num_noDup[BrProg_factor_filt_indices,init:final])

### Filter for markers with genetic map information

	NUST_Genotypes_Table_Mod_Filt <- as.matrix(NUST_Genotypes_Table_Mod_Filt1)
	NUST_Genotypes_Table_Mod_Num_Filt <- apply(as.matrix(NUST_Genotypes_Table_Mod_Num_Filt1),2,as.numeric)

	rownames(NUST_Genotypes_Table_Mod_Num_Filt) <- rownames(NUST_Genotypes_Table_Mod_Num_Filt1)
	
	dim(NUST_Genotypes_Table_Mod_Filt)
	dim(NUST_Genotypes_Table_Mod_Num_Filt)

######### Data Prep for GP model training

 StrainID_List <- strsplit(as.character(NUST_BLUEs[,1]),"_")

 StrainID <- unlist(lapply(StrainID_List,function(x) paste(x[[1]],x[[2]],sep="")))
 StrainID1 <- gsub("MGIII","",StrainID) 
 StrainID2 <- gsub("MGII","",StrainID1) 
 StrainID3 <- gsub("MGI","",StrainID2) 
 StrainID4 <- gsub("MG0","",StrainID3) 
 StrainID5 <- gsub("MG0","",StrainID4)
 StrainID6 <- gsub(" ","",StrainID5)
 StrainIDMod <- StrainID6
 #rownames(NUST_BLUEs) <-  StrainIDMod 
 length(which(StrainIDMod %in% rownames(NUST_Genotypes_Table_Mod_Filt)))
 
 commonStrainID <- StrainIDMod[(which(StrainIDMod %in% rownames(NUST_Genotypes_Table_Mod_Filt)))]
 
 Diff_StrainID <- setdiff(rownames(NUST_Genotypes_Table_Mod_Filt),commonStrainID)

 NUST_GenoTable_Filtered <- cbind(rownames(NUST_Genotypes_Table_Mod_Filt),NUST_Genotypes_Table_Mod_Filt)
 colnames(NUST_GenoTable_Filtered)[1] <- "StrainID" 
 
 
 NUST_BLUEs_Filt <- NUST_BLUEs[(which(StrainIDMod %in% rownames(NUST_Genotypes_Table_Mod_Filt))),]
 NUST_BLUEs_Filt[,1] <-  commonStrainID  
 colnames(NUST_BLUEs_Filt)[1] <- "StrainID" 
 
 NUST_Data_Table1 <- merge(NUST_GenoTable_Filtered,NUST_BLUEs_Filt,by="StrainID")
 
 dupIndices_Data <- which(duplicated(NUST_Data_Table1[,"StrainID"]))
 NUST_Data_Table <- NUST_Data_Table1[-dupIndices_Data,]
 
 NUST_GenoTable_Num_Filtered <- cbind(rownames(NUST_Genotypes_Table_Mod_Num_Filt),NUST_Genotypes_Table_Mod_Num_Filt)
 colnames(NUST_GenoTable_Num_Filtered)[1] <- "StrainID" 
 
   
####### Filtered Table in Numeric Format 
 
 NUST_Data_Table1_Num <- merge(NUST_GenoTable_Num_Filtered,NUST_BLUEs_Filt,by="StrainID")
 dupIndices_Data <- which(duplicated(NUST_Data_Table1_Num[,"StrainID"]))
 NUST_Data_Table_Num <- NUST_Data_Table1_Num[-dupIndices_Data,]

 NAIndices <- which(is.na(NUST_Data_Table_Num[,"YieldBuA"]))
 NUST_Data_Table_Num_Filt <- NUST_Data_Table_Num[-NAIndices,]
    
	
return(NUST_Data_Table_Num_Filt)

}


#' Calculate the U Statistic
#'
#' The `getU` function calculates the U statistic, which is a measure of the 
#' relative projection of a vector `v` transformed by a matrix `M.Pdt`.
#'
#' @param M.Pdt A square numeric matrix used to transform the vector `v`.
#' @param v A numeric vector to be projected and analyzed.
#'
#' @return A numeric value representing the U statistic, calculated as:
#'   \deqn{U = \frac{(v^\top M.Pdt^\top M.Pdt v)}{(v^\top v)}}
#'   where `v.hat = M.Pdt %*% v`.
#'
#' @details
#' The U statistic is used to assess the proportion of the vector's variance 
#' that is captured by its projection through the matrix `M.Pdt`.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' M <- matrix(c(2, 0, 0, 1), nrow = 2) # Example matrix
#' v <- c(1, 2) # Example vector
#' getU(M, v)
#' }
#'
#' @export
getU <- function(M.Pdt, v) {
  v.hat <- M.Pdt %*% v
  U <- (t(v.hat) %*% v.hat) / (t(v) %*% v)
  return(U)
}


########################

#' Title of getMergedData
#'
#' Description of what getMergedData does.
#'
#' @param gt2d Description of gt2d.
#' @param Pheno Description of Pheno.
#' @param testIDs Description of testIDs.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getMergedData
#' result <- getMergedData(...)
#' }
getMergedData <- function(gt2d, Pheno, testIDs){

    # Extract genotype columns and remove metadata columns (first 5)
    Genotypes_VCF <- gt2d[,-c(1:5)]
    markerID <- as.vector(unlist(gt2d[, 1]))
   
  # colnames(genoImpDF) <- paste("ss",as.vector(as.data.frame(genoImp)[,"SNPID"]),sep="-")
    markerIDMod <-  paste("ss",markerID,sep="-")
  
    # Clean Genotypes IDs
    Genotypes_VCF_ID <- gsub("[-._()]", "", colnames(Genotypes_VCF))
  
    # Clean test IDs
    if (!is.null(testIDs)) { 
        TestIDs <- gsub("[_-]", "", testIDs)
    } else {
        TestIDs <- testIDs
    }
    
    # Check for duplicates
    if (length(unique(Genotypes_VCF_ID)) != length(Genotypes_VCF_ID)) {
        warning("Duplicate Genotype IDs found!")
    }

    # Add cleaned IDs to genotype table
    Genotypes_Table_Mod <- cbind(Genotypes_VCF_ID, t(Genotypes_VCF))
    colnames(Genotypes_Table_Mod) <- c("Strain", markerIDMod)

     # !
    # Recoding genotypes from BB/AB/AA to -1/0/1 and convert to numeric
    Genotypes_Table_Mod_Num <- apply(Genotypes_Table_Mod[,-1], 2, function(x) {
        as.numeric(gsub("AA", "1", gsub("AB", "0", gsub("BB", "-1", x)))) + 1
    })
    
    Genotypes_Table_Mod_Num <- cbind.data.frame(Genotypes_Table_Mod[,1], Genotypes_Table_Mod_Num)
    colnames(Genotypes_Table_Mod_Num)[1] <- "Strain"
  
    # Remove duplicated IDs
    if (anyDuplicated(Genotypes_Table_Mod_Num[,1])) {
        dupIndices_Num <- which(duplicated(Genotypes_Table_Mod_Num[,1]))
		Genotype_Table_Num_noDup <- Genotypes_Table_Mod_Num[-dupIndices_Num,]
    } else{
        Genotype_Table_Num_noDup <- Genotypes_Table_Mod_Num
    }
    rownames(Genotype_Table_Num_noDup) <- Genotype_Table_Num_noDup[,1]

    # Filtered genotype table (numeric format)
    Genotypes_Table_Mod_Num_Filt <- apply(Genotype_Table_Num_noDup[,-1], 2, as.numeric)
    rownames(Genotypes_Table_Mod_Num_Filt) <- rownames(Genotype_Table_Num_noDup)

    # print(paste(dim(Genotypes_Table_Mod_Num_Filt)))
    # print(paste(rownames(Genotypes_Table_Mod_Num_Filt)[1:5])) 
	# print(paste(TestIDs))

    # Split train and test sets
    if (!is.null(TestIDs)) {
	    print("TestIn")
        #testIndices <- which(Genotypes_Table_Mod_Num_Filt[,1] %in% TestIDs)
        testIndices <- which(as.character(rownames(Genotypes_Table_Mod_Num_Filt)) %in% as.character(TestIDs))
		trainIndices <- setdiff(c(1:nrow(Genotypes_Table_Mod_Num_Filt)),testIndices)
        
        Test_Genotypes_Table_Mod_Num_Filt <- Genotypes_Table_Mod_Num_Filt[testIndices,]
        Train_Genotypes_Table_Mod_Num_Filt <- Genotypes_Table_Mod_Num_Filt[trainIndices,]
    }else{
        Train_Genotypes_Table_Mod_Num_Filt <- Genotypes_Table_Mod_Num_Filt
		Test_Genotypes_Table_Mod_Num_Filt <- NULL
    }

    # Process phenotype IDs
    Pheno[,1] <- gsub("[-_.()]", "", Pheno[,1])
    if (any(grepl("MG", Pheno[,1]))) {
        Pheno[,1] <- gsub("MG.*", "", Pheno[,1])
    }
	colnames(Pheno)[1] <- "StrainID"

    # Merge geno and pheno data
    trainStrainID <- gsub("[-_.()]", "", rownames(Train_Genotypes_Table_Mod_Num_Filt))
    commonStrainID <- Pheno[which(Pheno[, "StrainID"] %in% trainStrainID), "StrainID"]
    Diff_StrainID <- setdiff(trainStrainID, commonStrainID) # IDs not matching

    GenoTable_Filtered <- cbind.data.frame(rownames(Train_Genotypes_Table_Mod_Num_Filt),Train_Genotypes_Table_Mod_Num_Filt)
    colnames(GenoTable_Filtered)[1] <- "StrainID"
	
	PhenoTable_Filtered <- Pheno[which(Pheno[,"StrainID"] %in% trainStrainID),]
  
    Data_Table1_Num <- merge(GenoTable_Filtered, PhenoTable_Filtered,by="StrainID")

     print("Merge Done")
    # Check if merge was successful
    mergedIDs <- Data_Table1_Num$StrainID
    missingInGeno <- setdiff(PhenoTable_Filtered$StrainID, mergedIDs)
    missingInPheno <- setdiff(GenoTable_Filtered$StrainID, mergedIDs)

   
    # Handle merge failure gracefully
	if(length(mergedIDs) > 1){
     if (length(missingInGeno) > 0 | length(missingInPheno) > 0) {
	     mergeSuccess <- TRUE
		 if(length(missingInGeno)==0){missingInGeno="None"} 
		 if(length(missingInPheno)==0){missingInPheno="None"} 
		 warning("Some IDs did not match between genotype and phenotype tables.")
     }else if(length(missingInGeno)==0 && length(missingInPheno)== 0){
	      missingInGeno <- "None"
		  missingInPheno <- "None"
		  mergeSuccess <- TRUE
     }
    }else if(length(mergedIDs) <=1){
        mergeSuccess <- FALSE
    }
	
	mergeStatus <- list(
        missingInGeno = missingInGeno,  # IDs present in Pheno but not in Geno
        missingInPheno = missingInPheno # IDs present in Geno but not in Pheno
    )


    # Handle potential MG columns
    if ("MG" %in% colnames(Data_Table1_Num)) {
        positiveIndices <- which(as.numeric(as.character(Data_Table1_Num[,"MG"])) > 0)
        zeroIndices <- which(as.numeric(as.character(Data_Table1_Num[,"MG"])) == 0)
        negativeIndices <- which(as.numeric(as.character(Data_Table1_Num[,"MG"])) < 0)
      
        StrainIDModPhenoMG <- rep(0, nrow(Data_Table1_Num))
        StrainIDModPhenoMG[positiveIndices] <- paste(Data_Table1_Num[positiveIndices,"StrainID"], "MG", as.roman(Data_Table1_Num[positiveIndices,"MG"]), sep="")
        StrainIDModPhenoMG[zeroIndices] <- paste(Data_Table1_Num[zeroIndices,"StrainID"], "MG", Data_Table1_Num[zeroIndices,"MG"], sep="")
        StrainIDModPhenoMG[negativeIndices] <- paste(Data_Table1_Num[negativeIndices,"StrainID"], "MG00", sep="")
      
        Data_Table1_Num <- cbind.data.frame(Data_Table1_Num, StrainIDModPhenoMG)
        colnames(Data_Table1_Num)[ncol(Data_Table1_Num)] <- "StrainIDModPheno"
    } else {
        Data_Table1_Num <- cbind.data.frame(Data_Table1_Num,Data_Table1_Num[,"StrainID"])
        colnames(Data_Table1_Num)[ncol(Data_Table1_Num)] <- "StrainIDModPheno"
    }

    # Remove duplicate StrainID
    if (anyDuplicated(Data_Table1_Num[,"StrainIDModPheno"])) {
        dupIDIndices <- which(duplicated(Data_Table1_Num[,"StrainIDModPheno"]))
        Data_Table_Num <- Data_Table1_Num[-dupIDIndices,]
    } else {
        Data_Table_Num <- Data_Table1_Num
    }
  
    rownames(Data_Table_Num) <- Data_Table_Num[,"StrainIDModPheno"]
  
    return(list(
        Train_Data_Table_Num_Filt = Data_Table_Num, 
        Test_Genotypes_Table_Mod_Num_Filt = Test_Genotypes_Table_Mod_Num_Filt, 
        mergeStatus = mergeStatus,
        mergeSuccess = mergeSuccess
    ))

}


#' Title of getProcessedData
#'
#' Description of what getProcessedData does.
#'
#' @param Data_Table_Num_List Description of Data_Table_Num_List.
#' @param trait Description of trait.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getProcessedData
#' result <- getProcessedData(...)
#' }

getProcessedData <- function(Data_Table_Num_List,trait){

     TrainData_Table_Num <- Data_Table_Num_List[[1]]
	 TestData_Table_Num_Filt<- Data_Table_Num_List[[2]]
	 
	if(length(trait)>0){
	 ### Remove lines with 'NA' for trait values
	  if(length(trait)==1){
	   NAIndices <-  which(is.na(TrainData_Table_Num[,trait]))
	   if(length(NAIndices)>=1){
	    Train_Data_Table_Num_Filt <- TrainData_Table_Num[-NAIndices,]
	   }
	   if(length(NAIndices)<1){
	    Train_Data_Table_Num_Filt <- TrainData_Table_Num
	   }
	  }else if(length(trait)>1){
	   NAIndices <- c(unlist(apply(TrainData_Table_Num[,trait],2,function(x)which(is.na(x)))))
	   if(length(NAIndices)>=1){
	  
        Train_Data_Table_Num_Filt <- TrainData_Table_Num[-NAIndices,]
	   }
	   if(length(NAIndices)<1){
	  
        Train_Data_Table_Num_Filt <- TrainData_Table_Num
	   }
      }
	  
     return(list(Train_Data_Table_Num_Filt,TestData_Table_Num_Filt)) 
	 
	 }else if(length(trait)==0){
	 
	  return("Warning !.You've not selected any trait")
	  }
}



#' Title of getTasObj
#'
#' Description of what getTasObj does.
#'
#' @param infileVCF Description of infileVCF.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getTasObj
#' result <- getTasObj(...)
#' }
getTasObj <- function(infileVCF){
    tasGeno <- rTASSEL::readGenotypeTableFromPath(
		path = infileVCF,
		sortPositions=TRUE
	)
    return(tasGeno)
}



#' Title of getFilteredSitesGenoData
#'
#' Description of what getFilteredSitesGenoData does.
#'
#' @param tasGeno Description of tasGeno.
#' @param siteMinCnt Description of siteMinCnt.
#' @param MAF Description of MAF.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getFilteredSitesGenoData
#' result <- getFilteredSitesGenoData(...)
#' }

getFilteredSitesGenoData <- function(tasGeno,siteMinCnt,MAF){
   
	tasGenoFilt <- rTASSEL::filterGenotypeTableSites(
		tasObj = tasGeno,
		siteMinCount = siteMinCnt,
		siteMinAlleleFreq = MAF,
		siteMaxAlleleFreq = 1.0,
		siteRangeFilterType = "none"
	)

   
   return(tasGenoFilt)
}



#' Title of getFilteredTaxaGenoData
#'
#' Description of what getFilteredTaxaGenoData does.
#'
#' @param tasGeno Description of tasGeno.
#' @param MinNotMissing Description of MinNotMissing.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getFilteredTaxaGenoData
#' result <- getFilteredTaxaGenoData(...)
#' }
getFilteredTaxaGenoData <- function(tasGeno,MinNotMissing){

   	
    tasGenoFilt <- rTASSEL::filterGenotypeTableTaxa(
	   tasGeno,
	   minNotMissing = MinNotMissing,
	   minHeterozygous = 0,
	   maxHeterozygous = 1,
	   taxa = NULL
	)

   return(tasGenoFilt)
}





#' Title of getGenoTas_to_DF_Old
#'
#' Description of what getGenoTas_to_DF_Old does.
#'
#' @param tasGeno Description of tasGeno.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getGenoTas_to_DF_Old
#' result <- getGenoTas_to_DF_Old(...)
#' }

getGenoTas_to_DF_Old <- function(tasGeno){

    tasSumExp <- rTASSEL::getSumExpFromGenotypeTable(tasObj = tasGeno)
	
	tasGenoDF <- (SummarizedExperiment::assays(tasSumExp)[[1]])
	colnames(tasGenoDF) <- SummarizedExperiment::colData(tasSumExp)[,"Sample"]
	
   	
    ### Extract Table Report DF 
	
	tableReport <- rJava::new(
    rJava::J("net.maizegenetics.dna.map.PositionListTableReport"),
    tasGeno %>% rTASSEL:::getPositionList()) %>% 
    rTASSEL:::tableReportToDF() %>% as.data.frame()
      	
	varSplit <- strsplit(tableReport[,"VARIANT"],"/")
	varSplitTab <- cbind.data.frame(unlist(lapply(varSplit,function(x) x[1])),unlist(lapply(varSplit,function(x) x[2])))

    vcfIDTab <- cbind.data.frame(tableReport[,c("Name","Chromosome","Position")],varSplitTab)
	colnames(vcfIDTab) <- c("SNPID","Chrom","Position","REF","ALT")
		
	gt2d_tasGeno <-as_tibble(cbind.data.frame(vcfIDTab,tasGenoDF))
	
	return(gt2d_tasGeno)
  
}



#' Title of getGenoTas_to_DF
#'
#' Description of what getGenoTas_to_DF does.
#'
#' @param tasGeno Description of tasGeno.
#'
#' @importFrom rJava new J 
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getGenoTas_to_DF
#' result <- getGenoTas_to_DF(...)
#' }

getGenoTas_to_DF <- function(tasGeno){

    tasGenoMat <- as.matrix(tasGeno)
	
	tasGenoDF <- as.data.frame(t(tasGenoMat))
   	tasGenoDF$SNPID <- rownames(tasGenoDF)
	
    ### Extract Table Report DF 
	
	tableReport <- rJava::new(
		rJava::J("net.maizegenetics.dna.map.PositionListTableReport"),
		tasGeno %>% rTASSEL:::getPositionList()) %>% 
    rTASSEL:::tableReportToDF() %>% as.data.frame()
      	
	varSplit <- strsplit(tableReport[,"VARIANT"],"/")
	varSplitTab <- cbind.data.frame(unlist(lapply(varSplit,function(x) x[1])),unlist(lapply(varSplit,function(x) x[2])))

    vcfIDTab <- cbind.data.frame(tableReport[,c("Name","Chromosome","Position")],varSplitTab)
	colnames(vcfIDTab) <- c("SNPID","Chrom","Position","REF","ALT")
		
	# gt2d_tasGeno <-as_tibble(cbind.data.frame(vcfIDTab,tasGenoDF))
	gt2d_tasGeno <- as_tibble(merge(vcfIDTab,tasGenoDF,by="SNPID"))
	return(gt2d_tasGeno)
  
}



#' Title of getPredictionData
#'
#' Description of what getPredictionData does.
#'
#' @param Data_Table_Num_List Description of Data_Table_Num_List.
#' @param noCandidates Description of noCandidates.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getPredictionData
#' result <- getPredictionData(...)
#' }
getPredictionData <- function(Data_Table_Num_List,noCandidates){

     TrainData_Table_Num_Filt <- Data_Table_Num_List[[1]]
	 TestData_Table_Num_Filt <- Data_Table_Num_List[[2]]
	 set.seed(125)
	 CandidateIndices <- sample(c(1:nrow(TrainData_Table_Num_Filt)),noCandidates)
	 Candidates <- as.character(TrainData_Table_Num_Filt[CandidateIndices,1])
	 
	 if(length(Candidates)==nrow(TrainData_Table_Num_Filt)){
	  
      Candidate_Data_Table_Num_Filt <- TrainData_Table_Num_Filt
	 }
	 if(length(Candidates)< nrow(TrainData_Table_Num_Filt)){
	  
	  Candidate_Data_Table_Num_Filt <- TrainData_Table_Num_Filt[CandidateIndices,]
      
	 }
     
     return(list(Candidate_Data_Table_Num_Filt,TestData_Table_Num_Filt)) 
	 
} 


 
#' Title of getLibStats
#'
#' Description of what getLibStats does.
#'
#' @param buildLib Description of buildLib.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getLibStats
#' result <- getLibStats(...)
#' }

getLibStats <- function(buildLib){
    
    nInd <- nrow(buildLib)/2
    nMarkers <- ncol(buildLib[,-1]) 
    out <- paste("Uploaded haplotype library contains ",nInd," individuals and ",nMarkers," markers",sep="")
    return(out)
  }
  
  ##


  
 

	
#function to remove repeated genotypes

#' Title of cleanREPV2
#'
#' Description of what cleanREPV2 does.
#'
#' @param y Description of y.
#' @param gen Description of gen.
#' @param fam=NULL Description of fam=NULL.
#' @param thr=0.95 Description of thr=0.95.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of cleanREPV2
#' result <- cleanREPV2(...)
#' }

cleanREPV2 = function(y,gen,fam=NULL,thr=0.95){ 

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


###############################################################################


#' Title of getFixedData_List
#'
#' Description of what getFixedData_List does.
#'
#' @param Processed_Data_Table_Num_Split_Filt Description of Processed_Data_Table_Num_Split_Filt.
#' @param trait Description of trait.
#' @param fixedCoVar Description of fixedCoVar.
#' @param target_ID Description of target_ID.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getFixedData_List
#' result <- getFixedData_List(...)
#' }
getFixedData_List <- function(Processed_Data_Table_Num_Split_Filt,trait,fixedCoVar,target_ID){

  traitSel <- c(trait,fixedCoVar)
    
  Processed_Merged_BLUEs <- Processed_Data_Table_Num_Split_Filt[[1]][,c("StrainID",traitSel)]
  
  Processed_BLUEs_Test <- matrix(NA,nrow=nrow(Processed_Data_Table_Num_Split_Filt[[2]]), ncol=ncol(Processed_Merged_BLUEs))
  Processed_BLUEs_Test[,1]<- Processed_Data_Table_Num_Split_Filt[[2]][,1]
  Processed_BLUEs_Test[,ncol(Processed_BLUEs_Test)] <- rep("BD",nrow(Processed_BLUEs_Test))
  
  colnames(Processed_BLUEs_Test) <- colnames(Processed_Merged_BLUEs)
  Processed_BLUEs_Comb <- rbind(Processed_Merged_BLUEs,Processed_BLUEs_Test)
  
  
  testIndices <- which(Processed_BLUEs_Comb[,"StrainID"] %in% target_ID)
  trainIndices <- setdiff(c(1:nrow(Processed_BLUEs_Comb)),testIndices)
  
  X <- model.matrix(~as.factor(Processed_BLUEs_Comb[,fixedCoVar]))
  Fixed.X <- X[trainIndices,]
  Test.X <- X[testIndices,]

  return(list(Fixed.X,Test.X))
}



#' Title of getPhenoMEData
#'
#' Description of what getPhenoMEData does.
#'
#' @param PhenoME Description of PhenoME.
#' @param TraitME Description of TraitME.
#' @param nSelTraits Description of nSelTraits.
#' @param IDColsME Description of IDColsME.
#' @param StrainME Description of StrainME.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getPhenoMEData
#' result <- getPhenoMEData(...)
#' }
getPhenoMEData <- function(PhenoME,TraitME,nSelTraits,IDColsME,StrainME){

   print("Get Pheno ME Data")

	Data_Trt_List <- list()
	Data_Trt_Filt_List <- list()

	traits <- TraitME
	nTraits <- nSelTraits
	
	IDCols <- IDColsME

	for(nTrt in 1:nTraits){
	  
	  trait <- traits[nTrt]
	  selCols <- trait
	  Data_Trt <- PhenoME[,c(IDCols,selCols)]
	 
	  print(colnames(Data_Trt))
	  print(StrainME)
	 # colnames(Data_Trt)[which(colnames(Data_Trt) %in% StrainME)] <- "Strain"
	  Data_Trt$Strain <- as.factor(Data_Trt$Strain)
	  Data_Trt$Loc <- as.factor(Data_Trt$Loc)
	  
	  Loc <- levels(factor(Data_Trt$Loc))
	 
	  Data_Trt_Filt <-  Data_Trt
	  Data_Trt_Filt <- droplevels(Data_Trt_Filt)
	  
	  ####
	  Data_Trt_List[[nTrt]] <- Data_Trt
	  Data_Trt_Filt_List[[nTrt]] <- Data_Trt_Filt
	  
	}

## Pheno Mat List for traits
	Ph_Mat_list <- list()
	for(nTrt in 1:nTraits){
	  trait <- traits[nTrt]
	  Ph_Mat_0 <- reshape2::dcast(Strain ~ Loc, value.var = trait,fun.aggregate =mean,data=Data_Trt_List[[nTrt]])
	 ##Rm NA
	   naInd <- which(is.na(Ph_Mat_0[,1]) | is.nan (Ph_Mat_0[,1]))
	  if(length(naInd)>0){Ph_Mat_1 <- Ph_Mat_0[-naInd,] }else{Ph_Mat_1 <- Ph_Mat_0[,]}
	 ## Rm Dup  
	  dupInd <- which(duplicated(Ph_Mat_1[,1]))
	  if(length(dupInd)>0){Ph_Mat <- Ph_Mat_1[-dupInd,-1] }else{Ph_Mat <- Ph_Mat_1[,-1]}
	  rownames(Ph_Mat) <-  Ph_Mat_1[,1]
	  Ph_Mat_list[[nTrt]] <- Ph_Mat
	}

  return(list(Ph_Mat_list,Data_Trt_List,Data_Trt_Filt_List))
}




####################################################### 

#' Title of getMergedDataME
#'
#' Description of what getMergedDataME does.
#'
#' @param phData Description of phData.
#' @param genoImp Description of genoImp.
#' @param TargetIDs Description of TargetIDs.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getMergedDataME
#' result <- getMergedDataME(...)
#' }
getMergedDataME <- function(phData,genoImp,TargetIDs){

 print("MergingData")

	 DT_1_List <- list()
	 DT_1_Filt_List <- list()

     Ph_Mat_list <- phData[[1]]
	 Data_Trt_List <- phData[[2]]
	 Data_Trt_Filt_List <- phData[[3]]

	 genoDat_List <- list() 
	 phenoDat_List <- list()
	 
	 nTraits <- length(Ph_Mat_list)

	for(nTrt in 1:nTraits){
	  
	  phData <- Ph_Mat_list[[nTrt]]
	  
	  rmID <- which(colnames(genoImp) %in% c("SNPID","Chrom","Position","REF","ALT"))
	  genoImpDF <- as.data.frame(t(genoImp[,-c(1:5)]))
	  colnames(genoImpDF) <- paste("ss",as.vector(as.data.frame(genoImp)[,"SNPID"]),sep="-")
	  
	  ### Merge Data
	  
	  Data <- merge(phData,genoImpDF,by=0)
	  dim(Data)
	 
	 
	  #### GenoData
	  genoCols <- grep("ss",colnames(Data))
	  genoDat <- as.matrix(apply(Data[,genoCols],2,as.numeric))
	  rownames(genoDat) <- Data[,1]
	  is.matrix(genoDat)
	  
	  anyNA(genoDat)
	  table(genoDat)
	  genoDat_List[[nTrt]] <- genoDat
	  
	  #### PhenoData
	  phCols <- c(2:(genoCols[1]-1))
	  phenoDat <- as.matrix(apply(Data[,phCols],2,as.numeric))
	  rownames(phenoDat) <- Data[,1]
	  is.matrix(phenoDat)
	  
	  phenoDat_List[[nTrt]] <- phenoDat
	  
	  ###
	  
	  StrainGeno <- rownames(genoDat)
	  Trt_Data <- Data_Trt_List[[nTrt]]
	  
	  ###
	  selInd <- which(Trt_Data$Strain %in% StrainGeno)
	  DT_1 <- Trt_Data[selInd,]
	  DT_1_List[[nTrt]] <- DT_1
	  
	  #############
	  
	  Trt_Data_Filt <- Data_Trt_Filt_List[[nTrt]]
	  
	  selInd2 <- which(Trt_Data_Filt$Strain %in% StrainGeno)
	  DT_1_Filt <- Trt_Data_Filt[selInd2,]
	  DT_1_Filt <- droplevels(DT_1_Filt)
	  
	  DT_1_Filt_List[[nTrt]] <- DT_1_Filt
	  
	
	##
	
	 if(nTrt ==1){
		print("Merge Done")
		# Check if merge was successful
		mergedIDs <- Data[,1]
		oriPhIDs <- rownames(phData)
		oriGenoIDs <- rownames(genoImpDF)
		
		missingInGeno <- setdiff(oriPhIDs,mergedIDs)
		missingInPheno <- setdiff(oriGenoIDs,mergedIDs)

		mergeStatus <- list(
			missingInGeno = missingInGeno,  # IDs present in Pheno but not in Geno
			missingInPheno = missingInPheno # IDs present in Geno but not in Pheno
		)

		# Handle merge failure gracefully
		if(length(mergedIDs) >1){
		 if(length(missingInGeno) > 0 || length(missingInPheno) > 0){
		  mergeSuccess <- TRUE
		  warning("Some IDs did not match between genotype and phenotype tables.")
		 }else{ 
		  mergeSuccess <- TRUE
		 }
		}else if(length(mergedIDs) <=1){
		 mergeSuccess <- FALSE
		 warning("IDs did not match between genotype and phenotype tables.")
		}
	 }
	
	}
	
	return(list(
        DT_1_Filt_List = DT_1_Filt_List, 
        genoDat_List = genoDat_List, 
        mergeStatus = mergeStatus,
        mergeSuccess = mergeSuccess
    ))
	
 
}





#' Title of getPhenoDistbnPlots
#'
#' Description of what getPhenoDistbnPlots does.
#'
#' @param DT_1_Filt_List Description of DT_1_Filt_List.
#' @param TraitME Description of TraitME.
#' @param nTrt Description of nTrt.
#' @param outDirPath Description of outDirPath.
#'
#' @import ggplot2
#' @importFrom ggplot2 aes geom_boxplot position_dodge labs theme_minimal theme element_text element_blank
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getPhenoDistbnPlots
#' result <- getPhenoDistbnPlots(...)
#' }

getPhenoDistbnPlots <- function(DT_1_Filt_List,TraitME,nTrt,outDirPath){

###### Distribution of Trait BLUEs & BLUPs across locations 
     
    DT_1_Filt <- DT_1_Filt_List[[nTrt]]
	trait <- TraitME
#### Trait BLUPs
 
		 
	DT_2B <- DT_1_Filt
	DT_2B$Loc <- as.factor(DT_2B$Loc)
	  
	table(DT_2B$Loc)
	
	nanInd <-   which(is.na(DT_2B[,trait]) | is.nan(DT_2B[,trait]))
	if(length(nanInd)>0){DT_2B <- DT_2B[-nanInd,]}
	
	minY <- min(DT_2B[,trait])-4
	maxY <- max(DT_2B[,trait])+4
	  
	plot3 <- ggplot2::ggplot(DT_2B, aes(x = Loc, y = base::get(trait),fill=Loc)) +
		geom_boxplot(position = position_dodge(width = 0.8), color = "black") +
		labs(x = "Location", y = gsub("_"," ",trait))+
		#scale_fill_manual(values = c("North" = "gray","Central"="blue", "South" = "red")) +
		theme_minimal() +
		theme(axis.text.x = element_text(size=26,angle = 45, hjust = 1,face="bold",color = "black"),
			  axis.text.y= element_text(size=26,face="bold",color="black"),
			  panel.grid = element_blank(),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  axis.title.y = element_text(size = 26, face = "bold"),
			  axis.title.x = element_text(size = 26, face = "bold"),
			  legend.title = element_text(size = 26, face = "bold"),
			  legend.text = element_text(size = 26, face = "bold"))

     if(is.null(outDirPath) || outDirPath==""){
	   oFN <- paste(trait," Distribution.png",sep="")
	 }else{
       oFN <- paste(outDirPath,"/",trait," Distribution.png",sep="")
	 }
	 tryCatch({
             withCallingHandlers({
               png(oFN,width=1024,height=768,pointsize=20)
             }, warning = function(w){
               # This function will execute if a warning is issued
               # Capture the warning message
               warningMessage <- conditionMessage(w)
               print(paste("Warning captured:", warningMessage))
               
               # Invoke the default warning handler (optional)
               invokeRestart("muffleWarning")
             }, error=function(e) {
               print(paste("Error encountered:", e$message))
             })
          })
    
	 
	 print(plot3)
	 dev.off()

    return(plot3) 	 
	 
  }



 #' Title of writeCVOutTable
#'
#' Description of what writeCVOutTable does.
#'
#' @param CVROut Description of CVROut.
#' @param type Description of type.
#' @param outDirPath Description of outDirPath.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of writeCVOutTable
#' result <- writeCVOutTable(...)
#' }
writeCVOutTable <- function(CVROut,type,outDirPath){
    
    if(type=="ST"){
      oFN <-  paste(outDirPath,"/SingleTrait_CrossValidation_Ouput.txt",sep="")
    }else if(type=="MT"){
      oFN <-  paste(outDirPath,"/MultiTrait_CrossValidation_Ouput.txt",sep="") 
    }else if (type=="ME"){
      oFN <-  paste(outDirPath,"/MultiEnv_CrossValidation_Ouput.txt",sep="") 
    }
    
    write.table(CVROut,oFN,quote=F,sep="\t")
 }


 
 #' Title of writeGPOutTable
#'
#' Description of what writeGPOutTable does.
#'
#' @param GPOut Description of GPOut.
#' @param type Description of type.
#' @param outDirPath Description of outDirPath.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of writeGPOutTable
#' result <- writeGPOutTable(...)
#' }
writeGPOutTable <- function(GPOut,type,outDirPath){
    
    if(type=="ST"){
      oFNTr <-  paste(outDirPath,"/SingleTrait_GP_Ouput_TrainingData.txt",sep="")
	  oFNTgt <-  paste(outDirPath,"/SingleTrait_GP_Ouput_TargetData.txt",sep="")
    }else if(type=="MT"){
      oFNTr <-  paste(outDirPath,"/MultiTrait_GP_Ouput_TrainingData.txt",sep="") 
	  oFNTgt <- paste(outDirPath,"/MultiTrait_GP_Ouput_TargetData.txt",sep="") 
    }
	
    
    write.table(GPOut[[1]],oFNTr,quote=F,sep="\t",row.names=FALSE)
	
	if(!is.null(GPOut[[2]])){
	 write.table(GPOut[[2]],oFNTgt,quote=F,sep="\t",row.names=FALSE)
    }
 }
 
###

 
##### 
 
#' Title of getCombinedTab
#'
#' Description of what getCombinedTab does.
#'
#' @param outputListME Description of outputListME.
#' @param TraitME Description of TraitME.
#' @param IDColsME Description of IDColsME.
#' @param IDColME Description of IDColME.
#' @param fitEnvCov Description of fitEnvCov.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getCombinedTab
#' result <- getCombinedTab(...)
#' }
getCombinedTab <- function(outputListME,TraitME,IDColsME,IDColME,fitEnvCov){

 if(fitEnvCov==FALSE){
  Models <- c("MM","MDs","MDe")
 }else if(fitEnvCov == TRUE){
  Models <- c("EMM","EMDs","EMDe")
 }
  nTraits <- length(TraitME)
  
  for(nTrt in 1:nTraits){
  
    Pred_Out <- outputListME[[nTrt]]
	traitME <- TraitME[nTrt]
    Fit_Pred <- Pred_Out[[3]]
    nMod <- length(Fit_Pred)
	
    for(nM in 1:nMod){
	 if(nM==1){
	    Fit_Out <- Fit_Pred[[nM]]
		uniqIDCol <- which(colnames(Fit_Out) == IDColME)
		colnames(Fit_Out)[uniqIDCol] <- "UniqID"
		colnames(Fit_Out) <- gsub("Obs",paste("Obs",Models[nM],sep="-"),colnames(Fit_Out))
		colnames(Fit_Out) <- gsub("Pred",paste("Pred",Models[nM],sep="-"),colnames(Fit_Out))
	 }else if(nM>1){ 
	 
	    Fit_Pred_Tab <- Fit_Pred[[nM]]
		uniqIDCol <- which(colnames(Fit_Pred_Tab) == IDColME)
		colnames(Fit_Pred_Tab)[uniqIDCol] <- "UniqID"
		colnames(Fit_Pred_Tab) <- gsub("Obs",paste("Obs",Models[nM],sep="-"),colnames(Fit_Pred_Tab))
		colnames(Fit_Pred_Tab) <- gsub("Pred",paste("Pred",Models[nM],sep="-"),colnames(Fit_Pred_Tab))
		
	    Fit_Out <- merge(Fit_Out,Fit_Pred_Tab,by="UniqID",all=TRUE)
	 }
	}
	
	if(nTrt ==1){
	 colnames(Fit_Out) <- gsub("Obs",paste(traitME,"Obs",sep=""),colnames(Fit_Out))
	 colnames(Fit_Out) <- gsub("Pred",paste(traitME,"Pred",sep=""),colnames(Fit_Out))
	 Fit_Out_Tab <- Fit_Out
	 
	}else{
	 colnames(Fit_Out) <- gsub("Obs",paste(traitME,"Obs",sep=""),colnames(Fit_Out))
	 colnames(Fit_Out) <- gsub("Pred",paste(traitME,"Pred",sep=""),colnames(Fit_Out))
	 Fit_Out_Tab <- merge(Fit_Out_Tab,Fit_Out,by="UniqID",all=TRUE) 
	}

    	
   }
   
    valCols <- colnames(Fit_Out_Tab)[grep("Obs|Pred",colnames(Fit_Out_Tab))]
	selOutCols <- c("UniqID",setdiff(IDColsME,IDColME),valCols)
	Fit_Out_Tab_Sel <- Fit_Out_Tab[,selOutCols]
   
	
	return(Fit_Out_Tab_Sel)
	
  }



#' Title of getOutTab_ME_CV
#'
#' Description of what getOutTab_ME_CV does.
#'
#' @param ME_Out_CV_Trt Description of ME_Out_CV_Trt.
#' @param CVMet Description of CVMet.
#' @param Traits Description of Traits.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getOutTab_ME_CV
#' result <- getOutTab_ME_CV(...)
#' }
getOutTab_ME_CV <- function(ME_Out_CV_Trt,CVMet,Traits){
  
   nTraits <- length(ME_Out_CV_Trt) 
   
   if(CVMet=="CV_LOFO"){
	   Cor_LOF_CV_Tab_Out_Trt <- list()
	   
	   
	   for(nTrt in 1:nTraits){
		 ME_LOFCV <- ME_Out_CV_Trt[[nTrt]]
		 CorMM_LOF_CV <- unlist(lapply(ME_LOFCV,function(x) x[[1]][[1]]))
		 CorMDs_LOF_CV <- unlist(lapply(ME_LOFCV,function(x) x[[1]][[2]]))
		 CorMDe_LOF_CV <- unlist(lapply(ME_LOFCV,function(x) x[[1]][[3]]))
		 cor_LOF_CV_Tab <- cbind(CorMM_LOF_CV,CorMDs_LOF_CV,CorMDe_LOF_CV)
		 cor_LOF_CV_Tab_Out <- apply(cor_LOF_CV_Tab,2,summary)
		 colnames(cor_LOF_CV_Tab_Out) <- c("MM","MDs","MDe")
		 Cor_LOF_CV_Tab_Out_Trt[[nTrt]] <- cor_LOF_CV_Tab_Out
		}
	   
	   Cor_LOF_CV_Out_Trt_Tab <- do.call(cbind,lapply(c(1:nTraits),function(x) { 
				tab <- Cor_LOF_CV_Tab_Out_Trt[[x]]
				colnames(tab) <- paste(Traits[x],colnames(tab),sep="-")
				tab
		  })
		)
	
   	   Cor_CV_Out_Trt_Tab <- Cor_LOF_CV_Out_Trt_Tab

    }else if(CVMet!="CV_LOFO"){	
	
		for(nTrt in 1:nTraits){
		  ME_Out_CV <- ME_Out_CV_Trt[[nTrt]]
		  nReps <- length(ME_Out_CV[[3]])
		
		  for(nrep in 1:nReps){
			
			  fit_Tab <- ME_Out_CV[[3]][[nrep]]
			  nModels <- length(ME_Out_CV[[3]][[nrep]][[1]])
			  
			  corOut <- rep(0,nModels)
			  for(nModel in 1:nModels){ 	  
			   Out_Tab_Model <- do.call(rbind,lapply(fit_Tab,function(x) x[[nModel]]))
			   corOut[nModel] <- cor(Out_Tab_Model[,"Obs"],Out_Tab_Model[,"Pred"],use="pairwise.complete.obs")
			  }
			   names(corOut) <- c("MM","MDs","MDe")
			  
			  if(nrep==1){ 
				corOutReps <- corOut
			  }else if(nrep >1){ 
				corOutReps <- rbind(corOutReps,corOut)
			  }
		  } 
		
		if(nReps >1){
		  corOutTab <- apply(corOutReps,2,function(x) round(summary(x),digits=3))
		  colnames(corOutTab) <- paste(Traits[nTrt],c("MM","MDs","MDe"),sep="-")
		 }else if(nReps <=1){
		  corOutTab <- round(corOutReps,digits=3)
		  names(corOutTab) <- paste(Traits[nTrt],c("MM","MDs","MDe"),sep="-")
		 }
	
	    if(nTrt==1){
      	 Cor_CV_Out_Trt_Tab <- corOutTab
        }else if(nTrt >1){
		 
		 Cor_CV_Out_Trt_Tab <- rbind(Cor_CV_Out_Trt_Tab,corOutTab)
	    }
	  }
    }
  return(Cor_CV_Out_Trt_Tab)

}






#' Title of getFreq_Alleles
#'
#' Description of what getFreq_Alleles does.
#'
#' @param GenoTable Description of GenoTable.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getFreq_Alleles
#' result <- getFreq_Alleles(...)
#' }
getFreq_Alleles <- function(GenoTable){


	nMarkers <- dim(GenoTable)[1]
	nIndividuals <- dim(GenoTable)[2]
	n_alleles <- 2*nIndividuals

	FreqList <- lapply(1:ncol(GenoTable),function(j){ 
		p_2_P <- length(which(GenoTable[,j]==2))
		p_1_H <- length(which(GenoTable[,j]==1))
		p_0_Q <- length(which(GenoTable[,j]==0))

		count_1 <- ((2*p_2_P)+p_1_H )
		count_0 <- ((2*p_0_Q)+p_1_H )

		Freq_1 <- count_1/n_alleles
		Freq_0 <- count_0/n_alleles
		
		return(list(Freq_1,Freq_0))
	})
     
	 Freq_1_Vec <- unlist(lapply(FreqList,function(x) x[[1]]))
     Freq_0_Vec <- unlist(lapply(FreqList,function(x) x[[2]]))
	 Freq_Tab <- rbind.data.frame(Freq_1_Vec,Freq_0_Vec)
	
	return(Freq_Tab)

}



#' Title of numMAF
#'
#' Description of what numMAF does.
#'
#' @param Geno Description of Geno.
#' @param MAFTh Description of MAFTh.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of numMAF
#' result <- numMAF(...)
#' }
numMAF <- function(Geno,MAFTh){ 
   
   AF_Tab <- getFreq_Alleles(Geno)
   MAFVec <-  apply(AF_Tab,2,function(x) min(as.numeric(x)))
   numMAF_LT <- length(which(MAFVec <= MAFTh))
   return(numMAF_LT)
}

#' Title of numMissSites
#'
#' Description of what numMissSites does.
#'
#' @param GenoT Description of GenoT.
#' @param missSitesTH Description of missSitesTH.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of numMissSites
#' result <- numMissSites(...)
#' }
numMissSites <- function(GenoT,missSitesTH){ 

     Geno <- GenoT[,-c(1:5)]
     nInd <- ncol(Geno)
	 nSites <- nrow(Geno)
	 
	 NALen <- apply(Geno,1,function(x) length(which(is.na(x))))
	 numMissingSites <- length(which(NALen >= (missSitesTH*nInd)))
	 
	 return(numMissingSites)
 
 }
 
 
 
#' Title of numMissInd
#'
#' Description of what numMissInd does.
#'
#' @param Geno Description of Geno.
#' @param missIndTH Description of missIndTH.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of numMissInd
#' result <- numMissInd(...)
#' }
numMissInd <- function(Geno,missIndTH){ 

  nInd <- ncol(Geno)
  nSites <- nrow(Geno)
 
  NALen <- as.vector(apply(Geno,2,function(x) length(which(is.na(x)))))
  numMissingInd <- length(which(NALen >= missIndTH*nSites))
 
 return(numMissingInd)
 
} 



#' Title of getGenoDiff
#'
#' Description of what getGenoDiff does.
#'
#' @param Geno1 Description of Geno1.
#' @param Geno2 Description of Geno2.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getGenoDiff
#' result <- getGenoDiff(...)
#' }
getGenoDiff <- function(Geno1,Geno2){ 

  DiffInd <- abs(ncol(Geno1)-ncol(Geno2)) 
  DiffSites <- abs(nrow(Geno1)-nrow(Geno2))
  
  return(list(DiffSites,DiffInd))
  
 }


#' Title of getGenoQCStats
#'
#' Description of what getGenoQCStats does.
#'
#' @param GenoT Description of GenoT.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getGenoQCStats
#' result <- getGenoQCStats(...)
#' }
getGenoQCStats <- function(GenoT){ 
 
  Geno <- GenoT[,-c(1:5)]
  nInd <- ncol(Geno)
  nSites <- nrow(Geno)
  
  missNum <- sum(as.vector(apply(Geno,2,function(x) length(which(is.na(x))))))
  missFrac <- round(missNum/(nInd*nSites),digits=3)
  #missFrac <- (length(which(is.na(as.vector(Geno))))) / (length(as.vector(Geno)))
  
  code <- paste(names(table(apply(Geno,2,as.numeric))),sep="-",collapse="")
  
  genoLine <- paste("The genotype table has genotype scores for ", nInd," genotypes and ", nSites, " markers. \n",sep="")
  genoCodingLine <- paste("The genotype scores are coded in ",code," format. \n",sep="")
  missScoresLines <- paste(missFrac*100," % of the genotype scores in the table have missing values. \n",sep="") 
  outMsg <- paste(genoLine,genoCodingLine,missScoresLines,sep="")
  
  return(outMsg)
}
  
 
    
#' Title of getGenoQCStatsFilt1
#'
#' Description of what getGenoQCStatsFilt1 does.
#'
#' @param GenoT Description of GenoT.
#' @param GenoFilt1T Description of GenoFilt1T.
#' @param missSitesTH Description of missSitesTH.
#' @param MAFTH Description of MAFTH.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getGenoQCStatsFilt1
#' result <- getGenoQCStatsFilt1(...)
#' }
getGenoQCStatsFilt1 <- function(GenoT,GenoFilt1T,missSitesTH,MAFTH){ 

  Geno <- GenoT[,-c(1:5)]
  GenoFilt1 <- GenoFilt1T[,-c(1:5)]



  # missSitesGTH <- numMissSites(Geno,missSitesTH)
  # mafLTH <- numMAF(Geno,MAFTH)
  
  nInd <- ncol(GenoFilt1)
  nSites <- nrow(GenoFilt1)
  code <- paste(names(table(apply(GenoFilt1,2,as.numeric))),sep="-",collapse="")
  
  diffStats <- getGenoDiff(Geno,GenoFilt1)
  
  genoLine <- paste("The genotype table has genotype scores for ", nInd," genotypes and ", nSites, " markers. \n",sep="")
  genoCodingLine <- paste("The genotype scores are coded in ",code, " format. \n",sep="")
  # MAFLine <- paste(mafLTH," markers have MAF less than ",MAFTH, " threshold. \n",sep="\t")
  # missSiteLine <- paste(missSitesGTH," markers have missing values in more than ",missSitesTH*100," % of the genotypes. \n",sep="")
  filtLine <- paste(diffStats[[2]]," genotypes and ",diffStats[[1]]," markers have been removed in the filtered table. \n")
  #MAFLine,missSiteLine,
  outMsg <- paste(genoLine,genoCodingLine,filtLine,sep="")
  
  return(outMsg)
}
 


#' Title of getGenoQCStatsFilt2
#'
#' Description of what getGenoQCStatsFilt2 does.
#'
#' @param GenoFilt1T Description of GenoFilt1T.
#' @param GenoFilt2T Description of GenoFilt2T.
#' @param missIndTH Description of missIndTH.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getGenoQCStatsFilt2
#' result <- getGenoQCStatsFilt2(...)
#' }
getGenoQCStatsFilt2 <- function(GenoFilt1T,GenoFilt2T,missIndTH){ 

  GenoFilt1 <- GenoFilt1T[,-c(1:5)]
  GenoFilt2 <- GenoFilt2T[,-c(1:5)]
  
  nInd <- ncol(GenoFilt2)
  nSites <- nrow(GenoFilt2)
  diffStats <- getGenoDiff(GenoFilt1,GenoFilt2)
  code <- paste(names(table(apply(GenoFilt2,2,as.numeric))),sep="-",collapse="")
  
 # missIndGTH <- numMissInd(GenoFilt1,missIndTH)
  
  
  genoLine <- paste("The genotype table has genotype scores for ", nInd," genotypes and ", nSites, " markers. \n",sep="")
  genoCodingLine <- paste("The genotype scores are coded in ",code, " format. \n",sep="")
 # missIndLine <- paste(missIndGTH," genotypes had missing values in more than ",missIndTH*100," % of the markers. \n",sep="")
  filtLine <- paste(diffStats[[2]]," genotypes and ",diffStats[[1]]," markers have been removed in the filtered table. \n")
  #missIndLine,
  outMsg <- paste(genoLine,genoCodingLine,filtLine,sep="")
  
  return(outMsg)
}


####

#' Title of getGenoImp1Stats
#'
#' Description of what getGenoImp1Stats does.
#'
#' @param Geno_DF Description of Geno_DF.
#' @param GenoImp_DF1 Description of GenoImp_DF1.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getGenoImp1Stats
#' result <- getGenoImp1Stats(...)
#' }
getGenoImp1Stats <- function(Geno_DF,GenoImp_DF1){ 
   
  Geno <- Geno_DF[,-c(1:5)]
  nInd <- ncol(Geno)
  nSites <- nrow(Geno)
     
  GenoImp <- GenoImp_DF1[,-c(1:5)]
  nInd_I <- ncol(GenoImp)
  nSites_I <- nrow(GenoImp)
  
  missNum <- sum(as.vector(apply(Geno,2,function(x) length(which(is.na(x))))))
  missFrac <- round(missNum/(nInd*nSites),digits=3)
  code <- paste(names(table(apply(Geno,2,as.numeric))),sep="-",collapse="")

  
  
  missNumI <- sum(as.vector(apply(GenoImp,2,function(x) length(which(is.na(x))))))
  missFracI <- round(missNumI/(nInd_I*nSites_I),digits=3)
  codeI <- paste(names(table(apply(GenoImp,2,as.numeric))),sep="-",collapse="")

  diffStats <- getGenoDiff(Geno,GenoImp)
  
  breakLine0 <- "Input Genotype Table \n"
  genoLine <- paste("The input genotype table has genotype scores for ", nInd," genotypes and ", nSites, " markers. \n",sep="")
  genoCodingLine <- paste("The genotype scores are coded in ",code," format. \n",sep="")
  missScoresLines <- paste(missFrac*100," % of the genotype scores in the table have missing values. \n",sep="") 
  breakLine1 <- "\n"
  breakLine2 <- "Imputed Genotype Table \n"
  genoLineI <- paste("The imputed genotype table has genotype scores for ", nInd_I," genotypes and ", nSites_I, " markers. \n",sep="")
  genoCodingLineI <- paste("The genotype scores in the imputed table are coded in ",codeI," format. \n",sep="")
  missScoresLinesI <- paste(missFracI*100," % of the genotype scores in the table have missing values. \n",sep="") 
  
  filtLine <- paste(diffStats[[2]]," genotypes and ",diffStats[[1]]," markers have been removed in the imputed compared to the input table. \n")

  
  outMsg <- paste(breakLine0,genoLine,genoCodingLine,missScoresLines,breakLine1,breakLine2,genoLineI,genoCodingLineI,missScoresLinesI,filtLine,sep="")
 
}
   
