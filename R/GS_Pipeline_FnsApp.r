  
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
#' @export
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

  # gt2[(gt == "1/1")|(gt == "1|1")] <- 'BB'
  # gt2[(gt == "0/0")|(gt == "0|0")] <- 'AA'
  # gt2[(gt == "0/1")|(gt == "1/0")|(gt == "0|1")|(gt == "1|0")] <- 'AB'
  # gt2[(gt == "\\./\\.")|(gt == "\\.|\\.")] <- NA

  gt2 <- as_tibble(gt2d) %>%
    mutate(SNPID = fix_T$ID)

  gt.simple <- fix_T %>% dplyr::select("ID", "CHROM", POS, REF, ALT) %>%
    dplyr::rename(SNPID=ID, Chrom=CHROM, Position=POS) %>%
    dplyr::left_join(gt2, by = 'SNPID')

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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
	  genoImpDF <- as.data.frame(t(genoImp[,-rmID]))
	  colnames(genoImpDF) <- paste("ss",as.vector(as.data.frame(genoImp)[,"SNPID"]),sep="-")
	  
	  ### Merge Data
	  
	  Data <- merge(phData,genoImpDF,by=0)
	  dim(Data)
	 
	 
	  #### GenoData
	  genoCols <- grep("^ss",colnames(Data))
	  genoDat <- as.matrix(apply(Data[,genoCols],2,as.numeric))
	  rownames(genoDat) <- Data[,1]
	  is.matrix(genoDat)

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

getPhenoDistbnPlots <- function(DT_1_Filt_List,TraitME,nTrt,EnvVar,outDirPath){

###### Distribution of Trait BLUEs & BLUPs across locations

    DT_1_Filt <- DT_1_Filt_List[[nTrt]]
	trait <- TraitME
#### Trait BLUPs

	 DT_2B <- DT_1_Filt
	 envColInd <- which(colnames(DT_2B) %in% "Env")
	 colnames(DT_2B)[envColInd] <- "EnvVar"

	#DT_2B <- DT_1_Filt %>% rename(EnvVar = Env)
	DT_2B$EnvVar <- as.factor(DT_2B$EnvVar)

	print(table(DT_2B$EnvVar))

	nanInd <-   which(is.na(DT_2B[,trait]) | is.nan(DT_2B[,trait]))
	if(length(nanInd)>0){DT_2B <- DT_2B[-nanInd,]}

	minY <- min(DT_2B[,trait])-4
	maxY <- max(DT_2B[,trait])+4

	plot3 <- ggplot2::ggplot(DT_2B, aes(x = EnvVar, y = get(trait),fill=EnvVar)) +
		geom_boxplot(position = position_dodge(width = 0.8), color = "black") +
		labs(x = "Environment", y = gsub("_"," ",trait))+
		#scale_fill_manual(values = c("North" = "gray","Central"="blue", "South" = "red")) +
		theme_minimal() +
		theme(axis.text.x = element_text(size=20,angle = 45, hjust = 1,face="bold",color = "black"),
			  axis.text.y= element_text(size=20,face="bold",color="black"),
			  panel.grid = element_blank(),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  axis.title.y = element_text(size = 22, face = "bold"),
			  axis.title.x = element_text(size = 22, face = "bold"),
			  legend.position = "none",
			  legend.title = element_blank(), #text(size = 26, face = "bold"),
			  legend.text = element_blank()) #(size = 26, face = "bold"))

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
#' @export
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
#' @export
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
getCombinedTab <- function(outputListME, TraitME, IDColsME, IDColME,envVar,fitEnvCov,reaction) {

  # 1) Model labels
  if(fitEnvCov){ Models <-c("EMM","EMDs","EMDe")} else { Models <-c("MM","MDs","MDe")}
  if(reaction){ Models <-  c("RN-MM","RN-MDs","RN-MDe")}

  Pred_Out <- outputListME$preds
  nTraits  <- length(TraitME)

  # Helper: ensure a column exists, else stop with message
  .require_col <- function(df, col, who = "table") {
    if (!col %in% names(df)) stop(sprintf("Column '%s' not found in %s.", col, who), call. = FALSE)
  }

  # Helper: rename IDColME -> "UniqID" only if needed
  .rename_to_uniqid <- function(df) {
    if ("UniqID" %in% names(df)) return(df)
    .require_col(df, IDColME, "prediction table")
    names(df)[names(df) == IDColME] <- "UniqID"
    df
  }

  # Helper: suffix Obs/Pred with model label
  .tag_obs_pred <- function(df, label) {
    nn <- names(df)
    nn <- sub("^Obs$",  paste0("Obs-",  label), nn)
    nn <- sub("^Pred$", paste0("Pred-", label), nn)
    names(df) <- nn
    df
  }

  # Helper: left merge by UniqID and coalesce duplicate ID columns
  .safe_merge_by_uid <- function(x, y) {
    x <- .rename_to_uniqid(x)
    y <- .rename_to_uniqid(y)

    keep_idcols <- setdiff(IDColsME, IDColME) # other ID columns we care about

    out <- merge(x, y, by = "UniqID", all = TRUE, suffixes = c(".x", ".y"))

	#out <- apply(out,2,as.character)

    # Coalesce duplicated ID columns if present
    for (cid in keep_idcols) {
      cx <- paste0(cid, ".x")
      cy <- paste0(cid, ".y")
      if (cx %in% names(out) && cy %in% names(out)) {
        # prefer non-NA from .x, else .y
        out[[cid]] <- ifelse(!is.na(out[[cx]]), out[[cx]], out[[cy]])
        out[[cx]] <- NULL
        out[[cy]] <- NULL
      }
    }

    # Also drop any leftover ".x"/".y" duplicates for columns that are exactly identical
    dup_pairs <- grep("\\.(x|y)$", names(out), value = TRUE)
    if (length(dup_pairs)) {
      base_nms <- unique(sub("\\.(x|y)$", "", dup_pairs))
      for (bn in base_nms) {
        cx <- paste0(bn, ".x")
        cy <- paste0(bn, ".y")
        if (cx %in% names(out) && cy %in% names(out)) {
          if (identical(out[[cx]], out[[cy]])) {
            out[[bn]] <- out[[cx]]
            out[[cx]] <- NULL
            out[[cy]] <- NULL
          }
        }
      }
    }

    out
  }

  Fit_Out_Tab <- NULL

  for (nTrt in seq_len(nTraits)) {
    traitME  <- TraitME[nTrt]
    Fit_Pred <- Pred_Out[[nTrt]]

    # 2) Guard: models vs provided tables
    nMod <- length(Fit_Pred)
    if (nMod == 0L) stop(sprintf("No prediction tables for trait '%s'.", traitME), call. = FALSE)
    if (nMod != length(Models)) {
      warning(sprintf("Trait '%s': number of model tables (%d) != length(Models) (%d). Will label by index.",
                      traitME, nMod, length(Models)))
    }

    # Build per-trait combined table
    Fit_Out <- NULL
    for (nM in seq_len(nMod)) {
      lab <- if (nM <= length(Models)) Models[nM] else paste0("M", nM)

      tab <- Fit_Pred[[nM]]
      if (!is.data.frame(tab)) stop(sprintf("Trait '%s' model %d is not a data.frame.", traitME, nM), call. = FALSE)
      tab$Strain <- as.character(tab$Strain)
	  tab$Env <- as.character(tab$Env)

      tab <- .rename_to_uniqid(tab)
      tab <- .tag_obs_pred(tab, lab)

      if (is.null(Fit_Out)){
        Fit_Out <- tab
      } else {
        Fit_Out <- .safe_merge_by_uid(Fit_Out, tab)
      }
    }

    # Prefix Obs/Pred with trait name (as requested)
    names(Fit_Out) <- sub("^Obs-",  paste0(traitME, "Obs-"),  names(Fit_Out))
    names(Fit_Out) <- sub("^Pred-", paste0(traitME, "Pred-"), names(Fit_Out))

    # Merge across traits (by UniqID)
    if (is.null(Fit_Out_Tab)) {
      Fit_Out_Tab <- Fit_Out
    } else {
      Fit_Out_Tab <- .safe_merge_by_uid(Fit_Out_Tab, Fit_Out)
    }
  }

  envInd <- which(IDColsME %in% envVar)
  IDColsME[envInd] <- "Env"

  # 3) Final selection: ensure selOutCols all exist, else error
  valCols   <- grep("Obs-|Pred-", names(Fit_Out_Tab), value = TRUE)
  selOutCols <- c("UniqID", setdiff(IDColsME, IDColME), valCols)

  missing_cols <- setdiff(selOutCols, names(Fit_Out_Tab))


  if (length(missing_cols)) {

	selOutFiltCols <- setdiff(selOutCols,missing_cols)
	# stop(sprintf(
      # "Cannot subset final table: the following columns are missing: %s\nAvailable columns: %s",
      # paste(missing_cols, collapse = ", "),
      # paste(names(Fit_Out_Tab), collapse = ", ")
    # ), call. = FALSE)
  }else{ selOutFiltCols <- selOutCols }
  if(length(selOutFiltCols)){
     Fit_Out_Tab_Sel <- Fit_Out_Tab[, selOutFiltCols, drop = FALSE]
  }else {Fit_Out_Tab_Sel <- Fit_Out_Tab}

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
#' @export
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
#' @export
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
#' @export
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

####
getModelPredictiveAbilityPlots <- function(outputListME,TraitME,IDColsME,IDColME,fitEnvCov,nModel,outDirPath){

 print("PlotFnCall_In")
 library(dplyr)
 if(fitEnvCov==FALSE){
  Models <- c("MM","MDs","MDe")
 }else if(fitEnvCov == TRUE){
  Models <- c("EMM","EMDs","EMDe")
 }
  nTraits <- length(TraitME)

  nTrt <- 1
  Pred_Out <- outputListME$preds
  traitME <- TraitME[nTrt]
  Fit_Pred <- Pred_Out[[nTrt]]
  nMod <- length(Fit_Pred)
  nM <- nModel

  Fit_Out <- as.data.frame(Fit_Pred[[nM]])
  uniqIDCol <- which(colnames(Fit_Out) == IDColME)
  colnames(Fit_Out)[uniqIDCol] <- "UniqID"

  print(paste("Plot_Fit_Out-",dim(Fit_Out)))
  print(head(Fit_Out))

  FitValues_Summary_DF <- Fit_Out %>% dplyr::group_by(Strain) %>% dplyr::summarise(AxLoc_Avg.Obs = mean(Obs,na.rm=TRUE),AxLoc_Avg.Pred = mean(Pred,na.rm=TRUE),.groups = "drop")
  print(paste("Plot_DF_Dims-",dim(FitValues_Summary_DF)))
  print(colnames(FitValues_Summary_DF))
  print(head(FitValues_Summary_DF))

  mainTitle <- paste(Models[nM]," Model Predictive Ability",sep="")

  p <- plot_model_eval_strains_marginal(DF=FitValues_Summary_DF,
									title_size = 18,
									axis_title_size = 18,
									axis_text_size = 14,
									text_size = 5.0,
									point_size= 6,
									point_color = "red",
									line_color = "blue",
								    main_title =  mainTitle
								)
  OFN <- paste(outDirPath,"/",Models[nM],"_Model_PredAbility_PercentileGrid.png",sep="")
  png(OFN,w=760, h=760, pointsize=20)
  print(p)
  Sys.sleep(5)
  dev.off()
  return(p)

   }

####
### V4
plot_model_eval_strains_marginal <- function(DF,
                                             strain_col = Strain,
                                             obs_col = AxLoc_Avg.Obs,
                                             pred_col = AxLoc_Avg.Pred,
                                             probs = c(0, 0.20, 0.50, 0.80, 1),
                                             point_alpha = 0.5,
                                             point_size = 4,
                                             point_color = "black",
                                             line_color = "red",
                                             text_size = 4,
                                             title_size = 16,
										 main_title = "Model Predictive Ability" ,
                                             axis_title_size = 14,
                                             axis_text_size = 12) {
  library(dplyr)
  library(ggplot2)
  library(grid)

  # Rename columns
  DF <- DF %>% dplyr::rename(
      Strain = {{ strain_col }},
      observed = {{ obs_col }},
      predicted = {{ pred_col }}
    )

  # Percentile breaks
  pred_breaks <- quantile(DF$predicted, probs = probs, na.rm = TRUE)
  obs_breaks  <- quantile(DF$observed,  probs = probs, na.rm = TRUE)
  labels <- paste0(100 * probs[-length(probs)], "-", 100 * probs[-1], "%")

  # Bin assignment
  DF <- DF %>%
    dplyr::mutate(
      pred_bin = cut(predicted, breaks = pred_breaks, include.lowest = TRUE, labels = labels),
      obs_bin  = cut(observed,  breaks = obs_breaks, include.lowest = TRUE, labels = labels)
    )

  # Proportion within each pred_bin (conditional proportions)
  DF_summary <- DF %>%
    dplyr::group_by(pred_bin, obs_bin) %>%
    dplyr::summarise(unique_strains = n_distinct(Strain), .groups = "drop") %>%
    dplyr::group_by(pred_bin) %>%
    dplyr::mutate(proportion = unique_strains / sum(unique_strains)) %>%
    dplyr::ungroup()

  # Midpoints for placing text
  pred_mids <- sapply(1:(length(pred_breaks)-1), function(i) mean(pred_breaks[i:(i+1)]))
  obs_mids  <- sapply(1:(length(obs_breaks)-1),  function(i) mean(obs_breaks[i:(i+1)]))+2

  DF_summary <- DF_summary %>%
    dplyr::mutate(
      pred_mid = pred_mids[as.numeric(pred_bin)],
      obs_mid  = obs_mids[as.numeric(obs_bin)],
      x_min = pred_breaks[as.numeric(pred_bin)],
      x_max = pred_breaks[as.numeric(pred_bin)+1],
      y_min = obs_breaks[as.numeric(obs_bin)],
      y_max = obs_breaks[as.numeric(obs_bin)+1]
    )

  # Global statistics
  lm_fit <- lm(observed ~ predicted, data = DF)
  mse <- mean((DF$observed - DF$predicted)^2,na.rm=TRUE)
  cor_val <- cor(DF$predicted, DF$observed,use = "pairwise.complete")

  stats_text <- paste0(
    "Linear Fit: y = ", round(coef(lm_fit)[1], 2), " + ",
    round(coef(lm_fit)[2], 2), "x\n",
    "MSE = ", round(mse, 3), "\n",
    "Correlation = ", round(cor_val, 3)
  )

  # Plot
  p <- ggplot(DF, aes(x = predicted, y = observed)) +
    geom_rect(data = DF_summary,
              aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max,
                  fill = proportion),
              inherit.aes = FALSE, alpha = 0.5) +
    scale_fill_gradient(low = "white", high = "steelblue", name = "Proportion within \n Pred Bin") +
    geom_point(alpha = point_alpha,size=point_size,color=point_color) +
    geom_vline(xintercept = pred_breaks[c(1,2,4,5)], linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = pred_breaks[3], linetype = "solid", color = "black") +
    geom_hline(yintercept = obs_breaks[c(1,2,4,5)], linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = obs_breaks[3], linetype = "solid", color = "black") +
    geom_text(data = DF_summary,
              aes(x = pred_mid, y = obs_mid,
                  label = paste0(round(100 * proportion, 1), "%")),
              color = "black", size = text_size, inherit.aes = FALSE) +
    geom_smooth(method = "lm", color = line_color, se = FALSE) +
	annotate("text",
             x = min(DF$predicted) + 0.05 * diff(range(DF$predicted)),
             y = max(DF$observed,na.rm=TRUE) - 0.05 * diff(range(DF$observed,na.rm=TRUE)),
             label = stats_text, hjust = 0, vjust = 1, size = text_size, color = "black") +
    annotate("text", x = pred_mids, y = min(DF$observed,na.rm=TRUE) - 0.15 * diff(range(DF$observed,na.rm=TRUE)),
             label = labels, angle = 0, vjust = 1,size = text_size, fontface = "bold") +
    annotate("text", y = obs_mids, x = min(DF$predicted) - 0.1 * diff(range(DF$predicted)),
             label = labels, angle = 90, hjust = 1,size = text_size, fontface = "bold") +
    labs(
      title = main_title,
      x = "Predicted Values (percentile bins)",
      y = "Observed Values (percentile bins)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = title_size, face = "bold"),
      axis.title = element_text(size = axis_title_size, face = "bold"),
      axis.text = element_text(size = axis_text_size),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  return(p)
}


###########################################################################################
## Part B: New functions added from GS_Pipeline_Jan_2026_FnsApp.R
###########################################################################################


VCFtoDF_V2<- function(infile){


  vcf <- read.vcfR(infile, verbose = FALSE)
  gt <- extract.gt(vcf, element = "GT", as.numeric = F,IDtoRowNames = FALSE)
  fix_T <- as_tibble(getFIX(vcf))

  gt2 <- matrix(0, ncol = ncol(gt), nrow = nrow(gt))
  colnames(gt2) <- colnames(gt)

  gt2a <- apply(gt,2, function(x) gsub("1/1","BB",x))
  gt2b <- gsub("0[/|]0","AA",gt2a)
  gt2c <- gsub("[10][/|][10]","AB",gt2b)
  gt2d <- gsub("\\.[/|]\\.","NA",gt2c)

  # gt2[(gt == "1/1")|(gt == "1|1")] <- 'BB'
  # gt2[(gt == "0/0")|(gt == "0|0")] <- 'AA'
  # gt2[(gt == "0/1")|(gt == "1/0")|(gt == "0|1")|(gt == "1|0")] <- 'AB'
  # gt2[(gt == "\\./\\.")|(gt == "\\.|\\.")] <- NA

  gt2 <- as_tibble(gt2d) %>% mutate(SNPID = fix_T$ID)

  gt.simple <- fix_T %>% dplyr::select("ID", "CHROM", POS, REF, ALT) %>%
    dplyr::rename(SNPID=ID, Chrom=CHROM, Position=POS) %>%
    dplyr::left_join(gt2, by = 'SNPID')

  return(gt.simple)
}

##

VCFtoDF_NAM <- function(infile){

  vcf <- read.vcfR(infile, verbose = FALSE)
  gt <- extract.gt(vcf, element = "GT", as.numeric = F)
  fix_T <- as_tibble(getFIX(vcf))

  gt2 <- matrix(0, ncol = ncol(gt), nrow = nrow(gt))
  colnames(gt2) <- colnames(gt)


  gt2a <- apply(gt,2, function(x) gsub("1/1","BB",x))
  gt2b <- gsub("0[/|]0","AA",gt2a)
  gt2c <- gsub("[10][/|][10]","AB",gt2b)
  gt2d <- gsub("\\.[/|]\\.","NA",gt2c)




  gt2a <- apply(gt,2, function(x) gsub("1/1","BB",x))
  gt2b <- gsub("0[/|]0","AA",gt2a)
  gt2c <- gsub("[10][/|][10]","AB",gt2b)
  gt2d <- gsub("\\.[/|]\\.","NA",gt2c)

  # gt2[(gt == "1/1")|(gt == "1|1")] <- 'BB'
  # gt2[(gt == "0/0")|(gt == "0|0")] <- 'AA'
  # gt2[(gt == "0/1")|(gt == "1/0")|(gt == "0|1")|(gt == "1|0")] <- 'AB'
  # gt2[(gt == "\\./\\.")|(gt == "\\.|\\.")] <- NA

  gt2 <- as_tibble(gt2d) %>%
    mutate(SNPID = fix_T$ID)

  gt.simple <- fix_T %>% dplyr::select("ID", "CHROM", POS, REF, ALT) %>%
    dplyr::rename(SNPID=ID, Chrom=CHROM, Position=POS) %>%
    dplyr::left_join(gt2, by = 'SNPID')

  return(gt.simple)
}

####

getRankedPredictedValues_V2 <- function(Data_Table_Num_Filt_List,nTraits,trait,GPModel,optTS=NULL){

	 TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	 TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]


	 TestGenoTable <- apply(TestData_Table_Num_Filt[,-1],2,as.numeric)
	 initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]

	 if(is.null(optTS)){

		 trainPheno0 <- as.numeric(as.character(TrainData_Table_Num_Filt[,trait]))
		 geno_012 <- apply(TrainData_Table_Num_Filt[,c(initCol:(ncol(TrainData_Table_Num_Filt)-nTraits))],2,as.numeric)
		 trainGeno <- apply(geno_012,2,function(x) x-1)
		 Geno <- as.character(TrainData_Table_Num_Filt[,1])
	 }

	 if(!is.null(optTS)){
	     strainID <- as.character(TrainData_Table_Num_Filt[,1])
		 trainSetID <- as.character(as.vector(optTS))
	     trainIndices <- which(strainID %in% trainSetID)
		 trainPheno0 <- as.numeric(as.character(TrainData_Table_Num_Filt[trainIndices,trait]))
		 geno_012 <- apply(TrainData_Table_Num_Filt[trainIndices,c(initCol:(ncol(TrainData_Table_Num_Filt)-nTraits))],2,as.numeric)
		 trainGeno <- apply(geno_012,2,function(x) x-1)
		 Geno <- as.character(TrainData_Table_Num_Filt[trainIndices,1])
	 }


	 if(anyNA(trainGeno)){
		trainGeno_Imp0 <- snpQC(trainGeno,impute=TRUE,remove=FALSE)
	 }else{trainGeno_Imp0 <- trainGeno}


	  if(anyNA(TestGenoTable)){
	    testGeno_Imp0 <- snpQC(TestGenoTable,impute=TRUE,remove=FALSE)
		testGeno_Imp1 <- apply(testGeno_Imp0,2,function(x) as.numeric(x-1))
		testNA <- apply(testGeno_Imp,2,function(x) length(which(is.na(x))))
		testNAIndices <-(which(unlist(testNA) !=0))
	 }else{testGeno_Imp1 <- apply(TestGenoTable,2,function(x) as.numeric(x-1))
	    testNAIndices <- NULL
	 }

### Remove lines with NA
	 if(anyNA(trainPheno0)){
	   ph_NA_Indices <- which(is.na(trainPheno0))
	   trainPheno <- trainPheno0[-ph_NA_Indices]
	   trainGeno_Imp1 <- trainGeno_Imp0[-ph_NA_Indices,]
	   Geno <- as.character(TrainData_Table_Num_Filt[-ph_NA_Indices,1])

	 }
	 if(!anyNA(trainPheno0)){

	   trainPheno <- trainPheno0
	   trainGeno_Imp1 <- trainGeno_Imp0
	   Geno <- as.character(TrainData_Table_Num_Filt[,1])
	 }


### Check trainGeno and testGeno contain the same marker set
    if(!is.null(testNAIndices)){
	  trainssIDs <- colnames(trainGeno_Imp1)
	  testssIDs <- colnames(testGeno_Imp1)
	   if(length(trainssIDs) > length(testssIDs)){
		 NAssID <- setdiff(trainssIDs,testssIDs)
		 NAssIDIndex <- which(trainssIDs %in% NAssID)
		 trainGeno_Imp <-  trainGeno_Imp1[,-NAssIDIndex]
		 testGeno_Imp <-  testGeno_Imp1[,-testNAIndices]
	   }
	   if(length(trainssIDs) < length(testssIDs)){
		   NAssID <- setdiff(testssIDs,trainssIDs)
		   NAssIDIndex <- which(testssIDs %in% NAssID)
		   testGeno_Imp <-  testGeno_Imp1[,-c(NAssIDIndex,testNAIndices)]
		   trainGeno_Imp <-  trainGeno_Imp1
	   }
	   if(length(trainssIDs) == length(testssIDs)){
	       if(length(which(trainssIDs %in% testssIDs)) == length(trainssIDs)){
			   testGeno_Imp <-  testGeno_Imp1
			   trainGeno_Imp <-  trainGeno_Imp1
		   }
		   if(length(which(trainssIDs %in% testssIDs)) != length(trainssIDs)){
		       commonTestIndices <- which(testssIDs %in% trainssIDs)
			   commonTrainIndices <- which(trainssIDs %in% testssIDs)
			   testGeno_Imp <-  testGeno_Imp1[,commonTestIndices]
			   trainGeno_Imp <-  trainGeno_Imp1[,commonTrainIndices]
		   }
	   }

	}

	if(is.null(testNAIndices) | length(testNAIndices)==0){

	     trainssIDs <- colnames(trainGeno_Imp1)
	     testssIDs <- colnames(testGeno_Imp1)
	   if(length(trainssIDs) > length(testssIDs)){
		 NAssID <- setdiff(trainssIDs,testssIDs)
		 NAssIDIndex <- which(trainssIDs %in% NAssID)
		 trainGeno_Imp <-  trainGeno_Imp1[,-NAssIDIndex]
		 testGeno_Imp <-  testGeno_Imp1[,-testNAIndices]
	   }
	   if(length(trainssIDs) < length(testssIDs)){
		   NAssID <- setdiff(testssIDs,trainssIDs)
		   NAssIDIndex <- which(testssIDs %in% NAssID)
		   testGeno_Imp <-  testGeno_Imp1[,-NAssIDIndex]
		   trainGeno_Imp <-  trainGeno_Imp1
	   }
	   if(length(trainssIDs) == length(testssIDs)){
	       if(length(which(trainssIDs %in% testssIDs)) == length(trainssIDs)){
			   testGeno_Imp <-  testGeno_Imp1
			   trainGeno_Imp <-  trainGeno_Imp1
		   }
		   if(length(which(trainssIDs %in% testssIDs)) != length(trainssIDs)){
		       commonTestIndices <- which(testssIDs %in% trainssIDs)
			   commonTrainIndices <- which(trainssIDs %in% testssIDs)
			   testGeno_Imp <-  testGeno_Imp1[,commonTestIndices]
			   trainGeno_Imp <-  trainGeno_Imp1[,commonTrainIndices]
		   }
	   }

	 }



## Kinship Matrix
	 # A <- A.mat(trainGeno_Imp)
	 # colnames(A) <- Geno
	 # rownames(A) <- Geno

## Prepare Data Table for GP
	 Data <- cbind.data.frame(Geno,trainPheno)
	 colnames(Data) <- c("Geno","Pheno")
	 Geno <- "Geno"
	 Pheno <- "Pheno"




#### Impute trainGeno and train using mixed.solve


	if(GPModel == "rrBLUP (rrBLUP)"){

	 pred <- mixed.solve(trainPheno,Z=trainGeno_Imp,SE=FALSE,return.Hinv =FALSE)
	 Mean <- as.numeric(pred$beta)
	 Effects <- pred$u
     PredictedValues <- Mean + (trainGeno_Imp %*% Effects)
     SortedPredictedValues <- sort.int(PredictedValues,decreasing=TRUE,index.return=TRUE)

	 Test_PredictedValues <-  Mean + (testGeno_Imp %*% Effects)

	 Test_SortedPredictedValues <- sort.int(Test_PredictedValues,decreasing=TRUE,index.return=TRUE)
	 Test_StrainID <- rownames(TestData_Table_Num_Filt)

	 }


	 if(GPModel == "rrBLUP (bWGR)"){

		 pred <- emRR(trainPheno,trainGeno_Imp)
		 Mean <- as.numeric(pred$mu)
		 Effects <- pred$b


		 PredictedValues <- Mean + (trainGeno_Imp %*% Effects)
		 SortedPredictedValues <- sort.int(PredictedValues,decreasing=TRUE,index.return=TRUE)



		 Test_PredictedValues <-  Mean + (testGeno_Imp %*% Effects)


		 Test_SortedPredictedValues <- sort.int(Test_PredictedValues,decreasing=TRUE,index.return=TRUE)
		 Test_StrainID <- rownames(TestData_Table_Num_Filt)

	 }


	  if(GPModel == "BayesB (bWGR)"){

		 pred <- emBB(trainPheno,trainGeno_Imp)
		 Mean <- as.numeric(pred$mu)
		 Effects <- pred$b


		 PredictedValues <- Mean + (trainGeno_Imp %*% Effects)
		 SortedPredictedValues <- sort.int(PredictedValues,decreasing=TRUE,index.return=TRUE)
		 Test_PredictedValues <-  Mean + (testGeno_Imp %*% Effects)


		 Test_SortedPredictedValues <- sort.int(Test_PredictedValues,decreasing=TRUE,index.return=TRUE)
		 Test_StrainID <- rownames(TestData_Table_Num_Filt)
      }

	  if(GPModel == "BayesLASSO (bWGR)"){

		 pred <- emBL(trainPheno,trainGeno_Imp)
		 Mean <- as.numeric(pred$mu)
		 Effects <- pred$b

		 PredictedValues <- Mean + (trainGeno_Imp %*% Effects)
		 SortedPredictedValues <- sort.int(PredictedValues,decreasing=TRUE,index.return=TRUE)

		 Test_PredictedValues <-  Mean + (testGeno_Imp %*% Effects)
	     Test_SortedPredictedValues <- sort.int(Test_PredictedValues,decreasing=TRUE,index.return=TRUE)
		 Test_StrainID <- rownames(TestData_Table_Num_Filt)
      }


###
	# pred.kb <- kin.blup(as.data.frame(Data),Geno,Pheno,GAUSS=FALSE,K=A,covariate=NULL,PEV=TRUE,n.core=1,theta.seq=NULL)

### Upper bound of Reliability
    # if(length(testNAIndices)>0){
       # trainGeno_Imp2 <- apply(trainGeno_Imp[,-testNAIndices],2,function(x) x+1)
	# }
	# if(length(testNAIndices)==0){
	    # trainGeno_Imp2 <- apply(trainGeno_Imp,2,function(x) x+1)
	# }
    cleanData <- cleanREP(trainPheno,trainGeno_Imp)
    M <-  cleanData[[2]]
    M.Pdt <- t(M)%*% solve(M %*% t(M) + diag(nrow(M))) %*% M

    getU <- function(M.Pdt,v){
	   v.hat <- M.Pdt %*% v
	   U <- (t(v.hat)%*% v.hat)/ (t(v) %*% v)
	   return(U)
    }

	# if(length(testNAIndices)>0){
      # U <- apply(testGeno_Imp[,-testNAIndices],1,function(x) getU(M.Pdt,x))
	# }

	  U <- apply(testGeno_Imp,1,function(x) getU(M.Pdt,x))


	Test_SortedIndices <- Test_SortedPredictedValues[[2]]
	Sorted_Test_StrainID <- Test_StrainID[Test_SortedIndices]
	U_Sorted <- U[(Test_SortedIndices)]

	outputDF <- cbind.data.frame(Sorted_Test_StrainID ,round(Test_SortedPredictedValues[[1]],digits=2),round(U_Sorted,digits=5))
	colnames(outputDF) <- c("LineID",paste("Predicted Value for ",trait,sep=""),"Upper Bound of Reliability")


   return(outputDF)

 }


######

fixRefMismatches <- function(GenoTable_In, sampleInd){

	message("Input dimensions: ", paste(dim(GenoTable_In), collapse = " x "))

	## Get column indices for CHROM, POS, REF, ALT (both .x and .y)
	InfoInd <- unlist(lapply(c("CHROM","POS","REF","ALT"), function(x) grep(x, colnames(GenoTable_In))))
	RefxInd <- InfoInd[5]
	RefyInd <- InfoInd[6]
	AltxInd <- InfoInd[7]
	AltyInd <- InfoInd[8]
	colnames(GenoTable_In)[c(RefxInd,AltxInd)] <- c("REF.x","ALT.x")

	# Identify various match/mismatch types
	REFMatches_a1_a4_Ind <- which(GenoTable_In[,RefxInd] == GenoTable_In[,RefyInd])
	REFMatches_ALT_a1_a4_Ind <- which(GenoTable_In[,RefxInd] == GenoTable_In[,AltyInd])
	ALTMatches_REF_a1_a4_Ind <- which(GenoTable_In[,AltxInd] == GenoTable_In[,RefyInd])
	REF_ALT_Swap_a1_a4_Ind <- intersect(REFMatches_ALT_a1_a4_Ind, ALTMatches_REF_a1_a4_Ind)

	message("REF matches: ", length(REFMatches_a1_a4_Ind))
	message("REF == ALT.y matches: ", length(REFMatches_ALT_a1_a4_Ind))
	message("ALT == REF.y matches: ", length(ALTMatches_REF_a1_a4_Ind))
	message("REF/ALT swaps needed: ", length(REF_ALT_Swap_a1_a4_Ind))

	# Identify mismatches for filtering
	REFMisMatches_a1_a4_Ind0 <- which(GenoTable_In[,RefxInd] != GenoTable_In[,RefyInd])
	REFMatches_Other_a1_a4_Ind <- setdiff(REFMisMatches_a1_a4_Ind0, ALTMatches_REF_a1_a4_Ind)

	message("Total REF mismatches: ", length(REFMisMatches_a1_a4_Ind0))
	message("Non-swap mismatches: ", length(REFMatches_Other_a1_a4_Ind))

	## Fix REF/ALT mismatches
	GenoTable_filt_mod <- GenoTable_In
	swapCorrectedInd <- c()

    if(length(REF_ALT_Swap_a1_a4_Ind)>0){
	# Swap REF/ALT fields
	tmpRefToAlt <- GenoTable_filt_mod[REF_ALT_Swap_a1_a4_Ind, AltxInd]
	GenoTable_filt_mod[REF_ALT_Swap_a1_a4_Ind, AltxInd] <- GenoTable_filt_mod[REF_ALT_Swap_a1_a4_Ind, RefxInd]
	GenoTable_filt_mod[REF_ALT_Swap_a1_a4_Ind, RefxInd] <- tmpRefToAlt

	# Genotype re-encoding

	swap_genotypes <- function(gt) {
		gt <- gsub("0 / 0", "BB", gt)
		gt <- gsub("1 / 1", "0 / 0", gt)
		gt <- gsub("BB", "1 / 1", gt)
		gt <- gsub("1 / 0", "AB", gt)
		gt <- gsub("0 / 1", "1 / 0", gt)
		gt <- gsub("AB", "1 / 0", gt)
		return(gt)
	}

	Gt_for_Swap <- GenoTable_filt_mod[REF_ALT_Swap_a1_a4_Ind, sampleInd]
	Gt_Swapped <- apply(Gt_for_Swap, 2, swap_genotypes)
	GenoTable_filt_mod[REF_ALT_Swap_a1_a4_Ind, sampleInd] <- Gt_Swapped

    REF_ALT_Swap_a1_a4_Corr_Ind <- (which(GenoTable_filt_mod[REF_ALT_Swap_a1_a4_Ind,"REF.x"] == GenoTable_filt_mod[REF_ALT_Swap_a1_a4_Ind,"REF.y"]))

    message("REF/ALT swaps corrected: ", length(REF_ALT_Swap_a1_a4_Corr_Ind))

    swapCorrectedInd <- REF_ALT_Swap_a1_a4_Ind[REF_ALT_Swap_a1_a4_Corr_Ind]

	}

	## Strand flips (complementary base matching)
	flipInd <- which(
		(GenoTable_filt_mod[,RefxInd]=="A" & GenoTable_filt_mod[,RefyInd]=="T") |
		(GenoTable_filt_mod[,RefxInd]=="T" & GenoTable_filt_mod[,RefyInd]=="A") |
		(GenoTable_filt_mod[,RefxInd]=="C" & GenoTable_filt_mod[,RefyInd]=="G") |
		(GenoTable_filt_mod[,RefxInd]=="G" & GenoTable_filt_mod[,RefyInd]=="C")
	)
	message("Strand flip candidates: ", length(flipInd))

   flips_corrected_ind <- c()

   if(length(flipInd)>0){

	FlipAlleleMat_1 <- GenoTable_filt_mod[flipInd,]
	FlipAlleleMat_1[,"REF.x"] <- chartr("TAGC", "ATCG", FlipAlleleMat_1[,"REF.x"])
	FlipAlleleMat_1[,"ALT.x"] <- chartr("TAGC", "ATCG", FlipAlleleMat_1[,"ALT.x"])
	GenoTable_filt_mod[flipInd,] <- FlipAlleleMat_1


    flips_corr_Ind <- (which(GenoTable_filt_mod[flipInd,"REF.x"] == GenoTable_filt_mod[flipInd,"REF.y"]))
	message("Flips corrected: ", length(flips_corr_Ind))


	flips_remng_ind <- flipInd[setdiff(seq_len(length(flipInd)),flips_corr_Ind)]
    flips_corrected_ind <- flipInd[flips_corr_Ind]

   }

	## Flip-Swap
    RefMismatchInd_for_FlSw <- which(GenoTable_filt_mod[,RefxInd] != GenoTable_filt_mod[,RefyInd])
    flipSwapInd <- setdiff(RefMismatchInd_for_FlSw,unique(c(swapCorrectedInd,flips_corrected_ind)))

	message("Flip-swap candidates: ", length(flipSwapInd))

   refmatchIndFlipSwap <- c()
   if(length(flipSwapInd)>0){

	AltCompBaseMat_1 <- GenoTable_filt_mod[flipSwapInd,]

	TmpFlSwp_in_REF <- AltCompBaseMat_1[,"REF.x"]
	TmpFlSwp_in_ALT <- AltCompBaseMat_1[,"ALT.x"]


	AltCompBaseMat_1[,"REF.x"] <- chartr("TAGC", "ATCG", AltCompBaseMat_1[,"REF.x"])
	AltCompBaseMat_1[,"ALT.x"] <- chartr("TAGC", "ATCG", AltCompBaseMat_1[,"ALT.x"])

	###

	AltCompBaseMat_2 <- AltCompBaseMat_1

	##
	tmpNt <- AltCompBaseMat_2[,"REF.x"]
	AltCompBaseMat_2[,"REF.x"] <- AltCompBaseMat_2[,"ALT.x"]
	AltCompBaseMat_2[,"ALT.x"] <- tmpNt

    ## keep only those flip-swap markers that match the new ref
	refmatchIndFlipSwap <- which(AltCompBaseMat_2[,"REF.x"] == AltCompBaseMat_2[,"REF.y"])
	rmFlipSwapInd <- setdiff(seq_len(nrow(AltCompBaseMat_2)), refmatchIndFlipSwap)

	##
	AltCompBaseMat_2[rmFlipSwapInd,"REF.x"] <- TmpFlSwp_in_REF[rmFlipSwapInd]
	AltCompBaseMat_2[rmFlipSwapInd,"ALT.x"] <- TmpFlSwp_in_ALT[rmFlipSwapInd]

	Gt_FlipSwapped <- apply(AltCompBaseMat_2[refmatchIndFlipSwap,sampleInd], 2, swap_genotypes)
	AltCompBaseMat_2[refmatchIndFlipSwap,sampleInd] <- Gt_FlipSwapped

	GenoTable_filt_mod[flipSwapInd,] <- AltCompBaseMat_2

    flipswap_corr_Ind <- (which(GenoTable_filt_mod[flipSwapInd,"REF.x"] == GenoTable_filt_mod[flipSwapInd,"REF.y"]))
	message("Flipswaps corrected Ref: ", length(flipswap_corr_Ind))
    }

	# Remove mismatches still present after harmonization

	# Remove mismatches still present after harmonization
	flipswapMismatchesRmInd <- which(GenoTable_filt_mod[,"REF.x"] != GenoTable_filt_mod[,"REF.y"])
	message("Remaining REF mismatches after all corrections: ", length(flipswapMismatchesRmInd))
    if(length(flipswapMismatchesRmInd) >0){
	  GenoTable_Out <- GenoTable_filt_mod[-flipswapMismatchesRmInd,]
	}else{
	  GenoTable_Out <- GenoTable_filt_mod
	}

	# After fixing mismatches, report summary
	total_snps <- nrow(GenoTable_In)
	kept_snps <- nrow(GenoTable_Out)
	#removed_snps <- length(flipswapMismatchesRmInd)
	removed_snps <- (total_snps-kept_snps)
	n_swap <- length(REF_ALT_Swap_a1_a4_Ind)
	n_flip <- length(flipInd)
	#n_flipswap <- length(flipSwapInd)
	n_flipswap <- length(refmatchIndFlipSwap)
	n_match <- length(REFMatches_a1_a4_Ind)

	meta_info <- list(
		"Total SNPs input" = total_snps,
		"REF matches (kept as-is)" = n_match,
		"REF/ALT swapped" = n_swap,
		"Strand flipped" = n_flip,
		"Flip-swapped" = n_flipswap,
		"Remaining mismatches removed" = removed_snps,
		"SNPs in final output" = kept_snps
	)

	message("=== Harmonization Summary ===")
	for (key in names(meta_info)) {
	  message(sprintf("%-35s : %d", key, meta_info[[key]]))
	}

	return(list(
	  GenoTable_Harmonized = GenoTable_Out,
	  Harmonization_Report = meta_info
	))

}

######

getAGPVCF_Data <- function(AGP_Pop_Geno_Filt){

AGP_Pop_Geno_Filt$CHROM <- unlist(lapply(strsplit(rownames(AGP_Pop_Geno_Filt),"_"),function(x) x[4]))
AGP_Pop_Geno_Filt$POS <- unlist(lapply(strsplit(rownames(AGP_Pop_Geno_Filt),"_"),function(x) x[5]))
AGP_Pop_Geno_Filt$REF <- unlist(lapply(strsplit(rownames(AGP_Pop_Geno_Filt),"_"),function(x) x[6]))
AGP_Pop_Geno_Filt$ALT <- unlist(lapply(strsplit(rownames(AGP_Pop_Geno_Filt),"_"),function(x) x[7]))
AGP_Pop_Geno_Filt$ChrPosV1 <- paste(AGP_Pop_Geno_Filt$CHROM,AGP_Pop_Geno_Filt$POS,sep="-")

####

BARC_markerIDs <- grep("BARC",rownames(AGP_Pop_Geno_Filt),value=TRUE)
BARC_markerIDsInd <- grep("BARC",rownames(AGP_Pop_Geno_Filt))
trt_markerIDs <- rownames(AGP_Pop_Geno_Filt)[c((BARC_markerIDsInd[length(BARC_markerIDsInd)]+1):length(rownames(AGP_Pop_Geno_Filt)))]
trt_markerIDsInd <- c((BARC_markerIDsInd[length(BARC_markerIDsInd)]+1):length(rownames(AGP_Pop_Geno_Filt)))

### Extract INFO Cols in a1.v1

INFO_Cols1 <- as.data.frame(do.call(rbind,lapply(strsplit(BARC_markerIDs,"_"),function(x) x[c(4:7)])))
colnames(INFO_Cols1) <- c("CHROM","POS","REF","ALT")
INFO_Cols1$ChrPosV1 <- paste(INFO_Cols1$CHROM,INFO_Cols1$POS,sep="-")

INFO_Cols <- as.data.frame(do.call(rbind,lapply(strsplit(BARC_markerIDs,"_"),function(x) x[c(4:7)])))
colnames(INFO_Cols) <- c("CHROM","POS","REF","ALT")
INFO_Cols$ChrPosV1 <- paste(INFO_Cols$CHROM,INFO_Cols$POS,sep="-")


###
rownames(AGP_Pop_Geno_Filt)[c(length(rownames(AGP_Pop_Geno_Filt))-1,length(rownames(AGP_Pop_Geno_Filt)))]
AGP_Pop_Geno_Filt[trt_markerIDsInd,1:5]

AGP_Pop_Geno_Filt_Mod <- merge(INFO_Cols1,AGP_Pop_Geno_Filt,by="ChrPosV1")
colnames(AGP_Pop_Geno_Filt_Mod)[1] <- "CHROM.POS"

#### ssIDs are not provided in the AGP Genotype File. ssIDs are extracted from NUST BARCSoySNP6K
#### data.
#### Read InfoCols from NUST_Geno.vcf BARCSoySNP6K

infileVCF <- "NUST_Geno.vcf"
# Set the file path to the VCF file
vcf_file <- infileVCF


# Open a connection to the VCF file and skip the meta info lines
vcf_conn <- file(vcf_file, open = "r")
header <- c()
while (TRUE) {
  line <- readLines(vcf_conn, n = 1)
  if (startsWith(line, "##")) {
    header<- c(header,line)
    next
  } else {
    # Once we reach the information row, break out of the loop
    info_row <- line
    break
  }
}

# Read in the data in the VCF file, starting from the information row
NUST_Geno_vcf_data <- read.table(vcf_conn, header = FALSE)
close(vcf_conn)

# Print the information row and the first few rows of the VCF data
colnames(NUST_Geno_vcf_data) <- unlist(strsplit(info_row,"\t"))
#head(NUST_Geno_vcf_data)

colnames(NUST_Geno_vcf_data)[1]<- "CHROM"

#####

NUST_Geno_IDs <- NUST_Geno_vcf_data[,"ID"]

NUST_Geno_vcf_data <- droplevels(NUST_Geno_vcf_data)

chr9Ind <- which(as.numeric(NUST_Geno_vcf_data$CHROM) <10)
chr10Ind <- which(as.numeric(NUST_Geno_vcf_data$CHROM) >=10)

NUST_Geno_vcf_data$CHROM[chr9Ind] <- paste("Gm0",as.numeric(NUST_Geno_vcf_data$CHROM)[chr9Ind],sep="")
NUST_Geno_vcf_data$CHROM[chr10Ind] <- paste("Gm",NUST_Geno_vcf_data$CHROM[chr10Ind],sep="")
NUST_Geno_vcf_data$ChrPosV1 <- paste(NUST_Geno_vcf_data$CHROM,NUST_Geno_vcf_data$POS,sep="-")

### Merge NUSTGenoVCF, INFOCols

INFO_Cols_V1 <- merge(INFO_Cols,NUST_Geno_vcf_data,by="ChrPosV1",all.x=TRUE)

###

REF_Match <- unlist(lapply(c(1:nrow(INFO_Cols_V1)),function(x) if(INFO_Cols_V1[x,"REF.x"]== INFO_Cols_V1[x,"REF.y"]){1}else{0}))
ALT_Match <- unlist(lapply(c(1:nrow(INFO_Cols_V1)),function(x) if(INFO_Cols_V1[x,"ALT.x"]== INFO_Cols_V1[x,"ALT.y"]){1}else{0}))

length(REF_Match)
length(ALT_Match)

INFO_Cols_Mod <- INFO_Cols_V1[,c("CHROM.x","POS.x","ID","REF.x","ALT.x")]
colnames(INFO_Cols_Mod) <- c("CHROM","POS","ID","REF","ALT")
QUAL <- rep(".",nrow(INFO_Cols_Mod))
FILTER <- rep("PASS",nrow(INFO_Cols_Mod))
INFO <- rep(".",nrow(INFO_Cols_Mod))
FORMAT <- rep("GT",nrow(INFO_Cols_Mod))

INFO_Cols_Tab <- cbind.data.frame(INFO_Cols_Mod,QUAL,FILTER,INFO,FORMAT)
INFO_Cols_Tab$`CHROM.POS` <- paste(INFO_Cols_Tab$CHROM,INFO_Cols_Tab$POS,sep="-")

####
### Create INFO Cols Tab from

AGP_Pop_Geno_Filt_VCF <- merge(INFO_Cols_Tab,AGP_Pop_Geno_Filt_Mod,by="CHROM.POS")

rmColInd <- grep("\\.[xy]",colnames(AGP_Pop_Geno_Filt_VCF))
if(length(rmColInd) >0){ rmColIndx <- c(1,rmColInd)}else{rmColIndx <- 1}

AGP_Pop_Geno_Filt_VCF$POS <- as.numeric(AGP_Pop_Geno_Filt_VCF$POS)
AGP_Pop_Geno_Filt_VCF_Mod0 <- AGP_Pop_Geno_Filt_VCF[order(AGP_Pop_Geno_Filt_VCF$CHROM,AGP_Pop_Geno_Filt_VCF$POS),-rmColIndx]
print(dim(AGP_Pop_Geno_Filt_VCF_Mod0))


### Scan table for additional ALT alleles and add to ALT Field
InfoCols <- colnames(INFO_Cols_Tab)

scan_addALT <- function(AGP_Pop_Geno_Filt_VCF_Mod,InfoCols,nR){

	GenoVec <- unlist(AGP_Pop_Geno_Filt_VCF_Mod[nR,])
	GenoVecMod <- (GenoVec)
	samples <- setdiff(colnames(AGP_Pop_Geno_Filt_VCF_Mod),InfoCols)
	sampleColIndx2 <- which(colnames(AGP_Pop_Geno_Filt_VCF_Mod) %in% samples)

	markerVals0 <-levels(factor(GenoVec[sampleColIndx2]))
	markerVals0ModInd <- grep("/",markerVals0)

	if(length(markerVals0ModInd)>0){ markerVals0Mod_1 <- markerVals0[-markerVals0ModInd]}else{markerVals0Mod_1 <- markerVals0}

	markerVals0Mod0 <- grep("/",markerVals0,value=TRUE)
	if(length(markerVals0Mod0) >0){
	  markerVals0Mod <- gsub(" ","",markerVals0Mod0)
	  markerVals0Mod_2 <- unlist(lapply(markerVals0Mod,function(x) unlist(strsplit(x,"/"))))
	}else{markerVals0Mod_2 <- NULL}
	markerVals <- levels(factor(c(markerVals0Mod_1,markerVals0Mod_2)))

	markersTrueSet <- c(GenoVec[c("REF","ALT")],"FAIL")
	addtnMarkerStates <- setdiff(markerVals,markersTrueSet)
    if(length(addtnMarkerStates) >0 ){GenoVecMod["ALT"] <- paste(GenoVecMod["ALT"],paste(addtnMarkerStates,collapse=","),sep=",")}

	return(GenoVecMod)
}

 AGP_Pop_Geno_Filt_VCF_Mod <- do.call(rbind.data.frame,lapply(1:nrow(AGP_Pop_Geno_Filt_VCF_Mod0),function(x) scan_addALT(AGP_Pop_Geno_Filt_VCF_Mod0,InfoCols,x)))
 colnames(AGP_Pop_Geno_Filt_VCF_Mod) <- colnames(AGP_Pop_Geno_Filt_VCF_Mod0)

##### Translate allelic code to numeric code
InfoCols <- colnames(INFO_Cols_Tab)
samples <- setdiff(colnames(AGP_Pop_Geno_Filt_VCF_Mod),InfoCols)
sampleColIndx2 <- which(colnames(AGP_Pop_Geno_Filt_VCF_Mod) %in% samples)
RefIndx <- which(colnames(AGP_Pop_Geno_Filt_VCF_Mod) %in% "REF")
AltIndx <- which(colnames(AGP_Pop_Geno_Filt_VCF_Mod) %in% "ALT")

# Define a function to translate genotype codes
translate_genotype_vcf_num <- function(genotype,sampleColIndx2,RefIndx,AltIndx){
  genotypeMod <- genotype
  HomIndx <- grep("/",genotypeMod[sampleColIndx2],invert=TRUE)
  HetIndx <- grep("/",genotypeMod[sampleColIndx2],invert=FALSE)

  ## FAIL
  genotypeMod[sampleColIndx2] <- gsub("FAIL","\\./\\.",genotypeMod[sampleColIndx2])

  ## REF
  genotypeMod[sampleColIndx2] <- gsub(as.character(genotypeMod[RefIndx]),"0",genotypeMod[sampleColIndx2])
  genotypeMod[sampleColIndx2][HomIndx] <- gsub("0","0/0",genotypeMod[sampleColIndx2][HomIndx])


  ## ALT
  if(length(grep(",",genotypeMod[AltIndx]))>0){
     AltState <- unlist(strsplit(as.character(genotypeMod[AltIndx]),",")[[1]])
  }else{AltState <- genotypeMod[AltIndx]}

  AltState <- gsub(" ","",AltState)

  for(nAlt in 1:length(AltState)){
  #genotypeMod[sampleColIndx2] <- gsub(genotypeMod[AltIndx],"1",genotypeMod[sampleColIndx2])
   genotypeMod[sampleColIndx2] <- gsub(AltState[nAlt],as.character(nAlt),genotypeMod[sampleColIndx2])
  }
  for(nAlt in 1:length(AltState)){
   genotypeMod[sampleColIndx2][HomIndx] <- gsub(as.character(nAlt),paste(as.character(nAlt),"/",as.character(nAlt),sep=""),genotypeMod[sampleColIndx2][HomIndx])
  }

  genotypeMod[sampleColIndx2][HetIndx] <- gsub(" ","",genotypeMod[sampleColIndx2][HetIndx])
  genotypeMod
}

### Translate genotypes for all samples into numeric code

vcf_data_num <- do.call(rbind,lapply(c(1:nrow(AGP_Pop_Geno_Filt_VCF_Mod)),function(x)translate_genotype_vcf_num(AGP_Pop_Geno_Filt_VCF_Mod[x,],sampleColIndx2,RefIndx,AltIndx))) # Translate genotypes
vcf_data_num$CHROM <- as.factor(vcf_data_num$CHROM)
vcf_data_num$POS <- as.numeric(vcf_data_num$POS)

vcf_data_num_sort <- vcf_data_num[order(vcf_data_num$CHROM,vcf_data_num$POS),]
print(dim(vcf_data_num_sort))

print(table(apply(vcf_data_num_sort[,sampleColIndx2],2,as.character)))

### Write complete VCF
#### AGP Genotyping Plate


vcf_data <- vcf_data_num_sort
colCHR <- grep("CHR",colnames(vcf_data))
colFmt <- grep("FORMAT",colnames(vcf_data))
AGP_Geno_Data <- vcf_data

###Filter Markers
### Filt 1: Rm markers with NA IDs
Filt1_Ind <- which(is.na(AGP_Geno_Data[,"ID"]))
if(length(Filt1_Ind)>0){AGP_Geno_Data_Out_Filt1 <- AGP_Geno_Data[-Filt1_Ind,]}else{AGP_Geno_Data_Out_Filt1 <- AGP_Geno_Data}

### Filt2 Rm markers with duplicated IDs
Filt2_Ind <- which(duplicated(AGP_Geno_Data_Out_Filt1[,"ID"]))
if(length(Filt2_Ind)>0){AGP_Geno_Data_Out_Filt2 <- AGP_Geno_Data_Out_Filt1[-Filt2_Ind,]}else{AGP_Geno_Data_Out_Filt2 <- AGP_Geno_Data_Out_Filt1}
AGP_Geno_Data_Out_Filt2$CHROM <- gsub("Gm[0]*","",AGP_Geno_Data_Out_Filt2$CHROM)

scaffoldInd <- grep("scaffold",AGP_Geno_Data_Out_Filt2$CHROM)
if(length(scaffoldInd)>0){AGP_Geno_Data_Out_Filt3 <- AGP_Geno_Data_Out_Filt2[-scaffoldInd,]}else{AGP_Geno_Data_Out_Filt3 <- AGP_Geno_Data_Out_Filt2}

 print(length(samples))

 print(dim(AGP_Geno_Data_Out_Filt3))


 ###

 return(AGP_Geno_Data_Out_Filt3)
}


#######

scan_addALT <- function(AGP_Pop_Geno_Filt_VCF_Mod,InfoCols,nR){

	GenoVec <- unlist(AGP_Pop_Geno_Filt_VCF_Mod[nR,])
	GenoVecMod <- (GenoVec)
	samples <- setdiff(colnames(AGP_Pop_Geno_Filt_VCF_Mod),InfoCols)
	sampleColIndx2 <- which(colnames(AGP_Pop_Geno_Filt_VCF_Mod) %in% samples)

	markerVals0 <-levels(factor(GenoVec[sampleColIndx2]))
	markerVals0ModInd <- grep("/",markerVals0)

	if(length(markerVals0ModInd)>0){ markerVals0Mod_1 <- markerVals0[-markerVals0ModInd]}else{markerVals0Mod_1 <- markerVals0}

	markerVals0Mod0 <- grep("/",markerVals0,value=TRUE)
	if(length(markerVals0Mod0) >0){
	  markerVals0Mod <- gsub(" ","",markerVals0Mod0)
	  markerVals0Mod_2 <- unlist(lapply(markerVals0Mod,function(x) unlist(strsplit(x,"/"))))
	}else{markerVals0Mod_2 <- NULL}
	markerVals <- levels(factor(c(markerVals0Mod_1,markerVals0Mod_2)))

	markersTrueSet <- c(GenoVec[c("REF","ALT")],"FAIL")
	addtnMarkerStates <- setdiff(markerVals,markersTrueSet)
    if(length(addtnMarkerStates) >0 ){GenoVecMod["ALT"] <- paste(GenoVecMod["ALT"],paste(addtnMarkerStates,collapse=","),sep=",")}

	return(GenoVecMod)
}


translate_genotype_vcf_num <- function(genotype,sampleColIndx2,RefIndx,AltIndx){
  genotypeMod <- genotype
  HomIndx <- grep("/",genotypeMod[sampleColIndx2],invert=TRUE)
  HetIndx <- grep("/",genotypeMod[sampleColIndx2],invert=FALSE)

  ## FAIL
  genotypeMod[sampleColIndx2] <- gsub("FAIL","\\./\\.",genotypeMod[sampleColIndx2])

  ## REF
  genotypeMod[sampleColIndx2] <- gsub(as.character(genotypeMod[RefIndx]),"0",genotypeMod[sampleColIndx2])
  genotypeMod[sampleColIndx2][HomIndx] <- gsub("0","0/0",genotypeMod[sampleColIndx2][HomIndx])


  ## ALT
  if(length(grep(",",genotypeMod[AltIndx]))>0){
     AltState <- unlist(strsplit(as.character(genotypeMod[AltIndx]),",")[[1]])
  }else{AltState <- genotypeMod[AltIndx]}

  AltState <- gsub(" ","",AltState)

  for(nAlt in 1:length(AltState)){
  #genotypeMod[sampleColIndx2] <- gsub(genotypeMod[AltIndx],"1",genotypeMod[sampleColIndx2])
   genotypeMod[sampleColIndx2] <- gsub(AltState[nAlt],as.character(nAlt),genotypeMod[sampleColIndx2])
  }
  for(nAlt in 1:length(AltState)){
   genotypeMod[sampleColIndx2][HomIndx] <- gsub(as.character(nAlt),paste(as.character(nAlt),"/",as.character(nAlt),sep=""),genotypeMod[sampleColIndx2][HomIndx])
  }

  genotypeMod[sampleColIndx2][HetIndx] <- gsub(" ","",genotypeMod[sampleColIndx2][HetIndx])
  genotypeMod
}


#######

getCombinedTab_V0 <- function(outputListME,TraitME,IDColsME,IDColME,fitEnvCov){

 if(fitEnvCov==FALSE){
  Models <- c("MM","MDs","MDe")
 }else if(fitEnvCov == TRUE){
  Models <- c("EMM","EMDs","EMDe")
 }
  nTraits <- length(TraitME)
  Pred_Out <- outputListME$preds

  for(nTrt in 1:nTraits){

	traitME <- TraitME[nTrt]
    Fit_Pred <- Pred_Out[[nTrt]]
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
	print(colnames(Fit_Out_Tab))
	print(selOutCols)
	print(setdiff(selOutCols,colnames(Fit_Out_Tab)))
	Fit_Out_Tab_Sel <- Fit_Out_Tab[,selOutCols]


	return(Fit_Out_Tab_Sel)

  }


maskTargetSet <- function(DT_List,maskVar,maskVarLev,maskProp,LocME,YrME){

  nTrt <- 1
  DT_1_Filt_List <- DT_List
  DT_1 <- DT_1_Filt_List[[nTrt]]
  genoDat <- genoDat_List[[nTrt]]
  trait <- traits[nTrt]
  DT_Msk_List <- list()

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


  nanInd <-   which(is.nan(DT_2[,trait]))
  if(length(nanInd)>0){DT_2 <- DT_2[-nanInd,]}

  DT_2 <- droplevels(DT_2)

  maskFactCol <- which(colnames(DT_2) %in% maskFact)
  maskFactInd <- which(DT_2[,maskFactCol] %in% maskFactLev)
  nMaskFact <- length(maskFactInd)
  maskInd <- sample(c(1:nMaskFact),nMaskFact*maskProp,replace=FALSE)

  DT_2[maskInd,Trait] <- NA

  DT_Msk_List[[nTrt]] <- DT_2

  return(DT_Msk_List)
}

