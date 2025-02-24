##### SoygenGS App Functions
# renv::init()
# renv::use_python()


 # if(!require("dplyr", quietly = TRUE)){
    # install.packages("dplyr")
 # }
  # library(dplyr) 

 # if(!require("rrBLUP", quietly = TRUE)){
     # install.packages("rrBLUP")
 # }
  # library(rrBLUP)

# if(!require("vcfR", quietly = TRUE)){
    # install.packages("vcfR")
 # }
  # library(vcfR)

 
# if(!require("NAM", quietly = TRUE)){
     # install.packages("NAM")
 # }
  # library(NAM)

# if(!require("bWGR", quietly = TRUE)){
     # install.packages("bWGR")
 # }
 # library(bWGR)

 # if(!require("STPGA", quietly = TRUE)){
     # install.packages("STPGA")
 # }
  # library(STPGA)
# if(!require("BGLR", quietly = TRUE)){
     # install.packages("BGLR")
 # }
  # library(BGLR)

 
# if(!require("reshape2", quietly = TRUE)){
     # install.packages("reshape2")
 # }
  # library(reshape2)
  
 
 
# if(!require("gplots", quietly = TRUE)){
     # install.packages("gplots")
 # } 
  
  # library(gplots)

  
# #if (!require("BiocManager", quietly = TRUE))
 # # install.packages("BiocManager")
 # # BiocManager::install(version = "3.14")
 
 # # options(repos = BiocManager::repositories())
  # #suppressPackageStartupMessages(library(BiocManager))
   
 
  # options(rsconnect.http.trace = TRUE)
  
   # if(!require("devtools", quietly = TRUE)){
     # install.packages("devtools")
   # }
   # library(devtools)
   
   
   # if(!require("sommer", quietly = TRUE)){
     # library(devtools); install_github('covaruber/sommer')
   # }
   # suppressPackageStartupMessages(library(sommer))
   
   # if(!require("ggplot2", quietly = TRUE)){
     # install.packages("ggplot2")
   # }
   # library(ggplot2)
   
   
   # if(!require("reticulate", quietly = TRUE)){
     # install.packages("reticulate")
   # }
   # library(reticulate)

 
   
 # #   if(!require("rTASSEL", quietly = TRUE)){
# #	 devtools::install_bitbucket(
# #		repo = "bucklerlab/rTASSEL",
# #		host = "bitbucket.org",
# #		ref = "master",
# #		build_vignettes = FALSE,
# #		INSTALL_opts = "--no-multiarch"
# #	 ) 
# # } 
  
# #dyn.load('/usr/lib/jvm/java-17-openjdk-amd64/lib/server/libjvm.so')

 # if(!require("rJava", quietly = TRUE)){  
 
   # install.packages("rJava")
  # }


 # options(java.parameters = "-Xmx30g")


# if(file.exists("/usr/lib/jvm/java-17-openjdk-amd64/lib/server/libjvm.so")){
   
   # dyn.load('/usr/lib/jvm/java-17-openjdk-amd64/lib/server/libjvm.so')
# }




 # library(rJava)
# # if(!require("rTASSEL", quietly = TRUE)){
      # # devtools::install_bitbucket(
		# # repo = "bucklerlab/rTASSEL",
	 	# # ref = "master",
	 	# # build_vignettes = FALSE
      # # ) 
   # # } 
   
 # if(!require("rTASSEL", quietly = TRUE)){  
   # devtools::install_github("maize-genetics/rTASSEL")
  # }
 # library(rTASSEL)
   
  
# if(!require("qtl", quietly = TRUE)){
    # install.packages("qtl")
    
 # }
 # library(qtl)
 # if(!require("PopVar", quietly = TRUE)){
  # install.packages("PopVar")
  
 # }
 # library(PopVar)
 # if(!require("tibble", quietly = TRUE)){
  # install.packages("tibble")
  
 # }
  
 # library(tibble)
 
 # if(!require("EnvRtype", quietly = TRUE)){
  # devtools::install_github("allogamous/EnvRtype")
 # }
 # library(EnvRtype) 
 
# library(renv)
# # Automatically load all installed packages
# installed_packages <- renv::dependencies()$Package
# lapply(installed_packages, library, character.only = TRUE)
 
  
###########################################################################################  
  
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


######

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

	# markers_in_GenoTable <- colnames(NUST_Merged_Table_noDup[BrProg_factor_filt_indices,init:final]) 
	# markerFilter <- which(markers_in_GenoTable %in% markers_with_map)

	# length(markers_in_GenoTable[markerFilter])
	# length(which(markers_with_map %in% markers_in_GenoTable))
	# length(which(is.na(markers_with_map[which(markers_with_map %in% markers_in_GenoTable)])))

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




getMergedData <- function(gt2d,Pheno,testIDs){

   	
	Non_ZeroIndices <- NULL
	Genotypes_VCF <-   gt2d[,-c(1:5,Non_ZeroIndices)]
		
	Genotypes_VCF_ID <- colnames(Genotypes_VCF)
    markerID <- as.vector(unlist(gt2d[,1]))
  
### Remove special characters from strain ID in meta table and genotype table

	Genotypes_VCF_ID <- gsub("-","",Genotypes_VCF_ID) 
	Genotypes_VCF_ID1 <- gsub("\\.","",Genotypes_VCF_ID) 
	Genotypes_VCF_ID2 <- gsub("\\(","",Genotypes_VCF_ID1) 
	Genotypes_VCF_ID3 <- gsub("\\)","",Genotypes_VCF_ID2) 
	Genotypes_VCF_ID4 <- gsub("\\_","",Genotypes_VCF_ID3)


    if(!is.null(testIDs)){ TestIDs <- gsub("[_-]","",testIDs)
	}else if(is.null(testIDs)){TestIDs <- testIDs}
	
### Checks 1


	length(Genotypes_VCF_ID4)
	length(unique(Genotypes_VCF_ID4))

### Remove duplicated IDs from genotype table
	Genotypes_Table_Mod <- cbind(Genotypes_VCF_ID4,t(Genotypes_VCF)) 
	colnames(Genotypes_Table_Mod) <- c("Strain",markerID)

### Numeric coded genotype table from merged data table 

	Genotypes_Table_Mod_Merged <- Genotypes_Table_Mod
	Genotypes_Table_Mod_Num1 <- apply(Genotypes_Table_Mod_Merged[,-1],2,function(x) gsub("BB","-1",x)) 
	Genotypes_Table_Mod_Num2 <- apply(Genotypes_Table_Mod_Num1,2,function(x) gsub("AB","0",x)) 
	Genotypes_Table_Mod_Num3 <- apply(Genotypes_Table_Mod_Num2,2,function(x) gsub("AA","1",x)) 
	
	
	if(length(which(rownames(Genotypes_Table_Mod_Num3) %in% "SNPID")) >=1){
	Genotypes_Table_Mod_Num3A <- Genotypes_Table_Mod_Num3[-which(rownames(Genotypes_Table_Mod_Num3) %in% "SNPID"),]
    }
	if(length(which(rownames(Genotypes_Table_Mod_Num3) %in% "SNPID")) <1){
	Genotypes_Table_Mod_Num3A <- Genotypes_Table_Mod_Num3
    }
	
	Genotypes_Table_Mod_Num <- apply(Genotypes_Table_Mod_Num3A,2,function(x) as.numeric(x)+1)


	Genotypes_Table_Mod_Num_Comb <- cbind(Genotypes_Table_Mod_Merged[,1],Genotypes_Table_Mod_Num)
	colnames(Genotypes_Table_Mod_Num_Comb)[1] <- "Strain" 

	
	
	Genotype_Table_Num <- Genotypes_Table_Mod_Num_Comb
	dim(Genotype_Table_Num)
	colnames(Genotype_Table_Num)
	
# Filter Genotype Table and remove duplicate IDs
	
	dupIndices_Num <- which(duplicated(as.character(Genotype_Table_Num[,1])))
	dupStrain_Num <- as.character(Genotype_Table_Num[dupIndices_Num,1])
	dupIndices_in_Table_Num <- which(as.character(Genotype_Table_Num[,1]) %in% dupStrain_Num)

   if(length(dupIndices_Num) >=1) {
	Genotype_Table_Num_noDup <- Genotype_Table_Num[-dupIndices_Num,]
	rownames(Genotype_Table_Num_noDup) <- Genotype_Table_Num_noDup[,1]
	Genotype_Table_Num_Filt1 <- Genotype_Table_Num_noDup
	
	
	Genotypes_Table_Mod_Num_Filt0 <- apply(as.matrix(Genotype_Table_Num_Filt1[,-1]),2,as.numeric)
	rownames(Genotypes_Table_Mod_Num_Filt0) <- rownames(Genotype_Table_Num_Filt1)
	dim(Genotypes_Table_Mod_Num_Filt0)
	Genotypes_Table_Mod_Num_Filt <- cbind(rownames(Genotypes_Table_Mod_Num_Filt0),Genotypes_Table_Mod_Num_Filt0)
	
   }
   if(length(dupIndices_Num) <1) {
	Genotype_Table_Num_noDup <- Genotype_Table_Num
	rownames(Genotype_Table_Num_noDup) <- Genotype_Table_Num_noDup[,1]
	Genotype_Table_Num_Filt1 <- Genotype_Table_Num_noDup
		
	Genotypes_Table_Mod_Num_Filt0 <- apply(as.matrix(Genotype_Table_Num_Filt1[,-1]),2,as.numeric)
	rownames(Genotypes_Table_Mod_Num_Filt0) <- rownames(Genotype_Table_Num_Filt1)
	dim(Genotypes_Table_Mod_Num_Filt0)
	Genotypes_Table_Mod_Num_Filt<- cbind(rownames(Genotypes_Table_Mod_Num_Filt0),Genotypes_Table_Mod_Num_Filt0)
   }

### Filtered Genotype Table  
### Separate train and test sets
	
	Test_Genotypes_Table_Mod_Num_Filt <- NULL
#### 
    if(!is.null(TestIDs)){
	 
	  testIndices <- which(as.character(Genotypes_Table_Mod_Num_Filt[,1]) %in% TestIDs)
	  StrainIDs <- as.character(Genotypes_Table_Mod_Num_Filt[,1])
	  TrainIDs <- setdiff(StrainIDs,TestIDs)
	  trainIndices <- which(as.character(Genotypes_Table_Mod_Num_Filt[,1]) %in% TrainIDs)
	  Test_Genotypes_Table_Mod_Num_Filt <- Genotypes_Table_Mod_Num_Filt[testIndices,]
          Train_Genotypes_Table_Mod_Num_Filt <- Genotypes_Table_Mod_Num_Filt[trainIndices,]
    }else if(is.null(TestIDs)){ 
     
	  StrainIDs <- as.character(Genotypes_Table_Mod_Num_Filt[,1])
	  TrainIDs <- StrainIDs
	  trainIndices <- which(as.character(Genotypes_Table_Mod_Num_Filt[,1]) %in% TrainIDs)
	  Train_Genotypes_Table_Mod_Num_Filt <- Genotypes_Table_Mod_Num_Filt[trainIndices,]
    }
  
######### Data Prep for GP model training

############# Process IDs
### PhenoIDs  

	StrainID_List <- strsplit(as.character(Pheno[,1]),"[-_.()]")
	StrainID <- unlist(lapply(StrainID_List,function(x) paste(x,collapse="")))

## Remove MG from PhenoIDs
	if(length(grep("MG",as.character(StrainID))) >1){
		  
		  StrainIDMod <- gsub("MG.*","",as.character(StrainID))

	}else if(length(grep("MG",as.character(StrainID))) <1){
		StrainIDMod <- as.character(StrainID)
  	}	

   
        Pheno1 <- cbind(Pheno,StrainIDMod)
	Pheno1[,1] <- StrainID
	colnames(Pheno1)[ncol(Pheno1)] <- "StrainID"
				
### GenoIDs  
	
	# length(which(StrainIDMod %in% trainStrainID))
  
    # if(length(grep("MG",as.character(StrainID))) >1 & Train_Genotypes_Table_Mod_Num_Filt)[,"MG"] )  {
		  
		  # trainStrainIDMG <- paste(trainStrainID,"MG",
		  # length(which(StrainIDMod %in% rownames(Train_Genotypes_Table_Mod_Num_Filt)))
	# }

    # if(length(grep("MG",as.character(StrainID))) <1){
			# StrainIDMod <- as.character(StrainID)
    # }	
  
  
  trainStrainID <- gsub("[-_.()]","",(Train_Genotypes_Table_Mod_Num_Filt)[,1])
  commonStrainID <- Pheno1[(which(as.character(Pheno1[,"StrainID"]) %in% trainStrainID)),"StrainID"]
  Diff_StrainID <- setdiff(trainStrainID,commonStrainID)

  GenoTable_Filtered <- cbind(trainStrainID,Train_Genotypes_Table_Mod_Num_Filt)
  colnames(GenoTable_Filtered)[1] <- "StrainID" 
            
  PhenoTable_Filtered <- Pheno1[which(as.character(Pheno1[,"StrainID"]) %in% trainStrainID),]
	
  
### merge geno and pheno tables  
 # Data_Table1_Num <- merge(GenoTable_Filtered,PhenoTable_Filtered,by="StrainID")
  
  Data_Table1_Num <- merge(GenoTable_Filtered,PhenoTable_Filtered,by="StrainID",all=TRUE)
 
    
 
   if( length(grep("MG",colnames(Data_Table1_Num)))>0 ){ 
      
	   positiveIndices <- which(as.numeric(as.character(Data_Table1_Num[,"MG"])) >0) 
	   zeroIndices <-  which(as.numeric(as.character(Data_Table1_Num[,"MG"])) == 0) 
	   negativeIndices <- which(as.numeric(as.character(Data_Table1_Num[,"MG"])) < 0)
	   StrainIDModPhenoMG <- rep(0,nrow(Data_Table1_Num))
	   
       StrainIDModPhenoMG[positiveIndices] <- paste(as.character(Data_Table1_Num[positiveIndices,"StrainID"]),"MG",as.roman(as.numeric(as.character(Data_Table1_Num[positiveIndices,"MG"]))),sep="")
	   StrainIDModPhenoMG[zeroIndices] <- paste(as.character(Data_Table1_Num[zeroIndices,"StrainID"]),"MG",as.numeric(as.character(Data_Table1_Num[zeroIndices,"MG"])),sep="")
	   StrainIDModPhenoMG[negativeIndices] <- paste(as.character(Data_Table1_Num[negativeIndices,"StrainID"]),"MG00",sep="")
  
       Data_Table1_Num_Mod <- cbind(Data_Table1_Num,StrainIDModPhenoMG)
	   colnames(Data_Table1_Num_Mod)[ncol(Data_Table1_Num_Mod)] <- "StrainIDModPheno"
   
	}
    if(length(grep("MG",colnames(Data_Table1_Num)))==0 ){ 
 
        Data_Table1_Num_Mod <- cbind(Data_Table1_Num,Data_Table1_Num[,"StrainID"])
		colnames(Data_Table1_Num_Mod)[ncol(Data_Table1_Num_Mod)] <- "StrainIDModPheno"
	}
		 
 
### Dealing with duplicated data 
 
  # dupId_Data <- Data_Table1_Num[which(duplicated(Data_Table1_Num[,"StrainID"])),"StrainID"]
  # dupId_Indices <- lapply(dupId_Data,function(x) which(as.character(Data_Table1_Num[,"StrainID"]) %in% as.character(x)))
  # dupId_Indices_Len <- lapply(dupId_Indices,length)
  
  # table(unlist(dupId_Indices_Len)) 
  # dupData_List <- lapply(dupId_Indices,function(x) Data_Table1_Num[as.vector(x),])
  # dupData_NA_List <- lapply(dupData_List,function(x) apply(x,1,function(y) length(which(is.na(y)))))
  
  # StrainID_Table <- cbind(StrainIDMod,StrainID)  
  # colnames(StrainID_Table) <- c("StrainID","StrainIDPheno")
  # Data_Table1_Num_Mod <- merge(StrainID_Table,Data_Table1_Num,by="StrainID",all.y=TRUE)


  dupIDIndices <- which(duplicated(Data_Table1_Num_Mod[,"StrainIDModPheno"]))

    
 if(length(dupIDIndices) >0){
 
    Data_Table_Num <- Data_Table1_Num_Mod[-dupIDIndices,]
	rownames(Data_Table_Num) <- Data_Table_Num[,"StrainIDModPheno"]
  }
  
  if(length(dupIDIndices) ==0){
     Data_Table_Num <- Data_Table1_Num_Mod
	 rownames(Data_Table_Num) <- Data_Table_Num[,"StrainIDModPheno"]
  }
     
## Check difference of genotypic scores 
   
  
####### Filtered Table in Numeric Format 
####################################################################################
   #Data_Table_Num <- Data_Table1_Num
   #rownames(Data_Table_Num) <- Data_Table1_Num[,"StrainID"]
  
   Train_Data_Table_Num_Filt <- Data_Table_Num 
  	
   return(list(Train_Data_Table_Num_Filt,Test_Genotypes_Table_Mod_Num_Filt))

}


getProcessedData <- function(Data_Table_Num_List,trait){

     TrainData_Table_Num <- Data_Table_Num_List[[1]]
	 TestData_Table_Num_Filt<- Data_Table_Num_List[[2]]
	 
	 ### Remove lines with 'NA' for trait values
	 if(length(trait)==1){
	  NAIndices <-  which(is.na(TrainData_Table_Num[,trait]))
	  if(length(NAIndices)>1){
	  
        Train_Data_Table_Num_Filt <- TrainData_Table_Num[-NAIndices,]
	  }
	  if(length(NAIndices)<1){
	  
        Train_Data_Table_Num_Filt <- TrainData_Table_Num
	  }
	 }
	 
	 if(length(trait)>1){
	   NAIndices <- c(unlist(apply(TrainData_Table_Num[,trait],2,function(x)which(is.na(x)))))
	   if(length(NAIndices)>1){
	  
        Train_Data_Table_Num_Filt <- TrainData_Table_Num[-NAIndices,]
	   }
	   if(length(NAIndices)<1){
	  
        Train_Data_Table_Num_Filt <- TrainData_Table_Num
	   }
      
	  }
     return(list(Train_Data_Table_Num_Filt,TestData_Table_Num_Filt)) 
	 
} 



getTasObj <- function(infileVCF){
    tasGeno <- rTASSEL::readGenotypeTableFromPath(
		path = infileVCF,
		sortPositions=TRUE
	)
    return(tasGeno)
}



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




# getImputedData <- function(FiltGeno,l,k,impMethod){ 

  # if(impMethod=="LDKNNI"){
   # tasGenoImp <- imputeLDKNNi(FiltGeno, highLDSSites = l, knnTaxa = k, maxDistance = 1e+07)
  # } 
  # if(impMethod=="Numeric"){
   # tasGenoImp <- imputeNumeric(FiltGeno,byMean = TRUE,nearestNeighbors = 5, distance = c("Euclidean", "Manhattan", "Cosine")[1])
  # }
  
  # return(tasGenoImp)
# } 


getImputedData_LDKNNI <- function(FiltGeno,l,k){ 

   tasGenoImp <- rTASSEL::imputeLDKNNi(FiltGeno, highLDSSites = l, knnTaxa = k, maxDistance = 1e+07)
  
  return(tasGenoImp)
} 

getImputedData_Num <- function(FiltGeno,nN,Dist){
  
  
   tasGenoImp <- rTASSEL::imputeNumeric(FiltGeno,byMean = TRUE,nearestNeighbors = nN, distance = Dist)
  
  return(tasGenoImp)
}


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


### Step 2 : Get predicted genetic values usin kin.blup 

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



### Step 2 : Get predicted genetic values usin kin.blup 

 getRankedPredictedValues_Old <- function(Data_Table_Num_Filt_List,nTraits,trait,GPModel,fixedX=NULL,optTS=NULL){ 
    
     if(!is.null(fixedX) & fixedX !="NULL"){ 
	   GPModel <- "rrBLUP (rrBLUP)"
	   Fixed.X <- fixedX[[1]]
	   Test.X <- fixedX[[2]]
	   
	 }
	   
	 TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	 TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
	 
	
	 TestGenoTable <- apply(TestData_Table_Num_Filt[,-1],2,as.numeric)
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
		 trainSetID <- as.character(as.vector(optTS))
	     trainIndices <- which(strainID %in% trainSetID)
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
	 
	

 ### Remove lines with missing pheno	 
	  
	if(anyNA(trainPheno0)){
	   ph_NA_Indices <- which(is.na(trainPheno0))
	   trainPheno <- trainPheno0[-ph_NA_Indices] 
	   trainGeno_Imp <- trainGeno_Imp0[-ph_NA_Indices,]
	   Geno <- as.character(TrainData_Table_Num_Filt[-ph_NA_Indices,1])
	   if(!is.null(fixedX)  & fixedX !="NULL" ){
	     Fixed.X.Mod <- Fixed.X[-ph_NA_Indices,]
	   }
	   
	 }
	 if(!anyNA(trainPheno0)){
	  
	   trainPheno <- trainPheno0
	   trainGeno_Imp <- trainGeno_Imp0
	   Geno <- as.character(TrainData_Table_Num_Filt[,1])
	   if(!is.null(fixedX)  & fixedX !="NULL"){
	     Fixed.X.Mod <- Fixed.X
	   }
	  
	 } 
	 
####### Set the same set of markers for both train and test sets with the same order 
     ##Check with sorted markers and check with test set 
	  trainssIDs <- colnames(trainGeno_Imp)
	  testssIDs <- colnames(testGeno_Imp) 
	  testGeno_Imp_Ord <- testGeno_Imp[, match(trainssIDs,testssIDs)]
	  testGeno_Imp <- testGeno_Imp_Ord
	 
	  
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
	  
	
	 if(!is.null(fixedX) & fixedX !="NULL"){
	   pred <- mixed.solve(trainPheno,Z=trainGeno_Imp,X=Fixed.X.Mod,SE=FALSE,return.Hinv =FALSE) 
	   Mean <- Fixed.X.Mod %*% pred$beta
	 }
	 if(is.null(fixedX) | fixedX =="NULL"){
	     pred <- mixed.solve(trainPheno,Z=trainGeno_Imp,SE=FALSE,return.Hinv =FALSE) 
	     Mean <- as.numeric(pred$beta)
	 }
 	  Effects <- pred$u
   
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
 	  PredictedValues <- Mean + (trainGeno_Imp %*% Effects)
 	  SortedPredictedValues <- sort.int(PredictedValues,decreasing=TRUE,index.return=TRUE)
 	  
   
	   if(!is.null(fixedX) & fixedX !="NULL"){
	        Mean.tst <- Test.X %*% pred$beta
		  }
		  if(is.null(fixedX) | fixedX =="NULL" ){
		    Mean.tst <- as.numeric(pred$beta)
		  }
     
    
		 		  
	   Test_PredictedValues <-  Mean.tst + (testGeno_Imp %*% Effects)
	 
	   Test_SortedPredictedValues <- sort.int(Test_PredictedValues,decreasing=TRUE,index.return=TRUE) 
	   Test_StrainID <- rownames(TestData_Table_Num_Filt)
    
	 }
	 
	 
	 if(GPModel == "rrBLUP (bWGR)"){
	 
		 pred <- emRR(trainPheno,trainGeno_Imp) 
		 Mean <- as.numeric(pred$mu)
		 Effects <- pred$b
		 
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
	# if(length(testNAIndices)==0){ }
	 
	#print(anyNA(trainGeno_Imp2))
	#print(anyNA(trainPheno))
	#write.table(trainPheno,"trainPhenoCRep.csv")
	#write.table(trainGeno_Imp2,"trainGenoTab.txt",sep="\t",quote=FALSE,row.names=FALSE)
   	trainGeno_Imp2 <- apply(trainGeno_Imp,2,function(x) x+1) 
	cleanData <- cleanREPV2(trainPheno,trainGeno_Imp2)
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
	# if(length(testNAIndices)==0){ }
	  
	U <- apply(testGeno_Imp,1,function(x) getU(M.Pdt,x)) 
	
	
	Test_SortedIndices <- Test_SortedPredictedValues[[2]]
	Sorted_Test_StrainID <- Test_StrainID[Test_SortedIndices]
	U_Sorted <- U[(Test_SortedIndices)]
	
	outputDF <- cbind.data.frame(Sorted_Test_StrainID ,round(Test_SortedPredictedValues[[1]],digits=2),round(U_Sorted,digits=5))
	colnames(outputDF) <- c("LineID",paste("Predicted Value for ",trait,sep=""),"Upper Bound of Reliability")
	
	
   return(outputDF)
	
 }
 
 
	
### 


 getRankedPredictedValues <- function(Data_Table_Num_Filt_List,nTraits,trait,GPModel,fixedX=NULL,fixedData=NULL,optTS=NULL){ 
    
     if(!is.null(fixedX) & fixedX != "NULL" & fixedData !="NULL"){ 
	   GPModel <- "rrBLUP (rrBLUP)"
	   Fixed.X <- fixedData[[1]]
	   Test.X <- fixedData[[2]]
	 }
	   
	 TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	 TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
	 
	
	 TestGenoTable <- apply(TestData_Table_Num_Filt[,-1],2,as.numeric)
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
	 
####### Set the same set of markers for both train and test sets with the same order 
     ##Check with sorted markers and check with test set 
	  trainssIDs <- colnames(trainGeno_Imp)
	  testssIDs <- colnames(testGeno_Imp) 
	  testGeno_Imp_Ord <- testGeno_Imp[, match(trainssIDs,testssIDs)]
	  testGeno_Imp <- testGeno_Imp_Ord
	 
	  
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
	
 	  PredictedValues <- Mean + (trainGeno_Imp %*% Effects)
 	  SortedPredictedValues <- sort.int(PredictedValues,decreasing=TRUE,index.return=TRUE)
 	  
   
	  if(!is.null(fixedX) & fixedX !="NULL" & fixedData != "NULL" ){
	        Mean.tst <- Test.X %*% pred$beta
	  }
	  if(is.null(fixedX) | fixedX =="NULL" | fixedData == "NULL" ){
		    Mean.tst <- as.numeric(pred$beta)
	  }
     
    
		 		  
	   Test_PredictedValues <-  Mean.tst + (testGeno_Imp %*% Effects)
	 
	   Test_SortedPredictedValues <- sort.int(Test_PredictedValues,decreasing=TRUE,index.return=TRUE) 
	   Test_StrainID <- rownames(TestData_Table_Num_Filt)
    
	 }
	 
	 
	 if(GPModel == "rrBLUP (bWGR)"){
	 
		 pred <- emRR(trainPheno,trainGeno_Imp) 
		 Mean <- as.numeric(pred$mu)
		 Effects <- pred$b
		 
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
		 
		 	 
		 PredictedValues <- Mean + (trainGeno_Imp %*% Effects)
		 SortedPredictedValues <- sort.int(PredictedValues,decreasing=TRUE,index.return=TRUE)
	     Test_PredictedValues <-  Mean + (testGeno_Imp %*% Effects)
	     Test_SortedPredictedValues <- sort.int(Test_PredictedValues,decreasing=TRUE,index.return=TRUE) 
	     Test_StrainID <- rownames(TestData_Table_Num_Filt)
    }
	  
### 
	# pred.kb <- kin.blup(as.data.frame(Data),Geno,Pheno,GAUSS=FALSE,K=A,covariate=NULL,PEV=TRUE,n.core=1,theta.seq=NULL)
		 
### Upper bound of Reliability
   
	trainGeno_Imp2 <- apply(trainGeno_Imp,2,function(x) as.numeric(x)+1)
	
	rownames(trainGeno_Imp2) <- rownames(trainGeno_Imp)
	print(dim(trainGeno_Imp2))
	print(length(trainPheno))
	cleanData <- cleanREPV2(trainPheno,trainGeno_Imp)
	
	#cleanData <- cleanREP(trainPheno,trainGeno_Imp2)
    
	M <-  cleanData[[2]]
    M.Pdt <- t(M)%*% solve(M %*% t(M) + diag(nrow(M))) %*% M
      
    getU <- function(M.Pdt,v){
	   v.hat <- M.Pdt %*% v
	   U <- (t(v.hat)%*% v.hat)/ (t(v) %*% v)
	   return(U)
    }
	
		  
	U <- apply(testGeno_Imp,1,function(x) getU(M.Pdt,x)) 
		
	Test_SortedIndices <- Test_SortedPredictedValues[[2]]
	Sorted_Test_StrainID <- Test_StrainID[Test_SortedIndices]
	U_Sorted <- U[(Test_SortedIndices)]
	
	outputDF <- cbind.data.frame(Sorted_Test_StrainID ,round(Test_SortedPredictedValues[[1]],digits=2),round(U_Sorted,digits=5))
	colnames(outputDF) <- c("LineID",paste("Predicted Value for ",trait,sep=""),"Upper Bound of Reliability")
	
	
   return(outputDF)
   
    # if(length(testNAIndices)>0){
       # trainGeno_Imp2 <- apply(trainGeno_Imp[,-testNAIndices],2,function(x) x+1)  
	# }
	# if(length(testNAIndices)==0){ }
	
	#print(anyNA(trainGeno_Imp2))
	#print(anyNA(trainPheno))
	#write.table(trainPheno,"trainPhenoCRep.csv")
	#write.table(trainGeno_Imp2,"trainGenoTab.txt",sep="\t",quote=FALSE,row.names=FALSE)
   	
	# if(length(testNAIndices)>0){
      # U <- apply(testGeno_Imp[,-testNAIndices],1,function(x) getU(M.Pdt,x)) 
	# }
	# if(length(testNAIndices)==0){ }
	
 }
 
	
#### 


getRankedPredictedValuesMT <- function(Data_Table_Num_Filt_List,nTraits,trait,GPModelMT,optTS=NULL){ 
   
	 TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	 TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
	 			 
	 TestGenoTable <- apply(TestData_Table_Num_Filt[,-1],2,as.numeric)
	 
	
	 rownames(TestGenoTable) <- TestData_Table_Num_Filt[,1]
	 
	 set.seed(125)
	 strainID <- as.character(TrainData_Table_Num_Filt[,1])
	 
	 initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	 
	 # trainSetID <- as.character(as.vector( strainID[sample(c(1:length(strainID)),500)]))
	 # optTS <- cbind(trainSetID,c(1:500))
	 
	 trainSetID <- as.character(as.vector(strainID))
	 optTS <- cbind(trainSetID,c(1:length(trainSetID)))
	 
	  initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	  finalCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[length(grep("ss",colnames(TrainData_Table_Num_Filt)))]
	 
	
	 if(is.null(optTS)){
	
		 trainPheno0 <- TrainData_Table_Num_Filt[,trait]
		 geno_012 <- apply(TrainData_Table_Num_Filt[,c(initCol:finalCol)],2,as.numeric)
		 trainGeno <- apply(geno_012,2,function(x) x-1)
		 #Geno <- rownames(TrainData_Table_Num_Filt)
		 Geno <- as.character(TrainData_Table_Num_Filt[,1])
	 }
	 
	 if(!is.null(optTS)){
	     
	     strainID <- as.character(TrainData_Table_Num_Filt[,1])
		 trainSetID <- as.character(as.vector(optTS[,1]))
	     trainIndices <- which(strainID %in% trainSetID)
		 trainPheno0 <- TrainData_Table_Num_Filt[trainIndices,trait]
		 geno_012 <- apply(TrainData_Table_Num_Filt[trainIndices,c(initCol:finalCol)],2,as.numeric)
		 trainGeno <- apply(geno_012,2,function(x) x-1)
		 Geno <- as.character(TrainData_Table_Num_Filt[trainIndices,1])
	 }
	  
## Prepare Data Table for GP 
	 # Data <- cbind.data.frame(Geno,trainPheno)
	 # colnames(Data) <- c("Geno","Pheno")
	 # Geno <- "Geno"
	 # Pheno <- "Pheno"
	

#### Impute trainGeno and train using mixed.solve
	
	# if(anyNA(trainGeno)){	 
		# trainGeno_Imp <- snpQC(trainGeno,impute=TRUE,remove=FALSE)
	# }else{trainGeno_Imp <- trainGeno}
	
	
### Check if the markers are missing in > 99% of the lines 
	 
	  allNAMarkers <- apply(trainGeno,2,function(x) if(length(which(is.na(x)))==nrow(trainGeno)){1}else{0})
     
	  allNAMarkerIndices <- which(allNAMarkers!=0)
	  
	  if(length(allNAMarkerIndices)>=1){
		 trainGeno0 <- trainGeno[,-allNAMarkerIndices]
	  }else{ trainGeno0 <- trainGeno}
	 	 
	 
	 if(anyNA(trainGeno0)){	 
		trainGeno_Imp0 <- snpQC(trainGeno0,impute=TRUE,remove=FALSE)
	 }else{trainGeno_Imp0 <- trainGeno0}
	 
	 
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
	 
 ### Remove lines with missing pheno	 
	  
	if(anyNA(trainPheno0)){
	   ph_NA_Indices_Tab <- apply(trainPheno0,2,function(x) which(is.na(x)))
	   ph_NA_Indices <- unique(c(unlist(ph_NA_Indices_Tab)))
	   trainPheno <- trainPheno0[-ph_NA_Indices,]
	   trainGeno_Imp <- trainGeno_Imp0[-ph_NA_Indices,]
	   Geno <- as.character(TrainData_Table_Num_Filt[-ph_NA_Indices,1])
	}
	if(!anyNA(trainPheno0)){
	  
	   trainPheno <- trainPheno0
	   trainGeno_Imp <- trainGeno_Imp0
	   Geno <- as.character(TrainData_Table_Num_Filt[,1])
	}
	 
####### Set the same set of markers for both train and test sets in the same order 
     ##Check with sorted markers and check with test set 
	  trainssIDs <- colnames(trainGeno_Imp)
	  testssIDs <- colnames(testGeno_Imp) 
	  testGeno_Imp_Ord <- testGeno_Imp[, match(trainssIDs,testssIDs)]
	  testGeno_Imp <- testGeno_Imp_Ord

## Kinship Matrix 
    # rownames(trainGeno_Imp) <- Geno
	 #rownames(trainPheno) <- Geno
	 A <- A.mat(trainGeno_Imp)	
	 colnames(A) <- Geno
	 rownames(A) <- Geno

	 Y <- trainPheno
	 X <- rbind(trainGeno_Imp,testGeno_Imp)
	 #rownames(X) <- c(rownames(trainGeno_Imp),rownames(TestGenoTable))
	 n <- nrow(X)
	 p <- ncol(X)
	 
	 nTst <- nrow(testGeno_Imp)
	 Ytst <- matrix(rep(NA,nTst*ncol(Y)),nrow=nTst,ncol=ncol(Y))
	 colnames(Ytst) <- colnames(Y)
	 yNA <- as.matrix(rbind(Y,Ytst))
	 yNA_DF <- as.data.frame(yNA)
	 
	 #yNA_DF$id <- as.factor(rownames(X))
	 yNA_DF$id <- as.factor(paste(c(Geno,rownames(testGeno_Imp)),"_",yNA_DF[,1],sep=""))
	 
	 if(GPModelMT == "BRR (BGLR)"){ 
				  
			ETA <- list(list(X=X,model="BRR"))
			fm3 <- Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
	 }
	 if(GPModelMT == "RKHS (BGLR)"){ 
	        A.Tot <- A.mat(X)
			ETA <- list(list(K=A.Tot,model="RKHS"))
			fm3 <- Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
	 }
	 if(GPModelMT == "Spike-Slab (BGLR)"){ 
				  
			ETA <- list(list(X=X,model="SpikeSlab"))
			fm3 <- Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
	 }	
	 if(GPModelMT == "GBLUP (SOMMER)"){ 
	         
			A.Tot <- A.mat(X)
			# rownames(A.Tot) <- rownames(X) 
			# colnames(A.Tot) <- rownames(X)
			rownames(A.Tot) <- yNA_DF$id #c(Geno,rownames(TestGenoTable))
			colnames(A.Tot) <- yNA_DF$id #c(Geno,rownames(TestGenoTable))
			fm3 <- mmer(as.formula(paste("cbind(",paste(trait,collapse=","),")~1",sep="")),
            random=~vs(id,Gu=A.Tot),
            rcov=~units,data=yNA_DF,verbose = TRUE)
	 }
	 tst <- c((length(trainSetID)+1):nrow(X))
	 if(GPModelMT != "GBLUP (SOMMER)"){ 
	  Test_PredictedValues <- fm3$ETAHat[tst,]
	 }
	 if(GPModelMT == "GBLUP (SOMMER)"){  
	  UMat <- do.call(cbind,lapply(fm3$U$`u:id`,function(x) x))
	  Test_PredictedValues <- c()
	  for(nT in 1:length(trait)){
	    Test_PredictedValues <- cbind(Test_PredictedValues,(UMat[,nT]+fm3$Beta[nT,3])[tst])
	  }
	 }
	 
	Test_SortedPredictedValues <- sort.int(Test_PredictedValues[,1],decreasing=TRUE,index.return=TRUE) 
	Test_StrainID <- rownames(TestData_Table_Num_Filt)	 
### Upper bound of Reliability

   # trainGeno_Imp2 <- apply(trainGeno_Imp,2,function(x) x+1)  
   
    
    cleanData <- cleanREP(trainPheno[,1],apply(trainGeno_Imp,2,function(x) x))
    M <-  cleanData[[2]]
    M.Pdt <- t(M)%*% solve(M %*% t(M) + diag(nrow(M))) %*% M
      
    getU <- function(M.Pdt,v){
	   v.hat <- M.Pdt %*% v
	   U <- (t(v.hat)%*% v.hat)/ (t(v) %*% v)
	   return(U)
    }
    
	U <- apply(testGeno_Imp,1,function(x) getU(M.Pdt,x)) 

### Sort Output
	
	Test_SortedIndices <- Test_SortedPredictedValues[[2]]
	Sorted_Test_StrainID <- Test_StrainID[Test_SortedIndices]
	U_Sorted <- U[(Test_SortedIndices)]
	
	Test_SortedPredValuesTab <- apply(Test_PredictedValues[Test_SortedIndices,],2,function(x) round(x,digits=2))
	
### Output DF	
	outputDF <- cbind.data.frame(Sorted_Test_StrainID ,Test_SortedPredValuesTab,round(U_Sorted,digits=5))
	colnames(outputDF) <- c("LineID",colnames(Test_SortedPredValuesTab),"Upper Bound of Reliability")
	
	
   return(outputDF)
	
 }


  
	
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
 

# Data_Table_Num_Filt_List <- predictionData()
# trait<- unlist(Trait())
# nTraits <- nTraits()
# k<- k()
# nIter <- nIter()
  
 
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
		geno_Imp <- snpQC(geno,impute=TRUE,remove=FALSE)
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
	 }
	
	 if(GPModelMT == "RKHS (BGLR)"){ 
	        A.Tot <- rrBLUP::A.mat(X)
			ETA <- list(list(K=A.Tot,model="RKHS"))
			fm3 <- BGLR::Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
	 }
	 if(GPModelMT == "Spike-Slab (BGLR)"){ 
				  
			ETA <- list(list(X=X,model="SpikeSlab"))
			fm3 <- BGLR::Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
	 }	
	 if(GPModelMT == "GBLUP (SOMMER)"){ 
	         
			A.Tot <- A.mat(X)
			rownames(A.Tot) <- rownames(X) 
			colnames(A.Tot) <- rownames(X)
			fm3 <- sommer::mmer(as.formula(paste("cbind(",paste(trait,collapse=","),")~1",sep="")),
            random=~vs(id,Gu=A.Tot),
            rcov=~units,
            data=yNA_DF,verbose = TRUE)
	 }
	 
	 tst <- testIndices
	
	 if(GPModelMT != "GBLUP (SOMMER)"){ 
	  Test_PredictedValues <- fm3$ETAHat[tst,]
	 }
	 
	 if(GPModelMT == "GBLUP (SOMMER)"){  
	  UMat <- do.call(cbind,lapply(fm3$U$`u:id`,function(x) x))
	  Test_PredictedValues <- c()
	  for(nT in 1:length(trait)){
	    Test_PredictedValues <- cbind(Test_PredictedValues,(UMat[,nT]+fm3$Beta[nT,3])[tst])
		
	  }
	 }
	 	 
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

 getOptimalTS <- function(Data_Table_Num_Filt_List,trait,nTraits,noCandidates,nTrainToSelect,GAParameters){
  
    TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
	
	# NAIndices <-  which(is.na(TrainData_Table_Num[,trait]))
    # Train_Data_Table_Num_Filt <- NUST_Data_Table_Num[-NAIndices,]
		 
    initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	finalCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[length(grep("ss",colnames(TrainData_Table_Num_Filt)))]
	 
	initColTst <- grep("ss",colnames(TestData_Table_Num_Filt))[1]
	finalColTst <- grep("ss",colnames(TestData_Table_Num_Filt))[length(grep("ss",colnames(TestData_Table_Num_Filt)))]
	 
	
	trainPheno <- TrainData_Table_Num_Filt[,trait]
	trainGeno_012 <- apply(TrainData_Table_Num_Filt[,c(3:(ncol(TrainData_Table_Num_Filt)-nTraits-2))],2,as.numeric)
	testGeno_012 <- apply(TestData_Table_Num_Filt[,c(2:(ncol(TestData_Table_Num_Filt)))],2,as.numeric)
	
		
	# trainGeno_012 <- apply(TrainData_Table_Num_Filt[,c(initCol:finalCol)],2,as.numeric)
	# testGeno_012 <- apply(TestData_Table_Num_Filt[,c(initColTst:finalColTst)],2,as.numeric)
	
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
         Train_STPGA <- GenAlgForSubsetSelection(P=PC99,Candidates,Test,ntoselect=nTrainToSelect,InitPop=GAParameters$InitPop,
         npop=GAParameters$npop, nelite=GAParameters$nelite, mutprob=GAParameters$mutprob, mutintensity = GAParameters$mutintensity,
         niterations=GAParameters$niterations,minitbefstop=GAParameters$minitbefstop, tabu=GAParameters$tabu,
         tabumemsize = GAParameters$tabumemsize,plotiters=GAParameters$plotiters,errorstat=GAParameters$errorstat,
		 lambda=GAParameters$lambda, mc.cores=GAParameters$mc.cores)
		 
	   }
	   if(length(GAParameters$errorstat)>1){
         Train_STPGA <- GenAlgForSubsetSelectionMO(P=PC99,Candidates,Test,ntoselect=nTrainToSelect,InitPop=GAParameters$InitPop,
         npop=GAParameters$npop, nelite=GAParameters$nelite, mutprob=GAParameters$mutprob, mutintensity = GAParameters$mutintensity,
         niterations=GAParameters$niterations,minitbefstop=GAParameters$minitbefstop, tabu=GAParameters$tabu,
         tabumemsize = GAParameters$tabumemsize,plotiters=GAParameters$plotiters,errorstat=GAParameters$errorstat,
		 lambda=GAParameters$lambda, mc.cores=GAParameters$mc.cores)
		 
	   }
	  
    })
	
	Train_STPGA_Ids <- as.character(Train_STPGA$`Solution with rank 1`)
	Train_STPGA_Indices <-  which(Candidates %in% Train_STPGA_Ids)
    return(list(Train_STPGA,Train_STPGA_Ids,Train_STPGA_Indices))
 }
 
 
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
	 
	
	trainGeno_012 <- apply(TrainData_Table_Num_Filt[,c(3:(ncol(TrainData_Table_Num_Filt)-nTraits-2))],2,as.numeric)
	testGeno_012 <- apply(TestData_Table_Num_Filt[,c(2:(ncol(TestData_Table_Num_Filt)))],2,as.numeric)
	
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
  
   
 getTSComparisons <- function(Data_Table_Num_Filt_List,Train_STPGA,Train_Random,trait,nTraits,testIds){ 
  
    TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
	
	trainPheno <- TrainData_Table_Num_Filt[,trait]
	names(trainPheno) <- rownames(TrainData_Table_Num_Filt)
	
	initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	finalCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[length(grep("ss",colnames(TrainData_Table_Num_Filt)))]
	 
	initColTst <- grep("ss",colnames(TestData_Table_Num_Filt))[1]
	finalColTst <- grep("ss",colnames(TestData_Table_Num_Filt))[length(grep("ss",colnames(TestData_Table_Num_Filt)))]
	 
   	
	trainGeno_012 <- apply(TrainData_Table_Num_Filt[,c(3:(ncol(TrainData_Table_Num_Filt)-nTraits-2))],2,as.numeric)
	testGeno_012 <- apply(TestData_Table_Num_Filt[,c(2:(ncol(TestData_Table_Num_Filt)))],2,as.numeric)
	
		
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
	
	
	# if(is.null(optTS)){
	  
	   # train_STPGA_Ind <- c(1:length(Candidates))
	# } 
    # if(!is.null(optTS)){ 
  
       # train_STPGA_Ind <- which(Candidates %in% as.character(as.vector(optTS))) 
	# }
  
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
  
   
   
 getTSComparisonsMT <- function(Data_Table_Num_Filt_List,Train_STPGA,Train_Random,trait,nTraits,testIds){ 
  
    TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
		
	
	initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	finalCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[length(grep("ss",colnames(TrainData_Table_Num_Filt)))]
	 
	initColTst <- grep("ss",colnames(TestData_Table_Num_Filt))[1]
	finalColTst <- grep("ss",colnames(TestData_Table_Num_Filt))[length(grep("ss",colnames(TestData_Table_Num_Filt)))]
	 
	
	trainGeno_012 <- apply(TrainData_Table_Num_Filt[,c(3:(ncol(TrainData_Table_Num_Filt)-nTraits-2))],2,as.numeric)
	testGeno_012 <- apply(TestData_Table_Num_Filt[,c(2:(ncol(TestData_Table_Num_Filt)))],2,as.numeric)
		
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
		 }
		
		 if(GPModelMT == "RKHS (BGLR)"){ 
				A.Tot <- A.mat(X)
				ETA <- list(list(K=A.Tot,model="RKHS"))
				fm3 <- Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
		 }
		 if(GPModelMT == "Spike-Slab (BGLR)"){ 
					  
				ETA <- list(list(X=X,model="SpikeSlab"))
				fm3 <- Multitrait(y=yNA,ETA=ETA,nIter=1000,burnIn=500)
		 }	
		 if(GPModelMT == "GBLUP (SOMMER)"){ 
				 
				A.Tot <- A.mat(X)
				rownames(A.Tot) <- rownames(X) 
				colnames(A.Tot) <- rownames(X)
				fm3 <- mmer(as.formula(paste("cbind(",paste(trait,collapse=","),")~1",sep="")),
				random=~vs(id,Gu=A.Tot),
				rcov=~units,
				data=yNA_DF,verbose = TRUE)
		 }
		 
		tst <- testIndices
	
		if(GPModelMT != "GBLUP (SOMMER)"){ 
		  Test_PredictedValues <- fm3$ETAHat[tst,]
		}
		if(GPModelMT == "GBLUP (SOMMER)"){  
		  UMat <- do.call(cbind,lapply(fm3$U$`u:id`,function(x) x))
		  Test_PredictedValues <- c()
		  for(nT in 1:length(trait)){
		 
		    Test_PredictedValues <- cbind(Test_PredictedValues,(UMat[,nT]+fm3$Beta[nT,3])[tst])
		  }
		}
	 	 
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
 
 
 
  
	
#function to remove repeated genotypes

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



getRankedPredictedValuesME_SOMMER <- function(Data_Table_Num_Filt_List,nTraits,trait,GPModelVarCor,optTS=NULL){ 
   
	 TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	 TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
	 			 
	 TestGenoTable <- apply(TestData_Table_Num_Filt[,-1],2,as.numeric)
	 
	
	 rownames(TestGenoTable) <- TestData_Table_Num_Filt[,1]
	 
	 set.seed(125)
	 strainID <- as.character(TrainData_Table_Num_Filt[,1])
	 
	 initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	 
	 # trainSetID <- as.character(as.vector( strainID[sample(c(1:length(strainID)),500)]))
	 # optTS <- cbind(trainSetID,c(1:500))
	 
	 trainSetID <- as.character(as.vector(strainID))
	 optTS <- cbind(trainSetID,c(1:length(trainSetID)))
	 
	 initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	 finalCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[length(grep("ss",colnames(TrainData_Table_Num_Filt)))]
	 
	
	 if(is.null(optTS)){
	
		 trainPheno0 <- TrainData_Table_Num_Filt[,c(trait,"environ")]
		 geno_012 <- apply(TrainData_Table_Num_Filt[,c(initCol:finalCol)],2,as.numeric)
		 trainGeno <- apply(geno_012,2,function(x) x-1)
		 #Geno <- rownames(TrainData_Table_Num_Filt)
		 Geno <- as.character(TrainData_Table_Num_Filt[,1])
	 }
	 
	 if(!is.null(optTS)){
	     
	     strainID <- as.character(TrainData_Table_Num_Filt[,1])
		 trainSetID <- as.character(as.vector(optTS[,1]))
	     trainIndices <- which(strainID %in% trainSetID)
		 trainPheno0 <- TrainData_Table_Num_Filt[trainIndices,c(trait,"environ")]
		 geno_012 <- apply(TrainData_Table_Num_Filt[trainIndices,c(initCol:finalCol)],2,as.numeric)
		 trainGeno <- apply(geno_012,2,function(x) x-1)
		 Geno <- as.character(TrainData_Table_Num_Filt[trainIndices,1])
	 }
	  
## Prepare Data Table for GP 
	 # Data <- cbind.data.frame(Geno,trainPheno)
	 # colnames(Data) <- c("Geno","Pheno")
	 # Geno <- "Geno"
	 # Pheno <- "Pheno"
	

#### Impute trainGeno and train using mixed.solve
	
	# if(anyNA(trainGeno)){	 
		# trainGeno_Imp <- snpQC(trainGeno,impute=TRUE,remove=FALSE)
	# }else{trainGeno_Imp <- trainGeno}
	
	
### Check if the markers are missing in > 99% of the lines 
	 
	  allNAMarkers <- apply(trainGeno,2,function(x) if(length(which(is.na(x)))==nrow(trainGeno)){1}else{0})
     
	  allNAMarkerIndices <- which(allNAMarkers!=0)
	  
	  if(length(allNAMarkerIndices)>=1){
		 trainGeno0 <- trainGeno[,-allNAMarkerIndices]
	  }else{ trainGeno0 <- trainGeno}
	 	 
	 
	 if(anyNA(trainGeno0)){	 
		trainGeno_Imp0 <- snpQC(trainGeno0,impute=TRUE,remove=FALSE)
	 }else{trainGeno_Imp0 <- trainGeno0}
	 
	 
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
	 
	

 ### Remove lines with missing pheno	 
	  
	if(anyNA(trainPheno0)){
	   ph_NA_Indices_Tab <- apply(trainPheno0,2,function(x) which(is.na(x)))
	   ph_NA_Indices <- unique(c(unlist(ph_NA_Indices_Tab)))
	   trainPheno<- trainPheno0[-ph_NA_Indices,]
	   trainGeno_Imp <- trainGeno_Imp0[-ph_NA_Indices,]
	   Geno <- as.character(TrainData_Table_Num_Filt[-ph_NA_Indices,1])
	}
	if(!anyNA(trainPheno0)){
	  
	   trainPheno <- trainPheno0
	   trainGeno_Imp <- trainGeno_Imp0
	   Geno <- as.character(TrainData_Table_Num_Filt[,1])
	}
	 
####### Set the same set of markers for both train and test sets with the same order 
     ##Check with sorted markers and check with test set 
	  trainssIDs <- colnames(trainGeno_Imp)
	  testssIDs <- colnames(testGeno_Imp) 
	  testGeno_Imp_Ord <- testGeno_Imp[, match(trainssIDs,testssIDs)]
	  testGeno_Imp <- testGeno_Imp_Ord

## Kinship Matrix 
    # rownames(trainGeno_Imp) <- Geno
	 #rownames(trainPheno) <- Geno
	 A <- A.mat(trainGeno_Imp)	
	 colnames(A) <- Geno
	 rownames(A) <- Geno

	 Y <- trainPheno
	 X <- rbind(trainGeno_Imp,TestGenoTable)
	 rownames(X) <- c(Geno,rownames(TestGenoTable))
	 #rownames(X) <- c(rownames(trainGeno_Imp),rownames(TestGenoTable))
	 n <- nrow(X)
	 p <- ncol(X)
	 
	 A.Tot <- A.mat(X)	
	 colnames(A.Tot) <- rownames(X)
	 rownames(A.Tot) <- rownames(X)
	 
	 predEnviron <- PredEnviron
	 nTst <- nrow(TestGenoTable)
	 Ytst <- as.data.frame(matrix(rep(NA,nTst*ncol(Y)),nrow=nTst,ncol=ncol(Y)))
	 colnames(Ytst) <- colnames(Y)
	 Ytst$environ <- as.factor(rep(predEnviron,nTst))
	 yNA <- rbind.data.frame(Y,Ytst)
	 yNA_DF <- as.data.frame(yNA)
	 
	 #yNA_DF$id <- as.factor(rownames(X))
	 #yNA_DF$id <- as.factor(paste(c(Geno,rownames(TestGenoTable)),"_",yNA_DF[,1],sep=""))
	 yNA_DF$strain <- as.factor(c(Geno,rownames(TestGenoTable)))
     yNA_DF$environ <- as.factor(yNA_DF$environ)
	 
	 nTrt <- length(trait)
	 yNA_DF[,c(1:nTrt)] <- apply(yNA_DF[,c(1:nTrt)],2,as.numeric) 
	 
	 # tst <- c((length(trainSetID)+1):nrow(X))
	 # if(GPModelMT != "GBLUP (SOMMER)"){ 
	  # Test_PredictedValues <- fm3$ETAHat[tst,]
	 # }
	 # if(GPModelMT == "GBLUP (SOMMER)"){  
	  # UMat <- do.call(cbind,lapply(fm3$U$`u:id`,function(x) x))
	  # Test_PredictedValues <- c()
	  # for(nT in 1:length(trait)){
	    # Test_PredictedValues <- cbind(Test_PredictedValues,(UMat[,nT]+fm3$Beta[nT,3])[tst])
	  # }
	 # }
	 
	 
#### 
### a) Model with Main Effect

   if(GPModelVarCor == "Main Effect"){

       fitMain1 <- mmer(YieldBuA~factor(environ)-1,
                random=~vsr(strain,Gu=A.Tot),
                rcov=~units,
                data=yNA_DF,verbose=FALSE)
       summary(fitMain1)

	   m <- model.matrix(~factor(environ)-1 ,data=yNA_DF)
	   m1 <- model.matrix(~factor(strain)-1 ,data=yNA_DF)
	   m_beta <- m %*% as.numeric(fitMain1$Beta[,3]) 
	   PredMain1 <- m_beta + (m1%*% fitMain1$U$`u:strain`$YieldBuA)

	   #cor(PredMain1,yNA_DF[,"Yield_blup"]) 

      # plot(yNA_DF[,"Yield_blup"],PredMain,xlab="Observed",ylab="Predicted")

    }

### b) Model with Compound Symmetry Var-CoVar Structure 
  
   if(GPModelVarCor == "Homogeneous Variance (CS)"){
#E <- diag(length(unique(yNA_DF$environ)))

		E <- diag(length(levels(factor(yNA_DF$environ))))
		rownames(E) <- colnames(E) <- levels(factor(yNA_DF$environ))

		EA <- kronecker(E,A.Tot, make.dimnames = TRUE)
		yNA_DF$environ <- as.factor(yNA_DF$environ)
		yNA_DF$strain <- as.factor(yNA_DF$strain)

		fitCS <- mmer(Yield_blup~factor(environ)-1,
					  random= ~ vsr(strain, Gu=A.Tot) + vsr(factor(environ):strain, Gu=EA),
					  rcov= ~ units,tolPar=1e-03,tolParInv=10,
					  data=yNA_DF, verbose = TRUE)
					  

		summary(fitCS)

#####

		 m <- model.matrix(~ factor(environ)-1 ,data=yNA_DF)
		 m_beta <- m %*% as.numeric(fitCS$Beta[,3]) 

		 envNames <- levels(factor(yNA_DF$environ))
		 print(envNames)
		 strainSub <- lapply(envNames,function(x) yNA_DF[which((yNA_DF$environ) %in% x),"strain"])
		 strainEnv <- lapply(c(1:length(envNames)), function(x) paste(envNames[x],strainSub[[x]],sep=":"))
		 strainEnv_Unq <- lapply(strainEnv,unique)
		 
		 strSubInd <- lapply(strainEnv_Unq,function(x) (lapply(x, function(y) grep(y,names(fitCS$U$`u:environ:strain`$Yield_blup)))))
		 length(unlist(strSubInd))
		 unique(names(fitCS$U$`u:environ:strain`$Yield_blup))[unique(unlist(strSubInd))]

		 DT_Sub2 <- cbind.data.frame(yNA_DF,paste(yNA_DF$environ,yNA_DF$strain,sep=":"))
		 colnames(DT_Sub2)[ncol(DT_Sub2)] <- "environ-strain"
		 DT_Sub2A <- DT_Sub2[-which(duplicated(DT_Sub2$`environ-strain`)),]
		 dim(DT_Sub2A)
		 
		 length(which(names(PredCS) %in% DT_Sub2A$`environ-strain`))
		 
		 strEnvSubInd <- lapply(envNames, function(y) grep(y,names(fitCS$U$`u:environ:strain`$Yield_blup)))
		 
		 rmStrain <- setdiff(unique(names(fitCS$U$`u:environ:strain`$Yield_blup))[unique(unlist(strSubInd))],unlist(strainEnv_Unq))
		 # [1] "Wanatah_IN-2020:A14019141"     "Williamston_MI-2020:A14019141"
		 
		 strainEnv_out <- unique(names(fitCS$U$`u:environ:strain`$Yield_blup))[unique(unlist(strSubInd))]
		 strainEnv_unq_out <-  strainEnv_out[-which(unique(names(fitCS$U$`u:environ:strain`$Yield_blup))[unique(unlist(strSubInd))] %in% rmStrain)]

		 which(unique(names(fitCS$U$`u:environ:strain`$Yield_blup))[unique(unlist(strSubInd))] %in% rmStrain)

		 strainEnv_Unq_rm <-  strainEnv_Unq[which(unlist(strainEnv_Unq) %in% rmStrain)]
		 if(length(strainEnv_Unq_rm) ==0){ 
		 
				strainEnv_Unq_rm <- strainEnv_Unq
		 }

		 strSubInd <- lapply(strainEnv_unq_out,function(x) (lapply(x, function(y) grep(y,names(fitCS$U$`u:environ:strain`$Yield_blup)))))
		 length(unlist(strSubInd))
		 strSub<- lapply(strainEnv_unq_out,function(x) (lapply(x, function(y) grep(y,names(fitCS$U$`u:environ:strain`$Yield_blup),value=TRUE))))
		 
		 PredCS <- unlist(lapply(c(1:length(envNames)),function(x) as.numeric(fitCS$Beta[x,3]) + fitCS$U$`u:environ:strain`$Yield_blup[strEnvSubInd[[x]]]))

# strSubInd <- lapply(strainEnv_unq_out,function(y) grep(y,names(PredCS)))

		 
		 PredCS_DTInd <- (which(names(PredCS) %in% DT_Sub2A$`environ-strain`))
		 PredCS_DT_Sub2A <- PredCS[PredCS_DTInd]

		 PredCS_DT_Sub2A_Sort <- PredCS_DT_Sub2A[match(DT_Sub2A[,"environ-strain"],names(PredCS_DT_Sub2A))]
	     cor(PredCS_DT_Sub2A_Sort,DT_Sub2A[,"Yield_blup"])
# 0.9312556
          plot(PredCS,yNA_DF[,"Yield_blup"])
 }
 
#### c) Model with CS+DG (Heterogeneous Var-Covar 
 
   if(GPModelVarCor == "Heterogeneous Variance (CS+DG)"){


		fitCSDG <- mmer(YieldBuA~environ-1,
                random=~vsr(strain,Gu=A.Tot) +vsr(dsr(as.factor(environ)),strain,Gu=A.Tot),
                rcov=~units,tolPar=1e-03,tolParInv=1e-09,
                data=yNA_DF,verbose=TRUE) 
				

		
				
# Error in mmer(Yield_blup ~ environ - 1, random = ~vsr(strain, Gu = A.Tot) +  :
# Cube::init(): requested size is too large; suggest to enable ARMA_64BIT_WORD
# In addition: Warning message:
# fixed-effect model matrix is rank deficient so dropping 1262 columns / coefficients


# Error: Mat::init(): requested size is too large; suggest to enable
# | ARMA_64BIT_WORD or compiling usin C++11
# To enable ARMA_64BIT_WORD in RcppArmadillo on a locally installed version or RcppArmadillo we did

# wget https://cran.r-project.org/src/contrib/RcppArmadillo_0.11.4.0.1.tar.gz
# tar -xvzf RcppArmadillo_0.9.100.5.0.tar.gz
# vi RcppArmadillo/src/Makefile.in
# Add the line PKG_CPPFLAGS = -DARMA_64BIT_WORD=1 after the first line of the file. Then sta
# install.packages("RcppArmadillo", repos=NULL)



		geno_imp <- markov(apply(geno_num5[,-1],2,as.numeric))
		rownames(geno_imp) <- geno_num5[,"strain"]
		dim(geno_imp)
		env_geno_sub_indices <- which(rownames(geno_imp) %in% unique(DT[,"strain"]))
		geno_imp_sub <- geno_imp[env_geno_sub_indices,]

		dim(geno_imp_sub)

		Pheno_DT <- as.matrix(DT_Sub1[,"Yield_blup"]) 
		rownames(Pheno_DT) <- DT_Sub1[,"strain"] 

		genoQC1 <- cleanREPV2(Pheno_DT,geno_imp_sub)

		#### 

		K_rr <- A.mat(genoQC1[[2]])
		colnames(K_rr) <-rownames(genoQC1[[2]])
		rownames(K_rr) <- rownames(genoQC1[[2]])
		A <- K_rr
		dim(A)

		####


		K_rr <- A.mat(geno_imp_sub)
		colnames(K_rr) <-rownames(geno_imp_sub)
		rownames(K_rr) <- rownames(geno_imp_sub)
		A <- K_rr
		dim(A)

		#A_Sub <- A[1:500,1:500]

		#DT_Pheno <- Pheno_DT[which(rownames(Pheno_DT) %in% rownames(A_Sub)),]
		A_Sub <- A+diag(nrow(A))

		DT_Sub <- DT[which(DT[,"strain"] %in% rownames(A_Sub)),]

		E <- diag(length(unique(DT_Sub$environ)))
		rownames(E) <- colnames(E) <- unique(DT_Sub$environ)
		dim(E)

		### Same set of strains in each of the environments 
		# 
		# rmStrains <- names(which(table(DT_Sub[,"strain"]) <45))
		# DT_Sub1 <- DT_Sub[-which(DT_Sub[,"strain"] %in% rmStrains),]
		# 
		# A.Tot <- A_Sub[-which(rownames(A_Sub) %in% rmStrains),-which(rownames(A_Sub) %in% rmStrains)]
		# dim(A.Tot) 

		 EFactors <- levels(factor(DT_Sub$environ)) 
		 EFactors2019 <- EFactors[ grep("2019",levels(factor(DT_Sub$environ)))]
		 EFactors2020 <- EFactors[ grep("2020",levels(factor(DT_Sub$environ)))] 

		 E2019Indx <-  grep("2019",levels(factor(DT_Sub$environ)))
		 E2020Indx <-  grep("2020",levels(factor(DT_Sub$environ))) 
		  
		 Efactors_Sub1 <- c(EFactors2019,EFactors2020)
		 Sub1Indx <- which(as.character(DT_Sub$environ) %in% Efactors_Sub1)
		 #Sub1Indx <- which(as.character(DT_Sub$environ) %in% EFactors2020)


		### NUST - GxE 

		DT_Sub1 <- DT_Sub[Sub1Indx,]
		A_Sub_FiltIndx <- which(rownames(A_Sub) %in% as.character(DT_Sub1[,"strain"]))
		A.Tot <- A_Sub[A_Sub_FiltIndx,A_Sub_FiltIndx]

		fitCSDG1 <- mmer(Yield_blup~factor(environ)-1,
						random=~vsr((strain),Gu=A.Tot) +vsr(dsr(factor(environ)),(strain),Gu=A.Tot),
						rcov=~units,tolParInv=1e-03,
						data=DT_Sub1,verbose=TRUE) 

		summary(fitCSDG) 


  
		  
		m2 <- model.matrix(~factor(environ)-1 ,data=yNA_DF)
		m_beta <- m2 %*% as.numeric(fitCSDG$Beta[,3]) 
		length(m_beta)

		m_env_strain <- do.call(cbind,lapply(fitCSDG$U,function(x) x$Yield_blup))
		dim(m_env_strain)
		envStrain_blup <-c(m_env_strain[,2:4])                              
				

		m1 <- model.matrix(~factor(strain)-1,data=yNA_DF)

		YldBlp <- as.matrix(fitCSDG$U$`u:strain`$Yield_blup)

		YldBlp_Strain <- m1%*% YldB
		Pred <- m_beta+YldBlp_Strain[,1]
		
		cor(Pred,DT_Sub1[,"Yield_blup"])
### 
		summary(fitCSDG1)

		 
		m2 <- model.matrix(~factor(environ)-1 ,data=DT_Sub1)
		m_beta <- m2 %*% as.numeric(fitCSDG1$Beta[,3]) 
		length(m_beta)

		m_env_strain <- do.call(cbind,lapply(fitCSDG1$U,function(x) x$Yield_blup))
		dim(m_env_strain)
		envStrain_blup <-c(m_env_strain[,2:4])                              
				

		m1 <- model.matrix(~factor(strain)-1,data=DT_Sub1)

		YldBlp <- as.matrix(fitCSDG1$U$`u:strain`$Yield_blup)

		YldBlp_Strain <- m1%*% YldBlp
		Pred <- m_beta+YldBlp_Strain[,1]
		length(Pred)
		cor(Pred,DT_Sub1[,"Yield_blup"])

}

#### d) Model with unstructured Var-Covar

   if(GPModelVarCor == "Unstructured (US)"){
	
  	 fitUS <- mmer(YieldBuA~factor(environ)-1,
					random=~vsr(usr(factor(environ)),strain,Gu=A.Tot),
					rcov=~units,tolParInv=10,
					data=yNA_DF,verbose=TRUE) 
	 summary(fitUS) 
  
	 envNames <- levels(factor(yNA_DF$environ))
	 print(envNames)
	 print(names(fitUS$U))

	#env1Ind <- c(1,3,6)
	 EnvSplit <- strsplit(names(fitUS$U),":")
	 env1Ind <- which(unlist(lapply(EnvSplit,length))==2)

	 U_envStrain <- list()
	 PredUS <- list()
		for(i in 1:length(envNames)){
		   envInd <-  grep(envNames[i],names(fitUS$U))
		   envIndNames <-  grep(envNames[i],names(fitUS$U),value=TRUE)
		   strainSub <- yNA_DF[which(factor(yNA_DF$environ) %in% envNames[i]),"strain"]
		   strSubInd <- which(names(fitUS$U[[env1Ind[i]]]$Yield_blup) %in% strainSub)
		   U_envStrain[[i]] <-  as.numeric(fitUS$U[[env1Ind[i]]]$Yield_blup)[strSubInd]
		   
		   for(j  in 2:length(envInd)){ 
			 indJ <- envInd[j]
			 b <- cbind(names(fitUS$U[[indJ]]$Yield_blup),fitUS$U[[indJ]]$Yield_blup)
			 colnames(b) <- c("strain","Yield_blup")
			 b_sub <- b[which(b[,"strain"] %in% strainSub),]
			 b_group <- as_tibble(b_sub) %>% group_by(strain)
			 YldBlup_group <- b_group %>% summarise(Yield_blup = sum(as.numeric(Yield_blup)))
			 U_envStrain[[i]] <- U_envStrain[[i]] +YldBlup_group[,2]
			
		   }
		  # names(U_envStrain[[i]][[1]]) <- unlist(YldBlup_group[,1])
		   
		   PredUS[[i]] <- c(U_envStrain[[i]][[1]] + fitUS$Beta[i,3])
		   names(PredUS[[i]]) <-  paste(envNames[i],unlist(YldBlup_group[,1]),sep=":")
		}

	 PredUS1 <- unlist(PredUS)
	 
	 PredUS_DTInd <- (which(names(PredUS1) %in% DT_Sub2A$`environ-strain`))
	 PredUS_DT_Sub2A <- PredUS1[PredUS_DTInd]
	 PredUS_DT_Sub2A_Sort <- PredUS_DT_Sub2A[match(DT_Sub2A[,"environ-strain"],names(PredUS_DT_Sub2A))]
	 cor(unlist(PredUS_DT_Sub2A_Sort),DT_Sub2A[,"Yield_blup"]) 
	 
	}

	Test_SortedPredictedValues <- sort.int(Test_PredictedValues[,1],decreasing=TRUE,index.return=TRUE) 
	Test_StrainID <- rownames(TestData_Table_Num_Filt)	
	
### Upper bound of Reliability

   # trainGeno_Imp2 <- apply(trainGeno_Imp,2,function(x) x+1)  
   
    cleanData <- cleanREP(trainPheno,trainGeno_Imp)
    M <-  cleanData[[2]]
    M.Pdt <- t(M)%*% solve(M %*% t(M) + diag(nrow(M))) %*% M
      
    getU <- function(M.Pdt,v){
	   v.hat <- M.Pdt %*% v
	   U <- (t(v.hat)%*% v.hat)/ (t(v) %*% v)
	   return(U)
    }
    U <- apply(TestGenoTable,1,function(x) getU(M.Pdt,x)) 

### Sort Output
	
	Test_SortedIndices <- Test_SortedPredictedValues[[2]]
	Sorted_Test_StrainID <- Test_StrainID[Test_SortedIndices]
	U_Sorted <- U[(Test_SortedIndices)]
	
	Test_SortedPredValuesTab <- apply(Test_PredictedValues[Test_SortedIndices,],2,function(x) round(x,digits=2))
	
### Output DF	
	outputDF <- cbind.data.frame(Sorted_Test_StrainID ,Test_SortedPredValuesTab,round(U_Sorted,digits=5))
	colnames(outputDF) <- c("LineID",colnames(Test_SortedPredValuesTab),"Upper Bound of Reliability")
	
	
   return(outputDF)
	
 }



#### 



getRankedPredictedValuesME_BGGE <- function(Data_Table_Num_Filt_List,nTraits,trait,GPModelVarCor,optTS=NULL){ 
   
	 TrainData_Table_Num_Filt <- Data_Table_Num_Filt_List[[1]]
	 TestData_Table_Num_Filt <- Data_Table_Num_Filt_List[[2]]
	 			 
	 TestGenoTable <- apply(TestData_Table_Num_Filt[,-1],2,as.numeric)
	 rownames(TestGenoTable) <- TestData_Table_Num_Filt[,1]
	 
	 set.seed(125)
	 strainID <- as.character(TrainData_Table_Num_Filt[,1])
	 
	 initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	 
	 # trainSetID <- as.character(as.vector( strainID[sample(c(1:length(strainID)),500)]))
	 # optTS <- cbind(trainSetID,c(1:500))
	 
	 trainSetID <- as.character(as.vector(strainID))
	 optTS <- cbind(trainSetID,c(1:length(trainSetID)))
	 
	 initCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[1]
	 finalCol <- grep("ss",colnames(TrainData_Table_Num_Filt))[length(grep("ss",colnames(TrainData_Table_Num_Filt)))]
	 
	
	 if(is.null(optTS)){
	
		 trainPheno0 <- TrainData_Table_Num_Filt[,c(trait,"environ")]
		 geno_012 <- apply(TrainData_Table_Num_Filt[,c(initCol:finalCol)],2,as.numeric)
		 trainGeno <- apply(geno_012,2,function(x) x-1)
		 #Geno <- rownames(TrainData_Table_Num_Filt)
		 Geno <- as.character(TrainData_Table_Num_Filt[,1])
	 }
	 
	 if(!is.null(optTS)){
	     
	     strainID <- as.character(TrainData_Table_Num_Filt[,1])
		 trainSetID <- as.character(as.vector(optTS[,1]))
	     trainIndices <- which(strainID %in% trainSetID)
		 trainPheno0 <- TrainData_Table_Num_Filt[trainIndices,c(trait,"environ")]
		 geno_012 <- apply(TrainData_Table_Num_Filt[trainIndices,c(initCol:finalCol)],2,as.numeric)
		 trainGeno <- apply(geno_012,2,function(x) x-1)
		 Geno <- as.character(TrainData_Table_Num_Filt[trainIndices,1])
	 }
	  
## Prepare Data Table for GP 
	 # Data <- cbind.data.frame(Geno,trainPheno)
	 # colnames(Data) <- c("Geno","Pheno")
	 # Geno <- "Geno"
	 # Pheno <- "Pheno"
	

#### Impute trainGeno and train using mixed.solve
	
	# if(anyNA(trainGeno)){	 
		# trainGeno_Imp <- snpQC(trainGeno,impute=TRUE,remove=FALSE)
	# }else{trainGeno_Imp <- trainGeno}
	
	
### Check if the markers are missing in > 99% of the lines 
	 
	  allNAMarkers <- apply(trainGeno,2,function(x) if(length(which(is.na(x)))==nrow(trainGeno)){1}else{0})
     
	  allNAMarkerIndices <- which(allNAMarkers!=0)
	  
	  if(length(allNAMarkerIndices)>=1){
		 trainGeno0 <- trainGeno[,-allNAMarkerIndices]
	  }else{ trainGeno0 <- trainGeno}
	 	 
	 
	 if(anyNA(trainGeno0)){	 
		trainGeno_Imp0 <- snpQC(trainGeno0,impute=TRUE,remove=FALSE)
	 }else{trainGeno_Imp0 <- trainGeno0}
	 
	 
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
	 
	

 ### Remove lines with missing pheno	 
	  
	if(anyNA(trainPheno0)){
	   ph_NA_Indices_Tab <- apply(trainPheno0,2,function(x) which(is.na(x)))
	   ph_NA_Indices <- unique(c(unlist(ph_NA_Indices_Tab)))
	   trainPheno<- trainPheno0[-ph_NA_Indices,]
	   trainGeno_Imp <- trainGeno_Imp0[-ph_NA_Indices,]
	   Geno <- as.character(TrainData_Table_Num_Filt[-ph_NA_Indices,1])
	}
	if(!anyNA(trainPheno0)){
	  
	   trainPheno <- trainPheno0
	   trainGeno_Imp <- trainGeno_Imp0
	   Geno <- as.character(TrainData_Table_Num_Filt[,1])
	}
	 
####### Set the same set of markers for both train and test sets with the same order 
     ##Check with sorted markers and check with test set 
	  trainssIDs <- colnames(trainGeno_Imp)
	  testssIDs <- colnames(testGeno_Imp) 
	  testGeno_Imp_Ord <- testGeno_Imp[, match(trainssIDs,testssIDs)]
	  testGeno_Imp <- testGeno_Imp_Ord

## Kinship Matrix 
    # rownames(trainGeno_Imp) <- Geno
	 #rownames(trainPheno) <- Geno
	 A <- A.mat(trainGeno_Imp)	
	 colnames(A) <- Geno
	 rownames(A) <- Geno

	 Y <- trainPheno
	 X <- rbind(trainGeno_Imp,TestGenoTable)
	 rownames(X) <- c(Geno,rownames(TestGenoTable))
	 #rownames(X) <- c(rownames(trainGeno_Imp),rownames(TestGenoTable))
	 n <- nrow(X)
	 p <- ncol(X)
	 
	 A.Tot <- A.mat(X)	
	 colnames(A.Tot) <- rownames(X)
	 rownames(A.Tot) <- rownames(X)
	 
	 predEnviron <- PredEnviron
	 nTst <- nrow(TestGenoTable)
	 Ytst <- as.data.frame(matrix(rep(NA,nTst*ncol(Y)),nrow=nTst,ncol=ncol(Y)))
	 colnames(Ytst) <- colnames(Y)
	 Ytst$environ <- as.factor(rep(predEnviron,nTst))
	 yNA <- rbind.data.frame(Y,Ytst)
	 yNA_DF <- as.data.frame(yNA)
	 
	 #yNA_DF$id <- as.factor(rownames(X))
	 #yNA_DF$id <- as.factor(paste(c(Geno,rownames(TestGenoTable)),"_",yNA_DF[,1],sep=""))
	 yNA_DF$strain <- as.factor(c(Geno,rownames(TestGenoTable)))
     yNA_DF$environ <- as.factor(yNA_DF$environ)
	 
	 nTrt <- length(trait)
	 yNA_DF[,c(1:nTrt)] <- apply(yNA_DF[,c(1:nTrt)],2,as.numeric) 
	 
	 # tst <- c((length(trainSetID)+1):nrow(X))
	 # if(GPModelMT != "GBLUP (SOMMER)"){ 
	  # Test_PredictedValues <- fm3$ETAHat[tst,]
	 # }
	 # if(GPModelMT == "GBLUP (SOMMER)"){  
	  # UMat <- do.call(cbind,lapply(fm3$U$`u:id`,function(x) x))
	  # Test_PredictedValues <- c()
	  # for(nT in 1:length(trait)){
	    # Test_PredictedValues <- cbind(Test_PredictedValues,(UMat[,nT]+fm3$Beta[nT,3])[tst])
	  # }
	 # }
	 
	 
#### 
### a) Model with Main Effect

   if(GPModelVarCor == "Main Effect"){

       fitMain1 <- mmer(YieldBuA~factor(environ)-1,
                random=~vsr(strain,Gu=A.Tot),
                rcov=~units,
                data=yNA_DF,verbose=FALSE)
       summary(fitMain1)

	   m <- model.matrix(~factor(environ)-1 ,data=yNA_DF)
	   m1 <- model.matrix(~factor(strain)-1 ,data=yNA_DF)
	   m_beta <- m %*% as.numeric(fitMain1$Beta[,3]) 
	   PredMain1 <- m_beta + (m1%*% fitMain1$U$`u:strain`$YieldBuA)

	   #cor(PredMain1,yNA_DF[,"Yield_blup"]) 

      # plot(yNA_DF[,"Yield_blup"],PredMain,xlab="Observed",ylab="Predicted")

    }

### b) Model with Compound Symmetry Var-CoVar Structure 
  
   if(GPModelVarCor == "Homogeneous Variance (CS)"){
#E <- diag(length(unique(yNA_DF$environ)))

		E <- diag(length(levels(factor(yNA_DF$environ))))
		rownames(E) <- colnames(E) <- levels(factor(yNA_DF$environ))

		EA <- kronecker(E,A.Tot, make.dimnames = TRUE)
		yNA_DF$environ <- as.factor(yNA_DF$environ)
		yNA_DF$strain <- as.factor(yNA_DF$strain)

		fitCS <- mmer(Yield_blup~factor(environ)-1,
					  random= ~ vsr(strain, Gu=A.Tot) + vsr(factor(environ):strain, Gu=EA),
					  rcov= ~ units,tolPar=1e-03,tolParInv=10,
					  data=yNA_DF, verbose = TRUE)
					  

		summary(fitCS)

#####

		 m <- model.matrix(~ factor(environ)-1 ,data=yNA_DF)
		 m_beta <- m %*% as.numeric(fitCS$Beta[,3]) 

		 envNames <- levels(factor(yNA_DF$environ))
		 print(envNames)
		 strainSub <- lapply(envNames,function(x) yNA_DF[which((yNA_DF$environ) %in% x),"strain"])
		 strainEnv <- lapply(c(1:length(envNames)), function(x) paste(envNames[x],strainSub[[x]],sep=":"))
		 strainEnv_Unq <- lapply(strainEnv,unique)
		 
		 strSubInd <- lapply(strainEnv_Unq,function(x) (lapply(x, function(y) grep(y,names(fitCS$U$`u:environ:strain`$Yield_blup)))))
		 length(unlist(strSubInd))
		 unique(names(fitCS$U$`u:environ:strain`$Yield_blup))[unique(unlist(strSubInd))]

		 DT_Sub2 <- cbind.data.frame(yNA_DF,paste(yNA_DF$environ,yNA_DF$strain,sep=":"))
		 colnames(DT_Sub2)[ncol(DT_Sub2)] <- "environ-strain"
		 DT_Sub2A <- DT_Sub2[-which(duplicated(DT_Sub2$`environ-strain`)),]
		 dim(DT_Sub2A)
		 
		 length(which(names(PredCS) %in% DT_Sub2A$`environ-strain`))
		 
		 strEnvSubInd <- lapply(envNames, function(y) grep(y,names(fitCS$U$`u:environ:strain`$Yield_blup)))
		 
		 rmStrain <- setdiff(unique(names(fitCS$U$`u:environ:strain`$Yield_blup))[unique(unlist(strSubInd))],unlist(strainEnv_Unq))
		 # [1] "Wanatah_IN-2020:A14019141"     "Williamston_MI-2020:A14019141"
		 
		 strainEnv_out <- unique(names(fitCS$U$`u:environ:strain`$Yield_blup))[unique(unlist(strSubInd))]
		 strainEnv_unq_out <-  strainEnv_out[-which(unique(names(fitCS$U$`u:environ:strain`$Yield_blup))[unique(unlist(strSubInd))] %in% rmStrain)]

		 which(unique(names(fitCS$U$`u:environ:strain`$Yield_blup))[unique(unlist(strSubInd))] %in% rmStrain)

		 strainEnv_Unq_rm <-  strainEnv_Unq[which(unlist(strainEnv_Unq) %in% rmStrain)]
		 if(length(strainEnv_Unq_rm) ==0){ 
		 
				strainEnv_Unq_rm <- strainEnv_Unq
		 }

		 strSubInd <- lapply(strainEnv_unq_out,function(x) (lapply(x, function(y) grep(y,names(fitCS$U$`u:environ:strain`$Yield_blup)))))
		 length(unlist(strSubInd))
		 strSub<- lapply(strainEnv_unq_out,function(x) (lapply(x, function(y) grep(y,names(fitCS$U$`u:environ:strain`$Yield_blup),value=TRUE))))
		 
		 PredCS <- unlist(lapply(c(1:length(envNames)),function(x) as.numeric(fitCS$Beta[x,3]) + fitCS$U$`u:environ:strain`$Yield_blup[strEnvSubInd[[x]]]))

# strSubInd <- lapply(strainEnv_unq_out,function(y) grep(y,names(PredCS)))

		 
		 PredCS_DTInd <- (which(names(PredCS) %in% DT_Sub2A$`environ-strain`))
		 PredCS_DT_Sub2A <- PredCS[PredCS_DTInd]

		 PredCS_DT_Sub2A_Sort <- PredCS_DT_Sub2A[match(DT_Sub2A[,"environ-strain"],names(PredCS_DT_Sub2A))]
	     cor(PredCS_DT_Sub2A_Sort,DT_Sub2A[,"Yield_blup"])
# 0.9312556
          plot(PredCS,yNA_DF[,"Yield_blup"])
 }
 
#### c) Model with CS+DG (Heterogeneous Var-Covar 
 
   if(GPModelVarCor == "Heterogeneous Variance (CS+DG)"){


		fitCSDG <- mmer(YieldBuA~environ-1,
                random=~vsr(strain,Gu=A.Tot) +vsr(dsr(as.factor(environ)),strain,Gu=A.Tot),
                rcov=~units,tolPar=1e-03,tolParInv=1e-09,
                data=yNA_DF,verbose=TRUE) 
				

		
				
# Error in mmer(Yield_blup ~ environ - 1, random = ~vsr(strain, Gu = A.Tot) +  :
# Cube::init(): requested size is too large; suggest to enable ARMA_64BIT_WORD
# In addition: Warning message:
# fixed-effect model matrix is rank deficient so dropping 1262 columns / coefficients


# Error: Mat::init(): requested size is too large; suggest to enable
# | ARMA_64BIT_WORD or compiling usin C++11
# To enable ARMA_64BIT_WORD in RcppArmadillo on a locally installed version or RcppArmadillo we did

# wget https://cran.r-project.org/src/contrib/RcppArmadillo_0.11.4.0.1.tar.gz
# tar -xvzf RcppArmadillo_0.9.100.5.0.tar.gz
# vi RcppArmadillo/src/Makefile.in
# Add the line PKG_CPPFLAGS = -DARMA_64BIT_WORD=1 after the first line of the file. Then sta
# install.packages("RcppArmadillo", repos=NULL)



		geno_imp <- markov(apply(geno_num5[,-1],2,as.numeric))
		rownames(geno_imp) <- geno_num5[,"strain"]
		dim(geno_imp)
		env_geno_sub_indices <- which(rownames(geno_imp) %in% unique(DT[,"strain"]))
		geno_imp_sub <- geno_imp[env_geno_sub_indices,]

		dim(geno_imp_sub)

		Pheno_DT <- as.matrix(DT_Sub1[,"Yield_blup"]) 
		rownames(Pheno_DT) <- DT_Sub1[,"strain"] 

		genoQC1 <- cleanREPV2(Pheno_DT,geno_imp_sub)

		#### 

		K_rr <- A.mat(genoQC1[[2]])
		colnames(K_rr) <-rownames(genoQC1[[2]])
		rownames(K_rr) <- rownames(genoQC1[[2]])
		A <- K_rr
		dim(A)

		####


		K_rr <- A.mat(geno_imp_sub)
		colnames(K_rr) <-rownames(geno_imp_sub)
		rownames(K_rr) <- rownames(geno_imp_sub)
		A <- K_rr
		dim(A)

		#A_Sub <- A[1:500,1:500]

		#DT_Pheno <- Pheno_DT[which(rownames(Pheno_DT) %in% rownames(A_Sub)),]
		A_Sub <- A+diag(nrow(A))

		DT_Sub <- DT[which(DT[,"strain"] %in% rownames(A_Sub)),]

		E <- diag(length(unique(DT_Sub$environ)))
		rownames(E) <- colnames(E) <- unique(DT_Sub$environ)
		dim(E)

		### Same set of strains in each of the environments 
		# 
		# rmStrains <- names(which(table(DT_Sub[,"strain"]) <45))
		# DT_Sub1 <- DT_Sub[-which(DT_Sub[,"strain"] %in% rmStrains),]
		# 
		# A.Tot <- A_Sub[-which(rownames(A_Sub) %in% rmStrains),-which(rownames(A_Sub) %in% rmStrains)]
		# dim(A.Tot) 

		 EFactors <- levels(factor(DT_Sub$environ)) 
		 EFactors2019 <- EFactors[ grep("2019",levels(factor(DT_Sub$environ)))]
		 EFactors2020 <- EFactors[ grep("2020",levels(factor(DT_Sub$environ)))] 

		 E2019Indx <-  grep("2019",levels(factor(DT_Sub$environ)))
		 E2020Indx <-  grep("2020",levels(factor(DT_Sub$environ))) 
		  
		 Efactors_Sub1 <- c(EFactors2019,EFactors2020)
		 Sub1Indx <- which(as.character(DT_Sub$environ) %in% Efactors_Sub1)
		 #Sub1Indx <- which(as.character(DT_Sub$environ) %in% EFactors2020)


		### NUST - GxE 

		DT_Sub1 <- DT_Sub[Sub1Indx,]
		A_Sub_FiltIndx <- which(rownames(A_Sub) %in% as.character(DT_Sub1[,"strain"]))
		A.Tot <- A_Sub[A_Sub_FiltIndx,A_Sub_FiltIndx]

		fitCSDG1 <- mmer(Yield_blup~factor(environ)-1,
						random=~vsr((strain),Gu=A.Tot) +vsr(dsr(factor(environ)),(strain),Gu=A.Tot),
						rcov=~units,tolParInv=1e-03,
						data=DT_Sub1,verbose=TRUE) 

		summary(fitCSDG) 


  
		  
		m2 <- model.matrix(~factor(environ)-1 ,data=yNA_DF)
		m_beta <- m2 %*% as.numeric(fitCSDG$Beta[,3]) 
		length(m_beta)

		m_env_strain <- do.call(cbind,lapply(fitCSDG$U,function(x) x$Yield_blup))
		dim(m_env_strain)
		envStrain_blup <-c(m_env_strain[,2:4])                              
				

		m1 <- model.matrix(~factor(strain)-1,data=yNA_DF)

		YldBlp <- as.matrix(fitCSDG$U$`u:strain`$Yield_blup)

		YldBlp_Strain <- m1%*% YldB
		Pred <- m_beta+YldBlp_Strain[,1]
		
		cor(Pred,DT_Sub1[,"Yield_blup"])
### 
		summary(fitCSDG1)

		 
		m2 <- model.matrix(~factor(environ)-1 ,data=DT_Sub1)
		m_beta <- m2 %*% as.numeric(fitCSDG1$Beta[,3]) 
		length(m_beta)

		m_env_strain <- do.call(cbind,lapply(fitCSDG1$U,function(x) x$Yield_blup))
		dim(m_env_strain)
		envStrain_blup <-c(m_env_strain[,2:4])                              
				

		m1 <- model.matrix(~factor(strain)-1,data=DT_Sub1)

		YldBlp <- as.matrix(fitCSDG1$U$`u:strain`$Yield_blup)

		YldBlp_Strain <- m1%*% YldBlp
		Pred <- m_beta+YldBlp_Strain[,1]
		length(Pred)
		cor(Pred,DT_Sub1[,"Yield_blup"])

}

#### d) Model with unstructured Var-Covar

   if(GPModelVarCor == "Unstructured (US)"){
	
  	 fitUS <- mmer(YieldBuA~factor(environ)-1,
					random=~vsr(usr(factor(environ)),strain,Gu=A.Tot),
					rcov=~units,tolParInv=10,
					data=yNA_DF,verbose=TRUE) 
	 summary(fitUS) 
  
	 envNames <- levels(factor(yNA_DF$environ))
	 print(envNames)
	 print(names(fitUS$U))

	#env1Ind <- c(1,3,6)
	 EnvSplit <- strsplit(names(fitUS$U),":")
	 env1Ind <- which(unlist(lapply(EnvSplit,length))==2)

	 U_envStrain <- list()
	 PredUS <- list()
		for(i in 1:length(envNames)){
		   envInd <-  grep(envNames[i],names(fitUS$U))
		   envIndNames <-  grep(envNames[i],names(fitUS$U),value=TRUE)
		   strainSub <- yNA_DF[which(factor(yNA_DF$environ) %in% envNames[i]),"strain"]
		   strSubInd <- which(names(fitUS$U[[env1Ind[i]]]$Yield_blup) %in% strainSub)
		   U_envStrain[[i]] <-  as.numeric(fitUS$U[[env1Ind[i]]]$Yield_blup)[strSubInd]
		   
		   for(j  in 2:length(envInd)){ 
			 indJ <- envInd[j]
			 b <- cbind(names(fitUS$U[[indJ]]$Yield_blup),fitUS$U[[indJ]]$Yield_blup)
			 colnames(b) <- c("strain","Yield_blup")
			 b_sub <- b[which(b[,"strain"] %in% strainSub),]
			 b_group <- as_tibble(b_sub) %>% group_by(strain)
			 YldBlup_group <- b_group %>% summarise(Yield_blup = sum(as.numeric(Yield_blup)))
			 U_envStrain[[i]] <- U_envStrain[[i]] +YldBlup_group[,2]
			
		   }
		  # names(U_envStrain[[i]][[1]]) <- unlist(YldBlup_group[,1])
		   
		   PredUS[[i]] <- c(U_envStrain[[i]][[1]] + fitUS$Beta[i,3])
		   names(PredUS[[i]]) <-  paste(envNames[i],unlist(YldBlup_group[,1]),sep=":")
		}

	 PredUS1 <- unlist(PredUS)
	 
	 PredUS_DTInd <- (which(names(PredUS1) %in% DT_Sub2A$`environ-strain`))
	 PredUS_DT_Sub2A <- PredUS1[PredUS_DTInd]
	 PredUS_DT_Sub2A_Sort <- PredUS_DT_Sub2A[match(DT_Sub2A[,"environ-strain"],names(PredUS_DT_Sub2A))]
	 cor(unlist(PredUS_DT_Sub2A_Sort),DT_Sub2A[,"Yield_blup"]) 
	 
	}

	Test_SortedPredictedValues <- sort.int(Test_PredictedValues[,1],decreasing=TRUE,index.return=TRUE) 
	Test_StrainID <- rownames(TestData_Table_Num_Filt)	
	
### Upper bound of Reliability

   # trainGeno_Imp2 <- apply(trainGeno_Imp,2,function(x) x+1)  
   
    cleanData <- cleanREP(trainPheno,trainGeno_Imp)
    M <-  cleanData[[2]]
    M.Pdt <- t(M)%*% solve(M %*% t(M) + diag(nrow(M))) %*% M
      
    getU <- function(M.Pdt,v){
	   v.hat <- M.Pdt %*% v
	   U <- (t(v.hat)%*% v.hat)/ (t(v) %*% v)
	   return(U)
    }
    U <- apply(TestGenoTable,1,function(x) getU(M.Pdt,x)) 

### Sort Output
	
	Test_SortedIndices <- Test_SortedPredictedValues[[2]]
	Sorted_Test_StrainID <- Test_StrainID[Test_SortedIndices]
	U_Sorted <- U[(Test_SortedIndices)]
	
	Test_SortedPredValuesTab <- apply(Test_PredictedValues[Test_SortedIndices,],2,function(x) round(x,digits=2))
	
### Output DF	
	outputDF <- cbind.data.frame(Sorted_Test_StrainID ,Test_SortedPredValuesTab,round(U_Sorted,digits=5))
	colnames(outputDF) <- c("LineID",colnames(Test_SortedPredValuesTab),"Upper Bound of Reliability")
	
	
   return(outputDF)
	
 }

####

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
	 
	  colnames(Data_Trt)[which(colnames(Data_Trt) %in% StrainME)] <- "Strain"
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
	  Ph_Mat_0 <- dcast(Strain ~ Loc, value.var = trait,fun.aggregate =mean,data=Data_Trt_List[[nTrt]])
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
	  ####
	  
	  Data <- merge(phData,genoImpDF,by=0)
	  dim(Data)
	 
	  genoCols <- grep("ss",colnames(Data))
	  genoDat <- as.matrix(apply(Data[,genoCols],2,as.numeric))
	  rownames(genoDat) <- Data[,1]
	  is.matrix(genoDat)
	  
	  ####
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
	  
	}

  return(list(DT_1_Filt_List,genoDat_List))

}





getPhenoDistbnPlots <- function(DT_1_Filt_List,TraitME,nTrt){

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
	  
	plot3 <- ggplot(DT_2B, aes(x = Loc, y = get(trait),fill=Loc)) +
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

     png(paste(trait," Distribution.png",sep=""),width=1024,height=768,pointsize=20)
	 print(plot3)
	 dev.off()

    return(plot3) 	 
	 
  }



getEnvData <- function(Coord_Data_Tab,startDate,endDate){
    
	Dat_Loc_Coords_Filter <- Coord_Data_Tab

	Dat_Env <- as.character(Dat_Loc_Coords_Filter[,"Location"])
    Dat_Count <- as.character(Dat_Loc_Coords_Filter[,"Country"])
	lat <- as.numeric(as.character(Dat_Loc_Coords_Filter[,"Lat"]))
	long <- as.numeric(as.character(Dat_Loc_Coords_Filter[,"Long"]))
	startDay <- rep(as.character(startDate),length(Dat_Env)) #as.vector(paste(yrNamesSpl,"-05-01",sep=""))
	endDay <- rep(as.character(endDate),length(Dat_Env)) #paste(yrNamesSpl,"-11-20",sep="")
		
       
    #EnvRtype::
    ### Get weather data for the loc coords within the time interval 

   

	env.data <- get_weather_terra(env.id = Dat_Env,lat = lat,
								  lon =long,start.day = startDay,
								  end.day = endDay, country = Dat_Count,parallel=FALSE)


    return(env.data)

}

###

getEnvKernel <- function(Env.Data,process=FALSE,Gaussian=FALSE){
   
    env.Dat <- Env.Data
	
	
	if(process==FALSE){
	   
		var3 <- c("T2M","T2M_MAX","T2M_MIN","PRECTOT","WS2M","RH2M","T2MDEW","ALLSKY_SFC_LW_DWN","ALLSKY_SFC_SW_DWN")
		EC.Dat = EnvRtype::W_matrix(env.data = env.Dat,env.id="env", var.id = var3) 
		KE <- list(W = EnvRtype::env_kernel(env.data = EC.Dat,gaussian=Gaussian)[[2]])
    }else if(process==TRUE){
	
		prData <- processWTH(env.data = env.Dat)
		var4 <- c("T2M","T2M_MAX","T2M_MIN","PRECTOT","WS2M","RH2M","T2MDEW","ALLSKY_SFC_LW_DWN","ALLSKY_SFC_SW_DWN","RTA","VPD","SPV","ETP","PETP","GDD","FRUE","T2M_RANGE")
			
		EC.Dat.Pr = EnvRtype::W_matrix(env.data = prData,env.id="env", var.id = var4) 
		rmECVar <- which(apply(EC.Dat.Pr,2,function(x) length(which(x %in% NaN)))!=0)
		EC.Dat.Pr.Filt <- EC.Dat.Pr[,-rmECVar]

		KE = list(W = env_kernel(env.data = EC.Dat.Pr.Filt,gaussian=Gaussian)[[2]])
	}
	
	return(KE)
   }
	
# ### Based on four variables

	# var2 = c("PRECTOT","T2M","ALLSKY_SFC_SW_DWN","RH2M")
	# EC.Dat = W_matrix(env.data = env.data.Dat,env.id="env", var.id = var2) 
	# KE.Raw <- list(W = env_kernel(env.data = EC.Dat,gaussian=FALSE)[[2]])


plotEnvRel <- function(KE){
    #dev.off()
	#### Heatmap of Raw Weather Data  ##key.par=list(mar=c(1,4,2,1))
	par(oma = c(1, 1, 1, 1), mar = c(5, 4, 4, 2) + 0.1) # Adjust margins to avoid the 'pin' issue
	Hmp_Plot <- gplots::heatmap.2(KE$W,dendrogram="row",Colv=NA,col=bluered(75),scale="none",margins=c(10,10),key=TRUE,keysize=0.8,key.title="Env Relationship",
			  symm=TRUE,key.par=list(mar=c(5,1,1,1)),cexRow=1.5,cexCol=1.5,trace='none')

    return(Hmp_Plot)

}


syncEnvPhenoDat <- function(KE,LocCoord,OtherLoc){

 ke <- KE
 # Create a named vector to map abbreviations to full names
  if(OtherLoc==TRUE){
	  level_mapping <- setNames(as.character(LocCoord[,"OtherLocName"]),as.character(LocCoord[,"Location"]))
	  colnames(ke$W) <- (level_mapping[colnames(ke$W)])
	  rownames(ke$W) <- (level_mapping[rownames(ke$W)])
  }
  
  return(ke)
 }





#######

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

 # outputListME <- outputListSTME()
 # TraitME <- TraitME()
 # IDColME <- IDColME()
 # IDColsME <- IDColsME()
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

# resultsTab <- getOutTab_ME_CV(result,traits)
# resultsTab3 <- getOutTab_ME_CV(result3,CVMet,traits)
 # resultsTab4 <- getOutTab_ME_CV(result4,CVMet,traits)
 # resultsTab5 <- getOutTab_ME_CV(result5,CVMet,traits)

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

###

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


#### #foreach(nTst=1:nTsts,.packages=c("BGGE","EnvRtype")) %dopar%

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

get_kernel_MDe <- function(Y,KG,KE,intercept.random=FALSE,dimension_KE=NULL){
  
  ne = length(unique(Y$env))
  Zg <- stats::model.matrix(~0 + gid, Y)
  colnames(Zg) = gsub(colnames(Zg), pattern = "gid", replacement = "")
  ng <- length(unique(Y$gid))
  
  K_G = KG 
  K_E = KE
  
  model_b <- "MDe"
  
 # Kde = getK(Y = Y, setKernel = K_G, model = model_b, intercept.random = intercept.random) 
  
  Kde = BGGE::getK(Y = Y, X= X,kernel = "GK", model = model_b, intercept.random = intercept.random) 
  
  
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

get_geno_kernel<- function(Y,X,KG,KE,intercept.random=FALSE,model,method,dimension_KE=NULL){
  
  ne = length(unique(Y$env))
  Zg <- stats::model.matrix(~0 + gid, Y)
  colnames(Zg) = gsub(colnames(Zg), pattern = "gid", replacement = "")
  ng <- length(unique(Y$gid))
  
  K_G = KG 
  K_E = KE
  
  model_b <- model
  
  if(is.null(KG)){
    if(method=="GB"){
      K= BGGE:getK(Y = Y,X=X, kernel = "GB", model = model_b, intercept.random = intercept.random) 
    }else if(method=="GK"){
      K = BGGE::getK(Y = Y, X= X,kernel = "GK", model = model_b, intercept.random = intercept.random) 
    }
  }else if(!is.null(KG)){
    
    K= BGGE::getK(Y = Y, setKernel = K_G, model = model_b, intercept.random = intercept.random) 
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

Vcomp.BGGE <- function(model, env, gid, digits = digits, 
                       alfa = 0.1){
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





 
get_weather_terra <- function (env.id = NULL, lat = NULL, lon = NULL, start.day = NULL, 
    end.day = NULL, variables.names = NULL, dir.path = NULL, 
    save = FALSE, temporal.scale = "daily", country = NULL, parallel = TRUE, 
    workers = NULL, chunk_size = 29, sleep = 60) 
{
    if (!requireNamespace("doParallel", quietly = TRUE)) {
        utils::install.packages("doParallel")
    }
    if (!requireNamespace("parallel", quietly = TRUE)) {
        utils::install.packages("doParallel")
    }
    if (!requireNamespace("foreach", quietly = TRUE)) {
        utils::install.packages("foreach")
    }
    get_helper <- function(lon, lat, variables.names, start.day, 
        end.day, env.id, save) {
        CL = data.frame(nasapower::get_power(community = "ag", 
            lonlat = c(lon, lat), pars = variables.names, dates = c(start.day, 
                end.day), temporal_api = "daily"))
        cor_rain_name = which(names(CL) %in% "PRECTOTCORR")
        names(CL)[cor_rain_name] = "PRECTOT"
        CL$daysFromStart = 1:nrow(CL)
        CL$env <- env.id
        CL <- CL[, c(which(colnames(CL) == "env"), which(colnames(CL) != 
            "env"))]
        if (isTRUE(save)) {
            utils::write.csv(file = paste(env.id, ".csv", sep = ""), 
                row.names = F, x = CL)
        }
        return(CL)
    }
   
   # sec_to_hms <- function(t) {
        # paste(formatC(t%/%(60 * 60)%%24, width = 2, format = "d", 
            # flag = "0"), formatC(t%/%60%%60, width = 2, format = "d", 
            # flag = "0"), formatC(t%%60, width = 2, format = "d", 
            # flag = "0"), sep = ":")
    # }
    # progress <- function(min = 0, max = 100, leftd = "|", rightd = "|", 
        # char = "=", style = 2, width = getOption("width"), time = Sys.time()) {
        # return(list(min = min, max = max, leftd = leftd, rightd = rightd, 
            # char = char, style = style, width = width, time = time))
    # }
    # run_progress <- function(pb, actual, text = "", digits = 0, 
        # sleep = 0) {
        # Sys.sleep(sleep)
        # elapsed <- sec_to_hms(as.numeric(difftime(Sys.time(), 
            # pb$time, units = "secs")))
        # temp <- switch(pb$style, list(extra = nchar(text) + nchar(pb$leftd) + 
            # nchar(pb$rightd), text = paste(text, paste(pb$leftd, 
            # "%s%s", pb$right, sep = ""))), list(extra = nchar(text) + 
            # nchar(pb$leftd) + nchar(pb$rightd) + 6, text = paste(text, 
            # paste(pb$leftd, "%s%s", pb$right, sep = ""), "% s%%")), 
            # list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 
                # 9, text = paste(text, paste(pb$leftd, "%s%s", 
                # pb$rightd, sep = ""), elapsed)), list(extra = nchar(text) + 
                # nchar(pb$leftd) + nchar(pb$rightd) + 15, text = paste(text, 
                # paste(pb$leftd, "%s%s", pb$rightd, sep = ""), 
                # "% s%%", elapsed)))
        # step <- round(actual/pb$max * (pb$width - temp$extra))
	    # # Check if step and the width calculation are valid
        # if (step < 0) step <- 0
        # if (pb$width - step - temp$extra < 0) {
           # warning("Progress bar width is too small to display properly")
           # pb$width <- step + temp$extra + 1
         # }
        # cat(sprintf(temp$text, strrep(pb$char, step), strrep(" ", 
            # pb$width - step - temp$extra), round(actual/pb$max * 
            # 100, digits = digits)), "\r")
        # if (actual >= pb$max) {
            # cat("\n")
        # }
    # }
	
	# Function to initialize progress bar settings
	progress <- function(min = 0, max = 100, leftd = "|", rightd = "|", 
						 char = "=", style = 2, width = getOption("width"), time = Sys.time()) {
	  if (max <= 0) {
		stop("max must be a positive number")
	  }
	  if (width <= 0) {
		stop("width must be a positive number")
	  }
	  return(list(min = min, max = max, leftd = leftd, rightd = rightd, 
				  char = char, style = style, width = width, time = time))
	}

	# Function to run and update the progress bar
	run_progress <- function(pb, actual, text = "", digits = 0, sleep = 0){
	  Sys.sleep(sleep)
	  elapsed <- sec_to_hms(as.numeric(difftime(Sys.time(), pb$time, units = "secs")))
	  
	  temp <- switch(pb$style,
					 list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd),
						  text = paste(text, paste(pb$leftd, "%s%s", pb$rightd, sep = ""))),
					 list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 6,
						  text = paste(text, paste(pb$leftd, "%s%s", pb$rightd, sep = ""), "%s%%")),
					 list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 9,
						  text = paste(text, paste(pb$leftd, "%s%s", pb$rightd, sep = ""), elapsed)),
					 list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 15,
						  text = paste(text, paste(pb$leftd, "%s%s", pb$rightd, sep = ""), "%s%%", elapsed))
	  )
	  
	  step <- round(actual / pb$max * (pb$width - temp$extra))
	  
	  # Ensure step and width calculation are valid
	  if (step < 0) step <- 0
	  if (pb$width - step - temp$extra < 0) {
		# Adjust pb$width to a minimum valid size
		pb$width <- step + temp$extra + 1
	  }
	  
	  # Print the progress bar
	  cat(sprintf(temp$text, strrep(pb$char, step), strrep(" ", pb$width - step - temp$extra),
				  round(actual / pb$max * 100, digits = digits)), "\r")
	  
	  if (actual >= pb$max) {
		cat("\n")
	  }
	}

	# Helper function to convert seconds to HH:MM:SS format
	sec_to_hms <- function(seconds) {
	  h <- floor(seconds / 3600)
	  m <- floor((seconds %% 3600) / 60)
	  s <- seconds %% 60
	  sprintf("%02f:%02f:%02f", h, m, s)
	}


	
	
    split_chunk <- function(vec, length) {
        split(vec, ceiling(seq_along(vec)/length))
    }
    cat("------------------------------------------------ \n")
    cat("ATTENTION: This function requires internet access \n")
    cat("------------------------------------------------  \n")
    cat("Connecting to the NASA POWER API Client, Sparks et al 2018 \n")
    cat("https://docs.ropensci.org/nasapower \n")
    cat("------------------------------------------------  \n")
    if (is.null(env.id)) {
        env.id <- paste0("env", seq_along(lat))
    }
    if (!(is.character(env.id) || is.factor(env.id))) {
        stop("env.id should be a vector of characters (e.g. 'env1') or factors")
    }
    if (!requireNamespace("nasapower", quietly = TRUE)) {
        utils::install.packages("nasapower")
    }
    if (!requireNamespace("plyr", quietly = TRUE)) {
        utils::install.packages("plyr")
    }
    if (is.null(dir.path)) {
        dir.path = getwd()
    }
    if (is.null(start.day)) {
        start.day <- Sys.Date() - 1000
        cat(paste0("start.day is NULL", "\n"))
        cat(paste0("matched as ", start.day, "\n"))
    }
    if (is.null(end.day)) {
        end.day <- start.day + 30
        cat(paste0("end.day is NULL", "\n"))
        cat(paste0("matched as ", end.day, "\n"))
    }
    if (is.null(variables.names)) {
        variables.names = c("T2M", "T2M_MAX", "T2M_MIN", "PRECTOTCORR", 
            "WS2M", "RH2M", "T2MDEW", "ALLSKY_SFC_LW_DWN", "ALLSKY_SFC_SW_DWN")
    }
    variables.names[grepl(variables.names, pattern = "PRECTOT")] = "PRECTOTCORR"
    env.id = as.factor(env.id)
    if (parallel == FALSE) {
        results <- list()
        pb <- progress(max = length(env.id), style = 4)
        init_time <- Sys.time()
        iter <- 0
        for (i in 1:length(env.id)) {
            iter <- iter + 1
            query_issue <- as.numeric(difftime(Sys.time(), init_time, 
                units = "secs")) > 60
            if (iter >= 30 & query_issue) {
                message("Waiting ", sleep, "s for a new query to the API.")
                Sys.sleep(sleep)
                iter <- 0
                init_time <- Sys.time()
            }
            results[[i]] <- get_helper(lon = lon[i], lat = lat[i], 
                variables.names = variables.names, start.day = start.day[i], 
                end.day = end.day[i], env.id = env.id[i], save = save)
            msg <- paste0("Env ", env.id[i], " (", i, "/", length(env.id), 
                ") ", "downloaded")
            run_progress(pb, actual = i, text = msg)
        }
		cat("\nNASA POWER: Done!")
    }
    if (parallel == TRUE) {
        env.id_par = split_chunk(env.id, length = chunk_size)
        lat_par = split_chunk(lat, length = chunk_size)
        lon_par = split_chunk(lon, length = chunk_size)
        start.day_par = split_chunk(start.day, length = chunk_size)
        end.day_par = split_chunk(end.day, length = chunk_size)
        nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 
            0.9), workers)
        clust <- parallel::makeCluster(nworkers)
        on.exit(parallel::stopCluster(clust))
        results <- list()
        pb <- progress(max = length(env.id_par), style = 4)
        for (i in 1:length(env.id_par)) {
            env.id_par_tmp <- env.id_par[[i]]
            lat_par_tmp <- lat_par[[i]]
            lon_par_tmp <- lon_par[[i]]
            start.day_par_tmp <- start.day_par[[i]]
            end.day_par_tmp <- end.day_par[[i]]
            parallel::clusterExport(clust, varlist = c("get_helper", 
                "lat_par_tmp", "lon_par_tmp", "variables.names", 
                "start.day_par_tmp", "end.day_par_tmp", "env.id_par_tmp"), 
                envir = environment())
            length_chunk <- length(env.id_par[[i]])
            temp <- parallel::parLapply(clust, 1:length_chunk, 
                function(j) {
                  get_helper(lon = lon_par_tmp[j], lat = lat_par_tmp[j], 
                    variables.names = variables.names, start.day = start.day_par_tmp[j], 
                    end.day = end.day_par_tmp[j], env.id = env.id_par_tmp[j], 
                    save = save)
                })
            results[[i]] <- plyr::ldply(temp)
            if (i < length(env.id_par)) {
                message("Waiting ", sleep, "s for a new query to the API.")
                msg <- paste0("Chunk ", i, "/", length(env.id_par), 
                  " (", length_chunk, " points) downloaded")
                run_progress(pb, actual = i, text = msg)
                Sys.sleep(sleep)
            }
        }
        cat("\nNASA POWER: Done!")
    }
   
    df <- plyr::ldply(results)
    if (is.null(country)) {
        cat("\n")
        cat("country = NULL ; Not possible to download ALT data. Please provide a country ID (e.g., USA1, BRA, FRA)")
        variables.names[grepl(variables.names, pattern = "PRECTOTCORR")] = "PRECTOT"
        ids = which(names(df) %in% variables.names)
        df[, ids][df[, ids] == -999] = NA
        return(df)
    }
    if (!is.null(country)) {
        if (!requireNamespace("geodata")) 
            utils::install.packages("geodata")
        suppressMessages("geodata")
        suppressWarnings("geodata")
        cat("\n")
        cat("Connecting to https://geodata.ucdavis.edu/ using Hijmans 2021")
        cat("------------------------------------------------  \n")
        unique_country <- unique(country)
        geodata_alt <- vector(mode = "list", length = length(unique_country))
        names(geodata_alt) = unique_country
        for (n in 1:length(unique_country)) {
            if (isTRUE(grepl(x = unique_country[n], pattern = "USA"))) {
                id_region = as.numeric(gsub(x = unique_country[n], 
                  pattern = "USA", replacement = ""))
				geodata_alt[[n]] <- suppressMessages(geodata::elevation_30s("USA",getwd(),mask=TRUE)[[id_region]])
            }
            else {
                geodata_alt[[n]] <- suppressMessages(geodata::elevation_30s(country = unique_country[n], getwd(),mask = TRUE))
            }
        }
        df_no_alt <- df
        df <- c()
        for (i in 1:length(env.id)) {
		
            country_raster <- geodata_alt[[which(names(geodata_alt) %in% 
                country[i])]]
            id <- which(df_no_alt$env %in% env.id[i])
            df <- rbind(df, extract_GIS_terra(covraster = country_raster, 
                env.data = df_no_alt[id, ], name.out = "ALT"))
			
        }
		colnames(df) <- gsub("ALT.*","ALT",colnames(df))
        variables.names <- c(variables.names, "ALT")
        cat("\nSRTM: Done!")
        variables.names[grepl(variables.names, pattern = "PRECTOTCORR")] = "PRECTOT"
        ids = which(names(df) %in% variables.names)
        df[, ids][df[, ids] == -999] = NA
        return(df)
    }
}

####

extract_GIS_terra <- function(covraster = NULL, Latitude = NULL, Longitude = NULL, 
                              env.data = NULL, env.id = NULL, name.out = NULL) 
{
    # Setting default names if not provided
    if (is.null(name.out)) 
        name.out = "ALT"
    if (is.null(Latitude)) 
        Latitude <- "LAT"
    if (is.null(Longitude)) 
        Longitude <- "LON"
    if (is.null(env.id)) 
        env.id <- "env"
    
    # Creating a vector data frame from the environmental data for latitude and longitude
    loc <- data.frame(x = env.data[, Longitude], y = env.data[, Latitude], id = env.data[, env.id])
    points <- terra::vect(loc, geom=c("x", "y"), crs="+proj=longlat +datum=WGS84")
    
    # Extracting data from each raster in the list
    extracted_data <- lapply(covraster, function(rast) terra::extract(rast,points, fun=mean, na.rm=TRUE, df=TRUE))
    extracted_data_df <- do.call(rbind, extracted_data)
    names(extracted_data_df)[-1] <- sapply(names(covraster), function(nm) paste(name.out, nm, sep="_"))

    # Ensuring the ID column from points is included for merging
    extracted_data_df$id <- points$id
    
    # Merging extracted data with the environmental data using the unique identifier
    final_data <- merge(env.data, extracted_data_df, by.x = env.id, by.y = "id")
	
	idCol <- which(colnames(final_data) %in% "ID")
	if(length(idCol)>0) { final_data_V2 <- final_data[,-idCol]}else{final_data_V2 <- final_data} 
    
    return(final_data_V2)
}






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



numMAF <- function(Geno,MAFTh){ 
   
   AF_Tab <- getFreq_Alleles(Geno)
   MAFVec <-  apply(AF_Tab,2,function(x) min(as.numeric(x)))
   numMAF_LT <- length(which(MAFVec <= MAFTh))
   return(numMAF_LT)
}

numMissSites <- function(GenoT,missSitesTH){ 

     Geno <- GenoT[,-c(1:5)]
     nInd <- ncol(Geno)
	 nSites <- nrow(Geno)
	 
	 NALen <- apply(Geno,1,function(x) length(which(is.na(x))))
	 numMissingSites <- length(which(NALen >= (missSitesTH*nInd)))
	 
	 return(numMissingSites)
 
 }
 
 
 
numMissInd <- function(Geno,missIndTH){ 

  nInd <- ncol(Geno)
  nSites <- nrow(Geno)
 
  NALen <- as.vector(apply(Geno,2,function(x) length(which(is.na(x)))))
  numMissingInd <- length(which(NALen >= missIndTH*nSites))
 
 return(numMissingInd)
 
} 

getGenoDiff <- function(Geno1,Geno2){ 

  DiffInd <- abs(ncol(Geno1)-ncol(Geno2)) 
  DiffSites <- abs(nrow(Geno1)-nrow(Geno2))
  
  return(list(DiffSites,DiffInd))
  
 }


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
   
