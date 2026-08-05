
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
#' @export
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
#' @export
getImputedData_Num <- function(FiltGeno,nN,Dist){
  
  
   tasGenoImp <- rTASSEL::imputeNumeric(FiltGeno,byMean = TRUE,nearestNeighbors = nN, distance = Dist)
  
  return(tasGenoImp)
}


#' Paths of the AlphaPlantImpute intermediate files
#'
#' The genotype table, marker table, haplotype library and imputed output that
#' the AlphaPlantImpute step passes between its stages. These were formerly bare
#' filenames, i.e. resolved against the R working directory, which
#' \code{shiny::runApp()} pins to the app directory for the life of the session.
#' Under Docker/Singularity that directory lives in the container's writable
#' layer, so the files never appear in a bind-mounted results folder and are
#' discarded when the container is removed — taking \code{lib.phase} (expensive
#' to rebuild) with them. Addressing them relative to the run's output directory
#' keeps them with the rest of the results.
#'
#' @param outDir Directory the files live in. Defaults to the working directory,
#'   reproducing the historical behaviour.
#'
#' @return A named list of paths: \code{dir}, \code{genoFile}, \code{idTabFile},
#'   \code{libPrefix}, \code{libFile}, \code{impPrefix}, \code{impFile}.
#'   The two \code{*Prefix} entries are for AlphaPlantImpute's \code{-out}
#'   argument, which it interpolates directly into the output filename.
#' @examples
#' apiFilePaths(tempdir())
#' @export
apiFilePaths <- function(outDir = ".") {
  if (is.null(outDir) || !nzchar(trimws(outDir))) outDir <- "."
  list(
    dir       = outDir,
    genoFile  = file.path(outDir, "current_GenoTable.genotypes"),
    idTabFile = file.path(outDir, "currentVCFIDTab.txt"),
    libPrefix = file.path(outDir, "lib"),
    libFile   = file.path(outDir, "lib.phase"),
    impPrefix = file.path(outDir, "imputed_out"),
    impFile   = file.path(outDir, "imputed_out.genotypes")
  )
}

#' Run a fallible step from inside a Shiny observer
#'
#' An unhandled error in \code{observe()}/\code{observeEvent()} ends the whole
#' session — the user loses every dataset loaded so far, not just the step that
#' failed. This turns it into a notification and returns \code{NULL}, so the
#' observer stops there and everything already loaded survives.
#'
#' \code{req()} and \code{validate()} raise \code{shiny.silent.error} as
#' ordinary control flow ("not ready yet"), so that class is re-raised
#' untouched — catching it would turn every skipped observer into a spurious
#' error notification.
#'
#' @param what Human-readable name of the step, used in the notification.
#' @param expr Expression to evaluate. Lazily evaluated inside the handler.
#'
#' @return The value of \code{expr}, or \code{NULL} if it raised an error.
#' @examples
#' \dontrun{
#' observeEvent(input$RunGmat, withStepError("G matrix estimation", {
#'   kinGmatResult(getGmat(kinQCResult()$M.clean))
#' }))
#' }
#' @export
withStepError <- function(what, expr) {
  ## Many steps write their status through sink(temp_file); if the step throws,
  ## the matching sink() never runs and EVERY later print() in the session is
  ## silently diverted into that temp file. Restore the depth we started at.
  sinkDepth <- sink.number()

  ## One handler, and the decision made after tryCatch has returned. Re-raising
  ## a silent error from inside a `shiny.silent.error =` handler does NOT skip
  ## the sibling `error =` handler of the same tryCatch — it runs too, and a
  ## plain req() ends up reported as a step failure. Returning a tagged result
  ## and re-signalling out here keeps the two paths genuinely separate.
  res <- tryCatch(list(ok = TRUE, value = expr),
                  error = function(e) list(ok = FALSE, cond = e))

  while (sink.number() > sinkDepth) sink()

  if (isTRUE(res$ok)) return(res$value)
  ## req()/validate() signal shiny.silent.error as ordinary "not ready yet"
  ## control flow. Pass it on untouched: Shiny stops the observer quietly.
  if (inherits(res$cond, "shiny.silent.error")) stop(res$cond)

  shiny::showNotification(paste0(what, " failed: ", conditionMessage(res$cond)),
                          type = "error", duration = NULL)
  NULL
}

#' Banner for a background AlphaPlantImpute step that died
#'
#' A finished process is not a successful one — callr reports a non-zero exit
#' status when the worker failed. The useful explanation is nearly always the
#' last few lines of the run log the caller already has, so this only adds the
#' banner and decodes the one error shape that actively misleads: reticulate
#' re-raises AlphaPlantImpute's own \code{sys.exit()} as \code{"'SystemExit'
#' object does not support item assignment"}, which says nothing about the cause
#' printed just above it.
#'
#' @param what Name of the step, e.g. \code{"Imputation"}.
#' @param status Exit status from \code{$get_exit_status()}, or \code{NA}.
#' @param log Character vector of the worker's log lines, used to detect the
#'   \code{SystemExit} case.
#'
#' @return A character string to append to the status panel.
#' @examples
#' apiFailureNote("Imputation", 1L, "TypeError: 'SystemExit' object ...")
#' @export
apiFailureNote <- function(what, status = NA_integer_, log = NULL) {
  note <- paste0("=== ", what, " FAILED",
                 if (!is.na(status)) paste0(" (exit status ", status, ")") else "",
                 " ===")
  hint <- paste("The run log above stops at the point of failure;",
                "its last lines name the cause.")
  if (!is.null(log) && any(grepl("SystemExit", log, fixed = TRUE)))
    hint <- paste(hint,
                  paste("AlphaPlantImpute stopped itself; reticulate reports that as",
                        "\"'SystemExit' object does not support item assignment\".",
                        "Ignore that wording and read the ERROR line above it."),
                  sep = "\n")
  paste(note, hint, sep = "\n")
}

#' Title of getGenoData_API
#'
#' Description of what getGenoData_API does.
#'
#' @param FiltGeno Description of FiltGeno.
#' @param outDir Directory to write the genotype and marker tables to; see
#'   \code{\link{apiFilePaths}}.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getGenoData_API
#' result <- getGenoData_API(...)
#' }
#' @export
getGenoData_API <- function(FiltGeno, outDir = "."){

    apiPaths <- apiFilePaths(outDir)
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
    write.table(vcfIDTab,apiPaths$idTabFile,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

    write.table(currentGenoMat,apiPaths$genoFile,quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)

}

#####

#' Title of getImpGenoData_API
#'
#' Description of what getImpGenoData_API does.
#'
#' @param vcfIDTab Description of vcfIDTab.
#' @param outDir Directory holding the imputed genotypes and the marker table;
#'   see \code{\link{apiFilePaths}}.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getImpGenoData_API
#' result <- getImpGenoData_API(...)
#' }
#' @export
getImpGenoData_API <- function(vcfIDTab, outDir = "."){

	apiPaths <- apiFilePaths(outDir)
	# vcfIDTab_DF <- as.data.frame(vcfIDTab)
	if (!file.exists(apiPaths$impFile))
	  stop("No imputed genotypes at '", apiPaths$impFile,
	       "'. Run the AlphaPlantImpute step first.", call. = FALSE)
	if (!file.exists(apiPaths$idTabFile))
	  stop("No marker table at '", apiPaths$idTabFile, "'.", call. = FALSE)

    genoImp <-  read.table(apiPaths$impFile,sep=" ",header=FALSE)
	vcfIDTab_DF <- as.data.frame(read.table(apiPaths$idTabFile,sep="\t",header=TRUE))

	## AlphaGenes .genotypes carries no marker names, so column i is row i of the
	## marker table and a count mismatch misattributes every genotype. This also
	## catches the case that actually bit — an imputed_out left behind by an
	## earlier run on a different VCF, which looks valid until the column names
	## are assigned.
	n_marker_cols <- ncol(genoImp) - 1L
	if (n_marker_cols != nrow(vcfIDTab_DF))
	  stop(sprintf(paste0(
	    "Imputed genotypes have %d marker columns but the marker table has %d ",
	    "rows, so they are not from the same run. This usually means '%s' is ",
	    "left over from an earlier imputation on a different marker set. ",
	    "Re-run the AlphaPlantImpute build and impute steps to regenerate it."),
	    n_marker_cols, nrow(vcfIDTab_DF), apiPaths$impFile), call. = FALSE)

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
