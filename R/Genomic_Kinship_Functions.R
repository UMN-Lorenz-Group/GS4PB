
#' QC-filter a marker matrix using ASRgenomics::qc.filtering
#'
#' @param genoDF Data frame. Genotypic data frame produced by getGenoTas_to_DF.
#' @param maf Numeric. Minor allele frequency threshold (default 0.02).
#' @param marker.callrate Numeric. Marker call rate threshold (default 0.2).
#' @param ind.callrate Numeric. Individual call rate threshold (default 0.20).
#' @param heterozygosity Numeric. Heterozygosity threshold (default 0.7).
#' @param Fis Numeric. Inbreeding coefficient threshold (default 1).
#' @param impute Logical. Whether to impute missing values (default TRUE).
#'
#' @return A list with elements: M.clean, msg, plot.heteroz, plot.maf,
#'   plot.missing.ind, plot.missing.SNP, plot.Fis.
#' @examples
#' \dontrun{
#' # Example usage of getQCFilteredM
#' result <- getQCFilteredM(genoDF, maf = 0.05)
#' }
#' @export
getQCFilteredM <- function(genoDF, maf=0.02, marker.callrate=0.2,
                           ind.callrate=0.20, heterozygosity=0.7,
                           Fis=1, impute=TRUE) {
  M_raw <- .genoDF_to_M(genoDF)
  # Capture both stdout (cat) and stderr (message) — qc.filtering uses both
  msg_out <- character(0)
  msg_err <- character(0)
  result <- NULL
  msg_out <- capture.output({
    msg_err <- capture.output({
      result <- ASRgenomics::qc.filtering(
        M = M_raw, base = FALSE, ref = NULL,
        maf = maf, marker.callrate = marker.callrate,
        ind.callrate = ind.callrate, heterozygosity = heterozygosity,
        Fis = Fis, impute = impute, na.string = "-9", plots = TRUE
      )
    }, type = "message")
  }, type = "output")
  combined_msg <- paste(c(msg_out, msg_err), collapse = "\n")
  list(
    M.clean        = result$M.clean,
    msg            = combined_msg,
    plot.heteroz   = result$plot.heteroz,
    plot.maf       = result$plot.maf,
    plot.missing.ind = result$plot.missing.ind,
    plot.missing.SNP = result$plot.missing.SNP,
    plot.Fis       = result$plot.Fis
  )
}

#' Compute genomic relationship matrix (G matrix)
#' @param M_clean Numeric matrix. QC-filtered marker matrix from getQCFilteredM.
#' @param method Character. Estimation method (default "VanRaden").
#' @return A list with elements: G (matrix), msg (character).
#' @export
getGmat <- function(M_clean, method = "VanRaden") {
  G <- NULL
  msg <- capture.output({
    G <- ASRgenomics::G.matrix(M = M_clean, method = method, na.string = NA)$G
  })
  list(G = G, msg = paste(msg, collapse = "\n"))
}

#' Kinship diagnostics and duplicate detection
#' @param G Numeric matrix. Genomic relationship matrix.
#' @param duplicate.thr Numeric. Threshold for calling duplicates (default 0.98).
#' @return ASRgenomics kinship.diagnostics list with added msg element.
#' @export
getKinshipDiag <- function(G, duplicate.thr = 0.98) {
  check <- NULL
  msg <- capture.output({
    check <- ASRgenomics::kinship.diagnostics(K = G, duplicate.thr = duplicate.thr)
  })
  check$msg <- paste(msg, collapse = "\n")
  check
}

#' Remove duplicate individuals from G matrix (greedy)
#' @param G Numeric matrix. Genomic relationship matrix.
#' @param list.duplicate Data frame. Duplicate pairs from getKinshipDiag.
#' @return List with: G_clean, removed, n_pairs, msg.
#' @export
removeDuplicatesFromG <- function(G, list.duplicate) {
  if (is.null(list.duplicate) || nrow(list.duplicate) == 0)
    return(list(G_clean = G, removed = character(0),
                n_pairs = 0L, msg = "No duplicate pairs detected — G matrix unchanged."))

  to_remove <- character(0)
  for (i in seq_len(nrow(list.duplicate))) {
    s1 <- as.character(list.duplicate[i, 1])
    s2 <- as.character(list.duplicate[i, 2])
    if (!s1 %in% to_remove) {
      # s1 survives — mark s2 for removal (if not already marked)
      to_remove <- unique(c(to_remove, s2))
    }
    # if s1 is already removed, s2 is the survivor — no action
  }

  keep    <- setdiff(rownames(G), to_remove)
  G_clean <- G[keep, keep, drop = FALSE]

  msg <- paste0(
    "Duplicate pairs detected: ", nrow(list.duplicate), "\n",
    "Individuals removed (one per pair): ", length(to_remove), "\n",
    "  Removed: ", paste(to_remove, collapse = ", "), "\n",
    "Individuals retained in G: ", length(keep)
  )
  list(G_clean = G_clean, removed = to_remove, n_pairs = nrow(list.duplicate), msg = msg)
}

#' Fine-tune G matrix for positive-definiteness (blend or bend)
#' @param G Numeric matrix. Genomic relationship matrix.
#' @param method Character. "blend" or "bend" (default "blend").
#' @param pblend Numeric. Blend proportion for identity matrix (default 0.05).
#' @return List with: Gb (tuned matrix), msg, check, check_msg.
#' @export
getTunedG <- function(G, method = "blend", pblend = 0.05) {
  Gb <- NULL
  msg <- capture.output({
    if (method == "blend") {
      Gb <- ASRgenomics::G.tuneup(G = G, blend = TRUE, pblend = pblend)$Gb
    } else {
      Gb <- ASRgenomics::G.tuneup(G = G, bend = TRUE)$Gb
    }
  })
  check <- NULL
  check_msg <- capture.output({
    check <- ASRgenomics::kinship.diagnostics(K = Gb)
  })
  list(
    Gb        = Gb,
    msg       = paste(msg, collapse = "\n"),
    check     = check,
    check_msg = paste(check_msg, collapse = "\n")
  )
}

#' PCA on genomic relationship matrix
#' @param G Numeric matrix. Genomic relationship matrix.
#' @param ncp Integer. Number of principal components to retain (default 20).
#' @return FactoMineR PCA object.
#' @export
getKinshipPCA <- function(G, ncp = 20) {
  ASRgenomics::kinship.pca(K = G, ncp = ncp)
}
