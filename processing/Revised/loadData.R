load_files <- function() {

  cat("\nHow many files would you like to analyze?\n")
  
  num.files <- suppressWarnings(as.integer(readline()))
  
  if (is.na(num.files)) {
    cat("\nPlease provide answer as an integer!\n\n")
  } else {
    files <- c()
    sample_names <- c()
    for (i in 1:num.files) {
      cat(paste0("\nPlease input directory for file # ", i, "\n"))
      files[i] <- readline("Choice: ")
      
      cat("\nWhat name would you like to use for this sample?\n")
      sample_names[i] <- readline("Choice: ")
    }
    cat("\n")
    
    assign("files", files, envir = .GlobalEnv)
    assign("sample_names", sample_names, envir = .GlobalEnv)
  }
}

load_gff <- function() {
  
  cat("\nPlease provide directory to GFF file\n")
  annotation.gff <- readline("Choice: ")
  
  txdb <- GenomicFeatures::makeTxDbFromGFF(
    annotation.gff, format = "gff3", circ_seqs = character())
  txids <- AnnotationDbi::keys(txdb, keytype="TXNAME")
  cdsTransc <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = T)
  exonGRanges <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = T)
  cdsPosTransc <- RiboProfiling::orfRelativePos(cdsTransc, exonGRanges)
  
  tx <- list(
    "txdb" = txdb,
    "txids" = txids, 
    "cdsTransc" = cdsTransc,
    "exonGRanges" = exonGRanges,
    "cdsPosTransc" = cdsPosTransc)
  
  assign("tx", tx, envir = .GlobalEnv)
}