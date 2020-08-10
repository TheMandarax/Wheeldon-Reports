calc_enrich <- function(df, output, mem, sol) {
  output <- enquo(output)
  mem <- enquo(mem)
  sol <- enquo(sol)
  
  es_name <- quo_name(output)
  
  df %>% mutate(!! es_name := log2(((!! mem) / 1e6) / ((!! sol) / 1e6)))
}

calc_expression <- function(df, output, input, norm = NULL) {
  output <- enquo(output)
  input <- enquo(input)
  norm <- enquo(norm)
  
  output <- quo_name(output)
  
  if (grepl("tpm", output)) {
    df %>% 
      mutate(rpk = (!! input)/((!! norm)/1e3)) %>% 
      ungroup() %>% 
      mutate(!! output := rpk/(sum(rpk)/1e6)) %>% 
      select(-rpk)
  } else if (grepl("rpm", output)) {
    df %>% 
      group_by(.id) %>% 
      mutate(num = sum(smoothed)) %>% 
      group_by(mask, add = TRUE) %>% 
      mutate(den = sum(smoothed),
             den = first(den)) %>% 
      ungroup() %>% 
      mutate(rScaleFactor = num/den,
             rpm = nbrReads * rScaleFactor) %>% 
      filter(mask == FALSE) %>% 
      mutate(rpm = rpm/(sum(rpm)/1e6)) %>% 
      group_by(.id) %>% 
      mutate(rpm = sum(rpm)) %>% 
      select(-num, -den, -rScaleFactor)
  }
}

naTozeroRle <- function(rleObject) {
    ixNA <- which(is.na(S4Vectors::runValue(rleObject)))
  
    if(length(ixNA) > 0) {
        S4Vectors::runValue(rleObject)[ixNA] <- 0
        S4Vectors::runLength(rleObject)[ixNA] <- 0
    }
    
    return(rleObject)
}

applyShiftFeature <- function(transcGRangesList, shiftValue) {
    # Check parameter validity
    # If missing shift value or if shift value does not inherit from numeric class
    if(missing(shiftValue) || !inherits(shiftValue, "numeric")) {
        shiftValue <- 0
        warning("Incorrect shiftValue parameter! No shift is performed!\n")
    }
    
    # Check transcGRangesList if of class GRangesList
    if(!is(transcGRangesList, "GRangesList")) {
        stop(
            paste("transcGRangesList parameter is of class ", class(transcGRangesList), "instead of GRangesList!/n", sep=""))
    }
    
    # Determine transcript width from transcGRangesList
    transcWidth <- GenomicFeatures::transcriptWidths(start(transcGRangesList), end(transcGRangesList))
    
    # Takes absolute value of shift value
    absShiftVal <- abs(shiftValue)
    
    #if width of transcript is smaller than the absolute shiftValue eliminate the transcript
    ixSmallTransc <- which(transcWidth <= absShiftVal)
    if(length(ixSmallTransc) > 0) {
        transcGRangeList <- transGRangesList[-ixSmallTransc]
        transWidth <- GenomicFeatures::transcriptWidths(start(transcGRangesList), end(transcGRangesList))
    }
    
    # If the shiftValue is positive, the start of the transcript is shifted
    if(shiftValue > 0) {
        usefulRangeOnTransc <- cbind(startT = rep(absShiftVal + 1, length(transcGRangesList)), endT = transcWidth)
    }
    # Else it is the end of the transcript that we shift
    else {
        usefulRangeOnTransc <- cbind(startT=1, endT=transcWidth - absShiftVal)
    }
    
    # Make list of useful ranges using usefulRangeOnTrans (shifted values for transcript) across length of transcGRangesList
    listeUsefulRanges <- lapply(seq_len(length(transcGRangesList)), function(ixTransc){
        usefulRangeOnTransc[ixTransc, 1]:usefulRangeOnTransc[ixTransc, 2]
    })
    
    # For the remaining positions in the transcript, make 1bp bins of the genomic positions
    # TranscriptLocs2refLocs function for converting transcript-based locations in to reference- based locations.
    shiftedTransc <- GenomicFeatures::transcriptLocs2refLocs(listeUsefulRanges, start(transcGRangesList), end(transcGRangesList), as.character(S4Vectors::runValue(strand(transcGRangesList))), decreasing.rank.on.minus.strand=TRUE)
    
    # Give names for shiftedTransc the names from original transcGRangesList
    names(shiftedTransc) <- names(transcGRangesList)
    
    return(shiftedTransc)
}

# countShiftReads
countReads <- function (exonGRanges, cdsPosTransc, alnGRanges, originalAln, shiftValue, motifSize)  {
    
    if(missing(cdsPosTransc)) {
        stop("Missing cdsPosTransc parameter!\n")
    }
    
    if(length(exonGRanges) != length(cdsPosTransc)) {
        stop("Different lengths for exonGRanges and cdsPosTransc parameters!\n")
    }
    
    myCondNA <- which(is.na(unlist(cdsPosTransc)) | is.null(unlist(cdsPosTransc)))
    if(length(myCondNA) > 0) {
        stop("Non-null, non-NA values for the cdsPosTransc parameter!\n")
    }
    
    if(missing(shiftValue) || !inherits(shiftValue, "numeric")) {
        shiftValue <- 0
        warning("Incorrect shiftValue parameter! No shift is performed!\n")
    }
    
    if(!is(exonGRanges, "GRangesList")) {
        stop(paste("exonGRanges parameter is of class ", class(exonGRanges),
        " instead of GRangesList!\n", sep = ""))
    }
    
    if(!is(alnGRanges, "GRanges")) {
        stop(paste("alnGRanges parameter is of class ", class(alnGRanges),
        " instead of GRanges!\n", sep = ""))
    }
    
    if (missing(motifSize) || !is(motifSize, "numeric") || motifSize %% 1 != 0 || motifSize <= 0 || !(motifSize %in% c(3, 6, 9))) {
        warning("Param motifSize should be an integer! Accepted values 3, 6 or 9. Default value is 3.\n")
        motifSize <- 3
    }
    
    exonGRangesRestrict <- exonGRanges[names(cdsPosTransc)]
    if(length(exonGRangesRestrict) <= 5) {
        stop("Less than 5 common transcripts btw exonGRanges and cdsPosTransc!\n")
    }
    else {
        if (length(exonGRangesRestrict) <= 10) {
            warning("Less than 10 common transcripts between exonGRanges and cdsPosTransc!\n")
        }
    }
    
    transcWidth <- GenomicFeatures::transcriptWidths(start(exonGRangesRestrict), end(exonGRangesRestrict))
    
    absShiftVal <- abs(shiftValue)
    
    ixSmallTransc <- which(transcWidth <= absShiftVal)
    
    if(length(ixSmallTransc) > 0) {
        transcBig <- exonGRangesRestrict[-ixSmallTransc]
        cdsPosTRanscBig <- cdsPosTransc[-ixSmallTransc]
    }
    else {
        transcBig <- exonGRangesRestrict
    }
    
    overlapReads <- suppressWarnings(findOverlaps(originalAln, transcBig))
    
    startOverlapReads <- split(start(alnGRanges[queryHits(overlapReads)]), factor(subjectHits(overlapReads)))
    
    overlapReadsRle <- sapply(startOverlapReads, S4Vectors::Rle)
    
    transcWithReads <- transcBig[as.numeric(names(overlapReadsRle))]
    
    cdsPosTranscWithReads <- cdsPosTransc[as.numeric(names(overlapReadsRle))]
    
    cdslengthwithreads <- lapply(seq_len(NROW(cdsPosTranscWithReads)), function(ixTransc) {
        cdsPosTranscWithReads[[ixTransc]][2] - cdsPosTranscWithReads[[ixTransc]][1] + 1
    })
    
    newTranscWidth <- GenomicFeatures::transcriptWidths(start(transcWithReads), end(transcWithReads))
    
    cdsPosTranscShifted <- do.call(rbind, cdsPosTranscWithReads) + shiftValue
    
    listeRangesCDS <- lapply(seq_len(NROW(cdsPosTranscShifted)),
    
    function(ixTransc) {
        max(1, cdsPosTranscShifted[ixTransc, 1]):cdsPosTranscShifted[ixTransc, 2]
    })
    
    listeRanges5UTR <- lapply(seq_len(NROW(cdsPosTranscShifted)), function(ixTransc) {
        if((cdsPosTranscShifted[ixTransc, 1] - 1) < 1) {
            0
        }
        else {
            max(1, shiftValue):(cdsPosTranscShifted[ixTransc,
            1] - 1)
        }
    })
    
    listeRanges3UTR <- lapply(seq_len(NROW(cdsPosTranscShifted)),function(ixTransc) {
        (cdsPosTranscShifted[ixTransc, 2] + 1):min(newTranscWidth[ixTransc], newTranscWidth[ixTransc] + shiftValue)
    })
    
    binTransc <- applyShiftFeature(transcWithReads, 0)
    
    strandInfo <- S4Vectors::runValue(strand(transcWithReads))
    
    shiftedTranscMatches <- lapply(seq_len(length(binTransc)), function(ixTransc) {
        
        if(strandInfo[[ixTransc]] == "-") {
            binTranscVal <- sort(binTransc[[ixTransc]],
            decreasing = TRUE)
            txtail <- tail(binTranscVal, n = 1)
            binTranscVal <- c(binTranscVal, seq(txtail-1, txtail-25))
        }
        else {
            binTranscVal <- sort(binTransc[[ixTransc]])
            txtail <- tail(binTranscVal, n = 1)
            binTranscVal <- c(binTranscVal, seq(txtail+1, txtail+25))
        }
        
        matchedReadsTransc <- match(sort(overlapReadsRle[[ixTransc]]), binTranscVal)
        
        matchedReadsCDS <- naTozeroRle(match(matchedReadsTransc, listeRangesCDS[[ixTransc]]))
        
        matchedReads5UTR <- naTozeroRle(match(matchedReadsTransc, listeRanges5UTR[[ixTransc]]))
        
        matchedReads3UTR <- naTozeroRle(match(matchedReadsTransc, listeRanges3UTR[[ixTransc]]))
        
        if(length(matchedReadsCDS) > 0) {
            allCodonCounts <- aggregate(S4Vectors::runLength(matchedReadsCDS),
            by = list(ceiling(S4Vectors::runValue(matchedReadsCDS)/3)), FUN = sum)
            
            if(motifSize <= 3) {
                myCodonCounts <- allCodonCounts
            }
            else {
                if(motifSize == 6) {
                    myCodonCounts <- allCodonCounts[1:(nrow(allCodonCounts) - 1), ]
                }
                else {
                    if(motifSize == 9) {
                        myCodonCounts <- allCodonCounts[2:(nrow(allCodonCounts) - 1), ]
                    }
                }
            }
        }
        else {
            nbrCodons <- ceiling(length(listeRangesCDS[[ixTransc]])/motifSize)
            myCodonCounts <- data.frame(cbind(1:nbrCodons, rep(0, nbrCodons)))
        }
        
        nbrCodons <- ceiling(cdslengthwithreads[[ixTransc]]/motifSize)
        
        myCodonCounts2 <- data.frame(cbind(1:nbrCodons, rep(0, nbrCodons)))

        names(myCodonCounts) <- c("codonID", "nbrReads")
        names(myCodonCounts2) <- c("codonID", "nbrReads2")
        
        myCodonCounts3 <- merge.data.frame(myCodonCounts, myCodonCounts2, by = 'codonID', all = T)[,-3]
        
        names(myCodonCounts3) <- c("codonID", "nbrReads")
        
        list(c(sum(S4Vectors::runLength(matchedReadsCDS)), sum(S4Vectors::runLength(matchedReads5UTR)), sum(S4Vectors::runLength(matchedReads3UTR))), myCodonCounts3)
    })
    
    names(shiftedTranscMatches) <- names(transcWithReads)
    
    countsFeatures <- do.call(rbind, lapply(shiftedTranscMatches, `[[`, 1))
    
    colnames(countsFeatures) <- c("CDS_counts", "fiveUTR_counts", "threeUTR_counts")
    rownames(countsFeatures) <- names(shiftedTranscMatches)
    
    chrInfo <- S4Vectors::runValue(seqnames(transcWithReads))
    startInfo <- min(start(transcWithReads))
    endInfo <- max(end(transcWithReads))
    cdsInfo <- do.call(rbind, cdsPosTranscWithReads)
    cdsLength <- cdsInfo[, 2] - cdsInfo[, 1] + 1
    cdsStart <- cdsInfo[, 1]
    cdsEnd <- cdsInfo[, 2]
    
    countsData <- cbind(as.character(rownames(countsFeatures)), as.character(unlist(chrInfo)), as.character(unlist(strandInfo)), startInfo, endInfo, newTranscWidth, cdsStart, cdsEnd, cdsLength, countsFeatures)
    
    colnames(countsData) <- c("gene", "chr", "strand", "transc_genomic_start", "transc_genomic_end", "transc_length", "orf_start", "orf_end", "orf_length", colnames(countsFeatures))
    
    codonReadCoverage <- lapply(shiftedTranscMatches, `[[`, 2)
    names(codonReadCoverage) <- names(shiftedTranscMatches)
    
    return(list(as.data.frame(countsData), codonReadCoverage))
    
}

plotSummarizedCov <- function (covSummarized) 
{
  if (!inherits(covSummarized, "list")) {
    stop("The covSummarized object is not a list!\n")
  }
  else {
    if (!is(covSummarized[[1]], "GRanges")) {
      stop("The covSummarized object is not a list of GRanges objects!\n")
    }
  }
  listPlotSum <- lapply(covSummarized, function(iSumCov) {
    maxPeak <- max(iSumCov$values)
    maxPeakPos <- start(iSumCov)[which(iSumCov$values == 
                                         maxPeak)]
    if (maxPeak <= 100) {
      yLab <- "% of reads"
    }
    else {
      yLab <- "Number of Reads"
    }
    iPlot <- ggplot(iSumCov, ggplot2::aes(start, values)) + 
      geom_point(color = darkblue1) + geom_line(color = blue) + 
      xlab("Distance from start codon (nt)") + xlim(0, 30) +
      ylab(yLab) +
      scale_y_continuous(labels = comma_format()) +
      paper_theme
    return(iPlot)
  })
  
  return(listPlotSum)
}

