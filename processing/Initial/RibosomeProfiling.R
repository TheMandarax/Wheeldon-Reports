### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### Start Up
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
repeat {
  library(ggplot2)
  library(plotly)
  library(data.table)
  library(RiboProfiling)
  library(Biostrings)
  library(plyr)
  library(zoo)
  library(tools)
  library(methods)
  library(seqinr)
  library(clipr)
  library(rtracklayer)
  
  source("/Users/mandarax/Public/RStudio/scripts/chartron/justinFunctions.R")
  source("/Users/mandarax/Public/RStudio/scripts/avlayortFunctions.R")
  source("/Users/mandarax/Public/RStudio/scripts/Ontology.R")
  
  setwd("/Users/mandarax/Public/PyCharm/Sequencing")
  break }

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### Acquire files
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
load_data()
load_names(); load_annotations(cell.line)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### Artifical sequence generator to produce mask files
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (exists("annotation.cds")) {
  dna.fasta <- Biostrings::readDNAStringSet(annotation.cds)
  seq.names <- names(dna.fasta)
  gene.seq <- paste(dna.fasta)
  dna.dt <- data.table(seq.names, gene.seq)
}

# This will be changed depending on the way that Biostrings reads names
dna.dt[, id := sapply(seq.names, 
                      function(x) strpick(x, split = "|", pos = 1))]

dna.dt[, gene.name := sapply(seq.names, 
                             function(x) strpick(x, split = "|", pos = 2))]

dna.dt[, id := sub("_[0-9]*$", "", id)]

# Partition DNA sequences into 28mer fragments and generate DNA fasta file 
# representing pseudo Ribo-seq data. Use this dat to run Hisat2 on to get 
# sort.bam files needed to create mask
repeat {
  dna.dt[, seq.names := NULL]
  
  if ("gene.name" %in% names(dna.dt)) {
    setcolorder(dna.dt, c('id', 'gene.name', 'gene.seq'))
  } else { setcolorder(dna.dt, c('id', 'gene.seq'))}
  setkey(dna.dt, id)
  
  dna.dt[, gene.seq.mix := lapply(gene.seq, function(x) partition_sequence(x))]
  dsave(dnt.dt)
  file.txt <- paste0(save_directory, cell.line, "_mix.txt")
  write.table(dna.dt$gene.seq.mix, 
              file = file.txt, 
              row.names = FALSE, col.names = FALSE, sep = '\n')
  
  root <- "/Users/mandarax/Public/PyCharm/PhD/"
  
  output.sam <- unlist(strsplit(annotation.cds, "/"))
  output.sam[3:4] <- c("Sam", "mask.sam")
  output.sam <- paste(as.list(output.sam), collapse = "/")
  output.sam <- paste0(root, output.sam)
  
  annotation.file <- paste0(root, annotation.hisat)
  
  command <- paste0(
    "/Repository/miniconda3/envs/sequencing/bin/hisat2 ",
    "-x ", annotation.file, " ",
    "-U ", file.txt, " ",
    "-p 7 ",
    "-r ",
    "-S ", output.sam
  )
  system(command)
  break
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### ### Do Alignments
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
repeat {
  ## Original alignment using readGAlignments
  sing.aln <- lapply(sing.files, read_g_alignments)
  if (exists("mult.names")) {
    mult.aln <- lapply(mult.files, read_g_alignments)
  }
  
  ## Convert original alignments to 5'(start) or 3'(end) positions
  sing.alnGRanges <- lapply(
    sing.aln, RiboProfiling::readsToStartOrEnd, what = 'end')
  
  if (exists("mult.names")) {
    mult.alnGRanges <- lapply(
      mult.aln, RiboProfiling::readsToStartOrEnd, what = 'end')
  }
  break }

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### Determine p-site offsets
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
repeat {
  flank_width <- 30
  ## Calculate flank size around transcriptional start site
  sing.oneBinRanges <- lapply(
    seq(length(sing.names)), 
    function(x) RiboProfiling::aroundPromoter(
      txdb = tx$txdb,
      alnGRanges = sing.alnGRanges[[x]],
      percBestExpressed = 0.25,
      flankSize = flank_width)
  )
  
  ## Calculate shifts away from transcriptional start site per read length.
  if (!exists("mult.names")) {
    
    sing.matchLenDistr <- lapply(sing.aln, RiboProfiling::histMatchLength)
    
    sing.size_ranges <- list()
    for (i in 1) {
      sing.size_ranges[[i]] <- which(
        sing.matchLenDistr[[i]][[1]]$counts > 3000) + 14
    }
    
    sing.listPromoterCov <- lapply(
      seq(length(sing.names)), 
      function(x) RiboProfiling::readStartCov(
        alnGRanges = sing.alnGRanges[[x]],
        oneBinRanges = sing.oneBinRanges[[x]],
        matchSize = c(28,30),
        fixedInterval = c(-flank_width, flank_width),
        renameChr = "aroundTSS",
        charPerc = "sum")
    )
    
  } else { 
    sing.listPromoterCov <- lapply(
      seq(length(sing.names)),
      function(x) RiboProfiling::readStartCov(
        alnGRanges = sing.alnGRanges[[x]],
        oneBinRanges = sing.oneBinRanges[[x]],
        matchSize = "all",
        fixedInterval = c(-28, 28),
        renameChr = "aroundTSS",
        charPerc = "sum")
    )
  }
  
  ## Obtain shift values from listPromoterCov and put into list
  shift <- list()
  for (i in 1:length(sing.names)) {
    
    intnum <- which.max(
      sing.listPromoterCov[[i]]$sumUp@elementMetadata@listData[["values"]])
    
    shift[[i]] <- sing.listPromoterCov[[i]]$sumUp@ranges@start[intnum]
  }
break }

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### Get counts from alignments
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
repeat {
  sing.counts <- lapply(
    seq(length(sing.names)), 
    function(x) codonjustin(
      exonGRanges = tx$exonGRanges[names(tx$cdsPosTransc)],
      cdsPosTransc = tx$cdsPosTransc,
      alnGRanges = sing.alnGRanges[[x]],
      originalAln = sing.aln[[x]],
      shiftValue = as.numeric(shift[[x]]),
      motifSize = 3)
  )
  
  if (exists("mult.names")) {
    mult.counts <- lapply(
      seq(length(sing.names)), 
      function(x) codonjustin(
        exonGRanges = tx$exonGRanges[names(tx$cdsPosTransc)],
        cdsPosTransc = tx$cdsPosTransc,
        alnGRanges = mult.alnGRanges[[x]],
        originalAln = mult.aln[[x]],
        shiftValue = as.numeric(shift[[x]]),
        motifSize = 3)
    )
  }
  break }

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### Convert counts to individual sample data tables
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("mult.names")) {
  counts.dt <- list()
  for (i in 1:length(sing.names)) {
    ## Convert count list to individual data tables
    counts.dt[[i]] <- as.data.table(ldply(sing.counts[[i]][[2]]))
    
    ## Change data table variables names
    names(counts.dt[[i]]) <-
      c("gene", "codon", paste0(sing.names[[i]], ".codonReads"))
    
    ## Assign NEW variable, gene.codon
    #  Gene and codon number pasted
    counts.dt[[i]][, gene.codon := paste(gene,codon)]
    
    ## Set key for counts tables as gene.codon
    setkey(counts.dt[[i]], gene.codon)
  }
  names(counts.dt) <- sing.names
} else {
  ## Convert count list to individual data tables
  s.mask <- as.data.table(ldply(sing.counts[[1]][[2]]))
  m.mask <- as.data.table(ldply(mult.counts[[1]][[2]]))
  
  ## Change data table variable names
  names(s.mask) <- c("gene", "codon", "sing.codonReads")
  names(m.mask) <- c("gene", "codon", "mult.codonReads")
  
  ## Assign NEW variable, gene.codon
  #  Gene and codon number pasted together
  s.mask[, gene.codon := paste(gene, codon)]
  m.mask[, gene.codon := paste(gene, codon)]
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### Create MASK file: Merge single and multiple mapping data sets together
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (exists("mult.names")) {
  
  
  ## Merge data sets
  mask.dt <- merge(
    s.mask, m.mask, 
    by.x = "gene.codon", by.y = "gene.codon", all = TRUE)
  
  ## Ensure all gene information is merged between data sets
  #  Remove duplicate variables, and rename variables
  mask.dt[is.na(gene.x), gene.x := gene.y]
  mask.dt[is.na(codon.x), codon.x := codon.y]
  mask.dt[, c('gene.y', 'codon.y') := NULL]
  names(mask.dt)[1:3] <- c("gene.codon", "gene", "codon")
  
  ## Convert NA values to zero
  na_zero(mask.dt, c(4:5))
  
  ## Calculate read difference between sets
  mask.dt[, diff := mult.codonReads - sing.codonReads]
  
  ## Create NEW variable, theory.l
  #  Number of codons in each gene
  mask.dt[, theory.l := length(gene.codon), by = "gene"]
  
  ## Assign NEW variable, mask 
  #  TRUE and exlclude the first five and last five codons of the ORF
  mask.dt[codon <= 5 | theory.l - codon <= 5, mask := TRUE]
  
  ## Assign mask as TRUE 
  #  If there is no read difference between single and multiple assigned reads 
  mask.dt[diff != 0, mask := TRUE]
  
  ## Assign mask as FALSE
  #  If they are not true
  mask.dt[is.na(mask), mask := FALSE]
  
  ## Create NEW variable, pseudo.l
  #  Number of unmasked codons in each gene
  mask.dt[mask == FALSE, pseudo.l := length(gene.codon), by = "gene"]
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### Apply mask to get calculate read counts, RPK, and TPM scores
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("mult.names")) {
  load(paste0(
    "~/Public/RStudio/samples/", cell.line, "_files/mask/mask.dt.RData"))
  load(paste0(
    "~/Public/RStudio/samples/", cell.line, "_files/prot/prot.dt.RData"))
  
  tmp.dt <- list()
  pmil <- list()
  for (i in 1:length(sing.names)) {
    if (cell.line == "pp") {
      counts.dt[[i]][, id := gsub(".t01", "", gene)]
    } else if (cell.line == "sc") {
      counts.dt[[i]][, id := gsub("_mRNA", "", gene)]
    } else if (cell.line == "yl") {
      load("~/Public/RStudio/samples/yl_files/name.convert.dt.RData")
      counts.dt[[i]] <- merge(counts.dt[[i]], name.convert.dt, 
                              by.x = 'gene', by.y = 'gene', all = TRUE)
      names(counts.dt[[i]])[[5]] <- 'id'
    }
    
    ## Merge data sets with mask.dt and remove masked codons
    counts.dt[[i]] <- merge(
      mask.dt[mask == FALSE, .(gene.codon, gene, codon, pseudo.l)],
      counts.dt[[i]][, c(4,5,3)], by = 'gene.codon',
      all = TRUE)[!is.na(gene)]
    
    setcolorder(counts.dt[[i]], c(2,5,3,1,4,6))
    
    ## Convert NA values to zero in read count column
    na_zero(counts.dt[[i]], c(6L))
    
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Assign NEW variable, tmp.dt to filter to be used for metagene norm
    tmp.dt[[i]] <- copy(counts.dt[[i]])
    
    # Remove dubious reads from data sets
    tmp.dt[[i]] <- tmp.dt[[i]][tmp.dt[[i]]$id %in% prot.dt$id]
    
    cat(paste0(
      "\n", 
      round((nrow(counts.dt[[i]]) - nrow(tmp.dt[[i]])) / 
              nrow(counts.dt[[i]]) * 100, digits = 2), 
      "% of reads are dubious! Removing these reads...\n"))
    
    counts.dt[[i]] <- counts.dt[[i]][counts.dt[[i]]$id %in% prot.dt$id]
    
    ## Get total gene read counts by summing counts at each codon in each gene
    tmp.dt[[i]][, paste0(sing.names[[i]], ".reads") := sum(.SD),
                .SDcols = paste0(sing.names[[i]], ".codonReads"), by = gene]
    
    ## Truncate data set
    tmp.dt[[i]] <- unique(tmp.dt[[i]][, c(1, 2, 5, 7)])
    
    ## Calculate RPK scores for truncated data table
    tmp.dt[[i]][, paste0(sing.names[[i]], ".rpk") := .SD/(pseudo.l/1e3),
                .SDcols = paste0(sing.names[[i]], ".reads")]
    
    ## Calculate per million scaling factor to convert RPK values to TPM values
    pmil[[i]] <- sum(tmp.dt[[i]][, c(5L)])/1e6
    
    ## Use RPK and per million scaling factor to create TPM scores
    tmp.dt[[i]][, paste0(sing.names[[i]], ".tpm") := .SD/pmil[[i]],
                .SDcols = paste0(sing.names[[i]], ".rpk")]
  }
  rm(pmil)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### Merge sample sets together and calculate enrichment scores
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("mult.names")) {
  
  temp.dt <- Reduce(function(x, y) merge(x, y, all = TRUE), tmp.dt)
  
  if (cell.line == "pp") {
    if (exp.type == "triplicate") {
      temp.dt[, es.1 := log2((m1.tpm/1e6)/(s1.tpm/1e6))]
      temp.dt[, es.3 := log2((m3.tpm/1e6)/(s3.tpm/1e6))]
      temp.dt <- temp.dt[es.1 < 1 & es.3 < 1]
    } 
    else if (exp.type == "chx") {
      temp.dt[, es.chx.ts := log2((t.chx.tpm/1e6)/(s.chx.tpm/1e6))]
      temp.dt[, es.chx.ms := log2((m.chx.tpm/1e6)/(s.chx.tpm/1e6))]
      temp.dt[, es.neg.ts := log2((t.neg.tpm/1e6)/(s.neg.tpm/1e6))]
      temp.dt[, es.neg.ms := log2((m.neg.tpm/1e6)/(s.neg.tpm/1e6))]
      temp.dt <- 
        temp.dt[es.chx.ts < 1 & es.chx.ms < 1 & es.neg.ts < 1 & es.neg.ms < 1]
    }
  } 
  
  else if (cell.line == 'sc') {
    if (exp.type == "triplicate") {
      temp.dt[, es.1 := log2((m1.tpm/1e6)/(s1.tpm/1e6))]
      temp.dt[, es.2 := log2((m2.tpm/1e6)/(s2.tpm/1e6))]
      temp.dt <- temp.dt[es.1 < 1 & es.2 < 1]
    } 
    else if (exp.type == "chx") {
      temp.dt[, es.chx.ms := log2((mp.tpm/1e6)/(sp.tpm/1e6))]
      temp.dt <- temp.dt[es.chx.ms < 1]
    }
  }
  
  else if (cell.line == 'yl') {
    temp.dt[, es.chx.ts := log2((t.chx.tpm/1e6)/(s.chx.tpm/1e6))]
    temp.dt[, es.chx.ms := log2((m.chx.tpm/1e6)/(s.chx.tpm/1e6))]
    temp.dt[, es.neg.ts := log2((t.neg.tpm/1e6)/(s.neg.tpm/1e6))]
    temp.dt[, es.neg.ms := log2((m.neg.tpm/1e6)/(s.neg.tpm/1e6))]
    temp.dt <- temp.dt[es.chx.ts < 1 & es.chx.ms < 1 & 
                         es.neg.ts < 1 & es.neg.ms < 1]
  }
  
  rm(tmp.dt)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### Perform metagene normalization to get normalized count values
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
repeat {
  start.roll <- 10  # Width of rolling mean in first window
  start.end <- 1000  # End of first rolling mean window (mid.start)
  
  mid.roll <- 40   # Width of rolling mean in second window
  mid.end <- 2000  # End of second rolling mean window
  
  end.roll <- 200
  end.end <- 3400
  
  med.start <- 3400  # Starting place for reference median 
  med.ext <- 10000   # Extension after midroll to base median on in third window
  

  smoothing.parameters <- list("start.roll"=start.roll, "start.end"=start.end,
                       "mid.roll"=mid.roll, "mid.end"=mid.end,
                       "end.roll"=end.roll, "end.end"=end.end,
                       "med.start"=med.start, "med.ext"=med.ext)
break } 

if (!exists("mult.names")) {
  meta.save <- paste0(save_directory, "meta/")
  meta <- list()
  meta.fig <- list()
  load_colors()
  fixed.dt <- list()
  for (i in 1:length(sing.names)) {
    ## Assign NEW variable, rpc.100
    #  Reads per codon over the first 100 codons for each gene
    counts.dt[[i]][codon <= 100, rpc.100 := sapply(.SD, mean),
                   .SDcols = paste0(sing.names[[i]], ".codonReads"),
                   by = gene]
    
    ## Assign codons with NA rpc.100 the rpc.100 associated with its gene
    counts.dt[[i]][, rpc.100 := unique(rpc.100[!is.na(rpc.100)]), 
                   by = c('gene', 'id', 'pseudo.l')]
    
    ## Change rpc.100 zero values to NA so that no infinite values are produced
    zero_na(counts.dt[[i]], 7L)
    
    ## Assign NEW variable, rpc.norm (Divide by rpc.100 to normalize codonReads) 
    counts.dt[[i]][, norm.100 := .SD/rpc.100,
                   .SDcols = paste0(sing.names[[i]], ".codonReads")]
    
    ## Assign NEW variable, rpcpg (Reads per codon per gene)
    counts.dt[[i]][, rpcpg := sapply(.SD, mean),
                   .SDcols = paste0(sing.names[[i]], ".codonReads"),
                   by = gene]
    
    ## Assign NEW variable, meta
    #  Not enriched data
    #  Mean of codon reads normalized by reads per codon value of 1st 100 codons
    #  here is where you will filter out membrane enriched genes
    meta[[i]] <- counts.dt[[i]][
      rpcpg > 0.5 & rpc.100 > 0 & gene %in% temp.dt[, gene],
      mean(norm.100, trim = 0.05, na.rm = TRUE),
      by = codon]
    
    ## Set key for meta set to codon
    setkey(meta[[i]], codon)
    
    ## Set name for meta set
    names(meta[[i]])[[2]] <- "rpc"
    
    ## Assign NEW variable, smoothed
    #  Rolling mean of mean normalized reads at each codon for all genes
    meta[[i]][codon >= 0 & codon <= smoothing.parameters$start.end,
              smoothed := rollapply(
                rpc, smoothing.parameters$start.roll, mean, partial = TRUE)]
    
    #  Extend the width of rolling mean after 1000 codons
    meta[[i]][codon > smoothing.parameters$start.end &
                codon < smoothing.parameters$mid.end,
              smoothed := rollapply(
                rpc, smoothing.parameters$mid.roll, mean, partial = TRUE)]
    
    #  Smooth longer codon reads by taking median of initial longer codons
    meta[[i]][codon > smoothing.parameters$mid.end &
                codon < smoothing.parameters$end.end,
              smoothed := rollapply(
                rpc, smoothing.parameters$end.roll, median, partial = TRUE)]
    
    meta[[i]][codon >= smoothing.parameters$med.start, 
              smoothed := meta[[i]][codon >= smoothing.parameters$med.start &
                                      codon <= smoothing.parameters$med.start +
                                      smoothing.parameters$med.ext,
                                    median(rpc)]]

    ## Assign NEW variable, sample.smoothedFIG
    #  Plot of the smoothed normalized read counts per codon
    meta.fig[[i]] <- ggplot(data = meta[[i]]) +
      # geom_point(aes(x = codon, y = rpc), color = grey, na.rm = TRUE) +
      geom_point(aes(x = codon, y = smoothed), color = darkblue, na.rm = TRUE) +
      labs(title = paste0(sing.names[[i]], " Smoothed Metagene Plot"),
           x = "Codon Position", y = "Smoothed Read Count") + 
      theme_classic()
    
    ## Assign NEW variable, sample.adjusted
    fixed.dt[[i]] <- merge(counts.dt[[i]], meta[[i]], by = "codon", all = TRUE)
    
    ## Convert NA values to zero for scaling calculation
    na_zero(fixed.dt[[i]], c(6L, 7L, 9L, 10L))
    
    ## Assign NEW column, sample.scaled
    #  Scales unmasked codon reads by dividing by smoothed codon counts
    fixed.dt[[i]][, paste0(sing.names[[i]], ".scaled") :=  .SD /smoothed, 
                  .SDcols = paste0(sing.names[[i]], ".codonReads")]
    
    ## Convert NaN values to NA for scaling calculation
    nan_na(fixed.dt[[i]], 12L)
    
    ## Convert Inf values to NA for scaling calculation
    inf_na(fixed.dt[[i]], 12L)
    
    ## Convert NA values to zero for scaling calculation
    na_zero(fixed.dt[[i]], 12L)
  }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### Count reads and develope TPM scores for adjusted data sets
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("mult.names")) {
  pmil <- list()
  for (i in 1:length(sing.names)) {
    ## Get total gene read counts by summing counts at each codon in each gene
    fixed.dt[[i]][, paste0(sing.names[[i]], ".reads") := sum(.SD),
                  .SDcols = paste0(sing.names[[i]], ".scaled"),
                  by = gene]
    
    ## Truncate data set
    fixed.dt[[i]] <- unique(fixed.dt[[i]][, c(3, 5, 13)])
    
    ## Calculate RPK scores for truncated data table
    fixed.dt[[i]][, paste0(sing.names[[i]], ".rpk") := .SD/(pseudo.l/1e3),
                  .SDcols = paste0(sing.names[[i]], ".reads")]
    
    ## Calculate per million scaling factor to convert RPK values to TPM values
    pmil[[i]] <- sum(fixed.dt[[i]][, c(4L)])/1e6
    
    ## Use RPK and per million scaling factor to create TPM scores
    fixed.dt[[i]][, paste0(sing.names[[i]], ".tpm") := .SD/pmil[[i]], 
                  .SDcols = paste0(sing.names[[i]], ".rpk")]
  }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### Merge normalized sample sets together and calculate enrichment scores
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("mult.names")) {
  if (cell.line == "pp") {
    if (exp.type == "triplicate") {
      pp.trip.dt <- Reduce(function(x, y) merge(x, y, all = TRUE), fixed.dt)
      pp.trip.dt[, es.1 := log2((m1.tpm/1e6)/(s1.tpm/1e6))]
      pp.trip.dt[, es.3 := log2((m3.tpm/1e6)/(s3.tpm/1e6))]
    } 
    else if (exp.type == "chx") {
      pp.chx.dt <- Reduce(function(x, y) merge(x, y, all = TRUE), fixed.dt)
      pp.chx.dt[, es.chx.ts := log2((t.chx.tpm/1e6)/(s.chx.tpm/1e6))]
      pp.chx.dt[, es.chx.ms := log2((m.chx.tpm/1e6)/(s.chx.tpm/1e6))]
      pp.chx.dt[, es.neg.ts := log2((t.neg.tpm/1e6)/(s.neg.tpm/1e6))]
      pp.chx.dt[, es.neg.ms := log2((m.neg.tpm/1e6)/(s.neg.tpm/1e6))]
    }
  } 
  else if (cell.line == "sc") {
    if (exp.type == "triplicate") {
      sc.trip.dt <- Reduce(function(x, y) merge(x, y, all = TRUE), fixed.dt)
      sc.trip.dt[, es.1 := log2((m1.tpm/1e6)/(s1.tpm/1e6))]
      sc.trip.dt[, es.2 := log2((m2.tpm/1e6)/(s2.tpm/1e6))]
    } else if (exp.type == "chx") {
      sc.chx.dt <- Reduce(function(x, y) merge(x, y, all = TRUE), fixed.dt)
      sc.chx.dt[, es.chx.ms := log2((mp.tpm/1e6)/(sp.tpm/1e6))]
    } else if (exp.type == "wei") {
      sc.wei.dt <- Reduce(function(x, y) merge(x, y, all = TRUE), fixed.dt)
    }
  }
  else if (cell.line == "yl") {
    if (exp.type == "chx") {
      yl.chx.dt <- Reduce(function(x, y) merge(x, y, all = TRUE), fixed.dt)
      yl.chx.dt[, es.chx.ts := log2((t.chx.tpm/1e6)/(s.chx.tpm/1e6))]
      yl.chx.dt[, es.chx.ms := log2((m.chx.tpm/1e6)/(s.chx.tpm/1e6))]
      yl.chx.dt[, es.neg.ts := log2((t.neg.tpm/1e6)/(s.neg.tpm/1e6))]
      yl.chx.dt[, es.neg.ms := log2((m.neg.tpm/1e6)/(s.neg.tpm/1e6))]
    }
  }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### Protein sequence, SignalP, TopCons, TMHMM
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### Protein Sequence
if (exists("annotation.protein")) {
  protein.fasta <- Biostrings::readAAStringSet(annotation.protein)
  seq.names <- names(protein.fasta)
  protein.seq <- paste(protein.fasta)
  rm(protein.fasta)
  protein.dt <- data.table(seq.names, protein.seq)
  protein.dt[, id := as.character(lapply(
    seq.names, function(x) strsplit(x, " ")[[1]][[1]]
  ))]
  setkey(protein.dt, 'id')
  protein.dt[, protein.l := nchar(protein.seq)]
  protein.dt <- protein.dt[, .(id, protein.l, protein.seq)]
}

### SignalP_4.1
repeat {
  signalp.dt <- as.data.table(
    read.table(file = "Sequences/GS115/SignalP_4.1/cbs7435_PRT.fasta.out",
               col.names = c("id", "Cmax", "pos1", "Ymax", "sp.l", 
                             "Smax", "pos2", "Smean", "sp.score", "decision", 
                             "Dmaxcut", "Networks.used"))
  )
  signalp.dt[decision == "Y", sp.sp := TRUE]
  signalp.dt[decision == "N", sp.sp := FALSE]
  signalp.dt[, 
             c('Cmax', 'pos1', 'Ymax', 'Smax', 'pos2', 'Smean', 'Dmaxcut',
               'Networks.used', 'decision') := NULL]
  setkey(signalp.dt, "id")
  names(signalp.dt) <- c("id", "sp4.l", "sp4.score", "sp4.sp")
break }

### SignalP_5.0
repeat {
  signalp.dt <- as.data.table(
    read.delim2(
      file = "Sequences/GS115/SignalP_5.0/cbs7435_PRT_summary.signalp5",
      col.names = c("id", "prediction", "sp.perc", "other.perc", "position"), 
      as.is = c(1:5), header = FALSE
    ))
  signalp.dt[prediction == "SP(Sec/SPI)", sp := TRUE]
  signalp.dt[prediction == "OTHER", sp := FALSE]
  val_na(signalp.dt, "", c(5L))
  signalp.dt[!is.na(position),
             sp.l := sapply(position, function(x) unlist(strsplit(unlist(
               strsplit(x, "-"))[[2]], ".", fixed = TRUE))[[1]])]
  signalp.dt[, c("prediction", "other.perc", "position") := NULL]
  dsave(signalp.dt,
        directory = paste0(save_directory, "prot/"),
        file.name = "signalp_5.0.dt.RData")
  
break } 

### GPI Pred
repeat {
  gpipred.dt <- as.data.table(
    read.delim2(file = "Sequences/YEASX/GPIpred/yeast_protein_GPIpred.txt",
                col.names = c("id", "gpi.fpr", "gpi.omega"),
                as.is = c(1:3),
                header = FALSE))

  gpipred.dt[, id := gsub(">", "", id)]
  gpipred.dt[, gpi.fpr := gsub("FPrate:", "", gpi.fpr)]
  gpipred.dt[, gpi.omega := gsub("OMEGA:", "", gpi.omega)]
  gpipred.dt[, gpi.aa := sapply(
    omega, function(x) unlist(strsplit(x, "-", fixed = TRUE))[[1]])]
  gpipred.dt[, gpi.index := sapply(
    omega, function(x) unlist(strsplit(x, "-", fixed = TRUE))[[2]])]
  gpipred.dt[, gpi.omega := NULL]
  gpipred.dt[, gpi.specificity.index := (1 - gpi.fpr) * 100]
  prot.dt[, gpi.prediction := 0]
  prot.dt[gpi.specificity.index >= 99.0, gpi.prediction := 1]
  prot.dt[gpi.specificity.index >= 99.5, gpi.prediction := 2]
  prot.dt[gpi.specificity.index >= 99.9, gpi.prediction := 3]
break }

### DeepLoc
repeat {
  deeploc.dt <- as.data.table(
    read.delim2(file = "Sequences/YEASX/DeepLoc/yeast_DeepLoc.txt",
                col.names = c("id", "prediction"), 
                as.is = c(1:2),
                header = FALSE))

  deeploc.dt[, pred.loc := sapply(
    prediction, function(x) unlist(strsplit(x, ","))[[1]])]
  deeploc.dt[, pred.type := sapply(
    prediction, function(x) unlist(strsplit(x, ","))[[2]])]
  deeploc.dt[, pred.loc := gsub("Prediction: ", "", pred.loc)]
  deeploc.dt[, prediction := NULL]
break }

### TopCons
if (exists("annotation.TopCons")) {
  sys.cmd <- paste0(
    "awk '{print $7,$4,$3}' ", annotation.TopCons
  )
  topcons.dt <- as.data.table(system(sys.cmd, intern = TRUE))
  topcons.dt[, id := as.character(lapply(
    V1, function(x) strsplit(x, " ")[[1]][1]))]
  topcons.dt[, sp.tc := as.logical(lapply(
    V1, function(x) strsplit(x, " ")[[1]][2]))]
  topcons.dt[, tmd.tc.n := as.integer(lapply(
    V1, function(x) strsplit(x, " ")[[1]][3]))]
  topcons.dt[, V1 := NULL]
  setkey(topcons.dt, "id")
}

### TMHMM
if (exists("annotation.TMHMM")) {
  prob.cmd <- paste0(
    "grep \"prob\" ", annotation.TMHMM
  )
  probs <- as.data.table(system(prob.cmd, intern = TRUE))
  probs[, id := as.character(lapply(
    V1, function(x) strsplit(x, " ")[[1]][2]))]
  probs[, tmd.tmhmm.p := as.numeric(lapply(
    V1, function(x) strsplit(x, " ")[[1]][14]))]
  setkey(probs, "id")
  
  pred.cmd <- paste0(
    "grep \"predicted\" ", annotation.TMHMM
  )
  preds <- as.data.table(system(pred.cmd, intern = TRUE))
  preds[, id := as.character(lapply(
    V1, function(x) strsplit(x, " ")[[1]][2]))]
  preds[, tmd.tmhmm.n := as.numeric(lapply(
    V1, function(x) strsplit(x, " ")[[1]][8]))]
  setkey(preds, "id")
  
  tmhmm.dt <- merge(preds[, .(id, tmd.tmhmm.n)], probs[, .(id, tmd.tmhmm.p)],
                    all = TRUE)
  rm(preds, probs)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### Ontology
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### EggNOG
egg_file <-
  "Sequences/GS115/Annotations/cbs7435_PRT_2.fasta.emapper.annotations"
if (TRUE) {
  ## Read eggnog output file into R as a data table, setkey to 'id'
  eggnog.dt <- as.data.table(read.delim2(
    file =  egg_file,
    header = FALSE, as.is = c(1:13),
    na.strings = c("NA|NA|NA", ""),
    col.names = c("id",	"eggnog.match",	"eggnog.value",	"eggnog.score",
                  "eggnog.name",	"eggnog.GO",	"eggnog.KEGG",	"eggnog.BiGG",
                  "eggnog.taxon","eggnog.OG",	"eggnog.HMM",	"eggnog.cog",	
                  "eggnog.description")
  ))
  eggnog.dt <- eggnog.dt[,.(id, eggnog.name, eggnog.cog,
                            eggnog.description, eggnog.HMM)]
  setkey(eggnog.dt, 'id')
  
  ## Convert cog score output into list of cog scores
  eggnog.dt[, eggnog.cog := gsub(',', '', eggnog.cog)]
  eggnog.dt[, eggnog.cog := strsplit(eggnog.cog, ' ')]
  
  ## Create NEW data table by duplicating id n times, where n = length(cog.list)
  ## , and pairing with unique cog score from cog.list
  for (i in 1:nrow(eggnog.dt)) {
    
    dt <- as.data.table(eggnog.dt$id[[i]])
    
    for (j in 1:length(eggnog.dt[dt[[1]]]$eggnog.cog[[1]])) {
      eggnog.cog <- eggnog.dt[dt[[1]]]$eggnog.cog[[1]][j]
      tmp <- copy(dt)
      tmp[, V2 := eggnog.cog]
      if (j == 1) {
        dt.hold <- tmp
      } else {
        dt.hold <- merge(dt.hold, tmp, all = TRUE)
      }
    }
    if (i == 1) {
      dt.final <- dt.hold
    } else {
      dt.final <- merge(dt.final, dt.hold, all = TRUE)
    }
  }
  cog.dt <- dt.final; rm(dt.final)
  names(cog.dt) <- c('id', 'eggnog.cog')
  
  ## Create cog.table with relative cog frequency for hierarchichal clustering
  cog.table <- as.data.table(table(cog.dt$eggnog.cog))
  names(cog.table) <- c('eggnog.cog', 'cog.n')
  cog.table[, cog.freq := cog.n/(sum(cog.table$cog.n))]
  setkey(cog.table, "cog.freq")
  egg.hierarchy <- rev(cog.table$eggnog.cog)
  
  ## Cluster duplicated cog scores into highest frequency cog from cog.table
  while (nrow(cog.dt[duplicated(cog.dt$id)]) > 0) {
    dt <- cog.dt[duplicated(cog.dt$id)]
    
    id <- dt$id[[1]]
    
    obey(paste0("match.list <- dt[, id == '", id, "']"))
    potentials <- dt[match.list]$eggnog.cog
    
    tmp.best <- length(egg.hierarchy)
    for (i in 1:length(potentials)) {
      place <- lapply(egg.hierarchy, 
                      function(x) grep(potentials[[i]], x))
      place <- which(place == 1)
      if (place < tmp.best) {
        tmp.best  <- place
      }
    }
    best.fit <- egg.hierarchy[[tmp.best]]
    obey(paste0("cog.dt[id == '", id, "', eggnog.cog := '", best.fit, "']"))
    cog.dt <- unique(cog.dt)
    cat(paste0(
      "\nCHANGE: ", id, "'s COG value will be changed to ", best.fit, "\n"))
    cat(paste0(
      "Duplicates left: ",
      as.character(nrow(cog.dt[duplicated(cog.dt$id)])), "\n"))
    
  }
  
  eggnog.dt <- merge(
    eggnog.dt[,.(id, eggnog.name, eggnog.description, eggnog.HMM)],
    cog.dt,
    all = TRUE)
  
  ## Get COG tranlsations
  cog.def <- as.data.table(read.delim2(
    file = "/Users/mandarax/Public/RStudio/ontology/defCOG.txt",
    header = TRUE,
    as.is = (c(1:4)),
    col.names = c("eggnog.cog", "category", "subcategory", "color")
  ))
  
  eggnog.dt <- merge(eggnog.dt, cog.def,
                     by.x = "eggnog.cog", by.y = "eggnog.cog",
                     all = TRUE)
  setcolorder(eggnog.dt, c(2:5,1,6:8))
  setkey(eggnog.dt, id)
  eggnog.dt <- eggnog.dt[!is.na(id)]
  eggnog.dt[, color := sapply(color, function(x) rand_color(x))]
  
  rm(cog.def, cog.dt, cog.table, dt, dt.hold, tmp)
  
  eggnog.dt[is.na(eggnog.name), eggnog.name := 
              sapply(eggnog.HMM,
                     function(x) unlist(strsplit(x, "|", fixed = TRUE))[[1]])]
}

eggnog_1.dt <- eggnog.dt; rm(eggnog.dt)
eggnog_2.dt <- eggnog.dt; rm(eggnog.dt)

repeat { 
  eggnog.dt <- merge(eggnog_1.dt, eggnog_2.dt, all = TRUE)
  eggnog.dt[is.na(eggnog.name.x), eggnog.name.x := eggnog.name.y]
  eggnog.dt[is.na(eggnog.description.x), 
            eggnog.description.x := eggnog.description.y]
  eggnog.dt[is.na(eggnog.HMM.x), eggnog.HMM.x := eggnog.HMM.y]
  eggnog.dt[is.na(eggnog.cog.x), eggnog.cog.x := eggnog.cog.y]
  eggnog.dt[is.na(category.x), category.x := category.y]
  eggnog.dt[is.na(subcategory.x), subcategory.x := subcategory.y]
  eggnog.dt[is.na(color.x), color.x := color.y]
  eggnog.dt[, c("eggnog.name.y", "eggnog.description.y", "eggnog.HMM.y", 
                "eggnog.cog.y", "category.y", "subcategory.y",
                "color.y") := NULL]
  names(eggnog.dt) <- gsub(".x", "", names(eggnog.dt))
  load(
    paste0("~/Public/RStudio/samples/", cell.line, 
           "_files/prot/protein.dt.RData"))
  eggnog.dt <- merge(eggnog.dt, protein.dt[, .(id)], 
                     by.x = 'id', by.y = 'id', all = TRUE)
  eggnog.dt[is.na(eggnog.name), eggnog.name := id]
  eggnog.dt[is.na(eggnog.cog), eggnog.cog := "S"]
  eggnog.dt[is.na(color), color := "17,50,63"]
  eggnog.dt[is.na(category), category := "Poorly Characterized"]
  eggnog.dt[is.na(subcategory), subcategory := "Function Unknown"]
  eggnog.dt[is.na(eggnog.HMM), 
            color := sapply(color, function(x) rand_color(x))]
  
  voronoi.map <- na.omit(
    eggnog.dt[, .(category, subcategory, eggnog.name, id, color)])
  voronoi.map[, `protein:ID` := paste0(eggnog.name, ":", id)]
  voronoi.map[, c('eggnog.name', 'id') := NULL]
  setcolorder(voronoi.map, c(1:2,4))
  voronoi.map[, position := ""]
  
  write.table(voronoi.map, 
              file = paste0(save_directory, "tessellations/color_position.txt"),
              append = TRUE, 
              quote = FALSE, 
              sep = "\t", 
              row.names = FALSE,
              col.names = FALSE)
break }

### GO-Slim
repeat {
  ## GO-Slim prediction table from Yeast Genome's GO-Slim Wrapper
  ont.dt <- as.data.table(read.table(
    'Sequences/YEASX/Annotations/GO-Slim/GO_Slim_Process_Yeast_1.xls', 
    header = TRUE, sep = "\t", as.is = c(2:5),
    col.names = c("GO.id", "subsubcategory", "freq.x", "gfreq.x", "genes.x")))
  
  ont2.dt <- as.data.table(read.table(
    'Sequences/YEASX/Annotations/GO-Slim/GO_Slim_Process_Yeast_2.xls', 
    header = TRUE, sep = "\t", as.is = c(2:5),
    col.names = c("GO.id", "subsubcategory", "freq.y", "gfreq.y", "genes.y")))
  
  ont.dt <- merge(ont.dt[, .(subsubcategory, freq.x, gfreq.x, genes.x)], 
                  ont2.dt[, .(subsubcategory, freq.y, gfreq.y, genes.y)], 
                  by = "subsubcategory", all = TRUE); rm(ont2.dt)
  
  ## Develop subsubcategory hierarchy list
  if (!exists("go.hierarchy")) {
    ont.dt[, freq.x := as.numeric(gsub("%", "", as.character(lapply(
      freq.x, function(x) strsplit(x, " ")[[1]][6]))))]
    ont.dt[, freq.y := as.numeric(gsub("%", "", as.character(lapply(
      freq.y, function(x) strsplit(x, " ")[[1]][6]))))]
    ont.dt[, gfreq.x := as.numeric(gsub("%", "", as.character(lapply(
      gfreq.x, function(x) strsplit(x, " ")[[1]][5]))))]
    ont.dt[, gfreq.y := as.numeric(gsub("%", "", as.character(lapply(
      gfreq.y, function(x) strsplit(x, " ")[[1]][5]))))]
    na_zero(ont.dt, c(2L,3L,5L,6L))
    
    ont.dt[, hierarchy := (((freq.x + gfreq.x)/2) + (((freq.y + gfreq.y)/2)))]
    go.hierarchy <- ont.dt[, .(hierarchy, subsubcategory)]
    setkey(go.hierarchy, "hierarchy")
    go.hierarchy <- go.hierarchy$subsubcategory
    go.hierarchy <- rev(go.hierarchy)
  }
  
  ## Merge gene lists together by subsubcategory
  ont.dt[, genes := paste(gsub(',', ' ', genes.x), gsub(',', ' ', genes.y))]
  ont.dt[, c('freq.x', 'freq.y', 'gfreq.x', 'gfreq.y',
             'genes.x', 'genes.y') := NULL]
  
  ## Convert genes to list of genes from character
  ont.dt[, genes := strsplit(genes, split = ' ')]
  
  ## Convert genes list to individual gene/subsubcategory combination
  for (i in 1:nrow(ont.dt)) {
    dt <- as.data.table(ont.dt$genes[[i]])
    
    for (j in 1:nrow(dt)) {
      subsubcategory <- ont.dt$subsubcategory[[i]]
      dt[, V2 := subsubcategory]
    }
    if (i == 1) {
      dt.2 <- dt
    }
    dt.2 <- merge(dt.2, dt, all = T)
  }
  ont.dt <- dt.2; rm(dt.2)
  names(ont.dt) <- c('gene.name', 'subsubcategory')
  
  ## Remove gene.names that do not appear in sc.formals
  load("~/Public/RStudio/ontology/gene.names.RData")
  ont.dt <- ont.dt[gene.name %in% gene.names]
  
  ## Remove duplicated gene.name:subsubcategory combinations using GO.hierarchy
  while (nrow(ont.dt[duplicated(ont.dt$gene.name)]) > 0) {
    
    dt <- ont.dt[duplicated(ont.dt$gene.name)]
    
    name <- dt$gene.name[[1]]
    
    eval(parse(text = paste0(
      "match.list <- dt[, gene.name == '", name, "']"
    )))
    potentials <- dt[match.list]$subsubcategory
    
    tmp.best <- length(go.hierarchy)
    for (i in 1:length(potentials)) {
      place <- lapply(go.hierarchy, 
                      function(x) grep(paste0("^", potentials[[i]]), x))
      place <- which(place == 1)
      if (place < tmp.best) {
        tmp.best <- place
      }
    }
    best.fit <- go.hierarchy[[tmp.best]]
    eval(parse(text = paste0(
      "ont.dt[gene.name == '", name, "', subsubcategory := '", best.fit, "']"
    )))
    ont.dt <- unique(ont.dt)
    cat(paste0(
      '\nWARNING: ', name, ' category will be changed to ', best.fit, ':\n'))
    cat(paste0(
      'Duplicates left:', 
      as.character(nrow(ont.dt[duplicated(ont.dt$gene.name)])), '\n'))
  }
  
  ont.dt[, subcategory := lapply(subsubcategory, generalize)]
  ont.dt[, category := lapply(subcategory, summarize)]
  
  ## Change colors to subcategories
  ont.dt[, color := "rgb(r,g,b)"]
  load("~/Public/RStudio/ontology/subcategories.RData")
  load("~/Public/RStudio/ontology/subcategories_colors.RData")
  
  for (i in 1:length(subcategories)) {
    obey(paste0(
      "ont.dt[subcategory == '", subcategories[[i]], "', ",
      "color := as.character(lapply(color, function(x) rand_color(",
      subcategories_colors[[i]], ")))]"
    ))
  }
  
  ## Change characters to title case
  ont.dt[, subsubcategory := toTitleCase(subsubcategory)]
  break
}

