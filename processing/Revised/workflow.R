#### Package initiation ####
library(RiboProfiling)
library(GenomicAlignments)
library(GenomicRanges)
library(Rsamtools)
library(plotly)
library(zoo)
library(scales)
library(stats)
library(tidyr)
library(data.table)
library(plyr)
library(dplyr)

source("/Users/mandarax/Public/RStudio/CRG/Scripts/riboProfiling.R")
source("/Users/mandarax/Public/RStudio/CRG/Scripts/helperFunctions.R")
source("/Users/mandarax/Public/RStudio/CRG/Scripts/loadData.R")

#### Data initiation ####
# Load files and annotation file that will be used for analysis
load_files()
load_gff()

#### Riboseq workflow ####
# Convert reference to BAM files.
files <- lapply(files, Rsamtools::BamFile)

# Convert aligned reads from BAM file into GAlignments object.
aln <- lapply(files, GenomicAlignments::readGAlignments)

# Convert GAlignments object to End (3') positions.
alnGRanges <- lapply(aln, RiboProfiling::readsToStartOrEnd, what = "end")

# Returns a GRanges object containing the flank size around the transcriptional
# start site (TSS) for selected coding sequence (CDS).
flank_size <- 28
oneBinRanges <- lapply(seq(length(files)),
                       function(x) RiboProfiling::aroundPromoter(
                         txdb = tx$txdb,
                         alnGRanges = alnGRanges[[x]],
                         percBestExpressed = 0.25,
                         flankSize = flank_size))

# Create histogram of match length distribution of reads.
matchLenDistr <- lapply(aln, RiboProfiling::histMatchLength)

# Create vector of match lengths with read counts greater than 3000.
match_lengths <- lapply(seq(length(files)),
                      function(x) as.numeric(as.character(unlist(
                        matchLenDistr[[x]][[1]]$matchSize)))[
                          which(matchLenDistr[[x]][[1]]$counts > 3000)])

# Calculate summarized read coverages around TSS for specified match lengths.
listPromoterCov <- lapply(seq(length(files)),
                          function(x) RiboProfiling::readStartCov(
                            alnGRanges = alnGRanges[[x]],
                            oneBinRanges = oneBinRanges[[x]],
                            matchSize = match_lengths[[x]],
                            fixedInterval = c(-flank_size, flank_size),
                            renameChr = "aroundTSS",
                            charPerc = "sum"))

# Calculate psite offsets from listPromoterCov.
shift <- sapply(seq(length(files)),
                function(x) as.numeric(listPromoterCov[[x]]$sumUp@ranges@start[
                  which.max(listPromoterCov[[x]]$
                              sumUp@elementMetadata@listData$values)]))

# Applies psite offset on read start along trascript and returns
#   1. Information on ORF including names, position, lengths, and counts on 
#      5'UTR, CDS, and 3'UTR after offset is applied.
#   2. List of dataframes for each ORF containing read counts per codon.
counts <- lapply(seq(length(files)), 
                 function(x) countReads(
                   exonGRanges = tx$exonGRanges[names(tx$cdsPosTransc)],
                   cdsPosTransc = tx$cdsPosTransc,
                   alnGRanges = alnGRanges[[x]],
                   originalAln = aln[[x]],
                   shiftValue = shift[[x]],
                   motifSize = 3))

# Collapses list of dataframes for each ORF containing read counts per codon
# into one dataframe.
counts <- lapply(seq(length(files)), function(x) ldply(counts[[x]][[2]]))

# Apply masking to read count dataframe.
masked <- lapply(seq(length(files)), function(x) 
  full_join(mask, counts[[x]], by = c(".id", "codonID")) %>%
    filter(codonID > 5 & codonID <= theory.l - 5) %>% 
    replace_na(list(nbrReads = 0)))

# Calculate expression from masked read counts.
tmp <- lapply(seq(length(files)), function(x) 
  masked[[x]] %>% 
    filter(mask == FALSE) %>% 
    group_by(.id) %>%
    mutate(reads = sum(nbrReads)) %>% 
    select(-codonID, -nbrReads, -mask) %>% 
    distinct() %>%
    calc_expression(output = tpm, input = reads, norm = pseudo.l) %>% 
    rename(id = .id,
           !! paste0(sample_names[[x]], ".reads") := reads,
           !! paste0(sample_names[[x]], ".tpm") := tpm))

# Join samples together and calculate membrane enrichment to filter genes for 
# metagene analysis.
temp <- Reduce(function(x, y) full_join(x, y), tmp) %>% 
  calc_enrich(output = es.chx.ms, mem = chx.m.tpm, sol = chx.s.tpm) %>% 
  filter(es.chx.ms <= 1)

# Calculate average reads per codon for first 100 codons per gene to reduce 
# biases associated with increased read counts at beginning of transcripts.
rpc <- lapply(seq(length(files)), function(x)
  masked[[x]] %>% 
    filter(codonID <= 100 & mask == FALSE) %>% 
    group_by(.id) %>% 
    mutate(rpc.100 = mean(nbrReads),
           rpc.100 = na_if(rpc.100, 0)))

# Normalize reads per codon by average reads per codon for first 100 codons and
# determine reads per gene to establish cut off criterion for metagene analysis
norm <- lapply(seq(length(files)), function(x)
  full_join(masked[[x]], rpc[[x]]) %>% 
    filter(mask == FALSE) %>% 
    group_by(.id) %>% 
    fill(rpc.100) %>% 
    replace_na(list(rpc.100 = 0)) %>% 
    mutate(rpc.100 = na_if(rpc.100, 0),
           norm.100 = nbrReads/rpc.100,
           rpg = mean(nbrReads)) %>% 
    replace_na(list(norm.100 = 0)))

# Perform metagene analysis using non-enriched genes with adequate reads per 
# gene and average reads per codon for first 100 codons. Metagene analysis 
# calculates average normalized reads per codon for all genes. This average 
# normalized read per codon is smoothed using a rolling mean and rolling median. 
meta <- lapply(seq(length(files)), function(x)
  norm[[x]] %>% 
    filter(rpg > 0.5 & rpc.100 > 0 & .id %in% temp$id) %>% 
    group_by(codonID) %>% 
    mutate(meta = mean(norm.100, trim = 0.05, na.rm = TRUE)) %>% 
    select(codonID, meta) %>% 
    distinct() %>% 
    ungroup() %>% 
    full_join(meta[[x]] %>% 
                filter(codonID >= 0 & codonID <= 100) %>% 
                mutate(smoothed = rollapply(meta,
                                            10, mean, partial = TRUE))) %>% 
    full_join(meta[[x]] %>% 
                filter(codonID > 100 & codonID <= 1000) %>% 
                mutate(smoothed = rollapply(meta, 
                                            100, mean, partial = TRUE))) %>% 
    full_join(meta[[x]] %>% 
                filter(codonID > 1000) %>% 
                mutate(smoothed = rollapply(meta, 1000, 
                                            median, partial = TRUE))) %>% 
    filter(!is.na(smoothed)) %>% 
    mutate(smoothed = na_if(smoothed, 0)) %>% 
    fill(smoothed))

# Scale number of reads per codon per gene after normalizing reads per codon per 
# gene by smoothed metagene reads per codon. Use scaled reads per codon per gene
# to calculate expression in RPM and cTPM.
scaled <- lapply(seq(length(files)), function(x)
  full_join(masked[[x]], meta[[x]]) %>% 
    fill(smoothed) %>% 
    mutate(scaledReads = nbrReads/smoothed) %>% 
    group_by(.id) %>% 
    mutate(reads = sum(scaledReads)) %>% 
    calc_expression(output = rpm, input = nbrReads) %>% 
    select(-codonID, -nbrReads, -mask, -smoothed, -meta, -scaledReads) %>% 
    distinct() %>% 
    calc_expression(output = ctpm, input = reads, norm = pseudo.l) %>% 
    rename(id = .id,
           !! paste0(sample_names[[x]], ".reads") := reads,
           !! paste0(sample_names[[x]], ".rpm") := rpm,
           !! paste0(sample_names[[x]], ".ctpm") := ctpm))

# Collapse different samples' dataframes into one dataframe and calculate 
# membrane enrichment.
dt <- Reduce(function(x, y) full_join(x, y), scaled) %>% 
  calc_enrich(output = es.chx.ms, mem = chx.m.ctpm, sol = chx.s.ctpm)


