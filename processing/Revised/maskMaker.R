library(Biostrings)
library(data.table)

# Necessary functions
strpick <- function(string, split = " ", pos = NA) {
  nStr <- length(unlist(strsplit(string, split, fixed = TRUE)))
  str <- unlist(strsplit(string, split, fixed = TRUE))
  
  if (is.na(pos)) {
    return(str)
  }
  else {
    if (nStr == 1) {
      return(str[[nStr]])
    }
    else if (pos > nStr) {
      return(str[[nStr]])
    }
    else {
      return(str[[pos]])
    }
  }
}

partition_sequence <- function(string, footprint = 28) {
  gene.l <- length(unlist(strsplit(string, "", fixed = TRUE)))
  gene.seq <- unlist(strsplit(string, "", fixed = TRUE))
  
  if (gene.l <= footprint) {
    stop(paste0(
      "Gene length too low to partition for ribosome profiling (",
      footprint, " nt footprints)!\n")
    )
  }
  else {
    txt.fasta <- ""
    for (i in 1:(gene.l - footprint)) {
      start <- i
      start
      end <- i + (footprint - 1)
      end
      txt.sample <- paste0(gene.seq[start:end], collapse = "")
      if (i == 1) {
        txt.fasta <- txt.sample
      }
      else {
        txt.fasta <- paste0(txt.fasta, "\n", txt.sample)
      }
    }
    cmd <- paste0(
      "ORF scrambled... ", i, " unique fragments ", footprint,
      " nucleotides long produced. \n"
    )
    cat(cmd)
    return(txt.fasta)
  }
}

# Read CDS fasta file using Biostrings
annotation_cds <- 
  "/Users/mandarax/Public/PyCharm/Sequencing/Sequences/GS115/Annotations/GS115_CRG_cds.fa"

dna_fasta <- Biostrings::readDNAStringSet(annotation_cds)
seq.names <- names(dna_fasta)
gene.seq <- paste(dna_fasta)

# Make table of genes and their sequences
dna_dt <- data.table(seq.names, gene.seq)

# Get gene id (this changes depending on the way that Biostrings reads names)
dna_dt[, id := sapply(seq.names, 
                      function(x) strpick(x, split = " ", pos = 1))]

# Get gene name (not all files will produce a gene name)
dna_dt[, gene.name := sapply(seq.names, 
                             function(x) strpick(x, split = " ", pos = 2))]

# Tidy up data table
dna_dt[, seq.names := NULL]
setcolorder(dna_dt, c('id', 'gene.name', 'gene.seq'))

# Partition DNA sequences into 28mer fragments and generate DNA fasta file 
# representing pseudo Ribo-seq data. Use this dat to run Hisat2 on to get 
# sort.bam files needed to create mask
dna_dt[, gene.seq.mix := lapply(gene.seq, function(x) partition_sequence(x))]

output_txt <- "kphaffii_mix.txt"

write.table(dna_dt$gene.seq.mix, 
            file = output_txt, 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE, 
            sep = '\n')

output_sam <- "kphaffii_mix.sam"

annotation.hisat <- 
  "/Users/mandarax/Public/PyCharm/Sequencing/Sequences/GS115/Annotations/hisat2/love.fna"

command <- paste0(
  "/Users/Shared/Repository/miniconda3/envs/seq/bin/hisat2 ",
  "-x ", annotation.hisat, " ",
  "-U ", file.txt, " ",
  "-p 7 ",
  "-r ",
  "-S ", output.sam
)
system(command)