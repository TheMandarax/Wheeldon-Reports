#### Finish prot.dt with cell wall changes ####

prot.dt[, wall := sapply(description, function(x) grepl("cell wall", x))]
prot.dt[wall == FALSE, wall := sapply(description, function(x) grepl("Cell wall", x))]
prot.dt[wall == FALSE, wall := sapply(description, function(x) grepl("cell-wall", x))]
prot.dt[subcategory == "Function Unknown" & wall == TRUE, color := sapply(color, function(x) rand_color("169,111,164"))]
prot.dt[subcategory == "Function Unknown" & wall == TRUE, category := "Cellular Processes and Signaling"]
prot.dt[subcategory == "Function Unknown" & wall == TRUE, subcategory := "Cell Wall/Membrane/Envelope Biogenesis"]
prot.dt[, wall := NULL]

voronoi.dt<- na.omit(prot.dt[, .(category, subcategory, gene, id, color)])
voronoi.dt[, `protein:ID` := paste0(gene, ":", id)]
voronoi.dt[, c('gene', 'id') := NULL]
setcolorder(voronoi.dt, c(1:2,4))
voronoi.dt[, position := ""]

write.table(voronoi.dt,
            file = "/Users/mandarax/Public/RStudio/CRG/Data/tr_files/RUTC30/tessellations/vironoi_map.txt",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

#### Protein Sequence ####
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

#### Blastp ####
blast.file <- 
  "/Users/mandarax/Public/PyCharm/Sequencing/Sequences/YARLI/blast/blast_GCF_000002525.2_ASM252v1_protein_yeast.txt"

blast.dt <- as.data.table(read.delim(
  file = blast.file,
  header = FALSE,
  col.names = c("id", "match.id", "aln.percentage", "aln.length", "mismatches",
                "gapOpenings", "query.start", "query.end", "match.start", 
                "match.end", "evalue", "bitscore")
))
setkey(blast.dt, id)

unique.ids <- unique(blast.dt$id)

for (i in 1:length(unique.ids)) {
  cat(paste0("Generating best match for ", unique.ids[[i]], "...\n"))
  
  gene <- blast.dt[id == unique.ids[[i]]]
  
  tmp <- gene[evalue == min(gene$evalue)]
  if (nrow(tmp) > 1) {tmp <- tmp[bitscore == max(tmp$bitscore)]}
  if (nrow(tmp) > 1) {tmp <- tmp[1]}
  
  if (i == 1) {tmp.dt <- tmp} else {tmp.dt <- rbind(tmp.dt, tmp)}
}

blast.dt <- tmp.dt[,.(id, match.id)]

#### SignalP_4.1 ####
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

#### SignalP_5.0 ####
repeat {
  signalp.dt <- as.data.table(
    read.delim2(
      file = "/Users/mandarax/Public/PyCharm/Sequencing/Sequences/GS115/SignalP_5.0/Ppastoris_summary.signalp5",
      col.names = c("id", "prediction", "sp.perc", "other.perc", "position"), 
      as.is = c(1:5), header = FALSE
    ))
  signalp.dt[prediction == "SP(Sec/SPI)", sp := TRUE]
  signalp.dt[prediction == "OTHER", sp := FALSE]
  val_na(signalp.dt, "", c(5L))
  signalp.dt[!is.na(position),
             sp.l := sapply(position, function(x) unlist(strsplit(unlist(
               strsplit(x, "-"))[[2]], ".", fixed = TRUE))[[1]])]
  signalp.dt[, c("prediction", "other.perc") := NULL]
  dsave(signalp.dt,
        directory = paste0(save_directory, "prot/"),
        file.name = "signalp.dt")
  
  break } 

#### GPI Pred ####
repeat {
  gpipred.dt <- as.data.table(
    read.delim2(file = "/Users/mandarax/Public/PyCharm/Sequencing/Sequences/YEASX/GPIpred/yeast_protein_GPIpred.txt",
                col.names = c("id", "gpi.fpr", "gpi.omega"),
                as.is = c(1:3),
                header = FALSE))
  
  gpipred.dt[, id := gsub(">", "", id)]
  gpipred.dt[, gpi.fpr := gsub("FPrate:", "", gpi.fpr)]
  gpipred.dt[, gpi.omega := gsub("OMEGA:", "", gpi.omega)]
  gpipred.dt[, gpi.aa := sapply(
    gpi.omega, function(x) unlist(strsplit(x, "-", fixed = TRUE))[[1]])]
  gpipred.dt[, gpi.index := sapply(
    gpi.omega, function(x) unlist(strsplit(x, "-", fixed = TRUE))[[2]])]
  gpipred.dt[, gpi.omega := NULL]
  gpipred.dt[, gpi.specificity.index := (1 - as.numeric(gpi.fpr)) * 100]
  gpipred.dt[, gpi.prediction := 0]
  gpipred.dt[gpi.specificity.index >= 99.0, gpi.prediction := 1]
  gpipred.dt[gpi.specificity.index >= 99.5, gpi.prediction := 2]
  gpipred.dt[gpi.specificity.index >= 99.9, gpi.prediction := 3]
  break }

#### DeepLoc ####
repeat {
  deeploc.dt <- as.data.table(
    read.delim2(file = "/Users/mandarax/Public/PyCharm/Sequencing/Sequences/GS115/DeepLoc/deepLoc.txt",
                col.names = c("id", "dl.loc", "dl.type", "dl.nucleus", "dl.cytoplasm", "dl.extracellular", "dl.mitochondrion", "dl.membrane", "dl.er", "dl.plastid", "dl.golgi", "dl.lysosome", "dl.peroxisome"), 
                header = FALSE))
  
  deeploc.dt[, pred.loc := sapply(
    prediction, function(x) unlist(strsplit(x, ","))[[1]])]
  deeploc.dt[, pred.type := sapply(
    prediction, function(x) unlist(strsplit(x, ","))[[2]])]
  deeploc.dt[, pred.loc := gsub("Prediction: ", "", pred.loc)]
  deeploc.dt[, prediction := NULL]
  break }

#### TopCons ####
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

#### TMHMM ####
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

#### EggNOG 1.0 ####

egg_file <- "/Users/mandarax/Public/PyCharm/Sequencing/Sequences/RUTC30/eggnog/GCA_000513815.1_TrireRUTC30_1_protein.2.faa.emapper.annotations"

if (TRUE) {
  ## Read eggnog output file into R as a data table, setkey to 'id'
  eggnog.dt <- as.data.table(read.delim2(
    file =  egg_file,
    header = FALSE, as.is = c(1:13),
    na.strings = c("NA|NA|NA", ""),
    col.names = c("id",	"eggnog.match",	"eggnog.value",	"eggnog.score",
                  "gene.eggnog",	"eggnog.GO",	"eggnog.KEGG",	"eggnog.BiGG",
                  "eggnog.taxon","eggnog.OG",	"eggnog.HMM",	"eggnog.cog",	
                  "description.eggnog")
  ))
  eggnog.dt <- eggnog.dt[,.(id, gene.eggnog, eggnog.cog,
                            description.eggnog, eggnog.HMM)]
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
    eggnog.dt[,.(id, gene.eggnog, description.eggnog, eggnog.HMM)],
    cog.dt,
    all = TRUE)
  
  ## Get COG tranlsations
  cog.def <- as.data.table(read.delim2(
    file = "/Users/mandarax/Public/RStudio/CRG/Data/ontology/defCOG.txt",
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
  
  eggnog.dt[is.na(gene.eggnog), gene.eggnog := 
              sapply(eggnog.HMM,
                     function(x) unlist(strsplit(x, "|", fixed = TRUE))[[1]])]
}

eggnog_1.dt <- eggnog.dt; rm(eggnog.dt)
eggnog_2.dt <- eggnog.dt; rm(eggnog.dt)
eggnog_3.dt <- eggnog.dt; rm(eggnog.dt)

cell.line <- 'tr'

repeat { 
  eggnog.dt <- merge(eggnog.dt, eggnog_3.dt, all = TRUE)
  eggnog.dt[is.na(gene.eggnog.x), gene.eggnog.x := gene.eggnog.y]
  eggnog.dt[is.na(description.eggnog.x), 
            description.eggnog.x := description.eggnog.y]
  eggnog.dt[is.na(eggnog.HMM.x), eggnog.HMM.x := eggnog.HMM.y]
  eggnog.dt[is.na(eggnog.cog.x), eggnog.cog.x := eggnog.cog.y]
  eggnog.dt[is.na(category.x), category.x := category.y]
  eggnog.dt[is.na(subcategory.x), subcategory.x := subcategory.y]
  eggnog.dt[is.na(color.x), color.x := color.y]
  eggnog.dt[, c("gene.eggnog.y", "description.eggnog.y", "eggnog.HMM.y", 
                "eggnog.cog.y", "category.y", "subcategory.y",
                "color.y") := NULL]
  names(eggnog.dt) <- gsub(".x", "", names(eggnog.dt))
  
break }

# Load protein.dt to complete eggnog.dt
eggnog.dt <- merge(eggnog.dt, protein.dt[, .(id)],
                   by.x = 'id', by.y = 'id', all = TRUE)
eggnog.dt[is.na(gene.eggnog), gene.eggnog := id]
eggnog.dt[is.na(eggnog.cog), eggnog.cog := "S"]
eggnog.dt[is.na(color), color := "17,50,63"]
eggnog.dt[is.na(category), category := "Poorly Characterized"]
eggnog.dt[is.na(subcategory), subcategory := "Function Unknown"]
eggnog.dt[is.na(eggnog.HMM), 
          color := sapply(color, function(x) rand_color(x))]

#### EggNOG 2.0 ####
egg_file <-
  ""

if (TRUE) {
  ## Read eggnog output file into R as a data table, setkey to 'id'
  eggnog.dt <- as.data.table(read.delim2(
    file =  egg_file,
    header = FALSE, 
    na.strings = c("NA|NA|NA", ""),
    col.names = c("id",	"egg.ortholog",	"egg.value", "egg.score", "egg.tax",
                  "eggnog.name",	"egg.GO",	"egg.EC", "kegg.ko", "kegg.pathway",
                  "kegg.module",	"kegg.rxn",	"kegg.rclass", "egg.brite",
                  "kegg.tc", "egg.CAZy", "BiGG.rxn", "egg.annotlvl", "egg.og",
                  "egg.bestog", "eggnog.cog", "description.eggnog")
  ))
  eggnog.dt <- eggnog.dt[,.(id, eggnog.name, eggnog.cog,
                            description.eggnog)]
  setkey(eggnog.dt, 'id')
  
  ## Convert cog score output into list of cog scores
  eggnog.dt[, eggnog.cog := as.character(eggnog.cog)]
  eggnog.dt[is.na(eggnog.cog), eggnog.cog := "S"]
  eggnog.dt[, eggnog.cog := strsplit(eggnog.cog, '')]
  
  ## Create NEW data table by duplicating id n times, where n = length(cog.list)
  ## , and pairing with unique cog score from cog.list
  for (i in 1:nrow(eggnog.dt)) {
    
    dt <- as.data.table(eggnog.dt$id[[i]])
    
    for (j in 1:length(unlist(eggnog.dt[dt]$eggnog.cog))) {
      eggnog.cog <- unlist(eggnog.dt[dt]$eggnog.cog)[j]
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
    eggnog.dt[,.(id, eggnog.name, description.eggnog, eggnog.HMM)],
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