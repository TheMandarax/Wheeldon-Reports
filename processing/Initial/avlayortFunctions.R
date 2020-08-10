### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Clean system environment using MATLAB like calls
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
close.all <- structure(function(){}, class = "close.all")
print.close.all <- function(close.allObject) {
  if(!is.null(dev.list())) {invisible(dev.off())}
}

clear.all <- structure(function(){}, class = "clear.all")
print.clear.all <- function(clear.allObject) {
  rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)
}

clc <- structure(function(){}, class = "clc")
print.clc <- function(clcObject) {cat("\014")}

clear.vars <- structure(function(){}, class = "clear.vars")
print.clear.vars <- function(clear.varsObject) {
  lst <- c(ls.str(mode = "list", envir = .GlobalEnv),
           ls.str(mode = "character", envir = .GlobalEnv),
           ls.str(mode = "name", envir = .GlobalEnv),
           ls.str(mode = "call", envir = .GlobalEnv),
           ls.str(mode = "logical", envir = .GlobalEnv),
           ls.str(mode = "NULL", envir = .GlobalEnv),
           ls.str(mode = "numeric", envir = .GlobalEnv))
  
  rm(list = lst, envir = .GlobalEnv)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Save data object to directory
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsave <- function(object, file.name = NA, directory = NA) {
  if (!exists(as.character(substitute(object)))) {
    warning("Object does not exist!\n")
  }
  if (is.na(directory) & !exists("save_directory", envir = .GlobalEnv)) {
    warning("Make sure 'save_directory' is present in global environment or provide a directory in function call!\n")
  }
  if (is.na(directory)) {
    if (is.na(file.name)) {
      save(list = as.character(substitute(object)), 
           file = paste0(
             save_directory, as.character(substitute(object)), ".RData"))
    } else {
      save(list = as.character(substitute(object)), file = paste0(
        save_directory, file.name, ".RData"))
    }
  } else {
    if (is.na(file.name)) {
      save(list = as.character(substitute(object)),
           file = paste0(
             directory, as.character(substitute(object)), ".RData"))
    } else {
      save(list = as.character(substitute(object)), 
           file = paste0(directory, file.name, ".RData"))
    }
  }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Perform non standard evalutation for automation
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
obey <- function(string) {
  eval(parse(text = string), envir = .GlobalEnv)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Color palette for paper
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
assign("black", "#000000", envir = .GlobalEnv)
assign("white", "#ffffff", envir = .GlobalEnv)
assign("grey", "#d3d6d9", envir = .GlobalEnv)
assign("orange", "#ff6666", envir = .GlobalEnv)
assign("pink", "#ffb7a9", envir = .GlobalEnv)
assign("yellow", "#ffe47a", envir = .GlobalEnv)
assign("teal", "#33cccc", envir = .GlobalEnv)
assign("blue", "#4a708b", envir = .GlobalEnv)
assign("lightblue1", "#809aad", envir = .GlobalEnv)
assign("lightblue2", "#a4b7c5", envir = .GlobalEnv)
assign("lightblue3", "#c8d4dc", envir = .GlobalEnv)
assign("lightblue4", "#ecf0f3", envir = .GlobalEnv)
assign("darkblue1", "#334e61", envir = .GlobalEnv)
assign("darkblue2", "#253845", envir = .GlobalEnv)
assign("darkblue3", "#162129", envir = .GlobalEnv)
assign("darkblue4", "#070b0d", envir = .GlobalEnv)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Data table manipulations
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
zero_na <- function(DT, list = names(DT)) {
  invisible(lapply(
    list, function(.name) 
      set(DT, which(DT[[.name]] == 0), j = .name, value = NA)))
}

val_na <- function(DT, var,  list = names(DT)) {
  invisible(lapply(
    list, function(.name)
      set(DT, which(DT[[.name]] == var), j = .name, var = NA)
  ))
}

na_zero <- function(DT, list = names(DT)) {
  invisible(lapply(
    list,function(.name) 
      set(DT, which(is.na(DT[[.name]])), j = .name,value = 0)))
}

nan_na <- function(DT, list = names(DT)) {
  invisible(lapply(
    list,function(.name) 
      set(DT, which(is.nan(DT[[.name]])), j = .name,value = NA)))
}

nan_zero <- function(DT, list = names(DT)) {
  invisible(lapply(
    list,function(.name) 
      set(DT, which(is.nan(DT[[.name]])), j = .name,value = 0)))
}

inf_na <- function(DT, list = names(DT)) {
  invisible(lapply(
    list,function(.name) 
      set(DT, which(is.infinite(DT[[.name]])), j = .name,value = NA)))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Calculating standard deviation from mean for data curves
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
quant_sd <- function(data, predictor, response, quantiles = 5) {
  dt <- deparse(substitute(data))
  
  pred <- deparse(substitute(predictor))
  pred.dt <- paste0(dt, "$", pred)
  cat.sd <- paste0(pred, ".sd")
  obey(paste0(
    "pred.range <- as.character(floor(quantile(", pred.dt, ", probs = seq(0, 1, ", as.character(1/quantiles), "))))"
  ))
  
  response <- deparse(substitute(response))
  
  for (i in 1:(length(pred.range) - 1)) {

    obey(paste0(
      "pred.mean <- mean(", dt, "[", pred, " >= ", pred.range[[i]], " & ", pred, " <= ", pred.range[[i + 1]], "]$", response, ")"
    ))
    obey(paste0(
      "pred.sd <- sd(", dt, "[", pred, " >= ", pred.range[[i]], " & ", pred, " <= ", pred.range[[i + 1]], "]$", response, ")"
    ))
    obey(paste0(
      dt, "[", pred, " >= ", pred.range[[i]], " & ", pred, " <= ", pred.range[[i + 1]], " & ",
      response, " >= ", pred.mean, " - 3*", pred.sd, " & ",
      response, " <= ", pred.mean, " + 3*", pred.sd, ", ",
      cat.sd, " := 3]"
    ))
    obey(paste0(
      dt, "[", pred, " >= ", pred.range[[i]], " & ", pred, " <= ", pred.range[[i + 1]], " & ",
      response, " >= ", pred.mean, " - 2*", pred.sd, " & ",
      response, " <= ", pred.mean, " + 2*", pred.sd, ", ",
      cat.sd, " := 2]"
    ))
    obey(paste0(
      dt, "[", pred, " >= ", pred.range[[i]], " & ", pred, " <= ", pred.range[[i + 1]], " & ",
      response, " >= ", pred.mean, " - ", pred.sd, " & ",
      response, " <= ", pred.mean, " + ", pred.sd, ", ",
      cat.sd, " := 1]"
    ))
    obey(paste0(
      dt, "[", pred, " > ", pred.range[[i]], " & ",
      pred, " <= ", pred.range[[i + 1]], " & ",
      "is.na(", cat.sd, "), ", 
      cat.sd, " := 4]"
    ))
  }
}

plotSummarizedCov <- function (covSummarized) {
  listPlotSum <- lapply(covSummarized, function(iSumCov) {
    maxPeak <- max(iSumCov$values)
    
    iPlot <- ggplot(iSumCov, aes(start, values/1000)) + 
      geom_point(color = blue, na.rm = TRUE) + 
      geom_line(color = blue, na.rm = TRUE) + 
      labs(x = "Distance from P-site (nt)", y = bquote(bold("Reads ("~10^3~")"))) +
      xlim(0, 30) +
      paper_theme
    return(iPlot)
  })
  
  return(listPlotSum)
}

hydrophobize <- function(sp.seq) {
  seq.list <- unlist(strsplit(sp.seq, ""))
  score <- 0
  for (i in 1:length(seq.list)) {
    if (seq.list[[i]] == "A") {score <- score + 1.8}
    else if (seq.list[[i]] == "C") {score <- score + 2.5}
    else if (seq.list[[i]] == "D" | seq.list[[i]] == "E" | seq.list[[i]] == "N" | seq.list[[i]] == "Q") {score <- score - 3.5}
    else if (seq.list[[i]] == "F") {score <- score + 2.8}
    else if (seq.list[[i]] == "G") {score <- score - 0.4}
    else if (seq.list[[i]] == "H") {score <- score - 3.2}
    else if (seq.list[[i]] == "I") {score <- score + 4.5}
    else if (seq.list[[i]] == "K") {score <- score - 3.9}
    else if (seq.list[[i]] == "L") {score <- score + 3.8}
    else if (seq.list[[i]] == "M") {score <- score + 1.9}
    else if (seq.list[[i]] == "P") {score <- score - 1.6}
    else if (seq.list[[i]] == "R") {score <- score - 4.5}
    else if (seq.list[[i]] == "S") {score <- score - 0.8}
    else if (seq.list[[i]] == "T") {score <- score - 0.7}
    else if (seq.list[[i]] == "V") {score <- score + 4.2}
    else if (seq.list[[i]] == "W") {score <- score - 0.9}
    else if (seq.list[[i]] == "Y") {score <- score - 1.3}
  }
  return(score)
}
