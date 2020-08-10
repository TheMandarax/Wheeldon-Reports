### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Data read/write
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Clean system environment using MATLAB like calls
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

clear.functions <- structure(function(){}, class = "clear.functions")
print.clear.functions <- function(clear.functionsObject) {
  lst <- ls.str(mode = "function", envir = .GlobalEnv)
  
  rm(list = lst, envir = .GlobalEnv)
}

##  Save data object to directory or save_directory using file or object name
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

obey <- function(string) {
  eval(parse(text = string), envir = .GlobalEnv)
}

