
#' Get the base location of data for the SCS package
#'
#' @export
baseLoc <- function() {
  system.file(package='SCS')
}

#' Get the path of external files (tables, etc.) for SCS
#'
#' @export
extPath <- function() {
  file.path(baseLoc(), 'extdata')
}

#' Get the location of a data set file from a base name
#'
#' @export
data_set_file <- function (name, base.loc=baseLoc())  file.path(base.loc, "data_sets", paste0(name, ".rda"))

#' Load the data set into the global environment
#'
#' We use this instead of load() and use_data from usethis because large data sets take a
#' long time to load and save and must be reloaded every time there is a code change.
#'
#' @param name a string that represents the base name of the data set (both
#'              the file name and the contained object)
#' @param base.loc the base location of the files
#' @param bind if TRUE, bind the value to the global environment, with the given name, as is default with load()
#'
#' @return the object requested
#'
#' @export
load_data_set <- function(name, base.loc=baseLoc(), bind=TRUE) {
  loc = data_set_file(name, base.loc=base.loc)
  e <- if(bind) .GlobalEnv else new.env()
  load(loc, envir=e)
  get(as.character(name), envir=e)
}

#' Save a data set into the package data set space
#'
#' We use this instead of save() and use_data from usethis because large data sets take a
#' long time to load and save and must be reloaded every time there is a code change.
#'
#' @param ...  a symbol that represents the base name of the data set
#' @param name a character vector naming an object in the global environment to be saved.
#' @param compress by default, do not compress the file (greatly decreasing time and
#'   a minor increase disk space in most cases.)
#' @param base.loc the base location of the files
#' @param overwrite if set to T 2replace the old file if it exists
#'
#' @export
save_data_set <- function(..., name=NULL, compress=F, base.loc=baseLoc(), overwrite=F) {
  if(!missing(...) && !is.null(name)) {
    stop("object and name specification are mutually exclusive")
  }
  if(is.null(name)) {
    force(substitute(...))
    name = as.character(substitute(...))
  }

  data_file = data_set_file(name, base.loc=base.loc)
  if(file.exists(data_file) && !overwrite) {
    stop(paste("File", data_file, "exists. Use overwrite=T to override."))
  }
  
  data_dir = dirname(data_file)
  if(!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    if(!dir.exists(data_dir)) {
      stop(paste0("Could not create data directory ", data_dir, ". Check permissions."))
    }
  }
  save(list=name, file=data_file, compress=compress, envir=parent.frame())
}
