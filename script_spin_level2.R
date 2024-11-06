##### SCRIPT: SPectrometer for Ice Nucleation (SPIN) Level 2 #####
#' @author Christopher Rapp
#' @description
#'
#'
#'
# ---------------------------------------------------------------------------- #
##### SCRIPT: Libraries, Functions, Paths #####
#' These are directories that are in use for this code.
#' If files or scripts are renamed, moved, or changed, the code will NOT work!
#' Please keep updated for future readers of this code!
#'

{
  # Clears global environment and all saved variables
  rm(list = ls())

  # Data wrangling libraries
  library(data.table)
  library(lubridate)
  library(dplyr)
  library(stringr)
  library(tidyr)

  # Plotting libraries
  library(ggplot2)
  library(ggpubr)
  library(ggsci)

  setwd("~/Library/CloudStorage/Box-Box/SPIN1-FIN03")
  work.dir <- getwd()

  #' @import
  #' Specify import directory
  #' FIN03
  import.spin = paste0(work.dir, "/export/level1/")

  #' @export
  #' Specify where to send plots and text files
  export.data = paste0(work.dir, "/export/level2/")
  export.plot = paste0(work.dir, "/export/plots/")
}

# ---------------------------------------------------------------------------- #
##### SCRIPT: Combine Dataset #####
#'

{
  # List all SPIN files exported
  files.spin <- list.files(
    path = import.spin,
    recursive = TRUE,
    full.names = TRUE,
    pattern = '*.csv'
  )

  data.ls <- lapply(files.spin, function(x){

    tmp.df <- fread(x)

    if (nrow(tmp.df) > 0){
      return(tmp.df)
    }
  })

  data.df <- rbindlist(data.ls, fill = T)

  loop.time = format(Sys.time(), "%y%m%d%H%M%S")

  filename = paste0(export.data, "SPIN1_FIN03_Version_", loop.time, "_level2", ".csv")

  # Save data using data.table::fwrite
  data.table::fwrite(data.df, file = filename, showProgress = T, append = T)
}

