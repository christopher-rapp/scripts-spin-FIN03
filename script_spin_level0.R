##### SCRIPT: SPectrometer for Ice Nucleation (SPIN) Level 0 #####
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
  library(viridis)
  library(RColorBrewer)

  # Plotting parameters
  resolution.dpi = 400
  font.family = "Helvetica"

  setwd("~/Library/CloudStorage/Box-Box/SPIN1-FIN03")
  work.dir <- getwd()

  #' @import
  #' Specify import directory
  #' FIN03
  import.spin = paste0(work.dir, "/import/")

  #' @export
  #' Specify where to send plots and text files
  export.data = paste0(work.dir, "/export/level0/")
  export.plot = paste0(work.dir, "/export/plots/")

  #' @importFrom
  source(paste0("~/Documents/GitHub/functions/", "functions_microphysics.R"))
  source(paste0("~/Documents/GitHub/functions/", "functions_spin_plotting.R"))
}

# ---------------------------------------------------------------------------- #
##### SCRIPT: Start Loop #####
#' TZ LOCAL TIME DENVER STORM PEAK LAB

{
  # Directories
  spin.dirs <- list.dirs(path = import.spin,
                         recursive = FALSE,
                         full.names = FALSE)

  # Add year extension to path
  spin.path <- paste0(import.spin, "/", spin.dirs)

  # Pull POSIX style dates from directory names
  spin.dates <- as.Date(spin.dirs, format = '%Y%m%d')

  for (n in 1:length(spin.dates)){

    {
      print(paste0("Level 0 SPIN Data for ", spin.dates[n], " from directory ", spin.path[n]))

      # ---------------------------------------------------------------------- #
      ##### SCRIPT: List Files #####
      #'

      {
        # Code Timing
        pct <- proc.time()

        # List all SPIN files exported
        files.spin <- list.files(
          path = spin.path[n],
          recursive = TRUE,
          full.names = TRUE,
          pattern = '*.csv'
        )

        # List all log files
        files.log <- list.files(
          path = spin.path[n],
          recursive = TRUE,
          full.names = TRUE,
          pattern = '*.log'
        )

        # Filter files containing a possible string indicating low quality or not useable
        remove.filter <- c('temp', 'test', 'raw')

        # Apply word filter to remove files with key words
        files.spin <-
          files.spin[!grepl(paste(remove.filter, collapse = '|'),
                            ignore.case = T,
                            files.spin)]

        # Apply word filter to keep particle by particle data
        files.PbP <-
          files.spin[grepl(paste("PbP", collapse = '|'), ignore.case = T, files.spin)]

        # Only keep SPIN files, not PbP data or potentially saved pre-cooler files
        # Precooler data was originally saved into the SPIN directory but has since changed
        files.spin <-
          files.spin[!grepl(paste(c("PbP", "precooler"), collapse = '|'), ignore.case = T, files.spin)]

        rm(remove.filter)
      }

      # ---------------------------------------------------------------------- #
      ##### SECTION: Data Wrangling - SPIN #####
      #'

      {
        print("SECTION: Data Wrangling - SPIN")

        # Loop through list of files and create a dataframe of each
        data.ls <- lapply(files.spin, function(x) {
          # Use data.table's fread to read in data
          # Fastest method in reading CSV's and least memory consuming
          # fill MUST equal FALSE with SPIN files
          tmp.df <- data.table::fread(
            paste0(x),
            na.strings = c("", "NA", "NaN"),
            fill = FALSE,
            strip.white = TRUE,
            stringsAsFactors = FALSE
          )

          # Extract date and time from the filename using regular expression
          # (?<=...) is a positive-look behind assertion
          # SPIN files always follow this pattern
          tmp.date = str_extract(x, "(?<=SPIN\\d{3}\\S?)\\d{8}")
          tmp.time = str_extract(x, "(?<=SPIN\\d{3}\\S?\\d{8})\\d{6}")

          filename.c <- paste0(tmp.date, tmp.time)

          # Place the filename in the data as an identifier
          # When times are generated these should be strictly increasing
          tmp.df$Filename = filename.c

          # Identify the origin point of the file
          # Depending on when the software is restarted this may or may not be the same day for each file in a given directory
          # After 20150912 142002 the timezone of the SPIN laptop was changed
          if (filename.c < 20150912142002){

            tmp.df$Origin = first(lubridate::as_datetime(paste0(tmp.date, tmp.time), tz = "America/New_York") - lubridate::seconds(tmp.df$`Time (sec)`))
          } else {
            tmp.df$Origin = first(lubridate::as_datetime(paste0(tmp.date, tmp.time), tz = "America/Denver") - lubridate::seconds(tmp.df$`Time (sec)`))
          }

          # The time variable is seconds since midnight of the most recent day of the software started
          tmp.df <- tmp.df %>%
            mutate(`Local Time` = `Origin` + `Time (sec)`, .before = everything())

          # Identify the origin point of the file
          # Depending on when the software is restarted this may or may not be the same day for each file in a given directory
          # After 20150912 142002 the timezone of the SPIN laptop was changed
          if (filename.c < 20150912142002){

            tmp.df <- tmp.df %>%
              mutate(`Local Time` = force_tz(`Local Time`, tzone = "America/New_York")) %>%
              mutate(`Local Time` = with_tz(`Local Time`, tzone = "America/Denver"))

          } else {
            tmp.df <- tmp.df %>%
              mutate(`Local Time` = force_tz(`Local Time`, tzone = "America/Denver"))
          }

          return(tmp.df)
        })

        # Combine data.frames into one large data.frame
        # Column numbers must be the same length for each dataframe in order to merge using this function
        # Not sure why but there was a [-1] at the end of this function at one point
        rawSPIN.df <- data.table::rbindlist(data.ls, use.names = T, fill = T)

        # Round time variable to 0.1 seconds
        rawSPIN.df <- rawSPIN.df %>%
          mutate(`Time (sec)` = round(`Time (sec)`, 1)) %>%
          mutate(`Local Time` = sort(`Local Time`))

        duplicate.check <- length(which(duplicated(rawSPIN.df$`Local Time`)))

        if (duplicate.check != 0){
          stop("Times are not discrete")
        }

        # This is a checker to see if the time series is monotonically increasing
        if (all(diff(rawSPIN.df$`Local Time`) >= 0) != T){
          stop(paste0("Error - SPIN date times are not properly generated", " ", n))
        }

        # Remove potential missing time values
        rawSPIN.df <- rawSPIN.df[!is.na(rawSPIN.df$`Local Time`)]

        # Set any NaN's to NA
        rawSPIN.df[is.na(rawSPIN.df)] <- NA

        # Calculate the duration of the experiment from the first time logged
        # This is different than elapsed time since it is recorded differently somehow. DMT doesn't say
        rawSPIN.df <- rawSPIN.df %>%
          mutate(`Duration (s)` = round(as.numeric(`Local Time` - `Local Time`[1]))) %>%
          select(`Local Time`, `Duration (s)`, everything())

        #' Check if there is real data for the day
        #' If the compressor RPM is strictly zero, the instrument never entered
        #' the startup sequence
        if (max(rawSPIN.df$`Warm Cmpr (RPM)`) == 0) {
          # Print duration of loop and print date/status
          print(paste0(date.c, ' Empty, Elapsed Time (s) ', round((
            proc.time() - pct)[3], 2)))
          next
        }

        # Remove temporary variables
        rm(data.ls, files.spin)
      }

      # ---------------------------------------------------------------------- #
      ##### SECTION: Data Wrangling - Log File #####
      #' This is information generated by SPIN software to record status and
      #' any operations occurring with the instrument
      #' Some errors reported are ironically errors themselves in that they do not
      #' accurately reflect conditions of the instrument

      {
        if (length(files.log) != 0){

          print("SECTION: Data Wrangling - Log File")

          # Read in data using newline character as a separator
          # SPIN uses "\n" as a delimiter in log files
          tmp <- fread(files.log, fill = TRUE, sep = "\n")

          # Using the separate function to split the character column into two
          # This uses the regex pattern for a 24 time with a space after
          log.df <- tidyr::separate(tmp, col = 1, into = c("Time", "Message"),
                                    sep = "(?<=\\d{2}:\\d{2}:\\d{2}) "
          )

          # Extract date string from log file name
          log.df$Date <- as.Date(ymd(str_extract(files.log, "\\d{8}")))

          # Use run length encoding to find when a time shift across midnight occurs
          # Currently looks for times at which a 6 hour negative time shift happens
          switch.point = which(diff(rle(as.numeric(hms(log.df$Time)))$values) < -6*3600)
          if (length(switch.point) != 0){
            midnight.change = cumsum(rle(as.numeric(hms(log.df$Time)))$lengths)[switch.point] + 1

            log.df$Date[midnight.change:nrow(log.df)] <- as.Date(log.df$Date[1]) + 1
          }

          # Create a POSIX.ct object for further use
          log.df <- log.df %>%
            mutate(`Local Time` = as.POSIXct(paste0(log.df$Date, " ", log.df$Time),
                                             tz = "America/Denver"
            ))

          dataLOG.df <- log.df %>%
            select(-c(`Time`, `Date`)) %>%
            select(`Local Time`, `Message`) %>%
            arrange(`Local Time`) %>%
            mutate(`Date` = lubridate::as_date(`Local Time`), .before = everything())

          # Determine start times of specific sequences
          # If the names of the sequences have changed in the software this won't catch them correctly
          sequence.start.tm = dataLOG.df$`Local Time`[str_which(dataLOG.df$Message, "Starting Sequence Startup IV.")]
          sequence.evrmp.tm = dataLOG.df$`Local Time`[str_which(dataLOG.df$Message, "Starting Sequence Ramp with Evap Ramp II.")]
          sequence.icing.tm = dataLOG.df$`Local Time`[str_which(dataLOG.df$Message, "Starting Sequence Fill & Empty Chamber II.")]
          sequence.ramps.tm = dataLOG.df$`Local Time`[str_which(dataLOG.df$Message, "Starting Sequence Update SPs from Ramp.")]
          stopping.start.tm = dataLOG.df$`Local Time`[str_which(dataLOG.df$Message, "Stopping Sequence Startup IV.")]

          rm(tmp, log.df, files.log)
        } else {
          print(paste0("No log file detected"))

          dataLOG.df <- NULL
        }

        # This is a checker to see if the time series is monotonically increasing
        if (all(diff(dataLOG.df$`Local Time`) >= 0) != T){
          stop(paste0("Error - Log date times are not properly generated", " ", n))
        }
      }

      # ---------------------------------------------------------------------- #
      ##### SECTION: Data Wrangling - Particle by Particle #####
      #'

      {
        print("SECTION: Data Wrangling - Particle by Particle")

        # Loop through list of files and create a dataframe of each
        PbP.ls <- lapply(files.PbP, function(x) {
          # Use data.table's fread to read in data
          # Fastest method in reading CSV's and least memory consuming
          # fill MUST equal FALSE with SPIN files
          tmp.df <- data.table::fread(
            paste0(x),
            na.strings = c("", "NA", "NaN"),
            fill = TRUE,
            strip.white = TRUE,
            stringsAsFactors = FALSE
          )

          setnames(
            tmp.df,
            old = c("End ms", "S", "IPT (msec)"),
            new = c("End (ms)", "S1", "IPT (ms)")
          )

          # Extract date and time from the filename using regular expression
          # (?<=...) is a positive-look behind assertion
          # SPIN files always follow this pattern
          tmp.date = str_extract(x, "(?<=SPIN_PbP\\d{3}\\S?)\\d{8}")
          tmp.time = str_extract(x, "(?<=SPIN_PbP\\d{3}\\S?\\d{8})\\d{6}")

          # Place the filename in the data as an identifier
          # When times are generated these should be strictly increasing
          tmp.df <- tmp.df %>%
            mutate(`Filename` = paste0(tmp.date, tmp.time))

          return(tmp.df)
        })

        rawPBP.df <- data.table::rbindlist(PbP.ls, use.names = T, fill = T)

        # Only perform if there is actually data
        if (nrow(rawPBP.df) != 0){

          # Round time variable to 0.1 seconds
          rawPBP.df <- rawPBP.df %>%
            mutate(`Time (sec)` = round(`Time Stamp`, 1))

          # Retrieve SPIN times corresponding to time stamps
          # According to DMT, the PbP files should have matching timestamps as SPIN
          tmp.df <- rawSPIN.df %>%
            select(`Local Time`, `Time (sec)`)

          # Add POSIXct time to PBP data for splitting later
          dataPBP.df <- right_join(tmp.df, rawPBP.df, by = "Time (sec)") %>%
            mutate(`Time UTC` = lubridate::with_tz(`Local Time`, tzone = "UTC"), .after = `Local Time`) %>%
            mutate(`Date` = as.character(lubridate::as_date(`Local Time`)), .before = everything())

          # This is a checker to see if the time series is monotonically increasing
          if (all(diff(dataPBP.df$`Local Time`) >= 0) != T){
            stop(paste0("Error - Log date times are not properly generated", " ", n))
          }
        } else {
          dataPBP.df <- rawPBP.df
        }

        rm(PbP.ls, files.PbP)
      }

      # ---------------------------------------------------------------------- #
      ##### SECTION: Subset SPIN Data #####
      #' Combine all the SPIN data from the data list produced above
      #' Generate a secondary time column showing the duration of the experiment
      #' Subset data into groups that are later used in plotting
      #'

      {
        print("SECTION: Subset SPIN Data")

        # Create a string of the column names of the main dataframe
        tmp.nm = colnames(rawSPIN.df)

        # Indices of parameters that are not used
        indices.remove <- which(str_detect(tmp.nm, "Analog"))
        indices.remove <- append(indices.remove, which(str_detect(tmp.nm, "spare")))

        # Indices of binned particle data
        indices.bins <- which(str_detect(tmp.nm, "Bin\\s\\d+"))

        #  Indices of inlet valve status
        indices.inletvalve <- which(str_detect(tmp.nm, "Inlet Valve"))

        # Indices of relevant lamina variables
        indices.Lamina <- which(str_detect(tmp.nm, "SS%"))
        indices.Lamina <- append(indices.Lamina, which(str_detect(tmp.nm, "Lam")))
        indices.Lamina <- append(indices.Lamina, which(str_detect(tmp.nm, "\\(mm\\)")))

        # Indices of temperature setpoints
        indices.SP <- which(str_detect(tmp.nm, "SP"))

        # Indices of refrigerant variables
        indices.refrigerant <- which(str_detect(tmp.nm, "ColdTopRefr"))
        indices.refrigerant <- append(indices.refrigerant, which(str_detect(tmp.nm, "ColdMidRefr")))
        indices.refrigerant <- append(indices.refrigerant, which(str_detect(tmp.nm, "ColdBotRefr")))
        indices.refrigerant <- append(indices.refrigerant, which(str_detect(tmp.nm, "WarmTopRefr")))
        indices.refrigerant <- append(indices.refrigerant, which(str_detect(tmp.nm, "WarmBotRefr")))

        # Indices of flow rates
        indices.flows <- which(str_detect(tmp.nm, "Flow"))

        # Indices of pressure
        indices.pressure <- which(str_detect(tmp.nm, "Pressure"))

        # Indices of compressor values
        indices.cmpr <- which(str_detect(tmp.nm, "(?=.(psi))"))
        indices.cmpr <- append(indices.cmpr, which(str_detect(tmp.nm, ".*RPM")))
        indices.cmpr <- append(indices.cmpr, which(str_detect(tmp.nm, "Fault")))
        indices.cmpr <- append(indices.cmpr, which(str_detect(tmp.nm, "Cold.*(1st).*\\(C\\)")))
        indices.cmpr <- append(indices.cmpr, which(str_detect(tmp.nm, "Cold.*(2nd).*\\(C\\)")))
        indices.cmpr <- append(indices.cmpr, which(str_detect(tmp.nm, "Warm.*Cmpr.*\\(C\\)")))

        # Thermocouple Values
        indices.coldTCs <- which(str_detect(tmp.nm, "C\\d{0,2}\\(\\w\\)"))
        indices.warmTCs <- which(str_detect(tmp.nm, "W\\d{0,2}\\(\\w\\)"))

        # Heaters
        indices.ColdHeat <- which(str_detect(tmp.nm, "CH\\d+"))
        indices.ColdHeat <- append(indices.ColdHeat, which(str_detect(tmp.nm, "HeatCold_\\d+")))
        indices.WarmHeat <- which(str_detect(tmp.nm, "WH\\d+"))
        indices.WarmHeat <- append(indices.WarmHeat, which(str_detect(tmp.nm, "HeatWarm_\\d+")))

        # Evaporation

        # OPC Status
        indices.OPC <- which(str_detect(tmp.nm, "bandwidth"))
        indices.OPC <- append(indices.OPC, which(str_detect(tmp.nm, "baseline")))
        indices.OPC <- append(indices.OPC, which(str_detect(tmp.nm, "oversize")))
        indices.OPC <- append(indices.OPC, which(str_detect(tmp.nm, "Laser")))
        indices.OPC <- append(indices.OPC, which(str_detect(tmp.nm, "transit")))
      }

      # ---------------------------------------------------------------------- #
      ##### SUBSECTION: Flows and Chamber Pressure #####
      #'
      {
        dataFLOW.df <- rawSPIN.df %>%
          select(`Local Time`, all_of(c(indices.flows, indices.pressure)))

        # Rename the variables to a more friendly format
        tmp.nm <- str_squish(gsub('([[:upper:]])', ' \\1', colnames(dataFLOW.df)))
        tmp.nm = str_replace(tmp.nm, "S P", "SP")
        tmp.nm = str_replace(tmp.nm, " S L P M", "SLPM")
        tmp.nm = str_replace(tmp.nm, "slpm", "SLPM")
        tmp.nm = str_replace(tmp.nm, " L P M", "LPM")
        tmp.nm = str_replace(tmp.nm, "Flow\\(", "Flow (")
        tmp.nm = str_replace(tmp.nm, "\\(mbar", " (mbar")
        tmp.nm = str_replace(tmp.nm, "Flow%", "Flow %")
        tmp.nm = str_replace(tmp.nm, "Vol ", "Volumetric")
        tmp.nm = str_replace(tmp.nm, "cFlow", "c Flow")

        # Change column names
        setnames(dataFLOW.df, new = tmp.nm)
        setnames(dataFLOW.df, old = "Volumetric Flow (LPM)", new = "Sheath Volumetric Flow (LPM)")

        # Calculate total flow variable for volumetric flow
        dataFLOW.df <- dataFLOW.df %>%
          mutate(`Total Flow (LPM)` = `Sample Volumetric Flow (LPM)` + `Sheath Volumetric Flow (LPM)`)

        # Remove two variables that appear to be constants and reorder
        dataFLOW.df <- dataFLOW.df %>%
          select(c(colnames(dataFLOW.df)[1], sort(colnames(dataFLOW.df)[-1]))) %>%
          select(!c(`Closing Pressure`, `Sheath Flow Memory`))

        rm(indices.flows, indices.pressure, tmp.nm)
      }

      # ---------------------------------------------------------------------- #
      ##### SUBSECTION: Bin Data #####
      #' Binned aerosol data

      {
        # Midpoint diameter: NANOMETERS
        bins.nm = c(
          550,
          650,
          750,
          850,
          950,
          1500,
          2500,
          3500,
          4500,
          5500,
          6500,
          7500,
          8500,
          9500,
          10500,
          11500,
          12500,
          13500,
          14500
        )

        # Midpoint diameter: MICROMETERS
        bins.nm = bins.nm / 1000

        # Subset bins
        # Drop the first as it is not defined by a lower limit
        dataBINS.df <- rawSPIN.df %>%
          select(all_of(indices.bins)[-1]) %>%
          setnames(., new = paste0(bins.nm))

        # Append the time column to the front of the data.frame
        # Add a column indicated the units of the binned aerosol data
        dataBINS.df <- dataBINS.df %>%
          mutate(`Local Time` = rawSPIN.df$`Local Time`, .before = everything()) %>%
          mutate(`Units` = "#/s")

        rm(indices.bins, bins.nm)
      }

      # ---------------------------------------------------------------------- #
      ##### SUBSECTION: Inlet Filter Data #####
      #'

      {
        # 0 means closed
        # 1 means open
        dataFLTR.df <- rawSPIN.df %>%
          select(`Local Time`, all_of(indices.inletvalve)) %>%
          mutate(`Inlet Valve Inlet/Zero` = as.numeric(`Inlet Valve Inlet/Zero`))

        setnames(dataFLTR.df, old = "Inlet Valve Inlet/Zero", new = "Inlet Filter ON")

        rm(indices.inletvalve)
      }

      # ---------------------------------------------------------------------- #
      ##### SUBSECTION: Lamina Conditions #####
      #' Supersaturation and humidity
      #' Lamina temperature

      {
        dataLAM.df <- rawSPIN.df %>%
          select(`Local Time`, all_of(indices.Lamina))

        setnames(dataLAM.df,
                 old = colnames(dataLAM.df)[-1],
                 new = c("Lamina SS% Liquid",
                         "Lamina SS% Ice",
                         "Lamina Temp (C)",
                         "Lamina Centerline Ratio",
                         "Lamina Centerline Distance from Warm Wall (mm)",
                         "Lamina Thickness (mm)",
                         "Lamina Saturation Pressure Liquid/Ice",
                         "Nearside Distance (mm)",
                         "Farside Distance (mm)"
                 ))

        # Final processing
        # Convert to saturation ratio (S) rather than SS%
        # SS% in this case means percent above 100%
        # SS% = (S-1)*100
        # Also add dewpoint value
        dataLAM.df <- dataLAM.df %>%
          mutate(`Lamina S Liquid` = `Lamina SS% Liquid`/100 + 1) %>%
          mutate(`Lamina S Ice` = `Lamina SS% Ice`/100 + 1) %>%
          mutate("Dewpoint" = rawSPIN.df$`Dew Point (C)`)

        rm(indices.Lamina)
      }

      # ---------------------------------------------------------------------- #
      ##### SUBSECTION: Thermcouples #####
      #' Thermocouple variables

      {
        # Subset set points
        dataTCSP.df <- rawSPIN.df %>%
          select(all_of(c(indices.SP))) %>%
          select(`Warm_SP`, `Cold_SP`, `Evap_SP`) %>%
          setnames(., new = c("Warm Wall SP", "Cold Wall SP", "Evaporation SP")) %>%
          mutate(`CWSP Noise` = append(0, abs(diff(`Cold Wall SP`)))) %>%
          mutate(`WWSP Noise` = append(0, abs(diff(`Warm Wall SP`)))) %>%
          mutate("Local Time" = rawSPIN.df$`Local Time`, .before = everything())

        # Remove erroneous spikes from interuptions in the program
        dataTCSP.df$`Cold Wall SP`[which(dataTCSP.df$`CWSP Noise` > 22)] <- NA
        dataTCSP.df$`Warm Wall SP`[which(dataTCSP.df$`WWSP Noise` > 22)] <- NA

        # Drop temporary columns
        dataTCSP.df <- dataTCSP.df %>%
          select(-c(`CWSP Noise`, `WWSP Noise`))

        # Evaporation setpoint prior to the experiment time isn't a valid variable
        dataTCSP.df[which(rawSPIN.df$`Local Time` <= first(sequence.ramps.tm)), "Evaporation SP"] <- NA

        # Subset cold wall thermocouples
        dataCTC.df <- rawSPIN.df %>%
          select(all_of(c(indices.coldTCs))) %>%
          setnames(., new = c(str_remove(colnames(.), "\\(C\\)"))) %>%
          mutate("Local Time" = rawSPIN.df$`Local Time`, .before = everything())

        # Subset cold wall thermocouples
        dataWTC.df <- rawSPIN.df %>%
          select(all_of(c(indices.warmTCs))) %>%
          setnames(., new = c(str_remove(colnames(.), "\\(C\\)"))) %>%
          mutate("Local Time" = rawSPIN.df$`Local Time`,
                 .before = everything())

        rm(indices.coldTCs, indices.warmTCs)
      }

      # ---------------------------------------------------------------------- #
      ##### SUBSECTION: Compressor Data #####
      #' Compressor variables

      {
        #' Subset from rawSPIN.df
        #' Force to numeric
        dataCMPR.df <- rawSPIN.df %>%
          select(all_of(indices.cmpr)) %>%
          as.matrix() %>%
          as.data.frame() %>%
          mutate(`Local Time` = rawSPIN.df$`Local Time`, .before = everything()) %>%
          mutate(`Duration(s)` = rawSPIN.df$`Duration(s)`, .before = everything())

        # Rename the variables to a more friendly format
        tmp.nm = str_replace(colnames(dataCMPR.df), "Temp", "Temperature")
        tmp.nm = str_remove_all(tmp.nm, "Cmpr ")
        tmp.nm = str_remove_all(tmp.nm, "st")
        tmp.nm = str_remove_all(tmp.nm, "nd")
        tmp.nm = str_replace(tmp.nm, "Pres", "Pressure")

        setnames(dataCMPR.df, new = tmp.nm)

        # Add power supply current data
        dataCMPR.df$`Power Supply Current` <- rawSPIN.df$`Power Supply Curr (A)`

        # Remove unrealistic CMPR Out Temperature
        dataCMPR.df[dataCMPR.df < -273] <- NA
        dataCMPR.df$`Cold 2 Out Temperature (C)`[dataCMPR.df$`Cold 2 Out Temperature (C)` > 200] <- NA

        # Reorder columns
        dataCMPR.df <- dataCMPR.df %>%
          select(c(colnames(dataCMPR.df)[1], sort(colnames(dataCMPR.df)[-1])))

        rm(indices.cmpr, tmp.nm)
      }

      # ---------------------------------------------------------------------- #
      ##### SUBSECTION: Cooling System #####
      #'

      {
        dataREFR.df <- rawSPIN.df %>%
          select(`Local Time`, all_of(c(indices.refrigerant)))

        # Rename the variables to a more friendly format
        tmp.nm <- str_squish(gsub('([[:upper:]])', ' \\1', colnames(dataREFR.df)))
        tmp.nm = str_remove(tmp.nm, "_")
        tmp.nm = str_remove(tmp.nm, "Refr")
        tmp.nm = str_replace(tmp.nm, "S P", "SP")
        tmp.nm = str_replace(tmp.nm, "P V", "PV")
        tmp.nm = str_squish(tmp.nm)

        # Change column names
        setnames(dataREFR.df, new = tmp.nm)

        # Reorder dataframe
        dataREFR.df <- dataREFR.df %>%
          select(c(colnames(dataREFR.df)[1], sort(colnames(dataREFR.df)[-1])))

        dataCH.df <- rawSPIN.df %>%
          select(`Local Time`, all_of(c(indices.ColdHeat)))

        # Rename the variables to a more friendly format
        tmp.nm = str_replace(colnames(dataCH.df), "_", " ")
        tmp.nm = str_replace(tmp.nm, "OnTime", "On Time")
        tmp.nm = str_replace(tmp.nm, "OffTime", "Off Time")
        tmp.nm = str_replace(tmp.nm, "HeatCold", "Heat Cold")

        # Change column names
        setnames(dataCH.df, new = tmp.nm)

        # Reorder dataframe
        dataCH.df <- dataCH.df %>%
          select(c(colnames(dataCH.df)[1], sort(colnames(dataCH.df)[-1])))

        dataWH.df <- rawSPIN.df %>%
          select(`Local Time`, all_of(c(indices.WarmHeat)))

        # Rename the variables to a more friendly format
        tmp.nm = str_replace(colnames(dataWH.df), "_", " ")
        tmp.nm = str_replace(tmp.nm, "OnTime", "On Time")
        tmp.nm = str_replace(tmp.nm, "OffTime", "Off Time")
        tmp.nm = str_replace(tmp.nm, "HeatWarm", "Heat Warm")

        # Change column names
        setnames(dataWH.df, new = tmp.nm)

        # Reorder dataframe
        dataWH.df <- dataWH.df %>%
          select(c(colnames(dataWH.df)[1], sort(colnames(dataWH.df)[-1])))

        rm(indices.refrigerant, indices.ColdHeat, indices.WarmHeat)
      }

      # ---------------------------------------------------------------------- #
      ##### SUBSECTION: OPC System Variables #####
      #'

      {
        dataOPC.df <- rawSPIN.df %>%
          select(`Local Time`, all_of(c(indices.OPC)))

        # Rename the variables to a more friendly format
        tmp.nm = str_replace(colnames(dataOPC.df), "low", "Low")
        tmp.nm = str_replace(tmp.nm, "hi", "High")
        tmp.nm = str_replace(tmp.nm, "sizer", "Sizer")
        tmp.nm = str_replace(tmp.nm, "baseline", "Baseline")
        tmp.nm = str_replace(tmp.nm, "bandwidth", "Bandwidth")
        tmp.nm = str_replace(tmp.nm, "oversize", "Oversize")
        tmp.nm = str_replace(tmp.nm, "total transit time", "Total Transit Time")

        # Change column names
        setnames(dataOPC.df, new = tmp.nm)

        # Reorder dataframe
        dataOPC.df <- dataOPC.df %>%
          select(c(colnames(dataOPC.df)[1], sort(colnames(dataOPC.df)[-1])))

        rm(indices.OPC, tmp.nm)
      }

      # ---------------------------------------------------------------------- #
      ##### SECTION: Aggregate Data #####
      #' Combine datasets from above into one dataframe
      #' Create a finalized dataset for any experiments
      #'

      {
        print("SECTION: Aggregate Data")

        # Combine subsets
        dataSPIN.df <-
          cbind(
            dataBINS.df,
            dataFLTR.df,
            dataFLOW.df,
            dataLAM.df,
            dataOPC.df,
            dataTCSP.df,
            dataCTC.df,
            dataWTC.df,
            dataCMPR.df,
            dataREFR.df,
            dataCH.df,
            dataWH.df,
            rawSPIN.df[, 1:7]
          )

        dataSPIN.df <- subset(dataSPIN.df, select = !duplicated(names(dataSPIN.df))) %>%
          select(c(`Local Time`, `Time (sec)`, `Duration (s)`,
                   `Lamina S Ice`, `Lamina S Liquid`), everything())

        # Add a UTC time column
        dataSPIN.df <- dataSPIN.df %>%
          mutate(`Time UTC` = lubridate::with_tz(`Local Time`, tzone = "UTC"), .after = `Local Time`) %>%
          mutate(`Date` = as.character(lubridate::as_date(`Local Time`)), .before = everything()) %>%
          select(!`Timestamp`)

        rm(dataBINS.df, dataFLTR.df, dataFLOW.df, dataLAM.df, dataOPC.df,
           dataTCSP.df, dataCTC.df, dataWTC.df, dataCMPR.df, dataREFR.df,
           dataCH.df, dataWH.df)
      }

      # ---------------------------------------------------------------------- #
      ##### SECTION: Export Data #####
      #'
      {
        export.ls <- list()

        export.ls[[1]] <- dataSPIN.df %>%
          group_split(`Date`)

        if (nrow(dataPBP.df) != 0){
          export.ls[[2]] <- dataPBP.df %>%
            group_split(`Date`)
        } else {
          export.ls[[2]] <- list()
        }

        if (!is.null(dataLOG.df)){
          export.ls[[3]] <- dataLOG.df %>%
            group_split(`Date`)
        } else {
          export.ls[[3]] <- list()
        }

        rm(dataSPIN.df, dataPBP.df, dataLOG.df)

        for (i in 1:length(export.ls)){

          tmp.ls <- export.ls[[i]]

          if (length(tmp.ls) > 0){

            for (j in 1:length(tmp.ls)){

              # Merge back for exporting
              export.df <- tmp.ls[[j]]

              # Create a path to export plots with
              export.data.path = paste0(export.data, str_remove_all(unique(export.df$Date), '-'))

              # Check if export path exists
              # If it does not, create it
              if (!dir.exists(export.data.path)) {
                # Create a dated directory to send plots to
                dir.create(export.data.path, mode = "777")
              }

              if (i == 1){
                type = "SPIN001_SPIN"
              }

              if (i == 2){
                type = "SPIN001_PBP"
              }

              if (i == 3){
                type = "SPIN001_LOG"
              }

              filename = paste0(export.data.path, "/", type, "_Level0_", format(first(export.df$`Local Time`), format="%Y%m%d%H%M%S"), ".csv")

              print(paste0("Exporting Level 0 Data: ", type, "_Level0_", format(first(export.df$`Local Time`), format="%Y%m%d%H%M%S"), ".csv"))

              # Save data using data.table::fwrite
              data.table::fwrite(export.df, file = filename, showProgress = T)
            }
          }
        }
      }
    }
  }
}

