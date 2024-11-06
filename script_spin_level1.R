##### SCRIPT: SPectrometer for Ice Nucleation (SPIN) Level 1 #####
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
  library(latex2exp)

  # Plotting parameters
  resolution.dpi = 400
  font.family = "Helvetica"

  setwd("~/Library/CloudStorage/Box-Box/SPIN1-FIN03")
  work.dir <- getwd()

  #' @import
  #' Specify import directory
  #' FIN03
  import.spin = paste0(work.dir, "/export/level0/")

  #' @export
  #' Specify where to send plots and text files
  export.data = paste0(work.dir, "/export/level1/")
  export.plot = paste0(work.dir, "/export/plots/")

  #' @importFrom
  source(paste0("~/Documents/GitHub/functions/", "functions_microphysics.R"))
  source(paste0("~/Documents/GitHub/functions/", "functions_spin_plotting.R"))
}

# ---------------------------------------------------------------------------- #
##### SCRIPT: Start Loop #####
#'

{
  # Directories
  spin.dirs <- list.dirs(path = import.spin,
                         recursive = FALSE,
                         full.names = FALSE)

  # Pull POSIX style dates from directory names
  spin.dates <- as.Date(spin.dirs, format = '%Y%m%d')

  # Generate file paths
  spin.path <- paste0(import.spin, spin.dirs)

  for (n in 1:length(spin.path)){

    {
      # List all SPIN files exported
      files.all <- list.files(
        path = spin.path[n],
        recursive = TRUE,
        full.names = TRUE,
        pattern = '*.csv'
      )

      # Subset out the files that match PbP, SPIN, and LOG data
      files.spin <- files.all[str_which(files.all, pattern = "LOG", negate = T)]
      files.spin <- files.spin[str_which(files.spin, pattern = "PBP", negate = T)]
      files.pbp <- files.all[str_which(files.all, pattern = "PBP")]
      files.log <- files.all[str_which(files.all, pattern = "LOG")]

      # ---------------------------------------------------------------------- #
      ##### SECTION: SPIN Data #####
      #' Data for SPIN is split into three files
      #' The first is the main file that contains housekeeping data and binned counts
      #' The second is the log file which is used to determine when certain sequences have started
      #' The third is the particle by particle data recorded by the OPC
      #'

      {
        # Loop through all SPIN files and create a dataframe
        data.ls <- lapply(files.spin, function(x){

          # Read in data, no arguments needed as it is already formatted in Level 0
          tmp.df <- fread(x)

          # Correct the timezones. Fread reads times into UTC. with_tz
          tmp.df <- tmp.df  %>%
            mutate(`Local Time` = lubridate::with_tz(`Local Time`, tzone = "America/Denver"))

          return(tmp.df)
        })

        rawSPIN.df <- rbindlist(data.ls, use.names = T, fill = T) %>%
          relocate(`Lamina Temp (C)`, .before = `Lamina S Ice`)

        # For these older SPIN files `Inlet Filter ON` is inverted 0 means filter is ON
        rawSPIN.df$`Inlet Filter ON` <- abs(rawSPIN.df$`Inlet Filter ON` - 1)
      }

      # ------------------------------------------------------------------------ #
      ##### SECTION: LOG Data #####
      #' Data for SPIN is split into three files
      #' The first is the main file that contains housekeeping data and binned counts
      #' The second is the log file which is used to determine when certain sequences have started
      #' The third is the particle by particle data recorded by the OPC
      #'

      {
        data.ls <- lapply(files.log, function(x){

          tmp.df <- fread(x)

          return(tmp.df)
        })

        rawLOG.df <- rbindlist(data.ls, use.names = T, fill = T)

        LOG.break <<- FALSE

        tryCatch(if (length(rawLOG.df) < 1) {stop("Error")},
                 error = function(e) {LOG.break <<- TRUE})

        # Stop code from continuing if there are no ramps detected which means the dataset will be empty and continue to throw errors
        if (LOG.break) {next}

        dataLOG.df <- rawLOG.df %>%
          mutate(`Local Time` = lubridate::with_tz(`Local Time`, tzone = "America/Denver"))

        # Determine start times of specific sequences
        # If the names of the sequences have changed in the software this won't catch them correctly
        sequence.start.tm = dataLOG.df$`Local Time`[str_which(dataLOG.df$Message, "Starting Sequence Startup IV.")]
        sequence.evrmp.tm = dataLOG.df$`Local Time`[str_which(dataLOG.df$Message, "Starting Sequence Ramp with Evap Ramp II.")]
        sequence.icing.tm = dataLOG.df$`Local Time`[str_which(dataLOG.df$Message, "Starting Sequence Fill & Empty Chamber II.")]
        stopping.start.tm = dataLOG.df$`Local Time`[str_which(dataLOG.df$Message, "Stopping Sequence Startup IV.")]
      }

      # ---------------------------------------------------------------------- #
      ##### SECTION: PBP Data #####

      {
        # Loop through all SPIN files and create a dataframe
        data.ls <- lapply(files.pbp, function(x){

          # Read in data, no arguments needed as it is already formatted in Level 0
          tmp.df <- fread(x)

          # Correct the timezones. Fread reads times into UTC. with_tz
          tmp.df <- tmp.df %>%
            mutate(`Local Time` = lubridate::with_tz(`Local Time`, tzone = "America/Denver"))

          return(tmp.df)
        })

        tmp.df <- rbindlist(data.ls, use.names = T, fill = T)

        PbP.break <<- FALSE

        tryCatch(if (length(tmp.df) < 1) {stop("Error")},
                 error = function(e) {PbP.break <<- TRUE})

        # Stop code from continuing if there are no ramps detected which means the dataset will be empty and continue to throw errors
        if (PbP.break) {next}

        times.df <- tmp.df %>%
          distinct(`Time Stamp`, .keep_all = T) %>%
          select(c(`Date`:`Time Stamp`))

        # Save numeric values
        rawPBP.df <- tmp.df %>%
          select(!c(`Date`:`Time (sec)`))

        # Correct the stitching of the gain stages?
        # Not sure what this means... included in the original Matlab code
        rawPBP.df$S1[rawPBP.df$S1 > 14930] = rawPBP.df$S1[rawPBP.df$S1 > 14930] - 2000
        rawPBP.df$P1[rawPBP.df$P1 > 13570] = rawPBP.df$P1[rawPBP.df$P1 > 13570] - 2000
        rawPBP.df$P2[rawPBP.df$P2 > 14610] = rawPBP.df$P2[rawPBP.df$P2 > 14610] - 3000

        # Benchmarking
        x1 <- nrow(rawPBP.df)

        # Very very very fast method to do averaging for this dataset
        rawPBP.df <- rawPBP.df %>%
          collapse::fgroup_by(`Time Stamp`) %>%
          collapse::fmean()

        # Benchmarking
        x2 <- nrow(rawPBP.df)

        # Remove problematic values
        rawPBP.df <- na.omit(rawPBP.df)
        rawPBP.df <- rawPBP.df[!is.infinite(rowSums(rawPBP.df)), ]

        print(paste0(x1, " Observations averaged to ", x2))

        # Calculate log values
        rawPBP.df$`Log P1` <- log10(rawPBP.df$P1)
        rawPBP.df$`Log P2` <- log10(rawPBP.df$P2)
        rawPBP.df$`Log S1` <- log10(rawPBP.df$S1)
        rawPBP.df$`Log Size` <- log10(rawPBP.df$Size)
        rawPBP.df$`Log S1/P1` <- log10(rawPBP.df$S1/rawPBP.df$P1)

        # Remove any problematic values caused by taking the log
        rawPBP.df <- rawPBP.df[!is.nan(rowSums(rawPBP.df)), ]
        rawPBP.df <- rawPBP.df[!is.infinite(rowSums(rawPBP.df)), ]
        rawPBP.df <- na.omit(rawPBP.df)

        # Move the times back to the dataframe
        dataPBP.df <- right_join(times.df, rawPBP.df, by = "Time Stamp")

        # Time stamps for collapsing data to times where counts are being recorded
        starting.laser.tm = first(dataPBP.df$`Local Time`)
        stopping.laser.tm = last(dataPBP.df$`Local Time`)

        rm(times.df, x1, x2, tmp.df)
      }

      # ---------------------------------------------------------------------- #
      ##### SECTION: Aggregate Data #####

      {
        dataPBP.df <- rawPBP.df %>%
          select(`Time Stamp`, `Log S1/P1`, `Log P1`, `Log P2`, `Log S1`, `Log Size`,
                 `S1`, `P1`, `P2`, `Size`)

        setnames(dataPBP.df, old = "Time Stamp", new = "Time (sec)")

        # Time stamps need to be averaged to the nearest hundreth second to merge properly
        rawSPIN.df$`Time (sec)` <- round(rawSPIN.df$`Time (sec)`, 1)
        dataPBP.df$`Time (sec)` <- round(dataPBP.df$`Time (sec)`, 1)

        # Perform joins
        # Right join is focused on only allowing observations from the first join the second if there is a match
        # This inherently filters data so that only points where the OPC was on are included
        dataSPIN.df <- right_join(rawSPIN.df, dataPBP.df, by = "Time (sec)") %>%
          relocate(c(`Log S1/P1`, `Log P1`, `Log P2`, `Log S1`, `Log Size`, `S1`, `P1`, `P2`, `Size`), .after = `Inlet Filter ON`)

        # Date value to pass to plotting functions
        date.c <- unique(lubridate::as_date(dataSPIN.df$`Local Time`))

        print(paste0(n, " - ", date.c))
      }

      # ---------------------------------------------------------------------- #
      ##### PLOTTING: Overview Plot #####
      #'

      # Overview plot
      {
        #' Plot Cold Wall SP
        #' Plot Warm Wall SP
        #' Plot Lamina SSice
        #' Plot Lamina SSliquid
        #' Plot Lamina temperature
        #' Plot Inlet Valve Position
        #'
        title.main <- paste0("Spectrometer for Ice Nucleation (SPIN)")
        title.sub <- paste0("Experiment Overview: ", date.c)

        plot.filename = paste0(export.plot, date.c, '/', date.c, ' Overview', '.png')

        print(paste0("Plotting: ", plot.filename))

        overview.plot(title.main,
                      title.sub,
                      data = rawSPIN.df,
                      sequence.start.tm,
                      sequence.evrmp.tm,
                      sequence.icing.tm,
                      sequence.ramps.tm = NA,
                      stopping.start.tm,
                      plot.filename,
                      plot.width = 14,
                      plot.height = 8,
                      plot.resolution.dpi = 300)
      }

      # ---------------------------------------------------------------------- #
      ##### SECTION: Formatting Data #####

      # Filtering data to remove erroneous values
      # Calculating isothermal setpoints
      {
        #' Keep data when the sample flow is on
        #' Only keep positive sheath flow values
        #' Eliminate data that didn't have the cold2 compressor on
        #' Remove data without a lamina calculation
        #' Calculate rounded lamina temperatures to determine isothermal setpoints
        tmp.df <- dataSPIN.df %>%
          filter(`Sample Volumetric Flow (LPM)` > 0.5) %>%
          filter(`Sheath Volumetric Flow (LPM)` > 0) %>%
          filter(`Cold 2 (RPM)` != 0) %>%
          filter(!is.na(`Lamina Temp (C)`)) %>%
          mutate(`Mean Lamina` = zoo::rollapplyr(`Lamina Temp (C)`, 60*5, mean, fill = NA)) %>%
          mutate(`Lamina Breaks` = round(`Mean Lamina`/50, 1)*50, .after = `Lamina Temp (C)`)

        # This checks to see how many datapoints have a large difference in wall SP's
        # This corresponds to experiments
        exp.test = length(which(abs(tmp.df$`Cold Wall SP` - tmp.df$`Warm Wall SP`) > 5))

        exp.break <<- FALSE

        tryCatch(if (exp.test < 15*60) {
          stop("Did this work?")},
          error = function(e) {exp.break <<- TRUE})

        # Stop code from continuing if there are no ramps detected which means the dataset will be empty and continue to throw errors
        if (exp.break) {next}

        # 3 minute rolling standard deviation
        # Shorter windows will lead to more dataloss
        T.sd = zoo::rollapplyr(tmp.df$`Mean Lamina`, 60*3, sd, fill = NA)

        # This is threshold for outliers but we are using it to identify when the lamina is changing to a new isothermal setpoint
        threshold.nm <- as.numeric(quantile(T.sd, na.rm = T)[4] + IQR(T.sd, na.rm = T)*2)

        # Remove break data corresponding to a major shift due to lamina ramps
        tmp.df$`Lamina Breaks`[which(T.sd >= threshold.nm)] <- NA

        # Correct lamina breaks
        {
          # While loop to remove datapoints that are slightly above the "breaks"
          while(length(which(diff(tmp.df$`Lamina Breaks`) > 0)) != 0){

            upper = which(diff(tmp.df$`Lamina Breaks`) > 0) + 1
            target = which(diff(tmp.df$`Lamina Breaks`) > 0)

            tmp.df$`Lamina Breaks`[upper] <- tmp.df$`Lamina Breaks`[target]
          }

          # While loop to remove datapoints that are slightly below the "breaks"
          while(length(which(diff(tmp.df$`Lamina Breaks`) < 0)) != 0){

            lower = which(diff(tmp.df$`Lamina Breaks`) < 0) + 1
            target = which(diff(tmp.df$`Lamina Breaks`) < 0)

            tmp.df$`Lamina Breaks`[lower] <- tmp.df$`Lamina Breaks`[target]
          }
        }

        # Subset
        dataALL.df <- tmp.df

        tmp.df <- dataALL.df %>%
          filter(!is.na(`Lamina Breaks`)) %>%
          select(`Lamina Breaks`, `Lamina S Liquid`, `Lamina Temp (C)`, `Sheath Volumetric Flow (LPM)`, `Sample Volumetric Flow (LPM)`)

        # Convert to long format for easier plotting
        tmp.long <- tmp.df %>%
          pivot_longer(
            cols = !"Lamina Breaks",
            names_to = "Variable",
            values_to = "Value"
          )

        plot.ls <- NULL
        for (i in 1:length(unique(tmp.long$Variable))){

          tmp <- tmp.long %>%
            filter(`Variable` == unique(tmp.long$Variable)[i])

          plot.ls[[i]] <- ggplot(data = tmp, aes(x = `Variable`, y = `Value`, group = `Variable`, fill = `Variable`)) +
            facet_wrap(.~`Lamina Breaks`) +
            stat_boxplot(geom = "errorbar",
                         width = 0.15,
                         color = 1) +  # Error bar color
            geom_boxplot(fill = 2,           # Box color
                         alpha = 1,        # Transparency
                         color = 1,          # Border color
                         outlier.colour = 2) +
            xlab(label = paste0(tmp$Variable[1])) +
            scale_y_continuous(n.breaks = 10) +
            theme(
              plot.title = element_text(),
              plot.subtitle = element_text(color = "gray25"),
              panel.background = element_rect(fill = "white"),
              panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.1),
              panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.1),
              panel.grid.minor = element_line(colour = "grey80", linewidth = 0.1),
              panel.border = element_rect(colour = "black", fill = NA),
              axis.title.y = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_line(linewidth = 0.5),
              axis.ticks.y = element_line(linewidth = 0.5),
              axis.ticks.length = unit(1, "mm"),
              plot.margin = unit(c(0.2, 1, 0.2, 0.1), "cm"),
              legend.title = element_blank(),
              panel.spacing = unit(0.5, "cm")
            ) +
            guides(colour = guide_legend(override.aes = list(size=2)))
        }

        all.gg <- ggarrange(plotlist = plot.ls)
        all.gg <- annotate_figure(all.gg, top = text_grob(paste0(date.c), face = "bold", size = 14))

        plot.filename = paste0(export.plot, date.c, '/', date.c, ' Chamber Statistics', '.png')

        print(paste0("Plotting: ", plot.filename))
        ggsave(
          plot.filename,
          all.gg,
          width = 8,
          height = 8,
          dpi = 600,
          units = "in",
          bg = "#ffffff"
        )
      }

      # ---------------------------------------------------------------------- #
      ##### SECTION: Particle Analysis Data #####

      # Plot particle time series
      plot.filename = paste0(export.plot, date.c, '/', date.c, ' Particle Counts', '.png')

      particle.levelplot.png(data = dataALL.df,
                             col.palette = NULL,
                             title = date.c,
                             subtitle = "All Data",
                             time.type = "POSIX",
                             time.variable = dataALL.df$`Local Time`,
                             plot.filename = plot.filename,
                             plot.width = 12,
                             plot.height = 6,
                             plot.resolution.dpi = 400,
                             optional.y = NULL)

      # LOW PASS FILTER
      filter.nm = c(3)

      for (j in 1:length(filter.nm)){

        filter.counts = filter.nm[j]

        # Identify bins
        bins.ix <- suppressWarnings(which(!is.na(as.numeric(colnames(dataALL.df)))))
        bins.nm <- colnames(dataALL.df)[bins.ix]

        # Considering ice as being 5 um or larger
        cutoff.bins.nm <- bins.nm[which(as.numeric(bins.nm) >= 5)]

        # Select only the bins
        dataBINS.df <- dataALL.df %>%
          select(all_of(cutoff.bins.nm))

        remove.c <- which(rowSums(dataBINS.df) > filter.counts)

        dataBINS.df[remove.c, ] <- NA

        # Unit conversion from counts to concentration
        convert.nm = (dataALL.df$`Sample Volumetric Flow (LPM)` + dataALL.df$`Sheath Volumetric Flow (LPM)`)/60

        # Axis 1 indicates perform rowwise
        # Convert to number density rather than counts
        # Uses volumetric flow of the sample flow
        # Theoretically only the OPC sees 1 LPM at the lamina position
        dataBINS.df <- sweep(dataBINS.df, 1, convert.nm, "/")

        # Merge number density data back
        tmp.df <- dataALL.df %>%
          select(!all_of(bins.nm)) %>%
          cbind(., dataBINS.df) %>%
          relocate(all_of(cutoff.bins.nm), .after = `Lamina S Liquid`) %>%
          mutate(`Units` = "n/L")

        # Calculate and subtract backgrounds
        {
          # Is there enough data to accurately make any statistical inference?
          # 10 minutes minimum of data should be without a filter
          filter.test = length(which(tmp.df$`Inlet Filter ON` == 0))/60

          filter.break <<- FALSE

          tryCatch(if (filter.test < 10) {
            stop("Did this work?")},
            error = function(e) {filter.break <<- TRUE})

          # Stop code from continuing
          if (filter.break) {next}

          CF = 1.4

          # Calculate total particles larger than cutoff
          # Apply depolarization filter
          tmp.df <- tmp.df %>%
            mutate(`Total Particles` = rowSums(select(., .dots = all_of(cutoff.bins.nm))), .after = `Lamina S Liquid`) %>%
            filter(`Log Size` > 3.5) %>%
            filter(`Log S1/P1` > -0.25) %>%
            mutate(`Total Particles` = `Total Particles`*CF)

          # Calculate total number of particles when filter is on
          tmp.df <- tmp.df %>%
            mutate(`Total Background` = if_else(`Inlet Filter ON` == 1, `Total Particles`, NA), .after = `Total Particles`)

          # Using a linear interpolation to fill in the time series than averaging gives a higher background ice concentration
          # This gives a more conservative approach in declaring ice counts
          tmp.df <- tmp.df %>%
            mutate(`Background Particles` = zoo::na.approx(`Total Background`, na.rm = F), .after = `Total Background`)

          # Calculate background particles using a 5 minute window rolling mean
          tmp.df <- tmp.df %>%
            mutate(`Background Particles` = zoo::rollapplyr(`Background Particles`, 5*60, mean, fill = 0))

          # Subtract number density of backgrounds
          tmp.df <- tmp.df %>%
            mutate(`INP` = `Total Particles` - `Background Particles`, .after = `Total Background`) %>%
            mutate(`INP` = if_else(`INP` > 0, `INP`, 0))

          # Set values when the filter is on to NA
          # Using NA here as it will artifically lower the actual number due to weighting
          tmp.df <- tmp.df %>%
            mutate(`INP` = if_else(`Inlet Filter ON` == 0, `INP`, NA))

          # Remove data points that exceed water saturation
          tmp.df <- tmp.df %>%
            filter(`Lamina S Liquid` < 1)
        }

        # Apply outlier correction to INP concentration
        {
          # 3 minute rolling standard deviation
          # Shorter windows will lead to more dataloss
          N.sd = sd(tmp.df$`INP`, na.rm = T)

          # This is threshold for outliers but we are using it to identify when the lamina is changing to a new isothermal setpoint
          upper.threshold.nm <- as.numeric(N.sd + IQR(tmp.df$`INP`, na.rm = T)*1.5)

          # Remove data points entirely after the filter otherwise they will sway statistics
          # Perform this based on the lamina break (experiment setting)
          tmp.df <- tmp.df %>%
            group_by(`Lamina Breaks`) %>%
            mutate(`INP Corrected` = if_else(`INP` < upper.threshold.nm, `INP`, NA)) %>%
            ungroup()
        }

        # Convert to long format for easier plotting
        df1.long <- tmp.df %>%
          select(c(`Local Time`, `Total Background`, `Background Particles`, `INP`, `INP Corrected`, `Total Particles`)) %>%
          pivot_longer(
            cols = !`Local Time`,
            names_to = "Variable",
            values_to = "Value"
          ) %>%
          mutate(`Variable` = factor(`Variable`, levels = c("INP Corrected", "INP", "Total Particles", "Total Background", "Background Particles")))

        particle.gg <- ggplot(df1.long, aes(x = `Local Time`, y = `Value`, color = `Variable`)) +
          geom_point(shape = 1, size = 0.5) +
          theme_minimal() +
          scale_y_continuous(n.breaks = 10) +
          ylab(label = TeX("$N_{INP}\\  L^{-1}$")) +
          labs(title = paste0(date.c), subtitle = "Particle Filtering") +
          theme(
            plot.title = element_text(),
            plot.subtitle = element_text(color = "gray25"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.1),
            panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.1),
            panel.grid.minor = element_line(colour = "grey80", linewidth = 0.1),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.text.y.right = element_blank(),
            axis.title.y.right = element_blank(),
            axis.ticks.x = element_line(linewidth = 0.5),
            axis.ticks.y = element_line(linewidth = 0.5),
            axis.ticks.length = unit(1, "mm"),
            plot.margin = unit(c(0.2, 1, 0.2, 0.1), "cm"),
            legend.title = element_blank()
          ) +
          guides(colour = guide_legend(override.aes = list(size=2)))

        plot.filename = paste0(export.plot, date.c, '/', date.c, ' Particle Summary - Low Pass Filter ', filter.counts, ' Counts.png')

        print(paste0("Plotting: ", plot.filename))

        ggsave(
          plot.filename,
          particle.gg,
          width = 8,
          height = 6,
          dpi = 400,
          units = "in",
          bg = "#ffffff"
        )

        # Convert to long format for easier plotting
        df2.long <- tmp.df %>%
          select(c(`Lamina Breaks`, `Total Particles`, `Total Background`, `INP`, `INP Corrected`,)) %>%
          filter(!is.na(`Lamina Breaks`)) %>%
          pivot_longer(
            cols = !"Lamina Breaks",
            names_to = "Variable",
            values_to = "Value"
          ) %>%
          mutate(`Variable` = factor(`Variable`, levels = c("INP Corrected", "INP", "Total Particles", "Total Background")))

        particle2.gg <- ggplot(data = df2.long, aes(x = `Variable`, y = `Value`, group = `Variable`)) +
          labs(title = paste0(date.c), subtitle = "Particle Statistics") +
          facet_wrap(.~`Lamina Breaks`) +
          stat_boxplot(geom = "errorbar", na.rm = T) +
          geom_boxplot(aes(fill = `Variable`),
                       outlier.colour = "red",
                       na.rm = T) +
          scale_y_log10(limits = c(1E-3, 1E4), n.breaks = 7, expand = c(0,0),
                        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
          scale_x_discrete(expand = c(0.25, 0.25)) +
          ylab(label = TeX("$N_{INP}\\  L^{-1}$")) +
          annotation_logticks(sides = "l") +
          theme(
            plot.title = element_text(),
            plot.subtitle = element_text(color = "gray25"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major.y = element_line(colour = "gray90", linewidth = 0.1),
            panel.grid.minor.y = element_line(colour = "gray90", linewidth = 0.1),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(linewidth = 0.5),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.length = unit(1, "mm"),
            plot.margin = unit(c(0.2, 1, 0.2, 0.1), "cm"),
            legend.title = element_blank(),
            legend.position = "bottom",
            panel.spacing = unit(0.5, "cm")
          ) +
          guides(colour = guide_legend(override.aes = list(size=2)))

        plot.filename = paste0(export.plot, date.c, '/', date.c, ' Particle Statistics - Low Pass Filter ', filter.counts, ' Counts.png')

        print(paste0("Plotting: ", plot.filename))

        ggsave(
          plot.filename,
          particle2.gg,
          width = 8,
          height = 6,
          dpi = 400,
          units = "in",
          bg = "#ffffff"
        )

        # ------------------------------------------------------------------------ #
        ##### SECTION: Export and Summarize #####

        tmp <- tmp.df %>%
          filter(!is.na(`Lamina Breaks`))

        summary.df <- tmp %>%
          group_by(`Lamina Breaks`) %>%
          summarize(
            `Date` = unique(tmp$Date),
            `CF` = CF,
            `Outlier Correction` = F,
            `Low Pass Filter (counts/s)` = filter.counts,
            `INP Mean` = mean(`INP`, na.rm = T),
            `INP N Obs` = length(which(!is.na(`INP`))),
            `INP Total` = round(sum(`INP`, na.rm = T), 0),
            `INP Error` = sqrt(`INP Mean`),
            `Background Mean` = mean(`Background Particles`, na.rm = T),
            `Background N` = length(which(!is.na(`Background Particles`))),
            `Background Error` = sqrt(`Background Mean`),
            `Total Error` = sqrt((`INP Error`^2) + (`Background Error`^2)),
            `Test Statistic` = 1.64*`Total Error`,
            `Confidence Level` = 0.95,
            `INP Units` = "n/L",
            `T Mean` = mean(`Lamina Temp (C)`, na.rm = T),
            `T Error` = plotrix::std.error(`Lamina Temp (C)`, na.rm = T),
            `T SD` = sd(`Lamina Temp (C)`, na.rm = T),
            `RH Mean` = mean(100*`Lamina S Liquid`, na.rm = T),
            `RH SD` = sd(100*`Lamina S Liquid`, na.rm = T),
            `RH Max` = max(100*`Lamina S Liquid`, na.rm = T),
            `RH Error` = max(abs(quantile(`Lamina S Liquid`, probs = c(0.05, 0.95)) - `RH Mean`/100)*100)) %>%
          mutate(`Lamina T Target` = `Lamina Breaks`, .after = `Date`) %>%
          select(!`Lamina Breaks`)

        # Round values
        summary.df <- summary.df %>%
          mutate(across(where(is.numeric), ~round(., 4))) %>%
          mutate(`Pass/Fail` = if_else(`INP Mean` > `Test Statistic`, "Pass", "Fail"), .before = `INP Mean`)

        loop.time = format(Sys.time(), "%y%m%d%H%M%S")

        filename = paste0(export.data, "/", "SPIN1_FIN03_", date.c, "_", loop.time, "_level1", ".csv")

        print(paste0("Exporting Level 1 Data: ", date.c, " Level 1"))

        # Save data using data.table::fwrite
        data.table::fwrite(summary.df, file = filename, showProgress = T, append = T)
      }
    }
  }
}
