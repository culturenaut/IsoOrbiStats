####### Generalized script to asses isotope ratios of ions of interest from Orbi measurements #############

#------------------
# created by Merve Öztoprak
# m.oeztoprak@gmail.com
#------------------

# this code calculates isotope ratios of isotopically substituted vs. non substituted ion species recorded via the Orbitrap MS
# works for direct elutions ("bulk" measurements of intact molecular ions) and molecular fragments of interest (MS2 via HCD fragmentation)

# functions are intended to be used sequentially

# Input:
# .iosX file extracted from individual measurement runs or combined isoX file for multiple measurements, identified by their file name
# isotopologs.tsv file which was also used to generate the .isox file
# .txt file with start and end times for the peak of interest in each measurement run
# tolerance factor = the amount of deviation allowed from the exact mass of a molecule given in Daltons

# Output:
# Isotope ratios calculated for multiple ions/fragments of interest with and without background subtraction and other data culling options

# Isotope ratios are calculated three ways:
# 1) as described in Eiler et al., 2017 Nio = (intensity/peakNoise)*(4.4 /1)*((120000/resolution)^0.5)*(microscans^0.5) CN=4.4 ref.resolution=120000, z=1
# 2) weighted Eiler calculation Nio_weighted = (Nio*intensity)/sum(intensity)
# NOTE  this is unlike the ions.incremental output as calculated in IsoX ions. incremental = (Intensity/Peak.Noise) * 3 *sqrt(240000/FT.Resolution) *sqrt(Microscans)
# each isotope ratio is divided by the number of relevant atoms in the ion/fragment of interest as provided by the user in a vector (ie. (13C/12C)/6 for a C6 fragment)

# generates list with several data.frames each using different data culling and background correction methods, which the user can choose from based on preference.
# Table: Mean, median, stddev, stderr, RSE  for each ion species of interest X unfiltered and data filtered for <1/2/3 sigma from mean tic/it/ticxit, <10/20/30% Median average deviation (MAD) X baseline/background subtracted vs. non subtracted

# Plot RSE vs. Shotnoise for each ion species of interest

#----------------
#load required packages
packages = c("plyr","dplyr","ggplot2","ggpubr","slider","psych","tidyverse","scales","cowplot","tidyquant","gtable","ggExtra","ggstance","gridExtra","grid","lubridate","datasets","gapminder","ggsignif","reshape2", "viridis")
# Now load or install&load all
package.check <- lapply(packages,
                        FUN = function(x) {
                          if (!require(x, character.only = TRUE)) {
                            install.packages(x, dependencies = TRUE)
                            library(x, character.only = TRUE)
                          }
                        }
)

#-------------
#define function to import data set extracted from raw input using isoX. Use combined.isox file for multiple acquisitions

importIdata <- function(df,ranges,masses,ions,ilogs){

  #reformatting
  data= df
  colnames(masses)=c("compound","isotopolog","m/z", "tolerance [mmu]", "z")
  data= merge(data,masses[,c(1:3)], by=c("compound","isotopolog")) #append data.frame with exact masses for fragments
  data= merge (data, ranges, by= "filename", all.x=T, all.y=T) #append data.frame with start and end times for acquisition

  #add additional data quality columns and filter to asses incomplete/duplicate rows
  data= data %>%
    mutate("massdefect[da]" = abs(data$mzMeasured-data$`m/z`),#append column mass defect = absolute amount deviation from exact mass in daltons
           Nio = (intensity/peakNoise)*(4.4 /1)*((120000/resolution)^0.5)*(microscans^0.5) #append column with calculated Nio values (Eiler 2017) CN=4.4 ref.resolution=120000, z=1
    ) %>%
    filter(!is.na(data$start.min)& #only include acquisitions that are also defined in meta data (with start and end times)
             data$compound %in% ions & #only include scans with fragments/ions of interest i.e. C6H13 fragment
             data$isotopolog %in% ilogs #only include scans with isotopologs of interest i.e. C6H13 fragment with 12C (unsubstituted) and 1data13C
    ) %>%
    group_by(filename, compound, isotopolog) %>%
    mutate(Nio_weighted = (Nio*intensity)/sum(intensity),.after = Nio # weighted Nio according to Eiler lab meeting
    )%>%
    ungroup() %>%
    group_by(filename, scan.no, time.min, compound, isotopolog) %>% # rank duplicate isotopolog (assigned b/c mass match) in same scan based on signal intensity
    mutate(num_dups = n(), rank= rank(-intensity),
           "ticxit"= tic * it.ms
    ) %>%
    ungroup() %>%
    mutate(is_duplicated = num_dups > 1)%>%  # logical column identifying if there are duplicate isotopologs in single scan
    group_by(filename, scan.no, time.min, compound) %>% # append column to check if each compound has all isotope species of interest
    mutate(num_ilogs = length(unique(isotopolog))) %>%
    ungroup() %>%
    mutate(has_all_ilogs = num_ilogs == (length(ilogs)))# logical column identifying if all ilogs defined in ilogs are present in single scan
  data = data[order(data$filename,data$compound,data$isotopolog,data$scan.no),]
  data = data[,c("filename","compound","isotopolog","scan.no","time.min","intensity","ions.incremental",
                 "tic","it.ms","ticxit","resolution","peakResolution","peakNoise","basePeakIntensity","rawOvFtT",
                 "intensCompFactor","agc","agcTarget","microscans","numberLockmassesFound","analyzerTemperature",
                 "start.min","end.min","mzMeasured","m/z","massdefect[da]","Nio","Nio_weighted",
                 "num_dups","rank","is_duplicated","num_ilogs","has_all_ilogs")]#reorder columns
}

#output
#importIdata = dataframe only containing columns for further analysis i.e. time ranges for each acquisition, filter for fragments/ions of interest and isotopologs of interest and mass defect value.


#define function that will filter ilogs based on exact mass +- mass tolerance mmu, remove duplicate rows assigned to same ilog or rows that dont have all ilogs and add background level column per filename, compound, isotopolog

filter.Idata <- function (import.Idata_results, tolerance){

  #filter and remove columnes used to filter dataset
  mdf=importIdata_results[import.Idata_results$"massdefect[da]" <tolerance,] #removes rows with mass deferct larger then defined tolerance value in daltons
  is=mdf[mdf$rank == "1",] #removes rows assigned in duplicate to same isotopolog (in same scan) by keeping only the one with higher signal intensity
  filtered=is[is$has_all_ilogs == "TRUE",] #removes rows that do not contain the desired number of isotopologs in ilogs object
  filtered=filtered[,c(1:28)]#removes columns used for filtering

  #add background value columns
  b=filtered %>%
    group_by (filename)%>%
    arrange(desc(time.min)) %>%
    filter(time.min <= start.min)
  b=b%>% # background calculation method 1
    group_by(filename,compound,isotopolog)%>%
    top_n(4,time.min)%>%
    mutate(Nio_avrg4_background = mean(Nio),
           Nio_w_avrg4_background = mean(Nio_weighted))
  b=b%>% # background calculation method 2
    group_by(filename,compound,isotopolog)%>%
    top_n(20,time.min)%>%
    mutate(Nio_min20_background = min(Nio),
           Nio_w_min20_background = min(Nio_weighted))
  data_merged=merge(filtered,b[,c("filename","compound","isotopolog","Nio_avrg4_background","Nio_min20_background","Nio_w_avrg4_background","Nio_w_min20_background")], by=c("filename","compound","isotopolog"))

  #filter to only include rows within defined peak elution time
  data_out = data_merged%>%
    filter(time.min <= end.min & time.min >= start.min)
  data_out = distinct(data_out)

  #add columns with background subtracted
  data_out = data_out %>%
    mutate(Nio_b.avrg4 = Nio - Nio_avrg4_background,
           Nio_b.min20 = Nio - Nio_min20_background,
           Nio_w_b.min20 = Nio - Nio_w_min20_background,
           Nio_w_b.avrg4 = Nio - Nio_w_avrg4_background)

  return(data_out)
}

#output
#filtered = dataframe only containing scans with mass defect as defined by tolereance value, unique assignment of one isotopolog per defined exact mass,
#including scan that have all defined isotopologs of interest per scan within peak elution time and background values for each filename,compound, ilog of interest (2 background filter methods)


#calculate isotope-ratio values ; ions.incremental is not = Eiler 2017 NIO calculation!!! restructure data frame
#NL ratio, Nio ratio , Nio weighted ratio
#ONLY WORKS FOR TWO ILOGS!

I.ratios <- function(filter.Idata_results){
  data_ratios = filter.Idata_results
  data_ratios = data_ratios %>% #split dataframe by filename/compound/isotopologs and combine isotpolgs by matching filename and scan.no to get isotopologs in same row aligned by same scan number and filename
    group_by(filename, compound, isotopolog) %>%
    split(.$isotopolog)%>%
    lapply(function(x) x[(names(x) %in% c("filename","compound","scan.no","time.min","isotopolog","tic","it.ms","ticxit",
                                          "massdefect[da]","intensity","Nio","Nio_weighted","Nio_b.avrg4","Nio_b.min20",
                                          "Nio_w_b.avrg4","Nio_w_b.min20"))])%>%
    lapply(function(i) i[order(i$filename,i$scan.no),])%>%
    as.data.frame()%>% #this is the reason the script only works with 2 isotopologs
    mutate(NL_ratio = X1.intensity/X0.intensity, #NL ratio
           Nio_ratio = X1.Nio/X0.Nio, #Nio ratio
           Nio_w_ratio = X1.Nio_weighted/X0.Nio_weighted, #Nio weighted ratio
           Nio_b.avrg4_ratio = X1.Nio_b.avrg4/X0.Nio_b.avrg4, #Nio_b.avrg4 ratio
           Nio_w_b.avrg4_ratio = X1.Nio_w_b.avrg4/X0.Nio_w_b.avrg4, #Nio_b.avrg4 weighted ratio
           Nio_b.min20_ratio = X1.Nio_b.min20/X0.Nio_b.min20, #Nio_b.avrg4 ratio
           Nio_w_b.min20_ratio = X1.Nio_w_b.min20/X0.Nio_w_b.min20) #Nio_b.avrg4 weighted ratio
  data_out=data_ratios
  #background subtraction creates negative values (interestingly no neg values for background subtracted values using Nio_weighted) replace with NA
  is.na(data_out$Nio_b.avrg4_ratio) <- data_out$Nio_b.avrg4_ratio < 0
  is.na(data_out$Nio_b.min20_ratio) <- data_out$Nio_b.min20_ratio < 0
  data_out %>%
    select(X0.filename,X0.compound,X0.scan.no,X0.time.min,X0.tic,X0.it.ms,X0.ticxit,X0.intensity,X1.intensity,
           X0.Nio,X1.Nio,X0.Nio_weighted,X1.Nio_weighted,
           X0.Nio_b.avrg4,X1.Nio_b.avrg4,X0.Nio_b.min20,X1.Nio_b.min20,
           X0.Nio_w_b.avrg4,X1.Nio_w_b.avrg4,X0.Nio_w_b.min20,X1.Nio_w_b.min20,
           NL_ratio,Nio_ratio,Nio_w_ratio,Nio_b.avrg4_ratio,Nio_w_b.avrg4_ratio,Nio_b.min20_ratio,Nio_w_b.min20_ratio)%>%
    rename(filename=X0.filename,compound=X0.compound,scan.no=X0.scan.no,time.min=X0.time.min,tic=X0.tic,it.ms=X0.it.ms,ticxit=X0.ticxit,
           unsub.intensity=X0.intensity,sub.intensity=X1.intensity,unsub.Nio=X0.Nio,sub.Nio=X1.Nio,
           unsub.Nio_weighted=X0.Nio_weighted,sub.Nio_weighted=X1.Nio_weighted,
           unsub.Nio_b.avrg4=X0.Nio_b.avrg4,sub.Nio_b.avrg4=X1.Nio_b.avrg4,
           unsub.Nio_b.min20=X0.Nio_b.min20,sub.Nio_b.min20=X1.Nio_b.min20,
           unsub.Nio_w_b.avrg4=X0.Nio_w_b.avrg4, sub.Nio_w_b.avrg4=X1.Nio_w_b.avrg4,
           unsub.Nio_w_b.min20=X0.Nio_w_b.min20,sub.Nio_w_b.min20=X1.Nio_w_b.min20)#reduce and rename columns
}

#output
#data frame contains substituted and unsubstituted iostopolog data in columns next to each other with isotope ratios (NOT DEVIDED by #Cs)
#negative values after background subtraction are replaced with NA


#define culling function with options for different culling methods USING NIO values not weighted
#used to evaluate differnt culling options
#culled using ticxit

cull.Idata <- function(filter.Idata_results,ratio,it){ #better to do per filename and/or ion of interest b/c output becomes too large background subtracted columns becaome inaccurate because negative values have to be NA
  data= filter.Idata_results
  options(dplyr.summarise.inform = FALSE)  # Suppress summarise info
  data_0 = data %>%#----
  mutate(filter="no filter")
  data_0 = data_0
  summary_0 = data_0  %>%
    group_by(filename, compound) %>%
    summarise(
      filter="no filter",
      mean = mean(ratio),
      median = median(ratio),
      n = n(),
      sd = sd(ratio),
      se = sd / sqrt(n),
      RSD= sd(ratio)/mean(ratio)*100,
      RSE=se/mean(ratio)*100,
    )
  data_it = data %>% #----
  mutate(filter="it") %>%
    group_by(filename, compound) %>%
    filter(it.ms >= it) # filter dataset to only include rows that are <set "it" value to ensure constant AGC target
  summary_it = data_it  %>%
    group_by(filename, compound) %>%
    summarise(
      filter="it",
      mean = mean(ratio),
      median = median(ratio),
      n = n(),
      sd = sd(ratio),
      se = sd / sqrt(n),
      RSD= sd(ratio)/mean(ratio)*100,
      RSE=se/mean(ratio)*100,
    )
  data_sd1 = data %>% #----
  mutate(filter="sd1") %>%
    group_by(filename, compound) %>%
    mutate(plus_o=mean(ticxit)+(sd(ticxit)*1), #append group wise column with x+ sigma*level and x- sigma*level
           minus_o=mean(ticxit)-(sd(ticxit)*1)) %>%
    filter(ticxit <= plus_o & ticxit >= minus_o) # filter dataset to only include rows that are +- level sigma deviation of Nio
  summary_sd1 = data_sd1  %>%
    group_by(filename, compound) %>%
    summarise(
      filter="sd1",
      mean = mean(ratio),
      median = median(ratio),
      n = n(),
      sd = sd(ratio),
      se = sd / sqrt(n),
      RSD= sd(ratio)/mean(ratio)*100,
      RSE=se/mean(ratio)*100,
    )
  data_sd2 = data %>%#----
  mutate(filter="sd2") %>%
    group_by(filename, compound) %>%
    mutate(plus_o=mean(ticxit)+(sd(ticxit)*2), #append group wise column with x+ sigma*level and x- sigma*level
           minus_o=mean(ticxit)-(sd(ticxit)*2)) %>%
    filter(ticxit <= plus_o & ticxit >= minus_o) # filter dataset to only include rows that are +- level sigma deviation of Nio
  summary_sd2 = data_sd2  %>%
    group_by(filename, compound) %>%
    summarise(
      filter="sd2",
      mean = mean(ratio),
      median = median(ratio),
      n = n(),
      sd = sd(ratio),
      se = sd / sqrt(n),
      RSD= sd(ratio)/mean(ratio)*100,
      RSE=se/mean(ratio)*100,
    )

  data_sd3 = data %>%#----
  mutate(filter="sd3") %>%
    group_by(filename, compound) %>%
    mutate(plus_o=mean(ticxit)+(sd(ticxit)*3), #append group wise column with x+ sigma*level and x- sigma*level
           minus_o=mean(ticxit)-(sd(ticxit)*3)) %>%
    filter(ticxit <= plus_o & ticxit >= minus_o) # filter dataset to only include rows that are +- level sigma deviation of Nio
  summary_sd3 = data_sd3  %>%
    group_by(filename, compound) %>%
    summarise(
      filter="sd3",
      mean = mean(ratio),
      median = median(ratio),
      n = n(),
      sd = sd(ratio),
      se = sd / sqrt(n),
      RSD= sd(ratio)/mean(ratio)*100,
      RSE=se/mean(ratio)*100,
    )

  data_MAD10 = data %>%#----
  mutate(filter="MAD10") %>%
    group_by(filename, compound) %>%
    mutate(average= mean(ticxit), # append group wise column with +- MAD level
           MAD_insidef= abs(ticxit-average),
           plus_MAD= median(MAD_insidef)+(0.1*median(MAD_insidef)),
           minus_MAD= median(MAD_insidef)-(0.1*median(MAD_insidef)))%>%
    filter(ticxit <= plus_MAD & ticxit >= minus_MAD) # filter dataset to only include rows that are +- level MAD deviation of Nio
  summary_MAD10 = data_MAD10  %>%
    group_by(filename, compound) %>%
    summarise(
      filter="MAD10",
      mean = mean(ratio),
      median = median(ratio),
      n = n(),
      sd = sd(ratio),
      se = sd / sqrt(n),
      RSD= sd(ratio)/mean(ratio)*100,
      RSE=se/mean(ratio)*100,
    )

  data_MAD20 = data %>%#----
  mutate(filter="MAD20") %>%
    group_by(filename, compound) %>%
    mutate(average= mean(ticxit), # append group wise column with +- MAD level
           MAD_insidef= abs(ticxit-average),
           plus_MAD= median(MAD_insidef)+(0.2*median(MAD_insidef)),
           minus_MAD= median(MAD_insidef)-(0.2*median(MAD_insidef)))%>%
    filter(ticxit <= plus_MAD & ticxit >= minus_MAD) # filter dataset to only include rows that are +- level MAD deviation of Nio
  summary_MAD20 = data_MAD20  %>%
    group_by(filename, compound) %>%
    summarise(
      filter="MAD20",
      mean = mean(ratio),
      median = median(ratio),
      n = n(),
      sd = sd(ratio),
      se = sd / sqrt(n),
      RSD= sd(ratio)/mean(ratio)*100,
      RSE=se/mean(ratio)*100,
    )
  #----
  summary_data=rbind(summary_0,summary_it,summary_sd1,summary_sd2,summary_sd3,summary_MAD10,summary_MAD20)
  summary_data= summary_data %>%
    group_by(filename, compound) %>%
    slice(which.min(RSE))
  outList= list(data_0,data_it,data_sd1,data_sd2,data_sd3,data_MAD10,data_MAD20,summary_data)
  names(outList) <- c("no_filter","it","1sd","2sd","3sd","MAD10","MAD20","summary_data")
  return(outList)
}

#output
#data.frame of isotope ratios divided by number of atoms of interest


#summary data for different NL ratio calculation options per compound
#drop NA rows for background subtracted ratios
#this function is very repetitive have to rewrite summarize with loops

summary.ratios <- function(avrg.ratios_results){
  data = avrg.ratios_results
  NL_ratio = data %>%
    group_by (filename, compound) %>%
    summarise(
      ratio = "NL_ratio",
      mean = mean(NL_ratio),
      median = median(NL_ratio),
      n = n(),
      sd = sd(NL_ratio),
      se = sd / sqrt(n),
      RSD= sd(NL_ratio)/mean(NL_ratio)*100,
      RSE=se/mean(NL_ratio)*100
    )
  Nio_ratio = data %>%
    group_by (filename, compound) %>%
    summarise(
      ratio = "Nio_ratio",
      mean = mean(Nio_ratio),
      median = median(Nio_ratio),
      n = n(),
      sd = sd(Nio_ratio),
      se = sd / sqrt(n),
      RSD= sd(Nio_ratio)/mean(Nio_ratio)*100,
      RSE=se/mean(Nio_ratio)*100
    )
  Nio_w_ratio = data %>%
    group_by (filename, compound) %>%
    summarise(
      ratio = "Nio_w_ratio",
      mean = mean(Nio_w_ratio),
      median = median(Nio_w_ratio),
      n = n(),
      sd = sd(Nio_w_ratio),
      se = sd / sqrt(n),
      RSD= sd(Nio_w_ratio)/mean(Nio_w_ratio)*100,
      RSE=se/mean(Nio_w_ratio)*100
    )
  Nio_b.avrg4_ratio = data %>%
    group_by (filename, compound) %>%
    summarise(
      ratio = "Nio_b.avrg4_ratio",
      mean = mean(Nio_b.avrg4_ratio,na.rm = TRUE),
      median = median(Nio_b.avrg4_ratio,na.rm = TRUE),
      n = n(),
      sd = sd(Nio_b.avrg4_ratio,na.rm = TRUE),
      se = sd / sqrt(n),
      RSD= sd(Nio_b.avrg4_ratio,na.rm = TRUE)/mean(Nio_b.avrg4_ratio,na.rm = TRUE)*100,
      RSE=se/mean(Nio_b.avrg4_ratio,na.rm = TRUE)*100
    )
  Nio_w_b.avrg4_ratio = data %>%
    group_by (filename, compound) %>%
    summarise(
      ratio = "Nio_w_b.avrg4_ratio",
      mean = mean(Nio_w_b.avrg4_ratio),
      median = median(Nio_w_b.avrg4_ratio),
      n = n(),
      sd = sd(Nio_w_b.avrg4_ratio),
      se = sd / sqrt(n),
      RSD= sd(Nio_w_b.avrg4_ratio)/mean(Nio_w_b.avrg4_ratio)*100,
      RSE=se/mean(Nio_w_b.avrg4_ratio)*100
    )
  Nio_b.min20_ratio = data %>%
    group_by (filename, compound) %>%
    summarise(
      ratio = "Nio_b.min20_ratio",
      mean = mean(Nio_b.min20_ratio,na.rm = TRUE),
      median = median(Nio_b.min20_ratio,na.rm = TRUE),
      n = n(),
      sd = sd(Nio_b.min20_ratio,na.rm = TRUE),
      se = sd / sqrt(n),
      RSD= sd(Nio_b.min20_ratio,na.rm = TRUE)/mean(Nio_b.min20_ratio,na.rm = TRUE)*100,
      RSE=se/mean(Nio_b.min20_ratio,na.rm = TRUE)*100
    )
  Nio_w_b.min20_ratio = data %>%
    group_by (filename, compound) %>%
    summarise(
      ratio = "Nio_w_b.min20_ratio",
      mean = mean(Nio_w_b.min20_ratio),
      median = median(Nio_w_b.min20_ratio),
      n = n(),
      sd = sd(Nio_w_b.min20_ratio),
      se = sd / sqrt(n),
      RSD= sd(Nio_w_b.min20_ratio)/mean(Nio_w_b.min20_ratio)*100,
      RSE=se/mean(Nio_w_b.min20_ratio)*100
    )
  data_out = rbind(NL_ratio,Nio_ratio,Nio_w_ratio,Nio_b.avrg4_ratio,Nio_w_b.avrg4_ratio,Nio_b.min20_ratio,Nio_w_b.min20_ratio)
  return(data_out)
}

#output
#summary data table for each ratio calculation type per compound and filename


##############
#Shot-noise limit calcualtion per scan. modified after Cajetan Neubauers script
#have to subset data.frame to only contain one file of interest, ion of interest
#have to subset data.frame to conly contain unsub/sub. count to use for SNL and which ratio calculation method is preferred
#doesnt work for Nio_b.avrg4_ratio and Nio_b.min20_ratio b/c rows contain NA. Have to remove them for SNL calc without affecting the other ratios

calc.shotnoise <- function(d.neff){
  df.out <- data.frame(scans =1:length(d.neff$scan.no),
                       neff = 1:length(d.neff$scan.no),
                       rse= 1:length(d.neff$scan.no),
                       shotnoise = 1:length(d.neff$scan.no),
                       shotnoise.neff = 1:length(d.neff$scan.no))
  #Geometric mean, geometric SD and SE
  gmean <- function(x) exp(mean(log(x))) #define geometric mean
  gsd <- function(x) exp(mean(log(x)) + sd(log(x))) - exp(mean(log(x)))
  gse <- function(x) (exp(mean(log(x)) + sd(log(x))) - exp(mean(log(x)))) /sqrt(length(x))

  for (l in 1:length(d.neff$scan.no)) {
    df <- d.neff[1:l,]
    neff <- (sum(df$unsub.Nio) * sum(df$sub.Nio))/ (sum(df$unsub.Nio) + sum(df$sub.Nio))
    df.out$neff[l] <- neff
    shotnoise.neff <- neff^-0.5
    df.out$shotnoise.neff[l] <- shotnoise.neff
    shotnoise <- ((1/sum(df$sub.Nio) + (1/sum(df$unsub.Nio))))^0.5
    df.out$shotnoise[l] <- shotnoise
    rse <- gse(d.neff[1:l,]$Nio_ratio) / gmean(d.neff[1:l,]$Nio_ratio)
    df.out$rse[l] <- rse
  }
  df.out <- df.out %>% na.omit() %>% select(-shotnoise, -shotnoise.neff)
  df.out <- df.out %>% gather(key, value, -scans, -neff)
  return(df.out)
}

#output
#table with shot noise limits (value) and number of effective ion counts (n.eff) per scan. per ratio calculation method

#-----------------
