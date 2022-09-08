####### Generalized script to asses isotope ratios of ions of interest from Orbi measurements #############

#------------------
# created by Merve Öztoprak
# m.oeztoprak@gmail.com
#------------------

#'function that allows to import data set extracted from raw input using isoX. Use combined.isox file for multiple acquisitions
#'
#' @param df =data frame of isox data table output
#' @param ranges = data frame of .txt file with start and end times for the peak of interest in each measurement run
#' @param masses = dataframe of isotopologs.tsv file which was also used to generate the .isox file
#' @param ions = vector which contains all fragments/ions of interest
#' @param ilogs = vector of isotopologs of interest (has to be contained in isotopologs.tsv file)
#' @return = dataframe only containing columns for further analysis i.e. time ranges for each acquisition, filter for fragments/ions of interest and isotopologs of interest and mass defect value.
#' @examples

#-------------
#' @export
import.Idata <- function(df,ranges,masses,ions,ilogs){

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

