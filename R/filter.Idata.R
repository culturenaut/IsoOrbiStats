####### Generalized script to asses isotope ratios of ions of interest from Orbi measurements #############

#------------------
# created by Merve Öztoprak
# m.oeztoprak@gmail.com
#------------------

#'function that that will filter ilogs based on exact mass +- mass tolerance mmu, remove duplicate rows assigned to same ilog or rows that dont have all ilogs and add background level column per filename, compound, isotopolog
#' @param import.Idata_results  result from import.Idata function
#' @param tolerance = user defined amount of deviation from exact mass (mass defect) permitted in daltons
#' @return = dataframe only containing scans with mass defect as defined by tolereance value, unique assignment of one isotopolog per defined exact mass,
#'including only scan that have all defined isotopologs of interest per scan within peak elution time and background values for each filename,compound, ilog of interest (2 background filter methods)
#' @examples

#----------------
#' @export
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

