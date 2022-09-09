####### Generalized script to asses isotope {{ratio}}s of ions of interest from Orbi measurements #############

#------------------
# created by Merve Öztoprak
# m.oeztoprak@gmail.com
#------------------

#'Function that filters the dataframe using different culling methods
#'currently culling on preferred ratio can be switched to Nio_ratio, Nio_weighted etc. etc.
#'1/2/3 standard deviations from mean cullvalue i.e ticxit
#'10%/20% MAD median(abs((cullvalue)i-mean(cullvalue))
#' @param x results from I.ratios function
#' @param it <= it value defined by user
#' @param cullvalue column to use for culling (i.e. any ratio of preference, ticxit, etc.etc.)
#' @param ratio  isotope ratio calculation method preferred by user i.e. either "NL_{{ratio}}", "Nio_{{ratio}}", "Nio_w_{{ratio}}". Doesnt work for background subtracted {{ratio}}s since negative values have to be replaced by NA first
#' @return list containing dataframes that have been culled to include a subset of the input dataframe 1/2/3 standard deviations, 10%/20% MAD from the mean of Nio and below a certain it value. Output summary table can be used to evaluate different culling options for different fragments best hit for min. RSE
#' @examples

#----------------
#' @export
cull.Idata <- function(x,cullvalue,ratio,it){ #better to do per filename and/or ion of interest b/c output becomes too large background subtracted columns becaome inaccurate because negative values have to be NA
  data= x
  options(dplyr.summarise.inform = FALSE)  # Suppress summarise info
  data_0 = data %>%#----
  mutate(filter="no filter")
  data_0 = data_0
  summary_0 = data_0  %>%
    group_by(filename, compound) %>%
    summarise(
      filter="no filter",
      mean = mean({{ratio}}),
      median = median({{ratio}}),
      n = n(),
      sd = sd({{ratio}}),
      se = sd / sqrt(n),
      RSD= sd({{ratio}})/mean({{ratio}})*100,
      RSE=se/mean({{ratio}})*100,
    )
  data_it = data %>% #----
  mutate(filter="it") %>%
    group_by(filename, compound) %>%
    filter(it.ms >= it) # filter dataset to only include rows that are <set "it" value to ensure constant AGC target
  summary_it = data_it  %>%
    group_by(filename, compound) %>%
    summarise(
      filter="it",
      mean = mean({{ratio}}),
      median = median({{ratio}}),
      n = n(),
      sd = sd({{ratio}}),
      se = sd / sqrt(n),
      RSD= sd({{ratio}})/mean({{ratio}})*100,
      RSE=se/mean({{ratio}})*100,
    )
  data_sd1 = data %>% #----
  mutate(filter="sd1") %>%
    group_by(filename, compound) %>%
    mutate(plus_o=mean({{cullvalue}})+(sd({{cullvalue}})*1), #append group wise column with x+ sigma*level and x- sigma*level
           minus_o=mean({{cullvalue}})-(sd({{cullvalue}})*1)) %>%
    filter({{cullvalue}} <= plus_o & {{cullvalue}} >= minus_o) # filter dataset to only include rows that are +- level sigma deviation of Nio
  summary_sd1 = data_sd1  %>%
    group_by(filename, compound) %>%
    summarise(
      filter="sd1",
      mean = mean({{ratio}}),
      median = median({{ratio}}),
      n = n(),
      sd = sd({{ratio}}),
      se = sd / sqrt(n),
      RSD= sd({{ratio}})/mean({{ratio}})*100,
      RSE=se/mean({{ratio}})*100,
    )
  data_sd2 = data %>%#----
  mutate(filter="sd2") %>%
    group_by(filename, compound) %>%
    mutate(plus_o=mean({{cullvalue}})+(sd({{cullvalue}})*2), #append group wise column with x+ sigma*level and x- sigma*level
           minus_o=mean({{cullvalue}})-(sd({{cullvalue}})*2)) %>%
    filter({{cullvalue}} <= plus_o & {{cullvalue}} >= minus_o) # filter dataset to only include rows that are +- level sigma deviation of Nio
  summary_sd2 = data_sd2  %>%
    group_by(filename, compound) %>%
    summarise(
      filter="sd2",
      mean = mean({{ratio}}),
      median = median({{ratio}}),
      n = n(),
      sd = sd({{ratio}}),
      se = sd / sqrt(n),
      RSD= sd({{ratio}})/mean({{ratio}})*100,
      RSE=se/mean({{ratio}})*100,
    )

  data_sd3 = data %>%#----
  mutate(filter="sd3") %>%
    group_by(filename, compound) %>%
    mutate(plus_o=mean({{cullvalue}})+(sd({{cullvalue}})*3), #append group wise column with x+ sigma*level and x- sigma*level
           minus_o=mean({{cullvalue}})-(sd({{cullvalue}})*3)) %>%
    filter({{cullvalue}} <= plus_o & {{cullvalue}} >= minus_o) # filter dataset to only include rows that are +- level sigma deviation of Nio
  summary_sd3 = data_sd3  %>%
    group_by(filename, compound) %>%
    summarise(
      filter="sd3",
      mean = mean({{ratio}}),
      median = median({{ratio}}),
      n = n(),
      sd = sd({{ratio}}),
      se = sd / sqrt(n),
      RSD= sd({{ratio}})/mean({{ratio}})*100,
      RSE=se/mean({{ratio}})*100,
    )

  data_MAD10 = data %>%#----
  mutate(filter="MAD10") %>%
    group_by(filename, compound) %>%
    mutate(m= median({{cullvalue}}), # append group wise column with +- MAD level
           MAD_insidef= abs({{cullvalue}}-m),
           plus_MAD= median({{cullvalue}})+(0.1*median(MAD_insidef)),
           minus_MAD= median({{cullvalue}})-(0.1*median(MAD_insidef)))%>%
    filter({{cullvalue}} <= plus_MAD & {{cullvalue}} >= minus_MAD) # filter dataset to only include rows that are +- level MAD deviation of Nio
  summary_MAD10 = data_MAD10  %>%
    group_by(filename, compound) %>%
    summarise(
      filter="MAD10",
      mean = mean({{ratio}}),
      median = median({{ratio}}),
      n = n(),
      sd = sd({{ratio}}),
      se = sd / sqrt(n),
      RSD= sd({{ratio}})/mean({{ratio}})*100,
      RSE=se/mean({{ratio}})*100,
    )

  data_MAD20 = data %>%#----
  mutate(filter="MAD20") %>%
    group_by(filename, compound) %>%
    mutate(m= median({{cullvalue}}), # append group wise column with +- MAD level
           MAD_insidef= abs({{cullvalue}}-m),
           plus_MAD= median({{cullvalue}})+(0.2*median(MAD_insidef)),
           minus_MAD= median({{cullvalue}})-(0.2*median(MAD_insidef)))%>%
    filter({{cullvalue}} <= plus_MAD & {{cullvalue}} >= minus_MAD) # filter dataset to only include rows that are +- level MAD deviation of Nio
  summary_MAD20 = data_MAD20  %>%
    group_by(filename, compound) %>%
    summarise(
      filter="MAD20",
      mean = mean({{ratio}}),
      median = median({{ratio}}),
      n = n(),
      sd = sd({{ratio}}),
      se = sd / sqrt(n),
      RSD= sd({{ratio}})/mean({{ratio}})*100,
      RSE=se/mean({{ratio}})*100,
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
