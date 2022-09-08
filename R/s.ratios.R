####### Generalized script to asses isotope ratios of ions of interest from Orbi measurements #############

#------------------
# created by Merve Öztoprak
# m.oeztoprak@gmail.com
#------------------

#'summary data for different NL ratio calculation options per compound
#'drops NA rows for background subtracted ratios
#'this function is very repetitive have to rewrite summarize with loops
#' @param avrg.ratios_results results of avrg.ratios function
#' @return  = summary data table for each ratio calculation type per compound and filename
#' @examples

#----------------
#' @export
s.ratios <- function(avrg.ratios_results){
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

