####### Generalized script to asses isotope ratios of ions of interest from Orbi measurements #############

#------------------
# created by Merve Öztoprak
# m.oeztoprak@gmail.com
#------------------

#' function to calculate average isotope ratios of fragments/ions of interest. I.e. isotope ratio divided by the number of relevant atoms in fragment/ion of interest
#'
#' @param Iratios_results = data frame that results from applying the iratios() function
#' @param ions = vector which contains all fragments/ions of interest
#' @param atoms = vector wchih contains number of relevant atoms in fragment/ions of interest (ORDER MATTER!, has to be the same as in ions
#' @return data.frame of isotope ratios divided by number of atoms of interest
#' @examples
#'
#------------
#' @export
avrg.ratios <- function(Iratios_results,ions,n.atoms){
  data=Iratios_results
  div.factor=data.frame(ions,n.atoms)
  div.factor=div.factor%>%
    rename(compound=ions)
  data_out = data
  data_out = data_out %>% #divide ratios by number of ataoms of interest
    left_join(div.factor, by="compound")%>%
    mutate(across(c(22:28),.fns=~./n.atoms))%>%
    relocate(n.atoms,.after=compound)
}



