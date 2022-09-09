####### Generalized script to asses isotope ratios of ions of interest from Orbi measurements #############

#------------------
# created by Merve Öztoprak
# m.oeztoprak@gmail.com
#------------------

#' function to calculate average isotope ratios of fragments/ions of interest. I.e. isotope ratio divided by the number of relevant atoms in fragment/ion of interest
#'
#' @param x data frame that results from applying the iratios() function
#' @param ions vector which contains all fragments/ions of interest
#' @param atoms vector wchih contains number of relevant atoms in fragment/ions of interest (ORDER MATTER!, has to be the same as in ions
#' @return data.frame of isotope ratios divided by number of atoms of interest
#' @examples
#'
#------------
#' @export
avrg.ratios <- function(x,ions,n.atoms){
  data=x
  div.factor=data.frame(ions,n.atoms)
  div.factor=div.factor%>%
    rename(compound=ions)
  data_out = data
  data_out = data_out %>% #divide ratios by number of ataoms of interest
    left_join(div.factor, by="compound")%>%
    mutate(across(c("unsub.intensity","sub.intensity",
                    "unsub.Nio","sub.Nio",
                    "unsub.Nio_weighted","sub.Nio_weighted",
                    "unsub.Nio_b.avrg4","sub.Nio_b.avrg4",
                    "unsub.Nio_b.min20","unsub.Nio_b.min20",
                    "unsub.Nio_w_b.avrg4","sub.Nio_w_b.avrg4",
                    "unsub.Nio_w_b.min20","unsub.Nio_w_b.min20",
                    "NL_ratio","Nio_ratio","Nio_w_ratio",
                    "Nio_b.avrg4_ratio","Nio_w_b.avrg4_ratio",
                    "Nio_b.min20_ratio","Nio_w_b.min20_ratio"
                    ),.fns=~./n.atoms))%>%
    relocate(n.atoms,.after=compound)
}



