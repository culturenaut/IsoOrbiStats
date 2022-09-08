####### Generalized script to asses isotope ratios of ions of interest from Orbi measurements #############

#------------------
# created by Merve Öztoprak
# m.oeztoprak@gmail.com
#------------------

#'function that calautes isotope ratios from intensity values for each scan, also restructures the dataframe
#' ions.incremental is not = Eiler 2017 NIO calculation!!!
#'NL ratio, Nio ratio , Nio weighted ratio
#'ONLY WORKS FOR TWO ILOGS!
#' @param filter.Idata_results results from filter.Idata function
#' @return data frame contains substituted and unsubstituted iostopolog data in columns next to each other with isotope ratios (NOT DEVIDED by #Cs).Negative values after background subtraction are replaced with NA
#' @examples

#-------------
#' @export
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
