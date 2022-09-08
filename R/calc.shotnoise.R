####### Part of a generalized script to asses isotope ratios of ions of interest from Orbi measurements #############

#------------------
# created by Merve Öztoprak
# m.oeztoprak@gmail.com
#------------------

#'Shot-noise limit calcualtion per scan. modified after Cajetan Neubauers script
#'
#' @param d.neff = dataframe which is subset to to only contain one file of interest, ion of interest and one type of unsub/sub. count  and ratio calculation method to use for the SNL calculation. Doesnt work for Nio_b.avrg4_ratio and Nio_b.min20_ratio b/c rows contain NA. Have to remove them for SNL calc without affecting the other ratios
#' @return table with shot noise limits (value) and number of effective ion counts (n.eff) per scan.
#' @examples

#----------------
#' @export
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
