####### Generalized script to asses isotope ratios of ions of interest from Orbi measurements #############

#------------------
# created by Merve Öztoprak
# m.oeztoprak@gmail.com
#------------------

# script shows an example of how the IsoOrbiStats functions may be applied

# This example folder contains files which have been obtained form Orbitrap measurements conducted in John Eilers lab at Caltech in 2019
# Briefly Phytane was introduced to the GC via direct elution mode (as described in Eiler et al., 2017) to determine the elution time for our phytane peak of interest
# Subsequently, Phytane was injected again, but this time in "peak trapping" mode. Phytane was captured in a reservoir inside the GC oven as it eluded of the column.
# Phytane was then allowed to exit the reservior end elude into the Orbitrap MS where a 10amu Mass window around mass 85.5 amu was measured for a duration of 70min in order to capture all Phytane molecules eluding from the reservoir.

# raw files are provided in the RAW folder. The raw files have been extracted and combined via the application IsoX (contact Cajetan Neubauer for access).
# the .isox file as well as the .tsv file used by IsoX to extract the RAW files are used for this script and provided inside the Example 1 folder

#------------

# Input
df= read.table("combined.isox", header=TRUE)
#supply a table (read.csv...header=TRUE) with column names: filename, start.min, end.min
ranges= read.csv("peak_times.csv", header=TRUE)
# mass tolerance for mass assignment in daltons
tolerance= 0.001
#exact masses. use same file as IsoX
masses= read.table("isotopologs.tsv", header=F)
#supply list of ions of interest
ions= c("C6H13") # works with multiple fragments yay!!!!
#supply list of isotopologs of interest (ONLY WORKS WITH 2!!)
ilogs=c("0","1")
#supply number of relevant atoms in ions of interest (order matters!)
n.atoms=c(6)

#-------
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
#----------

#import
I_df = import.Idata(df,ranges,masses,ions,ilogs)
#filter for mass tolerance, fragments/ions of interest within defined peak time and only including scans with both all isotopologs
#and adds background value and background subtracted column (for Nio and weighted Nio) calculated 1) average of last4 scans before peak start and 2) minimum of last 20 scans before peak start
f_df = filter.Idata(I_df, tolerance)
#calculate isotope-ratio values and restructure data frame; ions.incremental is not = Eiler 2017 NIO calculation!!!
r_df = I.ratios(f_df)
#define culling function with options for different culling methods USING NIO values not weighted
#used to evaluate differnt culling options
#culled using ticxit
c_df = cull.Idata(r_df,Nio_ratio,Nio_ratio,300) # here I use Nio_ratio (only relevant when looking at summary data.frame), Nio_w_ratio or NL_ratio also works

#check summary file
check_up = c_df[["summary_data"]]
#grab whatever culled dataframe you like best based on evaluations using chek_up, a_df and s_df data.frames
final_df = c_df[["3sd"]] # here i exclude all datapoints larger then 3x std dev of mean Nio_ratio

#calculate average isotope ratios of fragments/ions of interest. I.e. isotope ratio divided by the number of relevant atoms in fragment/ion of interest
final_df = avrg.ratios(final_df,ions,n.atoms)
#summary data for different NL ratio calculation options per compound also drops NA rows for background subtracted ratios
final_df_summary =s.ratios_t(final_df)

plot(final_df$time.min[final_df$filename == "Cyano_Phytane_1"],final_df$Nio_ratio[final_df$filename == "Cyano_Phytane_1"])
plot(final_df$time.min[final_df$filename == "Halo_Phytane_1"],final_df$Nio_ratio[final_df$filename == "Halo_Phytane_1"])
plot(final_df$time.min[final_df$filename == "Std_Phytane_1"],final_df$Nio_ratio[final_df$filename == "Std_Phytane_1"])

