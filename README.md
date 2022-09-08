# IsoOrbiStats

Welcome to the GitHub page of IsoOrbiStats. Here I aim to develop a generalized script to evaluate intramolecular isotopic measurements on the Orbitrap platform for Position Specific Isotopic Analysis (PSIA).

Currently, the Master IsoOrbiStats Script contains all functions which are meant to be used consecutively and can be adapted to desired specifications.
The additonal files are required for use as a package within R. The script will be further expanded before submission to CRAN where IsoOrbiStats will be available for download directly within R. For now functions can be copied for use from the Master file or the folder may be downloaded directly into the working directory and installed using install("IsoOrbiStats") in R. 

The script has been tested on two datasets and is still in development.

The input for all functions are:

1. Dataframe of the combined.isox file that contains all files to be analyzed and is generated via isox.

2. Dataframe of isotopologs.tsv, which is also used to generate the .isox file

3. Dataframe that has 3 columns
  - filename (as in the isoX file)
  - start.min (start of the eluting peak)
  - end.min (end of the eluting peak)
  for each file to be analyzed

4. A factor named "tolerance" which is the amount of deviation allowed from the exact mass of a molecule given in Daltons

5. A vector containing ions of interest to be further analyzed (contained in isotopologs.tsv Compound column)

6. A vector containing isotopologs of interest to be further analyzed (contained in isotopologs.tsv Isotopolog column)

7. A vector supplying the number of relevant atoms in the ions of interest (ORDER MATTERS! Number of relevant atoms in ion has to fallow the same order as in ions vector)

There are two test folders that contain an applied script using the functions in IsoOrbiStats :

- An example of Phytane measured on a GC-Orbitrap platform where Phytane is captured in a reservoir and allowed to diffuse to the Orbitrap over the course of 70 minutes.

- An example of Phytol measured on a LC-Orbitrap platform via APCI soft ionization. Phytol is measured in both MS1 or in MS2 mode, where Phytol is fragmented in the HCD cell before introduction to the mass analyzer.


#####THINGS STILL IN DEVELOPMENT###### 

automatic peak picking --> replaces start.min and end.min 

function that outputs analytical plots: 
  - stats on background value choice 
  - stats on culling options choice 
  - isotopolog peak separation and stats on count of each isotopolog 
  - RSE plot for each fragment
