This readme will cover all of the folders in this repository.  Each folder should have its own readme for further direction, meaning that they are all self contained. Folders will contain intermediary files to document processing steps.  Most paths for reading in files in the R scripts will have to be changed, as your directory structure will be different than mine.

RAW sequence data is available on figshare:
	https://figshare.com/articles/GenomeSize_lowcoverage_Maizedata/5117827
Names of fastq files correspond to the names of individuals in the data sets.
Cell_num
	Jinliang's analysis of shoot apical meristem growth in maize inbreds

FISH_correlations
	Data for raw counts of knobs in mexicana populations by fish, and population averages by sequenceing
	
supp_tevar
	Data for top 15 RNA/DNA TE families in maize landraces, and the R code to generate the supplemental heatmap, showing variation across lines.
	
SelectionTests
	All resources necessary to examine tests for selection in both maize and mexicana.  Kinship matrices are included, along with R code (including customized EMMA script that prints out the beta value).  Also included is the R code to generate figure 1 of the manuscript. Here you will find the reports from BWA mapping of repeats - raw read counts, percentages, and megabase corrected values. GS_repeated measures are the raw data from flow cytometry of maize inbred lines where we repeated measures twice.
	
MakingUTE
	Folder for making the UTE database.  BLAST is necessary, but all the perl scripts to parse and create the new mapping reference are included.
	
BWA_mapping
	Contains the general scripts and parameters for mapping, though raw data are uploaded to figshare.  Requires bwa installed, but contains the scripts to add up the read counts from the sam outputs.  Also contains the repeat references!
	
Mexicana_GCexperiment
	Contains all the data and scripts that were used to develop the model for the growth chamber experiment.  To redo the final analysis by yourself, go to the Reproducible_growthchamber subfolder. The final analysis can be run using this R script Final_mexgc_analysis_cleaned.R, and all the data is contained in the csv grote_Indexmoms_groundedLL_d16NA. These files are contained in the Reproducible_growthchamber subfolder.  For further exploration, see the other folders and R scripts.  This directory lacks central readme's, as the R scripts are better commented, and many intermediary files are retained for the review process.    
