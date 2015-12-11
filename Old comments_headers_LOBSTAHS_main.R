# LOBSTAHS_main.R
#
# Created 12/4/2015 by Jamie Collins, james.r.collins@aya.yale.edu
#
# Purpose: Applies Lipid and Oxylipin Biomarker Screening through Adduct Hierarchy Sequences (LOBTAHS) screening to ID and annotate features in a multiple sample, HPLC-ESI-MS lipid dataset that has been processed using xcms and CAMERA. This dataset must contain peaks of a single ion mode/polarity (i.e., positive or negative). The script prepOrbidata.R can be used for the necessary processing in xcms and CAMERA. Tasks accomplished by prepOrbidata.R include peak picking, peak grouping, retention time correction, isotope identification, creation of pseudospectra using correlation of xcms peak groups between and within samples, and creation of the necessary xsAnnotate object.
#
# LOBSTAHS is under development by the Van Mooy Lab at Woods Hole Oceanographic Institution. Currently, this script is written to analyze lipid data from the experiment described in Graff van Creveld et al., 2015, "Early perturbation in mitochondria redox homeostasis in response to environmental stress predicts cell fate in diatoms," ISME Journal 9:385-395. This dataset is used to demonstrate the LOBSTAHS lipidomics pipeline in Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy, 2015, "Untargeted discovery and identification of oxidative stress biomarkers using a lipidomics pipeline for complex datasets."
#
# This script:
#
#  1. Conducts an in silico simulation to create the necessary lipid-oxylipin database from which assignments will be made. This is accomplished by calling the separate script LOBSTAHS_DB_genr.R using source(). LOBSTAHS_DB_genr.R will create a database according to chemical parameters specified by the user in several simple, editable input tables; the defaults will create entries for a wide range of intact polar lipids, oxidized intact lipids (ox-lipids), and oxylipins. Each entry in the database represents a different adduct ion of a compound of interest; fields include exact mass, empirical formula, and empirically determined adduct ion hierarchy data for the parent class of lipid. The database generated will be for the ion mode appropriate to the data. Users are cautioned that certain chemical species will only be included in the database for a single ion mode (for example, simulations for free fatty acids are only conducted in negative ion mode). LOBSTAHS_DB_genr.R can be used on its own to generate databases in .csv format, if this is desired. LOBSTAHS_DB_genr.R contains more specific, additional documentation about the various chemical properties that may be considered. If databases of the correct format have already been generated, the user is given the option to import directly from a .csv file.
#
#  2. The script then iterates through each pseudospectrum created by CAMERA (see script prepOrbidata.R) and performs the following tasks:
#
#    a. Eliminates of any secondary isotope peaks identifed by CAMERA. Since LOBSTAHS is designed to facilitate relative comparisons between samples in the same dataset, only monoisotopic peaks are retained.
#
#    b. Makes putative compound assignments to the dataset using exact mass. Assignments are made from the lipid-oxylipin database of appropriate ion mode created in LOBSTAHS_DB_genr.R.
#
#    c. Validation of these compound assignments is then accomplished by imposition of several orthogonal filters:
#
#      i. First, the script uses empirical observations of abundances of HPLC-ESI-MS adduct ions (given in LOBSTAHS_adduct_ion_hierarchies_pos.csv and LOBSTAHS_adduct_ion_hierarchies_neg.csv) to determine whether the features within each CAMERA pseudospectrum group satisfy the hierarchy for the parent lipid class of the putatively assigned compound. Only assignments for which the pseudospectrum components partially or fully satisfy the criteria are retained. Retained assignments are annotated for their degree of conformity with the relevant adduct ion hierarchy.
#
#      ii. If elected by the user, retention time window criteria are applied to exclude any assignments whose feature retention time lies outside of the range expected for the lipid class of the assigned compound. (These criteria will be particular to the chromatography used; must be specified in the .additional csv file, "VML_rt_windows.csv")
#
#      iii. If elected by user, the script will also eliminate any assignments containing an odd total number of acyl (fatty acid chain) carbon atoms. This criterion is useful if the user is evaluating biological data of eukaryotic origin, since we do not expect eukaryotes to synthesize fatty acids with odd numbers of carbon atoms.
#
#  3. In the course of making and evaluating these assignments, the script will annotate features to permit discovery of various types of isomer, indicate compliance with the relevant adduct ion hierarchy, and flag assignments which cannot be disambiguated from others within the ppm mass tolerance used in matching. Descriptions of all codes used in annotation are given in Collins et al. 2015.
#
#  4. Results are then exported to a .csv file to facilitate any follow-on analysis.
#
# As described in Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy, 2015, "Untargeted discovery and identification of oxidative stress biomarkers using a lipidomics pipeline for complex datasets"
#
# Obtain current versions of scripts and necessary dependencies at https://github.com/vanmooylipidomics/LOBSTAHS/
#
# Please direct questions/suggestions/comments to Jamie Collins, james.r.collins@aya.yale.edu, or Helen Fredricks, hfredricks@whoi.edu
#
# Revision history maintained on GitHub
#
################ Caveats and prerequisites #############
#
# Presumes user has installed the R packages "xcms", "CAMERA", "tools" with all required dependencies
#
# This script the following inputs:
#
#  1. A multiple sample, high-mass-accuracy HPLC-ESI-MS lipid dataset consisting of features in a single ion mode. This dataset should be in the form of an xsAnnotate object to which xcms and CAMERA pre-processing has already been applied. The script prepOrbidata.R can be used to perform the necessary processing. prepOrbidata.R is written to produce an object called xset_a, which is required by the code below. 
#
#  2. The .additional csv file, "VML_rt_windows.csv", containing the expected retention time windows for each class of lipid, if retention time screening is desired. This file is available at https://github.com/vanmooylipidomics/LOBSTAHS. The user should alter the retention times in this matrix according to his or her particular chromatography; the Van Mooy Lab data is distributed as a template.
#
# In addition, the user should verify that LOBSTAHS_DB_genr.R and its required .csv tables are in their default locations. This enables the database to be properly generated.

################ Initial setup and variable definition #############

# run this section before specifying user settings, below

# clear workspace
# 
# WS <- c(ls())
# rm(WS, list = WS)

# load required packages

library(tools) 

library(xcms)

library(CAMERA)

library(rsm)

# library(Rmpi)

library(snowfall) # if multicore tasking is desired

# ******************************************************************
################ Basic user begin editing here #############
# ******************************************************************

################ User: define locations of data files and database(s) #############

# *if* the in silico simulation results tables sim_results.neg and sim_results.pos are still in your R environment (because, say, you just ran "LOBSTAHS_DB_genr.R" and didn't clear the workspace), you don't have to specify the below; otherwise, you have to tell R where the database data is
DBfile_pos = "dependencies/LOBSTAHS_lip-oxy_DB_pos_2015-11-21T22-42-44-0300.csv" 
DBfile_neg = "dependencies/LOBSTAHS_lip-oxy_DB_neg_2015-11-21T22-42-45-0300.csv"

# ******************************************************************
################ Basic user stop editing here #############
# ******************************************************************