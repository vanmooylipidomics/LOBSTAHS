



# LOBSTAHS_DB_genr.R
#
# Created 11/19/2015 by Jamie Collins, james.r.collins@aya.yale.edu
#
# Purpose: Runs in silico simulations to generate exact masses for HPLC-ESI-MS adduct ions of a wide range of intact polar lipids, ox-lipids, and oxylipins, using input from user and data in a series of .csv tables. A component of the LOBSTAHS (Lipid and Oxylipin Biomarker Screening Through Adduct Hierarchy Sequences) lipidomics pipeline, developed in the Van Mooy Lab at Woods Hole Oceanographic Institution.
#
# As described in Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy, 2015, "Untargeted discovery and identification of oxidative stress biomarkers using a lipidomics pipeline for complex datasets"
#
# Obtain most recent version of script and accompanying databases at https://github.com/vanmooylipidomics/LOBSTAHS
#
# Please direct questions/suggestions/comments to Jamie Collins, james.r.collins@aya.yale.edu, or Helen Fredricks, hfredricks@whoi.edu
#
# Revision history maintained on GitHub
#
################ Caveats and prerequisites #############
#
# This script requires four tables as inputs; these are contained in .csv files maintained with this script at https://github.com/vanmooylipidomics/LOBSTAHS/dependencies. The user can also generate the necessary .csv files from Microsoft Excel (.xlsx) workbooks of the same name, after the latter have been customized. These Excel files are maintained in the same GitHub repository.
#
#  1. "LOBSTAHS_basic_mass_matrix.csv," containing elemental compositions of adduct ions, pigments, and "base" structures for several other lipid classes.
#
#  2. "LOBSTAHS_adduct_ion_hierarchies_pos.csv," adduct ion abundance hierarchies for target lipid classes in HPLC-ESI-MS positive ion mode; data from Table 2 in Collins et al., 2015
#
#  3. "LOBSTAHS_adduct_ion_hierarchies_neg.csv," adduct ion abundance hierarchies for target lipid classes in HPLC-ESI-MS negative ion mode; data from Table 2 in Collins et al., 2015
#
#  4. "LOBSTAHS_chem_prop_ranges.csv," containing ranges of three chemical properties (total no. of acyl carbon atoms, total no. of acyl C=C double bonds, and number of possible additional oxygen atoms) for which masses will be calculated for each class of lipid
#
#  5. "LOBSTAHS_addl_oxy_ranges.csv," specifying the number of possible additional oxygen atoms to be considered on species of each lipid class during the in silico simulation; default values distributed with package (0 to 4 additional oxygen atoms) are those used in Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy, 2015, "Untargeted discovery and identification of oxidative stress biomarkers using a lipidomics pipeline for complex datasets"
#
################ Initial setup and variable definition #############

# run this section before specifying user settings, below

# # ******************************************************************
# ################ Basic user begin editing here #############
# # ******************************************************************
# 
# ################ User: define locations of data files and database #############
# 
# setwd("/Users/jrcollins/Dropbox/code/LOBSTAHS/") # first, set working directory to location of LOBSTAHS files; dependencies should be in placed a subdirectory "dependencies"
# 
# LOB_massmatrix_file = "dependencies/LOBSTAHS_basic_mass_matrix.csv" # location of .csv file containing elemental composition matrix
# 
# LOB_AIH_file_pos = "dependencies/LOBSTAHS_adduct_ion_hierarchies_pos.csv" # location of .csv file containing position ion mode adduct hierarchy data for the various classes of lipid
# 
# LOB_AIH_file_neg = "dependencies/LOBSTAHS_adduct_ion_hierarchies_neg.csv" # location of .csv file containing negative ion mode adduct hierarchy data for the various classes of lipid
# 
# LOB_sim_ranges_file = "dependencies/LOBSTAHS_chem_prop_ranges.csv" # location of .csv file containing combinations of total acyl carbon atoms (FA_total_no_C) and acyl double bonds (FA_total_no_DB) for which masses of IP-DAG, FFA, PUA, and TAG will be computed during in silico simulation
# 
# LOB_oxyrange_file = "dependencies/LOBSTAHS_addl_oxy_ranges.csv" # location of .csv file specifying the number of possible additional oxygen atoms to be considered on species of each lipid class during the in silico simulation
# 
# LOB_DB_dest_pos = "dependencies/" # optional; directory to which results of positive mode in silico simulation are to be written, if desired
# 
# LOB_DB_dest_neg = "dependencies/" # optional; directory to which results of negative mode in silico simulation are to be written, if desired
# 
# 










