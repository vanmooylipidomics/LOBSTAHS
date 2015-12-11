# LOBSTAHS_main.R
#
# Created 11/18/2015 by Jamie Collins, james.r.collins@aya.yale.edu, using Apply_VML_lipidomics_rules_MAVEN_mzroll.R as a starting point
#
# Purpose: Core script of the Lipid and Oxylipin Biomarker Screening through Adduct Hierarchy Sequences (LOBTAHS) lipidomics screening pipeline, developed by the Van Mooy Lab at Woods Hole Oceanographic Institution. Identifies compounds across samples in HPLC-ESI-MS lipid data analyzed on an Exactive Plus Orbitrap mass spectrometer. Currently, the script is written to analyze lipid data from the experiment described in Graff van Creveld et al., 2015, "Early perturbation in mitochondria redox homeostasis in response to environmental stress predicts cell fate in diatoms," ISME Journal 9:385-395. This dataset is used to demonstrate the LOBSTAHS lipidomics pipeline in Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy, 2015, "Untargeted discovery and identification of oxidative stress biomarkers using a lipidomics pipeline for complex datasets."
#
# This script:
#
#  1. Uses xcms to perform (1) peak picking and integration, (2) chromatographic alignment, and (3) nonlinear grouping across samples. Requires package "IPO" for parameter optimization.
#
#  2. Uses CAMERA to perform (1) identification and elimination of secondary isotope peaks and (2) application of initial, putative compound assignments to the data based on exact mass; assignments are made from one of two lipid-oxylipin databases generated using the script "LOBSTAHS_DB_genr.R"
#
#  3. Finally, applies two orthogonal filters to screen, annotate, and winnow the list of putative assignments. Retention time window criteria are applied to exclude any features whose retention time lies outside of the range expected for the lipid class of its assigned compound. Second, the script uses observations of the relative abundances of HPLC-ESI-MS adduct ions in samples of known composition to screen the remaining features in the dataset for complicance. Application of these criteria generally significantly reduces the number of features present in the dataset; only high-confidence compound assignments are retained. All assignments retained through this step (including those that still cannot be disambiguated from others) are annotated with a number of codes to assist the user in further data analysis.  
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
# Presumes user has installed the R packages "xcms", "CAMERA", "tools", "IPO", "snowfall", "Rmpi" with all required dependencies
#
# If multicore tasking is desired, "snowfall" also required
#
# This script the following inputs:
#
#  1. A series of .mzXML files from the same dataset, containing centroided ms1 data of a single ion mode. File conversion from the Thermo .raw format, centroiding of data, and extraction of + and - mode scans into separate files can be accomplished in batch using the script "Exactive_full_scan_process_ms1+.r", available from https://github.com/vanmooylipidomics/LOBSTAHS. The .mzXML files should be placed together in a single directory, which can be specified by the user below.
#
#  2. A lipid-oxylipin database, containing exact masses and theoretical adduct hierarchy data for a wide range of lipids, ox-lipids, and oxylipins. The necessary database can be generated for the correct ion mode using the script "LOBSTAHS_DB_genr.R" and the three matrices LOBSTAHS_basic_mass_matrix.csv, LOBSTAHS_chem_prop_ranges.csv, and LOBSTAHS_adduct_ion_hierarchies_neg.csv (or LOBSTAHS_adduct_ion_hierarchies_pos.csv). You can run the database generation script just prior to executing this script, leaving the database tables sim_results.neg and sim_results.pos right in your R environment. (Or, you can import the databases from .csv files.) In generating the database, one can simply use the ranges of chemical properties for each lipid class from Collins et al., 2015; alternatively, these ranges may be specified according to the user's needs in the Microsoft Excel spreadsheet LOBSTAHS_chem_prop_ranges.xlsx, from which a custom "LOBSTAHS_chem_prop_ranges.csv" may be generated.
#
#  3. The .additional csv file, "VML_rt_windows.csv", containing the expected retention time windows for each class of lipid, if retention time screening is desired. This file is available at https://github.com/vanmooylipidomics/LOBSTAHS. The user should alter the retention times in this matrix according to his or her particular chromatography; the Van Mooy Lab data is distributed as a template.
#
#  4. If the package IPO was used to optimize xcms (centWave) peak-picking parameters AND automatic import of the optimized settings from a .csv file is desired: The file IPO_xcmsparamfits_ ... .csv, where ... is an ISO 8601 timestamp. This file can be generated using the helper script optim_xcmsParams.R, latest version at https://github.com/vanmooylipidomics/LOBSTAHS/blob/master/optim_xcmsParams.R

################ Initial setup and variable definition #############

# run this section before specifying user settings, below

# clear workspace
# you might *not* want to do this if the database matrices you generated (sim_results.neg, sim_results.pos, or both) are still in your R environment; otherwise, you'll have to load them in from external .csv files
# 
# WS <- c(ls())
# rm(WS, list = WS)

# load required packages

library(tools) 

library(xcms)

library(CAMERA)

library(rsm)

# library(Rmpi)

# run two lines below only if IPO hasn't been installed already

# library(devtools)
# install_github("glibiseller/IPO") 

library(IPO)

library(snowfall) # if multicore tasking is desired

# ******************************************************************
################ Basic user begin editing here #############
# ******************************************************************

################ User: define locations of data files and database(s) #############

working_dir = "/Users/jrcollins/Dropbox/code/LOBSTAHS/" # specify working directory
setwd(working_dir) # set working directory to working_dir

# specify directories subordinate to the working directory in which the .mzXML files for xcms can be found; per xcms documentation, use subdirectories within these to divide files according to treatment/primary environmental variable (e.g., station number along a cruise transect) and file names to indicate timepoint/secondary environmental variable (e.g., depth)
mzXMLfiles_folder_pos = "Pt_H2O2_mzXML_ms1_pos/" 
mzXMLfiles_folder_neg = "Pt_H2O2_mzXML_ms1_neg/"

# *if* the in silico simulation results tables sim_results.neg and sim_results.pos are still in your R environment (because, say, you just ran "LOBSTAHS_DB_genr.R" and didn't clear the workspace), you don't have to specify the below; otherwise, you have to tell R where the database data is
DBfile_pos = "dependencies/LOBSTAHS_lip-oxy_DB_pos_2015-11-21T22-42-44-0300.csv" 
DBfile_neg = "dependencies/LOBSTAHS_lip-oxy_DB_neg_2015-11-21T22-42-45-0300.csv"

# if you aren't planning on running IPO to optimize centWave and/or group/retcor parameters this session, but you have some parameter values from an earlier IPO run saved in a .csv file, you can specify the file paths below

# you will be given the option later to choose which parameters to use

saved_IPO_params_centW = "IPO_xcmsparamfits_2015-12-01T18-36-10-0300.csv"
saved_IPO_params_groupretcor = "IPO_retcorGroupparamfits_2015-12-01T18-36-10-0300.csv"

# ******************************************************************
################ Basic user stop editing here #############
# ******************************************************************

################# Load in mzXML files, get xcms settings from IPO or user input #############

mzXMLfiles = list.files(mzXMLfiles_folder_pos, recursive = TRUE, full.names = TRUE)

print(paste0("Loaded ",length(mzXMLfiles)," mzXML files. Dataset consists of:"))

print(mzXMLfiles)

# # exclude any files you don't want to push through xcms (e.g., blanks); note that the blanks for the Pt H2O2 dataset (Orbi_0481.mzXML and Orbi_0482.mzXML) have already been removed
# mzXMLfiles = mzXMLfiles[-c(1,2)]

################# Allow user to specify which centWave parameter values to use, and whether to run IPO #############

readinteger <- function()
{ 
  n = readline(prompt="Specify where to obtain parameters for peak picking using findPeaks.centWave. Enter '1' to run IPO for optimization now, then use those settings. Enter '2' to use default settings specified in script. Enter '3' to read in previously optimized parameter values from a .csv file specified in the file path definitions section of the script.")
  
  if(!grepl("^[0-9]+$",n))
  {
    return(readinteger())
  }
  
  return(as.integer(n))
  
}

centWparam.source = readinteger()

if (centWparam.source==1) { # user wants to run IPO now

#####################################################################################
################# Peak-picking & creation of xcmsSet using IPO and xcms #############
#####################################################################################

################# Use IPO to optimize some xcms peak-picking parameters #############

# IPO is described in Libiseller et al., 2015, "IPO: a tool for automated optimization of XCMS parameters," BMC Bioinformatics 16:118; see https://github.com/glibiseller/IPO/blob/master/vignettes/IPO.Rmd for installation instructions

# will use IPO to optimize settings for method = centWave

# define ranges of parameters to be tested
# if single value is specified for a parameter, or centWave default is used, that parameter will not be optimized

peakpickingParameters <- getDefaultXcmsSetStartingParams('centWave')
peakpickingParameters$min_peakwidth <- c(10,20) # centerpoint is 15
peakpickingParameters$max_peakwidth <- c(40,80) # centerpoint is 60
peakpickingParameters$ppm <- 2.5 # want to set this low to avoid peak data insertion errors from centWave; IPO wants to use something like 5.5 ppm if you allow it to "optimize," but this is too high
peakpickingParameters$prefilter <- 3 # a very long optimization routine settled on a value of 2.4
peakpickingParameters$value_of_prefilter <- c(1000,10000)
peakpickingParameters$snthresh <- c(10)
peakpickingParameters$noise <- c(500)

# only going to use 4 0 uM H2O2 treatment files from the dataset for optimization routine
# seems that mzXMLfiles[c(2,5)] throw up the error:
#
#       Error in checkForRemoteErrors(val) : 
#        ... nodes produced an error
#
# so, will use the other four files (mzXMLfiles[c(1,3,4,6)])

# # code to diagnose which files contain scans that are causing the package to give error

# xcmsraw = xcmsRaw(mzXMLfiles[5])
# scan121=getScan(xcmsraw,scan=107)
# scan121[2429,]

print(paste0("Using R package IPO to optimize centWave peak-picking settings with starting parameters user has specified in script. Using following subset of files for optimization:"))

print(mzXMLfiles[c(1,3,4,6)])

resultPeakpicking <- optimizeXcmsSet(files= mzXMLfiles[c(1,3,4,6)], 
                                     params=peakpickingParameters, nSlaves=4, subdir='rsmDirectory')
optimizedXcmsSetObject <- resultPeakpicking$best_settings$xset

################# Export IPO starting value(s) and optimal settings for each parameter to .csv #############

# generate unique timestamp for filename so we don't overwrite any existing output

output_DTG = format(Sys.time(), "%Y-%m-%dT%X%z") # return current time in a good format
output_DTG = gsub(" ", "_", output_DTG) # replace any spaces
output_DTG = gsub(":", "-", output_DTG) # replaces any colons with dashes (Mac compatibility)

# write 3-column table to .csv using write.table()

write.table(cbind(sort(rownames(as.matrix(peakpickingParameters))),
as.character(resultPeakpicking$best_settings$parameters[sort(names(resultPeakpicking$best_settings$parameters))]),
as.character(peakpickingParameters[sort(rownames(as.matrix(peakpickingParameters)))])),
file = paste("IPO_centWaveparamfits_",output_DTG,".csv",sep=""),
col.names = c("centWave_parameter","IPO_optim_value","Starting_value(s)"),
row.names = FALSE,
sep=",")

print(paste0("IPO optimization complete. Optimized and starting values for centWave parameters written to file IPO_centWaveparamfits_",output_DTG,".csv"))

if (length(resultPeakpicking$best_settings$parameters)==13) { # double-check that optimized settings exist, then use them
	
	print(paste0("Using IPO-optimized settings for findPeaks.centWave..."))
	
    centW.min_peakwidth = resultPeakpicking$best_settings$parameters$min_peakwidth
    centW.max_peakwidth = resultPeakpicking$best_settings$parameters$max_peakwidth
	centW.ppm = resultPeakpicking$best_settings$parameters$ppm
	centW.mzdiff = resultPeakpicking$best_settings$parameters$mzdiff
    centW.snthresh = resultPeakpicking$best_settings$parameters$snthresh
	centW.prefilter = c(resultPeakpicking$best_settings$parameters$prefilter,resultPeakpicking$best_settings$parameters$value_of_prefilter)
	centW.noise = resultPeakpicking$best_settings$parameters$noise
	
	# not using IPO settings from resultPeakpicking$best_settings$parameters for mzCenterFun, integrate, fitgauss, verbose.columns, nSlaves since those weren't targets of optimization
	
	} else if (centWparam.source==2) { # user wants to use settings specified below
		
	    print(paste0("Using values of centWave parameters specified in the script by user..."))

		# "non-optimized" settings listed here are based on recommended "HPLC/Orbitrap settings" from Table 1 of Patti et al., 2012, "Meta-analysis of untargeted metabolomic data from multiple profiling experiment," Nature Protocols 7: 508-516
		
	centW.min_peakwidth = 10 
	centW.max_peakwidth = 45 # lowered from Patti et al. recommended HPLC setting of 60 based on visual inspection of a single sample with plotPeaks 
	centW.ppm = 2.5
	centW.mzdiff = 0.005
	centW.snthresh = 10
	centW.prefilter = c(3,7500) # 3.5k recommended by Patti et al. appears to be too low
	centW.noise = 500

		} else if (centWparam.source==3) { # user wants to read in parameter values from file
			
	print(paste0("Loading values of centWave parameters from previous IPO optimization run in .csv file ", saved_IPO_params_centW))
			
	  centWprams.from.file = read.csv(saved_IPO_params_centW,colClasses = "character")
	  
	centW.min_peakwidth = as.numeric(centWprams.from.file[centWprams.from.file[,1]=="min_peakwidth",2])
	centW.max_peakwidth = as.numeric(centWprams.from.file[centWprams.from.file[,1]=="max_peakwidth",2])
	centW.ppm = as.numeric(centWprams.from.file[centWprams.from.file[,1]=="ppm",2])
	centW.mzdiff = as.numeric(centWprams.from.file[centWprams.from.file[,1]=="mzdiff",2])
	centW.snthresh = as.numeric(centWprams.from.file[centWprams.from.file[,1]=="snthresh",2])
	centW.prefilter = c(as.numeric(centWprams.from.file[centWprams.from.file[,1]=="prefilter",2]),as.numeric(centWprams.from.file[centWprams.from.file[,1]=="value_of_prefilter",2]))
	centW.noise = as.numeric(centWprams.from.file[centWprams.from.file[,1]=="noise",2])
	
	}
		
		# specify some additional settings we wish to keep constant, regardless of where the parameters above were obtained
		
		centW.fitgauss = TRUE
		centW.sleep = 1
		centW.mzCenterFun = c("wMean")
		centW.verbose.columns = TRUE
		centW.integrate = 1
		centW.profparam = list(step=0.001) # setting this very low, per Jan Stanstrup; low setting uses more memory but helps avoid the situation where mass accuracy eclipses the actual width of the m/z windows used to define each peak (a real possibility with Orbitrap data; see http://metabolomics-forum.com/viewtopic.php?f=8&t=598#p1853) 
		centW.nSlaves = 4 # if you have r package "snow" installed, can set to number of cores you wish to make use of

# ################# Peak visualization using individual sample files #############

# # optional section for method development

# # create xcmsRaw object from just a single sample (for method development)

# xfile_raw = xcmsRaw(mzXMLfiles[1], profparam = centW.profparam)
# profStep(xfile_raw) = 0.005

# rawpeaks = findPeaks.centWave(xfile_raw,
                 # ppm = centW.ppm,
                 # peakwidth = c(centW.min_peakwidth,centW.max_peakwidth),
                 # fitgauss = centW.fitgauss,
                 # noise = centW.noise,
                 # mzdiff = centW.mzdiff,
                 # verbose.columns = centW.verbose.columns,
                 # snthresh = centW.snthresh,
                 # integrate = centW.integrate,
                 # prefilter = centW.prefilter,
                 # mzCenterFun = centW.mzCenterFun
# #                 ,sleep = centW.sleep
# #                 nSlaves = centW.nSlaves
                 # ) 

# # # despite the good press, massifquant was picking some very bad looking features, using centWave for time being instead

# # # rawpeaks = findPeaks.massifquant(xfile_raw,
                # # criticalValue = 1,
                # # consecMissedLimit = 2, # supposedly optimal for Orbitrap data
                # # prefilter = c(3,10000), # documentation says the first argument only necessary if withWave = 1, but then the example shows it there with withWave = 0; using 10k rather than the 3.5k recommended by Patti et al.
                # # ppm = 2.5, # using recommended setting of Patti et al., 2012 (see below)
                # # unions = 1,
                # # profparam = centW.profparam,
                # # withWave = 1, # two arguments immediately below are if withWave = 1 only
                # # sleep = 1,
                # # peakwidth = c(10,45), # min. feature length in time scans, max chromatographic peak width
                # # snthresh = 10,
                # # integrate = 1,
                # # checkBack = 1,
                # # fitgauss = FALSE
# # #                nSlaves = 4 # if you have r package "snow" installed, can set to number of cores you wish to make use of
                # # ) 

# # plot some selected peaks

# plotPeaks(xfile_raw,rawpeaks[10150:10174,],figs = c(5,5),width = 100)
# plotPeaks(xfile_raw,rawpeaks[150:174,],figs = c(5,5),width = 100)
# plotPeaks(xfile_raw,rawpeaks[1:24,],figs = c(5,5),width = 100)

# # N.B., just because you can't see the full extent of the peaks in some of the subplots doesn't mean they're bad; appears to be something wonky with the ylim setting in plotPeaks; see www.metabolomics-forum.com/viewtopic.php?f=8&t=875

# # for example, can look at an individual peak this way:

# plotEIC(xfile_raw, mzrange = rawpeaks[10150,c("mzmin","mzmax")], rtrange = rawpeaks[10150,c("rtmin","rtmax")]   )

# # a diagnostic plot, showing necessity of setting profparam low enough

# mz_width <- rawpeaks@.Data[,"mzmax"] - rawpeaks@.Data[,"mzmin"]
# plot(density(mz_width,adjust=0.2))

# # some other plots

# plotChrom(xfile_raw)

################# Create xcmsSet using selected settings #############

print(paste0("Creating xcmsSet object from all ",length(mzXMLfiles)," mzXML files in dataset using specified settings..."))

# create xcms xset object; runs WAY faster with multicore tasking enabled; 

xset_centWave = xcmsSet(mzXMLfiles,
                 method = "centWave",
                 profparam = centW.profparam, 
                 ppm = centW.ppm,
                 peakwidth = c(centW.min_peakwidth,centW.max_peakwidth),
                 fitgauss = centW.fitgauss,
                 noise = centW.noise,
                 mzdiff = centW.mzdiff,
                 verbose.columns = centW.verbose.columns,
                 snthresh = centW.snthresh,
                 integrate = centW.integrate,
                 prefilter = centW.prefilter,
                 mzCenterFun = centW.mzCenterFun,
#                 sleep = centW.sleep
                 nSlaves = centW.nSlaves
                 )

print(paste0("xcmsSet object xset_centWave created:"))
             
print(xset_centWave)

# Some notes:
#
#  1. If using massifquant or centWave and you are sure your input data are centroided, can ignore warning message "It looks like this file is in profile mode. [method] can process only centroid mode data !" since this is just based on a heuristic. That is, you can ignore the message if you are certain data are in centroid mode. You can verify this by opening one of your converted .mzXML files in a text reader. You should see: <dataProcessing centroided="1"></dataProcessing> (a "0" is bad)
# 
#     For more on this error, see http://metabolomics-forum.com/viewtopic.php?f=8&t=267 or https://groups.google.com/forum/#!topic/xcms/xybDDQTaQiY
#
#  2. So long as the number of peak data insertion problems is relatively low (i.e., < 100), you can safely ignore the error. Otherwise, might try lowering the ppm
#
#  3. On-the-fly plotting features (i.e., with sleep ≥ 0.001 enabled) don't appear to function properly in Mac RStudio

#####################################################################################
############# Grouping and retention time correction using IPO and xcms #############
#####################################################################################

################# Use IPO to optimize some group.density, retcor.obiwarp parameters #############

# haven't figured out how to use IPO with retcor.loess yet

# get defaults for group.density and retcor.obiwarp

retcorGroupParameters <- getDefaultRetGroupStartingParams()

# specify some specific ranges for certain parameters

retcorGroupParameters$bw = c(3,15)
retcorGroupParameters$minfrac = c(0.2,0.5)
retcorGroupParameters$minsamp = 2
retcorGroupParameters$mzwid = c(0.001,0.035)
retcorGroupParameters$profStep = c(0.01,1)

# perform optimization

resultRetcorGroup <- optimizeRetGroup(xset=xset_centWave, params=retcorGroupParameters, 
                         nSlaves=4, subdir="rsmDirectory")

################# Export IPO starting value(s) and optimal settings for each parameter to .csv #############

# generate unique timestamp for filename so we don't overwrite any existing output

output_DTG = format(Sys.time(), "%Y-%m-%dT%X%z") # return current time in a good format
output_DTG = gsub(" ", "_", output_DTG) # replace any spaces
output_DTG = gsub(":", "-", output_DTG) # replaces any colons with dashes (Mac compatibility)

# write 3-column table to .csv using write.table()

# have to remove resultRetcorGroup$best_settings$center and append it to the end of the concatenated matrix, since there's no option to specify it in retcorGroupParameters

retcorGroup.exportmat = cbind(sort(rownames(as.matrix(retcorGroupParameters))),
as.character(resultRetcorGroup$best_settings[-length(resultRetcorGroup$best_settings)][sort(names(resultRetcorGroup$best_settings[-length(resultRetcorGroup$best_settings)]))]),
as.character(retcorGroupParameters[sort(rownames(as.matrix(retcorGroupParameters)))]))

retcorGroup.exportmat = rbind(retcorGroup.exportmat,c("center",resultRetcorGroup$best_settings$center,"NA"))

write.table(retcorGroup.exportmat,
file = paste("IPO_retcorGroupparamfits_",output_DTG,".csv",sep=""),
col.names = c("retcor_or_group_parameter","IPO_optim_value","Starting_value(s)"),
row.names = FALSE,
sep=",")

print(paste0("IPO optimization complete. Optimized and starting values for group.density and retcor.obiwarp parameters written to file IPO_retcorGroupparamfits_",output_DTG,".csv"))

################# Specify appropriate settings for group.density and retcor.obiwarp #############

# check for optimized group.density and retcor.obiwarp xcms settings, or give user a chance to specify them

if (length(resultRetcorGroup$best_settings)==16) { # optimized settings exist, let's use them
	
	print(paste0("IPO-optimized settings for group.density and retcor.obiwarp exist in curent R workspace. Loading..."))
	
	obiwarp.profStep = resultRetcorGroup$best_settings$profStep
	obiwarp.response = resultRetcorGroup$best_settings$response
	obiwarp.distFunc = resultRetcorGroup$best_settings$distFunc
	obiwarp.gapInit = resultRetcorGroup$best_settings$gapInit
	obiwarp.gapExtend = resultRetcorGroup$best_settings$gapExtend
	obiwarp.factorDiag = resultRetcorGroup$best_settings$factorDiag
	obiwarp.factorGap = resultRetcorGroup$best_settings$factorGap
	obiwarp.localAlignment = resultRetcorGroup$best_settings$localAlignment
	
	density.bw = resultRetcorGroup$best_settings$bw
	density.max = resultRetcorGroup$best_settings$max
	density.minfrac = resultRetcorGroup$best_settings$minfrac
	density.minsamp = resultRetcorGroup$best_settings$minsamp
	density.mzwid = resultRetcorGroup$best_settings$mzwid
		
	} else { # IPO-optimized settings likely don't exist, give user a chance to specify the parameters
		
	    print(paste0("Could not find IPO-optimized settings for group.density and retcor.obiwarp. Using values of centWave parameters specified by the user..."))

	# obiwarp settings below are the function defaults
	
	obiwarp.center = NULL
	obiwarp.profStep = 1
	obiwarp.response = 1
	obiwarp.distFunc = "cor_opt"
	obiwarp.gapInit = NULL
	obiwarp.gapExtend = NULL
	obiwarp.factorDiag = 2
	obiwarp.factorGap = 1
	obiwarp.localAlignment = 0
	
	# settings for group.density below are based on the recommended HPLC/Orbitrap settings from Table 1 of Patti et al., 2012, "Meta-analysis of untargeted metabolomic data from multiple profiling experiment," Nature Protocols 7: 508-516
	
	density.bw = 5 # 15?
	density.max = 50
	density.minfrac = 0.25
	density.minsamp = 2
	density.mzwid = 0.015 # 0.001?

		}
				
		# specify some additional settings that weren't optimized above, or for which we want to set at a specific value regardless of code above 
		
		obiwarp.center = NULL
		obiwarp.plottype = "deviation" # "none"
		density.sleep = 1

################# Perform grouping and retention time correction on dataset #############

print(paste0("Performing grouping and retention time correction on dataset..."))

# chromatographic alignment (retention time correction)

# # using method loess

# retcord_xset = retcor(xset_centWave,
                      # missing = 1,
                      # extra = 1,
                      # smooth = "loess",
                      # span = .2,
                      # family = "gaussian", # want to leave outliers in for the time being
			          # plottype = "mdevden"
                      # )

# using method obiwarp

retcord_xset = retcor(xset_centWave,
             method = "obiwarp",
             plottype = obiwarp.plottype,
             profStep = obiwarp.profStep,
             center = obiwarp.center,
             response = obiwarp.response,
             distFunc = obiwarp.distFunc,
             gapInit = obiwarp.gapInit,
             gapExtend = obiwarp.gapInit,
             factorDiag = obiwarp.factorDiag,
             factorGap = obiwarp.factorGap,
             localAlignment = obiwarp.localAlignment,
             initPenalty = 0
             )
             
# grouping

# # method "nearest" with settings below seems to work better than method = "density," but takes absolutely forever; however, it seems to take less time crunching centWave picked data than massifquant picked data

# xset_centWave = group(xset_centWave,
             # method = "nearest",
             # mzVsRTbalance=10,
             # mzCheck=0.2,
             # rtCheck=30,
             # kNN=10
             # ) 

# using method = "density" with settings from above

xset_centWave = group(xset_centWave,
			 method = "density",
             bw = density.bw,
             minfrac = density.minfrac,
             minsamp = density.minsamp,
             mzwid = density.mzwid,
             max = density.max,
             sleep = density.sleep
             )







# re-perform grouping

retcord_xset = group(retcord_xset,
             method = "nearest",
             mzVsRTbalance=10,
             mzCheck=0.2,
             rtCheck=30,
             kNN=10
             ) 
                     
retcord_xset2 = group(retcord_xset2,
                     method = "density",
                     bw = 5,
                     minfrac = .25,
                     minsamp = 2,
                     mzwid = 0.015,
                     max = 50
                     )


# fill missing peaks

xset_process = fillPeaks(retcord_xset,
                                 nSlaves = 4)
                                 
xset2_process = fillPeaks(retcord_xset2,
                                 nSlaves = 4)


################# Isotope peak annotation using CAMERA  #############

# first, a necessary workaround to avoid a import error; see https://support.bioconductor.org/p/69414/
imports = parent.env(getNamespace("CAMERA"))
unlockBinding("groups", imports)
imports[["groups"]] = xcms::groups
lockBinding("groups", imports)

# create annotated xset

xseta = xsAnnotate(xset_process)
xset2a = xsAnnotate(xset2_process)

# create CAMERA groups by RT; as of 11/26/15, not sure what settings to use (or if defaults are ok), or how this interacts with the groups already created in xcms

xsetaF = groupFWHM(xseta)
xsetaF = groupDen(xseta, bw = 1)


xset2aF = groupFWHM(xset2a) 

# find, annotate isotope peaks

xsetaFI = findIsotopes(xsetaF,
                       mzabs = 0.001,
                       maxcharge = 4,
                       maxiso = 4,
                       ppm = 2.5,
                       intval = c("maxo"),
                       minfrac = 0.5)
                       
xset2aFI = findIsotopes(xset2aF,
                       mzabs = 0.001,
                       maxcharge = 1,
                       ppm = 2.5,
                       intval = c("into"),
                       minfrac = 0.25)

write.csv(file="peaklist_with_isotopes.csv",getPeaklist(xsetaFI))
write.csv(file="isotopelist.csv",getIsotopeCluster(xsetaFI))


write.csv(file="peaklist_process.csv",getPeaklist(xseta))

xsetaC = groupCorr(xsetaF)
xset2aC = groupCorr(xset2aF) 

iso = getIsotopeCluster(xsetaFI)

xset_a_Fiso = groupFWHM(xset_a) 

################# Load in entries from relevant database #############



################# Load in entries from relevant database #############

print(paste0("Loading in database data..."))

lipid_DB <- read.csv(lipid_DB_file, skip = 0, col.names=c("frag_ID","mz_plus_H","mz","neutral_mass","lipid_class", "degree_oxidation","adduct","adduct_relative_percent","adduct_rank","FA_total_no_C","FA_total_no_DB","MAVEN_fullstring"), stringsAsFactors=FALSE)

lipid_DB <- lipid_DB[,1:11] # get rid of the MAVEN full string column since we don't need it in this script

# clean up any ?'s in the adduct_relative_percent or adduct_rank fields in the DB; these may exist in older versions of the VML databases there were cases (e.g., for some PUAs in - mode DB, for TAGs and pigments in + mode DB) where we hadn't yet populated this information for certain molecular species, either because we didn't have any standard data or we hadn't gotten to them yet 

lipid_DB$adduct_relative_percent[lipid_DB$adduct_relative_percent=="?"] <- NA  
lipid_DB$adduct_rank[lipid_DB$adduct_rank=="?"] <- NA  

# convert adduct_relative_percent and adduct_rank to numeric

lipid_DB$adduct_rank <- as.numeric(lipid_DB$adduct_rank)
lipid_DB$adduct_relative_percent <- as.numeric(lipid_DB$adduct_relative_percent)

# clean up known duplicate adduct levels in the DB

lipid_DB$adduct <- as.character(lipid_DB$adduct)
lipid_DB$adduct[lipid_DB$adduct==c("M+acetate")] <- c("M+Acetate") 
lipid_DB$adduct[lipid_DB$adduct==c("M+Chloride")] <- c("M+Cl") 

################# Import data from the MAVEN .mzroll file, pulling in necessary data from the lipid database as we go #############

print(paste0("Parsing XML in source .mzroll file..."))

parsedXML <- xmlTreeParse(MAVEN_mzroll_datafile) # parse the XML in the MAVEN mzroll file
XML_data <- xmlRoot(parsedXML) # extract data in XML format
XML_groupdata <- XML_data[[3]] # retain only the third node, which contains the actual data we want

numgroups <- xmlSize(XML_groupdata) # get the number of groups (i.e., adducts in the database) putatively assigned to various LC-MS features in the dataset by MAVEN; at this point, this will include overlapping/competing assignments; MAVEN functions by assigning any number of features to a specific group

numpeaks <- sum(xmlSApply(XML_groupdata, xmlSize)) # get the total number of mass spectral features ID'd by MAVEN using whatever peak picking/matching settings you specified

# create two destination arrays for the imported data: one array for group (i.e., compound assignment) data and one for peak (feature) data

MAVEN_groupdata <- matrix(NA,numgroups,27) # create matrix, 12 columns for data from the XML import plus 4 extras to keep track of how many peaks (features) we have attached to each group as we move throughout the script, plus 11 columns for data we extract for each group from the database
MAVEN_groupdata <- as.data.frame(MAVEN_groupdata) # convert to data frame so we can populate with both numeric and non-numeric data

colnames(MAVEN_groupdata) <- c("groupID","tagString","metaGroupId","expectedRtDiff","groupRank","label","type","changeFoldRatio","changePValue","compoundId","compoundDB","compoundName","peaksthisgroup_init","peaksthisgroup_quality","peaksthisgroup_rules","peaksthisgroup_final","frag_ID","DB_mz_plus_H","DB_mz","DB_neutral_mass","DB_lipid_class","DB_degree_oxidation","DB_adduct","DB_adduct_relative_percent","DB_adduct_rank","DB_FA_total_no_C","DB_FA_total_no_DB") # set column names

MAVEN_groupdata[,c("peaksthisgroup_init","peaksthisgroup_quality","peaksthisgroup_rules","peaksthisgroup_final")] <- 0 # set some fields to default value of 0

MAVEN_peakdata <- matrix(NA,numpeaks,37) # 36 columns of data per peak (feature) + 1 extra in which to note the groupID; the groupID will allow us to enter the first array with the second and retrieve group data as we need it
MAVEN_peakdata <- as.data.frame(MAVEN_peakdata) # convert to data frame so we can populate with both numeric and non-numeric data

colnames(MAVEN_peakdata) <- c("groupID","pos","minpos","maxpos","rt","rtmin","rtmax","mzmin","mzmax","scan","minscan","maxscan","peakArea","peakAreaCorrected","peakAreaTop","peakAreaFractional","peakRank","peakIntensity","peakBaseLineLevel","peakMz","medianMz","baseMz","quality","width","gaussFitSigma","gaussFitR2","groupNum","noNoiseObs","noNoiseFraction","symmetry","signalBaselineRatio","groupOverlap","groupOverlapFrac","localMaxFlag","fromBlankSample","label","sample") # set column names

# now, cycle through the XML structure, pull out data, format it, and populate our destination arrays

c_i_MAVEN_peakdata = 1 # set a counter that will allow us to insert peak (feature) data in the right place

for (i in 1:xmlSize(XML_groupdata)) { # cycle through assignment groups
  
  thisgroup <- XML_groupdata[[i]] # subset XML to grab group and peak data for just this group 
  
  # group data
  
  groupdata_thisgroup <- as.character(xmlAttrs(thisgroup, addNamespacePrefix = F)) # extract group data as character (otherwise will pull as factors, and we don't want this)
  thisgroupID <- as.numeric(groupdata_thisgroup[1]) # pull out the group ID for this group  
  thisfragID <- as.integer(gsub(" lipid.+","",groupdata_thisgroup[10])) # get an integer vector of the fragment ID for this adduct as it appears in the database, which we will then use to query the database and extract the information present in the DB's individual variables
  numpeaks_thisgroup <- xmlSize(thisgroup) # get number of peaks (features) associated with this group, i.e., the number of features that have this particular assignment attached to them
  MAVEN_groupdata[i,13] <- as.numeric(numpeaks_thisgroup) # record number of peaks initially associated with this group
  MAVEN_groupdata[i,1:9] <- as.numeric(groupdata_thisgroup[1:9]) # store numeric group data into the destination groupdata array
  MAVEN_groupdata[i,10:12] <- groupdata_thisgroup[10:12] # store non-numeric group data
  MAVEN_groupdata$frag_ID[i] <- thisfragID # store fragID
    
  # use fragID to retrieve species data from the database, then insert relevant database data
  
  DB_data_for_this_group <- lipid_DB[thisfragID==lipid_DB$frag_ID,c("mz_plus_H","mz","neutral_mass","lipid_class","degree_oxidation","adduct","adduct_relative_percent","adduct_rank","FA_total_no_C","FA_total_no_DB")] # extract DB data for this frag_ID
  # flow extracted data into the MAVEN_groupdata matrix:
  MAVEN_groupdata[i,c("DB_lipid_class","DB_degree_oxidation","DB_adduct","DB_adduct_relative_percent")] <- c(as.character(DB_data_for_this_group$lipid_class),as.character(DB_data_for_this_group$degree_oxidation),as.character(DB_data_for_this_group$adduct),as.character(DB_data_for_this_group$adduct_relative_percent)) # apparently, have to insert factor variables using as.character, then convert later back to factors or numbers
  MAVEN_groupdata[i,c("DB_mz_plus_H","DB_mz","DB_neutral_mass","DB_adduct_rank","DB_FA_total_no_C","DB_FA_total_no_DB")] <- DB_data_for_this_group[c("mz_plus_H","mz","neutral_mass","adduct_rank","FA_total_no_C","FA_total_no_DB")] # insert numerical data
  
  # peak (feature) data
  
  peakdata_thisgroup <- t(as.data.frame(xmlSApply(thisgroup,xmlToList,addAttributes = FALSE, simplify = TRUE))) # extract peak data for this group from XML structure
  endrow <- c_i_MAVEN_peakdata+numpeaks_thisgroup-1 # calculate the row in MAVEN_peakdata where we should finish inserting peakdata for this group
  MAVEN_peakdata[c_i_MAVEN_peakdata:endrow,1:36] <- cbind(thisgroupID,matrix(as.numeric(peakdata_thisgroup[,1:35]), ncol = 35)) # insert numeric peak data for these peaks (including the groupID) into our destination peakdata array
  MAVEN_peakdata[c_i_MAVEN_peakdata:endrow,37] <- peakdata_thisgroup[,36] # insert non-numeric peak data for these peak
  
  c_i_MAVEN_peakdata <- c_i_MAVEN_peakdata + numpeaks_thisgroup # increment our counter
  
}

rm(c_i_MAVEN_peakdata) # delete our counter

################# Some initial housekeeping after the initial (big) import #############

# consolidate some duplicate adduct levels we know can exist in MAVEN_groupdata because of a typo in early versions of the databases

MAVEN_groupdata$DB_adduct[MAVEN_groupdata$DB_adduct==c("M+acetate")] <- c("M+Acetate") 
MAVEN_groupdata$DB_adduct[MAVEN_groupdata$DB_adduct==c("M+Chloride")] <- c("M+Cl") 

# convert some fields to factors

MAVEN_groupdata$DB_adduct <- as.factor(MAVEN_groupdata$DB_adduct)
MAVEN_groupdata$DB_lipid_class <- as.factor(MAVEN_groupdata$DB_lipid_class)
MAVEN_groupdata$DB_degree_oxidation <- as.factor(MAVEN_groupdata$DB_degree_oxidation)

################# Quality and peak area (intensity) control #############

# This is where user-defined settings (specified above) are applied

# discard all peaks (features) with a MAVEN-assigned quality score less than user's min_qual cutoff (as specified above)

MAVEN_goodpeaks <- MAVEN_peakdata[MAVEN_peakdata$quality>=min_qual,]

# discard all peaks with a MAVEN-assigned peakAreaCorrected (intensity) less than user's min_peakarea cutoff (as specified above)

MAVEN_goodpeaks <- MAVEN_goodpeaks[MAVEN_goodpeaks$peakAreaCorrected>=min_peakarea,]

# now, we will probably need to eliminate some groups (assignments) and peaks (features) as a result of the criteria we've just applied:
# we will check to see if (1) there are any putative assignments (groups) for which < min_peakspergroup peaks remain in MAVEN_goodpeaks and (2) whether there are any assignments (groups) in MAVEN_groupdata for which application of our criteria eliminated all assigned features in MAVEN_peakdata

badgroups <- c(MAVEN_groupdata$groupID[MAVEN_groupdata$peaksthisgroup_init<min_peakspergroup],setdiff(MAVEN_groupdata$groupID,MAVEN_goodpeaks$groupID)) # get a combined list of bad groups

# delete these groups from MAVEN_goodpeaks and delete any remaining peaks associated with these groups from MAVEN_peakdata

MAVEN_goodpeaks <- MAVEN_goodpeaks[!(MAVEN_goodpeaks$groupID %in% badgroups),]
MAVEN_goodgroups <- MAVEN_groupdata[!(MAVEN_groupdata$groupID %in% badgroups),]

# calculate how many peaks still remain in the dataset for each group at this point (i.e., the number of features having a particular assignment from the databse), and record

num_peaks_per_group.now <- as.data.frame(table(MAVEN_goodpeaks$groupID)) # count number of peaks per group, right now

MAVEN_goodgroups[which(MAVEN_goodgroups$groupID %in% as.numeric(as.character(num_peaks_per_group.now[,1]))),c("peaksthisgroup_quality")] <- num_peaks_per_group.now[,2] # record this information to MAVEN_goodgroups

################# Other necessary housekeeping, and provide some output to user #############

print(paste0("After applying quality score cutoff of > ",min_qual,", a minimum peakAreaCorrected of ",min_peakarea,", and the requirement that at least ",min_peakspergroup," peaks (i.e., features) in the dataset were given the same assignment from the database, ",nrow(MAVEN_goodpeaks)," peaks were retained out of a total of ",nrow(MAVEN_peakdata)," features present in the dataset imported from MAVEN."))

print(paste0("These ",nrow(MAVEN_goodpeaks)," retained features represent ",nrow(MAVEN_goodgroups)," potential adducts of ",nrow(unique(MAVEN_goodgroups[c("DB_lipid_class","DB_degree_oxidation","DB_FA_total_no_C","DB_FA_total_no_DB")]))," unique parent compounds found in ",length(unique(as.character(MAVEN_goodpeaks$sample)))," samples in this experiment."))

print(paste0("This script will further winnow this list of ",nrow(MAVEN_goodpeaks)," features to retain only those meeting the adduct hierarchy rules. The script will also annotate peaks features which valid competing assignments were made, i.e., in cases where different adducts of two or more molecular species share the same m/z, or have different m/z's that fall within the ppm tolerance used for matching. Application of these rules may reduce the total number of unique potential parent compounds in the dataset to less than ",nrow(unique(MAVEN_goodgroups[c("DB_lipid_class","DB_degree_oxidation","DB_FA_total_no_C","DB_FA_total_no_DB")]))," if any the apparent adduct distributions for the parent compounds do not satisfy the hierarchy rules."))             

# also, at this point we will capture the meaning of the levels in our three factor variables, because the values of these variables will be converted to numbers as we load them into the array

lipid_class_levels <- levels(MAVEN_goodgroups$DB_lipid_class)
degree_oxidation_levels <- levels(MAVEN_goodgroups$DB_degree_oxidation)
adduct_levels <- levels(MAVEN_goodgroups$DB_adduct)

# define some lists

samples <- unique(MAVEN_peakdata$sample) # create a vector of all sample_IDs in the data
samples <- factor(samples, levels=as.character(samples)) # convert to factors

casecodes <- c("C1","C1x","C2a","C2b","C3c","C3f","C3r","C4","C5","C6a","C6b") # the possible case codes with which we'll be annotating each putative parent compound identification based on the data for its various adduct assignments

# N.B. on case codes used for annotation: All codes as defined in Collins et al., 2015 (Fig. 2), with addition of C1x, C6a and C6b, which are attached in instances where the database contains incomplete (C6a) or no (C6b) adduct hierarchy information for a given parent compound. C1x is assigned in addition to C6b when the DB contains no hierarchy information AND only one adduct is present, making it impossible to determine whether the assignment should be annotated as C1 or C4; this should be rare, since newest versions of the databases are nearly complete.

numcodes <- length(casecodes) # store number of case codes in a variable

################# Prep our peakdata output array screened_peakdata, create a vector current_casecodes we'll use to keep track of the casecodes for each assignment as we go along, and create an array case3s.bygroupID to hold groupIDs of case 3c and 3f compounds #############

# we will record most screened peak data to a single array screened_peakdata as we go; later, we'll go back and generate a results table by parent compound and sample ID

screened_peakdata <- array(NA,c(1,ncol(MAVEN_goodpeaks)+length(casecodes))) # create array to be populated with screened peakdata; setting first dimension size to nrow(MAVEN_goodpeaks) to give plenty of room for recording peak data (can easily truncate the extra rows of NA's at the end of the routine)

screened_peakdata <- as.data.frame(screened_peakdata) # convert to data frame

dimnames(screened_peakdata)[[2]] <- c(colnames(MAVEN_goodpeaks),casecodes) # label columns
screened_peakdata[,38:(38+numcodes-1)] <- 0 # set all case codes to 0 by default

current_casecodes <- array(NA,c(length = numcodes)) # create a vector to keep track of case codes
current_casecodes[1:numcodes] <- 0 # set all case codes to 0 by default
names(current_casecodes) <- casecodes # label columns

if (exists("peakdata_to_record")) {
rm(peakdata_to_record) # clear our peakdata to record vector, in case it's still populated with old data
}

# create a separate array case3s.bygroupID to hold the groupIDs of case 3c and 3f assignments (i.e., those for which there is a competing or overlapping assignment)
# assume we won't have more than 20 possible overlapping assignments in a given sample

case3s.bygroupID <- as.data.frame(array(NA,c(1,21))) # create array for case 3c/3f groupID data
colnames(case3s.bygroupID)[1] <- c("sample") # set field name of first column to sample

case3s.bygroupID_to_record <- as.data.frame(array(NA,c(1,21))) # create placeholder array for group 3c/3f data from each sample

################# Application of hierarchy rules #############

# at this point, we can finally apply our rules to the data

# we will use the adduct hierarchy information extracted from the database to evaluate each initial assignment MAVEN made to the data

# assuming we must observe the expected adduct hierarchy in each sample as it was run on the mass spectrometer, we'll apply the rules to the assignments in each sample independently

for (i in 1:length(samples)) { # cycle through each sample ID
  
  print(paste0("Now conducting initial application of adduct hierarchy rules to data for sample ",samples[i],"..."))
  
  MAVEN_peaks_this_sample <- MAVEN_goodpeaks[MAVEN_goodpeaks$sample==samples[i],] # get list of peaks (features) MAVEN found in this sample
  MAVEN_groups_this_sample <- MAVEN_goodgroups[MAVEN_goodgroups$groupID %in% MAVEN_peaks_this_sample$groupID,] # get list of groups (i.e., assignments from DB) associated with these peaks
  
  # in each sample, we will have two basic cases for assignments : (1) parent compounds for which MAVEN made only one assignment (i.e., case C1 or C4) and (2) parent compounds for which MAVEN made > 1 assignment, either by identifying a number of different adducts for this same compound (i.e., case C2a, C2b, or C5) and/or because MAVEN ID'd > 1 peak of the same adduct at different retention times (i.e., case 3r); since the second broad category will be trickier, we'll deal with them first  
  
  # we will also have some case 3c and 3f scenarios which can apply to compounds bearing any of the other annotations: the mz of an adduct assigned to one feature matches the mz of an adduct of a different parent compound assigned to some other feature; for these instances, we will want to note this situation but not eliminate these compounds from the dataset automatically (though annotating as we do here will allow a user to easily remove them later, if desired)
    
  # lastly we'll have a case 6 scenario that can occur in conjunction with any of the above: we do not have complete (or any) information in the DB on the relative abundances of the adducts for a particular species, but MAVEN has determined that one or more of them is present in the sample. we will have a "case 6a" scenario where rank information for some, but not all, possible adducts is present in the DB; and a "case 6b" scenario where we don't have rank information for any of the possible adducts. we will be checking for case 6 assignments at various points throughout the script. In a case 6b scenario, the script will report data for the adduct of the species which is most abundant (in terms of peak area) in that particular sample.
  
  ################# resolution of case 2a, case 2b, case 3r, and case 5 #############
  
  case_2_molec_species <- unique(MAVEN_groups_this_sample[duplicated(MAVEN_groups_this_sample[c("DB_lipid_class","DB_degree_oxidation","DB_FA_total_no_C","DB_FA_total_no_DB")]),c("DB_lipid_class","DB_degree_oxidation","DB_FA_total_no_C","DB_FA_total_no_DB")]) # first, we need a list of the case 2a/2b/3r/5 assignments in this sample (i.e., those parent compounds for which MAVEN ID'd more than one different adduct)
  
  if (nrow(case_2_molec_species) > 0) {  # if we don't have any case 2a/2b/3r/5 species, skip this part
  
  print(paste0("   Applying hierarchy rules to ",nrow(case_2_molec_species)," case 2a/2b/3r/5 parent compounds putatively present in this sample, i.e., those compounds for which more than one adduct was identified..."))
  
  retcount = 0 # create a counter to keep track of how many species of each kind we end up retaining

  for (j in 1:nrow(case_2_molec_species)) { # cycle through each of the case 2a/2b/3r/5 compounds represented in this sample
           
   groups_this_molec_species <- MAVEN_groups_this_sample[which(MAVEN_groups_this_sample$DB_lipid_class==case_2_molec_species$DB_lipid_class[j] & MAVEN_groups_this_sample$DB_degree_oxidation==case_2_molec_species$DB_degree_oxidation[j] & MAVEN_groups_this_sample$DB_FA_total_no_C==case_2_molec_species$DB_FA_total_no_C[j] & MAVEN_groups_this_sample$DB_FA_total_no_DB==case_2_molec_species$DB_FA_total_no_DB[j]),] # return group data corresponding to assignments in this sample that involve this parent compound
   
   peaks_this_molec_species <- MAVEN_peaks_this_sample[MAVEN_peaks_this_sample$groupID %in% groups_this_molec_species$groupID,] # return peak data for all assignments in this sample that involve this parent compound
  
  adduct_distribution <- table(groups_this_molec_species$DB_adduct) # use the assignments involving this parent compound to obtain the distribution of its adducts in this sample
  
  possible_adducts <- lipid_DB[lipid_DB$lipid_class==as.character(case_2_molec_species$DB_lipid_class[j]) & lipid_DB$degree_oxidation==as.character(case_2_molec_species$DB_degree_oxidation[j]) & lipid_DB$FA_total_no_C==as.character(case_2_molec_species$DB_FA_total_no_C[j]) & lipid_DB$FA_total_no_DB==as.character(case_2_molec_species$DB_FA_total_no_DB[j]),c("adduct","adduct_rank")] # for comparison purposes, get a list of all possible adducts for this parent compound from the DB
  
  possible_adducts <- possible_adducts[order(possible_adducts$adduct_rank, decreasing = TRUE, na.last = FALSE),] # reorder the list of possible adducts to iterate through them; must return any NA's (i.e., any adducts for which ranking was not given in DB) first so that rest of code works ok
  
  for (k in 1:nrow(possible_adducts)) { # iterate through the list of all possible adducts from least to most abundant, checking to see whether MAVEN made an assignment for each one
    
    if (adduct_distribution[as.character(possible_adducts$adduct[k])]>=1) { # at least one peak was identified in the data for this adduct
      
      possible_adducts$present_in_sample[k]=1 # mark this adduct as being present in the data
      possible_adducts$num_times_present[k]=adduct_distribution[as.character(possible_adducts$adduct[k])] # indicate number of different peaks MAVEN identified as this adduct
      
    } else { # assume no peak was identified for this adduct
      
      possible_adducts$present_in_sample[k]=0 # mark this adduct as being absent from the data
      possible_adducts$num_times_present[k]=0 # indicate number of different peaks MAVEN identified as this adduct (0) 
      
    }
    
  }
  
  ################# case 6 check #############
  
  # before beginning evaluation of case 2a/2b/3r/5 parent compounds, first check to see whether they happen to be case 6's, then proceed accordingly
  
  if (!all(is.na(possible_adducts$adduct_rank))) { # either we have a case 6a situation, or no case 6 at all; either way, proceed with evaluation of case 2a/2b/3r/5 parent compounds
    
    if (any(is.na(possible_adducts$adduct_rank))) { # we have a case 6a situation
      
      current_casecodes["C6a"] <- 1 # note that case 6a is true
      
    }
  
  # apply subrules to the case 2a/2b/3r/5 assignments for this parent compound
  # Note that the way i've coded this, using independent "if" statments, some subrule assignments will be nonexclusive; i.e., the script can tag a particular parent compound as both case 3r and 2a or 2b if it meets the minimum criteria for each
  
  # case 3r: we have multiple assignments for the same putative parent compound, and this is because MAVEN identified > 1 peak representing the same adduct at two or more different retention times
    
  if (any(adduct_distribution>1)) { # parent compound is case 3r
    
    current_casecodes["C3r"] <- 1 # note that case 3r is true
        
  }

  # case 2a, 2b, or 5: we have multiple assignments for the same putative parent compound, and this is because MAVEN identified multiple peaks representing different adducts
  # note that for case 2a/2b determination, we will only consider adduct_rank; the logic in the script doesn't at this point consider the relative abundances (%'s) of the different adducts 

  if (any(adduct_distribution==1)) { # putative parent compound is case 2a/2b/5... but, have to determine which subtybe

  # now, use adduct information to determine whether 2a/2b/5; easiest to begin with 5 and then address 2a and 2b

  # case 5: the different adducts identified by MAVEN do NOT satisfy the hypothesized hierarchy; adducts of hypothetically lower rank are present while the abundant adduct which should be most abundant is not present
  
  if ((sum(possible_adducts$present_in_sample[1:nrow(possible_adducts)-1]) >= 1) & (possible_adducts$present_in_sample[nrow(possible_adducts)]==0)) { # putative parent compound is case 5 --> bad assignment
    
    current_casecodes["C5"] <- 1 # note that case 5 is true
    
    } else { # we should have case 2a or 2b
      
    # case 2a: the different adducts identified by MAVEN strictly satisfy the hypothesized adduct hierarchy for this parent compound, i.e., MAVEN must have identified a peak in the sample for every adduct more abundant than the least abundant adduct in the sample (however, this does not require MAVEN to have identified a peak in the sample for EVERY possible adduct of this parent compound)
      
    adduct_theoretically_least <- match(1,possible_adducts$present_in_sample) # of the adducts present in the sample, this one should be the least abundant according to the rules; for case 2a to be satisfied, we must then have in the sample all adducts which are supposed to be more abundant than this one
    
    if (sum(possible_adducts$present_in_sample[adduct_theoretically_least:nrow(possible_adducts)]) == nrow(possible_adducts) - adduct_theoretically_least + 1 ) { # parent compound is case 2a
      
      current_casecodes["C2a"] <- 1 # note that case 2a is true
      
    }
    
    # case 2b: the different adducts identified by MAVEN do not perfectly satisfy the hypothesized adduct hierarchy for this putative parent compound, but we consider it a good parent compound match because the adduct which should be most abundant according to the rules is present in the sample; in this case, some adducts of intermediate hypothesized abundance may be absent from the sample while some other less abundant adducts are present, however the adduct which should be most abundant is definitely present
    
    if ((sum(possible_adducts$present_in_sample[adduct_theoretically_least:nrow(possible_adducts)]) < nrow(possible_adducts) - adduct_theoretically_least + 1 ) & possible_adducts$present_in_sample[nrow(possible_adducts)]==1) { # parent compound is case 2b
      
      current_casecodes["C2b"] <- 1 # note that case 2b is true
      
    }
    
    }
  
  }
  
  # last order of business before moving onto the next case 2a/2b/3r/5 species in this sample: select and then write to the results array the appropriate data for the current parent compound
  # method: we will use the rules to determine the particular adduct of this parent compound from which we will pull the data (logic as agreed on by HFF and JRC in discussion on 1/13/2015)
  # if the parent compound is case 3r but the adduct which should be most abundant is not present, or the parent compound is case 5, assume the assignment was wrong to begin with --> record no data
  # if the parent compound is case 2a or 2b, the adduct which should be most abundant is present in sample --> record data for this adduct, regardless of whether it is actually the most abundant in this particular sample; as long as we are consistent in applying this throughout the experiment, we should be ok since we're only really concerned with relative changes between samples
  # if the parent compound is case 3r and the adduct which should be most abundant is present, we have a slightly more complicated situation (could be because >= 2 regioisomers of the species were simultaneously identified in the sample); we will record data for all of these assignments by inserting as many additional row(s) as is necessary
  # lastly, if the parent compound is case 6b, we'll record data for the adduct of the parent compound that is present in the sample in highest actual abundance
  
  if ((current_casecodes["C3r"]==1 & possible_adducts$present_in_sample[nrow(possible_adducts)]==0) | current_casecodes["C5"]==1) { # case 3r but the adduct which should be most abundant is not present, or case C5 --> record no peak data
  
    # do nothing!
    # should return NULL to console
    
  } else if ((current_casecodes["C2a"]==1 | current_casecodes["C2b"]==1) & possible_adducts$num_times_present[nrow(possible_adducts)]==1) { # 2a or 2b, and MAVEN only identified a single peak for the adduct of theoretical greatest abundance  --> use data for this single peak assignment
        
    peakdata_to_record <- peaks_this_molec_species[peaks_this_molec_species$groupID==groups_this_molec_species[as.character(possible_adducts$adduct[nrow(possible_adducts)])==as.character(groups_this_molec_species$DB_adduct),c("groupID")],] # select peak data to record
            
  } else if ((current_casecodes["C3r"]==1 & possible_adducts$num_times_present[nrow(possible_adducts)]>1)) { # C3r, with more than one assignment at different retention times for the adduct that should most abundant
   
    peakdata_to_record <- peaks_this_molec_species[peaks_this_molec_species$groupID %in% groups_this_molec_species[as.character(possible_adducts$adduct[nrow(possible_adducts)])==as.character(groups_this_molec_species$DB_adduct),c("groupID")],] # select peak data to record
    
  }
  
  } else { # we have a case 6b situation: this will require some evaluation to determine which data we'll record
    
    current_casecodes["C6b"] <- 1 # note that case 6b is true
    
    adduct_maxabund <- groups_this_molec_species[groups_this_molec_species$groupID==peaks_this_molec_species$groupID[peaks_this_molec_species$peakIntensity==max(peaks_this_molec_species$peakIntensity)],c("DB_adduct")] # return the adduct of this parent compound which is most abundant in the sample
    
    # now, check to see how many peaks we have for this particular adduct in this sample; this will determine what data we record
    
    if (adduct_distribution[adduct_maxabund]==1) { # only one peak for this adduct
      
      peakdata_to_record <- peaks_this_molec_species[peaks_this_molec_species$groupID==groups_this_molec_species[as.character(adduct_maxabund)==as.character(groups_this_molec_species$DB_adduct),c("groupID")],] # select peak data to record
              
    } else if (adduct_distribution[adduct_maxabund]>1) { # we have more than one peak of this particular adduct in this sample, essentially a case 3r-case 6b situation
      
      current_casecodes["C3r"] <- 1 # at this point, we should note that case 3r is true since it wouldn't have been captured earlier
      
      peakdata_to_record <- peaks_this_molec_species[peaks_this_molec_species$groupID %in% groups_this_molec_species[as.character(adduct_maxabund)==as.character(groups_this_molec_species$DB_adduct),c("groupID")],] # select peak data to record
    
    }
          
  }
  
  # now, insert/append the peakdata and case codes for this parent compound into our peakdata results array
  
  if (exists("peakdata_to_record")) { # only force insertion if there's data
    
    retcount = retcount + nrow(peakdata_to_record) # increment our counter
    
  if (all(is.na(screened_peakdata[1,1:36]))) { # it's the first entry
    
    screened_peakdata <- as.data.frame(c(peakdata_to_record[,1:37],current_casecodes))
    
  } else { # it's not the first sample, so rbind
    
    screened_peakdata <- rbind(screened_peakdata,as.data.frame(c(peakdata_to_record[,1:37],current_casecodes)))
    
  }
  
  # finally, reset/remove our casecode and peakdata_to_record vectors before proceeding
  
  rm(peakdata_to_record)
  
  }
  
  current_casecodes[1:numcodes] <- 0
  
  }
  
  print(paste0("       Retained ",retcount,"."))
  
  }
    
  ################# resolution of case 1 and case 4 #############
  
  case_1_molec_species <- rbind(MAVEN_groups_this_sample[c("DB_lipid_class","DB_degree_oxidation","DB_FA_total_no_C","DB_FA_total_no_DB")],case_2_molec_species[,1:4])[!duplicated(rbind(MAVEN_groups_this_sample[c("DB_lipid_class","DB_degree_oxidation","DB_FA_total_no_C","DB_FA_total_no_DB")],case_2_molec_species[,1:4]),fromLast = FALSE)&!duplicated(rbind(MAVEN_groups_this_sample[c("DB_lipid_class","DB_degree_oxidation","DB_FA_total_no_C","DB_FA_total_no_DB")],case_2_molec_species[,1:4]),fromLast = TRUE),]  # first, we need a list of the case 1/4 putative parent compounds in this sample (i.e., those parent compounds for which MAVEN ID'd only one adduct; should be all those parent compounds that we did not ID as case 2a/2b/3r/5 above)
  
  if (nrow(case_1_molec_species) > 0) {  # if we don't have any case 1 or 4 species, skip this part
    
  print(paste0("   Applying hierarchy rules to ",nrow(case_1_molec_species)," case 1 or case 4 parent compounds putatively present in this sample, i.e., those compounds for which only one adduct was identified..."))
  
  retcount = 0 # create a counter to keep track of how many parent compounds of each kind we end up retaining
  
  for (l in 1:nrow(case_1_molec_species)) { # cycle through each of the case 1/4 putative compounds in this sample  
    
    groups_this_molec_species <- MAVEN_groups_this_sample[which(MAVEN_groups_this_sample$DB_lipid_class==case_1_molec_species$DB_lipid_class[l] & MAVEN_groups_this_sample$DB_degree_oxidation==case_1_molec_species$DB_degree_oxidation[l] & MAVEN_groups_this_sample$DB_FA_total_no_C==case_1_molec_species$DB_FA_total_no_C[l] & MAVEN_groups_this_sample$DB_FA_total_no_DB==case_1_molec_species$DB_FA_total_no_DB[l]),] # return group data corresponding to assignments in this sample that involve this putative parent compound (should be only one row in this case)
    
    peaks_this_molec_species <- MAVEN_peaks_this_sample[MAVEN_peaks_this_sample$groupID %in% groups_this_molec_species$groupID,] # return peak data for all assignments in this sample that involve this parent compound (should be only one row in this case)
            
    if (nrow(groups_this_molec_species)!=1) { # check to make sure there's only one assignment for this parent compound in this particular sample (otherwise something's wrong with the code)
      stop("More than one assignment for a case 1 or 4 compound. Something is wrong with the code. Not sure how it snuck by!") # stop script if this is the case
    }
    
    adduct_this_assignment <- groups_this_molec_species$DB_adduct # pull out the adduct type for which this assignment was made
    
    possible_adducts <- lipid_DB[lipid_DB$lipid_class==as.character(case_1_molec_species$DB_lipid_class[l]) & lipid_DB$degree_oxidation==as.character(case_1_molec_species$DB_degree_oxidation[l]) & lipid_DB$FA_total_no_C==as.character(case_1_molec_species$DB_FA_total_no_C[l]) & lipid_DB$FA_total_no_DB==as.character(case_1_molec_species$DB_FA_total_no_DB[l]),c("adduct","adduct_rank")] # for comparison purposes, get a list of all possible adducts for this parent compound from the DB
        
    possible_adducts <- possible_adducts[order(possible_adducts$adduct_rank, decreasing = TRUE, na.last = FALSE),] # reorder the list of possible adducts; must return any NA's (i.e., any adducts for which ranking was not given in DB) first so that rest of code works ok
  
    # case 1: the lone adduct of this parent compound in this particular sample is the theoretically most abundant adduct for the parent compound; hypothesized hierarchy is satisfied
    # and
    # case 4: the lone adduct is some adduct of lesser abundance, and MAVEN did not find a peak matching the adduct of this parent compound we expect to be most abundant; hierarchy is NOT satisfied
    
    ################# case 6 check #############
    
    # before continuing, check again to see whether we have a case 6 situation, then proceed accordingly
    
    if (!all(is.na(possible_adducts$adduct_rank))) { # either we have a case 6a situation, or no case 6 at all; either way, proceed with evaluation of case 1/4 compounds
      
      if (any(is.na(possible_adducts$adduct_rank))) { # we have a case 6a situation
        
        current_casecodes["C6a"] <- 1 # note that case 6a is true
        
      }
      
    if (as.character(possible_adducts$adduct[nrow(possible_adducts)])==as.character(adduct_this_assignment) & !is.na(possible_adducts$adduct[nrow(possible_adducts)])) { # parent compound is case 1
      
      current_casecodes["C1"] <- 1 # note that case 1 is true
      
    } else if (as.character(possible_adducts$adduct[nrow(possible_adducts)])!=as.character(adduct_this_assignment)) { # parent compound is case 4
      
      current_casecodes["C4"] <- 1 # note that case 4 is true
      
    }
    
    # last order of business before moving onto the next case 1/4 parent compound in this sample: select and then write to the peakdata results array the data for the current parent compound
    # method: we will use the rules (pretty simple for case 1 IDs) to determine whether data gets written or not (logic as agreed on by HFF and JRC in discussion on 1/13/2015)
    # if the ID is case 1 --> record data for this adduct
    # if the ID is case 4, assume misidentification since adduct which should be most abundant is not present --> write no data
    # if the species is case 6b and only one adduct was identified, record data for that adduct 
    
    if (current_casecodes["C1"]==1) { # case 1 --> record data
            
      peakdata_to_record <- peaks_this_molec_species[peaks_this_molec_species$groupID==groups_this_molec_species[as.character(possible_adducts$adduct[nrow(possible_adducts)])==as.character(groups_this_molec_species$DB_adduct),c("groupID")],] # select peak data to record
      
    } else if (current_casecodes["C4"]==1) { # case 4 --> record no peak data
      
      # do nothing!
      # should return NULL to console
      
    }
        
    } else { # we have a case 1-case 6b scenario: MAVEN ID'd a single peak representing only one good adduct for this parent compound, but we don't know how it ranks since there's no ranking data in the DB; we will call this case C1x
      
      current_casecodes["C1x"] <- 1  # note that case 1x is true; we can't really say whether it's a case 1 or 4 since we don't know
      
      current_casecodes["C6b"] <- 1  # note that case 6b is true
      
      peakdata_to_record <- peaks_this_molec_species[peaks_this_molec_species$groupID==groups_this_molec_species[as.character(adduct_this_assignment)==as.character(groups_this_molec_species$DB_adduct),c("groupID")],] # select peak data to record
      
    }
  
  # now, insert/append the peakdata and case codes for this parent compound into our peakdata results array
  
  if (exists("peakdata_to_record")) { # only force insertion if there's data
    
    retcount = retcount + nrow(peakdata_to_record) # increment our counter
    
    if (all(is.na(screened_peakdata[1,1:36]))) { # it's the first entry
      
      screened_peakdata <- as.data.frame(c(peakdata_to_record[,1:37],current_casecodes))
      
    } else { # it's not the first sample, so rbind
      
      screened_peakdata <- rbind(screened_peakdata,as.data.frame(c(peakdata_to_record[,1:37],current_casecodes)))
      
    }
    
    # finally, reset/remove our casecode and peakdata_to_record vectors before proceeding
    
    rm(peakdata_to_record)
    
  }
  
  current_casecodes[1:numcodes] <- 0
  
}

print(paste0("       Retained ",retcount,"."))

}

}

print(paste0("Done applying hierarchy rules. Groups representing parent compounds for which the adduct of greatest theoretical abundance was not present (i.e., C4 and C5 IDs) have been eliminated."))

# before proceeding, calculate how many peaks still remain in the dataset for each group at this point, and record

num_peaks_per_group.now <- as.data.frame(table(screened_peakdata$groupID)) # count number of peaks per group, right now

MAVEN_goodgroups[which(MAVEN_goodgroups$groupID %in% as.numeric(as.character(num_peaks_per_group.now[,1]))),c("peaksthisgroup_rules")] <- num_peaks_per_group.now[,2] # record this information to MAVEN_goodgroups

################# Final cleanup of screened_peakdata before flowing into results tables, below #############

# get rid of any groups for which we now have less than some number of peaks (min_peakspergroup_final, as specified by user at outset; doesn't have to be the same as min_peakspergroup used at the outset of analysis)

# this should very substantially reduce the number of groups still in the dataset

goodgroups_after_screening <- MAVEN_goodgroups[MAVEN_goodgroups$peaksthisgroup_rules>min_peakspergroup_final,]

################# Collate results of hierarchy rules application into a single matrix organized by group and sample #############

print(paste0("Now extracting appropriate peak data for all groups (assignments) still in dataset. In this step, the script will calculate mean peak quality, retention time, and m/z for each group. In addition, the script will screen check to see whether we might have missed any case 3r assignments where the same adduct was observed at different retention times but in different samples..."))

# also, add a column of overall case codes for each group to goodgroups_after_screening, and columns for casecode summary information

goodgroups_after_screening$casecodes_thisgroup.summary <- NA # create a new column for the code summary in goodgroups_after_screening
goodgroups_after_screening[,c("casecodes_thisgroup.C1","casecodes_thisgroup.C1x","casecodes_thisgroup.C2a","casecodes_thisgroup.C2b","casecodes_thisgroup.C3c","casecodes_thisgroup.C3f","casecodes_thisgroup.C3r","casecodes_thisgroup.C4","casecodes_thisgroup.C5","casecodes_thisgroup.C6a","casecodes_thisgroup.C6b","casecodes_thisgroup.C3f.rt.and.ignore.oddC","casecodes_thisgroup.C3c.rt.and.ignore.oddC")] <- NA # create new columns for the group code information in goodgroups_after_screening

# "casecodes_thisgroup.C3f.rt.and.ignore.oddC" and "casecodes_thisgroup.C3c.rt.and.ignore.oddC" annotate C3f and C3c assignments while excluding any potential competitors with an # of odd fatty acid carbon atoms, an additional rule invoked in Collins et al., 2015, for data of biological origin

species_info <- c("lipid_class","FA_total_no_DB","FA_total_no_C","degree_oxidation","DB_adduct","DB_adduct_mz","DB_parent_exactneutralmass") # the basic chemical properties of each unique parent compound for which we'll be reporting results: IPL type, the number of total fatty acid chain carbon atoms ("FA_total_no_C"), number of total fatty acid chain double bonds ("FA_total_no_DB"), number of oxygen atoms in excess of those which are normally part of the headgroup structure ("degree_oxidation"), the adduct of the compound for which the peak information in each sample is reported, the m/z of that particular adduct, and the exact mass of the unionized parent compound

addl_group_stats <- c("mean_qual","stddev_qual","mean_rt","stddev_rt","mean_groupmz","stddev_groupmz","mean_match_delta_ppm","peaksthisgroup") # additional statistics we'll calculate and record for the peakdata that remain in each group 

# create array screenedpeaks_bygroup.full into which we'll flow results

screenedpeaks_bygroup.full <- matrix(NA,nrow(goodgroups_after_screening),2+length(species_info)+length(addl_group_stats)+length(samples)*5) # create summary array
screenedpeaks_bygroup.full <- as.data.frame(screenedpeaks_bygroup.full) # convert to data frame
colnames(screenedpeaks_bygroup.full) <- c("MAVEN_groupID","VML_DB_fragID",species_info,addl_group_stats,apply(expand.grid(samples,c("quality","rt","peakAreaCorrected","observed_mz","casecodes_sample"))[,c(1,2)],1,paste,collapse=".")) # create and apply column names
screenedpeaks_bygroup.full[1,c(1,2,4,5,8,9,10:(10+length(samples)*4+length(addl_group_stats)-1))] <- 0 # force data type num on certain columns
screenedpeaks_bygroup.full[1,10:(10+length(samples)*4+length(addl_group_stats)-1)] <- NA # reset values to NA
screenedpeaks_bygroup.full$casecodes_thisgroup.summary <- NA # create a new column for the code summary in screenedpeaks_bygroup.full

screenedpeaks_bygroup.full[,c("casecodes_thisgroup.C1","casecodes_thisgroup.C1x","casecodes_thisgroup.C2a","casecodes_thisgroup.C2b","casecodes_thisgroup.C3c","casecodes_thisgroup.C3f","casecodes_thisgroup.C3r","casecodes_thisgroup.C4","casecodes_thisgroup.C5","casecodes_thisgroup.C6a","casecodes_thisgroup.C6b","casecodes_thisgroup.C3f.rt.and.ignore.oddC","casecodes_thisgroup.C3c.rt.and.ignore.oddC")] <- 0 # create new columns for the group code information in screenedpeaks_bygroup.full

# "casecodes_thisgroup.C3f.rt.and.ignore.oddC" and "casecodes_thisgroup.C3c.rt.and.ignore.oddC" annotate C3f and C3c assignments while excluding any potential competitors with an # of odd fatty acid carbon atoms, an additional rule invoked in Collins et al., 2015, for data of biological origin

# cycle through each row in goodgroups_after_screening and flow data into screenedpeaks_bygroup.full as appropriate

for (n in 1:nrow(goodgroups_after_screening)) {
  
  group.casecodes <- vector(mode="double",length=numcodes) # create or reset the vector to keep track of codes for this group
  
  # convert some numerics back to plain language by way of original factor levels
  thislipid_class <- lipid_class_levels[goodgroups_after_screening$DB_lipid_class[n]]
  thisoxy <- degree_oxidation_levels[goodgroups_after_screening$DB_degree_oxidation[n]]
  thisadduct <- adduct_levels[goodgroups_after_screening$DB_adduct[n]]
  
  # insert basic chemical information for this parent compound (elements of species_info)
  screenedpeaks_bygroup.full$lipid_class[n] <- thislipid_class
  screenedpeaks_bygroup.full$FA_total_no_C[n]  <- goodgroups_after_screening$DB_FA_total_no_C[n]
  screenedpeaks_bygroup.full$FA_total_no_DB[n] <- goodgroups_after_screening$DB_FA_total_no_DB[n]
  screenedpeaks_bygroup.full$degree_oxidation[n] <- thisoxy
  screenedpeaks_bygroup.full$DB_adduct[n] <- thisadduct
  screenedpeaks_bygroup.full$DB_adduct_mz[n] <- goodgroups_after_screening$DB_mz[n]
  screenedpeaks_bygroup.full$DB_parent_exactneutralmass[n] <- goodgroups_after_screening$DB_neutral_mass[n]
  
  # insert MAVEN_groupID and VML_DB_fragID for this parent compound
  screenedpeaks_bygroup.full$MAVEN_groupID[n]  <- goodgroups_after_screening$groupID[n]
  screenedpeaks_bygroup.full$VML_DB_fragID[n]  <- goodgroups_after_screening$frag_ID[n]
      
  for (o in 1:length(samples)) { # now, cycle through sample IDs, extract peak data for this compound in that sample (if it exists), and insert it in the right place
    
    peakdata_this_sample <- screened_peakdata[screened_peakdata$groupID==goodgroups_after_screening$groupID[n] & as.character(screened_peakdata$sample)==as.character(samples[o]),] # extract matching peak data
                          
    if (nrow(peakdata_this_sample)==1) { # there is data for this parent compound in this sample; let's insert it
      
      codes_this_feature <- casecodes[peakdata_this_sample[38:(38+numcodes-1)]==1] # pull case codes for this feature
      codestring_this_feature <- paste(codes_this_feature,collapse="; ") # build a single string of case codes
      group.casecodes <- group.casecodes+peakdata_this_sample[38:(38+numcodes-1)] # while we're at it, add these codes to the group code vector
      
      loc_end_spectdat <- (2+length(species_info)+length(addl_group_stats)) # location of last column of parent compound data in the output arrays
      
      insertcols_full <- c(loc_end_spectdat+o,loc_end_spectdat+length(samples)+o,loc_end_spectdat+length(samples)*2+o,loc_end_spectdat+length(samples)*3+o,loc_end_spectdat+length(samples)*4+o) # calculate the columns into which we should insert the gathered peakdata in screenedpeaks_bygroup.full
                  
      screenedpeaks_bygroup.full[n,insertcols_full[1:4]] <- peakdata_this_sample[c("quality","rt","peakAreaCorrected","peakMz")] # insert peak data
      screenedpeaks_bygroup.full[n,insertcols_full[5]] <- codestring_this_feature # insert case codes
      
      
    }
    
  }
    
  # now, calculate summary LCMS statistics for this group and insert into screenedpeaks_bygroup.full
  
  # calculate values
  
  thisgroup_mean_qual <- mean(as.numeric(screenedpeaks_bygroup.full[n,(loc_end_spectdat+1):(loc_end_spectdat+length(samples))]), na.rm = T)
  thisgroup_stddev_qual <- sd(as.numeric(screenedpeaks_bygroup.full[n,(loc_end_spectdat+1):(loc_end_spectdat+length(samples))]), na.rm = T)
  thisgroup_mean_rt <- mean(as.numeric(screenedpeaks_bygroup.full[n,(loc_end_spectdat+1+length(samples)):(loc_end_spectdat+2*length(samples))]), na.rm = T)
  thisgroup_stddev_rt <- sd(as.numeric(screenedpeaks_bygroup.full[n,(loc_end_spectdat+1+length(samples)):(loc_end_spectdat+2*length(samples))]), na.rm = T)
  thisgroup_mean_groupmz <- mean(as.numeric(screenedpeaks_bygroup.full[n,(loc_end_spectdat+1+3*length(samples)):(loc_end_spectdat+4*length(samples))]), na.rm = T)
  thisgroup_stddev_groupmz <- sd(as.numeric(screenedpeaks_bygroup.full[n,(loc_end_spectdat+1+3*length(samples)):(loc_end_spectdat+4*length(samples))]), na.rm = T)
     
  peaksthisgroup <- goodgroups_after_screening$peaksthisgroup_rules[n]
  
  # insert
  
  screenedpeaks_bygroup.full[n,c("mean_qual","stddev_qual","mean_rt","stddev_rt","mean_groupmz","stddev_groupmz","peaksthisgroup")] <- c(thisgroup_mean_qual,thisgroup_stddev_qual,thisgroup_mean_rt,thisgroup_stddev_rt,thisgroup_mean_groupmz,thisgroup_stddev_groupmz,peaksthisgroup)  
  
  # check to see whether we missed C3r designation for this group
  
  # inspection of data from earlier versions of the script showed that some case C3r's may sneak through to this point unidentified; because we examined data by sample, we wouldn't catch any instances where the same adduct appeared at different r/ts in *different* samples
    
  # can use presence of duplicate fragment IDs (from the database) to diagnose this situation
  if (length(which(goodgroups_after_screening$frag_ID[n]==goodgroups_after_screening$frag_ID))>1) { # this fragID appears more than once in the groups we still have    
    
    group.casecodes[c("C3r")] <- group.casecodes[c("C3r")] + 1 # update the case codes for this current group accordingly 
    
  }
  
  # insert the case code summary information for this group into goodgroups_after_screening and screenedpeaks_bygroup.full
  
  codes_this_group <- casecodes[group.casecodes>=1] # pull case codes for this feature
  codestring_this_group <- paste(codes_this_group,collapse="; ") # build a single string of case codes to insert into casecodes_thisgroup.summary
  
  # insert information in casecodes_thisgroup.summary
  goodgroups_after_screening$casecodes_thisgroup.summary[n] <- codestring_this_group 
  screenedpeaks_bygroup.full$casecodes_thisgroup.summary[n] <- codestring_this_group
  
  # insert case code tabulation for this group into appropriate fields in each array
  
  goodgroups_after_screening[n,c("casecodes_thisgroup.C1","casecodes_thisgroup.C1x","casecodes_thisgroup.C2a","casecodes_thisgroup.C2b","casecodes_thisgroup.C3c","casecodes_thisgroup.C3f","casecodes_thisgroup.C3r","casecodes_thisgroup.C4","casecodes_thisgroup.C5","casecodes_thisgroup.C6a","casecodes_thisgroup.C6b")] <- as.numeric(group.casecodes)
  screenedpeaks_bygroup.full[n,c("casecodes_thisgroup.C1","casecodes_thisgroup.C1x","casecodes_thisgroup.C2a","casecodes_thisgroup.C2b","casecodes_thisgroup.C3c","casecodes_thisgroup.C3f","casecodes_thisgroup.C3r","casecodes_thisgroup.C4","casecodes_thisgroup.C5","casecodes_thisgroup.C6a","casecodes_thisgroup.C6b")] <- as.numeric(group.casecodes)
  
}

print(paste0("Extraction and matching of peak data complete. Output results table screenedpeaks_bygroup.full has been created."))

print(paste0("Now calculating delta ppm between the mean m/z of the observed features in each group and the exact m/z of the assigned adduct as it appears in the database..."))

# finally, calculate & store delta ppm for each group

screenedpeaks_bygroup.full[,c("mean_match_delta_ppm")] <- ((screenedpeaks_bygroup.full[,c("DB_adduct_mz")]-screenedpeaks_bygroup.full[,c("mean_groupmz")])/screenedpeaks_bygroup.full[,c("DB_adduct_mz")])*1000000

################# Further annotation using crude retention time window rules #############

# will annotate each group in screenedpeaks_bygroup.full to indicate whether it satisfies the retention time window criteria for IPL species (from data in rt_window_file); this information (as 0 or 1) will get stored in an additional field "rt_window.passfail" 

rt_windata <- read.csv(rt_window_file, skip = 0, header = T, stringsAsFactors=FALSE) # read in rt window data

# preallocate the destination field in each of the two output arrays with a "fail" value of 0

screenedpeaks_bygroup.full$rt_window.passfail <- 0

# cycle through each class of lipid in the retention time window data and apply it to the screened data in screenedpeaks_bygroup.full and screenedpeaks_bygroup.sum
# using the less conservative rt_win_min for lower bound right now, but could easily change here to use rt2_win_min, which would accomodate some additional potential IDs in each class that are highly unsaturated with low number of C

print(paste0("Annotating data for compliance with crude retention time windows established for each IPL type..."))

for (s in 1:nrow(rt_windata)) {
  
  if (!is.na(rt_windata$rt_win_min[s]) & !is.na(rt_windata$rt_win_max[s])) { # we have an upper and lower bound for this IPL
    
    screenedpeaks_bygroup.full[(screenedpeaks_bygroup.full$lipid_class==rt_windata$lipid_class[s]) & ((as.numeric(screenedpeaks_bygroup.full$mean_rt)+as.numeric(screenedpeaks_bygroup.full$stddev_rt))>=rt_windata$rt_win_min[s]) & ((as.numeric(screenedpeaks_bygroup.full$mean_rt)-as.numeric(screenedpeaks_bygroup.full$stddev_rt))<=rt_windata$rt_win_max[s]),c("rt_window.passfail")] <- 1
    
  } else if (!is.na(rt_windata$rt_win_min[s]) & is.na(rt_windata$rt_win_max[s])) { # we only have a lower bound
    
    screenedpeaks_bygroup.full[(screenedpeaks_bygroup.full$lipid_class==rt_windata$lipid_class[s]) & ((as.numeric(screenedpeaks_bygroup.full$mean_rt)+as.numeric(screenedpeaks_bygroup.full$stddev_rt))>=rt_windata$rt_win_min[s]),c("rt_window.passfail")] <- 1
    
  } else if (is.na(rt_windata$rt_win_min[s]) & !is.na(rt_windata$rt_win_max[s])) { # we only have an upper bound
    
    screenedpeaks_bygroup.full[(screenedpeaks_bygroup.full$lipid_class==rt_windata$lipid_class[s]) & ((as.numeric(screenedpeaks_bygroup.full$mean_rt)-as.numeric(screenedpeaks_bygroup.full$stddev_rt))<=rt_windata$rt_win_max[s]),c("rt_window.passfail")] <- 1
    
  } else if (is.na(rt_windata$rt_win_min[s]) & is.na(rt_windata$rt_win_max[s])) { # we have no rt time data for this IPL type --> record NA
    
    screenedpeaks_bygroup.full[(screenedpeaks_bygroup.full$lipid_class==rt_windata$lipid_class[s]),c("rt_window.passfail")] <- NA
    
  }
  
}

print(paste0("Retention time annotation complete."))

################# Finally, resolution of case 3c and 3f compounds: check to see whether groups (i.e., putative IDs) still in the dataset at this point overlap or compete with any others #############

# Case 3c and/or 3f annotation is indicated when MAVEN has made "overlapping" or "competing" assignments, i.e., multiple groups with different database assignments have been created for the same set of peaks. There are two cases: adducts of two or more different parent compounds share the exact same m/z (case 3f, functional structural isomers; easy to identify), or adducts of two or more different parent compounds are so close in their m/z (i.e., not exactly the same, but within the ppm matching window specified by the user in MAVEN) that MAVEN could not discriminate between them and therefore created multiple groups for a single set of peaks (case 3c; more difficult to diagnose)

print(paste0("Finally, perform screening to identify competing/overlapping (i.e., case 3c or 3f) assignments. Note that unlike earlier annotation of compounds using the hierarchy rules, where we ID'd cases by sample, we will perform the search for case 3c/3f compounds by group (i.e., by assignment). Thus, do not expect to find any case 3c/3f annotation codes on individual sample data -- the case 3c/3f designation will be assigned to groups as a whole."))

print(paste0("The script will also perform an additional, more restrictive case 3c/3f screening while ignoring assignments whose parent compound has an odd number of fatty acid carbon atoms, and those which have failed the RT window screening critera. Results of this screening will be dumped to casecodes_thisgroup.C3f.rt.and.ignore.oddC and casecodes_thisgroup.C3c.rt.and.ignore.oddC"))

# first, create an alternate dataset screenedpeaks_bygroup.full.evenC.rtpass that includes only assignments with and even # of FA C atoms, and which have passed the RT rules (need this for the additional, more restrictive case 3 screening)

screenedpeaks_bygroup.full.evenC.rtpass <- screenedpeaks_bygroup.full[screenedpeaks_bygroup.full$rt_window.passfail==1 & screenedpeaks_bygroup.full$FA_total_no_C%%2==0,]

print(paste0("   Now searching for any case 3c/3f molecular species..."))

retcount.3f = 0 # create some counters to keep track of how many species of each kind we end up retaining
retcount.3f.rt.and.ignore.oddC = 0
retcount.3c = 0
retcount.3c.rt.and.ignore.oddC = 0

for (t in 1:nrow(screenedpeaks_bygroup.full)) { # now that groups have peak data attached to them, cycle through each and check to see whether there's any overlap 
  
  # case 3f, first
  
  # including adducts that failed RT rules and those with odd # of FA C
  
  if (length(which(screenedpeaks_bygroup.full$DB_adduct_mz[t]==screenedpeaks_bygroup.full$DB_adduct_mz))>1) { # this exact m/z appears more than once in the groups we still have --> case 3f
    
    screenedpeaks_bygroup.full$casecodes_thisgroup.C3f[t] <- 1 # update casecodes_thisgroup.C3f for the current group accordingly
    screenedpeaks_bygroup.full$casecodes_thisgroup.summary[t] <- paste0(screenedpeaks_bygroup.full$casecodes_thisgroup.summary[t],"; C3f") # append "C3f" to casecodes_thisgroup.summary for the current group
    
    # ***** future plan: record the groupIDs of the assignments that compete with this one in another field
    
    retcount.3f <- retcount.3f + 1 # update our counter
    
  }
  
  # excluding adducts that failed RT rules and those with odd # of FA C
  
  if (length(which(screenedpeaks_bygroup.full$DB_adduct_mz[t]==screenedpeaks_bygroup.full.evenC.rtpass$DB_adduct_mz))>1) { # this exact m/z appears more than once in the groups we still have --> case 3f  
    
    screenedpeaks_bygroup.full$casecodes_thisgroup.C3f.rt.and.ignore.oddC[t] <- 1 # update casecodes_thisgroup.C3f.rt.and.ignore.oddC for the current group accordingly
    
    # ***** future plan: record the groupIDs of the assignments that compete with this one in another field
    
    retcount.3f.rt.and.ignore.oddC <- retcount.3f.rt.and.ignore.oddC + 1 # update our counter
    
  }
  
  # now, case 3c
  
  # screen all other groups against this one using the ppm tolerance that user applied initially in MAVEN; other groups with a DB_adduct_mz which are within the DB_adduct_mz of this group ± the ppm tolerance may represent potential competing assignments
  
  potential_case3cs.thisgroup <- screenedpeaks_bygroup.full[which((abs(screenedpeaks_bygroup.full$DB_adduct_mz-screenedpeaks_bygroup.full$DB_adduct_mz[t])<=(MAVEN_ppm_tolerance/1000)) & (abs(screenedpeaks_bygroup.full$DB_adduct_mz-screenedpeaks_bygroup.full$DB_adduct_mz[t])>0)),] # pull out all other groups that fall within DB_adduct_mz ± ppm tolerance of this group, but don't share the *exact* same DB_adduct_mz (since we've presumably already captured those as case 3fs, above)
  
  potential_case3cs.thisgroup.rt.and.ignore.oddC <- screenedpeaks_bygroup.full[which((abs(screenedpeaks_bygroup.full$DB_adduct_mz-screenedpeaks_bygroup.full$DB_adduct_mz[t])<=(MAVEN_ppm_tolerance/1000)) & (abs(screenedpeaks_bygroup.full$DB_adduct_mz-screenedpeaks_bygroup.full$DB_adduct_mz[t])>0) & (screenedpeaks_bygroup.full$FA_total_no_C%%2==0) & (screenedpeaks_bygroup.full$rt_window.passfail==1)),] # do the same as above, but with the additional RT and C # restrictions
  
  # including adducts that failed RT rules and those with odd # of FA C
  
  if (nrow(potential_case3cs.thisgroup)>0) { # it's possible we have a case 3c... but it's not enough that there's another assignment that falls within DB_adduct_mz ± ppm tolerance of this one --> we have to observe actual overlap of some or all actual peak data
    # we will diagnose this overlap of peak data by testing for commonality in the peakAreaCorrected data, assuming it's extremely improbable that any of the features attached to the group we're considering would share by way of pure chance any of the exact same peakAreaCorrected as a feature within the subset we've just drawn out
    
    peakareadata.start <- 2+length(species_info)+length(addl_group_stats)+length(samples)*2+1 # column at which individual peakAreaCorrected entries start in screenedpeaks_bygroup.full
    peakareadata.end <- 2+length(species_info)+length(addl_group_stats)+length(samples)*3 # column at which individual peakAreaCorrected entries end in screenedpeaks_bygroup.full
    peakAreaCorrected.data.from.potential_case3cs <- c(t(potential_case3cs.thisgroup[,peakareadata.start:peakareadata.end])) # pull out peakAreaCorrected data for all the potential case 3cs that we suspect overlap with this group
    
    if (length(which(screenedpeaks_bygroup.full[t,peakareadata.start:peakareadata.end] %in% peakAreaCorrected.data.from.potential_case3cs))>0) { # this is a bona fide case 3c situation -- we have some overlap of peakAreaCorrected data --> record as case 3c
      
      screenedpeaks_bygroup.full$casecodes_thisgroup.C3c[t] <- 1 # update casecodes_thisgroup.C3c for the current group accordingly
      screenedpeaks_bygroup.full$casecodes_thisgroup.summary[t] <- paste0(screenedpeaks_bygroup.full$casecodes_thisgroup.summary[t],"; C3c") # append "C3c" to casecodes_thisgroup.summary for the current group
      
      # ***** future plan: record the groupIDs of the assignments that compete with this one in another field
      
      retcount.3c <- retcount.3c + 1 # update our counter
      
    }
    
  }
  
  # excluding adducts that failed RT rules and those with odd # of FA C
  
  if (nrow(potential_case3cs.thisgroup.rt.and.ignore.oddC)>0) { # it's possible we have a case 3c... but it's not enough that there's another assignment that falls within DB_adduct_mz ± ppm tolerance of this one --> we have to observe actual overlap of some or all actual peak data
    # we will diagnose this overlap of peak data by testing for commonality in the peakAreaCorrected data, assuming it's extremely improbable that any of the features attached to the group we're considering would share by way of pure chance any of the exact same peakAreaCorrected as a feature within the subset we've just drawn out
    
    peakareadata.start <- 2+length(species_info)+length(addl_group_stats)+length(samples)*2+1 # column at which individual peakAreaCorrected entries start in screenedpeaks_bygroup.full
    peakareadata.end <- 2+length(species_info)+length(addl_group_stats)+length(samples)*3 # column at which individual peakAreaCorrected entries end in screenedpeaks_bygroup.full
    peakAreaCorrected.data.from.potential_case3cs <- c(t(potential_case3cs.thisgroup.rt.and.ignore.oddC[,peakareadata.start:peakareadata.end])) # pull out peakAreaCorrected data for all the potential case 3bs that we suspect overlap with this group
    
    if (length(which(screenedpeaks_bygroup.full.evenC.rtpass[screenedpeaks_bygroup.full$MAVEN_groupID[t]==screenedpeaks_bygroup.full.evenC.rtpass$MAVEN_groupID,peakareadata.start:peakareadata.end] %in% potential_case3cs.thisgroup.rt.and.ignore.oddC))>0) { # this is a bona fide case 3c situation -- we have some overlap of peakAreaCorrected data --> record as case 3c
      
      screenedpeaks_bygroup.full$casecodes_thisgroup.C3c.rt.and.ignore.oddC[t] <- 1 # update casecodes_thisgroup.C3c.rt.and.ignore.oddC for the current group accordingly
      
      # ***** future plan: record the groupIDs of the assignments that compete with this one in another field
      
      retcount.3c.rt.and.ignore.oddC <- retcount.3c.rt.and.ignore.oddC + 1 # update our counter
      
    }
    
  }
  
}

print(paste0("       Found ",retcount.3f," case 3f species while including potential assignments with an odd number of fatty acid carbon atoms and those which failed the retention time window criteria; ",retcount.3f.rt.and.ignore.oddC," while ignoring."))
print(paste0("       Found ",retcount.3c," case 3c species while including potential assignments with an odd number of fatty acid carbon atoms and those which failed the retention time window criteria; ",retcount.3c.rt.and.ignore.oddC," while ignoring."))

################# Set correct data type on fields in screenedpeaks_bygroup.full #############

# as of 5/15/15, appears that some fields which contain numeric data are set to chr at this point; resetting manually, though it would of course be good to (eventually) sleuth out the origin of the issue

screenedpeaks_bygroup.full[,-c(3,6,7,(loc_end_spectdat+length(samples)*4+1):(loc_end_spectdat+length(samples)*5+1))] <- lapply(screenedpeaks_bygroup.full[,-c(3,6,7,(loc_end_spectdat+length(samples)*4+1):(loc_end_spectdat+length(samples)*5+1))],as.numeric)

################# File export, if desired #############

# sort the output table by m/z

screenedpeaks_bygroup.full <- screenedpeaks_bygroup.full[order(screenedpeaks_bygroup.full$DB_adduct_mz),]

# write data to .csv file, making sure we add a unique timestamp to the filename so we don't overwrite any existing output

output_DTG <- format(Sys.time(), "%Y-%m-%dT%X%z") # return current time in a good format
output_DTG <- gsub(" ", "_", output_DTG) # replace any spaces
output_DTG <- gsub(":", "-", output_DTG) # replaces any colons with dashes (Mac compatibility)

print(paste0("Exporting output .csv file..."))

write.csv(screenedpeaks_bygroup.full, paste("screenedpeaks_bygroup.full_",output_DTG,".csv",sep=""))

print(paste0("File saved as screenedpeaks_bygroup.full_",output_DTG))