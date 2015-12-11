################ Set classes, methods #############

# create a class for our screened xsAnnotate object to which compound assignments and confidence codes have been applied

setClass("xsALOB", contains="xsAnnotate", representation(confcode = "factor"), prototype = prototype(confcode = factor(character(0)))) # maybe don't want to use "contains" ? can decide later

################ Wrapper function #############

# makeLOBAssignments: Wrapper function for screening & annotation of an xsAnnotate object

makeLOBAssignments = function(xsA, polarity = NULL, database = NULL, remove.iso = TRUE, rt.restrict =  TRUE, rt.windows = NULL, exclude.oddFA = TRUE, match.ppm = 2.5) { # add nSlaves option at some point?
  
################ Check arguments, load necessary data #############
  
  # verify that input xsA is of correct class, has pseudospectra, and that rt data are in seconds
  
  if (!class(xsA)=="xsAnnotate") {
    
    stop("Input 'xsA' is not an 'xsAnnotate' object.\n")
    
  }
  
  if (length(xsA@pspectra)<1) {
    
    stop("Input 'xsA' does not appear to have any CAMERA pseudospectra. Re-run CAMERA or examine the input object!\n")
    
  }
  
  if (sum(xsA@groupInfo[colnames(xsA@groupInfo)=="rt"]>100)==0) { # likely that rt data are in minutes; they need to be in seconds
    
    stop("Retention time data in input object 'xsA' appear to be in minutes. Please ensure data are recorded in seconds during conversion of manufacturer data files to mzXML.\n")
    
  }
  
  # check for consistent polarity between xsA and value of argument 'polarity' provided by user
  
  if (is.null(polarity)) { # user didn't provide value for argument, try to automatically detect polarity from xsA@polarity
    
    cat("User did not specify a value for argument 'polarity.' Attempting to detect current polarity from property of input 'xsA'...\n")
    
    if (!is.null(xsA@polarity)) {
      
      if (is.element(xsA@polarity,c("positive","negative"))) {
        
        polarity = xsA@polarity
        cat("Input 'xsA' appears to be of polarity '",polarity,"'\n")
        
      } else {
        
        stop("Polarity of input 'xsA' is not 'positive' or 'negative.'\nIf value for argument 'polarity' is not specified by user in call to makeLOBAssignments, input object must have a valid polarity.\n")
        
      }
      
    } else {
      
      stop("Input 'xsA' has no assigned polarity. Re-create xsAnnotate object with CAMERA, or specify value for argument 'polarity' in call to makeLOBAssignments.\n")
      
    }
    
  } else { # user provided some value for 'polarity'
    
    if (!is.element(polarity,c("positive","negative"))) {
      
      stop("Value for argument 'polarity' must be 'positive' or 'negative.'\n")
      
    } else {
      
      if (!is.null(xsA@polarity)) {
        
        if (xsA@polarity!=polarity) {
          
          stop("Value '",polarity,"' given for argument 'polarity' does not match apparent polarity '",xsA@polarity,"' of input 'xsA.'\n")
          
        }
        
      }
      
    }
    
  }
  
  # check for database input, use correct polarity default DB if no input
  
  if (!is.null(database)) {
    
    if (!class(database)=="LOBdbase") {
      
      stop("Input 'database' is not a 'LOBdbase' object. Use loadLOBdbase() to read a user-supplied lipid-ox-lipid-oxylipin database from a .csv file, or use generateLOBdbase() to generate a new LOBSTAHS database.\n")
      
    } else { # make sure it's the correct polarity
      
      if (polarity!=database@polarity) {
        
        stop("Polarity '",database@polarity,"' of database does not match polarity '",polarity,"' of input 'xsA' object. Specify a 'LOBdbase' object of the appropriate polarity.\n")
        
      }
      
    }
    
  } else { # user didn't specify a database, use the default DB of the correct polarity
    
    cat("User did not specify an external database. Loading default LOBSTAHS database for polarity '",polarity,"' ...\n")
    
    load("dependencies/default.LOBdbase.RData")
    database = LOBdbase[[polarity]]
    
  }
  
  # set things up for rt window screening
  
  if (rt.restrict = TRUE) {
    
    if (rt.windows = NULL) { # use defaults
      
      load("dependencies/default.rt.windows.RData")
      
      warning("User did not specify an external source of retention time window data. Loading default rt windows used for the Van Mooy Lab Orbitrap. Non-VML users are strongly urged to use argument 'rt.wintable' to load retention time data specific to their chromatography...\n")
      
    } else { # load and perform basic check of user-specified rt window data, assuming they pointed rt.wintable to a valid R matrix with reasonbly named column headings
      
      if (sum(grepl("[Ll][Ii][Pp][Ii][Dd][ |\\.|_]*[Cc][Ll][Aa][Ss][Ss]",colnames(rt.windows)))!=1) {
        
        stop("Could not find unique field 'lipid_class' in user's retention time data.\n")
        
      } else {
        
        rt.windows$lipid_class = as.character(rt.windows[,grepl("[Ll][Ii][Pp][Ii][Dd][ |\\.|_]*[Cc][Ll][Aa][Ss][Ss]",colnames(rt.windows))])
        
      }
      
      if (sum(grepl("[Rr][Tt][ |\\.|_]*[Ww][Ii][Nn][ |\\.|_]*[Mm][Ii][Nn]",colnames(rt.windows)))!=1) {
        
        stop("Could not find unique field 'rt_win_min' in user's retention time data.\n")
        
      } else {
        
        rt.windows$rt_win_min = as.numeric(rt.windows[,grepl("[Rr][Tt][ |\\.|_]*[Ww][Ii][Nn][ |\\.|_]*[Mm][Ii][Nn]",colnames(rt.windows))])
        
      }
      
      if (sum(grepl("[Rr][Tt][ |\\.|_]*[Ww][Ii][Nn][ |\\.|_]*[Mm][Aa][Xx]",colnames(rt.windows)))!=1) {
        
        stop("Could not find unique field 'rt_win_max' in user's retention time data.\n")
        
      } else {
        
        rt.windows$rt_win_max = as.numeric(rt.windows[,grepl("[Rr][Tt][ |\\.|_]*[Ww][Ii][Nn][ |\\.|_]*[Mm][Aa][Xx]",colnames(rt.windows))])
        
      }
      
      if (any((rt.windows$rt_win_min|rt.windows$rt_win_max)>100, na.rm = TRUE)) { # likely that user's data are in seconds, not minutes
        
        warning("Retention time data should be specified in minutes, not seconds. User's data appear to be in seconds.")
        
      }
      
    }
    
  }
  
  # lastly, if user has elected remove.iso, make sure there are actually isotopes ID'd in the xsAnnotate object
  
  if ((remove.iso = TRUE) & (length(xsA@isotopes)<1)) {
    
    stop("Input 'xsA' does not appear to contain any identified isotopes. Re-run CAMERA and use findIsotopes!")
    
  }
  
}

################ Perform screening #############

pspectra = 1:length(xsA@pspectra)

screened.pspecdata = lapply(pspectra, screenPSpectrum, xsA, polarity = polarity, database = database, remove.iso = remove.iso, rt.restrict = rt.restrict, rt.windows = rt.windows, exclude.oddFA = exclude.oddFA, match.ppm = match.ppm)



################ Helper functions (some private) #############

# loadLOBdbase: loads and creates a LOBdbase object from a properly formatted .csv file
# presumes file has correct header names somewhere in the file, but order and spelling/capitalization shouldn't matter

loadLOBdbase = function(file, polarity, num.compounds = NULL) {
  
  polarity = match.arg(polarity, c("positive","negative"), several.ok = FALSE)
  
  db.readraw = read.table(file, sep = ",", skip = 0, header = FALSE, colClasses = "character") # read in raw data first to determine where headers are (in case there are any junk lines at top of file)
  
  ind.header = which(matrix(grepl("[Mm]/*[Zz]", db.readraw), ncol=ncol(db.readraw)), arr.ind=TRUE) # assuming position of m/z column header is a reasonable indicator of where data starts
  
  header.row = ind.header[1]
  
  if (length(header.row)>1) {
    
    stop("Found multiple possible header rows in .csv file. Please check the format of your database.")
    
  }
  
  db.rawdata = read.table(file, sep = ",", skip = header.row-1, header = TRUE, colClasses = "character") # read in file, for real
  
  object = new("LOBdbase") # create new LOBdbase object
  
  # determine which fields are where, being as forgiving as possible regarding capitalization and punctuation of column headers, then load data from that column into appropriate slot in the LOBdbase object
  
  object@mz = as.numeric(db.rawdata[,grepl("[Mm]/*[Zz]",colnames(db.rawdata))])
  
  # only optional field is frag_ID
  
  if (sum(grepl("[Ff][Rr][Aa][Gg][ |\\.|_]*[Ii][Dd]",colnames(db.rawdata)))==1) { # user's DB appears to have frag_IDs
    
    object@frag_ID = as.integer(db.rawdata[,grepl("[Ff][Rr][Aa][Gg][ |\\.|_]*[Ii][Dd]",colnames(db.rawdata))])
    
  } else {
    
    object@frag_ID = 1:length(object@mz)
    
    cat("Database being imported doesn't appear to have field 'frag_ID'; frag_ID will be assigned sequentially based on order of entries in .csv file.\n")
  
    }
  
  object@exact_parent_neutral_mass = as.numeric(db.rawdata[,rowSums(sapply(c("[Ee][Xx][Aa][Cc][Tt][ |\\.|_]*[Pp][Aa][Rr][Ee][Nn][Tt][ |\\.|_]*[Nn][Ee][Uu][Tt][Rr][Aa][Ll][ |\\.|_]*[Mm][Aa][Ss][Ss]","[Ee][Xx][Aa][Cc][Tt][ |\\.|_]*[Nn][Ee][Uu][Tt][Rr][Aa][Ll][ |\\.|_]*[Mm][Aa][Ss][Ss]","[Ee][Xx][Aa][Cc][Tt][ |\\.|_]*[Mm][Aa][Ss][Ss]","[Pp][Aa][Rr][Ee][Nn][Tt][ |\\.|_]*[Nn][Ee][Uu][Tt][Rr][Aa][Ll][ |\\.|_]*[Mm][Aa][Ss][Ss]"),grepl,colnames(db.rawdata),any))>0])
  
  object@lipid_class = as.factor(db.rawdata[,grepl("[Ll][Ii][Pp][Ii][Dd][ |\\.|_]*[Cc][Ll][Aa][Ss][Ss]",colnames(db.rawdata))])
  
  object@species = as.character(db.rawdata[,grepl("[Ss][Pp][Ee][Cc][Ii][Ee][Ss]",colnames(db.rawdata))])

  object@adduct = as.factor(db.rawdata[,grepl("^[Aa][Dd][Dd][Uu][Cc][Tt]$",colnames(db.rawdata))])
  
  object@adduct_rank = as.integer(db.rawdata[,grepl("[Aa][Dd][Dd][Uu][Cc][Tt][ |\\.|_]*[Rr][Aa][Nn][Kk]",colnames(db.rawdata))])
  
  object@FA_total_no_C = as.integer(db.rawdata[,grepl("[Ff][Aa][ |\\.|_]*[Tt][Oo][Tt][Aa][Ll][ |\\.|_]*([Nn][Oo]|[Nn][Uu][Mm]|#)[ |\\.|_]*([Cc]|[Cc][Aa][Rr][Bb][Oo][Nn])",colnames(db.rawdata))])
  
  object@FA_total_no_DB = as.integer(db.rawdata[,grepl("[Ff][Aa][ |\\.|_]*[Tt][Oo][Tt][Aa][Ll][ |\\.|_]*([Nn][Oo]|[Nn][Uu][Mm]|#)[ |\\.|_]*([Dd][Bb]|[Dd][Oo][Uu][Bb][Ll][Ee][ |\\.|_]*[Bb][Oo][Nn][Dd][Ss])",colnames(db.rawdata))])
  
  object@degree_oxidation = as.integer(db.rawdata[,grepl("([Dd][Ee][Gg][Rr][Ee][Ee]|[Dd][Ee][Gg])[ |\\.|_]*([Oo][Xx][Ii][Dd][Aa][Tt][Ii][Oo][Nn]|[Oo][Xx])",colnames(db.rawdata))])
  
  object@parent_elem_formula = as.character(db.rawdata[,grepl("[Pp][Aa][Rr][Ee][Nn][Tt][ |\\.|_]*([Ee][Ll][Ee][Mm][Ee][Nn][Tt][Aa][Ll]|[Ee][Ll][Ee][Mm])[ |\\.|_]*[Ff][Oo][Rr][Mm][Uu][Ll][Aa]",colnames(db.rawdata))])
  
  object@parent_compound_name = as.character(db.rawdata[,grepl("[Pp][Aa][Rr][Ee][Nn][Tt][ |\\.|_]*[Cc][Oo][Mm][Pp][Oo][Uu][Nn][Dd][ |\\.|_]*[Nn][Aa][Mm][Ee]",colnames(db.rawdata))])
  
  # try to verify user's indication of polarity
  
  if (polarity=="positive") {
    
    if (sum(grepl("^\\[.*\\]\\-{1,}$",object@adduct))>0) {
      
      stop("At least one of the adducts in field 'adduct' of the database being imported appears to be of ion mode opposite that of indicated polarity '",polarity,".' Check the ion mode specified. Aborting...")
      
    }
    
    if (sum(grepl("^\\[.*\\]\\+{1,}$",object@adduct))!=length(object@mz)) {
      
      warning(paste0("Could not determine that all adducts in field 'adduct' of the database being imported are of indicated polarity '",polarity,".' User may experience unexpected performance.\n"))
      
    }
    
  }
  
  if (polarity=="negative") {
    
    if (sum(grepl("^\\[.*\\]\\+{1,}$",object@adduct))>0) {
      
      stop("At least one of the adducts in field 'adduct' of the database being imported appears to be of ion mode opposite that of indicated polarity '",polarity,".' Check the ion mode specified. Aborting...")
      
    }
    
    if (sum(grepl("^\\[.*\\]\\-{1,}$",object@adduct))!=length(object@mz)) {
      
      warning(paste0("Could not determine that all adducts in field 'adduct' of the database being imported are of indicated polarity '",polarity,".' User may experience unexpected performance.\n"))
      
    }
    
  }
  
  # assign object@polarity
  
  object@polarity = as.factor(polarity)
  
  # assign object@num.entries
  
  object@num.entries = as.integer(length(object@mz))
  
  # assign object@num.compounds, if user has provided it
  
  if (!is.null(num.compounds)) {
    
    object@num.compounds = as.integer(num.compounds)
  
  }
  
  object
  
}

# getformattedIsoData: generate plain-text labels for each peakgroup in an xsAnnotate object for which findIsotopes returned data; should be called with sapply or lapply, e.g., sapply(isodata, getformattedIsoData, polarity = "positive")

# DISCLAIMER: Some code in this function was borrowed and then modified from source code for the CAMERA function getPeaklist. This code was borrowed under the GPL 2.0 license from the CAMERA source R file "xsAnnotate.R" on 12/10/2015. The modified code is provided without warranty. CAMERA is described in Kuhl, C., Tautenhahn, R., Boettcher, C., Larson, T.R., and Neumann, S., 2012, "CAMERA: an integrated strategy for compound spectra extraction and annotation of liquid chromatography/mass spectrometry data sets," Anal. Chem. 84:283–289. See http://pubs.acs.org/doi/abs/10.1021/ac202450g and https://bioconductor.org/packages/release/bioc/html/CAMERA.html

getformattedIsoData = function(isodata, polarity) { # input "isodata" is all or part of the @isotopes slot from a CAMERA "xsAnnotate" object
  
  polarity = match.arg(polarity, choices = c("positive","negative"), several.ok = FALSE)

  polsym = switch(
    polarity,
    positive="+",
    negative="-")
  
  if (!is.null(isodata)) {
    
    num.iso = isodata$iso
    
    # begin building isotope data string
    
    if (num.iso == 0) {
      
      str.iso = "[M]"
      
    } else { 
      
      str.iso = paste("[M+", num.iso, "]", sep="")
    }
    # if multiply charged
    
    if (isodata$charge > 1) {
      
      isotopes = paste("[", isodata$y, "]", str.iso, isodata$charge, polsym, sep="")
      
    } else {
      
      isotopes = paste("[", isodata$y, "]", str.iso, polsym, sep="")
      
    }
    
  } else { 
    
    # no isotope information for this peakgroup
    
    isotopes = ""
    
  }
  
}

# matchbyPPM: assigns compound identities to a dataset from a database based on a ppm restriction

matchbyPPM = function(mz, database, ppm) { # mz is list of m/z ratios from a dataset; database is a "LOBdbase" object; ppm inherited from wrapper function 
  
  matches = database@frag_ID[abs((database@mz-mz)/database@mz*1000000)<=ppm]
  
}

# evalFeatureRT: evaluates retention times of xcms peakgroups to which putative database assignments have been made against retention time window restrictions supplied by user 

evalFeatureRT = function(matched.frag_IDs, assignment.rt, rt.windows, database) { # assignment.rt is the mean retention time of a peakgroup; matched.frag_IDs are the database fragment IDs of the putative database assignment(s) matched to that peakgroup
  
  assignment.rt = assignment.rt/60 # convert observed feature rt from seconds to minutes
  
  ID.eval = rep(TRUE, length(matched.frag_IDs))
  
  if (length(matched.frag_IDs)>0) { # matches were made to this feature 
    
    for (i in 1:length(matched.frag_IDs)) {
      
      rt.windowdata = rt.windows[database@species[matched.frag_IDs[i]]==rt.windows$lipid_class,] # get any window data for this feature
      
      if (nrow(rt.windowdata)==1) { # there is rt window data for this lipid class
        
        if (!is.na(rt.windowdata$rt_win_min) & !is.na(rt.windowdata$rt_win_max)) { # we have an upper and lower bound for this class
          
          if (!(assignment.rt>rt.windowdata$rt_win_min & assignment.rt<rt.windowdata$rt_win_max)) {
            
            ID.eval[i] = FALSE
            
          }
          
        } else if (!is.na(rt.windowdata$rt_win_min) & is.na(rt.windowdata$rt_win_max)) { # we only have a lower bound
          
          if (!(assignment.rt>rt.windowdata$rt_win_min)) {
            
            ID.eval[i] = FALSE
            
          }        
          
        } else if (is.na(rt.windowdata$rt_win_min) & !is.na(rt.windowdata$rt_win_max)) { # we only have an upper bound
          
          if (!(assignment.rt<rt.windowdata$rt_win_max)) {
            
            ID.eval[i] = FALSE
            
          }        
          
        }
        
      }
      
    }
    
  }
  
  return(matched.frag_IDs[ID.eval])
  
}

# screenPSpectrum: applies screening & annotation routine to peakgroups that belong to a given CAMERA pseudospectrum 

screenPSpectrum = function(pseudospectrum, xsA, polarity, database, remove.iso, rt.restrict, rt.windows, exclude.oddFA, match.ppm) { # first argument is a CAMERA pseudospectrum; second argument is the xsA object from wrapper function; others self-explanatory or passed down from wrapper function
  
  # first, get all peakgroup and isotope data associated with the pseudospectrum, appending the xcms peakgroup number, isotope data, and the pseudospectrum number
  
  pgdata = xsA@groupInfo[xsA@pspectra[[pseudospectrum]],] # get pg data
  xcms_peakgroup = xsA@pspectra[[pseudospectrum]] # get pg number, store to a vector
  CAMERA_pseudospectrum = rep(pseudospectrum,length(pgdata)/ncol(xsA@groupInfo)) # create a vector that captures this ps number
  isodata = xsA@isotopes[xcms_peakgroup] # extract isotope data
  isotopes = as.character(sapply(isodata, getformattedIsoData, polarity = polarity)) # generate isotope strings
  
  pgdata = data.frame(pgdata, xcms_peakgroup, isotopes, CAMERA_pseudospectrum, stringsAsFactors = FALSE) # combine into data frame
  
  # define a matrix for diagnostic data
  
  diagnostic.data = matrix(data = integer(0), nrow = 1, ncol = 14)
  colnames(diagnostic.data) = c("pg.init","peaks.init","postiso.pg","postiso.peaks","pg.with.initial.assign","peaks.with.initial.assign","num.initial.assigned.adducts","num.initial.assigned.parents","postrt.pg","postrt.peaks","postAIH.pg","postAIH.peaks","postevenFA.pg","postevenFA.peaks")
  
  diagnostic.data[c("pg.init","peaks.init")] = c(nrow(pgdata),sum(pgdata[,11:(10+length(sampnames(xsA@xcmsSet)))]>0))

  # take care of secondary isotopes, if user wants
  
  if (remove.iso = TRUE) {
    
    pgdata = pgdata[(pgdata$isotopes=="") | (grepl("\\[M\\]\\+$",pgdata$isotopes)),]
    
    diagnostic.data[c("postiso.pg",
                      "postiso.peaks")] = 
      c(nrow(pgdata),
        sum(pgdata[,11:(10+length(sampnames(xsA@xcmsSet)))]>0))
    
  }
  
  # get initial matches, record results
  
  init.matches = sapply(pgdata$mz, matchbyPPM, database = database, ppm = ppm) # returns list of frag_IDs of potential matches to the mz of each peakgroup in this pseudospectrum
  
  diagnostic.data[c("pg.with.initial.assign",
                    "peaks.with.initial.assign",
                    "num.initial.assigned.adducts",
                    "num.initial.assigned.parents")] = 
    c(sum(sapply(init.matches, function(x) length(x)>0)),
      sum(pgdata[sapply(init.matches, function(x) length(x)>0),11:(10+length(sampnames(xsA@xcmsSet)))]>0),
      sum(sapply(init.matches, length)),
      
      
      
      unique(database@ unlist(init.matches)
      
                                                                                                                                             )

  # apply rt restrictions, if user asked for it
  
  if (rt.restrict = TRUE) {
    
    evalFeatureRT = function(assignment.rt, matched.frag_IDs, rt.windows, database)
      
      rt.matches = mapply(evalFeatureRT, matched.frag_IDs = init.matches, assignment.rt = pgdata$rt, MoreArgs =  list(database = database, rt.windows = rt.windows))
      
      num.rt.matches = sum(sapply(rt.matches,length))
      
  }
    
    
    
  as.matrix(init.matches)
    
  }
}


  
}