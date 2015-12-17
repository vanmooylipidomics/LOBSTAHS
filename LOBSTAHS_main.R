################ Set classes, methods #############

# create a class for our screened xsAnnotate object to which compound assignments and confidence codes have been applied

xsALOB = setClass("xsALOB", contains="xsAnnotate",
                  representation(confcodes = "factor",
                                 iso.C3r = "list",
                                 iso.C3f = "list",
                                 iso.C3c = "list",
                                 LOBdbase.frag_ID = "integer",
                                 LOBdbase.exact_parent_neutral_mass = "numeric",
                                 LOBdbase.ppm.match = "numeric",
                                 lipid_class = "factor",
                                 species = "character",
                                 major.adduct = "factor",
                                 FA_total_no_C = "integer",
                                 FA_total_no_DB = "integer",
                                 degree_oxidation = "integer",
                                 elem_formula = "character",
                                 compound_name = "character"
                                 ),
                  
                  prototype(confcodes = factor(character(0)),
                            iso.C3r = list(),
                            iso.C3f = list(),
                            iso.C3c = list(),
                            frag_ID = integer(0),
                            LOBdbase.frag_ID = integer(0),
                            LOBdbase.exact_parent_neutral_mass = numeric(0),
                            LOBdbase.ppm.match = numeric(0),
                            lipid_class = factor(character(0)),
                            species = character(0),
                            major.adduct = factor(character(0)),
                            FA_total_no_C = integer(0),
                            FA_total_no_DB = integer(0),
                            degree_oxidation = integer(0),
                            elem_formula = character(0),
                            compound_name = character(0)
                  ))

################ Wrapper function #############

# makeLOBAssignments: Wrapper function for screening & annotation of an xsAnnotate object

makeLOBAssignments = function(xsA, polarity = NULL, database = NULL, remove.iso = TRUE, rt.restrict =  TRUE, rt.windows = NULL, exclude.oddFA = TRUE, match.ppm = 2.5) { # add nSlaves option at some point?
  
  ################ Define some global constants #############
  
  casecodes <<- c("C1","C1x","C2a","C2b","C3c","C3f","C3r","C4","C5","C6a","C6b") # the possible case codes with which we'll be annotating each putative parent compound identification based on the data for its various adduct assignments
  
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
  
  if (rt.restrict==TRUE) {
    
    if (rt.windows==NULL) { # use defaults
      
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
  
  if ((remove.iso==TRUE) & (length(xsA@isotopes)<1)) {
    
    stop("Input 'xsA' does not appear to contain any identified isotopes. Re-run CAMERA and use findIsotopes!")
    
  }
  
}

################ Perform screening #############

pspectra = 1:length(xsA@pspectra)

screened.pspecdata = lapply(pspectra, screenPSpectrum, xsA, polarity = polarity, database = database, remove.iso = remove.iso, rt.restrict = rt.restrict, rt.windows = rt.windows, exclude.oddFA = exclude.oddFA, match.ppm = match.ppm)



################ Helper functions (some will be private) #############

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

matchbyPPM = function(mz, database, match.ppm) { # mz is list of m/z ratios from a dataset; database is a "LOBdbase" object; match.ppm inherited from wrapper function 
  
  matches = database@frag_ID[abs((database@mz-mz)/database@mz*1000000)<=match.ppm]
  
}

# extractLOBdbasedata: extracts data from a LOBdbase object for a given frag_ID

extractLOBdbasedata = function(frag_ID, database) {
  
  DBdata = c(as.integer(database@frag_ID[database@frag_ID==frag_ID]),
             as.numeric(database@exact_parent_neutral_mass[database@frag_ID==frag_ID]),
             as.character(database@lipid_class[database@frag_ID==frag_ID]),
             as.character(database@species[database@frag_ID==frag_ID]),
             as.character(database@adduct[database@frag_ID==frag_ID]),
             as.integer(database@FA_total_no_C[database@frag_ID==frag_ID]),
             as.integer(database@FA_total_no_DB[database@frag_ID==frag_ID]),
             as.integer(database@degree_oxidation[database@frag_ID==frag_ID]),
             as.character(database@parent_elem_formula[database@frag_ID==frag_ID]),
             as.character(database@parent_compound_name[database@frag_ID==frag_ID]))
  
  names(DBdata) = c("LOBdbase.frag_ID","LOBdbase.exact_parent_neutral_mass","lipid_class","species","major.adduct","FA_total_no_C","FA_total_no_DB","degree_oxidation","elem_formula","compound_name")
  
  return(DBdata)

}

# evalFeatureRT: evaluates retention times of xcms peakgroups to which putative database assignments have been made against retention time window restrictions supplied by user 

evalFeatureRT = function(matched.frag_IDs, assignment.rt, rt.windows, database) { # assignment.rt is the mean retention time of a peakgroup; matched.frag_IDs are the database fragment IDs of the putative database assignment(s) matched to that peakgroup
  
  assignment.rt = assignment.rt/60 # convert observed feature rt from seconds to minutes
  
  ID.eval = rep(TRUE, length(matched.frag_IDs)) # vector of logicals to record of the compliance of each assignment
  
  if (length(matched.frag_IDs)>0) { # matches were made to this feature 
    
    for (i in 1:length(matched.frag_IDs)) {
      
      rt.windowdata = rt.windows[database@species[database@frag_ID==matched.frag_IDs[i]]==rt.windows$lipid_class,] # get any window data for this feature
      
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

# excludeoddFAlength: removes IP-DAG, TAG, PUA, FFA assignments which have an odd total number of fatty acid C atoms

excludeoddFAlength = function(matched.frag_IDs, database) {
  
  ID.eval = rep(TRUE, length(matched.frag_IDs)) # vector of logicals to record of the compliance of each assignment
  
  if (length(matched.frag_IDs)>0) { # matches were made to this feature 
    
    for (i in 1:length(matched.frag_IDs)) {
      
      if (database@FA_total_no_C[database@frag_ID==matched.frag_IDs[i]]%%2!=0) {
        
        ID.eval[i] = FALSE
        
      }
      
    }
    
  }
  
  return(matched.frag_IDs[ID.eval])
  
}

# trimAssignments: for a given compound, removes any C4/C5 assignments, i.e., those representing adducts of lesser theoretical abundance in cases where the pseudospectrum doesn't also contain an assignment representing the adduct of greatest theoretical abundance

trimAssignments = function(matched.frag_IDs, database, pspectrum.frag_IDs, pspectrum.compound_names) { # matched.frag_IDs are the database fragment IDs of the putative database assignment(s) matched to a given peakgroup; pspectrum.frag_IDs is a list of all the fragment IDs of assignments in the entire pseudospectrum; pspectrum.compound_names is a list of the parent compound names in the pseudospectrum (with duplicates)
  
  ID.eval = rep(TRUE, length(matched.frag_IDs)) # vector of logicals to record of the compliance of each assignment
  
  if (length(matched.frag_IDs)>0) { # matches were made to this feature 
    
    for (i in 1:length(matched.frag_IDs)) {
      
      if (database@adduct_rank[database@frag_ID==matched.frag_IDs[i]]==1) {
        
        ID.eval[i] = TRUE
        
      } else {
        
        other.assignments.same.parent = pspectrum.frag_IDs[pspectrum.compound_names==database@parent_compound_name[database@frag_ID==matched.frag_IDs[i]]]
        other.adduct.ranks.this.parent = database@adduct_rank[sapply(other.assignments.same.parent, match, database@frag_ID)]
        
        if (any(other.adduct.ranks.this.parent==1)==TRUE) {
          
          ID.eval[i] = TRUE
          
        } else {
          
          ID.eval[i] = FALSE
          
        }
        
      }
      
    }
    
  }
  
  return(matched.frag_IDs[ID.eval])
  
}

# screenPSpectrum: applies screening & annotation routine to peakgroups that belong to a given CAMERA pseudospectrum 

screenPSpectrum = function(pseudospectrum, xsA, polarity, database, remove.iso, rt.restrict, rt.windows, exclude.oddFA, match.ppm) { # first argument is a CAMERA pseudospectrum; second argument is the xsA object from wrapper function; others self-explanatory or passed down from wrapper function
  
  # get all peakgroup and isotope data associated with the pseudospectrum, appending the xcms peakgroup number, isotope data, and the pseudospectrum number
  
  pgdata = xsA@groupInfo[xsA@pspectra[[pseudospectrum]],] # get pg data
  xcms_peakgroup = xsA@pspectra[[pseudospectrum]] # get pg number, store to a vector
  CAMERA_pseudospectrum = rep(pseudospectrum,length(pgdata)/ncol(xsA@groupInfo)) # create a vector that captures this ps number
  isodata = xsA@isotopes[xcms_peakgroup] # extract isotope data
  isotopes = as.character(sapply(isodata, getformattedIsoData, polarity = polarity)) # generate isotope strings
  
  pgdata = data.frame(pgdata, xcms_peakgroup, isotopes, CAMERA_pseudospectrum, stringsAsFactors = FALSE) # combine into data frame
  
  # define a matrix to hold diagnostic data elements
  
  diagnostic.data = data.frame(matrix(data = NA, nrow = 6, ncol = 4))
  colnames(diagnostic.data) = c("pg","peaks","adducts","parent_compounds")
  rownames(diagnostic.data) = c("initial","post.remove.iso","initial.assignments","post.rt.restrict","post.exclude.oddFA","post.AIHscreen")
  
  diagnostic.data[c("initial"),c("pg","peaks")] = c(nrow(pgdata),sum(pgdata[,11:(10+length(sampnames(xsA@xcmsSet)))]>0))
  
  ################# begin pre-screening ################# 
  
  # take care of secondary isotopes, if user wants
  
  if (remove.iso==TRUE) {
    
    pgdata = pgdata[(pgdata$isotopes=="") | (grepl("\\[M\\]\\+$",pgdata$isotopes)),]
    
    diagnostic.data[c("post.remove.iso"),c("pg","peaks")] = 
      c(nrow(pgdata),
        sum(pgdata[,11:(10+length(sampnames(xsA@xcmsSet)))]>0))
    
  }
  
  # get initial matches, record results
  
  init.matches = sapply(pgdata$mz, matchbyPPM, database = database, match.ppm = match.ppm) # returns list of frag_IDs of potential matches to the mz of each peakgroup in this pseudospectrum
  
  diagnostic.data[c("initial.assignments"),c("pg","peaks","adducts","parent_compounds")] = 
    c(sum(sapply(init.matches, function(x) length(x)>0)),
      sum(pgdata[sapply(init.matches, function(x) length(x)>0),11:(10+length(sampnames(xsA@xcmsSet)))]>0),
      sum(sapply(init.matches, length)),
      length(unique(database@parent_compound_name[unlist(init.matches)])))
  
  current.matches = init.matches
  
  # apply rt restrictions, if user asked for it
  
  if (rt.restrict==TRUE) {
    
    rt.matches = mapply(evalFeatureRT, matched.frag_IDs = current.matches, assignment.rt = pgdata$rt, MoreArgs =  list(database = database, rt.windows = rt.windows))
    
    diagnostic.data[c("post.rt.restrict"),c("pg","peaks","adducts","parent_compounds")] = 
      c(sum(sapply(rt.matches, function(x) length(x)>0)),
        sum(pgdata[sapply(rt.matches, function(x) length(x)>0),11:(10+length(sampnames(xsA@xcmsSet)))]>0),
        sum(sapply(rt.matches, length)),
        length(unique(database@parent_compound_name[unlist(rt.matches)])))
    
    current.matches = rt.matches
    
  }
  
  # apply even-odd FA chain length criteria, if user asked for it
  
  if (exclude.oddFA==TRUE) {
    
    evenFA.matches = lapply(current.matches, excludeoddFAlength, database = database)
    
    diagnostic.data[c("post.exclude.oddFA"),c("pg","peaks","adducts","parent_compounds")] = 
      c(sum(sapply(evenFA.matches, function(x) length(x)>0)),
        sum(pgdata[sapply(evenFA.matches, function(x) length(x)>0),11:(10+length(sampnames(xsA@xcmsSet)))]>0),
        sum(sapply(evenFA.matches, length)),
        length(unique(database@parent_compound_name[unlist(evenFA.matches)])))
    
    current.matches = evenFA.matches
    
  }
  
  # eliminate case C4/C5 assignments; i.e., those assignments representing adducts of lesser theoretical abundance where the most abundant adduct of the same parent is not present; should leave us with only adducts of abundance rank = 1, and those lesser adducts for which the most abundant adduct is also present
  
  pspectrum.compound_names = database@parent_compound_name[sapply(unlist(current.matches), match, database@frag_ID)]
  
  current.matches = lapply(current.matches, trimAssignments, database = database, pspectrum.frag_IDs = unlist(current.matches), pspectrum.compound_names = pspectrum.compound_names)
  
  ################# begin adduct ion hierarchy screening ################# 
  
  current_casecodes <- array(NA,c(length = length(casecodes))) # create a vector to keep track of case codes
  current_casecodes[1:length(casecodes)] <- 0 # set all case codes to 0 by default
  names(current_casecodes) <- casecodes # label columns
  
  row.counter = 0 # set a counter so we know how to build screened.peaktable
  
  ################# case 2a/2b/3r ################# 
  
  parent.compounds.current.matches = lapply(current.matches, function(x) database@parent_compound_name[database@frag_ID %in% x]) # obtain list of parent compounds of all matches still remaining in this pseudospectrum
  multiadduct.parent.compounds = unique(unlist(parent.compounds.current.matches)[duplicated(unlist(parent.compounds.current.matches))]) # get list of compounds still remaining in the pseudospectrum which appear in assignments to more than one peakgroup (these will be case 2a, 2b, and 3r assignments; note that we will also have to check for case 3r assignments again later, in case there are some which gap across pseudospectra) 
  adducts.current.matches = lapply(current.matches, function(x) database@adduct[database@frag_ID %in% x])
  ranks.current..matches = lapply(current.matches, function(x) database@adduct_rank[database@frag_ID %in% x])
  
  if (length(multiadduct.parent.compounds) > 0) {  # if we don't have any potential case 2a/2b/3r species, skip this part
    
    for (i in 1:length(multiadduct.parent.compounds)) {
      
      observed.adducts.this.parent = mapply(function(x,y) as.character(x[y==multiadduct.parent.compounds[i]]), adducts.current.matches, parent.compounds.current.matches) # return all observed adducts of this case 2/3r parent compound appearing in this pseudospectrum
      
      print(multiadduct.parent.compounds[i])
      print(unlist(observed.adducts.this.parent))
      
      frag_IDs.this.parent = mapply(function(x,y) as.character(x[y==multiadduct.parent.compounds[i]]), current.matches, parent.compounds.current.matches) # return all frag_IDs for this case 2/3r parent compound appearing in this pseudospectrum
      
      adduct_distribution.this.pspectrum <- table(unlist(observed.adducts.this.parent)) # use the assignments involving this parent compound to obtain the distribution of its adducts in this pseudospectrum
      
      possible.adducts.this.parent = data.frame(as.character(database@adduct[database@parent_compound_name==multiadduct.parent.compounds[i]]),database@adduct_rank[database@parent_compound_name==multiadduct.parent.compounds[i]])
      colnames(possible.adducts.this.parent) = c("adduct","adduct_rank") # return list of possible adducts for this parent compound from the database
      
      possible.adducts.this.parent <- possible.adducts.this.parent[order(possible.adducts.this.parent$adduct_rank, decreasing = TRUE, na.last = FALSE),] # reorder the list of possible adducts to iterate through them; must return any NA's (i.e., any adducts for which ranking was not given in DB) first so that rest of code works ok
      
      for (j in 1:nrow(possible.adducts.this.parent)) { # iterate through the list of all possible adducts from least to most abundant, checking to see whether we made an assignment for each one
        
        if (!is.na(adduct_distribution.this.pspectrum[as.character(possible.adducts.this.parent$adduct[j])]) && adduct_distribution.this.pspectrum[as.character(possible.adducts.this.parent$adduct[j])]>=1) { # at least one peak was identified in the data for this adduct
          
          possible.adducts.this.parent$present_in_pspectrum[j]=1 # mark this adduct as being present in the pseudospectrum
          possible.adducts.this.parent$num_times_present[j]=adduct_distribution.this.pspectrum[as.character(possible.adducts.this.parent$adduct[j])] # indicate number of different peakgroups identified as this adduct
          
        } else { # assume no peak was identified for this adduct
          
          possible.adducts.this.parent$present_in_pspectrum[j]=0 # mark this adduct as being absent from the data
          possible.adducts.this.parent$num_times_present[j]=0 # indicate number of different peaks identified as this adduct (0)
          
        }
        
      }
      
      print(possible.adducts.this.parent)
      
      ################# case 6 check #############
      
      # before continuing evaluation of case 2a/2b/3r parent compounds, first check to see whether they happen to be case 6's, then proceed accordingly
      
      if (!all(is.na(possible.adducts.this.parent$adduct_rank))) { # either we have a case 6a situation, or no case 6 at all; either way, proceed with evaluation of case 2a/2b/3r parent compounds
        
        if (any(is.na(possible.adducts.this.parent$adduct_rank))) { # we have a case 6a situation
          
          current_casecodes["C6a"] <- 1 # note that case 6a is true
          
        }
        
        # apply subrules to the case 2a/2b/3r assignments for this parent compound
        # Note that the way i've coded this, using independent "if" statments, some subrule assignments will be nonexclusive; i.e., the script can tag a particular parent compound as both case 3r and 2a or 2b if it meets the minimum criteria for each
        
        # case 3r: we have multiple assignments of the dominant adduct of the same putative parent compound, and this is because we've identified > 1 peakgroup in this pseudospectrum representing the same dominant adduct 
        
        if (possible.adducts.this.parent$num_times_present[possible.adducts.this.parent$adduct_rank==1]>1) { # parent compound is case 3r
          
          current_casecodes["C3r"] <- 1 # note that case 3r is true
          
          if (sum((adduct_distribution.this.pspectrum>=1))==1) { # this is a case 1-case 3r compound
            
            current_casecodes["C1"] <- 1 # note that case 1 is true
            
          }
          
        }
        
        # case 2a or 2b: we have multiple assignments for the same putative parent compound, and this is because we've identified multiple peaks representing different adducts
        
        if (sum((adduct_distribution.this.pspectrum>=1))>1) { # putative parent compound is case 2a/2b... but, have to determine which subtybe
          
          # case 2a: the different adducts identified strictly satisfy the hypothesized adduct hierarchy for this parent compound, i.e., we have identified a peakgroup in this pseudospectrum for every adduct more abundant than the least abundant adduct in the pseudospectrum (however, this does not mean we have to have identified a peakgroup in the pseudospectrum for EVERY possible adduct of this parent compound)
          
          adduct_theoretically_least <- match(1,possible.adducts.this.parent$present_in_pspectrum) # of the adducts present in the pseudospectrum, this one should be the least abundant according to the rules; for case 2a to be satisfied, we must then have in the pseudospectrum all adducts which are supposed to be more abundant than this one
          
          if (sum(possible.adducts.this.parent$present_in_pspectrum[adduct_theoretically_least:nrow(possible.adducts.this.parent)]) == nrow(possible.adducts.this.parent) - adduct_theoretically_least + 1 ) { # parent compound is case 2a
            
            current_casecodes["C2a"] <- 1 # note that case 2a is true
            
          }
          
          # case 2b: the different adducts identified do not perfectly satisfy the hypothesized adduct hierarchy for this putative parent compound, but we consider it a good parent compound match because the adduct which should be most abundant according to the rules is present in the pseudospectrum; in this case, some adducts of intermediate hypothesized abundance may be absent from the pseudospectrum while some other less abundant adducts are present, however the adduct which should be most abundant is definitely present
          
          if ((sum(possible.adducts.this.parent$present_in_pspectrum[adduct_theoretically_least:nrow(possible.adducts.this.parent)]) < nrow(possible.adducts.this.parent) - adduct_theoretically_least + 1 ) & possible.adducts.this.parent$present_in_pspectrum[nrow(possible.adducts.this.parent)]==1) { # parent compound is case 2b
            
            current_casecodes["C2b"] <- 1 # note that case 2b is true
            
          }
          
        }
        
        # last order of business before moving onto the next case 2a/2b/3r species in this pseudospectrum: select and then write to the results array the appropriate data for the current parent compound
        # method: we will use the rules to determine the particular adduct of this parent compound from which we will pull the data
        # if the parent compound is case 2a or 2b, the adduct which should be most abundant is present in pseudospectrum --> record data for this adduct, regardless of whether it is actually the most abundant in this particular pseudospectrum; as long as we are consistent in applying this throughout the experiment, we should be ok since we're only really concerned with relative changes between samples
        # if the parent compound is case 3r, we have a slightly more complicated situation (because >= 2 regioisomers of the species were simultaneously identified in the pseudospectrum); we will record data for all of these assignments by inserting as many additional row(s) as is necessary
        # lastly, if the parent compound is case 6b, we'll record data for the adduct of the parent compound that is present in the pseudospectrum in highest actual abundance
        
        if ((current_casecodes["C2a"]==1 | current_casecodes["C2b"]==1) & possible.adducts.this.parent$num_times_present[nrow(possible.adducts.this.parent)]==1) { # 2a or 2b, and we identified only a single peak for the adduct of theoretical greatest abundance  --> use data for this single peak assignment
          
          peakdata_to_record = cbind(pgdata[as.character(possible.adducts.this.parent$adduct[nrow(possible.adducts.this.parent)])==observed.adducts.this.parent,],t(extractLOBdbasedata(as.integer(frag_IDs.this.parent[as.character(possible.adducts.this.parent$adduct[nrow(possible.adducts.this.parent)])==observed.adducts.this.parent]),database)))
          isodata.C3r_to_record = NULL
          
        } else if ((current_casecodes["C3r"]==1 & possible.adducts.this.parent$num_times_present[nrow(possible.adducts.this.parent)]>1)) { # C3r, with more than one assignment at different retention times for the adduct that should most abundant
          
          peakdata_to_record = cbind(pgdata[as.character(possible.adducts.this.parent$adduct[nrow(possible.adducts.this.parent)])==observed.adducts.this.parent,],t(sapply(as.integer(frag_IDs.this.parent[as.character(possible.adducts.this.parent$adduct[nrow(possible.adducts.this.parent)])==observed.adducts.this.parent]), extractLOBdbasedata, database = database)))
          isodata.C3r_to_record = as.integer(peakdata_to_record[[c("xcms_peakgroup")]])
          
        }
        
      } else { # we have a case 6b situation: unable to determine what data to record
        
        warning(paste0("Pseudospectrum ",pseudospectrum," contains an assignment for which no adduct hierarchy data are given in the database. Cannot determine what peak data to report for this pseudospectrum.\n"))
        
        current_casecodes["C6b"] <- 1 # note that case 6b is true
        
        LOBdbdata.thisparent = t(sapply(as.integer(frag_IDs.this.parent[as.character(possible.adducts.this.parent$adduct[nrow(possible.adducts.this.parent)])==observed.adducts.this.parent]), extractLOBdbasedata, database = database))
        
        peakdata_to_record = cbind(matrix(nrow = nrow(LOBdbdata.thisparent), ncol = 29),LOBdbdata.thisparent)
        isodata.C3r_to_record = NULL
        
      }
      
      # now, insert/append the peakdata and case codes for this parent compound into screened.peaktable and isodata.C3r
      
      if (exists("peakdata_to_record")) { # only force insertion if there's data
        
        if (!exists("screened.peaktable")) { # it's the first peaktable entry
          
          # create a matrix for our results, and a list for storing any C3r isodata
          
          screened.peaktable = data.frame(matrix(data = NA, nrow = nrow(peakdata_to_record), ncol = (ncol(pgdata)+12+length(casecodes))))
          colnames(screened.peaktable) = c(colnames(pgdata),"LOBdbase.frag_ID","LOBdbase.exact_parent_neutral_mass","lipid_class","species","major.adduct","FA_total_no_C","FA_total_no_DB","degree_oxidation","elem_formula","compound_name",casecodes,"confcodes","LOBdbase.ppm.match")
          
          # record data
          
          screened.peaktable[,1:50] = cbind(peakdata_to_record,t(replicate(nrow(peakdata_to_record),current_casecodes)))
          
        } else { # it's not the first sample, so rbind
          
          peakdata_to_record = cbind(peakdata_to_record,
                                     t(replicate(nrow(peakdata_to_record),current_casecodes)),
                                     matrix(data = NA, nrow = nrow(peakdata_to_record), ncol = 2))
          
          colnames(peakdata_to_record) = colnames(screened.peaktable)
          
          screened.peaktable <- rbind(screened.peaktable,
                                      peakdata_to_record)
          
        }
        
        if (!exists("isodata.C3r")) { # it's the first isodata.C3r entry
          
          isodata.C3r = vector("list", nrow(peakdata_to_record))
          isodata.C3r[1:nrow(peakdata_to_record)] = rep(list(isodata.C3r_to_record),nrow(peakdata_to_record))
          
        } else { # it's not the first entry
          
          isodata.C3r = append(isodata.C3r, rep(list(isodata.C3r_to_record),nrow(peakdata_to_record)))
          
        }
        
        # finally, reset/remove our casecode and peakdata_to_record vectors before proceeding
        
        rm(peakdata_to_record)
        
      }
      
      current_casecodes[1:length(current_casecodes)] <- 0
      
    }
    
  }
  