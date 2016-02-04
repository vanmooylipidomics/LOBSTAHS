################ Wrapper function #############

# doLOBscreen: Wrapper function for LOBSTAHS screening & annotation of an xsAnnotate object; returns a LOBSet object

doLOBscreen = function(xsA, polarity = NULL, database = NULL, remove.iso = TRUE, rt.restrict =  TRUE, rt.windows = NULL, exclude.oddFA = TRUE, match.ppm = NULL) { # planning to add nSlaves option at some point

  cat("\n")

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

  # check for match.ppm input

  if (is.null(match.ppm)) {

    stop("User did not specify a ppm matching tolerance.")

  }

  # check for consistent polarity between xsA and value of argument 'polarity' provided by user

  if (is.null(polarity)) { # user didn't provide value for argument, try to automatically detect polarity from xsA@polarity

    cat("User did not specify a value for argument 'polarity.' Attempting to detect current polarity from property of input 'xsA'...\n")

    if (!is.null(xsA@polarity)) {

      if (is.element(xsA@polarity,c("positive","negative"))) {

        polarity = xsA@polarity
        cat("Input 'xsA' appears to be of polarity '",polarity,"'\n\n")

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

    cat("User did not specify an external database. Using default LOBSTAHS database for polarity '",polarity,"'\n\n")

    defDB = 1
    
    default.LOBdbase = NULL # to satisfy R CMD CHECK
    data(default.LOBdbase, envir = environment())
    LOBdbase = default.LOBdbase
    database = LOBdbase[[polarity]]

  }

  # set things up for rt window screening

  if (rt.restrict==TRUE) {
    
    warning("User elected screening based on retention time. If any flavor of xcms retention time correction was applied to the dataset prior to analysis with LOBSTAHS, user is advised to consider whether RT screening was a wise choice. Although LOBSTAHS adds a 10% buffer to all retention times to account for small shifts that may occur during RT correction, poor results may arise if the original data contained wide variance in retention time.\n")

    if (is.null(rt.windows)) { # use defaults
      
      default.rt.windows = NULL # to satisfy R CMD CHECK
      data(default.rt.windows, envir = environment())
      rt.windows = default.rt.windows
      defRTwin = 1

      warning("User did not specify an external source of retention time window data. Default rt windows for the Van Mooy Lab Orbitrap were used. Non-VML users are strongly urged to use argument 'rt.wintable' to load retention time data specific to their chromatography.\n")

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

        warning("Retention time data should be specified in minutes, not seconds. User's data appear to be in seconds.\n\n")

      }

    }

  }

  # lastly, if user has elected remove.iso, make sure there are actually isotopes ID'd in the xsAnnotate object

  if ((remove.iso==TRUE) & (length(xsA@isotopes)<1)) {

    stop("Input 'xsA' does not appear to contain any identified isotopes. Re-run CAMERA and use findIsotopes!")

  }

  # # remove any holdovers from a previous run, if they exist
  #
  # if (exists("screened.peaktable")) {
  #
  #   rm(screened.peaktable)
  #
  # }
  #
  # if (exists("isodata.C3r")) {
  #
  #   rm(isodata.C3r)
  #
  # }
  #
  # if (exists("diagnostic.data")) {
  #
  #   rm(diagnostic.data)
  #
  # }

  ################ Perform screening #############

  cat("Performing initial screening and annotation of dataset. Compound assignments will be made from database with",match.ppm,"ppm match tolerance...\n")

  pspectra = 1:length(xsA@pspectra)

  casecodes = c("C1","C1x","C2a","C2b","C3c","C3f","C3r","C4","C5","C6a","C6b") # the possible case codes with which we'll be annotating each putative parent compound identification based on the data for its various adduct assignments

  screened.pspecdata = lapply(pspectra, screenPSpectrum, xsA = xsA, polarity = polarity, database = database, remove.iso = remove.iso, rt.restrict = rt.restrict, rt.windows = rt.windows, exclude.oddFA = exclude.oddFA, match.ppm = match.ppm, casecodes = casecodes)

  # extract & reformat results, plus extract some aggregate diagnostics

  screenedpeaks = as.data.frame(do.call("rbind", lapply(screened.pspecdata, function(x) x[["screened.peaktable"]])))
  colnames(screenedpeaks)[1:6] = c("peakgroup.mz","peakgroup.mzmin","peakgroup.mzmax","peakgroup.rt","peakgroup.rtmin","peakgroup.rtmax") # to clarify distinction between DB values and observed data
  isodata.C3r = unlist(lapply(screened.pspecdata, function(x) x[["isodata.C3r"]]), recursive = FALSE)
  LOBscreen.diagnostics = Reduce("+",lapply(screened.pspecdata, function(x) x[["diagnostic.data"]]))
  
  # calculate the match delta ppm for each item in screenedpeaks

  screenedpeaks$LOBdbase.ppm.match = (screenedpeaks$LOBdbase.mz-screenedpeaks$peakgroup.mz)/screenedpeaks$LOBdbase.mz*1000000

  # append a match_ID

  screenedpeaks$match_ID = 1:nrow(screenedpeaks)

  cat("Initial screening and annotation complete.",LOBscreen.diagnostics[c("post.AIHscreen"),c("peakgroups")],"peakgroups remain in dataset, to which",LOBscreen.diagnostics[c("post.AIHscreen"),c("parent_compounds")],"parent compound identities have been assigned from database.\n\n")

  # further screening to identify isomers
  
  num.treatments = length(levels(xsA@xcmsSet@phenoData$class)) # first, get number of sample treatments in original dataset (i.e., xcms "classes")

  # create a matrix to keep track of isomer data

  LOBisoID.diagnostics = data.frame(matrix(data = NA, nrow = 3, ncol = 4))
  colnames(LOBisoID.diagnostics) = c("peakgroups","parent_compounds","assignments","features")
  rownames(LOBisoID.diagnostics) = c("C3r_regio.iso","C3f_funct.struct.iso","C3c_isobars")

  # check for & tag any additional "case 3r" regioisomers

  cat("Identifying and annotating possible regioisomers...\n")

  screenedpeaks$C3r[screenedpeaks$compound_name %in% screenedpeaks$compound_name[duplicated(screenedpeaks$compound_name)]] = 1 # set "C3r" for all compounds with regioisomers to 1

  C3r.peakdata = screenedpeaks[screenedpeaks$C3r==1,]
  C3r.parent_compounds = unique(C3r.peakdata$compound_name)

  # update dataset.isodata.C3r with correct cross-references

  for (i in 1:length(C3r.parent_compounds)) {

    pg.this.parent = screenedpeaks[screenedpeaks$compound_name==C3r.parent_compounds[i],]

    isodata.C3r[screenedpeaks$compound_name==C3r.parent_compounds[i]] = rep(list(pg.this.parent$match_ID),nrow(pg.this.parent))

  }

  LOBisoID.diagnostics[c("C3r_regio.iso"),c("peakgroups","parent_compounds","assignments","features")] =
    c(length(unique(screenedpeaks[screenedpeaks$C3r==1,c("xcms_peakgroup")])),
      length(unique(screenedpeaks$compound_name[screenedpeaks$C3r==1])),
      nrow(screenedpeaks[screenedpeaks$C3r==1,]),
      sum(apply(screenedpeaks[(!duplicated(screenedpeaks$xcms_peakgroup) & screenedpeaks$C3r==1),(12+num.treatments):(11+num.treatments+length(sampnames(xsA@xcmsSet)))], c(1,2), function(x) x>0))
      )

  cat("Found",LOBisoID.diagnostics$peakgroups[1],"possible regioisomers of",LOBisoID.diagnostics$parent_compounds[1],"parent compounds.\n\n")

  # check for & tag isobars (case C3c) and functional structural isomers (case C3f)

  cat("Identifying and annotating isobars and possible functional structural isomers...\n")

  # create some lists to hold the C3f isomer/C3c isobar information
  isodata.C3c = vector(mode = "list", length = nrow(screenedpeaks))
  isodata.C3f = vector(mode = "list", length = nrow(screenedpeaks))

  C3f_C3c.peakdata = screenedpeaks[screenedpeaks$xcms_peakgroup %in% screenedpeaks$xcms_peakgroup[duplicated(screenedpeaks$xcms_peakgroup)],]
  C3f_C3c.peakgroups = unique(C3f_C3c.peakdata$xcms_peakgroup)

  for (i in 1:length(C3f_C3c.peakgroups)) {

    IDs.this.pg = screenedpeaks[screenedpeaks$xcms_peakgroup==C3f_C3c.peakgroups[i],]

    parent.mzs.this.pg = unique(IDs.this.pg$LOBdbase.mz)

    if (length(parent.mzs.this.pg)>1) { # we have a C3c scenario

      screenedpeaks$C3c[screenedpeaks$xcms_peakgroup==C3f_C3c.peakgroups[i]] = 1
      isodata.C3c[screenedpeaks$xcms_peakgroup==C3f_C3c.peakgroups[i]] = rep(list(IDs.this.pg$match_ID),nrow(IDs.this.pg))

    }

    for (j in 1:length(parent.mzs.this.pg)) {

      if (nrow(IDs.this.pg[IDs.this.pg$LOBdbase.mz==parent.mzs.this.pg[j],])>1) { # we have at least 2 functional structural isomers that go together

        screenedpeaks$C3f[screenedpeaks$match_ID %in% IDs.this.pg$match_ID[IDs.this.pg$LOBdbase.mz==parent.mzs.this.pg[j]]] = 1
        isodata.C3f[screenedpeaks$match_ID %in% IDs.this.pg$match_ID[IDs.this.pg$LOBdbase.mz==parent.mzs.this.pg[j]]] = rep(list(IDs.this.pg$match_ID[IDs.this.pg$LOBdbase.mz==parent.mzs.this.pg[j]]),nrow(IDs.this.pg[IDs.this.pg$LOBdbase.mz==parent.mzs.this.pg[j],]))

      }

    }

  }

  LOBisoID.diagnostics[c("C3f_funct.struct.iso"),c("peakgroups","parent_compounds","assignments","features")] =
    c(length(unique(screenedpeaks[screenedpeaks$C3f==1,c("xcms_peakgroup")])),
      length(unique(screenedpeaks$compound_name[screenedpeaks$C3f==1])),
      nrow(screenedpeaks[screenedpeaks$C3f==1,]),
      sum(apply(screenedpeaks[(!duplicated(screenedpeaks$xcms_peakgroup) & screenedpeaks$C3f==1),(12+num.treatments):(11+num.treatments+length(sampnames(xsA@xcmsSet)))], c(1,2), function(x) x>0))
    )

  LOBisoID.diagnostics[c("C3c_isobars"),c("peakgroups","parent_compounds","assignments","features")] =
    c(length(unique(screenedpeaks[screenedpeaks$C3c==1,c("xcms_peakgroup")])),
      length(unique(screenedpeaks$compound_name[screenedpeaks$C3c==1])),
      nrow(screenedpeaks[screenedpeaks$C3c==1,]),
      sum(apply(screenedpeaks[(!duplicated(screenedpeaks$xcms_peakgroup) & screenedpeaks$C3c==1),(12+num.treatments):(11+num.treatments+length(sampnames(xsA@xcmsSet)))], c(1,2), function(x) x>0))
    )
  
  cat("Found",LOBisoID.diagnostics$peakgroups[2],"functional structural isomers and",LOBisoID.diagnostics$peakgroups[3],"isobars, representing",sum(LOBisoID.diagnostics$parent_compounds[c(2,3)]),"parent compounds.\n")

#   cat("Found",LOBisoID.diagnostics$peakgroups[2],"functional structural isomers and",LOBisoID.diagnostics$peakgroups[3],"isobars. These isobars are compound assignments that differ in m/z so narrowly that they cannot be resolved from each other at",match.ppm,"ppm. Together, these assignments represent",sum(LOBisoID.diagnostics$parent_compounds[c(2,3)]),"different parent compounds.\n")

  # populate screenedpeaks.casecodes

  codestrings = apply(screenedpeaks[,casecodes], c(1), function (x) casecodes[x>=1])
  screenedpeaks$casecodes = unlist(lapply(codestrings, function(x) paste(x,collapse="; ")))

  # create LOBSet object for return to user

  object = new("LOBSet")

  object@peakdata = screenedpeaks
  object@polarity = as.factor(polarity)
  object@sampnames = sampnames(xsA@xcmsSet)
  object@iso.C3r = isodata.C3r
  object@iso.C3c = isodata.C3c
  object@iso.C3f = isodata.C3f
  object@LOBscreen.diagnostics = LOBscreen.diagnostics
  object@LOBisoID.diagnostics = LOBisoID.diagnostics

  if (defDB==1) {

    database.used = "default"

  } else {

    database.used = "user-supplied database"

  }

  if (rt.restrict==TRUE) {

    if (defRTwin==1) {

      rt.windows.used = "default"

    } else {

      rt.windows.used = "user-supplied rt window data"

    }

  } else {

    rt.windows.used = NULL

  }

  object@LOBscreen.settings = list(database = database.used,
                                   remove.iso = remove.iso,
                                   rt.restrict = rt.restrict,
                                   rt.windows = rt.windows.used,
                                   exclude.oddFA = exclude.oddFA,
                                   match.ppm = match.ppm)

  return(object)

}

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

      stop("At least one of the adducts in field 'adduct' of the database being imported appears to be of ionization mode opposite that of indicated polarity '",polarity,".' Check the ionization mode specified. Aborting...\n")

    }

    if (sum(grepl("^\\[.*\\]\\+{1,}$",object@adduct))!=length(object@mz)) {

      warning(paste0("Could not determine that all adducts in field 'adduct' of the database being imported are of indicated polarity '",polarity,".' User may experience unexpected performance.\n"))

    }

  }

  if (polarity=="negative") {

    if (sum(grepl("^\\[.*\\]\\+{1,}$",object@adduct))>0) {

      stop("At least one of the adducts in field 'adduct' of the database being imported appears to be of ionization mode opposite that of indicated polarity '",polarity,".' Check the ionization mode specified. Aborting...\n")

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

  DBdata = data.frame(as.integer(database@frag_ID[database@frag_ID==frag_ID]),
                      as.numeric(database@exact_parent_neutral_mass[database@frag_ID==frag_ID]),
                      as.numeric(database@mz[database@frag_ID==frag_ID]),
                      as.character(database@lipid_class[database@frag_ID==frag_ID]),
                      as.character(database@species[database@frag_ID==frag_ID]),
                      as.character(database@adduct[database@frag_ID==frag_ID]),
                      as.integer(database@FA_total_no_C[database@frag_ID==frag_ID]),
                      as.integer(database@FA_total_no_DB[database@frag_ID==frag_ID]),
                      as.integer(database@degree_oxidation[database@frag_ID==frag_ID]),
                      as.character(database@parent_elem_formula[database@frag_ID==frag_ID]),
                      as.character(database@parent_compound_name[database@frag_ID==frag_ID]),
                      stringsAsFactors = FALSE)

  names(DBdata) = c("LOBdbase.frag_ID","LOBdbase.exact_parent_neutral_mass","LOBdbase.mz","lipid_class","species","major.adduct","FA_total_no_C","FA_total_no_DB","degree_oxidation","elem_formula","compound_name")

  DBdata

}

# evalFeatureRT: evaluates retention times of xcms peakgroups to which putative database assignments have been made against retention time window restrictions supplied by user

evalFeatureRT = function(matched.frag_IDs, assignment.rtmin, assignment.rtmax, rt.windows, database) { # assignment.rtmin and max are the min and max retention times of a peakgroup; matched.frag_IDs are the database fragment IDs of the putative database assignment(s) matched to that peakgroup
  
  if (is.null(rt.windows)) { # use defaults
    
    default.rt.windows = NULL # to satisfy R CMD CHECK
    data(default.rt.windows, envir = environment())
    rt.windows = default.rt.windows
    
  }

  # convert observed feature rts from seconds to minutes
  
  assignment.rtmin = assignment.rtmin/60
  assignment.rtmax = assignment.rtmax/60
  
  # incorporate a conservative (10%) error into each value, to account for the fact that rt's can be altered from observed values during xcms retention time correction; the extra 10% should be sufficient in most cases, unless the user fed xcms a low-quality dataset with very large variance in rt 
  
  assignment.rtmin = assignment.rtmin-assignment.rtmin*.1
  assignment.rtmax = assignment.rtmax+assignment.rtmax*.1
  
  ID.eval = rep(TRUE, length(matched.frag_IDs)) # vector of logicals to record of the compliance of each assignment

  if (length(matched.frag_IDs)>0) { # matches were made to this feature

    for (i in 1:length(matched.frag_IDs)) {

      rt.windowdata = rt.windows[database@species[database@frag_ID==matched.frag_IDs[i]]==rt.windows$lipid_class,] # get any window data for this feature

      if (nrow(rt.windowdata)==1) { # there is rt window data for this lipid class

        if (!is.na(rt.windowdata$rt_win_min) & !is.na(rt.windowdata$rt_win_max)) { # we have an upper and lower bound for this class

          if (!(assignment.rtmax>rt.windowdata$rt_win_min & assignment.rtmin<rt.windowdata$rt_win_max)) {

            ID.eval[i] = FALSE

          }

        } else if (!is.na(rt.windowdata$rt_win_min) & is.na(rt.windowdata$rt_win_max)) { # we only have a lower bound

          if (!(assignment.rtmax>rt.windowdata$rt_win_min)) {

            ID.eval[i] = FALSE

          }

        } else if (is.na(rt.windowdata$rt_win_min) & !is.na(rt.windowdata$rt_win_max)) { # we only have an upper bound

          if (!(assignment.rtmin<rt.windowdata$rt_win_max)) {

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

      if (database@lipid_class[database@frag_ID==matched.frag_IDs[i]] %in% c("IP_DAG","PUA","FFA","TAG")) {
        
        if (database@FA_total_no_C[database@frag_ID==matched.frag_IDs[i]]%%2!=0) {

          ID.eval[i] = FALSE

        }

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

# getLOBpeaklist: generates a peaklist from a screened & annotated LOBSet object, with options to include isotope cross-references and export to .csv

getLOBpeaklist = function(LOBSet, include.iso = TRUE, gen.csv = FALSE) {

  if (!class(LOBSet)=="LOBSet") {

    stop("Input 'LOBSet' is not a 'LOBSet' object.\n")

  }

  export.df = LOBSet@peakdata

  # get rid of some junk, reorder some columns

  export.df = export.df[,-c(which(colnames(export.df) %in% c("npeaks","isotopes")))]
  leadcols = export.df[,c("match_ID","compound_name","elem_formula","LOBdbase.mz","peakgroup.mz","LOBdbase.ppm.match","peakgroup.rt")]
  export.df = export.df[,-c(which(colnames(export.df) %in% c("match_ID","compound_name","elem_formula","LOBdbase.mz","peakgroup.mz","LOBdbase.ppm.match","peakgroup.rt")))]
  export.df = data.frame(leadcols,export.df)

  # argument-dependent options

  if (include.iso==TRUE) {

    iso.C3r.match_ID = sapply(LOBSet@iso.C3r, paste, collapse = ", ")
    iso.C3f.match_ID = sapply(LOBSet@iso.C3f, paste, collapse = ", ")
    iso.C3c.match_ID = sapply(LOBSet@iso.C3c, paste, collapse = ", ")

    export.df = data.frame(export.df,iso.C3r.match_ID,iso.C3f.match_ID,iso.C3c.match_ID, stringsAsFactors = FALSE)

  }

  if (gen.csv==TRUE) {

    output_DTG = genTimeStamp()

    fname = paste0("LOBSTAHS_screened_peakdata_",output_DTG,".csv")

    write.csv(export.df, fname, row.names = FALSE)

    cat("Peak data exported to:",fname,"\n")

  }

  return(export.df)

}

# screenPSpectrum: applies screening & annotation routine to peakgroups that belong to a given CAMERA pseudospectrum

screenPSpectrum = function(pseudospectrum, xsA, polarity, database, remove.iso, rt.restrict, rt.windows, exclude.oddFA, match.ppm, casecodes) { # first argument is a CAMERA pseudospectrum; second argument is the xsA object from wrapper function; others self-explanatory or passed down from wrapper function

  cat("Pseudospectrum:",pseudospectrum,"\n")

  # get all peakgroup and isotope data associated with the pseudospectrum, appending the xcms peakgroup number, isotope data, and the pseudospectrum number

  pgdata = xsA@groupInfo[xsA@pspectra[[pseudospectrum]],] # get pg data
  xcms_peakgroup = xsA@pspectra[[pseudospectrum]] # get pg number, store to a vector
  CAMERA_pseudospectrum = rep(pseudospectrum,length(pgdata)/ncol(xsA@groupInfo)) # create a vector that captures this ps number
  isodata = xsA@isotopes[xcms_peakgroup] # extract isotope data
  isotopes = as.character(sapply(isodata, getformattedIsoData, polarity = polarity)) # generate isotope strings
  num.treatments = length(levels(xsA@xcmsSet@phenoData$class)) # get number of sample treatments in original dataset (i.e., xcms "classes")

  if (length(pgdata)<=(7+num.treatments+length(sampnames(xsA@xcmsSet)))) {

    pgdata = t(pgdata)

  }
  
  pgdata = data.frame(pgdata, xcms_peakgroup, isotopes, CAMERA_pseudospectrum, stringsAsFactors = FALSE) # combine into data frame

  # define a matrix to hold diagnostic data elements

  diagnostic.data = data.frame(matrix(data = NA, nrow = 6, ncol = 4))
  colnames(diagnostic.data) = c("peakgroups","peaks","assignments","parent_compounds")
  rownames(diagnostic.data) = c("initial","post.remove.iso","initial.assignments","post.rt.restrict","post.exclude.oddFA","post.AIHscreen")

  diagnostic.data[c("initial"),c("peakgroups","peaks")] = c(nrow(pgdata),sum(pgdata[,(8+num.treatments):(7+num.treatments+length(sampnames(xsA@xcmsSet)))]>0))

  ################# begin pre-screening #################

  # take care of secondary isotopes, if user wants

  if (remove.iso==TRUE) {

    pgdata = pgdata[(pgdata$isotopes=="") | (grepl("\\[M\\]\\+$",pgdata$isotopes)),]

    if (nrow(pgdata) > 0) { # we still have matches

      diagnostic.data[c("post.remove.iso"),c("peakgroups","peaks")] =
        c(nrow(pgdata),
          sum(pgdata[,(8+num.treatments):(7+num.treatments+length(sampnames(xsA@xcmsSet)))]>0))

    } else {

      diagnostic.data[c("post.remove.iso"),c("peakgroups","peaks")] = c(0,0)

    }

  }

  # get initial matches, record results

  init.matches = lapply(pgdata$mz, matchbyPPM, database = database, match.ppm = match.ppm) # returns list of frag_IDs of potential matches to the mz of each peakgroup in this pseudospectrum

  if (length(unlist(init.matches)) > 0) { # we still have matches

    diagnostic.data[c("initial.assignments"),c("peakgroups","peaks","assignments","parent_compounds")] =
      c(sum(sapply(init.matches, function(x) length(x)>0)),
        sum(pgdata[sapply(init.matches, function(x) length(x)>0),(8+num.treatments):(7+num.treatments+length(sampnames(xsA@xcmsSet)))]>0),
        sum(sapply(init.matches, length)),
        length(unique(database@parent_compound_name[unlist(init.matches)])))

  } else {

    diagnostic.data[c("initial.assignments"),c("peakgroups","peaks","assignments","parent_compounds")] = c(rep(0,4))

  }

  current.matches = init.matches

  # apply rt restrictions, if user asked for it

  if (rt.restrict==TRUE) {

    rt.matches = mapply(evalFeatureRT, matched.frag_IDs = current.matches, assignment.rtmin = pgdata$rtmin, assignment.rtmax = pgdata$rtmax, MoreArgs =  list(database = database, rt.windows = rt.windows), SIMPLIFY = FALSE)

    if (length(unlist(rt.matches)) > 0) { # we still have matches

      diagnostic.data[c("post.rt.restrict"),c("peakgroups","peaks","assignments","parent_compounds")] =
        c(sum(sapply(rt.matches, function(x) length(x)>0)),
          sum(pgdata[sapply(rt.matches, function(x) length(x)>0),(8+num.treatments):(7+num.treatments+length(sampnames(xsA@xcmsSet)))]>0),
          sum(sapply(rt.matches, length)),
          length(unique(database@parent_compound_name[unlist(rt.matches)])))

    } else {

      diagnostic.data[c("post.rt.restrict"),c("peakgroups","peaks","assignments","parent_compounds")] = c(rep(0,4))

    }

    current.matches = rt.matches

  }

  # apply even-odd FA chain length criteria, if user asked for it

  if (exclude.oddFA==TRUE) {

    evenFA.matches = lapply(current.matches, excludeoddFAlength, database = database)

    if (length(unlist(evenFA.matches)) > 0) { # we still have matches

      diagnostic.data[c("post.exclude.oddFA"),c("peakgroups","peaks","assignments","parent_compounds")] =
        c(sum(sapply(evenFA.matches, function(x) length(x)>0)),
          sum(pgdata[sapply(evenFA.matches, function(x) length(x)>0),(8+num.treatments):(7+num.treatments+length(sampnames(xsA@xcmsSet)))]>0),
          sum(sapply(evenFA.matches, length)),
          length(unique(database@parent_compound_name[unlist(evenFA.matches)])))

    } else {

      diagnostic.data[c("post.exclude.oddFA"),c("peakgroups","peaks","assignments","parent_compounds")] = c(rep(0,4))

    }

    current.matches = evenFA.matches

  }

  # eliminate case C4/C5 assignments; i.e., those assignments representing adducts of lesser theoretical abundance where the most abundant adduct of the same parent is not present; should leave us with only adducts of abundance rank = 1, and those lesser adducts for which the most abundant adduct is also present

  if (length(unlist(current.matches)) > 0) { # we still have matches

    pspectrum.compound_names = database@parent_compound_name[sapply(unlist(current.matches), match, database@frag_ID)]

    current.matches = lapply(current.matches, trimAssignments, database = database, pspectrum.frag_IDs = unlist(current.matches), pspectrum.compound_names = pspectrum.compound_names)

    if (length(unlist(current.matches)) > 0) { # we still have matches

      ################# begin adduct ion hierarchy screening #################

      current_casecodes = array(NA,c(length = length(casecodes))) # create a vector to keep track of case codes
      current_casecodes[1:length(casecodes)] = 0 # set all case codes to 0 by default
      names(current_casecodes) = casecodes # label columns

      ################# resolution of case 2a/2b/3r #################

      parent.compounds.current.matches = lapply(current.matches, function(x) database@parent_compound_name[database@frag_ID %in% x]) # obtain list of parent compounds of all matches still remaining in this pseudospectrum
      multiadduct.parent.compounds = unique(unlist(parent.compounds.current.matches)[duplicated(unlist(parent.compounds.current.matches))]) # get list of compounds in the pseudospectrum for which we ID'd more than one adduct (these will be case 2a, 2b, and 3r assignments; note that we will also have to check for case 3r assignments again later, in case there are some which gap across pseudospectra)
      adducts.current.matches = lapply(current.matches, function(x) database@adduct[database@frag_ID %in% x])
      ranks.current..matches = lapply(current.matches, function(x) database@adduct_rank[database@frag_ID %in% x])

      if (length(multiadduct.parent.compounds) > 0) {  # if we don't have any potential case 2a/2b/3r species, skip this part

        for (i in 1:length(multiadduct.parent.compounds)) {
          
          observed.adducts.this.parent = mapply(function(x,y) as.character(x[y==multiadduct.parent.compounds[i]]), adducts.current.matches, parent.compounds.current.matches) # return all observed adducts of this case 2/3r parent compound appearing in this pseudospectrum

          frag_IDs.this.parent = mapply(function(x,y) as.character(x[y==multiadduct.parent.compounds[i]]), current.matches, parent.compounds.current.matches) # return all frag_IDs for this case 2/3r parent compound appearing in this pseudospectrum

          adduct_distribution.this.pspectrum = table(unlist(observed.adducts.this.parent)) # use the assignments involving this parent compound to obtain the distribution of its adducts in this pseudospectrum

          possible.adducts.this.parent = data.frame(as.character(database@adduct[database@parent_compound_name==multiadduct.parent.compounds[i]]),database@adduct_rank[database@parent_compound_name==multiadduct.parent.compounds[i]])
          colnames(possible.adducts.this.parent) = c("adduct","adduct_rank") # return list of possible adducts for this parent compound from the database

          possible.adducts.this.parent = possible.adducts.this.parent[order(possible.adducts.this.parent$adduct_rank, decreasing = TRUE, na.last = FALSE),] # reorder the list of possible adducts to iterate through them; must return any NA's (i.e., any adducts for which ranking was not given in DB) first so that rest of code works ok

          for (j in 1:nrow(possible.adducts.this.parent)) { # iterate through the list of all possible adducts from least to most abundant, checking to see whether we made an assignment for each one

            if (!is.na(adduct_distribution.this.pspectrum[as.character(possible.adducts.this.parent$adduct[j])]) && adduct_distribution.this.pspectrum[as.character(possible.adducts.this.parent$adduct[j])]>=1) { # at least one peak was identified in the data for this adduct

              possible.adducts.this.parent$present_in_pspectrum[j]=1 # mark this adduct as being present in the pseudospectrum
              possible.adducts.this.parent$num_times_present[j]=adduct_distribution.this.pspectrum[as.character(possible.adducts.this.parent$adduct[j])] # indicate number of different peakgroups identified as this adduct

            } else { # assume no peak was identified for this adduct

              possible.adducts.this.parent$present_in_pspectrum[j]=0 # mark this adduct as being absent from the data
              possible.adducts.this.parent$num_times_present[j]=0 # indicate number of different peaks identified as this adduct (0)

            }

          }

          ################# case 6 check #############

          # before continuing evaluation of case 2a/2b/3r parent compounds, first check to see whether they happen to be case 6's, then proceed accordingly

          if (!all(is.na(possible.adducts.this.parent$adduct_rank))) { # either we have a case 6a situation, or no case 6 at all; either way, proceed with evaluation of case 2a/2b/3r parent compounds

            if (any(is.na(possible.adducts.this.parent$adduct_rank))) { # we have a case 6a situation

              current_casecodes["C6a"] = 1 # note that case 6a is true

            }

            # apply subrules to the case 2a/2b/3r assignments for this parent compound
            # Note that the way i've coded this, using independent "if" statments, some subrule assignments will be nonexclusive; i.e., the script can tag a particular parent compound as both case 3r and 2a or 2b if it meets the minimum criteria for each

            # case 3r: we have multiple assignments of the dominant adduct of the same putative parent compound, and this is because we've identified > 1 peakgroup in this pseudospectrum representing the same dominant adduct

            if (possible.adducts.this.parent$num_times_present[possible.adducts.this.parent$adduct_rank==1]>1) { # parent compound is case 3r

              current_casecodes["C3r"] = 1 # note that case 3r is true

              if (sum((adduct_distribution.this.pspectrum>=1))==1) { # this is a case 1-case 3r compound

                current_casecodes["C1"] = 1 # note that case 1 is true

              }

            }

            # case 2a or 2b: we have multiple assignments for the same putative parent compound, and this is because we've identified multiple peaks representing different adducts

            if (sum((adduct_distribution.this.pspectrum>=1))>1) { # putative parent compound is case 2a/2b... but, have to determine which subtybe

              # case 2a: the different adducts identified strictly satisfy the hypothesized adduct hierarchy for this parent compound, i.e., we have identified a peakgroup in this pseudospectrum for every adduct more abundant than the least abundant adduct in the pseudospectrum (however, this does not mean we have to have identified a peakgroup in the pseudospectrum for EVERY possible adduct of this parent compound)

              adduct_theoretically_least = match(1,possible.adducts.this.parent$present_in_pspectrum) # of the adducts present in the pseudospectrum, this one should be the least abundant according to the rules; for case 2a to be satisfied, we must then have in the pseudospectrum all adducts which are supposed to be more abundant than this one

              if (sum(possible.adducts.this.parent$present_in_pspectrum[adduct_theoretically_least:nrow(possible.adducts.this.parent)]) == nrow(possible.adducts.this.parent) - adduct_theoretically_least + 1 ) { # parent compound is case 2a

                current_casecodes["C2a"] = 1 # note that case 2a is true

              }

              # case 2b: the different adducts identified do not perfectly satisfy the hypothesized adduct hierarchy for this putative parent compound, but we consider it a good parent compound match because the adduct which should be most abundant according to the rules is present in the pseudospectrum; in this case, some adducts of intermediate hypothesized abundance may be absent from the pseudospectrum while some other less abundant adducts are present, however the adduct which should be most abundant is definitely present

              if ((sum(possible.adducts.this.parent$present_in_pspectrum[adduct_theoretically_least:nrow(possible.adducts.this.parent)]) < nrow(possible.adducts.this.parent) - adduct_theoretically_least + 1 ) & possible.adducts.this.parent$present_in_pspectrum[nrow(possible.adducts.this.parent)]==1) { # parent compound is case 2b

                current_casecodes["C2b"] = 1 # note that case 2b is true

              }

            }

            # last order of business before moving onto the next case 2a/2b/3r species in this pseudospectrum: select and then write to the results array the appropriate data for the current parent compound
            # method: we will use the rules to determine the particular adduct of this parent compound from which we will pull the data
            # if the parent compound is case 2a or 2b, the adduct which should be most abundant is present in pseudospectrum --> record data for this adduct, regardless of whether it is actually the most abundant in this particular pseudospectrum; as long as we are consistent in applying this throughout the experiment, we should be ok since we're only really concerned with relative changes between samples
            # if the parent compound is case 3r, we have a slightly more complicated situation (because >= 2 regioisomers of the species were simultaneously identified in the pseudospectrum); we will record data for all of these assignments by inserting as many additional row(s) as is necessary
            # lastly, if the parent compound is case 6b, warn user

            if ((current_casecodes["C2a"]==1 | current_casecodes["C2b"]==1) & possible.adducts.this.parent$num_times_present[nrow(possible.adducts.this.parent)]==1) { # 2a or 2b, and we identified only a single peak for the adduct of theoretical greatest abundance  --> use data for this single peak assignment

              peakdata_to_record = cbind(pgdata[as.character(possible.adducts.this.parent$adduct[nrow(possible.adducts.this.parent)])==observed.adducts.this.parent,],extractLOBdbasedata(as.integer(frag_IDs.this.parent[as.character(possible.adducts.this.parent$adduct[nrow(possible.adducts.this.parent)])==observed.adducts.this.parent]),database))
              isodata.C3r_to_record = NULL

            } else if ((current_casecodes["C3r"]==1 & possible.adducts.this.parent$num_times_present[nrow(possible.adducts.this.parent)]>1)) { # C3r, with more than one assignment at different retention times for the adduct that should most abundant

              LOBdbdata.thisparent = t(sapply(as.integer(frag_IDs.this.parent[as.character(possible.adducts.this.parent$adduct[nrow(possible.adducts.this.parent)])==observed.adducts.this.parent]), extractLOBdbasedata, database = database, simplify = FALSE))
              LOBdbdata.thisparent = as.data.frame(do.call("rbind", lapply(LOBdbdata.thisparent, function(x) x)))

              peakdata_to_record = cbind(pgdata[as.character(possible.adducts.this.parent$adduct[nrow(possible.adducts.this.parent)])==observed.adducts.this.parent,],LOBdbdata.thisparent)
              isodata.C3r_to_record = as.integer(peakdata_to_record[[c("xcms_peakgroup")]])

            }

          } else { # we have a case 6b situation: unable to determine what data to record

            warning(paste0("Pseudospectrum ",pseudospectrum," contains an assignment for which no adduct hierarchy data are given in the database. Cannot determine what peak data to report for this pseudospectrum.\n"))

            current_casecodes["C6b"] = 1 # note that case 6b is true

            LOBdbdata.thisparent = t(sapply(as.integer(frag_IDs.this.parent[as.character(possible.adducts.this.parent$adduct[nrow(possible.adducts.this.parent)])==observed.adducts.this.parent]), extractLOBdbasedata, database = database, simplify = FALSE))
            LOBdbdata.thisparent = as.data.frame(do.call("rbind", lapply(LOBdbdata.thisparent, function(x) x)))

            peakdata_to_record = cbind(matrix(nrow = nrow(LOBdbdata.thisparent), ncol = 29),LOBdbdata.thisparent)
            isodata.C3r_to_record = NULL

          }

          # now, insert/append the peakdata and case codes for this parent compound into screened.peaktable and isodata.C3r

          if (exists("peakdata_to_record")) { # only force insertion if there's data

            if (!exists("screened.peaktable")) { # it's the first peaktable entry

              # create a matrix for our results, and a list for storing any C3r isodata
              
              peaktable.ncols = ncol(pgdata)+13+length(casecodes)
              
              screened.peaktable = data.frame(matrix(data = NA, nrow = nrow(peakdata_to_record), ncol = peaktable.ncols))
              
              colnames(screened.peaktable) = c(colnames(pgdata),"LOBdbase.frag_ID","LOBdbase.exact_parent_neutral_mass","LOBdbase.mz","lipid_class","species","major.adduct","FA_total_no_C","FA_total_no_DB","degree_oxidation","elem_formula","compound_name",casecodes,"casecodes","LOBdbase.ppm.match")

              # record data

              screened.peaktable[,1:(peaktable.ncols-2)] = cbind(peakdata_to_record,t(replicate(nrow(peakdata_to_record),current_casecodes)))

            } else { # it's not the first sample, so rbind

              peakdata_to_record = cbind(peakdata_to_record,
                                         t(replicate(nrow(peakdata_to_record),current_casecodes)),
                                         matrix(data = NA, nrow = nrow(peakdata_to_record), ncol = 2))

              colnames(peakdata_to_record) = colnames(screened.peaktable)

              screened.peaktable = rbind(screened.peaktable,
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

          current_casecodes[1:length(current_casecodes)] = 0

        }

      }

      ################# resolution of case 1 species #############

      single_adduct.parent.compounds = unlist(parent.compounds.current.matches)[!duplicated(unlist(parent.compounds.current.matches),fromLast = FALSE)&!duplicated(unlist(parent.compounds.current.matches),fromLast = TRUE)] # get list of compounds for which we ID'd only one adduct in this pseudospectrum (case 1 compounds; should be all those parent compounds that we did not ID as case 2a/2b/3r above)

      if (length(single_adduct.parent.compounds) > 0) {  # if we don't have any case 1 species, skip this part

        for (l in 1:length(single_adduct.parent.compounds)) { # cycle through each of the case 1 putative compounds in this sample

          frag_IDs.this.parent = mapply(function(x,y) as.character(x[y==single_adduct.parent.compounds[l]]), current.matches, parent.compounds.current.matches) # return frag_ID (there should be only one) for this case 1 parent compound appearing in this pseudospectrum

          if (length(unlist(frag_IDs.this.parent))!=1) { # check to make sure there's only one assignment for this parent compound in this particular pseudospectrum (otherwise something's wrong with the code)

            stop(paste0("More than one assignment remaining in pseudospectrum ",pseudospectrum," for a case 1 compound. Something is wrong! Please capture conditions leading to this error and send to developer.\n")) # stop script if this is the case

          }

          observed.adducts.this.parent = mapply(function(x,y) as.character(x[y==single_adduct.parent.compounds[l]]), adducts.current.matches, parent.compounds.current.matches) # return observed adduct (there should be only one) of this case 1 parent compound appearing in this pseudospectrum

          possible.adducts.this.parent = data.frame(as.character(database@adduct[database@parent_compound_name==single_adduct.parent.compounds[l]]),database@adduct_rank[database@parent_compound_name==single_adduct.parent.compounds[l]])
          colnames(possible.adducts.this.parent) = c("adduct","adduct_rank") # return list of possible adducts for this parent compound from the database
          possible.adducts.this.parent = possible.adducts.this.parent[order(possible.adducts.this.parent$adduct_rank, decreasing = TRUE, na.last = FALSE),] # reorder the list of possible adducts

          # another final check to make sure this is a case 1 species

          ################# case 6 check #############

          # before continuing, check again to see whether we have a case 6 situation, then proceed accordingly

          if (is.na(possible.adducts.this.parent$adduct_rank[possible.adducts.this.parent$adduct==unlist(observed.adducts.this.parent)])) { # we have a case 6b-1x situation; we ID'd a single peak representing only one good adduct for this parent compound, but we don't know how it ranks since there's no ranking data in the DB

            warning(paste0("Pseudospectrum ",pseudospectrum," contains an assignment for which no adduct hierarchy data are given in the database. Since only one adduct of this parent compound appears in this pseudospectrum, peak data will be reported for the lone adduct.\n"))

            current_casecodes["C1x"] = 1  # note that case 1x is true; we can't really say whether it's a case 1 or 4 since we don't know

            current_casecodes["C6b"] = 1  # note that case 6b is also true

          } else if (possible.adducts.this.parent$adduct_rank[possible.adducts.this.parent$adduct==unlist(observed.adducts.this.parent)]==1) { # this is definitely a case 1 situation, proceed

            current_casecodes["C1"] = 1 # note that case 1 is true

          } else {

            stop("Adduct rank data for this parent compound appear to be incorrect. Aborting...\n")

          }

          # last order of business before moving onto the next case 1 parent compound in this sample: select and then write to the screened.peakdata results array the data for the current parent compound
          # method: we will use the rules (pretty simple for case 1 IDs) to determine whether data gets written or not
          # if the ID is case 1 --> record data for this adduct
          # if the species is case 6b and only one adduct was identified, record data for that adduct

          if (current_casecodes["C1"]==1) { # case 1

            peakdata_to_record = cbind(pgdata[as.character(possible.adducts.this.parent$adduct[nrow(possible.adducts.this.parent)])==observed.adducts.this.parent,],extractLOBdbasedata(as.integer(frag_IDs.this.parent[as.character(possible.adducts.this.parent$adduct[nrow(possible.adducts.this.parent)])==observed.adducts.this.parent]),database))
            isodata.C3r_to_record = NULL

          } else { # we have a case 1-case 6b scenario:

            current_casecodes["C1x"] = 1  # note that case 1x is true; we can't really say whether it's a case 1 or 4 since we don't know

            current_casecodes["C6b"] = 1  # note that case 6b is true

            peakdata_to_record = cbind(pgdata[sapply(frag_IDs.this.parent, function(x) length(x)>0),],extractLOBdbasedata(as.integer(frag_IDs.this.parent[sapply(frag_IDs.this.parent, function(x) length(x)>0)]),database))
            isodata.C3r_to_record = NULL

          }

          # now, insert/append the peakdata and case codes for this parent compound into our peakdata results array

          if (exists("peakdata_to_record")) { # only force insertion if there's data

            if (!exists("screened.peaktable")) { # it's the first peaktable entry

              # create a matrix for our results, and a list for storing any C3r isodata

              peaktable.ncols = ncol(pgdata)+13+length(casecodes)
              
              screened.peaktable = data.frame(matrix(data = NA, nrow = nrow(peakdata_to_record), ncol = peaktable.ncols))
              
              colnames(screened.peaktable) = c(colnames(pgdata),"LOBdbase.frag_ID","LOBdbase.exact_parent_neutral_mass","LOBdbase.mz","lipid_class","species","major.adduct","FA_total_no_C","FA_total_no_DB","degree_oxidation","elem_formula","compound_name",casecodes,"casecodes","LOBdbase.ppm.match")
              
              # record data
              
              screened.peaktable[,1:(peaktable.ncols-2)] = cbind(peakdata_to_record,t(replicate(nrow(peakdata_to_record),current_casecodes)))
              
            } else { # it's not the first sample, so rbind

              peakdata_to_record = cbind(peakdata_to_record,
                                         t(replicate(nrow(peakdata_to_record),current_casecodes)),
                                         matrix(data = NA, nrow = nrow(peakdata_to_record), ncol = 2))

              colnames(peakdata_to_record) = colnames(screened.peaktable)

              screened.peaktable = rbind(screened.peaktable,
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

          current_casecodes[1:length(current_casecodes)] = 0

        }

      }

    }

  }

  # update our diagnostics

  if (length(unlist(current.matches)) > 0) {

    diagnostic.data[c("post.AIHscreen"),c("peakgroups","peaks","assignments","parent_compounds")] =
      c(length(unique(screened.peaktable$xcms_peakgroup)),
        sum(apply(screened.peaktable[!duplicated(screened.peaktable$xcms_peakgroup),(8+num.treatments):(7+num.treatments+length(sampnames(xsA@xcmsSet)))], c(1,2), function(x) x>0)),
        nrow(screened.peaktable),
        length(unique(screened.peaktable$compound_name)))

  } else {

    diagnostic.data[c("post.AIHscreen"),c("peakgroups","peaks","assignments","parent_compounds")] = c(0,0,NA,0)
    screened.peaktable = NULL
    isodata.C3r = NULL

  }

  # return our results

  list(diagnostic.data = diagnostic.data, screened.peaktable = screened.peaktable, isodata.C3r = isodata.C3r)

}
