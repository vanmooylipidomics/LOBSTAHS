################ Set classes, methods for LOBSet #############

# define class representation

LOBSet = setClass("LOBSet",
                  representation(peakdata = "data.frame",
                                 iso_C3r = "list",
                                 iso_C3f = "list",
                                 iso_C3c = "list",
                                 LOBscreen_diagnostics = "data.frame",
                                 LOBisoID_diagnostics = "data.frame",
                                 LOBscreen_settings = "list",
                                 polarity = "factor",
                                 sampnames = "character"
                  ),
                  
                  prototype(peakdata = data.frame(),
                            iso_C3r = list(),
                            iso_C3f = list(),
                            iso_C3c = list(),
                            LOBscreen_diagnostics = data.frame(),
                            LOBisoID_diagnostics = data.frame(),
                            LOBscreen_settings = list(),
                            polarity = factor(character(0)),
                            sampnames = character(0)
                  ))

# create constructor

LOBSet = function(peakdata = NULL, iso_C3r = NULL, iso_C3f = NULL,
                  iso_C3c = NULL, LOBscreen_diagnostics = NULL,
                  LOBisoID_diagnostics = NULL, LOBscreen_settings = NULL,
                  polarity = c("positive","negative"), sampnames = NULL) {
  
  # first, warn user that it's better to use doLOBscreen
  
  warning("The LOBSet constructor should only be used for manual construction ",
          "of a LOBSet. Under almost all circumstances, the function ",
          "doLOBscreen should be used to create a LOBSet by apply the ",
          "LOBSTAHS screening criteria to an xsAnnotate object.")
  
  # create new (empty) LOBSet object
  
  object = new("LOBSet")
  
  # check arguments as best as possible 
  
  # peakdata
  
  if (!is.null(peakdata)) {
    
    if (class(peakdata)=="data.frame") {
      
      # spot check of the minimum necessary column names; leaving this non-
      # all-inclusive to account for future development of new classes or data
      # fields
      
      if (!all(c("peakgroup_mz","peakgroup_rt","LOBdbase_frag_ID",
            "LOBdbase_exact_parent_neutral_mass","LOBdbase_mz","lipid_class",
            "LOBdbase_ppm_match","species","major_adduct","FA_total_no_C",
            "FA_total_no_DB","degree_oxidation","elem_formula","compound_name",
            "C1","casecodes","LOBdbase_ppm_match","match_ID") %in% 
          colnames(peakdata))) {
        
        stop("Input 'peakdata' does not appear to have the minimum necessary ",
             "elements.\n")
        
      }
      
    } else {
      
      stop("Input 'peakdata' is not of class 'data.frame'\n")
      
    }
    
  }
        
  # iso_C3r
  
  if (!is.null(iso_C3r)) {
    
    if (!(class(iso_C3r)=="list")) {
      
      stop("Input 'iso_C3r' is not of class 'list'\n")
      
    }
    
  }
  
  # iso_C3f
  
  if (!is.null(iso_C3f)) {
    
    if (!(class(iso_C3f)=="list")) {
      
      stop("Input 'iso_C3f' is not of class 'list'\n")
      
    }
    
  }
  
  # iso_C3c
  
  if (!is.null(iso_C3c)) {
    
    if (!(class(iso_C3c)=="list")) {
      
      stop("Input 'iso_C3c' is not of class 'list'\n")
      
    }
    
  }
  
  # LOBscreen_diagnostics
  
  if (!is.null(LOBscreen_diagnostics)) {
    
    if (class(LOBscreen_diagnostics)=="data.frame") {
      
      # spot check of the minimum necessary column names; leaving this non-
      # all-inclusive to account for future development of new classes or data
      # fields
      
      if ((!all(c("peakgroups","peaks","assignments",
                 "parent_compounds") %in% colnames(LOBscreen_diagnostics))) |
          (!all(apply(LOBscreen_diagnostics,2,function(x) class(x)=="numeric")))
      ) {
        
        stop("Input 'LOBscreen_diagnostics' does not appear to contain data ",
             "of the correct type.\n")
        
        }
      
      } else {
      
      stop("Input 'LOBscreen_diagnostics' is not of class 'data.frame'\n")
        
      }
    
  }
  
  # LOBisoID_diagnostics
  
  if (!is.null(LOBisoID_diagnostics)) {
    
    if (class(LOBisoID_diagnostics)=="data.frame") {
      
      # spot check of the minimum necessary column names; leaving this non-
      # all-inclusive to account for future development of new classes or data
      # fields
      
      if ((!all(c("peakgroups","parent_compounds","assignments",
                  "features") %in% colnames(LOBisoID_diagnostics))) |
          (!all(apply(LOBisoID_diagnostics,2,function(x) class(x)=="numeric")))
      ) {
        
        stop("Input 'LOBisoID_diagnostics' does not appear to contain data of ",
             "the correct type.\n")
        
      }
      
    } else {
      
      stop("Input 'LOBisoID_diagnostics' is not of class 'data.frame'\n")
      
    }
    
  }

  # LOBscreen_settings
  
  if (!is.null(LOBscreen_settings)) {
    
    if (class(LOBscreen_settings)=="list") {
      
      # spot check of the minimum information that should have been recorded if
      # this is a LOBSet
      
      if (!all(c("database","remove.iso","rt.restrict",
                  "rt.windows","exclude.oddFA","match.ppm") %in%
               names(LOBscreen_settings))) {
        
        stop("Input 'LOBscreen_settings' does not appear to contain data of ",
             "the correct type.\n")
        
      }
      
    } else {
        
        stop("Input 'LOBscreen_settings' is not of class 'list'\n")
      
    }
    
  }
  
  polarity = match.arg(polarity, several.ok = FALSE)
  
  return(object)
  
}
  
# set generics and accessors

setMethod("show", "LOBSet", function(object) {
  
  if (nrow(peakdata(object))>0) {
    
    cat("A",as.character(polarity(object)),
        "polarity \"LOBSet\" containing LC-MS peak data. Compound",
        "assignments\n",
        "and adduct ion hierarchy screening annotations applied to",
        length(sampnames(object)),"samples using the\n",
        "\"LOBSTAHS\" package.\n\n")
    
    if (!is.null(LOBscreen_settings(object)$retain.unidentified) &&
        LOBscreen_settings(object)$retain.unidentified==TRUE) {
      
      # unidentified/discarded features *were* retained, provide additional
      # diagnostics
      
      cat("Unidentified features and those discarded during the LOBSTAHS",
          "screening process\n",
          "have been retained in this LOBSet. These features",
          "will show a value of 'NA' in\n",
          "the 'match_ID' field.\n\n")
      
      cat("Total no. individual peaks in this LOBSet:",
          LOBscreen_diagnostics(object)$peaks[1],"\n")
      cat("Total no. peak groups in this LOBSet:",
          length(unique(peakdata(object)$xcms_peakgroup)),"\n")
      
      if (.hasSlot(object, "LOBisoID_diagnostics")) {
        
        # can assume it is an newer LOBSet that has underscores for column names
        # in peakdata
        
        cat("m/z range of all features:",
            paste(min(peakdata(object)$peakgroup_mz, na.rm = TRUE),
                  max(peakdata(object)$peakgroup_mz, na.rm = TRUE),
                  sep = "-"),"\n\n")
        
      } else if (.hasSlot(object, "LOBisoID.diagnostics")) {
        
        # can assume it is an older LOBSet that has periods for column names in
        # peakdata instead of underscores
        
        cat("m/z range of all features:",
            paste(min(peakdata(object)$peakgroup.mz, na.rm = TRUE),
                  max(peakdata(object)$peakgroup.mz, na.rm = TRUE),
                  sep = "-"),"\n\n")
        
      }
      
    }
    
    cat("No. individual peaks with LOBSTAHS compound assignments:",
        LOBscreen_diagnostics(object)$peaks[6],"\n")
    cat("No. peak groups with LOBSTAHS compound assignments:",
        LOBscreen_diagnostics(object)$peakgroups[6],"\n")
    cat("No. LOBSTAHS compound assignments:",
        LOBscreen_diagnostics(object)$parent_compounds[6],
        "\n")
    
    if (.hasSlot(object, "LOBisoID_diagnostics")) {
      
      # can assume it is an newer LOBSet that has underscores for column names
      # in peakdata
      
      cat("m/z range of features identified using LOBSTAHS:",
          paste(min(peakdata(object)$peakgroup_mz[
            !is.na(peakdata(object)$match_ID)], na.rm = TRUE),
            max(peakdata(object)$peakgroup_mz[
              !is.na(peakdata(object)$match_ID)], na.rm = TRUE),
            sep = "-"),"\n\n")
      
    } else if (.hasSlot(object, "LOBisoID.diagnostics")) {
      
      # can assume it is an older LOBSet that has periods for column names in
      # peakdata instead of underscores
      
      cat("m/z range of features identified using LOBSTAHS:",
          paste(min(peakdata(object)$peakgroup.mz[
            !is.na(peakdata(object)$match.ID)], na.rm = TRUE),
            max(peakdata(object)$peakgroup.mz[
              !is.na(peakdata(object)$match.ID)], na.rm = TRUE),
            sep = "-"),"\n\n")
      
    }
      
    cat("Identified peak groups having possible regisomers:",
        paste(LOBisoID_diagnostics(object)$peakgroups[1],"\n"))
    cat("Identified peak groups having possible structural functional isomers:",
        paste(LOBisoID_diagnostics(object)$peakgroups[2],"\n"))
    cat("Identified peak groups having isobars indistinguishable within ppm",
        "matching\n",
        "tolerance:",
        paste(LOBisoID_diagnostics(object)$peakgroups[3],"\n\n"))
    
    cat("Restrictions applied prior to conducting adduct ion hierarchy",
        "screening:\n",
        paste(c("remove.iso","rt.restrict","exclude.oddFA")[
          unlist(LOBscreen_settings(object)[
            c("remove.iso","rt.restrict","exclude.oddFA")])], 
          collapse = ", "),"\n\n")
    
    cat("Match tolerance used for LOBSTAHS database assignments:",
        LOBscreen_settings(object)$match.ppm,"ppm\n\n")
    #   cat("Ranges of chemical parameters represented in molecules other than",
    #       "pigments:\n\n")
    #   cat("Total number of acyl carbon atoms:",
    #       paste(min(peakdata(object)$FA_total_no_C,
    #                 na.rm = TRUE),
    #             max(peakdata(object)$FA_total_no_C,
    #                 na.rm = TRUE), sep = "-"),"\n")
    #   cat("Total number of acyl carbon-carbon double bonds:",
    #       paste(min(peakdata(object)$FA_total_no_DB, na.rm = TRUE),
    #             max(peakdata(object)$FA_total_no_DB, na.rm = TRUE), 
    #                   sep = "-"),"\n")
    #   cat("Number of additional oxygen atoms:",
    #       paste(min(peakdata(object)$degree_oxidation, na.rm = TRUE),
    #             max(peakdata(object)$degree_oxidation, na.rm = TRUE), 
    #                   sep = "-"),"\n\n")
    
  } else {
    
    cat("\nA \"LOBSet\" object that appears to be empty.\n\n")
    
  }
  
  memsize = object.size(object)
  cat("Memory usage:", signif(memsize/2^20, 3), "MB\n")
  
})

setGeneric("LOBscreen_diagnostics",
           function(object) standardGeneric("LOBscreen_diagnostics"))

setMethod("LOBscreen_diagnostics", "LOBSet",
          function(object) {
            
            if (.hasSlot(object, "LOBscreen_diagnostics")) {
              
              object@LOBscreen_diagnostics

            } else if (.hasSlot(object, "LOBscreen.diagnostics")) {
              
              # likely a LOBSet created before slot names were changed from
              # dots to underscores; adjust accordingly

              object@LOBscreen.diagnostics
              
            }
            
          })

setGeneric(
  "LOBscreen_diagnostics<-", 
  function(object, value) standardGeneric("LOBscreen_diagnostics<-"))

setReplaceMethod("LOBscreen_diagnostics", "LOBSet", 
                 function(object, value) {
                   
                   object@LOBscreen_diagnostics <- value
                   
                   object
                 })

setGeneric("LOBisoID_diagnostics",
           function(object) standardGeneric("LOBisoID_diagnostics"))

setMethod("LOBisoID_diagnostics", "LOBSet",
          function(object) {
            
            if (.hasSlot(object, "LOBisoID_diagnostics")) {
              
              object@LOBisoID_diagnostics
              
            } else if (.hasSlot(object, "LOBisoID.diagnostics")) {
              
              # likely a LOBSet created before slot names were changed from
              # dots to underscores; adjust accordingly
              
              object@LOBisoID.diagnostics
              
            }
            
          })

setGeneric(
  "LOBisoID_diagnostics<-", 
  function(object, value) standardGeneric("LOBisoID_diagnostics<-"))

setReplaceMethod("LOBisoID_diagnostics", "LOBSet", 
                 function(object, value) {
                   
                   object@LOBisoID_diagnostics <- value
                   
                   object
                 })


setGeneric("polarity", function(object) standardGeneric("polarity"))

setMethod("polarity", "LOBSet", function(object) object@polarity)

setGeneric(
  "polarity<-", 
  function(object, value) standardGeneric("polarity<-"))

setReplaceMethod("polarity", "LOBSet", 
                 function(object, value) {
                   
                   object@polarity <- value
                   
                   object
                 })

setGeneric("peakdata", function(object) standardGeneric("peakdata"))

setMethod("peakdata", "LOBSet", function(object) object@peakdata)

setGeneric(
  "peakdata<-", 
  function(object, value) standardGeneric("peakdata<-"))

setReplaceMethod("peakdata", "LOBSet", 
                 function(object, value) {
                   
                   object@peakdata <- value
                   
                   object
                 })

setGeneric("sampnames", function(object) standardGeneric("sampnames"))

setMethod("sampnames", "LOBSet", function(object) object@sampnames)

setGeneric(
  "sampnames<-", 
  function(object, value) standardGeneric("sampnames<-"))

setReplaceMethod("sampnames", "LOBSet", 
                 function(object, value) {
                   
                   object@sampnames <- value
                   
                   object
                 })

setGeneric("LOBscreen_settings",
           function(object) standardGeneric("LOBscreen_settings"))

setMethod("LOBscreen_settings", "LOBSet",
          function(object) {
            
            if (.hasSlot(object, "LOBscreen_settings")) {
              
              object@LOBscreen_settings
              
            } else if (.hasSlot(object, "LOBscreen.settings")) {
              
              # likely a LOBSet created before slot names were changed from
              # dots to underscores; adjust accordingly
              
              object@LOBscreen.settings
              
            }
            
          })

setGeneric(
  "LOBscreen_settings<-", 
  function(object, value) standardGeneric("LOBscreen_settings<-"))

setReplaceMethod("LOBscreen_settings", "LOBSet", 
                 function(object, value) {
                   
                   object@LOBscreen_settings <- value
                   
                   object
                 })