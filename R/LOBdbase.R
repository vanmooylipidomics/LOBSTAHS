################ Set classes, methods for LOBdbase #############

# define class representation

LOBdbase = setClass("LOBdbase", representation(frag_ID = "integer",
                                               mz = "numeric",
                                               exact_parent_neutral_mass = 
                                                 "numeric",
                                               lipid_class = "factor",
                                               species = "character",
                                               adduct = "factor",
                                               adduct_rank = "integer",
                                               FA_total_no_C = "integer",
                                               FA_total_no_DB = "integer",
                                               degree_oxidation = "integer",
                                               parent_elem_formula = 
                                                 "character",
                                               parent_compound_name = 
                                                 "character",
                                               polarity = "factor",
                                               num_entries = "integer",
                                               num_compounds = "integer"
),

prototype(frag_ID = integer(0),
          mz = numeric(0),
          exact_parent_neutral_mass = numeric(0),
          lipid_class = factor(character(0)),
          species = character(0),
          adduct = factor(character(0)),
          adduct_rank = integer(0),
          FA_total_no_C = integer(0),
          FA_total_no_DB = integer(0),
          degree_oxidation = integer(0),
          parent_elem_formula = character(0),
          parent_compound_name = character(0),
          polarity = factor(character(0)),
          num_entries = integer(0),
          num_compounds = integer(0)
))

# create constructor

LOBdbase = function(frag_ID = NULL, mz = NULL, 
                    exact_parent_neutral_mass = NULL, lipid_class = NULL,
                    species = NULL, adduct = NULL, adduct_rank = NULL,
                    FA_total_no_C = NULL, FA_total_no_DB = NULL,
                    degree_oxidation = NULL, parent_elem_formula = NULL,
                    parent_compound_name = NULL, 
                    polarity = NULL, num_entries = NULL,
                    num_compounds = NULL) {
  
  # first, warn user that it's better to use generateLOBdbase or loadLOBdbase
  
  warning("The LOBdbase constructor should only be used for manual ",
          "construction of a LOBdbase. Under almost all circumstances, the ",
          "function generateLOBdbase should be used to create a LOBdbase via ",
          "an in silico simulation. If necessary, loadLOBdbase provides a ",
          "more robust way of recreating a database from a properly ",
          "formatted text file. Use the LOBdbase constructor only as a last ",
          "resort.")
  
  # define a whole number checker
  
  is.wholenumber = function(x, tol = .Machine$
                              double.eps^0.5)  abs(x - round(x)) < tol
  
  # create new (empty) LOBdbase object
  
  object = new("LOBdbase")
  
  # check arguments as best as possible 
  
  # frag_ID
  
  if (!is.null(frag_ID)) {
    
    if (!(class(frag_ID)=="integer")) {
      
      stop("Input 'frag_ID' is not of class 'integer'\n")
      
    }
    
  }
  
  # mz
  
  if (!is.null(mz)) {
    
    if (!(class(mz)=="numeric")) {
      
      stop("Input 'mz' is not of class 'numeric'\n")
      
    }
    
  }
  
  # exact_parent_neutral_mass
  
  if (!is.null(exact_parent_neutral_mass)) {
    
    if (!(class(exact_parent_neutral_mass)=="numeric")) {
      
      stop("Input 'exact_parent_neutral_mass' is not of class 'numeric'\n")
      
    }
    
  }
  
  # lipid_class
  
  if (!is.null(lipid_class)) {
    
    if (!(class(lipid_class)=="factor")) {
      
      stop("Input 'lipid_class' is not of class 'factor'\n")
      
    }
    
  }
  
  # species
  
  if (!is.null(species)) {
    
    if (!(class(species)=="character")) {
      
      stop("Input 'species' is not of class 'character'\n")
      
    }
    
  }
  
  # adduct
  
  if (!is.null(adduct)) {
    
    if (!(class(adduct)=="factor")) {
      
      stop("Input 'adduct' is not of class 'factor'\n")
      
    }
    
  }
  
  # adduct_rank
  
  if (!is.null(adduct_rank)) {
    
    if (!(class(adduct_rank)=="integer")) {
      
      stop("Input 'adduct_rank' is not of class 'integer'\n")
      
    }
    
  }
  
  # FA_total_no_C
  
  if (!is.null(FA_total_no_C)) {
    
    if (!(class(FA_total_no_C)=="integer")) {
      
      stop("Input 'FA_total_no_C' is not of class 'integer'\n")
      
    }
    
  }
  
  # FA_total_no_DB
  
  if (!is.null(FA_total_no_DB)) {
    
    if (!(class(FA_total_no_DB)=="integer")) {
      
      stop("Input 'FA_total_no_DB' is not of class 'integer'\n")
      
    }
    
  }
  
  # degree_oxidation
  
  if (!is.null(degree_oxidation)) {
    
    if (!(class(degree_oxidation)=="integer")) {
      
      stop("Input 'degree_oxidation' is not of class 'integer'\n")
      
    }
    
  }
  
  # parent_elem_formula
  
  if (!is.null(parent_elem_formula)) {
    
    if (!(class(parent_elem_formula)=="character")) {
      
      stop("Input 'parent_elem_formula' is not of class 'character'\n")
      
    }
    
  }
  
  # parent_compound_name
  
  if (!is.null(parent_compound_name)) {
    
    if (!(class(parent_compound_name)=="character")) {
      
      stop("Input 'parent_compound_name' is not of class 'character'\n")
      
    }
    
  }
  
  # polarity
  
  if (!is.null(polarity)) {
    
    if (!(length(polarity)==1 && (polarity %in% c("positive","negative")))) {
      
      stop("Input 'polarity' must be either 'negative' or 'positive'\n")
      
    } else {
      
      polarity = factor(as.character(polarity), c("positive","negative"))
      
    }
    
  }
  
  # num_entries
  
  if (!is.null(num_entries)) {

    if (!(length(num_entries)==1 && is.wholenumber(num_entries))) {
      
      stop("Input 'num_entries' is an integer of length = 1.'\n")
      
    } else {
      
      num_entries = as.integer(num_entries)
      
    }
    
  }
  
  # num_compounds
  
  if (!is.null(num_compounds)) {
    
    if (!(length(num_compounds)==1 && is.wholenumber(num_compounds))) {
      
      stop("Input 'num_compounds' is an integer of length = 1.'\n")
      
    } else {
      
      num_compounds = as.integer(num_compounds)
      
    }
    
  }
  
  # now, check to make sure we all supplied arguments are the same length
  
  # obtain list of values passed to the function, including defaults if they are
  # NULL
  
  argList =  mget(names(formals()),sys.frame(sys.nframe()))
  
  # pull out arguments that contain data
  
  dataArgs = argList[1:(length(argList)-3)]
  
  dataArgs.lengths = lengths(dataArgs)
  
  # now, check for compliance
  
  if (length(unique(dataArgs.lengths))==1) {
    
    if (!is.null(num_entries)) {
      
      if (num_entries!=dataArgs.lengths[1]) {
        
        stop("Lengths of supplied data elements do not match the value ",
             "provided for the metadata field 'num_entries'\n")
        
      }
      
    }
    
  } else {
    
    stop("Supplied data elements are not of equal length.")
    
  }
  
  return(object)

}
  
# set generics and accessors

setMethod("show", "LOBdbase", function(object) {
  
  if (length(frag_ID(object))>0) {
    
    cat("\nA",as.character(polarity(object)),
        "polarity \"LOBdbase\" object.\n\n")
    cat("Contains entries for",num_entries(object), "possible adduct ions of",
        num_compounds(object),"unique parent compounds.\n\n")
    
    cat("Parent lipid types (",length(levels(lipid_class(object))),"):",
        paste(levels(lipid_class(object)), collapse = ", "),"\n")
    cat("IP-DAG classes (",length(as.character(unique(species(object)[
      lipid_class(object)=="IP_DAG"]))),"):",
        paste(as.character(unique(species(object)[
          lipid_class(object)=="IP_DAG"])), collapse = ", "),"\n")
    cat("IP-MAG classes (",length(as.character(unique(species(object)[
      lipid_class(object)=="IP_MAG"]))),"):",
      paste(as.character(unique(species(object)[
        lipid_class(object)=="IP_MAG"])), collapse = ", "),"\n")
    cat("Photosynthetic pigments (",length(as.character(unique(species(object)[
      lipid_class(object)=="pigment"]))),"):",
        paste(as.character(unique(species(object)[
          lipid_class(object)=="pigment"])), collapse = ", "),"\n")
    cat("Adducts (",length(levels(adduct(object))),"):",
        paste(levels(adduct(object)), collapse = ", "),"\n\n")
    cat("m/z range:",paste(min(mz(object), na.rm = TRUE),
                           max(mz(object), na.rm = TRUE), sep = "-"),"\n\n")
    
    cat("Ranges of chemical parameters represented in molecules with acyl",
        "moieties:\n\n")
    cat("Total number of acyl carbon atoms:",
        paste(min(FA_total_no_C(object), na.rm = TRUE),
              max(FA_total_no_C(object), na.rm = TRUE), sep = "-"),"\n")
    cat("Total number of acyl carbon-carbon double bonds:",
        paste(min(FA_total_no_DB(object), na.rm = TRUE),
              max(FA_total_no_DB(object), na.rm = TRUE), sep = "-"),"\n")
    cat("Number of additional oxygen atoms:",
        paste(min(degree_oxidation(object), na.rm = TRUE),
              max(degree_oxidation(object), na.rm = TRUE), sep = "-"),"\n\n")
    
  } else {
    
    cat("\nA \"LOBdbase\" object that appears to be empty.\n\n")
    
  }
  
  memsize = object.size(object)
  cat("Memory usage:", signif(memsize/2^20, 3), "MB\n")
  
})

setGeneric("frag_ID", function(object) standardGeneric("frag_ID"))

setMethod("frag_ID", "LOBdbase", function(object) object@frag_ID)

setGeneric("frag_ID<-", function(object, value) standardGeneric("frag_ID<-"))

setReplaceMethod("frag_ID", "LOBdbase", function(object, value) {
  
  object@frag_ID <- value
  
  object
})

setGeneric("mz", function(object) standardGeneric("mz"))

setMethod("mz", "LOBdbase", function(object) object@mz)

setGeneric("mz<-", function(object, value) standardGeneric("mz<-"))

setReplaceMethod("mz", "LOBdbase", function(object, value) {
  
  object@mz <- value
  
  object
})

setGeneric("exact_parent_neutral_mass", 
           function(object) standardGeneric("exact_parent_neutral_mass"))

setMethod("exact_parent_neutral_mass", "LOBdbase", 
          function(object) object@exact_parent_neutral_mass)

setGeneric(
  "exact_parent_neutral_mass<-", 
  function(object, value) standardGeneric("exact_parent_neutral_mass<-"))

setReplaceMethod("exact_parent_neutral_mass", "LOBdbase", 
                 function(object, value) {
  
  object@exact_parent_neutral_mass <- value
  
  object
})

setGeneric("lipid_class", function(object) standardGeneric("lipid_class"))

setMethod("lipid_class", "LOBdbase", function(object) object@lipid_class)

setMethod("lipid_class", "LOBdbase", 
          function(object) object@lipid_class)

setGeneric(
  "lipid_class<-", 
  function(object, value) standardGeneric("lipid_class<-"))

setReplaceMethod("lipid_class", "LOBdbase", 
                 function(object, value) {
                   
                   object@lipid_class <- value
                   
                   object
                 })

setGeneric("species", function(object) standardGeneric("species"))

setMethod("species", "LOBdbase", function(object) object@species)

setGeneric(
  "species<-", 
  function(object, value) standardGeneric("species<-"))

setReplaceMethod("species", "LOBdbase", 
                 function(object, value) {
                   
                   object@species <- value
                   
                   object
                 })

setGeneric("adduct", function(object) standardGeneric("adduct"))

setMethod("adduct", "LOBdbase", function(object) object@adduct)

setGeneric(
  "adduct<-", 
  function(object, value) standardGeneric("adduct<-"))

setReplaceMethod("adduct", "LOBdbase", 
                 function(object, value) {
                   
                   object@adduct <- value
                   
                   object
                 })

setGeneric("adduct_rank", function(object) standardGeneric("adduct_rank"))

setMethod("adduct_rank", "LOBdbase", function(object) object@adduct_rank)

setGeneric(
  "adduct_rank<-", 
  function(object, value) standardGeneric("adduct_rank<-"))

setReplaceMethod("adduct_rank", "LOBdbase", 
                 function(object, value) {
                   
                   object@adduct_rank <- value
                   
                   object
                 })

setGeneric("FA_total_no_C", function(object) standardGeneric("FA_total_no_C"))

setMethod("FA_total_no_C", "LOBdbase", function(object) object@FA_total_no_C)

setGeneric(
  "FA_total_no_C<-", 
  function(object, value) standardGeneric("FA_total_no_C<-"))

setReplaceMethod("FA_total_no_C", "LOBdbase", 
                 function(object, value) {
                   
                   object@FA_total_no_C <- value
                   
                   object
                 })


setGeneric("FA_total_no_DB", function(object) standardGeneric("FA_total_no_DB"))

setMethod("FA_total_no_DB", "LOBdbase", function(object) object@FA_total_no_DB)

setGeneric(
  "FA_total_no_DB<-", 
  function(object, value) standardGeneric("FA_total_no_DB<-"))

setReplaceMethod("FA_total_no_DB", "LOBdbase", 
                 function(object, value) {
                   
                   object@FA_total_no_DB <- value
                   
                   object
                 })

setGeneric("degree_oxidation", 
           function(object) standardGeneric("degree_oxidation"))

setMethod("degree_oxidation", "LOBdbase", 
          function(object) object@degree_oxidation)

setGeneric(
  "degree_oxidation<-", 
  function(object, value) standardGeneric("degree_oxidation<-"))

setReplaceMethod("degree_oxidation", "LOBdbase", 
                 function(object, value) {
                   
                   object@degree_oxidation <- value
                   
                   object
                 })

setGeneric("parent_elem_formula", 
           function(object) standardGeneric("parent_elem_formula"))

setMethod("parent_elem_formula", "LOBdbase", 
          function(object) object@parent_elem_formula)

setGeneric(
  "parent_elem_formula<-", 
  function(object, value) standardGeneric("parent_elem_formula<-"))

setReplaceMethod("parent_elem_formula", "LOBdbase", 
                 function(object, value) {
                   
                   object@parent_elem_formula <- value
                   
                   object
                 })

setGeneric("parent_compound_name", 
           function(object) standardGeneric("parent_compound_name"))

setMethod("parent_compound_name", "LOBdbase", 
          function(object) object@parent_compound_name)

setGeneric(
  "parent_compound_name<-", 
  function(object, value) standardGeneric("parent_compound_name<-"))

setReplaceMethod("parent_compound_name", "LOBdbase", 
                 function(object, value) {
                   
                   object@parent_compound_name <- value
                   
                   object
                 })

# setGeneric("polarity", function(object) standardGeneric("polarity"))

setMethod("polarity", "LOBdbase", function(object) object@polarity)

# setGeneric(
#   "polarity<-", 
#   function(object, value) standardGeneric("polarity<-"))

setReplaceMethod("polarity", "LOBdbase", 
                 function(object, value) {
                   
                   object@polarity <- value
                   
                   object
                 })

setGeneric("num_entries", function(object) standardGeneric("num_entries"))

setMethod("num_entries", "LOBdbase", function(object) object@num_entries)

setGeneric(
  "num_entries<-", 
  function(object, value) standardGeneric("num_entries<-"))

setReplaceMethod("num_entries", "LOBdbase", 
                 function(object, value) {
                   
                   object@num_entries <- value
                   
                   object
                 })

setGeneric("num_compounds", function(object) standardGeneric("num_compounds"))

setMethod("num_compounds", "LOBdbase", function(object) object@num_compounds)

setGeneric(
  "num_compounds<-", 
  function(object, value) standardGeneric("num_compounds<-"))

setReplaceMethod("num_compounds", "LOBdbase", 
                 function(object, value) {
                   
                   object@num_compounds <- value
                   
                   object
                 })

## to generate .Rd file:
# library(methods)
# promptClass(LOBdbase)