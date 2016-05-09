################ Set classes, methods for LOBdbase #############

# constructor

LOBdbase = setClass("LOBdbase", representation(frag_ID = "integer",
                                               mz = "numeric",
                                               exact_parent_neutral_mass = "numeric",
                                               lipid_class = "factor",
                                               species = "character",
                                               adduct = "factor",
                                               adduct_rank = "integer",
                                               FA_total_no_C = "integer",
                                               FA_total_no_DB = "integer",
                                               degree_oxidation = "integer",
                                               parent_elem_formula = "character",
                                               parent_compound_name = "character",
                                               polarity = "factor",
                                               num.entries = "integer",
                                               num.compounds = "integer"
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
          num.entries = integer(0),
          num.compounds = integer(0)
))

# set generics and accessors

setMethod("show", "LOBdbase", function(object) {
  
  cat("\nA",as.character(object@polarity),"polarity \"LOBdbase\" object.\n\n")
  cat("Contains entries for",object@num.entries, "possible adduct ions of",object@num.compounds,"unique parent compounds.\n\n")
  
  cat("Parent lipid types (",length(levels(object@lipid_class)),"):",paste(levels(object@lipid_class), collapse = ", "),"\n")
  cat("IP-DAG classes (",length(as.character(unique(object@species[object@lipid_class=="IP_DAG"]))),"):",paste(as.character(unique(object@species[object@lipid_class=="IP_DAG"])), collapse = ", "),"\n")
  cat("Pigments (",length(as.character(unique(object@species[object@lipid_class=="pigment"]))),"):",paste(as.character(unique(object@species[object@lipid_class=="pigment"])), collapse = ", "),"\n")
  cat("Adducts (",length(levels(object@adduct)),"):",paste(levels(object@adduct), collapse = ", "),"\n\n")
  cat("m/z range:",paste(min(object@mz, na.rm = TRUE),max(object@mz, na.rm = TRUE), sep = "-"),"\n\n")
  
  cat("Ranges of chemical parameters represented in molecules other than pigments:\n\n")
  cat("Total number of acyl carbon atoms:",paste(min(object@FA_total_no_C, na.rm = TRUE),max(object@FA_total_no_C, na.rm = TRUE), sep = "-"),"\n")
  cat("Total number of acyl carbon-carbon double bonds:",paste(min(object@FA_total_no_DB, na.rm = TRUE),max(object@FA_total_no_DB, na.rm = TRUE), sep = "-"),"\n")
  cat("Number of additional oxygen atoms:",paste(min(object@degree_oxidation, na.rm = TRUE),max(object@degree_oxidation, na.rm = TRUE), sep = "-"),"\n\n")
  
  memsize = object.size(object)
  cat("Memory usage:", signif(memsize/2^20, 3), "MB\n")
  
})

# setMethod("[", signature("LOBdbase", "ANY", "ANY", "missing"),
#           function(x, i, j, ..., drop=TRUE) {
#             
#             .frag_ID = x@frag_ID[i]
#             .mz = x@mz[i]
#             .exact_parent_neutral_mass = x@exact_parent_neutral_mass[i]
#             .lipid_class = x@lipid_class[i]
#             .species = x@species[i]
#             .adduct = x@adduct[i]
#             .adduct_rank = x@adduct_rank[i]
#             .FA_total_no_C = x@FA_total_no_C[i]
#             .FA_total_no_DB = x@FA_total_no_DB[i]
#             .degree_oxidation = x@degree_oxidation[i]
#             .parent_elem_formula = x@parent_elem_formula[i]
#             .parent_compound_name = x@parent_compound_name[i]
#             .polarity = x@polarity
#             
#             LOBdbase(frag_ID = .frag_ID,
#                      mz = .mz,
#                      exact_parent_neutral_mass = .exact_parent_neutral_mass,
#                      lipid_class = .lipid_class,
#                      species = .species,
#                      adduct = .adduct,
#                      adduct_rank = .adduct_rank,
#                      FA_total_no_C = .FA_total_no_C,
#                      FA_total_no_DB = .FA_total_no_DB,
#                      degree_oxidation = .degree_oxidation,
#                      parent_elem_formula = .parent_elem_formula,
#                      parent_compound_name = .parent_compound_name,
#                      polarity = .polarity)
#             
#           })

setGeneric("frag_ID", function(object) standardGeneric("frag_ID"))

setMethod("frag_ID", "LOBdbase", function(object) object@frag_ID)

setGeneric("mz", function(object) standardGeneric("mz"))

setMethod("mz", "LOBdbase", function(object) object@mz)

setGeneric("exact_parent_neutral_mass", function(object) standardGeneric("exact_parent_neutral_mass"))

setMethod("exact_parent_neutral_mass", "LOBdbase", function(object) object@exact_parent_neutral_mass)

setGeneric("lipid_class", function(object) standardGeneric("lipid_class"))

setMethod("lipid_class", "LOBdbase", function(object) object@lipid_class)

setGeneric("species", function(object) standardGeneric("species"))

setMethod("species", "LOBdbase", function(object) object@species)

setGeneric("adduct", function(object) standardGeneric("adduct"))

setMethod("adduct", "LOBdbase", function(object) object@adduct)

setGeneric("adduct_rank", function(object) standardGeneric("adduct_rank"))

setMethod("adduct_rank", "LOBdbase", function(object) object@adduct_rank)

setGeneric("FA_total_no_C", function(object) standardGeneric("FA_total_no_C"))

setMethod("FA_total_no_C", "LOBdbase", function(object) object@FA_total_no_C)

setGeneric("FA_total_no_DB", function(object) standardGeneric("FA_total_no_DB"))

setMethod("FA_total_no_DB", "LOBdbase", function(object) object@FA_total_no_DB)

setGeneric("degree_oxidation", function(object) standardGeneric("degree_oxidation"))

setMethod("degree_oxidation", "LOBdbase", function(object) object@degree_oxidation)

setGeneric("parent_elem_formula", function(object) standardGeneric("parent_elem_formula"))

setMethod("parent_elem_formula", "LOBdbase", function(object) object@parent_elem_formula)

setGeneric("parent_compound_name", function(object) standardGeneric("parent_compound_name"))

setMethod("parent_compound_name", "LOBdbase", function(object) object@parent_compound_name)

setGeneric("polarity", function(object) standardGeneric("polarity"))

setMethod("polarity", "LOBdbase", function(object) object@polarity)

setGeneric("num.entries", function(object) standardGeneric("num.entries"))

setMethod("num.entries", "LOBdbase", function(object) object@num.entries)

setGeneric("num.compounds", function(object) standardGeneric("num.compounds"))

setMethod("num.compounds", "LOBdbase", function(object) object@num.compounds)

## to generate .Rd file:
# library(methods)
# promptClass(LOBdbase)