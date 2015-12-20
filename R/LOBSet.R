################ Set classes, methods for LOBSet #############

# create a class "LOBSet" for the results of LOBSTAHS screening & compound assignments

LOBSet = setClass("LOBSet",
                  representation(LOBpeakdata = "data.frame",
                                 iso.C3r = "list",
                                 iso.C3f = "list",
                                 iso.C3c = "list",
                                 LOBscreen.diagnostics = "data.frame",
                                 LOBisoID.diagnostics = "data.frame",
                                 LOBscreen.settings = "list",
                                 parent.xsAnnotate = "xsAnnotate"
                  ),
                  
                  prototype(LOBpeakdata = data.frame(),
                            iso.C3r = list(),
                            iso.C3f = list(),
                            iso.C3c = list(),
                            LOBscreen.diagnostics = data.frame(),
                            LOBisoID.diagnostics = data.frame(),
                            LOBscreen.settings = list(),
                            parent.xsAnnotate = new("xsAnnotate")
                  ))

setMethod("show", "LOBSet", function(object) {
  
  cat("A",as.character(object@parent.xsAnnotate@polarity),"polarity \"LOBSet\" containing LC-MS peak data. Compound assignments and adduct ion hierarchy screening annotations applied to",length(sampnames(object@parent.xsAnnotate@xcmsSet)),"samples using the \"LOBSTAHS\" package.\n\n")
  cat("Individual peaks:",object@LOBscreen.diagnostics$peaks[6],"\n")
  cat("Peak groups:",object@LOBscreen.diagnostics$peakgroups[6],"\n")
  cat("Compound assignments:",object@LOBscreen.diagnostics$parent_compounds[6],"\n")
  cat("m/z range:",paste(min(object@LOBpeakdata$peakgroup.mz, na.rm = TRUE),max(object@LOBpeakdata$peakgroup.mz, na.rm = TRUE), sep = "-"),"\n\n")
  
  cat("Possible regisomers:",paste(object@LOBisoID.diagnostics$peakgroups[1],"\n"))
  cat("Possible structural functional isomers:",paste(object@LOBisoID.diagnostics$peakgroups[2],"\n"))
  cat("Isobars indistinguishable within ppm matching tolerance:",paste(object@LOBisoID.diagnostics$peakgroups[3],"\n\n"))
  
  cat("Restrictions applied prior to conducting adduct ion hierarchy screening:",paste(c("remove.iso","rt.restrict","exclude.oddFA")[unlist(object@LOBscreen.settings[c("remove.iso","rt.restrict","exclude.oddFA")])], collapse = ", "),"\n\n")
  
  cat("Match tolerance used for database assignments:",object@LOBscreen.settings$match.ppm,"ppm\n\n")
  #   cat("Ranges of chemical parameters represented in molecules other than pigments:\n\n")
  #   cat("Total number of acyl carbon atoms:",paste(min(object@LOBpeakdata$FA_total_no_C, na.rm = TRUE),max(object@LOBpeakdata$FA_total_no_C, na.rm = TRUE), sep = "-"),"\n")
  #   cat("Total number of acyl carbon-carbon double bonds:",paste(min(object@LOBpeakdata$FA_total_no_DB, na.rm = TRUE),max(object@LOBpeakdata$FA_total_no_DB, na.rm = TRUE), sep = "-"),"\n")
  #   cat("Number of additional oxygen atoms:",paste(min(object@LOBpeakdata$degree_oxidation, na.rm = TRUE),max(object@LOBpeakdata$degree_oxidation, na.rm = TRUE), sep = "-"),"\n\n")
  
  cat("The CAMERA \"xsAnnotate\" object used for screening was retained in slot @parent.xsAnnotate.\n\n")
  
  memsize = object.size(object)
  cat("Memory usage:", signif(memsize/2^20, 3), "MB\n")
  
})