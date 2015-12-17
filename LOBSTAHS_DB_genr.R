################ Set classes, methods #############

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
  
  memsize <- object.size(object)
  cat("Memory usage:", signif(memsize/2^20, 3), "MB\n")
  
})

setMethod("[", "LOBdbase",
          function(x,i,drop="missing") {
            
            .frag_ID = x@frag_ID[i]
            .mz = x@mz[i]
            .exact_parent_neutral_mass = x@exact_parent_neutral_mass[i]
            .lipid_class = x@lipid_class[i]
            .species = x@species[i]
            .adduct = x@adduct[i]
            .adduct_rank = x@adduct_rank[i]
            .FA_total_no_C = x@FA_total_no_C[i]
            .FA_total_no_DB = x@FA_total_no_DB[i]
            .degree_oxidation = x@degree_oxidation[i]
            .parent_elem_formula = x@parent_elem_formula[i]
            .parent_compound_name = x@parent_compound_name[i]
            .polarity = x@polarity
            
            LOBdbase(frag_ID = .frag_ID,
                     mz = .mz,
                     exact_parent_neutral_mass = .exact_parent_neutral_mass,
                     lipid_class = .lipid_class,
                     species = .species,
                     adduct = .adduct,
                     adduct_rank = .adduct_rank,
                     FA_total_no_C = .FA_total_no_C,
                     FA_total_no_DB = .FA_total_no_DB,
                     degree_oxidation = .degree_oxidation,
                     parent_elem_formula = .parent_elem_formula,
                     parent_compound_name = .parent_compound_name,
                     polarity = .polarity)
            
          })

################ Wrapper function #############

# generateLOBdbase: wrapper function for lipid-ox-lipid-oxylipin database generation

generateLOBdbase = function(polarity = c("positive","negative"), gen.csv = FALSE, output.dir = "LOBDB", component.defs = NULL, AIH.defs = NULL, acyl.ranges = NULL, oxy.ranges = NULL) {
  
  polarity = match.arg(polarity, several.ok = TRUE)
  
  # determine whether user has specified external source files for any input parameters, get correct file paths, and warn user of consequences of using improperly formatted or incomplete external data
  # perhaps later, can write a function to check the formatting of the data in each file; not a high priority right now since we presume user would have followed the layout of the default tables if creating their own
  
  if (is.null(component.defs)) { # user didn't specify external component definitions table, use defaults
    
    componentTable.loc = "dependencies/LOBSTAHS_basic_component_matrix.csv"
    
  } else { # user specified something
    
    componentTable.loc = component.defs
    
    # let user know he/she provided external input, and the possible consequences
    
    cat("User specified external source for basic component composition matrix.\n")
    cat("Ensure .csv file is properly formatted and sufficient entries exist in\n")
    cat("other external files to support any additional adducts, lipid classes, or molecules.\n")
    
  }
  
  if (is.null(AIH.defs)) { # user didn't specify external AIH table, use defaults
    
    AIHtable.loc = "dependencies/LOBSTAHS_adduct_ion_hierarchies.csv"
    
  } else { # user specified something
    
    AIHtable.loc = AIH.defs
    
    # let user know he/she provided external input, and the possible consequences
    
    cat("User specified external source for adduct ion hierarchy matrix.\n")
    cat("Ensure .csv file is properly formatted and sufficient entries exist in\n")
    cat("other external files to support any additional adducts, lipid classes, or molecules.\n")
    
  }
  
  if (is.null(acyl.ranges)) { # user didn't specify external in silico acyl property range table, use defaults
    
    acylRanges.loc = "dependencies/LOBSTAHS_acyl_prop_ranges.csv"
    
  } else { # user specified something
    
    acylRanges.loc = acyl.ranges
    
    # let user know he/she provided external input, and the possible consequences
    
    cat("User specified external source for in silico simulation acyl property ranges.\n")
    cat("Ensure .csv file is properly formatted and sufficient entries exist in\n")
    cat("other external files to support any additional adducts, lipid classes, or molecules.\n")
    
  }
  
  if (is.null(oxy.ranges)) { # user didn't specify external AIH table, use defaults
    
    oxyRanges.loc = "dependencies/LOBSTAHS_addl_oxy_ranges.csv"
    
  } else { # user specified something
    
    oxyRanges.loc = oxy.ranges
    
    # let user know he/she provided external input, and the possible consequences
    
    cat("User specified external source for numbers of additional oxygen atoms to be considered.\n")
    cat("Ensure .csv file is properly formatted and sufficient entries exist in\n")
    cat("other external files to support any additional adducts, lipid classes, or molecules.\n")
    
  }
  
  # load in silico simulation parameters using helper functions
  
  masses = calcComponentMasses(componentTable.loc)
  ranges = loadSimRanges(acylRanges.loc,oxyRanges.loc)
  adductHierarchies = loadAIH(AIHtable.loc)
  
  # run simulation(s)
  
  sapply(polarity, runSim, acylRanges = ranges$acylC_DB, oxyRanges = ranges$addl_oxy, adductHierarchies = adductHierarchies, baseComponent.masses = masses$baseComponents, adduct.masses = masses$adducts, gen.csv = gen.csv, output.dir = output.dir, simplify = FALSE,USE.NAMES = TRUE)
  
}

################ Helper functions (some will be private) #############

# defineElemExactMasses: creates table of exact masses of the basic chemical building blocks used in DB generation (elements, monatomic ions, and 2 species that are adduct components)

defineElemExactMasses = function() { 
  
  # specify exact masses of some necessary chemical species
  
  # values from de Laeter et al., 2003, "Atomic weights of the elements," Pure Appl. Chem. 75(6): 683–800; C. Amsler et al., 2008, "Review of particle physics," Physics Letters B667: 1
  
  m_C = 12
  m_H = 1.00782504
  m_H_plus = 1.007276467
  m_N = 14.00307401
  m_O = 15.99491462
  m_P = 30.97376151
  m_S = 31.97207069
  m_Na = 22.98977
  m_Cl = 34.96885271
  m_K = 38.96370686
  m_e_minus = 0.00054858
  m_Mg = 23.985045
  
  # calculate exact masses of acetonitrile and acetate using data we just specified
  
  m_ACN = 2*m_C + 3*m_H + 1*m_N
  m_Ac_minus = 2*m_C + 3*m_H + 2*m_O
  
  # create an exact-mass lookup table and put it in alphabetical order
  
  exact.masses = c(m_C,m_H,m_H_plus,m_N,m_O,m_P,m_S,m_Na,m_Cl,m_K,m_e_minus,m_Mg,m_ACN,m_Ac_minus)
  exact.masses = as.table(exact.masses)
  
  names(exact.masses) = c("m_C","m_H","m_H_plus","m_N","m_O","m_P","m_S","m_Na","m_Cl","m_K","m_e_minus","m_Mg","m_ACN","m_Ac_minus")
  
  exact.masses = exact.masses[order(names(exact.masses))]
  
  return(exact.masses)
  
}

# calcComponentMasses: calculates exact masses of components specified in the component table, using the given compositions and the exact masses produced by defineElemExactMasses()
# creates two tables: (1) baseComponents, with compositions and exact masses of basic, non-adduct components and (2) adducts, containing exact masses of adducts

calcComponentMasses = function(componentTableLoc) { # input should be the component composition file location obtained from getTableLocs, or other
  
  # get exact masses using defineElemExactMasses()
  
  exact.masses = defineElemExactMasses()
  
  # load in component composition matrix
  
  componentCompTable.raw = read.table(componentTableLoc, sep=",", header = TRUE, row.names = 1)
  
  # put columns in alphabetical order
  
  componentCompTable = componentCompTable.raw[order(colnames(componentCompTable.raw))]
  
  # calculate exact masses of basic components and extract into a few separate tables
  # note: this will calculate full exact masses of pigments and DNP-PE
  
  # check to make sure we have same number of elemental building blocks in our exact.masses table and along the second dimension of the composition table
  
  if (length(exact.masses)!=(ncol(componentCompTable)-1)) {
    
    stop("Different number of chemical building blocks in the component composition matrix and in the onboard list of exact masses. Check your composition matrix carefully. Aborting...") # stop script if this is the case
    
  }
  
  # assuming same no. of building blocks, calculate exact masses & store as additional column in componentCompTable
  
  componentCompTable[,ncol(componentCompTable)+1] = apply(componentCompTable[,1:ncol(componentCompTable)-1], 1, function(x) sum(x*exact.masses,na.rm = TRUE))
  colnames(componentCompTable)[ncol(componentCompTable)] = c("Exact_mass")
  
  # extract masses of adducts into separate table (we'll need these later); create basecomponent.masses by removing adduct data from componentCompTable
  
  adduct.masses = componentCompTable[componentCompTable$Species_class %in% c("adduct_neg","adduct_pos"),]
  
  baseComponent.masses = componentCompTable[!componentCompTable$Species_class  %in% c("adduct_neg","adduct_pos"),]
  
  # remove ACN and Ac- from baseComponent.masses
  
  baseComponent.masses = baseComponent.masses[-match(c("ACN","Ac-"),rownames(baseComponent.masses)),]
  
  # return adduct.masses, baseComponent.masses as list
  
  list(adducts=adduct.masses,baseComponents=baseComponent.masses)
  
}

# loadSimRanges: loads the necessary in silico simulation range data from the acylRanges and oxyRanges file locations. creates two data frames: (1) acylC_DB and (2) addl_oxy

loadSimRanges = function(acylRangeTableLoc,oxyRangeTableLoc) { # inputs should be the two file locations obtained from getTableLocs, or other
 
  # load in ranges of total acyl C atoms and double bonds to be considered during simulations
  
  acylRanges = read.table(acylRangeTableLoc, sep=",", skip = 1, header = TRUE)
  
  # load in ranges of additional oxygen atoms to be considered during simulations
  
  oxyRanges = read.table(oxyRangeTableLoc, sep=",", skip = 1, header = TRUE)
 
  list(acylC_DB=acylRanges,addl_oxy=oxyRanges)
  
}

# loadAIH: loads the adduct ion hierarchy data from the AIHfile location. creates one data frame called adductHierarchies

loadAIH = function(AIHTableLoc) { # input should be file location obtained from getTableLocs, or other
  
  adductHierarchies = read.table(AIHTableLoc, sep=",", skip = 1, header = TRUE)
  
  # for compatibility, also assign values in "Adduct" as row names
  
  row.names(adductHierarchies) = adductHierarchies$Adduct
  
  return(adductHierarchies)
  
}

# calcNumCombs: calculates the number of parent compounds and adduct ions for which masses are to be generated in a given ion mode, based on the user-specified ranges of lipid classes and chemical properties

calcNumCombs = function(polarity, acylRanges, oxyRanges, adductHierarchies, baseComponent.masses, adduct.masses) {
  
  # extract adduct hierarchies for this mode, compare the number to those in adduct.masses to ensure the mass of each of them is defined in the adducts mass table
  
  AIHs.thismode = adductHierarchies[adductHierarchies$Adduct_ion_mode==polarity,]
  
  if (!all(is.element(AIHs.thismode$Adduct,rownames(adduct.masses)))) {
    
    stop("Not all adduct ions given in the hierarchy table for this mode are defined in the component composition table.")
    
  }
  
  # now, proceed with obtaining number of combinations
  
  # initialize counters
  numAddIons = 0
  numCompounds = 0
  
  for (i in 1:nrow(baseComponent.masses)) {
    
    # retrieve necessary data for this class
    
    this.class = baseComponent.masses$Species_class[i]
    
    if (this.class %in% c("IP_DAG","FFA","TAG","PUA")) {
      
      this.oxymin = as.numeric(oxyRanges[grep(paste0(as.character(baseComponent.masses$Species_class[i]),"_min"),colnames(oxyRanges))])
      this.oxymax = as.numeric(oxyRanges[grep(paste0(as.character(baseComponent.masses$Species_class[i]),"_max"),colnames(oxyRanges))])
      this.C_DBmindata = acylRanges[,grep(paste0(as.character(baseComponent.masses$Species_class[i]),"_min"),colnames(acylRanges))]
      this.C_DBmaxdata = acylRanges[,grep(paste0(as.character(baseComponent.masses$Species_class[i]),"_max"),colnames(acylRanges))]
      
      num.adducts = sum(!is.na(AIHs.thismode[,colnames(AIHs.thismode)[colnames(AIHs.thismode)==rownames(baseComponent.masses)[i]]]))
      
      if (num.adducts>0) {
        
        num_compounds.this_species = (this.oxymax-this.oxymin+1)*sum(this.C_DBmaxdata-this.C_DBmindata+1,na.rm = TRUE)
      } else {
        
        num_compounds.this_species = 0
        
      }
      
      num_ions.this_species = num.adducts*num_compounds.this_species
      
    } else if (this.class %in% c("DNPPE","pigment")) {
      
      if (this.class=="pigment") {
        
        num.adducts = sum(!is.na(AIHs.thismode[,colnames(AIHs.thismode)[colnames(AIHs.thismode)==baseComponent.masses$Species_class[i]]]))
        
      } else if (this.class=="DNPPE") {
        
        num.adducts = sum(!is.na(AIHs.thismode[,colnames(AIHs.thismode)[colnames(AIHs.thismode)==rownames(baseComponent.masses)[i]]]))
        
      }
      
      if (num.adducts>0) {
        
        num_compounds.this_species = 1 # because we considered DNPPE and each pigment individually
        
      } else {
        
        num_compounds.this_species = 0
        
      }
      
      num_ions.this_species = num.adducts
      
    }
    
    numCompounds = numCompounds + num_compounds.this_species
    numAddIons = numAddIons + num_ions.this_species
    
  }
  
  return(list(numCompounds = as.integer(numCompounds), numAddIons = as.integer(numAddIons)))
  
}

# genTimeStamp: generates a timestamp string based on the current system time

genTimeStamp = function () {
  
  output_DTG = format(Sys.time(), "%Y-%m-%dT%X%z") # return current time in a good format
  output_DTG = gsub(" ", "_", output_DTG) # replace any spaces
  output_DTG = gsub(":", "-", output_DTG) # replaces any colons with dashes (Mac compatibility)
  
}

# exportDBtoCSV: writes a LOBdbase object to file

exportDBtoCSV = function(output.dir = NULL, LOBdbase) {
  
  if (is.null(output.dir)) {
    
    export.filepath = NULL
    
  } else {
    
    if (!file.exists(output.dir)) {
      
      dir.create(file.path(output.dir))
      
    }
    
    export.filepath = paste0(output.dir,"/")
    
  }
  
  
  output_DTG = genTimeStamp()
  
  cat("Exporting .csv file containing",as.character(LOBdbase@polarity),"mode simulation output...\n")
  
  fname = paste0(export.filepath,"LOBSTAHS_lipid-oxy_DB_",strtrim(as.character(LOBdbase@polarity),3),"_",output_DTG,".csv")
  
  exportmat = data.frame(LOBdbase@frag_ID,
                         LOBdbase@mz,
                         LOBdbase@exact_parent_neutral_mass,
                         as.character(LOBdbase@lipid_class),
                         as.character(LOBdbase@species),
                         as.character(LOBdbase@adduct),
                         as.character(LOBdbase@adduct_rank),
                         LOBdbase@FA_total_no_C,
                         LOBdbase@FA_total_no_DB,
                         LOBdbase@degree_oxidation,
                         LOBdbase@parent_elem_formula,
                         LOBdbase@parent_compound_name,
                         stringsAsFactors = FALSE)
  
  colnames(exportmat) = c("frag_ID","mz","exact_parent_neutral_mass","lipid_class","species","adduct","adduct_rank","FA_total_no_C","FA_total_no_DB","degree_oxidation","parent_elem_formula","parent_compound_name")
  
  write.csv(exportmat, fname)
  
  cat(as.character(LOBdbase@polarity),"mode database exported to:",fname,"\n")
  
}

# runSim: runs the in silico simulation for a given ion mode, returns a 

runSim = function(polarity, acylRanges, oxyRanges, adductHierarchies, baseComponent.masses, adduct.masses, output.dir, gen.csv) {
  
  # calculate number of combinations for which data are to be calculated in this mode
  
  numCombs = calcNumCombs(polarity, acylRanges, oxyRanges, adductHierarchies, baseComponent.masses, adduct.masses)
  
  # provide feedback to user
  
  cat("\nPerforming in silico simulation to generate data for",polarity,"mode species...\n\n")
  
  cat("Based on user settings, LOBSTAHS will generate exact masses and adduct hierarchy data for",numCombs$numAddIons,"possible",polarity,"mode adduct ions, representing",numCombs$numCompounds,"unique parent compounds.\n\n")
  
  # preallocate matrix for results
  
  sim_results = matrix(data = NA, nrow = numCombs$numAddIons, ncol = 12)
  sim_results = as.data.frame(sim_results)
  
  # extract adduct hierarchies for this ion mode
  
  AIHs.thismode = adductHierarchies[adductHierarchies$Adduct_ion_mode==polarity,]
  
  # get exact masses using defineElemExactMasses()
  
  exact.masses = defineElemExactMasses()
  
  # now, perform simulation by species
  
  ins.row = 1 # define variable to keep track of insertion point in results table
  
  for (i in 1:nrow(baseComponent.masses)) {
    
    # retrieve, store this.lipid_class and this.species
    this.lipid_class = as.character(baseComponent.masses$Species_class[i])
    this.species = rownames(baseComponent.masses)[i]
    
    # provide sensible feedback to user
    
    if (this.lipid_class==this.species) {
      
      cat("Calculating data for lipid class:",this.species,"...\n")
      
    } else if (this.lipid_class=="pigment") {
      
      cat("Calculating data for",this.lipid_class,":",this.species,"...\n")
      
    } else if (this.lipid_class=="IP_DAG") {
      
      cat("Calculating data for",this.lipid_class,"lipid class:",this.species,"...\n")
      
    }
    
    # first, need to define an "adduct lookup class" since pigment adduct hierarchy data are currently defined for all pigments, whereas hierarchy data for other classes are class specific
    
    if (baseComponent.masses$Species_class[i]==c("pigment")) { # it's a pigment
      
      adduct.lookup.class = baseComponent.masses$Species_class[i]
      
    } else {
      
      adduct.lookup.class = rownames(baseComponent.masses)[i]
      
    }
    
    # now, check to see whether hierarchy data exists for this species or lipid class in this ionization mode
    
    if (sum(AIHs.thismode[,colnames(AIHs.thismode)[colnames(AIHs.thismode)==adduct.lookup.class]],na.rm = TRUE)>=1) { # hierarchy data does exist for this species or class in this mode; proceed
      
      # get number of adducts formed by this species/class
      
      adducts.thisclass = rownames(AIHs.thismode)[!is.na(AIHs.thismode[,as.character(adduct.lookup.class)])]
      
      # check to see whether this is a pigment or DNP-PE, since we will treat these differently in the simulation (they're the only species for which we don't consider ranges in # of acyl C, double bonds, etc.)
      
      if (adduct.lookup.class %in% c("pigment","DNPPE")) { # this element is a pigment or DNPPE; we only need to drill down to the adduct level
        
        for (j in 1:length(adducts.thisclass)) { # cycle thru adducts
          
          this.mz = adduct.masses$Exact_mass[rownames(adduct.masses) %in% adducts.thisclass[j]]+baseComponent.masses$Exact_mass[i] # calculate mz for this adduct ion of this species
          
          # ascertain parent compound elemental formula
          
          these.base_elements = subset(baseComponent.masses[i,], select=-c(Species_class,Exact_mass))          
          this.parent_formula = paste0(apply(rbind(colnames(these.base_elements)[these.base_elements>=1],as.numeric(these.base_elements[these.base_elements>=1])),2,paste,collapse=""),collapse="")
          
          # get name
          
          this.parent_compound_name = this.species
          
          # now, record data for this adduct ion
          
          this.frag_ID = ins.row
          this.parent_exactneutralmass = baseComponent.masses$Exact_mass[i]
          this.degree_oxidation = 0
          this.adduct = adducts.thisclass[j]
          this.adduct_rank = AIHs.thismode[rownames(AIHs.thismode)==this.adduct,colnames(AIHs.thismode)==adduct.lookup.class]
          this.FA_total_no_C = NA
          this.FA_total_no_DB = NA
          
          # integers
          sim_results[ins.row,c(1,7,8,9,10)] = c(this.frag_ID,this.adduct_rank,this.FA_total_no_C,this.FA_total_no_DB,this.degree_oxidation)
          
          # doubles
          sim_results[ins.row,c(2,3)] = c(this.mz,this.parent_exactneutralmass)
          
          # text fields
          sim_results[ins.row,c(4,5,6,11,12)] = c(this.lipid_class,this.species,this.adduct,this.parent_formula,this.parent_compound_name)
          
          ins.row = ins.row + 1 # advance our insertion point
          
        }
        
        rm(j)
        
      } else { # this species is not a pigment, or DNPPE --> requires more involved simulation
        
        # retrieve, store "base" exact mass for this lipid class
        
        this.base_exactmass = baseComponent.masses$Exact_mass[i]
        
        # load acyl C - double bond distributions from sim.ranges; oxidation state parameters from user-specified oxy_range variable; set variable for number of carboxyl groups (for calculation IP-DAG, FFA, TAG, PUA masses)
        
        these.sim.ranges = acylRanges[!is.na(acylRanges[,grep(paste0(as.character(baseComponent.masses$Species_class[i]),"_min"),colnames(acylRanges))]),c("FA_total_no_C",paste0(as.character(baseComponent.masses$Species_class[i]),"_min"),paste0(as.character(baseComponent.masses$Species_class[i]),"_max"))]
        this.oxymin = as.numeric(oxyRanges[grep(paste0(as.character(baseComponent.masses$Species_class[i]),"_min"),colnames(oxyRanges))])
        this.oxymax = as.numeric(oxyRanges[grep(paste0(as.character(baseComponent.masses$Species_class[i]),"_max"),colnames(oxyRanges))])
        these.oxystates = this.oxymin:this.oxymax
        
        num.carboxyl = switch(
          as.character(baseComponent.masses$Species_class[i]),
          IP_DAG=2,
          FFA=1,
          PUA=0,
          TAG=3)
        
        for (j in 1:nrow(these.sim.ranges)) { # cycle thru number of allowable acyl C atoms
          
          # retrieve, store this.FA_total_no_C
          this.FA_total_no_C = these.sim.ranges[j,1]
          
          # generate vector of allowable acyl C=C double bonds for this FA_total_no_C
          these.allowable.DBs = these.sim.ranges[j,2]:these.sim.ranges[j,3]
          
          for (k in 1:length(these.allowable.DBs)) { # cycle thru allowable no. of double bonds
            
            # retrieve, store this.FA_total_no_DB
            this.FA_total_no_DB = these.allowable.DBs[k]
            
            for (l in 1:length(these.oxystates)) { # cycle thru allowable additional oxygen atoms
              
              # retrieve, store this.degree_oxidation
              this.degree_oxidation = these.oxystates[l]
              
              for (m in 1:length(adducts.thisclass)) { # cycle thru adducts
                
                # retrieve, store this.adduct, this.adduct_exact_mass, this.adduct_rank
                this.adduct_exact_mass = adduct.masses$Exact_mass[rownames(adduct.masses) %in% adducts.thisclass[m]]
                this.adduct = adducts.thisclass[m]
                this.adduct_rank = AIHs.thismode[rownames(AIHs.thismode)==this.adduct,colnames(AIHs.thismode)==adduct.lookup.class]
                
                # now, finally, can assemble and record data for this combination
                
                # calculate exact neutral mass of parent compound that reflects this combination of properties, and mz for this adduct ion of that compound 
                
                this.parent_exactneutralmass = this.base_exactmass+
                  exact.masses[c("m_C")]*this.FA_total_no_C+
                  exact.masses[c("m_H")]*(this.FA_total_no_C*2)-
                  exact.masses[c("m_H")]*num.carboxyl+
                  exact.masses[c("m_O")]*this.degree_oxidation-
                  exact.masses[c("m_H")]*(this.FA_total_no_DB*2)
                
                this.mz = this.parent_exactneutralmass+this.adduct_exact_mass
                
                # ascertain parent compound elemental formula
                
                these.base_elements = subset(baseComponent.masses[i,], select=-c(Species_class,Exact_mass)) 
                these.base_elements["C"] = these.base_elements["C"]+this.FA_total_no_C
                these.base_elements["H"] = these.base_elements["H"]+this.FA_total_no_C*2-num.carboxyl-this.FA_total_no_DB*2
                these.base_elements["O"] = these.base_elements["O"]+this.degree_oxidation
                this.parent_formula = paste0(apply(rbind(colnames(these.base_elements)[these.base_elements>=1],as.numeric(these.base_elements[these.base_elements>=1])),2,paste,collapse=""),collapse="")
                
                # get name
                
                if (this.degree_oxidation>0) {
                  
                  oxystring = paste0(" +",this.degree_oxidation,"O")
                  
                } else {
                  
                  oxystring = NULL
                  
                }
     
                this.parent_compound_name = paste0(this.species," ",this.FA_total_no_C,":",this.FA_total_no_DB,oxystring)
                
                # now, record data
                
                this.frag_ID = ins.row
                
                # integers
                sim_results[ins.row,c(1,7,8,9,10)] = c(this.frag_ID,this.adduct_rank,this.FA_total_no_C,this.FA_total_no_DB,this.degree_oxidation)
                
                # doubles
                sim_results[ins.row,c(2,3)] = c(this.mz,this.parent_exactneutralmass)
                
                # text fields
                sim_results[ins.row,c(4,5,6,11,12)] = c(this.lipid_class,this.species,this.adduct,this.parent_formula,this.parent_compound_name)
                
                ins.row = ins.row + 1 # advance our insertion point
                
              }
              
              rm(m)
              
            }
            
            rm(l)
            
          }
          
          rm(k)
          
        }
        
        rm(j)
        
      }
      
    }
    
  }
  
  rm(i)
  
  # add column headings; store simulation results for this mode to a matrix of the appropriate name
  
  colnames(sim_results) = c("frag_ID","mz","exact_parent_neutral_mass","lipid_class","species","adduct","adduct_rank","FA_total_no_C","FA_total_no_DB","degree_oxidation","parent_elem_formula","parent_compound_name")
  
  # create new LOBdbase object
  
  object = new("LOBdbase")
  
  object@frag_ID = as.integer(sim_results$frag_ID)
  object@mz = as.numeric(sim_results$mz)
  object@exact_parent_neutral_mass = as.numeric(sim_results$exact_parent_neutral_mass)
  object@lipid_class = as.factor(as.character(sim_results$lipid_class))
  object@species = as.character(sim_results$species)
  object@adduct = as.factor(as.character(sim_results$adduct))
  object@adduct_rank = as.integer(sim_results$adduct_rank)
  object@FA_total_no_C = as.integer(sim_results$FA_total_no_C)
  object@FA_total_no_DB = as.integer(sim_results$FA_total_no_DB)
  object@degree_oxidation = as.integer(sim_results$degree_oxidation)
  object@parent_elem_formula = as.character(sim_results$parent_elem_formula)
  object@parent_compound_name = as.character(sim_results$parent_compound_name)
  object@polarity = as.factor(polarity)
  object@num.entries = as.integer(numCombs$numAddIons)
  object@num.compounds = as.integer(numCombs$numCompounds)
  
  # export to .csv, if user has elected the option
  
  if (gen.csv==TRUE) {
    
    exportDBtoCSV(output.dir = output.dir, LOBdbase = object)
      
  }
  
  # return the LOBdbase object, invisibly
  
  invisible(assign(paste0(as.character(object@polarity)),object))
  
}