################ Wrapper function #############

# generateLOBdbase: wrapper function for lipid-ox-lipid-oxylipin database
# generation

generateLOBdbase = function(polarity = c("positive","negative"), 
                            gen.csv = FALSE, component.defs = NULL, 
                            AIH.defs = NULL, acyl.ranges = NULL, 
                            oxy.ranges = NULL) {
  
  polarity = match.arg(polarity, several.ok = TRUE)
  
  # determine whether user has specified external source files for any input 
  # parameters, get correct file paths, and warn user of consequences of using 
  # improperly formatted or incomplete external data
  
  # perhaps later, can write a function to check the formatting of the data in 
  # each file; not a high priority right now since we presume user would have 
  # followed the layout of the default tables if creating their own
  
  if (is.null(component.defs)) { # user didn't specify external component 
    # definitions table, use defaults
    
    componentTable.loc = NULL
    use.default.componentTable = TRUE
    
  } else { # user specified something
    
    componentTable.loc = component.defs
    use.default.componentTable = FALSE
    
    # let user know he/she provided external input, and the possible 
    # consequences
    
    warning("User specified external source for basic component composition ",
            "matrix. Ensure .csv file is properly formatted and sufficient ",
            "entries exist in other external files to support any additional ",
            "adducts, lipid classes, or molecules.\n")
    
  }
  
  if (is.null(AIH.defs)) { # user didn't specify external AIH table, use 
    # defaults
    
    AIHtable.loc = NULL
    use.default.AIHtable = TRUE
    
  } else { # user specified something
    
    AIHtable.loc = AIH.defs
    use.default.AIHtable = FALSE
    
    # let user know he/she provided external input, and the possible 
    # consequences
    
    warning("User specified external source for adduct ion hierarchy ",
            "matrix. Ensure .csv file is properly formatted and sufficient ",
            "entries exist in other external files to support any additional ",
            "adducts, lipid classes, or molecules.\n")
    
  }
  
  if (is.null(acyl.ranges)) { # user didn't specify external in silico acyl 
    # property range table, use defaults
    
    acylRanges.loc = NULL
    use.default.acylRanges = TRUE
    
  } else { # user specified something
    
    acylRanges.loc = acyl.ranges
    use.default.acylRanges = FALSE
    
    # let user know he/she provided external input, and the possible 
    # consequences
    
    warning("User specified external source for in silico simulation acyl ",
            "property ranges. Ensure .csv file is properly formatted and ",
            "sufficient entries exist in other external files to support any ",
            "additional adducts, lipid classes, or molecules.\n")
    
  }
  
  if (is.null(oxy.ranges)) { # user didn't specify external AIH table, use 
    # defaults
    
    oxyRanges.loc = NULL
    use.default.oxyRanges = TRUE
    
  } else { # user specified something
    
    oxyRanges.loc = oxy.ranges
    use.default.oxyRanges = FALSE
    
    # let user know he/she provided external input, and the possible 
    # consequences
    
    warning("User specified external source for additional oxygen atoms to be ",
            "considered. Ensure .csv file is properly formatted and ",
            "sufficient entries exist in other external files to support any ",
            "additional adducts, lipid classes, or molecules.\n")
    
  }
  
  # load in silico simulation parameters using helper functions
  
  masses = calcComponentMasses(componentTable.loc,use.default.componentTable)
  ranges = loadSimRanges(acylRanges.loc,oxyRanges.loc,use.default.acylRanges,
                         use.default.oxyRanges)
  adductHierarchies = loadAIH(AIHtable.loc,use.default.AIHtable)
  
  # run simulation(s)
  
  sapply(polarity, runSim, acylRanges = ranges$acylC_DB, 
         oxyRanges = ranges$addl_oxy, adductHierarchies = adductHierarchies, 
         baseComponent.masses = masses$baseComponents, 
         adduct.masses = masses$adducts, gen.csv = gen.csv, simplify = FALSE, 
         USE.NAMES = TRUE)
  
}

################ Helper functions (some will be private) #############

# defineElemExactMasses: creates table of exact masses of the basic chemical 
# building blocks used in DB generation (elements, monatomic ions, and 2 species
# that are adduct components)

defineElemExactMasses = function() { 
  
  # specify exact masses of some necessary chemical species
  
  # values from de Laeter et al., 2003, "Atomic weights of the elements," Pure 
  # Appl. Chem. 75(6): 683-800; C. Amsler et al., 2008, "Review of particle 
  # physics," Physics Letters B667: 1
  
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
  m_Si = 27.97692649
  
  # calculate exact masses of acetonitrile and acetate using data we just 
  # specified
  
  m_ACN = 2*m_C + 3*m_H + 1*m_N
  m_Ac_minus = 2*m_C + 3*m_H + 2*m_O
  
  # create an exact-mass lookup table and put it in alphabetical order
  
  exact.masses = c(m_C,m_H,m_H_plus,m_N,m_O,m_P,m_S,m_Na,m_Cl,m_K,m_e_minus,
                   m_Mg,m_ACN,m_Ac_minus,m_Si)
  exact.masses = as.table(exact.masses)
  
  names(exact.masses) = c("m_C","m_H","m_H_plus","m_N","m_O","m_P","m_S","m_Na",
                          "m_Cl","m_K","m_e_minus","m_Mg","m_ACN","m_Ac_minus","m_Si")
  
  exact.masses = exact.masses[order(names(exact.masses))]
  
  return(exact.masses)
  
}

# calcComponentMasses: calculates exact masses of components specified in the 
# component table, using the given compositions and the exact masses produced by
# defineElemExactMasses()

# creates two tables: (1) baseComponents, with compositions and exact masses of 
# basic, non-adduct components and (2) adducts, containing exact masses of 
# adducts

calcComponentMasses = function(componentTableLoc,use.default.componentTable) {
  
  # input to calcComponentMasses should be the component composition file 
  # location obtained from getTableLocs, or other
  
  # get exact masses using defineElemExactMasses()
  
  exact.masses = defineElemExactMasses()
  
  # load in component composition matrix, or default
  
  if (use.default.componentTable==TRUE) {
    
    default.componentCompTable = NULL # to satisfy R CMD CHECK
    data(default.componentCompTable, envir = environment())
    componentCompTable = default.componentCompTable
    
  } else {
    
    componentCompTable = read.table(componentTableLoc, sep=",", header = TRUE, 
                                    row.names = 1)
    
    # # in the future, may want to consider using this alternative code instead;
    # # could help avoid problems with .csv file formatting, which can be
    # # affected by the platform (Linux, Mac, PC) on which they were created 
    # 
    # raw.componentCompTable = readChar(componentTableLoc,
    #                                   file.info(componentTableLoc)$size)
    # componentCompTable = read.table(text=gsub("\\r\\n","\\\n",
    #                                           raw.componentCompTable,
    #                                           perl=T),
    #                                 sep=",", header = TRUE, row.names = 1)

  }
  
  # first, a check to make sure user has supplied a "new" style
  # componentCompTable (i.e., after renaming and addition of new fields on or
  # about January 2017)
  
  if (sum(c("DB_gen_compound_type","Adduct_hierarchy_lookup_class") %in% 
          colnames(componentCompTable))!=2) {
    
    stop("User-supplied componentCompTable does not appear to have the ",
         "correct fields. In LOBSTAHS v1.1.3 and later, componentCompTable ",
         "must include the fields DB_gen_compound_type and ",
         "Adduct_hierarchy_lookup_class. See package documentation for ",
         "details. Aborting...\n")
    
    # stop script if this is the case
    
  }
  
  # also, a check to make sure the componentCompTable contains a column for
  # silicon atoms; must be the case for v.1.3.0 and later, but earlier 
  # versions (prior to addition of the PDMS contaminant series) did not
  # contain the Si column
  
  if (!(c("Si") %in% colnames(componentCompTable))) {
    
    stop("User-supplied componentCompTable does not appear to have the ",
         "correct fields. In LOBSTAHS v1.3.0 and later, componentCompTable ",
         "must include the field Si, for specifying the number of silicon ",
         "atoms in a given molecule. See package documentation for ",
         "details. Aborting...\n")
    # stop script if this is the case
    
  }
  
  # put columns in alphabetical order, with text fields at end
  
  componentCompTable.text = componentCompTable[,((ncol(componentCompTable)-2):
                                                   ncol(componentCompTable))]
  componentCompTable = componentCompTable[,-c((ncol(componentCompTable)-2):
                                                ncol(componentCompTable))]
  
  componentCompTable = componentCompTable[order(colnames(componentCompTable))]
  componentCompTable = cbind(componentCompTable,componentCompTable.text)

  # calculate exact masses of basic components and extract into a few separate 
  # tables
  # note: this will calculate full exact masses of all species in the component
  # composition table for which DB_gen_compound_type = "DB_unique_species" 
  
  # check to make sure we have same number of elemental building blocks in our 
  # exact.masses table and along the second dimension of the composition table
  
  if (length(exact.masses)!=(ncol(componentCompTable)-3)) {
    
    stop("Different number of chemical building blocks in the component ",
         "composition matrix and in the onboard list of exact masses. Check ",
         " your composition matrix carefully. Aborting...\n")
    # stop script if this is the case
    
  }
  
  # assuming same no. of building blocks, calculate exact masses & store as 
  # additional column in componentCompTable
  
  componentCompTable = cbind(componentCompTable,
                             matrix(data = NA, 
                                    nrow = nrow(componentCompTable), 
                                    ncol = 1)) # preallocate
  
  componentCompTable[,ncol(componentCompTable)] = 
    apply(as.matrix(sapply(componentCompTable[
      ,1:(ncol(componentCompTable)-4)], as.numeric)), 1, 
      function(x) sum(x*exact.masses,na.rm = TRUE))
  colnames(componentCompTable)[ncol(componentCompTable)] = c("Exact_mass")
  
  # extract masses of adducts into separate table (we'll need these later); 
  # create basecomponent.masses by removing adduct data from componentCompTable
  
  adduct.masses = componentCompTable[componentCompTable$Species_class %in% 
                                       c("adduct_neg","adduct_pos"),]
  
  baseComponent.masses = componentCompTable[
    !componentCompTable$Species_class  %in% c("adduct_neg","adduct_pos"),]
  
  # remove ACN and Ac- from baseComponent.masses
  
  baseComponent.masses = baseComponent.masses[
    -match(c("ACN","Ac-"),rownames(baseComponent.masses)),]
  
  # return adduct.masses, baseComponent.masses as list
  
  list(adducts=adduct.masses,baseComponents=baseComponent.masses)
  
}

# loadSimRanges: loads the necessary in silico simulation range data from the 
# acylRanges and oxyRanges file locations. creates two data frames: (1) acylC_DB
# and (2) addl_oxy

loadSimRanges = function(acylRangeTableLoc,oxyRangeTableLoc,
                         use.default.acylRanges,use.default.oxyRanges) {
  
  # inputs to loadSimRanges should be the two file locations obtained from 
  # getTableLocs, or other
  
  # load in ranges of total acyl C atoms and double bonds to be considered 
  # during simulations
  
  if (use.default.acylRanges==TRUE) {
    
    default.acylRanges = NULL # to satisfy R CMD CHECK
    data(default.acylRanges, envir = environment())
    acylRanges = default.acylRanges
    
  } else {
    
    acylRanges = read.table(acylRangeTableLoc, sep=",", skip = 1, header = TRUE)
    
  }
  
  # load in ranges of additional oxygen atoms to be considered during 
  # simulations
  
  if (use.default.oxyRanges==TRUE) {
    
    default.oxyRanges = NULL # to satisfy R CMD CHECK
    data(default.oxyRanges, envir = environment())
    oxyRanges = default.oxyRanges
    
  } else {
    
    oxyRanges = read.table(oxyRangeTableLoc, sep=",", skip = 1, header = TRUE)
    
  }
  
  list(acylC_DB=acylRanges,addl_oxy=oxyRanges)
  
}

# loadAIH: loads the adduct ion hierarchy data from the AIHfile location. 
# creates one data frame called adductHierarchies

loadAIH = function(AIHTableLoc,use.default.AIHtable) {
  
  # input to loadAIH should be file location obtained from getTableLocs, or 
  # other
  
  if (use.default.AIHtable==TRUE) {
    
    default.adductHierarchies = NULL # to satisfy R CMD CHECK
    data(default.adductHierarchies, envir = environment())
    adductHierarchies = default.adductHierarchies
    
  } else {
    
    adductHierarchies = read.table(AIHTableLoc, sep=",", skip = 1, 
                                   header = TRUE, stringsAsFactors = FALSE)
    
  }
  
  # for compatibility, also assign values in "Adduct" as row names
  
  row.names(adductHierarchies) = adductHierarchies$Adduct
  
  return(adductHierarchies)
  
}

# calcNumCombs: calculates the number of parent compounds and adduct ions for 
# which masses are to be generated in a given ionization mode, based on the 
# user-specified ranges of lipid classes and chemical properties

# this is a wrapper function; combCalc() performs the actual calculation for 
# each lipid class

calcNumCombs = function(polarity, acylRanges, oxyRanges, adductHierarchies, 
                        baseComponent.masses, adduct.masses) {
  
  # extract adduct hierarchies for this mode, compare the number to those in 
  # adduct.masses to ensure the mass of each of them is defined in the adducts 
  # mass table
  
  AIHs.thismode = adductHierarchies[adductHierarchies$Adduct_ion_mode==
                                      polarity,]
  
  if (!all(is.element(AIHs.thismode$Adduct,rownames(adduct.masses)))) {
    
    stop("Not all adduct ions given in the hierarchy table for this mode are ",
         "defined in the component composition table.\n")
    
  }
  
  # now, proceed with obtaining number of combinations
  
  combSums = apply(
    apply(
      data.frame(
        as.character(baseComponent.masses$Adduct_hierarchy_lookup_class),
        as.character(baseComponent.masses$Species_class),
        as.character(baseComponent.masses$DB_gen_compound_type)),
      1, 
      combCalc, 
      AIHs.thismode = AIHs.thismode,
      acylRanges = acylRanges,
      oxyRanges = oxyRanges),
    1,
    sum)
  
  
  return(list(numCompounds = as.integer(combSums[1]), 
              numAddIons = as.integer(combSums[2])))
  
}

# combCalc: calculates number of possible DB entires for a given lipid class;
# designed to work with the wrapper calcNumCombs, which calculates total no.
# of combinations for an entire theoretical database

combCalc = function(classInfo, AIHs.thismode, acylRanges, oxyRanges) {
  
  # classInfo: matrix consisting of the adduct hierarchy lookup classes,
  # species class name, and DB generation type
  #
  # e.g., classInfo = data.frame(
  #       as.character(baseComponent.masses$Adduct_hierarchy_lookup_class),
  #       as.character(baseComponent.masses$Species_class),
  #       as.character(baseComponent.masses$DB_gen_compound_type))
  
  # first, check to make sure baseComponent.masses$DB_gen_compound_type are of
  # acceptable kind
  
  if (!(classInfo[3] %in% c("DB_unique_species","DB_acyl_iteration",
                            "basic_component","adduct_neg","adduct_pos"))) {
    
  stop("The database generation type must be either DB_acyl_iteration, ",
       "DB_unique_species, basic_component, adduct_neg, or adduct_pos. Check ",
       "your composition matrix carefully. Aborting...\n")
  # stop script if this is the case
    
  }
  
  # retrieve necessary data for this class
  
  if (classInfo[3]=="DB_acyl_iteration") {
    
    # i.e., if this is a lipid class for which we will be generating entries for
    # different molecules with various numbers of acyl C, DB, and additional
    # oxygen atoms, e.g., "IP_DAG","IP_MAG","FFA","TAG","PUA"
    
    this.oxymin = as.numeric(oxyRanges[
      grep(paste0(as.character(classInfo[2]),"_min"),
           colnames(oxyRanges))])
    this.oxymax = as.numeric(oxyRanges[
      grep(paste0(as.character(classInfo[2]),"_max"),
           colnames(oxyRanges))])
    this.C_DBmindata = acylRanges[
      ,grep(paste0(as.character(classInfo[2]),
                   "_min"),colnames(acylRanges))]
    this.C_DBmaxdata = acylRanges[
      ,grep(paste0(as.character(classInfo[2]),
                   "_max"),colnames(acylRanges))]
    
    num.adducts = sum(!is.na(AIHs.thismode[,colnames(AIHs.thismode)[
      colnames(AIHs.thismode)==classInfo[1]]]))
    
    if (num.adducts>0) {
      
      num_compounds.this_species = (this.oxymax-this.oxymin+1)*
        sum(this.C_DBmaxdata-this.C_DBmindata+1,na.rm = TRUE)
      
    } else {
      
      num_compounds.this_species = 0
      
    }
    
    num_ions.this_species = num.adducts*num_compounds.this_species
    
  } else if (classInfo[3]=="DB_unique_species") {
    
    # this is a unique molecular species for which there will be no iteration
    # (according to user's specifications in the "DB_gen_compound_type" field of
    # the basic component matrix)
    
    num.adducts = sum(!is.na(AIHs.thismode[
      ,colnames(AIHs.thismode)[colnames(AIHs.thismode)==
                                 classInfo[1]]]))
    
    if (num.adducts>0) {
      
      num_compounds.this_species = 1 # because each of these "DB_unique_species"
      # entries represents only one compound
      
    } else {
      
      num_compounds.this_species = 0
      
    }
    
    num_ions.this_species = num.adducts
    
  }
  
  # return number of compounds and adduct ions for this species
  
  c(num_compounds.this_species, num_ions.this_species)
  
}

# genTimeStamp: generates a timestamp string based on the current system time

genTimeStamp = function () {
  
  output_DTG = format(Sys.time(), "%Y-%m-%dT%X%z") # return current time in a 
  # good format
  output_DTG = gsub(" ", "_", output_DTG) # replace any spaces
  output_DTG = gsub(":", "-", output_DTG) # replaces any colons with dashes (Mac
  # compatibility)
  
}

# exportDBtoCSV: writes a LOBdbase object to file

exportDBtoCSV = function(LOBdbase) {
  
  output_DTG = genTimeStamp()
  
  cat("Exporting .csv file containing",as.character(polarity(LOBdbase)),
      "mode simulation output...\n")
  
  fname = paste0("LOBSTAHS_lipid-oxy_DB_",strtrim(as.character(
    polarity(LOBdbase)),3),"_",output_DTG,".csv")
  
  exportmat = data.frame(frag_ID(LOBdbase),
                         mz(LOBdbase),
                         exact_parent_neutral_mass(LOBdbase),
                         as.character(lipid_class(LOBdbase)),
                         as.character(species(LOBdbase)),
                         as.character(adduct(LOBdbase)),
                         as.character(adduct_rank(LOBdbase)),
                         FA_total_no_C(LOBdbase),
                         FA_total_no_DB(LOBdbase),
                         degree_oxidation(LOBdbase),
                         parent_elem_formula(LOBdbase),
                         parent_compound_name(LOBdbase),
                         stringsAsFactors = FALSE)
  
  colnames(exportmat) = c("frag_ID","mz","exact_parent_neutral_mass",
                          "lipid_class","species","adduct","adduct_rank",
                          "FA_total_no_C","FA_total_no_DB","degree_oxidation",
                          "parent_elem_formula","parent_compound_name")
  
  write.csv(exportmat, fname)
  
  cat(as.character(polarity(LOBdbase)),"mode database exported to:",fname,"\n")
  
}

# runSim: runs the in silico simulation for a given ionization mode, returns a 
# LOBdbase object

runSim = function(polarity, acylRanges, oxyRanges, adductHierarchies, 
                  baseComponent.masses, adduct.masses, gen.csv) {
  
  # calculate number of combinations for which data are to be calculated in this
  # mode
  
  numCombs = calcNumCombs(polarity, acylRanges, oxyRanges, adductHierarchies, 
                          baseComponent.masses, adduct.masses)
  
  # provide feedback to user
  
  cat("\nPerforming in silico simulation to generate data for",polarity,
      "mode species...\n\n")
  
  cat("Based on user settings, LOBSTAHS will generate exact masses and adduct",
      "hierarchy data for",numCombs$numAddIons,"possible",polarity,"mode",
      "adduct ions, representing",numCombs$numCompounds,"unique parent",
      "compounds.\n\n")
  
  # preallocate matrix for results
  
  sim_results = matrix(data = NA, nrow = numCombs$numAddIons, ncol = 12)
  sim_results = as.data.frame(sim_results)
  
  # extract adduct hierarchies for this ionization mode
  
  AIHs.thismode = adductHierarchies[adductHierarchies$Adduct_ion_mode==
                                      polarity,]
  
  # get exact masses using defineElemExactMasses()
  
  exact.masses = defineElemExactMasses()
  
  # now, perform simulation by species
  
  ins.row = 1 # define variable to keep track of insertion point in results 
  # table
  
  for (i in 1:nrow(baseComponent.masses)) {
    
    # retrieve, store this.lipid_class, this.species, this.DB_gen_compound_type,
    # this.adduct_lkup_class
    this.lipid_class = as.character(baseComponent.masses$Species_class[i])
    this.species = rownames(baseComponent.masses)[i]
    this.DB_gen_compound_type = 
      as.character(baseComponent.masses$DB_gen_compound_type[i])
    this.adduct_lkup_class = 
      as.character(baseComponent.masses$Adduct_hierarchy_lookup_class[i])

    # provide sensible feedback to user
    
    # first, check to make sure this.DB_gen_compound_type is of an
    # acceptable type (at least in this version of LOBSTAHS)
    
    if (!(this.DB_gen_compound_type %in% c("DB_unique_species",
                                           "DB_acyl_iteration",
                                           "basic_component","adduct_neg",
                                           "adduct_pos"))) {
      
      stop("The database generation type must be either DB_acyl_iteration, ",
           "DB_unique_species, basic_component, adduct_neg, or adduct_pos. ",
           "Check your composition matrix carefully. Aborting...\n")
      # stop script if this is the case
      
    }
    
    # now, generate a logical feedback string
    
    if (this.DB_gen_compound_type=="DB_acyl_iteration") {
      
      if (this.lipid_class==this.species) {
        
        cat("Calculating data for lipid class:",this.species,"...\n")
        
      } else {
        
        cat("Calculating data for",this.lipid_class,"lipid class:",this.species,
            "...\n")
        
      }
      
    } else if (this.DB_gen_compound_type=="DB_unique_species") {
      
      cat("Calculating data for",this.lipid_class,":",this.species,"...\n")
      
    }
    
    # now, check to see whether hierarchy data exists for this species or lipid 
    # class in this ionization mode
    
    if (sum(AIHs.thismode[,colnames(AIHs.thismode)[
      colnames(AIHs.thismode)==this.adduct_lkup_class]],na.rm = TRUE)>=1) { 
      
      # hierarchy data exists for this species or class in this mode; proceed
      
      # get number of adducts formed by this species/class
      
      adducts.thisclass = rownames(AIHs.thismode)[
        !is.na(AIHs.thismode[,as.character(this.adduct_lkup_class)])]
      
      # check to see whether this is a "one-off", i.e., a unique species for
      # which we aren't considering ranges # of acyl C, double bonds, etc.
      
      if (this.DB_gen_compound_type=="DB_unique_species") {
        
        # this element is "one-off"; we only need to drill down to the
        # adduct level
        
        for (j in 1:length(adducts.thisclass)) { # cycle thru adducts
          
          # calculate mz for this adduct ion of this species
          
          this.mz = adduct.masses$Exact_mass[
            rownames(adduct.masses) %in% adducts.thisclass[j]]+
            baseComponent.masses$Exact_mass[i]
          
          # ascertain parent compound elemental formula
          
          these.base_elements = subset(baseComponent.masses[i,], 
                                       select=-c(Species_class,Exact_mass,
                                                 Adduct_hierarchy_lookup_class,
                                                 DB_gen_compound_type))          
          this.parent_formula = paste0(apply(rbind(colnames(
            these.base_elements)[these.base_elements>=1],
            as.numeric(these.base_elements[these.base_elements>=1])),
            2,paste,collapse=""),collapse="")
          
          # get name
          
          this.parent_compound_name = this.species
          
          # now, record data for this adduct ion
          
          this.frag_ID = ins.row
          this.parent_exactneutralmass = baseComponent.masses$Exact_mass[i]
          this.degree_oxidation = 0
          this.adduct = adducts.thisclass[j]
          this.adduct_rank = AIHs.thismode[
            rownames(AIHs.thismode)==this.adduct,
            colnames(AIHs.thismode)==this.adduct_lkup_class]
          this.FA_total_no_C = NA
          this.FA_total_no_DB = NA
          
          # integers
          sim_results[ins.row,c(1,7,8,9,10)] = c(this.frag_ID,this.adduct_rank,
                                                 this.FA_total_no_C,
                                                 this.FA_total_no_DB,
                                                 this.degree_oxidation)
          
          # doubles
          sim_results[ins.row,c(2,3)] = c(this.mz,this.parent_exactneutralmass)
          
          # text fields
          sim_results[ins.row,c(4,5,6,11,12)] = c(this.lipid_class,this.species,
                                                  this.adduct,
                                                  this.parent_formula,
                                                  this.parent_compound_name)
          
          ins.row = ins.row + 1 # advance our insertion point
          
        }
        
        rm(j)
        
      } else if (this.DB_gen_compound_type=="DB_acyl_iteration") {
        # this species requires more involved simulation/iteration
        
        # retrieve, store "base" exact mass for this lipid class
        
        this.base_exactmass = baseComponent.masses$Exact_mass[i]
        
        # load acyl C - double bond distributions from sim.ranges; oxidation 
        # state parameters from user-specified oxy_range variable; set variable
        # for number of carboxyl groups (for calculation IP-DAG, IP-MAG, FFA,
        # TAG, PUA masses)
        
        these.sim.ranges = acylRanges[!is.na(acylRanges[,grep(paste0(
          as.character(baseComponent.masses$Species_class[i]),"_min"),
          colnames(acylRanges))]),c("FA_total_no_C",paste0(as.character(
            baseComponent.masses$Species_class[i]),"_min"),paste0(
              as.character(baseComponent.masses$Species_class[i]),"_max"))]
        this.oxymin = as.numeric(oxyRanges[grep(paste0(
          as.character(baseComponent.masses$Species_class[i]),"_min"),
          colnames(oxyRanges))])
        this.oxymax = as.numeric(oxyRanges[grep(paste0(
          as.character(baseComponent.masses$Species_class[i]),"_max"),
          colnames(oxyRanges))])
        these.oxystates = this.oxymin:this.oxymax
        
        num.carboxyl = switch(
          as.character(baseComponent.masses$Species_class[i]),
          IP_DAG=2,
          IP_MAG=1,
          FFA=1,
          PUA=0,
          TAG=3)
        
        for (j in 1:nrow(these.sim.ranges)) { # cycle thru number of allowable 
          # acyl C atoms
          
          # retrieve, store this.FA_total_no_C
          this.FA_total_no_C = these.sim.ranges[j,1]
          
          # generate vector of allowable acyl C=C double bonds for this 
          # FA_total_no_C
          these.allowable.DBs = these.sim.ranges[j,2]:these.sim.ranges[j,3]
          
          for (k in 1:length(these.allowable.DBs)) { # cycle thru allowable no. 
            # of double bonds
            
            # retrieve, store this.FA_total_no_DB
            this.FA_total_no_DB = these.allowable.DBs[k]
            
            for (l in 1:length(these.oxystates)) { # cycle thru allowable 
              # additional oxygen atoms
              
              # retrieve, store this.degree_oxidation
              this.degree_oxidation = these.oxystates[l]
              
              for (m in 1:length(adducts.thisclass)) { # cycle thru adducts
                
                # retrieve, store this.adduct, this.adduct_exact_mass, 
                #this.adduct_rank
                this.adduct_exact_mass = adduct.masses$Exact_mass[
                  rownames(adduct.masses) %in% adducts.thisclass[m]]
                this.adduct = adducts.thisclass[m]
                this.adduct_rank = AIHs.thismode[
                  rownames(AIHs.thismode)==this.adduct,
                  colnames(AIHs.thismode)==this.adduct_lkup_class]
                
                # now, finally, can assemble and record data for this 
                # combination
                
                # calculate exact neutral mass of parent compound that reflects 
                # this combination of properties, and mz for this adduct ion of 
                # that compound 
                
                this.parent_exactneutralmass = this.base_exactmass+
                  exact.masses[c("m_C")]*this.FA_total_no_C+
                  exact.masses[c("m_H")]*(this.FA_total_no_C*2)-
                  exact.masses[c("m_H")]*num.carboxyl+
                  exact.masses[c("m_O")]*this.degree_oxidation-
                  exact.masses[c("m_H")]*(this.FA_total_no_DB*2)
                
                this.mz = this.parent_exactneutralmass+this.adduct_exact_mass
                
                # ascertain parent compound elemental formula
                
                Species_class = NULL # to satisfy R CMD CHECK
                Exact_mass = NULL # to satisfy R CMD CHECK
                Adduct_hierarchy_lookup_class = NULL # to satisfy R CMD CHECK
                DB_gen_compound_type = NULL # to satisfy R CMD CHECK
                
                these.base_elements = 
                  subset(baseComponent.masses[i,], 
                         select=-c(Species_class,Exact_mass,
                                   Adduct_hierarchy_lookup_class,
                                   DB_gen_compound_type)) 
                these.base_elements["C"] = 
                  these.base_elements["C"]+
                  this.FA_total_no_C
                these.base_elements["H"] = 
                  these.base_elements["H"]+
                  this.FA_total_no_C*2-
                  num.carboxyl-
                  this.FA_total_no_DB*2
                these.base_elements["O"] = 
                  these.base_elements["O"]+
                  this.degree_oxidation
                this.parent_formula = 
                  paste0(apply(rbind(colnames(these.base_elements)[
                    these.base_elements>=1],as.numeric(these.base_elements[
                      these.base_elements>=1])),2,paste,collapse=""),
                    collapse="")
                
                # get name
                
                if (this.degree_oxidation>0) {
                  
                  oxystring = paste0(" +",this.degree_oxidation,"O")
                  
                } else {
                  
                  oxystring = NULL
                  
                }
                
                this.parent_compound_name = paste0(this.species," ",
                                                   this.FA_total_no_C,":",
                                                   this.FA_total_no_DB,
                                                   oxystring)
                
                # now, record data
                
                this.frag_ID = ins.row
                
                # integers
                sim_results[ins.row,c(1,7,8,9,10)] = c(this.frag_ID,
                                                       this.adduct_rank,
                                                       this.FA_total_no_C,
                                                       this.FA_total_no_DB,
                                                       this.degree_oxidation)
                
                # doubles
                sim_results[ins.row,c(2,3)] = c(this.mz,
                                                this.parent_exactneutralmass)
                
                # text fields
                sim_results[ins.row,c(4,5,6,11,12)] = c(this.lipid_class,
                                                        this.species,
                                                        this.adduct,
                                                        this.parent_formula,
                                                     this.parent_compound_name)
                
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
  
  # add column headings; store simulation results for this mode to a matrix of 
  # the appropriate name
  
  colnames(sim_results) = c("frag_ID","mz","exact_parent_neutral_mass",
                            "lipid_class","species","adduct","adduct_rank",
                            "FA_total_no_C","FA_total_no_DB","degree_oxidation",
                            "parent_elem_formula","parent_compound_name")
  
  # create new LOBdbase object
  
  object = new("LOBdbase")
  
  object@frag_ID = as.integer(sim_results$frag_ID)
  mz(object) = round(as.numeric(sim_results$mz),8)
  # 8 decimal places is far more than enough precision for what we need; also, 
  # this is currently the limit of precision of our input exact mass data -- so 
  # we shouldn't exceed it. if we don't round here, we'll get some very small 
  # arithmetic errors beyond the 14th decimal place that will throw off 
  # distinction between isobars and true isomers
  exact_parent_neutral_mass(object) = 
    as.numeric(sim_results$exact_parent_neutral_mass)
  lipid_class(object) = as.factor(as.character(sim_results$lipid_class))
  species(object) = as.character(sim_results$species)
  adduct(object) = as.factor(as.character(sim_results$adduct))
  adduct_rank(object) = as.integer(sim_results$adduct_rank)
  FA_total_no_C(object) = as.integer(sim_results$FA_total_no_C)
  FA_total_no_DB(object) = as.integer(sim_results$FA_total_no_DB)
  degree_oxidation(object) = as.integer(sim_results$degree_oxidation)
  parent_elem_formula(object) = as.character(sim_results$parent_elem_formula)
  parent_compound_name(object) = as.character(sim_results$parent_compound_name)
  polarity(object) = as.factor(polarity)
  num_entries(object) = as.integer(numCombs$numAddIons)
  num_compounds(object) = as.integer(numCombs$numCompounds)
  
  # export to .csv, if user has elected the option
  
  if (gen.csv==TRUE) {
    
    exportDBtoCSV(LOBdbase = object)
    
  }
  
  # return the LOBdbase object, invisibly
  
  invisible(assign(paste0(as.character(polarity(object))),object))
  
}
