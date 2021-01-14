############################################################
#   INRA - ANR project Pseudorasbora
#
#       genind2structure
#
############################################################
# by Thomas Brazier
# brazier.thomas@gmail.com
# MSc internship

# Supervisor: Scott McCairns
#             INRA UMR ESE

# Lindsay V. Clark, 26 July 2015, initiated a genind2structure function, basic yet efficient, with some features lacking
# https://github.com/lvclark/R_genetics_conv/blob/master/genind2structure.R
# Here is some improvements of this first function...

#==========================================================
# FUNCTION genind2str
#==========================================================


#---------------------------------------------------------
# DESCRIPTION

# Function to export to STRUCTURE format from genind object.
# genind objects are created in the R package adegenet. Function below is an R function.

# STRUCTURE input files are described as following by the user manual:
# Each row of individual data contains the following elements. These form columns in the data file.
# 1. Label (Optional; string) A string of integers or characters used to designate each individual in the sample.
# 2. PopData (Optional; integer) An integer designating a user-defined population
# 3. PopFlag (Optional; 0 or 1) A Boolean flag which indicates whether to use the PopData
#       when using learning samples (Note: A Boolean variable (flag) is a variable which takes the values TRUE or FALSE, which are designated here by the integers
#       1 (use PopData) and 0 (don’t use PopData), respectively.)
# 4. LocData (Optional; integer) An integer designating a user-defined sampling location
# 5. Phenotype (Optional; integer) An integer designating the value of a phenotype of interest
# 6. Extra Columns (Optional; string) It may be convenient for the user to include additional data in the input file which are ignored by the program.
# 7. Genotype Data (Required; integer) Each allele at a given locus should be coded by a unique integer (eg microsatellite repeat score).

#---------------------------------------------------------
# USAGE
# genind2structure(obj, file="", pops=FALSE, missingdata.label=-9, popflag=NA, locdata=NA, phenotype=NA, extra=NA)

#---------------------------------------------------------
# ARGUMENTS
# obj: a genind object to convert
# file: file name to write
# pops: whether to include population info in the file, population names are taken in the genind object
# missingdata.label: value is the expected value for missing data
# popflag: if popflag is not NA, then it is a vector of popflags (numeric, 0 or 1) of the same size as the number of individuals
# locdata: if locdata is not NA, then it is a a numeric vector of user-defined sampling locations
# phenotype: if phenotype is not NA, then it is a numeric vector of phenotypes
# extra: if extra is not NA, then it is a character vector with extra informations

# example use: 
# data(nancycats)
# genind2structure(nancycats, file="nancy_structure.txt", pops=TRUE)

#---------------------------------------------------------
# VALUE

# Return either a data.frame if file="", or write the data.frame into the file given in argument.


#---------------------------------------------------------
genind2structure = function(obj, file="", pops=TRUE, missingdata.label=-9, popflag=NA, locdata=NA,
                            phenotype=NA, extra=NA){
  # Dependencies
  if ('adegenet' %in% installed.packages()) {
    library('adegenet')
  } else {
    warning("adegenet is not installed.")
    ans = readline("Do you want to install it? (y/n)")
    if (ans=="y" | ans=="yes") {
      install.packages('adegenet')
    } else {
      stop("Function ended because adegenet is required and not installed.")
    }
  }
  
  # Check first that the object is a genind object; if not, stop the function
  if(!"genind" %in% class(obj)){
    stop("Function was designed for genind objects.")
  }
  
  # get the max ploidy of the dataset
  ploidy = max(obj@ploidy)
  # get the number of individuals
  Nind = nInd(obj)
  
  #-------------------------
  # 1. ADD IND. LABEL
  # column of individual names to write; set up data.frame
  tab = data.frame(ind=rep(indNames(obj), each=ploidy))
  
  # 2. ADD POP. LABEL
  # column of pop ids to write
  if(pops){
    popnums = 1:nPop(obj)
    names(popnums) = as.character(unique(pop(obj)))
    popcol = rep(popnums[as.character(pop(obj))], each=ploidy)
    tab = cbind(tab, data.frame(pop=popcol))
  }
  
  # 3. Add POPFLAG column if needed
  # Consider the 2 cases where popflag is the same length as structure data frame (replicated values for ploidy)
  # or same length as number of individuals (just one value per individual that we need to replicate given the ploidy level)
  if (all(!is.na(popflag))) {
    if (is.numeric(popflag)) {
      if (length(popflag)==Nind) {
        tab=cbind(tab, data.frame(popflag=rep(popflag, each=ploidy))) # popflags need to be replicated given the ploidy level
      } else {
        if (length(popflag)==nrow(tab)) { # popflags do not need to be replicated given the ploidy
          tab=cbind(tab, data.frame(popflag=popflag))
        }else {
          warning("Not the correct number of popflags in the list")
        }
      }
    } else {
      warning("Popflag is not a numeric vector")
    }
  }
  
  # 4. ADD LOCDATA
  if (all(!is.na(locdata))) {
    if (is.numeric(locdata)) {
      if (length(locdata)==Nind) {
        tab=cbind(tab, data.frame(locdata=rep(locdata, each=ploidy))) # locdatas need to be replicated given the ploidy level
      } else {
        if (length(locdata)==nrow(tab)) { # locdatas do not need to be replicated given the ploidy
          tab=cbind(tab, data.frame(locdata=locdata))
        }else {
          warning("Not the correct number of locdatas in the list")
        }
      }
    } else {
      warning("Locdata is not a numeric vector")
    }
  }
  
  # 5. ADD PHENOTYPE
  if (all(!is.na(phenotype))) {
    if (is.numeric(phenotype)) {
      if (length(phenotype)==Nind) {
        tab=cbind(tab, data.frame(phenotype=rep(phenotype, each=ploidy))) # phenotypes need to be replicated given the ploidy level
      } else {
        if (length(phenotype)==nrow(tab)) { # phenotypes do not need to be replicated given the ploidy
          tab=cbind(tab, data.frame(phenotype=phenotype))
        }else {
          warning("Not the correct number of phenotypes in the list")
        }
      }
    } else {
      warning("Phenotype is not a numeric vector")
    }
  }
  
  
  # 6. ADD EXTRA COLUMN
  if (all(!is.na(extra))) {
    if (is.character(extra)) {
      if (length(extra)==nrow(tab)) {
        tab=cbind(tab, data.frame(extra=rep(extra, each=ploidy)))
      } else {
        warning("Not as many extra informations as individuals")
      }
    } else {
      warning("Extra information is not a character vector")
    }
  }
  
  if (all(!is.na(extra))) {
    if (is.character(extra)) {
      if (length(extra)==Nind) {
        tab=cbind(tab, data.frame(extra=rep(extra, each=ploidy))) # extras need to be replicated given the ploidy level
      } else {
        if (length(extra)==nrow(tab)) { # extras do not need to be replicated given the ploidy
          tab=cbind(tab, data.frame(extra=extra))
        }else {
          warning("Not the correct number of extras in the list")
        }
      }
    } else {
      warning("Extra is not a character vector")
    }
  }
  
  # 7. ADD GENOTYPES
  # get loci names
  loci = locNames(obj) 
  # add columns for genotypes
  tab = cbind(tab, matrix(missingdata.label, nrow=dim(tab)[1], ncol=nLoc(obj), dimnames=list(NULL,loci)))
  # begin going through loci
  for(L in loci){
    gen = obj@tab[,grep(paste("^", L, "\\.", sep=""), 
                              dimnames(obj@tab)[[2]]), 
                        drop = FALSE] # genotypes by locus
    al = 1:dim(gen)[2] # numbered alleles
    for(s in 1:Nind){
      if(all(!is.na(gen[s,]))){
        tabrows = (1:dim(tab)[1])[tab[[1]] == indNames(obj)[s]] # index of rows in output to write to
        tabrows = tabrows[1:sum(gen[s,])] # subset if this is lower ploidy than max ploidy
        tab[tabrows,L] = rep(al, times = gen[s,])
      }
    }
  }

  #-------------------------------
  # export table as a file (if no filename specified, then tab is returned at the end of the function)
  if (file != "") {
    write.table(tab, file=file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  } else {
    return(tab)
  }
}

#==========================================================
# TESTS
#==========================================================
# data(nancycats)

# Test with a genind object and no option, no output file
# genind2structure(nancycats, file="", pops=TRUE)
# Test with a genind object and no option, with output file
# genind2structure(nancycats, file="nancy_structure.txt", pops=TRUE)
# Test with a genind object and a correct POPFLAG option, with output file
# popflag=rep(1,2*nrow(nancycats@tab)) # create a correct vector of popflag (length = nrow of structure file given ploidy)
# or
# popflag=rep(1,nrow(nancycats@tab)) # create a correct vector of popflag (length = number of individuals)
# genind2structure(nancycats, file="nancy_structure.txt", pops=TRUE, popflag=popflag)


#==========================================================
# THE END
#==========================================================