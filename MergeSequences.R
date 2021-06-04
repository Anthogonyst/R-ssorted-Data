
##########
### Summary: Pipeline for merging fasta files and exporting to our collaborators
###
### Input: Filepaths of fasta folders
### Output: Sanitized copies of merged and unmerged fasta files, with metadata
###
###
### Author: Anthony A
### Date created: 3/22/2021
### Date updated: 5/24/2021
##########


### Packages
source("~/GISAIDWorkplace/PathUtils.R")
packages = c("seqinr", "tidyr", "stringr", "testthat", "tictoc")
LoadPackages(packages)
rm(packages)


### Pipeline for exporting the GISAID fasta files
PipelineGISAID <- function(path, excludes = NA, PathFUN = .PathsDev, pipeOutFastas = FALSE) {
  ### Initializing
  tic("Gisaid Pipeline finished!")
  prefix = "gisaid"
  root = PathFUN()[1]
  .DateAndTime(TRUE)
  
  ### Write destinations
  outputFolder = .PathsOut(prefix)
  mergeOutputFolder = .PathsMergeOut()
  
  ### Saves a table for the sample's year of collection
  yearsTable = FetchYears()
  
  ### The pipeline returns a list of fasta files
  theFastaFiles = HandleSubpaths(path) %>%
    ParseNestedFolder(root, .) %>%
      HandleBadFiles(.) %>%
        LoadFiles(.) %>%
          ExcludeFiles(., excludes) %>%
            ModifyNames(., yearsTable)  %>%
              AddFastaAttributes(.)
  
  ### Writes the files onto the drive
  WriteOrganizationTaggedFasta(theFastaFiles, outputFolder)
  WriteMergedFasta(theFastaFiles, mergeOutputFolder, prefix)
  toc("Gisaid Pipeline finished!")
  
  ### Returns the filepaths of the resulting folders unless fastas are desired
  if (pipeOutFastas)
    return(theFastaFiles)
  else return(c(outputFolder, mergeOutputFolder))
}


### Pipeline for exporting the University fasta files
PipelineUniversity <- function(path, excludes = NA, PathFUN = .PathsArchive, pipeOutFastas = FALSE) {
  ### Initializing
  tic("University Pipeline finished!")
  analysisFolder = PathFUN()[1]
  prefix = "University"
  .DateAndTime(TRUE)
  
  ### Write destinations
  outputFolder = .PathsOut(prefix)
  mergeOutputFolder = .PathsMergeOut()
  
  ### Saves a table for the sample's year of collection
  yearsTable = FetchYears()
  
  ### The pipeline returns a list of fasta files
  theFastaFiles = ParseFolder(analysisFolder) %>%
    HandleBadFiles(.) %>%
      LoadFiles(.) %>%
        ExcludeFiles(., excludes) %>%
          ModifyNames(., yearsTable)  %>%
            AddFastaAttributes(.) %>%
              UniversityAdditionalErrorChecking(., yearsTable)
  
  ### Writes the files onto the drive
  WriteOrganizationTaggedFasta(theFastaFiles, outputFolder)
  WriteMergedFasta(theFastaFiles, mergeOutputFolder, prefix)
  toc("Gisaid Pipeline finished!")
  
  
  ### Returns the filepaths of the resulting folders unless fastas are desired
  if (pipeOutFastas)
    return(theFastaFiles)
  else return(c(outputFolder, mergeOutputFolder))
}


### Pipeline for merging files with the expectation of being similar to GISAID
PipelineSemiGISAID <- function(path, filesArray = NA, excludes = NA,
                               MetadataFUN = FetchYears, prefix = "explicit",
                               ExtraFUN = .Empty, pipeOutFastas = FALSE) {
  ### Initializing
  tic("Alternate Gisaid Pipeline finished!")
  .DateAndTime(TRUE)
  
  ### Write destinations
  outputFolder = .PathsOut(prefix)
  mergeOutputFolder = .PathsMergeOut()
  
  ### Saves a table for the sample's year of collection
  yearsTable = MetadataFUN()
  
  ### The pipeline returns a list of fasta files
  theFastaFiles = ParseFlexibleFileFolder(path) %>%
    LoadFiles(.) %>%
      FilterFor(., filesArray) %>%
        ExcludeFiles(., excludes) %>%
          ModifyNames(., yearsTable) %>%
            ExtraFUN(.) %>%
              AddFastaAttributes(.)
  
  ### Writes the files onto the drive
  WriteOrganizationTaggedFasta(theFastaFiles, outputFolder)
  WriteMergedFasta(theFastaFiles, mergeOutputFolder, prefix)
  toc("Alternate Gisaid Pipeline finished!")
  
  ### Returns the filepaths of the resulting folders unless fastas are desired
  if (pipeOutFastas)
    return(theFastaFiles)
  else return(c(outputFolder, mergeOutputFolder))
}


### Merges only the specified files
PipelineGenericCalls <- function(path, filesArray = NA, prefix = "generic", pipeOutFastas = FALSE) {
  tic("Pipeline finished!", quiet = FALSE)
  .DateAndTime(TRUE)
  
  ### Parse anything to files
  files = ParseFlexibleFileFolder(path)
  
  ### Write destinations
  outputFolder = .PathsOut(prefix)
  mergeOutputFolder = .PathsMergeOut()
  
  ### Saves a table for the sample's year of collection
  yearsTable = FetchYears()
  
  ### The pipeline returns a list of fasta files
  theFastaFiles = LoadFiles(files) %>%
    FilterFor(., filesArray) %>%
      AddFastaAttributes(.)

  ### Writes the files onto the drive
  WriteSanitizedFasta(theFastaFiles, outputFolder)
  WriteMergedFasta(theFastaFiles, mergeOutputFolder, prefix)
  toc("Pipeline finished!")
  
  
  ### Returns the filepaths of the resulting folders unless fastas are desired
  if (pipeOutFastas)
    return(theFastaFiles)
  else return(c(outputFolder, mergeOutputFolder))
}


### Easy way to split an already clean, merged fasta file
SplitMergedFile <- function(path, outputPath, appendNumber = FALSE) {
  ### If it detects a folder instead of file, reads all files in folder
  if (!file.exists(path)) {
    ParseFolder(path) %>%
      LoadFiles(.) %>%
        WriteUntaggedFasta(., outputPath, appendNumber)
  } else {
    LoadFiles(path) %>%
      WriteUntaggedFasta(., outputPath, appendNumber)
  }
}

### Easier for function calls
ParseFolder <- function(folders, childFolder = "", root = "", pattern = ".*\\.fa(sta)?$", hardStop = TRUE) {
  ### Gets a list of filepaths for each folder conforming to the regex
  listFilePaths = lapply(folders, function(x) {
    list.files(path = paste0(root, x, childFolder),
               pattern = pattern, all.files = TRUE, full.names = TRUE)
  }) %>%
    unlist(., recursive = FALSE)
  
  ### Failsafe with optional stop
  if (length(listFilePaths) > 0)
    print("Found files!")
  else if (hardStop) {
    stop("No files found. Please check your path.")
  } else print("No files found. Please check your path.")
  
  return(listFilePaths)
}


### More human friendly version of ParseFolder
ParseNestedFolder <- function(root = "", folders, childFolder = "", pattern = ".*\\.fa(sta)?$", hardStop = TRUE) {
  ParseFolder(folders, childFolder, root, pattern, hardStop)
}


### Determines and parses anything (folders, files) to its file contents.
###   If a folder is provided, it returns every file, whereas it will only
###   take one file if a full filepath is provided for that particular one.
ParseFlexibleFileFolder <- function(direction) {
  ### Parse anything to file contents
  files = lapply(direction, function(x) {
    if (dir.exists(x)) {
      list.files(x, full.names = TRUE)
    } else {
      x
    }
  }) %>%
    unlist(., recursive = FALSE)
  
  return(files)
}


### Appends the midpath to CloudRuns
HandleSubpaths <- function(listFolderPaths, overrideToIgnore = FALSE) {
  # TODO: Merge this correctly with path functions
  # Perhaps: Get a separate path function depending on filename
  
  ### Just in case
  if (overrideToIgnore) {
    return(listFolderPaths)
  }
  
  ### The expected subpaths
  modifiedPaths = lapply(listFolderPaths, function(x) { 
    if (grepl("*CloudRun[^/]*(/)?$", x)) {
      paste0(x, "/filtered/samtools/out/")
    }
    else if (grepl("*NanoporeRun[^/]*(/)?$", x)) {
      paste0(x, "/nanoFilter/samtools/out/")
    }
    else x
  })
  
  return(modifiedPaths)
}


### Filters out files that match a regex
.HandleBadSimple <- function(listFilePaths, regex) {
  filteredFiles = lapply(listFilePaths, function(x) {
    if (grepl(regex, x))
      FALSE
    else TRUE
  }) %>%
    unlist(., recursive = FALSE)
  
  return(listFilePaths[filteredFiles])
}


### Filter CloudRuns for the consensus fasta in *filepaths*
.HandleCloudFiles <- function(listFilePaths) {
  ### For Cloud runs, filters for consensus.fasta files exclusively
  filteredFiles = lapply(listFilePaths, function(x) {
    if (grepl("*CloudRun*", x)) {
      if (grepl("*\\.consensus\\.*", x))
        TRUE
      else FALSE
    }
    else TRUE
  }) %>%
    unlist(., recursive = FALSE)
  
  return(listFilePaths[filteredFiles])
}


### Filter based on criteria that may exist in a folder
HandleBadFiles <- function(listFilePaths) {
  ### Filters out isolate and tmp files
  relevantFilePaths = .HandleCloudFiles(listFilePaths) %>%
    .HandleBadSimple(., "*isolate*") %>%
      .HandleBadSimple(., "*\\.tmp\\.fa(sta)?") %>%
        .HandleBadSimple(., "*isolate\\.fa(sta)?") %>%
          .HandleBadSimple(., "*aligned.fa(sta)?") %>%
            .HandleBadSimple(., "*merged.fa(sta)?") %>%
              .HandleBadSimple(., "*aligned_all.fa(sta)?") %>%
                .HandleBadSimple(., "*merged_all\\.fa(sta)?")
  
  ### Logs the excluded files
  .LogEvent(listFilePaths %>% .[. %notin% relevantFilePaths], "bad_filepaths_excluded")
  return(relevantFilePaths)
}


### Loads all of the fasta files specified
LoadFiles <- function(listFilePaths) {
  ### Loads all of the fasta files into memory
  fastaFiles = lapply(listFilePaths, function(x) {
    seqinr::read.fasta(file = x, as.string = TRUE, forceDNAtolower = FALSE, set.attributes = TRUE)
  }) %>% 
    unlist(., recursive = FALSE)
  
  ### Console print and recording data stream
  print("Done reading fasta files!")
  .LogEvent(listFilePaths, "loaded_filepaths")
  .LogEvent(names(fastaFiles), "loaded_files")
  
  return(fastaFiles)
}


### Returns data table for samples and their respective year of collection
#   The format is an M×2 matrix:  virus   date
FetchYears <- function() {
  ### Preparing the key-values map
  yearsForSamples = LoadMetadata("University_export") %>% 
    dplyr::transmute(virus = GISAID.name, date = Collection.Date)

  ### Process data to separate sample name and collection year
  yearsForSamples$virus = yearsForSamples$virus %>%
    substring(0, nchar(.) - 4)
  yearsForSamples$date = yearsForSamples$date %>%
    substring(0, 4)
  
  return(yearsForSamples)
}


### Reads a generic table in case you're using samples not on GISAID_Export
YearsFromText <- function(path = "~/GISAIDWorkplace/allDataPlusYears-20210505.tsv",
                          hasColnames = TRUE, delim = "\t") {
  yearsForSamples = read.table(path, header = hasColnames, sep = delim, stringsAsFactors = FALSE)
  colnames(yearsForSamples) = c("virus", "date")
  
  return(yearsForSamples)
}

### Helper function to remove files
.ExcludeFromSet <- function(x, y) {
  x[!(names(x) %in% y)]
}


### Removes files specified from the set and returns the new set
ExcludeFiles <- function(fastaFiles, exclude, isFile = "automatic") {
  ### Checks if there is a list of excluded files
  if (!is.na(exclude[1])) {
    ### Tracks before/after running
    beforeExcludes = names(fastaFiles)
    
    ### Assumes it is a file if there is a path attached with a / in it
    if (typeof(isFile) == "logical") {
      parseAsFile = isFile[1]
    } else {
      if (grepl("/", exclude))
        parseAsFile = TRUE
      else parseAsFile = FALSE
    }
    
    ### Works to convert the file to a vector if not already, and then filters
    if (parseAsFile) {
      excludedFilesList = .ReadTextFile(exclude)
      fastaFiles = .ExcludeFromSet(fastaFiles, excludedFilesList)
    } else {
      excludedFilesList = exclude
      fastaFiles = .ExcludeFromSet(fastaFiles, excludedFilesList)
    }
    
    ### Logs the excluded files
    .LogEvent(beforeExcludes %>% .[. %notin% names(fastaFiles)], "post_exclusions")
    .LogEvent(names(fastaFiles), "post_exclude_files_loaded")
  }
  
  return(fastaFiles)
}


### Parses a set of data for specific files
FilterFor <- function(fastaFiles, desiredFiles) {
  if (is.na(desiredFiles)) {
    print("Received NA when filtering files; assuming no filter.")
    return(fastaFiles)
  }
  
  found = fastaFiles %>% .[names(.) %in% desiredFiles]
  notfound = fastaFiles %>% .[names(.) %notin% desiredFiles]
  
  .LogEvent(names(found), "files_kept_after_filter")
  .LogEvent(names(notfound), "files_filtered_out")
  
  return(found)
}


### Renames legacy files that are tagged with consensus, threshold, and quality
RemoveConsensus <- function(fastaFiles) {
  unLegacy = stringr::str_replace_all(names(fastaFiles), pattern = "Consensus_19COV", replacement = "19COV") %>%
    stringr::str_replace_all(., pattern = "_threshold_\\d.\\d_quality_\\d\\d", replacement = "")
  names(fastaFiles) = unLegacy
  
  return(fastaFiles)
}

### Changes the names of all the fasta files
ModifyNames <- function(fastaFiles, yearsForSamples) {
  ### Just in case it catches Consensus
  unlegacyNames = RemoveConsensus(fastaFiles)
  
  ### Replace header and file name of fasta files with GISAID format
  revisedNames = stringr::str_replace_all(names(unlegacyNames), pattern = ".*19COV", replacement = "hCoV-19/USA/NY-NYCOrganization-") %>%
    stringr::str_replace_all(pattern = "-NYC-.*", replacement = "/") %>%
      stringr::str_replace_all(pattern = "_barcode.*", replacement = "")
  
  logs = {}
  ### Applying year to filenames
  finishedNames = lapply(revisedNames, function(x) {
    year = yearsForSamples[grep(x, yearsForSamples$virus), "date"]
    
    if (length(year) == 0) {
      ### Parses errors individually
      err = paste("Virus not logged, using current year:", x)
      warning(err)
      logs <<- c(logs, err)
      
      return(paste0(x, .DateAndTime() %>% substring(0, 4)))
    }
    
    return(paste0(x, year))
  })
  
  ### Writes the error file if there exists any
  if (length(logs) > 0) {
    .LogEvent(logs, "years_missing")
  }
  
  names(fastaFiles) = finishedNames
  return(fastaFiles)
}


### Adds attributes to each of the fasta files
AddFastaAttributes <- function(fastaFiles) {
  ### Tagging all of the data based on attributes tags so that we can write filenames easier
  names = names(fastaFiles)
  
  fastaFiles = mapply(fastaFiles, names, SIMPLIFY = FALSE,
                      # TODO: In the event that the collection date is empty in University_export, the resulting
                      #       file will be called NY-NYCOrganization-####.fa with two missing digits and the fasta
                      #       will have the year "NA" in the virus name. This is not the desired result.
                      # NOTE: This is actually two separate bugs; the first bug is caused by noMod = TRUE
                      FUN = function(x, y) {
                        if (grepl("hCoV-19/USA/NY-NYCOrganization-*", y)) {
                          if (is.null(attributes(x)$sample)) { attr(x, "sample") <- substring(y, 23, nchar(y) - 5) }
                        } else  {
                          if (is.null(attributes(x)$sample)) { attr(x, "sample") <- .SafeFilename(y, "-") }
                        }
                        if (is.null(attributes(x)$filename)) { attr(x, "filename") <- y }
                        return(x)
                      })
  
  return(fastaFiles)
}


### Double checks that all of the University files were successfully done
UniversityAdditionalErrorChecking <- function(fastaFiles, yearsForSamples) {
  # TODO: See if it is better to merge this with another function FunA() { ... FuncB() }
  
  ### Excluded files
  outta = .ExcludeFromSet(fastaFiles, yearsForSamples$virus)
  outta = fastaFiles[substring(names(fastaFiles), 0, nchar(names(fastaFiles)) - 4) %in% yearsForSamples$virus]
  whatsmissing = (yearsForSamples$virus)[!(paste0(yearsForSamples$virus, yearsForSamples$date) %in% names(outta))]
  
  if (length(whatsmissing) > 0) {
    lapply(whatsmissing, function(x) { print(paste("These files are missing:", x))})
    .LogEvent(whatsmissing, "missing_years_in_University_spreadsheet")
  }
  
  return(fastaFiles)
}


### Fasta write function for files that passed through AddFastaAttributes
### to clean GISAID names
### Thus, this function is meant for fasta produced by NYC Organization
WriteOrganizationTaggedFasta <- function(fastaFiles, outputPath, forceWritePath = TRUE) {
  .CheckOutputPath(outputPath, forcePath = forceWritePath)
  
  len = lapply(fastaFiles, function(x) {
    ### For our official uploads, we have a strict naming convention that
    ### requires the files to be piped first
    if (is.null(attributes(x)$filename) | is.null(attributes(x)$sample)) {
      # Note: this acts like "continue" where the lapply keeps iterating but skips one file
      warning(paste("Please tag fasta as per Organization standard with AddFastaAttributes() first:", x))
      return(1)
    }
    
    ### If the file is polished correctly, we write it
    seqinr::write.fasta(sequences = x,
                        names = attributes(x)$filename[1], 
                        file.out = paste0(outputPath, "/NY-NYCOrganization-", attributes(x)$sample[1], ".fa"),
                        open = "w")
  })
  
  ### The return code for a successful I/O operation is 0 on most computers
  paste("Done writing", length(len %>% .[. == 0]), "files!") %>% 
    print(.)
  .LogEvent(attributes(fastaFiles)$filename, "finished_written_Organization_filenames")
  .LogEvent(attributes(fastaFiles)$sample, "finished_written_Organization_sampleIds")
  
  return(outputPath)
}


### Fasta write function for files that passed through AddFastaAttributes
### to clean any fasta passed through it
WriteSanitizedFasta <- function(fastaFiles, outputPath, forceWritePath = TRUE) {
  .CheckOutputPath(outputPath, forcePath = forceWritePath)
  
  ### Intended for external files, it will simply name the file based on the fasta name
  len = lapply(fastaFiles, function(x) {
    seqinr::write.fasta(sequences = x,
                        names = if (!is.null(attributes(x)$filename))
                                  attributes(x)$filename[1]
                                else attributes(x)$name, 
                        file.out = if (!is.null(attributes(x)$sample)) 
                                      paste0(outputPath, "/", attributes(x)$sample[1], ".fa")
                                   else attributes(x)$name %>% .SafeFilename(.) %>%
                                      paste0(outputPath, "/", .),
                        open = "w")
  })
  
  ### The return code for a successful I/O operation is 0 on most computers
  paste("Done writing", length(len %>% .[. == 0]), "files!") %>%
    print(.)
  .LogEvent(attributes(fastaFiles)$filename, "finished_written_tagged_filenames")
  .LogEvent(attributes(fastaFiles)$sample, "finished_written_tagged_sampleIds")
  
  return(outputPath)
}


### Fasta write function for files that use the generic names()
WriteUntaggedFasta <- function(fastaFiles, outputPath, appendVal = FALSE, forceWritePath = TRUE) {
  .CheckOutputPath(outputPath, forcePath = forceWritePath)
  
  ### Just in case of filename issues, appends a value since untagged files are "unclean"
  if (appendVal) {
    i = 0
    len = lapply(fastaFiles, function(x) {
      i <<- i + 1
      seqinr::write.fasta(sequences = x,
                          names = attributes(x)$name, 
                          file.out = paste0(outputPath, "/", attributes(x)$name %>% 
                                              gsub("[\\/\\:\\*\\?\"\\<\\>\\|]", "", .), i, ".fa"),
                          open = "w")
    })
  } else {
    len = lapply(fastaFiles, function(x) { 
      seqinr::write.fasta(sequences = x,
                          names = attributes(x)$name, 
                          file.out = paste0(outputPath, "/", attributes(x)$name %>% 
                                              gsub("[\\/\\:\\*\\?\"\\<\\>\\|]", "", .), ".fa"),
                          open = "w")
    })
  }
  
  ### The return code for a successful I/O operation is 0 on most computers
  paste("Done writing", length(len %>% .[. == 0]), "files!") %>%
    print(.)
  .LogEvent(attributes(fastaFiles)$name, "finished_written_untagged")
  
  return(outputPath)
}
  

### Fasta write function that creates one large merged file
WriteMergedFasta <- function(fastaFiles, outputPath, prefix, forceWritePath = TRUE) {
  ### Checks that the output path is good
  .CheckOutputPath(outputPath, forcePath = forceWritePath)
  mergeFilename = paste0(outputPath, "/", prefix, "_upload_", .DateAndTime(), ".fa")
  
  ### And writes one single merged file consisting of many sequences
  seqinr::write.fasta(sequences = fastaFiles, 
                      names = names(fastaFiles),
                      file.out = mergeFilename, 
                      open = "w")
  
  ### The return code for a successful I/O operation is 0 on most computers
  paste("Done merging", length(fastaFiles), "files!") %>%
    print(.)
  .LogEvent(names(fastaFiles), "finished_merge")
  
  return(outputPath)
}


### Another empty function meant to be substituted if necessary
.Empty <- function(x) {
  return(x)
}

### Just an empty function so that pipes don't print too much
.NoPrint <- function(x, ...) {
  print("All Done!")
  
  # Just food for thought
  DoSomethingInteresting <- function(y) { return(y) }
  
  inputList <- list(...)
  outputList <- lapply(inputList, DoSomethingInteresting)
  return(outputList)
}

