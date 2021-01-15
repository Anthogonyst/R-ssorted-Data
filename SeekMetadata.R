### Packages
packages = c("seqinr", "readxl", "dplyr", "plyr", "WriteXLS", "stringr", "openxlsx")

### Install packages not yet installed
installedPackages = packages %in% rownames(installed.packages())
if (any(installedPackages == FALSE)) {
  install.packages(packages[!installedPackages])
}

### Load packages
invisible(lapply(packages, library, character.only = TRUE))


### Default paths if script is not called from command line, but from here instead
Paths = function() {
  # Change these as necessary
  defaultMetadataPath = "Workspace/metadata/"
  defaultFastaPath =    "Workspace/sequences/"
  defaultOutputPath =   "Workspace/upload/"

  # Don't change this line
  c(defaultMetadataPath, defaultFastaPath, defaultOutputPath)
}

### Summary: Given a metadata file and a folder containing FASTA files,
###   this seeks and filters out the metadata for every FASTA record in
###   each FASTA file to reduce file size and better aggregate the data.
###
###
### Input: .xls File, .fa Folder, output Folder
### Output: .xls File
###
###
### Author: Anthony A
SeekMetadata <- function(metadataPath,
                       fastaPath,
                       outputPath,
                       verbose = FALSE,
                       forcePath = FALSE) {

### Use this if you need to debug or rip the source for some reason
#if (TRUE) { metadataPath = defaultMetadataPath; fastaPath = defaultFastaPath; outputPath = defaultOutputPath; verbose = FALSE; forcePath = FALSE; }
#fastaPath = listFiles[1]

  ### Define working directory and constants
  date = gsub("-", "", Sys.Date())
  time = gsub(":", "", substring(Sys.time(), 12))
  outFile = gsub(".fa", ".xls", fastaPath)
  failure = FALSE

  if (outFile == fastaPath) {
    outFile = paste0("gisaid_upload_", date, "_", time)
  }

  ### Initial error checking
  if (is.na(metadataPath) | is.na(fastaPath) |
      is.na(outputPath) | !file.exists(metadataPath) |
      !file.exists(fastaPath)) {
      stop("File does not exist.")
  }

  ### Make new directory if necessary
  if (!dir.exists(outputPath)) {
    if (forcePath) {
      dir.create(outputPath)
    } else {
      question = readline("Output path does not exist, create new folder? (Y/N)")
      if (regexpr(question, 'y', ignore.case = TRUE)) {
        dir.create(outputPath)
      } else {
        stop(paste0("Output path does not exist at \"", ))
      }
    }
  }

  ### Load files in, with multiple fail checks on the excel file
  fastaInput = read.fasta(file = fastaPath)

  ### Regex: "^metadata_.*\\.xlsx$"
  # ^                     Start of filename
  # metadata_             Direct match
  # .*                    Any data allowed, for the date/time
  # \\.                   Period
  # xlsx                  Direct match
  # $                     End of filename
  #######

  ## Original six+ sheet xlsx
  excelFile = tail(list.files(path = metadataPath, pattern = "^metadata_.*\\.xlsx$", all.files = TRUE, full.names = TRUE, recursive = FALSE), 1)
  if (!is.na(excelFile)) excelInput = read.xlsx(excelFile, sheet = "GISAID_export")

  ## One sheet xlsx
  if (is.na(excelFile)) {
    excelFile = tail(list.files(path = metadataPath, pattern = "^gisaid_.*\\.xlsx$", all.files = TRUE, full.names = TRUE, recursive = FALSE), 1)
    if (!is.na(excelFile)) excelInput = read.xlsx(excelFile, sheet = 1)
  }

  ## One sheet xls
  if (is.na(excelFile)) {
    excelFile = tail(list.files(path = metadataPath, pattern = "^gisaid.*\\.xls$", all.files = TRUE, full.names = TRUE, recursive = FALSE), 1)
    if (!is.na(excelFile)) excelInput = read_xls(excelFile)
  }

  ## Fail condition
  if (is.na(excelFile)) {
    failure = TRUE
    stop("Could not find excel file...")
  }

  ### Print header if necessary
  if (verbose) {
    head(excelInput)
  }

  ### Define header for spreadsheet and some error checking
  excelHeader = excelInput[1,]
  if (length(excelHeader) != 29) {
    stop("Metadata does not match expected # of columns (29).")
  }
  excelInput = excelInput[-1,]

  ### Acquire all metadata related to FASTA files
  woahbro = excelInput[excelInput$covv_virus_name %in% names(fastaInput), ]
  `%notin%` <- Negate(`%in%`)
  eybro = fastaInput[names(fastaInput) %notin% excelInput$covv_virus_name]

  ### Check if all FASTA files have metadata and writes an error log if not
  if (length(woahbro$covv_virus_name) != length(fastaInput)) {
    error = "Critical warning: Metadata entries are missing and do not correspond to FASTA files.
      The missing files are the following:\n"
    missingFiles = paste(unlist(names(eybro)), collapse = "\n")
    logName = paste0(outputPath, "outputFailure_", date, "_", time, ".log")
    logMsg = paste0(error, missingFiles)
    write(logMsg, file = logName, append = TRUE)
    error = paste0(error, missingFiles)
    failure = TRUE
    stop(error)
  }

  ### Attach header, rows, and sanitize data
  woahbro$fn = basename(fastaPath)
  woahbro = woahbro %>% arrange(covv_virus_name)
  finalExcel = rbind(excelHeader, woahbro)

  ### Write each data table to excel file
  writeSuccess = WriteXLS(finalExcel, ExcelFileName = paste0(outputPath, "/", basename(outFile)), Encoding = c("UTF-8"), AllText = TRUE, SheetNames = "Submissions")
  # Encoding = c("UTF-8", "latin1", "cp1252")


  ### CSV files require extra scrubbing;
  ### Removes last two columns, first row, and replaces empty with "unknown" data.
  ### Quotes all fields that potentially contain commas
  finalCSV = finalExcel[, -28:-29]
  finalCSV$covv_orig_lab[-1] = paste0("\"", finalCSV$covv_orig_lab[-1], "\"")
  finalCSV$covv_orig_lab_addr[-1] = paste0("\"", finalCSV$covv_orig_lab_addr[-1], "\"")
  finalCSV$covv_subm_lab_addr[-1] = paste0("\"", finalCSV$covv_subm_lab_addr[-1], "\"")
  finalCSV$covv_authors[-1] = paste0("\"", finalCSV$covv_authors[-1], "\"")
  finalCSV = finalCSV[-1,]
  csvOutFile = gsub(".xls", ".csv", outFile)


  ### Write each data table to csv file
  writeSuccess = write.csv(finalCSV, file = paste0(outputPath, "/", basename(csvOutFile)), na = "", row.names = FALSE, quote = FALSE)

  ### Write log if no errors encountered
  if (!failure & verbose) {
    logName = paste0(outputPath, "outputSuccess_", date, "_", time, ".log")
    logMsg = "Successfully minimized metadata file!"
    write(logMsg, file = logName, append = TRUE)
  }

  ### Console print
  if (verbose) {
    print(logMsg)
  }
}


### Command line argument with default parameters
args = commandArgs(trailingOnly=TRUE)

if (is.na(args[1])) {
  args = Paths()
}


### Irregular length of args is expected from malformed user input
if (length(args) < 2 & length(args) > 0) {
  warning("Incorrect number of arguments, expected 3: metadata, fasta, output path.")
}


### To make command line args easier for common file access
# Regex checks for either "def" or "default"
if (grepl("^def(ault)?$", args[1], ignore.case = TRUE))
  args[1] = Paths()[1]

if (grepl("^def(ault)?$", args[2], ignore.case = TRUE))
  args[2] = Paths()[2]

if (grepl("^def(ault)?$", args[3], ignore.case = TRUE))
  args[3] = Paths()[3]


### Read all FASTA files in provided directory and run the script
listFiles = list.files(path = args[2], pattern = ".*\\.fa$", all.files = TRUE, full.names = TRUE, recursive = FALSE)
done = lapply(listFiles, function(x) {
  SeekMetadata(metadataPath = args[1], fastaPath = x, outputPath = args[3])
})


### Clear variables from global environment
suppressWarnings(rm(packages, installedPackages, defaultFastaPath, defaultMetadataPath, defaultOutputPath, listFiles, args, SeekMetadata, done))

