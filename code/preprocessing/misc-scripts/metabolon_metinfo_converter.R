###########################################
# Helper script to convert metabolon files
# adapted from krumsieklab
###########################################
zap()
#### init ----
setwd("/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/metMetabolon/")
library(readxl)
source(codes.makepath("snippets/data/parseMetabolonFile.R"))
source(codes.makepath('snippets/data/writexlsx_append.R'))

onlyFirstSheet = T # restrict to first valid sheet?
info = T # do xlsx files already have info sheet (first sheet)


#### read files from all over my HDD, listed in text file ----

# load list of files
#files <- read.csv("metinfo_list.txt", header=F, sep = "\t", stringsAsFactors=F)[,1]
files <- list.files(".",full.names = T,recursive = T)
# for each file, for each sheet
for (i in 1:length(files)) {
  # skip if there is a #
  if (!(substr(files[i],0,1)=="#" || nchar(files[i])==0)) {
    # construct file name
    newfile <- paste0(
      tools::file_path_sans_ext(files[i]),
      "_metinfo.",
      tools::file_ext(basename(files[i]))
    )
    unlink(newfile) # ensure that file does not exist
    # loop over sheets
    sheets <- excel_sheets(files[i])
    success=F
    if(info){m=2}else{m=1}
    for (j in m:length(sheets)) {
      # skip if we already successfully loaded one and are supposed to only read one
      if (!(success && onlyFirstSheet)) {
        s <- NA
        tryCatch(
          s <- parseMetabolonFile(files[i], sheet=sheets[j]),
          error = function(e) {
            print(sprintf("Warning, could not identify Metabolon format in sheet '%s', file '%s'", sheets[j], basename(files[i])))
          }
        )
        if (!all(is.na(s$metinfo))) {
          # write out metinfo, if there is one
          
          if(info){
            x <- read_excel(path = files[i],sheet = 1)
            writexlsx_append(
              x = x,
              file = newfile,
              sheetName = sheets[1],
              rowNames = F
            )
            
          }
          
          writexlsx_append(
            x = s$metinfo,
            file = newfile,
            sheetName = sheets[j],
            rowNames = F
          )
          
          success=T  
        }
      }
    }
  }
}



#### parse Karsten's files ----

xlsfile = "data/Karsten MetabolonHeaders 2018-07-10.xlsx"

# get name map
dfnames <- read_excel(xlsfile, "LEGEND")

# read one by one
for (i in dfnames$SHEET) {
  # construct file name
  newfile <- paste0(
    "metinfos/automatic/",
    tools::file_path_sans_ext(basename(dfnames$File[dfnames$SHEET==i])),
    "_metinfo.xlsx"
  )
  # read actual data
  df <- read_excel(xlsfile, sprintf("Sheet%d",i), col_names=F)
  # find first row with header
  ifirst <- min(which(!is.na(df[[1]])))
  # extract data, set header
  dfout <- df[(ifirst+1):nrow(df),]
  colnames(dfout) <- df[ifirst,]
  # write
  stopifnot(!file.exists(newfile))
  write.xlsx(dfout, file=newfile)
}




       