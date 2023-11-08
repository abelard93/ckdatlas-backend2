# read in metabolon file

#library(xlsx)
library(readxl)
library(stringr)

parseMetabolonFile <- function(file, sheet, copynansheet='') {
  
  # using xlsx package:
  #raw = read.xlsx(file=file, sheetName=sheet, header=F, stringsAsFactors=F)
  
  # using readxl package:
  raw = read_excel(path=file, sheet=sheet, col_names = F)
  
  result=list()
  
  # find metabolite header row and last metabolite row
  imetheader = min(which(!is.na(raw[,1])))
  imetlast = max(which(apply(is.na(raw),1,sum)<dim(raw)[2]))
  # find sample header column and last sample row
  isampheader = min(which(!is.na(raw[1,])))
  isamplast = max(which(apply(is.na(raw),2,sum)<dim(raw)[1]))
  
  # fix overlapping cell
  overl=str_replace(gsub("\\s+", " ", str_trim(raw[imetheader,isampheader])), "B", "b")
  overl=strsplit(overl, " ")[[1]]
  overlmet = overl[2]
  overlsamp = overl[1]
  
  
  # extract metabolite information
  result$metinfo = raw[(imetheader+1):imetlast, 1:isampheader]
  colnames(result$metinfo) = raw[imetheader, 1:isampheader]
  rownames(result$metinfo) = c()
  # fix last one
  colnames(result$metinfo)[dim(result$metinfo)[2]] = overlmet
  # also add Metabolon ID in M00000 format
  if("COMP_ID"%in%colnames(result$metinfo)){
    result$metinfo$COMP_IDstr = sapply(sprintf('%05d',as.numeric(result$metinfo$COMP_ID)), function(x)paste0('M',x))
  }
  
  
  # extract sample information
  result$sampleinfo = data.frame(t(raw[1:imetheader-1,(isampheader+1):isamplast]))
  colnames(result$sampleinfo) = as.list(raw[1:imetheader-1,isampheader])[[1]] # dirty hack, had something to do with the output format of read_excel
  rownames(result$sampleinfo) = c()
  
  
  # extract data
  d=t(raw[(imetheader+1):imetlast, (isampheader+1):isamplast])
  result$data = as.data.frame(apply(d,2, as.numeric))
  
  # as column names, use "BIOCHEMICAL", if available
  if ("BIOCHEMICAL" %in% colnames(result$metinfo)) {
    colnames(result$data) = result$metinfo$BIOCHEMICAL
  } else {
    colnames(result$data) = c()
  }
  
  # as row names, use "SAMPLE_NAME", if available 
  if ("SAMPLE_NAME" %in% colnames(result$sampleinfo)) {
    rownames(result$data) = result$sampleinfo$SAMPLE_NAME
  } else {
    rownames(result$data) = c()
  }
  
  # set info flags
  result$info$file = file
  result$info$sheet = sheet
  
  # copy NanN from another sheet?
  if (nchar(copynansheet)>0) {
  # recursive call
    nandf = parseMetabolonFile(file=file, sheet=copynansheet)
    # sanity checks
    if (!(all.equal(colnames(result$data), colnames(nandf$data)))==T)
      stop('some metabolites are different between data sheet and NaN sheet');
    if (!all.equal(rownames(result$data), rownames(nandf$data)))
      stop('some sample names are different between data sheet and NaN sheet');
    # copy over NaNs
    result$data[is.na(nandf$data)]=NA
    # set info flag
    result$info$copynansheet = copynansheet

  }
  # return
  result
}


# test code from jan
if (F) { # never execute automatically
  file = 'C:/Users/Jan/work/metabo/highfatgut/Test samples non-targeted metabolome/NonTargeted_Test sample_results_HELM-13-14ML+ CDT_150414.xlsx'
  sheet = 'OrigScale'
  
  D = parseMetabolonFile(file,sheet)
}
