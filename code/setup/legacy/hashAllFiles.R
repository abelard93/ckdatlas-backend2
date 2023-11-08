# Generate a combined md5 hash sum of all files in all subfolders of a list of directories

require(tools)
require(digest)

hashAllFiles <- function(lst) { # lst is list of directories, where all subfolders will be parsed each
  # for each directory, list all files, hash, and merge into one large string
  allhashes <- sapply(lst, function(d){
    paste0(md5sum(list.files(d, recursive = TRUE, full.names=T)),collapse='')
  })
  # again paste into one and hash again
  digest(paste0(allhashes,collapse=''),algo="md5")
}

# example call: 
# hashAllFiles( c(sysdiab.makepath("packages"), sysdiab.makepath("R")) )
