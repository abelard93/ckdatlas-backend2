git <- "~/Documents/ICB/neo4j/" # path to git repository 
log_path <- "~/Desktop/"
db_addr <- "http://dzdcon1:9494/db/data/"
## only required if Neo4j is running on different server
neo4j_server <- "ssh alex@dzdcon1"
wd_path_server <- "/home/alex/sysdiab/"

#################
## SETUP
#################

## SET WORKING DIRECTORY
wd_path <- paste0(git,"data") 
setwd(wd_path)


checkSource <- function(source,db,log){
  # Check of the source node already exists and create it if not. 
  #
  # Args:
  #   source: List of source information.
  #   db: Graph object
  #   log: Log file.
  query <-  paste("UNWIND {properties} AS source MERGE (n:source {sid:source.sid}) set n=source",sep = "")
  logit(cypher(db,query,properties=source),log)
}
getInfo <- function(file, needed, log) {
  # Reads in the info sheet for given file and checks
  # if all required fields are present.
  #
  # Args:
  #   file: File path.
  #   needed: Information "fields" that are required for this node type
  #           as specified in requirements.txt. 
  #   log: Error log.
  #
  # Returns:
  #   A list containing all the information.
  
  source <- NA
  file.type <- tools::file_ext(file)
  
  infoWarning <- function(file, log) {
    # Print error message to log if info sheet is missing.
    msg.txt <-
      paste0(
        "Info sheet for file ",
        basename(file),
        " does not have a matching info file/sheet as specified in requirements.txt. File was skipped."
      )
    msg <- paste(Sys.time(), "NA", "INFO", msg.txt, sep = "\t")
    cat(msg, file = log, sep = "\n", append = T)
  }
  
  ## CSV FILE
  if (file.type == "csv") {
    file.info <- gsub(".csv", "_info.csv", file)
    if (file.exists(file.info)) { # does info file exist?
      tmp <- read.csv(file.info, header = F)
      source <- setNames(as.list(tmp$V2), tmp$V1)
    } else{
      # error - no info file
      infoWarning(file, log)
    }
    ## XML
  } else if (file.type=="xml"){
    file=gsub(".xml","_info.csv",file)
    if (file.exists(file)){
      tmp <- read.csv(file,header = F)
      source <- setNames(as.list(tmp$V2), tmp$V1)
    } else {
      # error - no info file
      infoWarning(file, log)
    }
    ## EXCEL FILE
  } else if (file.type == "xlsx" | file.type == "xls") {
    if ("info" %in% excel_sheets(path = file)) { # does info sheet exist?
      tmp <- read_excel(path = file, sheet = "info", col_names = F)
      source <- setNames(as.list(tmp$X__2), tmp$X__1)
    } else{
      # error - no info file
      infoWarning(file, log)
    }
  }
  
  if (!all(needed %in% names(source))) {
    # error - not all required fields 
    msg.txt <-
      paste0(
        "Info sheet for file ",
        source["file_name"],
        " not as specified in requirements.txt. File was skipped."
      )
    msg <- paste(Sys.time(), source["node_type"], "INFO", msg.txt, sep = "\t")
    cat(msg, file = log, sep = "\n", append = T)
    # set source to NA
    source <-  NA
  }
  
  # convert to character
  source <- lapply(source, as.character)
  
  return(source)
  
}

getFiles <- function(dir,pat=""){
  # Get all files in directory excluding _info files.
  #
  # Args:
  #   file: Directory to search.
  #   pat: Regular expression for file names to search for.
  #
  # Returns: 
  #   List of files,
  if(pat!=""){ # only find files according to "pat"
    tmp <- grep(list.files(dir,full.names = T,recursive = T,pattern = pat), pattern='_info.csv', inv=T, value=T)
  }else{ # find all files
    tmp <- grep(list.files(dir,full.names = T,recursive = T), pattern='_info.csv', inv=T, value=T)
  }
  
  return(tmp)
}

create_sid_reBiGG <- function(react,source){
  
  react$sid <- with(react, paste0( rxns, "_",  source))
  return(react)
  
}

create_sid_metBiGG <- function(met,source){
  
  met$sid <- with(met, paste0( mets, "_",  source))
  names(met)[which(names(met)=="id")] <- "mid"
  names(met)[which(names(met)=="mets")] <- "id"
  return(met)
  
}



rename_sid <- function(data,col){
  colnames(data)[which(colnames(data)==col)]  <-  "sid"
  return(data)
}

readBiGGModel <- function(file){
  
  needed=c("info","metabolites","S","reactions","rxnGeneMat","genes")
  sheetnames <- readxl::excel_sheets(file)
  
  if(all(needed %in% sheetnames)){
    
    M = suppressWarnings(read_excel(path=file, sheet="metabolites")) #metabolites
    M = read_excel(path=file, sheet="metabolites",col_types = c(rep("text",ncol(M))))
    if(any(grepl("\\[", M$mets))|any(grepl("\\]", M$mets))){
      M$mets <- gsub('\\[', "_", M$mets)
      M$mets <- gsub('\\]', "", M$mets)
    }
    M <- separate(M,mets,c("id","comp"),sep = -2, remove=F) 
    M <- separate(M,comp,c("no","comp"),sep = 1)
    M <- select(M,-no)
    M$id <- gsub('_\\>', "", M$id)
    
    
    
    rawS = read_excel(path=file, sheet="S")
    R = suppressWarnings(read_excel(path=file, sheet="reactions"))
    R = read_excel(path=file, sheet="reactions",col_types = c(rep("text",ncol(R))))
    
    rawG = read_excel(path=file, sheet="rxnGeneMat")
    G = suppressWarnings(read_excel(path=file, sheet="genes")) #genes
    if(is.numeric(G$genes)){
      if(ncol(G)>1){
        G = read_excel(path=file, sheet="genes",col_types = c("numeric","text")) #genes
      }
      G$genes=unlist(lapply(G$genes,as.integer))
    }else if(ncol(G)>1){
      G = read_excel(path=file, sheet="genes",col_types = c("text","text")) #genes
    }  
    
    return(list(M,G,R,rawS,rawG))
    
  }else{
    
    msg=paste(Sys.time(),"bigModel","INFO",paste("Info sheet for file ",basename(file)," not as specified in requirements.txt. File was skipped.",sep = ""),sep = "\t")
   #cat(msg,file = log,sep="\n",append = T)
    
    return(NA)
  }
}


removeBiGGModel <- function(){
  for (x in getFiles("BiGGModel",".csv")){
    system(paste0("git rm ",x))
  }
  system("git commit -m 'deleted tmp files'")
  system("git push")
}

checkBiGGModel <- function(path){
  
  files=getFiles("BiGGModel",".csv")
  if(all(c(paste0(path,"_genes.csv"),paste0(path,"_genes.csv"),paste0(path,"_reactions.csv"))%in%files)){
    return(TRUE)
  }else{
    return(FALSE)
  }
  
}

write_met_reac <- function(M,R,rawS){
  
  metabo=c()
  reac=c()
  stoech=c()
  
  for (i in 1:dim(rawS)[1]){
    x <- rawS[i,]
    metabo=c(metabo,M[x$row])
    reac=c(reac,R[x$col])
    stoech=c(stoech,x$val)
  }
  
  result <- cbind(metBiGG=metabo,reBiGG=reac,stochiometry=stoech)
  result <- as.data.frame(result) %>% group_by(metBiGG,stochiometry) %>% summarise(reBiGG  = paste(reBiGG, collapse ="|"))
  write.csv(result,file = paste(tmp_path,"_met_reac.csv",sep = ""),na = "",row.names = F)
  
  
}


write_gene_reac <- function(G,R,rawG){
  
  reac=c()
  gene=c()
  
  for (i in 1:dim(rawG)[1]){
    x <- rawG[i,]
    #print(G[x$col])
    if(!is.na(G[x$col])){
      gene=c(gene,G[x$col])
      reac=c(reac,R[x$row])
    }
  }
  
  
  result <- cbind(gene=gene,reBiGG=reac)
  result <- aggregate(reBiGG ~ gene, data = result, FUN = paste, collapse = "|")
  
  write.csv(result,file = paste(tmp_path,"_gene_reac.csv",sep = ""),na = "",row.names = F)
  
  
}


tmp_path <- tools::file_path_sans_ext(file)
source <- getInfo(file, c("sid", "tax", "taxID", "file_name"), err)
#checkSource(source, db, log)
if(!checkBiGGModel(tmp_path)){
  
  tmp <- readBiGGModel(file) ##error
  
  if(!all(is.na(tmp))){
    
    M <- tmp[[1]]
    M <- create_sid_metBiGG(M,source["sid"])
    G <- data.frame(tmp[[2]])
    G <- rename_sid(G,"genes")
    R <- tmp[[3]]
    R <- create_sid_reBiGG(R,source["sid"])
    rawS <- tmp[[4]]
    rawG <- tmp[[5]]
    
    
    write.csv(M,file=paste(tmp_path,"_metabolites.csv",sep = ""),na = "",row.names = F)
    
    if(source["gene"]=="symbol"){
      func <- function(x){ ifelse(is.na(as.numeric(x)), x, NA)}
      G$sid <- suppressWarnings(unlist(lapply(G$sid,func)))
      write.csv(data.frame(sid=G[!is.na(G$sid),]),file=paste(tmp_path,"_genes.csv",sep = ""),na = "",row.names = F)
    }else if(source["gene"]=="ncbi"){
      write.csv(G,file=paste(tmp_path,"_genes.csv",sep = ""),na = "",row.names = F)
      
    }
    
    write.csv(R,file=paste(tmp_path,"_reactions.csv",sep = ""),na = "",row.names = F)
    
    write_met_reac(M$sid,R$sid,rawS)
    write_gene_reac(G$sid,R$sid,rawG)
    
    system(paste("git add ",wd_path,"/",dirname(file),"/*.csv",sep = ""))
    system("git commit -m 'added tmp files'")
    system("git push")
   # system(paste0(neo4j_server," '(cd ",wd_path_server,"; git pull)'"))
  }
}
