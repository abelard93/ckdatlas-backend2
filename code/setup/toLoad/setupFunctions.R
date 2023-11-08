#  Helper functions  ----------------------------------------

#' Create log file
#'
#' Create a new logger object with `create.logger()`.
#'
#' @param log_path Path to log File

create_log <- function(log_path) {
  log <- create.logger()
  # Set the logger's file output.
  logfile(log) <- log_path
  # Set logger level.
  level(log) <- 'INFO'
  
  return(log)
  
}

#' Print to log
#'
#' Catch warnings/errors and print to log.
#'
#' @param operation Command to perform
#' @param log Path to log File

logit <- function(operation, log) {
  tryCatch(
    operation,
    warning = function(w) {
      warn(log, conditionMessage(w))
      
    },
    error = function(e) {
      error(log, conditionMessage(e))
      
      print("FATAL ERROR: please consult neo4j.log for more info")
      quit()
    }
  )
  
}

#' Add constraints
#'
#' Add unique constraint to node type in neo4j DB.
#'
#' @param types Node type or list of node types
#' @param unique Name of unique constraint

add_constr <-  function(types, unique) {
  for (t in types) {
    runQuery(paste0("CREATE CONSTRAINT FOR (n:",
             t,
             ") REQUIRE n.",
             unique,
             " IS UNIQUE")
             )
  }
}

#' List files
#'
#' AGet all files in directory excluding _info files.
#'
#' @param dir Directory to search
#' @param pat Regular expression for file names to search for
#'
#' @return List of files

getFiles <- function(dir, pat = "") {
  if (pat != "") {
    # only find files according to "pat"
    tmp <-
      grep(
        list.files(
          dir,
          full.names = T,
          recursive = T,
          pattern = pat
        ),
        pattern = '_info.csv',
        inv = T,
        value = T
      )
  } else{
    # find all files
    tmp <-
      grep(
        list.files(dir, full.names = T, recursive = T),
        pattern = '_info.csv',
        inv = T,
        value = T
      )
  }
  return(tmp)
}

#' Get source information
#'
#' Reads in the info sheet for given file and checks if all required fields are present.
#'
#' @param file File path
#' @param needed Information "fields" that are required for this node type as specified in documentation.
#' @param log Error log
#'
#' @return A list containing all the information from the file, `NA` if requirements not met.

getInfo <- function(file, needed, log) {
  source <- NA
  file.type <- tools::file_ext(file)
  
  infoWarning <- function(file, log) {
    # Print error message to log if info sheet is missing.
    msg.txt <-
      paste0(
        "Info sheet for file ",
        basename(file),
        " does not have a matching info file/sheet as specified in documentation. File was skipped."
      )
    msg <- paste(Sys.time(), "NA", "INFO", msg.txt, sep = "\t")
    cat(msg,
        file = log,
        sep = "\n",
        append = T)
  }
  
  ## CSV FILE
  if (file.type == "csv" || file.type == "gz") {
    if (file.type == "gz") {
      file.info <- gsub(".csv.gz", "_info.csv", file)
    } else{
      file.info <- gsub(".csv", "_info.csv", file)
    }
    if (file.exists(file.info)) {
      # does info file exist?
      tmp <- read.csv(file.info, header = F)
      source <- setNames(as.list(tmp$V2), tmp$V1)
    } else{
      # error - no info file
      infoWarning(file, log)
    }
    ## XML
  } else if (file.type == "xml") {
    file = gsub(".xml", "_info.csv", file)
    if (file.exists(file)) {
      tmp <- read.csv(file, header = F)
      source <- setNames(as.list(tmp$V2), tmp$V1)
    } else {
      # error - no info file
      infoWarning(file, log)
    }
    ## EXCEL FILE
  } else if (file.type == "xlsx" | file.type == "xls") {
    if ("info" %in% excel_sheets(path = file)) {
      # does info sheet exist?
      tmp <-
        suppressMessages(read_excel(
          path = file,
          sheet = "info",
          col_names = F
        ))
      colnames(tmp) <- c("key", "value")
      source <- setNames(tmp$value, tmp$key)
    } else{
      # error - no info file
      infoWarning(file, log)
    }
  }
  
  if (!all(needed %in% names(source))) {
    # error - not all required fields
    msg.txt <-
      paste0("Info sheet for file ",
             source["file_name"],
             " not as specified in documentation. File was skipped.")
    msg <-
      paste(Sys.time(), source["node_type"], "INFO", msg.txt, sep = "\t")
    cat(msg,
        file = log,
        sep = "\n",
        append = T)
    # set source to NA
    source <-  NA
  }
  # convert to character
  source <- lapply(source, as.character)
  
  return(source)
}

#' Check if source node exists and create
#' 
#' Checks if the source node already exists and create it if not. 
#' Further checks if following are given and creates/connects to node:
#' * organism
#' * cohort - multiple entries separated by "|"
#' * sample_type - multiple entries separated by "|"
#' * publication - multiple entries separated by "|" 
#'
#' @param source List of source information - output of `getInfo()`
#' @param log Log file

checkSource <- function(source,log){
  
  ## SOURCE
  query <-  paste0("UNWIND $properties AS prop MERGE (n:source {sid:prop.sid}) ON CREATE SET n=prop")
  runQuery(query, param=TRUE, param_list=source)

  ## ORGANISM 
  if("organism"%in%names(source)){
    query <-  paste0("UNWIND $properties AS prop 
                    MATCH (s:source {sid:prop.sid})
                    MERGE (n:organism {sid:prop.organism})
                    CREATE (s)-[:FROM_ORGANISM]->(n)
                    ")
    runQuery(query,param=TRUE, param_list=source)
  }
  ## COHORT 
  if("cohort"%in%names(source)){
    query <-  paste0("UNWIND $properties AS prop 
                    with split(prop.cohort,'|') as co, prop
                    MATCH (s:source {sid:prop.sid})
                    UNWIND co as cohort
                    MERGE (n:cohort {sid:cohort})
                    CREATE (s)-[:FROM_COHORT]->(n)
                    ")
    runQuery(query,param=TRUE, param_list=source)
  }
  ## SAMPLE TYPE: E.G. SERUM, LIVER, BRAIN etc 
  if("sample_type"%in%names(source)){
    query <-  paste0("UNWIND $properties AS prop 
                    with split(prop.sample_type,'|') as st, prop
                    MATCH (s:source {sid:prop.sid})
                    UNWIND st as type
                    MERGE (n:sampleType {sid:type})
                    CREATE (s)-[:FROM_SAMPLE_TYPE]->(n)
                    ")
    runQuery(query,param=TRUE, param_list=source)
  }
  ## PUBLICATION
  if("publication"%in%names(source)){
    query <-  paste0("UNWIND $properties AS prop 
                    with split(prop.publication,'|') as pb, prop
                    MATCH (s:source {sid:prop.sid})
                    UNWIND pb as pub
                    MERGE (n:publication {sid:pub}) 
                    CREATE (s)-[:FROM_PUBLICATION]->(n)
                    ")
    runQuery(query,param=TRUE, param_list=source)
  }
}



#' Read Metabolon annotation file.
#' 
#' Parse annotation of Metabolon metabolites as sent by vendor. 
#' File format is excel and file needs to have 2 sheets, the first being named 'info' containing meta data.
#' Attention: Annotation df eeds 'COMP_ID' column or error is thrown.
#'
#' @param file Path to file
#' @param log Logger object
#'
#' @return Dataframe of Metabolon data

readMetaboFile <- function(file, log) {

  metabo <- NA
  sheetnr <- length(excel_sheets(file))
  # check if correct number of sheets and if info sheet is first.
  if (sheetnr == 2 & excel_sheets(file)[1] == "info") {
    # read metabolon file
    metabo <-  read_excel(file, sheet = 2, col_types ="text")
    # remove uninformative column
    if ("PATHWAY_SORTORDER" %in% colnames(metabo)) {
      metabo <- metabo[, -which(colnames(metabo) == "PATHWAY_SORTORDER")]
    } 
    # skip file and issue error if comp_id not given
    if (!"COMP_ID" %in% colnames(metabo)) {
      metabo <- NA
      msg <-
        paste(
          Sys.time(), "metMetabolon", "FORMAT",
          paste0("Internal id coloumn (COMP_ID) missing: ", basename(file), " skipped, please check file"),
          sep = "\t"
        )
      cat(msg,file = err,sep="\n",append = T)
    }
    # error if format requirements not met   
  } else{
    msg <-
      paste(
        Sys.time(), "metMetabolon", "FORMAT",
        paste0( basename(file)," must have more than one sheet (first must be named info). 
                Please check requirements.txt for appropriate data conventions. File was skipped.\n"),
        sep = "\t"
      )
    cat(msg,file = log,sep="\n",append = T)
  }
  return(metabo)
}


#' Check metabolite ID mapping.
#' 
#' Checks for 1:n metabolite ID mappings. Results are written to the log file.
#' 
#' @param log Log file

checkMapping  <-  function(log) {

  # 1. CHECK IF MULTIPLE METABOLITE MAP TO ONE HMDB
  # metmetabolon
  query = "
  match (h:metHMDB)-[:MAP]-(m:metOntology)-[:MAP]-(mm:metMeasured)
  with h.sid as HMDB, collect(distinct h.name) as hName, collect(distinct m.CHEMICAL_ID) as multipleID, collect(distinct mm.name) as mName
  where size(multipleID)>1
  return HMDB,hName, multipleID,mName
  "
  multipleHMDB <- runQuery(query, read = TRUE)
  # write to log
  if (length(multipleHMDB)!=0) {
    apply(multipleHMDB, 1, function(x) {
      msg <- paste(Sys.time(),"metHMDB","MAPPING",paste0(x[1]," (",x[2],")"," mapped to metOntology: ",x[3]," name: ",x[4]),sep = "\t")
      cat(msg,file = log,sep="\n",append = T)
    })
  }

  # 2. CHECK IF MULTIPLE HMDB MAP TO ONE METABOLITE
  query = "
    match (mm:metMeasured)-[:MAP]-(m:metOntology)-[:MAP]-(h:metHMDB)
    with m.sid as metabolite,collect(distinct mm.name) as mName,collect(distinct h.sid) as multipleID, collect(h.name) as hName
    where size(multipleID)>1
    return metabolite, mName, multipleID, hName
    "
  multipleHMDB <- runQuery(query, read = TRUE)
  
  # write to log
  if (length(multipleHMDB)!=0) {
    apply(multipleHMDB, 1, function(x) {
      msg <- paste(Sys.time(),"metOntology","MAPPING",paste0(x[1]," (",x[2],")"," mapped to HMDB: ",x[3]," name: ",x[4]),sep = "\t")
      cat(msg,file = log,sep="\n",append = T)
    })
  }
  # 3. CHECK IF MULTIPLE METABOLON MAP TO ONE METABOLITE
  query = "
    match (m:metMetabolon)-[:MAP]-(n:metOntology)
    with distinct m.COMP_ID as metMetabolon,collect(distinct m.name) as mName, collect(distinct n.sid) as multipleID
    where size(multipleID)>1
    return metMetabolon, mName, multipleID
    "
  multipleMet <- runQuery(query, read = TRUE)
  
  # write to log
  if (length(multipleMet)!=0) {
    apply(multipleMet, 1, function(x) {
      msg <- paste(Sys.time(),"metMetabolon","MAPPING",paste0(x[1]," (",x[2],")"," mapped to metOntology: ",x[3]),sep = "\t")
      cat(msg,file = log,sep="\n",append = T)
    })
  }
  
  # 4. CHECK IF MULTIPLE metBiocrates MAP TO ONE METABOLITE
  query = "
    match (m:metBiocrates)-[r:MAP]-(n:metOntology)
    with distinct m.sid as metBiocrates,collect(distinct n.sid) as multipleID
    where size(multipleID)>1
    return metBiocrates,multipleID
    "
  multipleMet <- runQuery(query)
  
  #write to log
  if (length(multipleMet)!=0) {
    apply(multipleMet, 1, function(x) {
      msg <- paste(Sys.time(), "metOntology", "MAPPING", paste0(x[1]," mapped to metBiocrates: ", x[2]), sep = "\t")
      cat(msg,file = log,sep="\n",append = T)
    })
  }
}

checkMappingGwas  <-  function(log){
  # Checks GWAS to metMetabolon mapping. 
  #
  # Args:
  #   db: Graph object
  #   log: Log file.  
  ## CHECK IF METABOLON NOT ATTACHED TO SOURCE NODE
  query = "
  match (m:metMetabolon:metMeasured)
  where not (m)-[:FROM]-()
  return m.sid
  "
  multipleMet <- readSimpleQuery(query)
  ## WRITE TO LOG
  if(!is.null(multipleMet)){
      apply(multipleMet,1,function(x){
        tmp="gwasHit associated with this metMetabolon but no further information availible (no source/not measured)" 
        msg=paste(Sys.time(),"metMetabolon","MAPPING",x[1],tmp,sep = "\t")
        cat(msg,file = log,sep="\n",append = T)
      })
  }
}

#' Check HMDB to metabolite mapping.
#' 
#' Checks for multiple metabolite to HMDB mappings and vice versa. Results are written to the log file.
#' 
#' @param log Log file

checkMappingHmdb  <-  function(log){

  # CHECK IF SA IS LINKED TO MULTIPLE PA
  query = "
  match (h:metHMDB)-[:PA]-(a:accHMDB)-[:SA]-(m:metHMDB)
  where h <> m
  return h.sid, m.sid,a.sid
  "
  multipleHMDB <- runQuery(query, read = TRUE)

  # write to log
  if(length(multipleHMDB)!=0){
    apply(multipleHMDB,1,function(x){
      tmp <- paste(unlist(x[2]),collapse="|")
      msg=paste(Sys.time(),"accHMDB","MAPPING",x[3],paste("PA:",x[1],"|SA:",x[2],sep = ""),sep = "\t")
      cat(msg,file = log,sep="\n",append = T)
    })
  }
}

#' Check ChEBI to metabolite mapping.
#' 
#' Checks for multiple metabolite to ChEBI mappings and vice versa. Results are written to the log file.
#' 
#' @param log Log file

checkMappingChEBI  <-  function(log){
  
  # CHECK IF SA IS LINKED TO MULTIPLE PA
  query = "
  match (h:metChEBI)-[:PA]-(a:accChEBI)-[:SA]-(m:metChEBI)
  where h <> m
  return h.sid, m.sid,a.sid
  "
  multipleChEBI <- runQuery(query, read = TRUE)
  
  # write to log
  if(length(multipleChEBI)!=0){
    apply(multipleChEBI,1,function(x){
      tmp <- paste(unlist(x[2]),collapse="|")
      msg=paste(Sys.time(),"accChEBI","MAPPING",x[3],paste("PA:",x[1],"|SA:",x[2],sep = ""),sep = "\t")
      cat(msg,file = log,sep="\n",append = T)
    })
  }
}



# runSimpleQuery <- function(query){
#  
#   session <- driver$session()
#   session$write_transaction(submitQuery,query)
#   session$close()
# }
# 
# readSimpleQuery <- function(query){
#   
#   session <- driver$session()
#   result <- session$read_transaction(submitReadQuery,query)
#   session$close()
#   result %>% lapply(function(x) lapply(x,function(y) paste0(y,collapse = "|"))) %>% rbindlist() %>% data.frame
#   
# }
# 
# submitReadQuery <- function(tx,query){
#   
#     result <- tx$run(query)
#     result$data()
#   
# }
# 
# 
# readQueryWithParam <- function(query,param){
#   
#   session <- driver$session()
#   result <- session$read_transaction(submitReadQueryWithParam,query,param)
#   session$close()
#   result %>% lapply(function(x) lapply(x,function(y) paste0(y,collapse = "|"))) %>% rbindlist() %>% data.frame
#   
# }
# 
# submitReadQueryWithParam <- function(tx,query,param){
#   
#   result <- tx$run(query, properties=param)
#   result$data()
#   
# }
# 
# submitQuery <- function(tx,query){
# 
#   logit(
#     tx$run(query),
#     log
#   )
#   
# }
# 
# runQueryWithParam <- function(query,param){
#   
#   session <- driver$session()
#   session$write_transaction(submitQueryWithParam,query,param)
#   session$close()
# }
# 
# submitQueryWithParam <- function(tx,query,param){
#   
#   logit(
#     tx$run(query, properties=param),
#     log
#   )
#   
# }
# 
# runQueryPeriodicWithParam <- function(query,param){
#   
#   session <- driver$session()
#   logit(session$run(query,properties=param),log)
#   session$close()
# }
# 
# runQueryPeriodic <- function(query){
#   
#   session <- driver$session()
#   logit(session$run(query),log)
#   session$close()
# }