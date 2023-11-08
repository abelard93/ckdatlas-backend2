#' Add sources to neo4j DB
#'
#' @description
#' Function loops through files in 'data > sources' and creates node for each file
#' Folder structure corresponds to node types: 
#' * sampleType
#' * organism
#' * cohort
#' Files are .csv tables which must contain `sid`,unique source id, and `node_type` as row name

addSources <- function(){
  
  for (file in getFiles("source")){
    # read file with info on dataset
    source <- read.csv(file,header = F) 
    source <- setNames(as.character(source$V2), source$V1)
    source["file_name"] <- basename(file)
    source <- source %>% t %>% data.frame(.,stringsAsFactors = F)
    # skip file if no/incomplete info
    if (all(c("sid","node_type")%in%names(source))){
      tic(paste0("Added source file ", source["file_name"]))
      # add as node to db (sid)
      query <- paste("UNWIND $properties AS props 
                     MERGE (n:",source$node_type," {sid:props.sid}) 
                     SET n += props 
                     ",sep="")
      
      runQuery(query, param=T, param_list=source)
      info(log,capture.output(toc()))
     }
}
}
