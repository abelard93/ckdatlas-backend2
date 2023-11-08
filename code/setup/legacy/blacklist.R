#################
## Blacklist
#################

doBlacklisting <- function(file){
  if(missing(file)){
    warning("No blacklisting peformed. Please provide file - csv format with one sid per line")
  }else{
    
    tic(paste0("Blacklisting (marking) nodes from file ",basename(file)))
    query <- paste0("
                      LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS props
                      MATCH (n:recon3D) where n.sid=props.sid 
                      call apoc.create.addLabels(id(n), ['blacklist']) yield node RETURN \"\"
                        ")
    tmp <- cypher(db,query)
    
    info(log,capture.output(toc()))
  }
  
}
