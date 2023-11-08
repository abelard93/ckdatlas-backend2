#################
## WCM
#################

addWCM <- function(){

    for (file in getFiles("metWCM")){
      #read info on dataset
      source <- getInfo(file,c("sid","tax","taxID","file_name","node_type"),err)
      #skip file if no/incomplete info
      if (!all(is.na(source))){
        tic(paste0("Added metWCM file ", source["file_name"]))
        #check/add source node for dataset
        checkSource(source,db,log)
        
        ## ADD metWCM NODE
        # sid is source_id
        query <- paste("
                       USING PERIODIC COMMIT 1000
                       LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS props 
                       MERGE (n:metWCM:metMeasured {sid:props.id+\"_\"+\"",source["sid"],"\"}) 
                       ",sep="")
        
        cypher(db,query)
        
        query <- paste("
                       USING PERIODIC COMMIT 1000
                       LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS props 
                       MATCH (n:metWCM:metMeasured {sid:props.id+\"_\"+\"",source["sid"],"\"}) 
                       MATCH (s:source {sid:\"",source["sid"],"\"}) 
                       MERGE (n)-[:FROM]->(s) 
                       ",sep="")
        
        cypher(db,query)
        
        ## CREATE AND MAP TO DUMMY metOntology NODE [to do]
        query <- paste("
                       USING PERIODIC COMMIT 1000
                       LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS props 
                       MERGE (m:metOntology {sid:props.id})
                       SET m.sid=props.id, m.hmdb=props.hmdb with m,props
                       MATCH (n:metWCM {sid:props.id+\"_\"+\"",source["sid"],"\"})
                       MERGE (n)-[:MAP]->(m) 
                       ",sep="")
        
        cypher(db,query)
        
        ## metWCM to HMDB
        query <- paste("
                       USING PERIODIC COMMIT 1000
                       LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS props 
                       with props where exists(props.hmdb)
                       MATCH (m:metOntology {sid:props.id}) 
                       MATCH (h1:accHMDB)-[]-(h:metHMDB) where h1.sid=m.hmdb 
                       MERGE (m)-[:MAP]->(h) 
                       ",sep="")
        
        cypher(db,query) # send query
        info(log,capture.output(toc()))
                       }
    }
}