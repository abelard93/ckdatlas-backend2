#################
## gwasHit
#################

addGwasHit <- function(){
  
  # create index on Symbol for faster search in functions
  query <- "CREATE INDEX ON :gwasHit(rsID)"
  cypher(db,query)
  
  for(file in getFiles("gwasHit")){
    # read info on dataset
    source <- getInfo(file,c("sid","tax","taxID","file_name","node_type"),err)
    # skip file if no/incomplete info
    if (!all(is.na(source))){
      # check if node for source is present, if not create
      checkSource(source,db,log)
      tic(paste0("Added gwasHits ",source["file_name"]))
      ## ADD GWAS HIT AND CONNECT TO METABOLON and SNP
      query <- paste("
                    USING PERIODIC COMMIT 1000 
                    LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS props 
                    CREATE (n:gwasHit {pvalue:props.pvalue, rsID:props.rsID, metabolite:replace(props.metabolite,\"M\",\"\")}) with n
                    MATCH (s:source {sid:\"",source["sid"],"\"}) with n,s
                    MERGE (m:metMetabolon:metMeasured {COMP_ID:n.metabolite}) 
                    CREATE (n)-[:FROM]->(s)
                    CREATE (m)-[:HIT]->(n)
                    ",sep="")
      cypher(db,query)
      
      ## TODOOOO
      ## ADD SNP AND CONNECT TO GWAS HIT
      query <- paste("USING PERIODIC COMMIT 1000 
                     LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS props  
                     CREATE (n:SNP {sid:props.sid, chr:props.chr, position:props.position}) with n 
                     MERGE (g:gwasHit {rsID:n.sid})
                     CREATE (n)-[:HIT]->(g)
                     ",sep="")
      cypher(db,query)
      
      
      info(log,capture.output(toc()))
    }
  }
}