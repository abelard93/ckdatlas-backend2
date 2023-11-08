#' Add ChEBI - Chemical Entities of Biological Interest  
#'
#' @description
#' Reads .csv file in 'data > ChEBI' containing HMDB db and do:
#' * create ChEBI node
#' * connect node to source
#' * create SA and PA mapping/nodes (PA: primary accession, SA: secondary accession)
#' * add ChEBI ontology (hierarchical relationships between nodes)
#' * add ChEBI to hmdb mapping

addChEBI <- function(){
  
  for (file in getFiles("metChEBI")){
    # read info on dataset
    source <- getInfo(file,c("sid","file_name","node_type"),err)
    # skip file if no/incomplete info
    if (!all(is.na(source))){
      
      # check if node for source is present, if not create
      checkSource(source,log)
      
      # create ChEBI node
      tic("Added ChEBI info")
      query <- paste0("USING PERIODIC COMMIT 1000
                        LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS props 
                        MERGE (n:metChEBI {sid:props.sid})
                        on create SET n += props
                        ")
      runQuery(query, param = TRUE, periodic = TRUE, param_list=list(file=file))
      
      # connect node to source
      query <- paste0("USING PERIODIC COMMIT 1000
                        LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS props 
                        MATCH (n:metChEBI {sid:props.sid})
                        MATCH (s:source {sid:$properties.source}) 
                        MERGE (n)-[:FROM]->(s)
                        ") 
      runQuery(query, param = TRUE, periodic = TRUE, param_list=list(file=file,source=source["sid"]$sid))
      info(log,capture.output(toc()))
      
      # create SA and PA mapping/nodes
      tic("Added ChEBI sa and pa info")
      query <- paste0(" 
                      MATCH (n:metChEBI) 
                      with n.sid as pa
                      MERGE (:accChEBI {sid:pa})
                      ")
      runQuery(query)

      query <- paste0(" 
                        MATCH (n:metChEBI) 
                        with n as ChEBI, n.sid as pa
                        MATCH (p:accChEBI {sid:pa})
                        MERGE (ChEBI)-[:PA]->(p)
                        ")
      runQuery(query)
      
      query <- paste0(" 
                      MATCH (n:metChEBI) 
                      with n as ChEBI, split(n.alt_id,'|') as sa
                      UNWIND sa as secAcc
                      MERGE (s:accChEBI {sid:secAcc}) 
                      CREATE (ChEBI)-[:SA]->(s) 
                      ")
      runQuery(query)
      info(log,capture.output(toc()))
      
      # check which secondary accessions (SA) are mapped to multiple primary accessions (PA)
      checkMappingChEBI(err) 
      
      # add ChEBI ontology
      tic("Added ChEBI relationships")
      relationships <-
        c(
          "is_a",
          "is_conjugate_base_of",
          "has_part",
          "has_role",
          "is_enantiomer_of",
          "has_functional_parent",
          "has_parent_hydride",
          "is_conjugate_acid_of",
          "is_tautomer_of",
          "is_substituent_group_from"
        )
      # loop through relations
      for(rel_type in relationships){
        query <- paste0("
                          MATCH (n:metChEBI)
                          with n as ChEBI, split(n.",rel_type,",'|') as rel
                          UNWIND rel as relNode
                          MERGE (r:metChEBI {sid:relNode})
                          MERGE (ChEBI)-[:",rel_type," {source:\"",source["file_name"],"\"}]->(r)
                          ")
        runQuery(query)
      }
      info(log,capture.output(toc()))
    }
  }
  
  # add ChEBI to hmdb mapping
  tic("Added ChEBI-HMDB info")
  query <- paste0("MATCH (n:metHMDB)
                    where exists(n.chebi_id)
                    with n,n.chebi_id as id
                    MATCH (:accChEBI {sid:\"CHEBI:\"+id})-[:PA|SA]-(c:metChEBI)
                    MERGE (n)-[:MAP]->(c)
                    ")
  runQuery(query)
  info(log,capture.output(toc()))
}
