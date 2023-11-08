#' Add metabolites measured on a Metabolon platform  
#'
#' @description
#' 1 - Reads Metabolon specific results files in 'data > metMetabolon' in excel format:
#' * create :metMetabolon:metMeasured node
#' * connect node to source
#' * create metaInfo nodes for specific source properties
#' * add hmdb to smpdb mapping
#' 2 - Create pseudo ontology layer (using chemical id)
#' 3 - Reads .csv file with meta Information on Metabolon metabolites.

addMetabolon <- function(){
  
  for (file in getFiles("metMetabolon","_metinfo")) {
    # read file
    metabo <- readMetaboFile(file, err)
    # check that nothing went wrong/format errors
    if (!all(is.na(metabo))) {
      # read info on dataset
      source <- getInfo(file, c("sid", "tax", "taxID", "file_name", "node_type"), err)
      
      # skip file if no/incomplete info
      if (!all(is.na(source))) {
        
        tic(paste0("Added metMetabolon file ", source["file_name"]))
        
        # generate dataset specific sid and id
        metabo$sid <- with(metabo, paste0(COMP_ID, "_", tolower( gsub("[^A-Za-z]|[[:space:]]", "", PLATFORM) )))
      
        # check if node for source is present, if not create
        checkSource(source, log)
        
        metabo$source <- source["sid"]$sid
        
        # create :metMetabolon:metMeasured node
        metabo_neo <- unname(split(metabo, 1:nrow(metabo))) %>% lapply(function(x){x %>% select_if(~ !any(is.na(.)))})
        query <-  paste0("
                           UNWIND $properties AS prop 
                           MERGE (m:metMetabolon:metMeasured {sid:prop.sid}) 
                           ON CREATE SET m.COMP_ID=prop.COMP_ID, m.PLATFORM=prop.PLATFORM, m.name=prop.BIOCHEMICAL
                           ")
        
        runQuery(query,param = TRUE, param_list = metabo_neo)
        
        # connect node to source
        query <-  paste0("
                         UNWIND $properties AS prop 
                         MATCH (m:metMetabolon:metMeasured {sid:prop.sid}) 
                         MATCH (s:source {sid:prop.source}) 
                         MERGE (m)-[:FROM]->(s)")
        runQuery(query,param = TRUE, param_list = metabo_neo)
        
        # create metaInfo nodes for following properties:
        # BIOCHEMICHAL, SUPER_PATHWAY, SUB_PATHWAY, PLATFORM, MASS, CAS, PUBCHEM, KEGG, HMDB, RI
     
        ## BIOCHEMICAL
        query <-  paste0("
                         UNWIND $properties AS prop with prop where exists(prop.BIOCHEMICAL)
                         MATCH (m:metMetabolon:metMeasured {sid:prop.sid}) 
                         MATCH (s:source {sid:prop.source}) 
                         MERGE (n:metName {sid:apoc.text.replace(toLower(prop.BIOCHEMICAL), '[^A-Za-z0-9]', '')}) with m,n,s
                         MERGE (m)-[:HAS_NAME]->(n) with n,s
                         MERGE (n)-[:FROM]->(s)
                         ")
        runQuery(query,param = TRUE, param_list = metabo_neo)
        
        query <-  paste0("
                         UNWIND $properties AS prop with prop where exists(prop.BIOCHEMICAL)
                         MATCH (m:metMetabolon:metMeasured {sid:prop.sid}) 
                         MATCH (s:source {sid:prop.source}) 
                         MERGE (n:metName {sid:prop.BIOCHEMICAL}) with m,n,s
                         MERGE (m)-[:HAS_NAME]->(n) with n,s
                         MERGE (n)-[:FROM]->(s)
                         ")
        runQuery(query,param = TRUE, param_list = metabo_neo)
        
        ## SUPER_PATHWAY
        query <-  paste0("
                         UNWIND $properties AS prop with prop where exists(prop.SUPER_PATHWAY)
                         MATCH (m:metMetabolon:metMeasured {sid:prop.sid}) 
                         MATCH (s:source {sid:prop.source}) 
                         MERGE (n:metSuperPathway {sid:prop.SUPER_PATHWAY}) with m,n,s
                         MERGE (m)-[:IN_PATHWAY]->(n) with n,s
                         MERGE (n)-[:FROM]->(s)
                         ")
        runQuery(query,param = TRUE, param_list = metabo_neo)
      
        ## SUB_PATHWAY
        query <-  paste0("
                         UNWIND $properties AS prop with prop where exists(prop.SUB_PATHWAY)
                         MATCH (m:metMetabolon:metMeasured {sid:prop.sid}) 
                         MATCH (s:source {sid:prop.source}) 
                         MERGE (n:metSubPathway {sid:prop.SUB_PATHWAY}) with m,n,s
                         MERGE (m)-[:IN_PATHWAY]->(n) with n,s
                         MERGE (n)-[:FROM]->(s)
                         ")
        runQuery(query,param = TRUE, param_list = metabo_neo)
       
        # Platform
        query <-  paste0("
                         UNWIND $properties AS prop with prop where exists(prop.PLATFORM)
                         MATCH (m:metMetabolon:metMeasured {sid:prop.sid})
                         MATCH (s:source {sid:prop.source})
                         MERGE (n:metPlatform {sid:apoc.text.replace(toLower(prop.PLATFORM), '[^A-Za-z0-9]', '')}) with m,n,s
                         MERGE (m)-[:MEASURED_WITH]->(n) with n,s
                         MERGE (n)-[:FROM]->(s)
                         ")
        runQuery(query,param = TRUE, param_list = metabo_neo)
        
        ## CAS
        query <-  paste0("
                         UNWIND $properties AS prop with prop where exists(prop.CAS)
                         MATCH (m:metMetabolon:metMeasured {sid:prop.sid}) 
                         MATCH (s:source {sid:prop.source})
                         MERGE (n:metCAS {sid:prop.CAS}) with m,n,s
                         MERGE (m)-[:HAS_ID]->(n) with n,s
                         MERGE (n)-[:FROM]->(s)
                         ")
        runQuery(query,param = TRUE, param_list = metabo_neo)
      
        ## MASS
        query <-  paste0("
                         UNWIND $properties AS prop with prop where exists(prop.MASS)
                         MATCH (m:metMetabolon:metMeasured {sid:prop.sid}) 
                         MATCH (s:source {sid:prop.source})
                         MERGE (n:metWeight {sid:prop.MASS}) with m,n,s
                         MERGE (m)-[:HAS_WEIGHT]->(n) with n,s
                         MERGE (n)-[:FROM]->(s)
                         ")
        runQuery(query,param = TRUE, param_list = metabo_neo)
       
        ## RI
        query <-  paste0("
                         UNWIND $properties AS prop with prop where exists(prop.RI)
                         MATCH (m:metMetabolon:metMeasured {sid:prop.sid}) 
                         MATCH (s:source {sid:prop.source})
                         MERGE (n:metRI {sid:prop.RI}) with m,n,s
                         MERGE (m)-[:HAS_RI]->(n) with n,s
                         MERGE (n)-[:FROM]->(s)
                         ")
        runQuery(query,param = TRUE, param_list = metabo_neo)
      
        ## CHEMSPIDER
        query <-  paste0("
                         UNWIND $properties AS prop with prop where exists(prop.CHEMSPIDER)
                         MATCH (m:metMetabolon:metMeasured {sid:prop.sid}) 
                         MATCH (s:source {sid:prop.source})
                         MERGE (n:metChemSpider {sid:prop.CHEMSPIDER}) with m,n,s
                         MERGE (m)-[:HAS_ID]->(n) with n,s
                         MERGE (n)-[:FROM]->(s)
                         ")
        runQuery(query,param = TRUE, param_list = metabo_neo)
               
        ## PubChem
        query <-  paste0("
                         UNWIND $properties AS prop with prop where exists(prop.PUBCHEM)
                         MATCH (m:metMetabolon:metMeasured {sid:prop.sid}) 
                         MATCH (s:source {sid:prop.source})
                         MERGE (n:metPubChem {sid:prop.PUBCHEM}) with m,n,s
                         MERGE (m)-[:HAS_ID]->(n) with n,s
                         MERGE (n)-[:FROM]->(s)
                         ")
        runQuery(query,param = TRUE, param_list = metabo_neo)
    
        
     
        ## KEGG
        query <-  paste0("
                         UNWIND $properties AS prop with prop where exists(prop.KEGG)
                         MATCH (m:metMetabolon:metMeasured {sid:prop.sid}) 
                         MATCH (s:source {sid:prop.source})
                         MERGE (n:metKEGG {sid:prop.KEGG}) with m,n,s
                         MERGE (m)-[:HAS_ID]->(n) with n,s
                         MERGE (n)-[:FROM]->(s)
                         ")
        runQuery(query,param = TRUE, param_list = metabo_neo)
       
     
        ## HMDB
        query <-  paste0("
                         UNWIND $properties AS prop with prop where exists(prop.HMDb_ID)
                         MATCH (m:metMetabolon:metMeasured {sid:prop.sid}) 
                         MATCH (h1:accHMDB)-[]-(h:metHMDB) where h1.sid=prop.HMDb_ID with m,h
                         MERGE (m)-[:HAS_ID]->(h) 
                         ")
        runQuery(query,param = TRUE, param_list = metabo_neo)
        
        # 2 - Create pseudo ontology layer using Chemical_id of metbabolon metabolites
        
        if ("CHEMICAL_ID" %in% colnames(metabo)) {
          
          metabo$CHEMICAL_ID[metabo$'CHEMICAL_ID'=="0"] <- NA
          data <- metabo[, c("CHEMICAL_ID", "BIOCHEMICAL", "HMDb_ID")]
          # use biochemical name if no chemical_id given
          data$sid <-
            apply(data, 1, function(x)
              if (is.na(x['CHEMICAL_ID']))
                x['BIOCHEMICAL']
              else
                x['CHEMICAL_ID'])
          
          data <- data %>% select(-BIOCHEMICAL)
          
          ## CREATE DUMMY metOntology
          # if a property value is different to the one already existent then save
          # the property name in "diff" and overwrite property value
         
          query <- paste0("
                            UNWIND $properties AS prop with prop 
                            where exists(prop.sid)
                            MERGE (m:metOntology {sid:prop.sid}) 
                            ON CREATE SET m = prop 
                            ")
          runQuery(query,param = TRUE, param_list = unname(split(data, 1:nrow(data))) %>% lapply(function(x){x %>% select_if(~ !any(is.na(.)))}))
          
          
          ## MAP metMetabolon TO metOntology
          query <-  paste0("
                             UNWIND $properties AS prop 
                             MATCH (m:metMetabolon:metMeasured {sid:prop.sid}) 
                             MATCH (h:metOntology {sid:COALESCE(prop.CHEMICAL_ID,prop.BIOCHEMICAL)}) 
                             MERGE (m)-[:MAP]->(h)")
          runQuery(query,param = TRUE, param_list = metabo_neo)
           
          
          ## MAP metOntology TO HMDB
          query <- paste0("
                            UNWIND $properties AS prop with prop 
                            where exists(prop.CHEMICAL_ID) and exists(prop.HMDb_ID) 
                            MATCH (m:metOntology {sid:COALESCE(prop.CHEMICAL_ID,prop.BIOCHEMICAL)}) 
                            MATCH (h1:accHMDB)-[]-(h:metHMDB) where h1.sid=prop.HMDb_ID 
                            MERGE (m)-[:MAP]->(h)")
          runQuery(query,param = TRUE, param_list = metabo_neo)
          
        }else{ 
          msg=paste(Sys.time(),"metMetabolon","FORMAT",paste0("Chemical id coloumn (CHEMICAL_ID) missing: for ",source["file_name"]," no links to metOntology and metHMDB generated"),sep = "\t")
          cat(msg,file = err,sep="\n",append = T)
        }
        info(log,capture.output(toc()))
      }
    }
  }
  
  
  # 3 - Add info from metafile
  if(length(getFiles("metMetabolon","_meta"))!=0){
    
    query <- "
    LOAD CSV WITH HEADERS FROM 'file:///metMetabolon/metMetabolon_meta.csv' AS line with line 
    where exists(line.chemical_id) and exists(line.metabolonID)
    MERGE (n:metMetabolon:metMeasured  {COMP_ID:line.metabolonID})
    ON CREATE SET n.sid = line.metabolonID + '_unknown'
    MERGE (m:metOntology {sid:line.chemical_id,CHEMICAL_ID:line.chemical_id})
    ON MATCH SET m.HMDb_ID = apoc.coll.toSet(m.HMDb_ID + split(line.metHMDB,'|'))
    ON CREATE SET m.HMDb_ID = apoc.coll.toSet(split(line.metHMDB,'|'))
    MERGE (n)-[:MAP]->(m) 
    "
    runQuery(query)
    
    ## biochemical names (lower case - no special characters)
    query <- "
    LOAD CSV WITH HEADERS FROM 'file:///metMetabolon/metMetabolon_meta.csv' AS line with line 
    where exists(line.metabolonName) 
    MERGE (n:metMetabolon:metMeasured  {COMP_ID:line.metabolonID})
    ON CREATE SET n.sid = line.metabolonID + '_unknown'
    with n, split(line.metabolonName,'|') as names unwind names as name
    MERGE (m:metName {sid:apoc.text.replace(toLower(name), '[^A-Za-z0-9]', '')})
    MERGE (n)-[:HAS_NAME]->(m)
    "
    runQuery(query)
    
    # set name property of metMetabolon node
    query <- "
    LOAD CSV WITH HEADERS FROM 'file:///metMetabolon/metMetabolon_meta.csv' AS line with line 
    where exists(line.metabolonName) 
    MERGE (n:metMetabolon:metMeasured  {COMP_ID:line.metabolonID})
    with n, split(line.metabolonName,'|') as names 
    where not exists(n.name)
    SET n.name=names[0]
    "
    runQuery(query)
    
    ## biochemical names 
    query <- "
    LOAD CSV WITH HEADERS FROM 'file:///metMetabolon/metMetabolon_meta.csv' AS line with line 
    where exists(line.metabolonName) 
    MERGE (n:metMetabolon:metMeasured  {COMP_ID:line.metabolonID})
    with n, split(line.metabolonName,'|') as names unwind names as name
    MERGE (m:metName {sid:name})
    MERGE (n)-[:HAS_NAME]->(m)
    "
    runQuery(query)
    
    ##superpathway
    query <- "
    LOAD CSV WITH HEADERS FROM 'file:///metMetabolon/metMetabolon_meta.csv' AS line with line 
    where exists(line.super_pathway) 
    MERGE (n:metMetabolon:metMeasured  {COMP_ID:line.metabolonID})
    with n, split(line.super_pathway, '|') as ptwys unwind ptwys as pt
    MERGE (m:metSuperPathway {sid:pt})
    MERGE (n)-[:IN_PATHWAY]->(m) 
    "
    runQuery(query)
    
    ## subpathway
    query <- "
    LOAD CSV WITH HEADERS FROM 'file:///metMetabolon/metMetabolon_meta.csv' AS line with line 
    where exists(line.sub_pathway)
    MERGE (n:metMetabolon:metMeasured  {COMP_ID:line.metabolonID})
    with n, split(line.sub_pathway, '|') as ptwys unwind ptwys as pt
    MERGE (m:metSubPathway {sid:pt})
    MERGE (n)-[:IN_PATHWAY]->(m) 
    "
    runQuery(query)
    
    ##BiGG
    query <- "
    LOAD CSV WITH HEADERS FROM 'file:///metMetabolon/metMetabolon_meta.csv' AS line with line 
    where exists(line.metBiGG)
    MERGE (n:metMetabolon:metMeasured  {COMP_ID:line.metabolonID})
    with n, split(line.metBiGG, '|') as ids unwind ids as id
    MERGE (m:metBiGG {mid:id})
    MERGE (n)-[:HAS_ID]->(m) 
    "
    runQuery(query)
    
    #KEGG
    query <- "
    LOAD CSV WITH HEADERS FROM 'file:///metMetabolon/metMetabolon_meta.csv' AS line with line 
    where exists(line.metKEGG)
    MERGE (n:metMetabolon:metMeasured  {COMP_ID:line.metabolonID})
    with n, split(line.metKEGG, '|') as ids unwind ids as id 
    MERGE (m:metKEGG {sid:id})
    MERGE (n)-[:HAS_ID]->(m) 
    "
    runQuery(query)
    
    #PUBCHEM
    query <- "
    LOAD CSV WITH HEADERS FROM 'file:///metMetabolon/metMetabolon_meta.csv' AS line with line 
    where exists(line.metPUBCHEM)
    MERGE (n:metMetabolon:metMeasured  {COMP_ID:line.metabolonID})
    with n, split(line.metPUBCHEM, '|') as ids unwind ids as id 
    MERGE (m:metPubChem {sid:id})
    MERGE (n)-[:HAS_ID]->(m) 
    "
    runQuery(query)
    
    #HMDB
    query <- "
    LOAD CSV WITH HEADERS FROM 'file:///metMetabolon/metMetabolon_meta.csv' AS line with line 
    where exists(line.metHMDB)
    MERGE (n:metMetabolon:metMeasured  {COMP_ID:line.metabolonID})
    with n, split(line.metHMDB, '|') as ids unwind ids as id
    MATCH (h1:accHMDB)-[]-(h:metHMDB) where h1.sid=id 
    MERGE (n)-[:HAS_ID]->(h) 
    "
    runQuery(query)
  }
}
