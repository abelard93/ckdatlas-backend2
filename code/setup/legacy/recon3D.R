#################
## Recon3D
#################

addRecon3D <- function(){

  file <- "recon3d/Recon3D_301.xml"
  # read info on dataset
  source <- getInfo(file,c("sid","tax","taxID","file_name"),err)
  # skip file if no/incomplete info
  if(!all(is.na(source))){
    tic(paste0("Added Recon3D ",source["sid"]))
    
    # check if node for source is present, if not create
    checkSource(source,db,log)
    ## ADD RECON GENES: ENTREZ TO SYMBOL MAPPING
    query <- paste0("
                        USING PERIODIC COMMIT 1000 
                        LOAD CSV WITH HEADERS FROM 'file:///recon3d/recon3d_genes.csv' AS props 
                        MERGE (s:geneSymbol {sid:props.symbol})
                        MERGE (g:gene {sid:props.entrez_id})
                        set g.recon3d_id=props.recon3d_id
                        MERGE (s)-[:SYNONYM]-> (g)
                        ")
    cypher(db,query)
    ## ADD RECON METABOLITES
    query <- paste0("
                        USING PERIODIC COMMIT 1000
                        LOAD CSV WITH HEADERS FROM 'file:///recon3d/recon3d_metabolites.csv' AS props
                        MERGE (n:metBiGG:recon3D {sid:props.sid})
                        on create SET n += props with n
                        MATCH (s:source {sid:\"",source["sid"],"\"})
                        MERGE (n)-[:FROM]->(s)
                        ")
    cypher(db,query)
    ## MAP METABOLITES TO HMDB
    query <- paste0("
                        USING PERIODIC COMMIT 1000
                        LOAD CSV WITH HEADERS FROM 'file:///recon3d/recon3d_metabolites.csv' AS props
                        MATCH (n:metBiGG:recon3D {sid:props.sid}) WHERE exists(n.hmdb) with n, split(n.hmdb,'|') as id
                        UNWIND id as hmdb
                        MATCH (h1:accHMDB)-[]-(h:metHMDB) where h1.sid=hmdb 
                        with n,h
                        MERGE (n)-[:MAP]->(h)
                        ")
    cypher(db,query)
    ## ADD RECON REACTIONS
    query <- paste0("
                        USING PERIODIC COMMIT 1000
                        LOAD CSV WITH HEADERS FROM 'file:///recon3d/recon3d_reactions.csv' AS props
                        MERGE (n:reBiGG:recon3D {sid:props.sid})
                        ON CREATE SET n += props with n
                        MATCH (s:source {sid:\"",source["sid"],"\"})
                        MERGE (n)-[:FROM]->(s)
                        ")
    cypher(db,query)
    ## ADD "CATALYZES" REACTION
    query <- paste0("
                        USING PERIODIC COMMIT 1000
                        LOAD CSV WITH HEADERS FROM 'file:///recon3d/recon3d_re_gene.csv' AS props 
                        MATCH (g:gene) where g.sid=toUpper(props.gene)
                        MATCH (r:reBiGG:recon3D {sid:props.reBiGG})
                        MERGE (g)-[:CATALYZES {source:\"",source["sid"],"\",taxID:\"",source["taxID"],"\"}]->(r) 
                        ")
    cypher(db,query)
    ## ADD DIRECTED "PARTICIPATES" RELATIONSHIP
    # positiv stochiometry: reaction --> metabolite 
    query <- paste0("
                        USING PERIODIC COMMIT 1000 
                        LOAD CSV WITH HEADERS FROM 'file:///recon3d/recon3d_met_re.csv' AS props 
                        with props where toFloat(props.stoichiometry) >= 0
                        MATCH (m:metBiGG:recon3D {sid:props.metBiGG}) 
                        MATCH (r:reBiGG:recon3D {sid:props.reBiGG})
                        MERGE (r)-[:PARTICIPATES {source:\"",source["sid"],"\",taxID:\"",source["taxID"],"\", stoichiometry:props.stoichiometry}]->(m)
                        ")
    cypher(db,query)
    # negative stochiometry: metabolite --> reaction
    query <- paste0("
                        USING PERIODIC COMMIT 1000 
                        LOAD CSV WITH HEADERS FROM 'file:///recon3d/recon3d_met_re.csv' AS props 
                        with props where toFloat(props.stoichiometry) < 0
                        MATCH (m:metBiGG:recon3D {sid:props.metBiGG}) 
                        MATCH (r:reBiGG:recon3D {sid:props.reBiGG})
                        MERGE (m)-[:PARTICIPATES {source:\"",source["sid"],"\",taxID:\"",source["taxID"],"\", stoichiometry:props.stoichiometry}]->(r)
                        ")
    cypher(db,query)
    
    query <- "CREATE INDEX ON :metBiGG(mid)"
    cypher(db,query)
    
    query <- "CREATE INDEX ON :metBiGG(id)"
    cypher(db,query)
    info(log,capture.output(toc()))
  }
}