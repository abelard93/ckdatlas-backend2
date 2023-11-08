#################
## genes
#################

addGenes <- function(){
    # create index on Symbol for faster search in functions
    query <- "CREATE INDEX ON :gene(Symbol)"
    cypher(db,query)
    
    tic("Added gene info")
    
    ## ADD AND MAP GENE AND GENE SYMBOL NODES
    # ncbi genes
    query <- "
    USING PERIODIC COMMIT 5000 
    LOAD CSV FROM 'file:///genes/ncbigenes.csv' AS line 
    CREATE (n:gene) 
    SET n.sid = line[0], n.taxID = line[1], n.GeneID = line[2], n.Symbol = line[3], n.LocusTag = line[4], n.Synonyms = line[5], n.dbXrefs = line[6], n.chromosome = line[7], n.map_location = line[8], n.description = line[9], n.type_of_gene = line[10], n.Symbol_from_nomenclature_authority = line[11], n.Full_name_from_nomenclature_authority = line[12], n.Nomenclature_status = line[13], n.Other_designations = line[14], n.Modification_date = line[15], n.Feature_type = line[16],n.source='ncbi'
    with n.Symbol as symbol,n
    MERGE (g:geneSymbol {sid:toUPPER(symbol)}) 
    MERGE (g)-[:SYNONYM]->(n)"
    cypher(db,query)
    # add synonyms as gene symbol and map to gene node
    query <- "
    USING PERIODIC COMMIT 5000 
    LOAD CSV FROM 'file:///genes/ncbigenes.csv' AS line 
    MATCH (n:gene {sid:line[0]})
    with split(n.Synonyms,'|') as names,n
    UNWIND names as symbol
    MERGE (g:geneSymbol {sid:toUPPER(symbol)}) 
    MERGE (g)-[:SYNONYM]->(n)"
    cypher(db,query)
    # ensemble genes
    query <- "
    USING PERIODIC COMMIT 5000 
    LOAD CSV FROM 'file:///genes/ensemblgenes.csv' AS line 
    CREATE (n:gene) 
    SET n.sid = line[0], n.Symbol = line[1], n.taxID = line[2], n.source = line[3]
    with n.Symbol as symbol,n
    MERGE (g:geneSymbol {sid:toUPPER(symbol)}) 
    MERGE (g)-[:SYNONYM]->(n)"
    cypher(db,query)
    
    ## MAP NCBI TO ENSEMBLE
    query <- "
    USING PERIODIC COMMIT 5000 
    LOAD CSV FROM 'file:///genes/ncbi_ensembl_mappings.csv' AS line 
    MATCH (a:gene), (b:gene) 
    WHERE a.sid = line[0] AND b.sid = line[1] 
    CREATE (a)-[r:MAP]->(b) 
    SET r.source = line[2]"
    cypher(db,query)
    
    info(log,capture.output(toc()))
}
