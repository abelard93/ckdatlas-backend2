#' Add proteins defined by ensembl and add mapping to uniprot

addProtein <- function(){
  
  tic("Added protein info")
  
  # add ensembl proteins
  query <- "
  USING PERIODIC COMMIT 1000 
  LOAD CSV WITH HEADERS FROM 'file:///protein/ensembl_protein.csv' AS line 
  MERGE (g:protein:ensembl {sid:line.ensembl_peptide_id, ensembl_id:line.ensembl_peptide_id}) 
  SET g.uniprot_swissprot_id = line.uniprotswissprot, g.uniprot_trembl_id = line.uniprotsptrembl, g.version = 'current'
  "
  runQuery(query,periodic = TRUE)

  
  # add proteiens from SNiPA (legacy)
  query <- "
  USING PERIODIC COMMIT 1000 
  LOAD CSV WITH HEADERS FROM 'file:///protein/legacy/ProteinID2UniProt.csv' AS line 
  MERGE (g:protein:ensembl {sid:line.Ensembl_Protein_ID, ensembl_id:line.Ensembl_Protein_ID}) 
  ON CREATE SET g.uniprot_swissprot_id = line.UniProt_SwissProt, g.uniprot_gene_name = line.Uniprot_Gene_Name, g.version = 'legacy'
  "
  runQuery(query,periodic = TRUE)
  
  # make uniprot node and link to ensembl where uniprot id is given (ensembl based)
  # swissprot
  query <- "
  MATCH (g:protein:ensembl) where exists(g.uniprot_swissprot_id) with g
  MERGE (u:uniprot {sid:g.uniprot_swissprot_id,uniprot_id:g.uniprot_swissprot_id})-[:MAP]->(g)
  ON CREATE SET u.source = 'ensembl' 
  "
  runQuery(query,periodic = TRUE)
  
  #trembl
  query <- "
  MATCH (g:protein:ensembl) where exists(g.uniprot_trembl_id) with g
  MERGE (u:uniprot {sid:g.uniprot_trembl_id, uniprot_id:g.uniprot_trembl_id})-[:MAP]->(g)
  ON CREATE SET u.source = 'ensembl'
  "
  runQuery(query,periodic = TRUE)
  
  ## MATCH!! on ensembl to ensure connectivity
  # make uniprot node and link to ensembl where uniprot id is given (uniprot based)
  query <- "
  USING PERIODIC COMMIT 1000 
  LOAD CSV WITH HEADERS FROM 'file:///protein/uniprot.csv' AS line with line where exists(line.ensembl_id)
  with split(line.ensembl_id,'|') as ens, line.uniprot_id as id
  UNWIND ens as ensemblID
  MATCH (g:protein:ensembl {sid:ensemblID, ensembl_id:ensemblID}) 
  MERGE (u:uniprot {sid:id, uniprot_id:id})-[:MAP]->(g)
  ON CREATE SET u.source='uniprot'
  "
  runQuery(query,periodic = TRUE)
  
  #"MATCH (n:protein) where apoc.meta.type(n.uniprot_id)='STRING' set n.uniprot_id=[n.uniprot_id]" %>% cypher(db,.)
  #"MATCH (n:protein) where apoc.meta.type(n.uniprot_swissprot_id)='STRING' set n.uniprot_swissprot_id=[n.uniprot_swissprot_id]" %>% cypher(db,.)
  
  info(log,capture.output(toc()))
}

