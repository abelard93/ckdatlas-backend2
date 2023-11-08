###
## ADD SIGN. LINKS
###
# - all associations with p-/q-value below 
#   individual threshold
# - will help speedup and simplify analysis
# - for now only regard protein coding 
###

addSignRels <- function(){
  
  tic("Added significant relationships info")
  
  # BUILD TABIX FILE -----------------
  # for HTML label lookup 
  files <- c("CKDatlas/coabundance/ggm_protein.csv",
             #"CKDatlas/coexpression/coexpr.csv",
             "CKDatlas/genetic_association/association_gene_trait.csv",
             "CKDatlas/genetic_association/association_metabo_gene.csv",
             "CKDatlas/metabolic_association/association_metabo_trait.csv",
             "CKDatlas/partial_correlation/ggm.csv",
             "CKDatlas/coregulation/significant_eQTL.csv",
             "CKDatlas/coregulation/significant_eQTL_all.csv")
  count = 0
  
  # add ID to relationships 
  for(file in files){
    
    system(paste0("awk -F',' -v OFS=',' 'NR == 1 {print \"\\\"ID\\\"\", $0; next} {print (NR-1+",count,"), $0}' ",file," > ",str_replace(file,".csv","_indexed.csv")," && mv ",str_replace(file,".csv","_indexed.csv")," ",file))
    count = count + (read.table(pipe(paste0("wc -l ", file)))[[1]] - 1)
  }
  
  
  # build tabix
 # system("cd CKDatlas && bash buildTabix.sh")
  
  ## eQTL gene-gene links -------
  
  
  ## COMMENTED FOR NOW TO KEEP DIRECTION !! RELATIONSHIP NOT SYMMETRIC
  ## PROTEIN A -> SNP -> eQTL -> PROTEIN B
  
  query <- paste0("CALL apoc.periodic.iterate(\"
  LOAD CSV WITH HEADERS FROM 'file:///CKDatlas/coregulation/significant_eQTL_all.csv' AS line RETURN line\", 
  \"MATCH (s:proteinCoding {ensembl_id:line.gene_snp})
  MATCH (g:proteinCoding {ensembl_id:line.gene})
  CREATE (s)-[:COREGULATION {pvalue: toFloat(line.pvalue), index: toInteger(line.ID), source: split(line.source,'|')}]->(g)
  \",{batchSize:100, parallel:false})")
  
  runQuery(query, periodic = TRUE)
  
  
  
  
  query <-  paste0("CALL apoc.periodic.iterate(\"
      LOAD CSV WITH HEADERS FROM 'file:///CKDatlas/coregulation/significant_eQTL.csv' AS line with line return line\", 
      \"MATCH (s:proteinCoding {ensembl_id:line.gene_snp})
      MATCH (g:proteinCoding {ensembl_id:line.gene})
      CALL apoc.create.relationship(s, line.relationship, {pvalue: toFloat(line.pvalue), type: line.type, tissue: line.tissue, index: toInteger(line.ID), source: split(line.source,'|')}, g)
      YIELD rel
      RETURN count(*)
       \",{batchSize:500, parallel:false})")
  runQuery(query, periodic = TRUE)
  
  
  
  
  ## COEXPR gene-gene links -------
  
  ## For the moment this comment is commented out because there is no data on it. Uncomment it if you have data
  
 # query <-  paste0("CALL apoc.periodic.iterate(\"
  #LOAD CSV WITH HEADERS FROM 'file:///CKDatlas/coexpression/coexpr.csv' AS line RETURN line\", 
  #\"MATCH (g1:proteinCoding {ensembl_id:line.gene_a})
  #MATCH (g2:proteinCoding {ensembl_id:line.gene_b})
  #CREATE (g1)-[:COEXPRESSION {tissue: line.tissue, index: toInteger(line.ID)}]->(g2)
  #CREATE (g2)-[:COEXPRESSION {tissue: line.tissue, index: toInteger(line.ID)}]->(g1)\",{batchSize:500, parallel:false})")
  
  #runQuery(query, periodic = TRUE)
  
  
  
  
  ## GGM metabolite-metabolite links -------
  
  
  query <-  paste0("CALL apoc.periodic.iterate(\"
  LOAD CSV WITH HEADERS FROM 'file:///CKDatlas/partial_correlation/ggm.csv' AS line RETURN line\", 
  \"MATCH (m1:metMeta) where line.metabolite_a in m1.sid
  MATCH (m2:metMeta) where line.metabolite_b in m2.sid
  CREATE (m1)-[:PARTIAL_CORRELATION {pvalue: toFloat(line.pvalue), index: toInteger(line.ID), sample_type: line.sample_type}]->(m2)
  CREATE (m2)-[:PARTIAL_CORRELATION {pvalue: toFloat(line.pvalue), index: toInteger(line.ID), sample_type: line.sample_type}]->(m1)\",{batchSize:500, parallel:false})")
  
  runQuery(query, periodic = TRUE)
  
  
  
  ## GENETIC_ASSOCIATION gene-metabolite links -------
  
  
  query <-  paste0("CALL apoc.periodic.iterate(\"
  LOAD CSV WITH HEADERS FROM 'file:///CKDatlas/genetic_association/association_metabo_gene.csv' AS line RETURN line\", 
  \"MATCH (s:proteinCoding {ensembl_id:line.gene})
  MATCH (m:metMeta) where line.metabolite in m.sid
  CREATE (s)-[:GENETIC_ASSOCIATION {pvalue: toFloat(line.pvalue), index: toInteger(line.ID), sample_type: line.sample_type}]->(m)\",{batchSize:500, parallel:false})")
  
  runQuery(query, periodic = TRUE)
  
  
  
  
  
  ## METABOLIC_ASSOCIATION trait-metabolite links -------
  
  
  query <-  paste0("CALL apoc.periodic.iterate(\"
  LOAD CSV WITH HEADERS FROM 'file:///CKDatlas/metabolic_association/association_metabo_trait.csv' AS line RETURN line\", 
  \"MATCH (m:metMeta) where line.metabolite in m.sid
  MATCH (t:trait) where line.trait IN t.sid
  CREATE (t)<-[:METABOLIC_ASSOCIATION {pvalue: toFloat(line.pvalue), index: toInteger(line.ID), sample_type: line.sample_type}]-(m)\",{batchSize:500, parallel:false})")
  
  runQuery(query, periodic = TRUE)
  
  
  
  ## GENETIC_ASSOCIATION trait-gene links -------
  
  
  query <-  paste0("CALL apoc.periodic.iterate(\"
  LOAD CSV WITH HEADERS FROM 'file:///CKDatlas/genetic_association/association_gene_trait.csv' AS line RETURN line\", 
  \"MATCH (s:proteinCoding {ensembl_id:line.gene})
  MATCH (t:trait) where line.trait IN t.sid
  CREATE (t)<-[:GENETIC_ASSOCIATION {pvalue: toFloat(line.pvalue), index: toInteger(line.ID)}]-(s)\",{batchSize:500, parallel:false})")
  
  runQuery(query, periodic = TRUE)
  
  
  ## COABUNDANCE gene-gene (protein-protein) links -------
  
  query <-  paste0("CALL apoc.periodic.iterate(\"
  LOAD CSV WITH HEADERS FROM 'file:///CKDatlas/coabundance/ggm_protein.csv' AS line RETURN line\", 
  \"MATCH (a:proteinCoding {ensembl_id:line.gene_a})
  MATCH (b:proteinCoding {ensembl_id:line.gene_b})
  CREATE (a)-[:COABUNDANCE {pvalue: toFloat(line.pvalue), index: toInteger(line.ID), sample_type: line.sample_type}]->(b)
  CREATE (b)-[:COABUNDANCE {pvalue: toFloat(line.pvalue), index: toInteger(line.ID), sample_type: line.sample_type}]->(a)\",{batchSize:500, parallel:false})")
  
  runQuery(query, periodic = TRUE)
  
  info(log,capture.output(toc()))
}


# #needed for queries
# query <- 'CALL apoc.periodic.iterate(
#   "MATCH ()-[r:GENETIC_ASSOCIATION]-() RETURN r",
#   "SET r.pvalue = toFloat(r.pvalue)",
#   {batchSize:100, parallel:true})'
# 
# n <- cypher(db,query)


# 
# query <- paste0("CALL apoc.export.csv.query(\"MATCH (g:proteinCoding)-[:MAP]-(s:SNP)-[:HIT]-(q:mQTL) RETURN distinct g.ensembl_id as gene,s.sid as snp, g.cutoff, 'snipa' as tissue \",\"mgwas_snipa.csv\", {})")
# cypher(db,query)
# 
# query <- paste0("CALL apoc.export.csv.query(\"MATCH (g:proteinCoding)-[:HIT]-(e:eQTL)-[:HIT]-(s:SNP)-[:HIT]-(q:mQTL) MATCH (e)-[:FROM]-(su:source) return distinct g.ensembl_id as gene,s.sid as snp, g.cutoff, collect(distinct su.sample_type) as tissue \",\"mgwas_coreg.csv\", {})")
# cypher(db,query)