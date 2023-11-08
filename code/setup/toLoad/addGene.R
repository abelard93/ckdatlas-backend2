#' Add gene nodes from Ensembl.
#' 
#' Adds gene nodes defined by ensembl ids and adds/creates gene symbols from both ensembl and SNiPA.

addGene <- function(){

    tic("Added gene info")
    
    # add ensembl genes
    query <- "
    USING PERIODIC COMMIT 2000 
    LOAD CSV WITH HEADERS FROM 'file:///gene/ensembl_gene.csv' AS line 
    MERGE (g:gene:ensembl {sid:line.ensembl_gene_id, ensembl_id:line.ensembl_gene_id}) 
    SET g.description = line.description, g.source = 'ensembl', g.band = line.band, g.chromosome = line.chromosome_name, g.start_position = line.start_position, g.end_position = line.end_position, g.gene_type = line.gene_biotype, g.symbol = toUPPER(line.hgnc_symbol), g.version = 'current'
    with g where exists(g.symbol)
    MERGE (g)-[:SYNONYM]->(n:geneSymbol {sid:g.symbol, symbol_id:g.symbol})
    ON CREATE SET n.symbol_source='ensembl'
    "
    runQuery(query,periodic = TRUE)
    
    # add gene nodes from SNiPA (legacy) and add gene names/symbol
    query <- "
    USING PERIODIC COMMIT 2000 
    LOAD CSV WITH HEADERS FROM 'file:///gene/legacy/Gene2Symbol.csv' AS line 
    MERGE (g:gene:ensembl {sid:line.Ensembl_Gene_ID, ensembl_id:line.Ensembl_Gene_ID}) 
    ON CREATE SET g.symbol = line.HGNC_symbol, g.name_other = line.Associated_Gene_Name, g.chromosome = line.Chromosome_Name, g.start_position = line.Gene_Start_bp, g.end_position = line.Gene_End_bp, g.strand = line.Strand, g.version = 'legacy'
    with g where exists(g.symbol) 
    MERGE (g)-[:SYNONYM]->(n:geneSymbol {sid:g.symbol, symbol_id:g.symbol})
    ON CREATE SET n.symbol_source='legacy'
    "
    runQuery(query,periodic = TRUE)
    
    # add gene names/symbol specified in column 'name_other'
    query <- "
    USING PERIODIC COMMIT 2000 
    LOAD CSV WITH HEADERS FROM 'file:///gene/legacy/Gene2Symbol.csv' AS line 
    MATCH (g:gene:ensembl {sid:line.Ensembl_Gene_ID, ensembl_id:line.Ensembl_Gene_ID}) 
    with g where exists(g.name_other) 
    MERGE (g)-[:SYNONYM]->(n:geneSymbol {sid:g.name_other, symbol_id:g.name_other})
    ON CREATE SET n.symbol_source='legacy'
    "
    runQuery(query,periodic = TRUE)
    
    info(log,capture.output(toc()))
}


#' Add gene-wise cutoff to gene nodes
#' 
#' Adds the gene-wise cutoff, derived from the SNP to gene annotation from SNiPA and defined as 0.05/#SNPs, to gene nodes.

addSNiPACutoff <- function(){
  
  tic("Added SNiPA cutoff info")

  query <- "
  USING PERIODIC COMMIT 2000 
  LOAD CSV WITH HEADERS FROM 'file:///gene/numsnps_per_gene_cutoff.csv' AS line 
  MERGE (g:gene:ensembl {sid:line.sid}) 
  SET g.NUMSNPS_SNiPA=line.NUMSNPS_SNiPA, g.cutoff=toFloat(line.cutoff), g.cutoff_info=line.cutoff_info
  "
  runQuery(query,periodic = TRUE)
  
  info(log,capture.output(toc()))
}
