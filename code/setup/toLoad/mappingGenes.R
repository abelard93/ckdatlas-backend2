#' Add gene - transcript - protein links

mappingGenes <- function(){
    
    tic("Add gene 2 transcript 2 protein info")
    
    # add gene - CODES -> transcript [from ensembl]
    query <- "
    USING PERIODIC COMMIT 1000 
    LOAD CSV WITH HEADERS FROM 'file:///gene/mapping/ensembl_gene2transcript2protein.csv' AS line 
    MERGE (n:gene:ensembl {sid:line.ensembl_gene_id, ensembl_id:line.ensembl_gene_id}) 
    ON CREATE SET n.source = 'mapping'
    with n, line 
    where exists(line.ensembl_transcript_id)
    MERGE (t:transcript:ensembl {sid:line.ensembl_transcript_id, ensembl_id:line.ensembl_transcript_id})
    ON CREATE SET t.source = 'mapping' with n, line, t
    MERGE (n)-[c:CODES]->(t)
    ON CREATE SET c.source='ensembl' 
    "
    runQuery(query,periodic = TRUE)
    
    # add transcript - CODES -> protein [from ensembl]
    query <- "
    USING PERIODIC COMMIT 1000 
    LOAD CSV WITH HEADERS FROM 'file:///gene/mapping/ensembl_gene2transcript2protein.csv' AS line 
    MERGE (t:transcript:ensembl {sid:line.ensembl_transcript_id, ensembl_id:line.ensembl_transcript_id})
    ON CREATE SET t.source = 'mapping' with line, t
    where exists(line.ensembl_peptide_id)
    MERGE (p:protein:ensembl {sid:line.ensembl_peptide_id, ensembl_id:line.ensembl_peptide_id})
    ON CREATE SET p.source = 'mapping' with p, line, t
    MERGE (t)-[s:CODES]->(p)
    ON CREATE SET s.source='ensembl'
    "
    runQuery(query,periodic = TRUE)
    
    # add gene - CODES -> transcript [from SNiPA (legacy)]
    query <- "
    USING PERIODIC COMMIT 1000 
    LOAD CSV WITH HEADERS FROM 'file:///gene/mapping/legacy/GeneID2TranscriptID2ProteinID.csv' AS line 
    MERGE (n:gene:ensembl {sid:line.Ensembl_Gene_ID, ensembl_id:line.Ensembl_Gene_ID}) 
    ON CREATE SET n.source = 'legacy'
    with n, line 
    where exists(line.Ensembl_Transcript_ID)
    MERGE (t:transcript:ensembl {sid:line.Ensembl_Transcript_ID, ensembl_id:line.Ensembl_Transcript_ID})
    ON CREATE SET t.source = 'legacy' with n, line, t
    MERGE (n)-[c:CODES]->(t)
    ON CREATE SET c.source='legacy' 
    "
    runQuery(query,periodic = TRUE)
    
    # add transcript - CODES -> protein [from SNiPA (legacy)]
    query <- "
    USING PERIODIC COMMIT 1000 
    LOAD CSV WITH HEADERS FROM 'file:///gene/mapping/legacy/GeneID2TranscriptID2ProteinID.csv' AS line 
    MERGE (t:transcript:ensembl {sid:line.Ensembl_Transcript_ID, ensembl_id:line.Ensembl_Transcript_ID})
    ON CREATE SET t.source = 'legacy' with line, t
    where exists(line.Ensembl_Protein_ID)
    MERGE (p:protein:ensembl {sid:line.Ensembl_Protein_ID, ensembl_id:line.Ensembl_Protein_ID})
    ON CREATE SET p.source = 'legacy' with p, line, t
    MERGE (t)-[s:CODES]->(p)
    ON CREATE SET s.source='legacy'
    "
    runQuery(query,periodic = TRUE)

    info(log,capture.output(toc()))
}
