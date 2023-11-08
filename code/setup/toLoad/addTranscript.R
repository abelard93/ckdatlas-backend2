#' Add transcripts defined by ensembl

addTranscript <- function(){

    tic("Added transcript info")
    
    # add ensembl transcripts
    query <- "
    USING PERIODIC COMMIT 2000 
    LOAD CSV WITH HEADERS FROM 'file:///transcript/ensembl_transcript.csv' AS line 
    MERGE (g:transcript:ensembl {sid:line.ensembl_transcript_id, ensembl_id:line.ensembl_transcript_id}) 
    SET g.band = line.band, g.chromosome = line.chromosome_name, g.transcript_start = line.transcript_start, g.transcript_end = line.transcript_end, g.transcript_type = line.transcript_biotype, g.transcription_start_site = line.transcription_start_site, g.version = 'current'
    "
    runQuery(query,periodic = TRUE)
    
    
    # add transcripts from SNiPA (legacy)
    query <- "
    USING PERIODIC COMMIT 2000 
    LOAD CSV WITH HEADERS FROM 'file:///transcript/legacy/Transcript2Position.csv' AS line 
    MERGE (g:transcript:ensembl {sid:line.Ensembl_Transcript_ID, ensembl_id:line.Ensembl_Transcript_ID}) 
    ON CREATE SET g.strand = line.Strand, g.chromosome = line.Chromosome_Name, g.transcript_start = line.Transcript_Start_bp, g.transcript_end = line.Transcript_End_bp, g.transcript_type = line.transcript_biotype, g.transcription_start_site = line.transcription_start_site, g.version = 'legacy'
    "
    runQuery(query,periodic = TRUE)
    
    info(log,capture.output(toc()))
}
