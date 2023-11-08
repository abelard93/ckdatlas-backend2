USING PERIODIC COMMIT 5000 
LOAD CSV FROM 'file:///ncbi_ensembl_mappings.csv' AS line 
MATCH (a:Gene), (b:Gene) 
WHERE a.sid = line[0] AND b.sid = line[1] 
CREATE (a)-[r:MAPS]->(b) 
SET r.source = line[2] 
RETURN count(DISTINCT r)