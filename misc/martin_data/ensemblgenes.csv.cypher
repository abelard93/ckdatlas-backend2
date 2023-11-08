USING PERIODIC COMMIT 5000 
LOAD CSV FROM 'file:///ensemblgenes.csv' AS line 
CREATE (n:Gene) 
SET n.sid = line[0], n.name = line[1], n.taxid = line[2], n.source = line[3]