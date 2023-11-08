USING PERIODIC COMMIT 5000 
LOAD CSV FROM 'file:///ncbigenes.csv' AS line 
CREATE (n:Gene) 
SET n.sid = line[0], n.tax_id = line[1], n.GeneID = line[2], n.Symbol = line[3], n.LocusTag = line[4], n.Synonyms = line[5], n.dbXrefs = line[6], n.chromosome = line[7], n.map_location = line[8], n.description = line[9], n.type_of_gene = line[10], n.Symbol_from_nomenclature_authority = line[11], n.Full_name_from_nomenclature_authority = line[12], n.Nomenclature_status = line[13], n.Other_designations = line[14], n.Modification_date = line[15], n.Feature_type = line[16]