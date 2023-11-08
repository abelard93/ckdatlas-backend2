// get all measured metabolites (metabolon) in liver samples, return biochemical and hmdb
MATCH (m:metMetabolon)-[:FROM]-(s:source)
WHERE s.sample_type='liver'
RETURN distinct(m.sid) as sid, m.BIOCHEMICAL as name, m.HMDb_ID as hmdb
