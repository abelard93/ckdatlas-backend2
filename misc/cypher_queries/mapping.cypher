//get HMDB ids (primary and secondary) for all measured metabolon metabolites
MATCH (m:metMetabolon)-[:MAP*2]-(h:metHMDB)
return m.sid as sid, m.BIOCHEMICAL as name, h.sid as primary_HMDB, collect(h.sa) as secondary_HMDB
