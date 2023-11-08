//for all WCM metabolites return hmdb id and annotated SMPDB pathways
match(m:metWCM)-[:MAP*2]-(h:metHMDB)-[:ANNOTATE]-(s:SMPDB)
return m.sid as metWCM, h.sid as hmdb ,collect(s.sid) as SMPDB
