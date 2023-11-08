//get measured metabolon metabolites (human) that participate in multiple reactions (>1) and belong to the super pathway 'Amino Acid'
MATCH (:reBiGG)-[p:PARTICIPATES]-(m:metBiGG) where (m)-[:FROM]-(:source {sid:"Recon2"}) WITH count(p) AS nr, m AS metabo
WHERE nr>1
MATCH (metabo)-[:MAP*3]-(m:metMetabolon)
WHERE m.SUPER_PATHWAY='Amino Acid' RETURN distinct(m.sid) as sid, m.BIOCHEMICAL as name, m.HMDb_ID as hmdb
