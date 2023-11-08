//get measured metabolon metabolites that participate in 5-10 reactions (in human) and belong to the super pathway 'Amino Acid'
MATCH (:reBiGG)-[p:PARTICIPATES {source:"Recon2"}]-(m:metBiGG) WITH count(p) AS nr, m AS metabo
WHERE 5<nr<10
MATCH p=(:reBiGG)-[:PARTICIPATES {source:"Recon2"}]-(metabo)-[:MAP*3]-(m:metMetabolon {SUPER_PATHWAY:'Amino Acid'}) RETURN p
