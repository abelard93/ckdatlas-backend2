//for (100) measured metabolites get genes that catalyse reactions they participate in (for Human)
MATCH (n:metMetabolon) with n.COMP_ID as metabolite LIMIT 100
MATCH (g:gene)-[:CATALYZES]-(:reBiGG)-[:PARTICIPATES]-(n:metBiGG)-[:MAP*]-(m:metMetabolon)
WHERE (n)-[:FROM]-(:source {sid:'Recon2'}) and m.COMP_ID = metabolite
return m.sid as Metabolite, collect(distinct g.Symbol) as Genes
