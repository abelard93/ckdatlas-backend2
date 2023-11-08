//show 4 measured metabolites and the genes that catalyse reactions they participate in (for Human)
MATCH (n:metMetabolon) with n.COMP_ID as metabolite LIMIT 37
MATCH p=(g:gene)-[:CATALYZES]-(:reBiGG)-[:PARTICIPATES]-(n:metBiGG)-[:MAP*]-(m:metMetabolon)
WHERE (n)-[:FROM]-(:source {sid:'Recon2'}) and m.COMP_ID = metabolite
return p
