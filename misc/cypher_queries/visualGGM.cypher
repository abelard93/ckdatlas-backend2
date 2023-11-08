//visualization of two weakly correlated (GGM) metabolites and their connectivity to gene ‘NDUFAB1’
match p=(n:metMetabolon {sid:"59_lcmsneg"})-[:MAP*3]-(:metBiGG)-[:PARTICIPATES]-(:reBiGG)-[:CATALYZES]-(:gene {Symbol:"NDUFAB1"})-[:CATALYZES]-(:reBiGG)-[:PARTICIPATES]-(:metBiGG)-[:MAP*3]-(m:metMetabolon {sid:"1302_lcmspos"})-[:GGM*2]-(n) return p
