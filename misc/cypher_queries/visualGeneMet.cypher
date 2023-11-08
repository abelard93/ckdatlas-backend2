//visualization of measured 'aspartate' metabolites to 'SLC1A6' gene pathway
match p=(n:metMetabolon {BIOCHEMICAL:"aspartate"})-[:MAP*3]-(metBiGG)-[:PARTICIPATES]-(:reBiGG)-[:CATALYZES]-(:gene {Symbol:"SLC1A6"}) return p
