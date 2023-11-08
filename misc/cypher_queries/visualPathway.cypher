//visualization of measured ALL ’adenosine’ metabolite to gene pathways
match p=(n:metMetabolon {BIOCHEMICAL:"adenosine"})-[:MAP*3]-(metBiGG)-[:PARTICIPATES]-(:reBiGG)-[:CATALYZES]-(:gene) return p
