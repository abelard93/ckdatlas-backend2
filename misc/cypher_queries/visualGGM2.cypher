//retrieve (visual) GGM edges between metabolites from the super_pathway 'Amino Acid' and 'Lipid'
MATCH (m:metMetabolon) where m.SUPER_PATHWAY="Amino acid"
MATCH (n:metMetabolon) where n.SUPER_PATHWAY="Lipid"
MATCH p=(m)-[]-(e:edgeGGM)-[]-(n)
return p
