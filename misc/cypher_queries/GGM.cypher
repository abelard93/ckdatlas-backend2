//retrieve GGM edges between metabolites from the super_pathway 'Amino Acid' and 'Lipid'
MATCH (m:metMetabolon) where m.SUPER_PATHWAY="Amino acid"
MATCH (n:metMetabolon) where n.SUPER_PATHWAY="Lipid"
MATCH (m)-[]-(e:edgeGGM)-[]-(n) with m,n,e
return m.sid as sid_a, m.BIOCHEMICAL as name_a, n.sid as sid_b, n.BIOCHEMICAL as name_b, e.partialCorr as partialCorr, e.pvalue as pvaue
