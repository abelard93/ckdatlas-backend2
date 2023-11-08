MATCH (a:metMeta) where (a)-[:GENETIC_ASSOCIATION]-() set a :atlas
MATCH (a:proteinCoding) where (a)-[:GENETIC_ASSOCIATION]-() set a :atlas

MATCH (a:metMeta) where (a)-[:METABOLIC_ASSOCIATION]-() set a :atla

MATCH (a:metMeta) where (a)-[:PARTIAL_CORRELATION]-() set a :atlas

MATCH (a:proteinCoding) where (a)-[:COEXPRESSION]-() set a :atlas

MATCH (a:proteinCoding) where (a)-[:COREGULATION]-() set a :atlas


CALL algo.louvain.stream("atlas", "MATCH (:atlas)-[r]-(:atlas) where type(r) IN ['GENETIC_ASSOCIATION','METABOLIC_ASSOCIATION','PARTIAL_CORRELATION','COREGULATION','COEXPRESSION'] ", {})
YIELD nodeId, community with distinct community, collect(distinct nodeId) as n return community, size(n) as size

CALL algo.louvain.stream("MATCH (a:atlas) return ID(a) as id", "MATCH (:atlas)-[r]-(:atlas) where type(r) IN ['GENETIC_ASSOCIATION','METABOLIC_ASSOCIATION','PARTIAL_CORRELATION','COREGULATION','COEXPRESSION'] RETURN ID(STARTNODE(r)) as source, ID(ENDNODE(r)) as target", {graph:"cypher"})
YIELD nodeId, community with distinct community, collect(distinct nodeId) as n return community, size(n) as size