# uniformize metabolite ids according to list 
# list provided by script written by Parviz


query <- paste0("
                LOAD CSV WITH HEADERS FROM 'file:////fixes/uniformize_ids.csv' AS fix
                MATCH (n:metBiGG:recon3D {mid:fix.mid})
                SET n.id=fix.replacement+'_'+n.`compartment.id`"
)
cypher(db,query)

query <- paste0("
                LOAD CSV WITH HEADERS FROM 'file:////fixes/uniformize_ids.csv' AS fix
                MATCH (n:metBiGG:recon3D {mid:fix.replacement})
                MATCH (m:metBiGG:recon3D {mid:fix.mid}) with m,n
                where n.id=m.id with collect(m)+[n] as nodes
                call apoc.refactor.mergeNodes(nodes,{mergeRels:true}) yield node
                return count(*)"
)
cypher(db,query)

query <- paste0("
                MATCH (n:metBiGG:recon3D)-[r:FROM]-(s:source {sid:'Recon3D'})
                WITH n, type(r) as type, collect(r) as rels, s
                WHERE size(rels) > 1
                UNWIND tail(rels) as rel
                DELETE rel"
)
cypher(db,query)
