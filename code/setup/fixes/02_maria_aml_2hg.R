
########################################
# Add D/L for 
# Mitochondria, Cytosol, Extracellular
########################################

query <- paste0("
                MATCH (n:metBiGG:recon3D {id:'2hydog_c'})
                return n"
)
node <- as.data.frame(t(cbind(unlist(cypherToList(db,query)))))
node <- node %>% setNames(sub('..', '', names(node))) %>% 
  mutate(freq = 6) %>% 
  uncount(freq) %>% 
  mutate(sid = paste0(rep(c("D","L"),3),"2hydog_",c("c","c","m","m","e","e"),"_Recon3D"),
         id = paste0(rep(c("D","L"),3),"2hydog_",c("c","c","m","m","e","e")),
         mid = paste0(rep(c("D","L"),3),"2hydog"),
         compartment = rep(c("Cytoplasm","Mitochondrion","Extracellular"),c(2,2,2)),
         compartment.id = rep(c("c","m","e"),c(2,2,2)),
         enantiomere = rep(c("D","L"),3))

query <- paste0("UNWIND {node} as n
                MERGE (m:metBiGG:recon3D {sid:n.sid})
                SET m += n"
)
cypher(db,query,node=node)


# add connection to hmdb id
query <- paste0("
                MATCH (a:metBiGG:recon3D) where a.sid starts with 'D2hydog'
                MATCH (b:metHMDB {sid:'HMDB0000606'})
                MERGE (a)-[:MAP]->(b)"
)
cypher(db,query)

query <- paste0("
                MATCH (a:metBiGG:recon3D) where a.sid starts with 'L2hydog'
                MATCH (b:metHMDB {sid:'HMDB0000694'})
                MERGE (a)-[:MAP]->(b)"
)
cypher(db,query)

query <- paste0("
                MATCH (a:metBiGG:recon3D) where a.name='2-Hydroxyglutarate'
                MATCH (b:metHMDB {sid:'HMDB0059655'})
                MERGE (a)-[:MAP]->(b)"
)
cypher(db,query)

query <- paste0("
                MATCH (a:metBiGG:recon3D) where a.name='2-Hydroxyglutarate'
                MATCH (s:source {sid:'Recon3D'})
                MERGE (a)-[:FROM]->(s)"
)
cypher(db,query)

#################
# connect fad

query <- paste0("
                match(a:metBiGG:recon3D) where a.fullName='FAD'
                MATCH (b:metHMDB {sid:'HMDB0001248'})
                MERGE (a)-[:MAP]->(b)"
)
cypher(db,query)


#################################################################################
# rename 2HG_e to 3hydog

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'2hydog_e'})
                SET a.sid='3hydog_e_Recon3D', a.fullName='3-hydroxyglutaric acid', a.mid='3hydog', a.id='3hydog_e', a.inchi='1S/C5H8O5/c6-3(1-4(7)8)2-5(9)10/h3,6H,1-2H2,(H,7,8)(H,9,10)', a.avgmolweight='148.114',a.monoisotopicweight='148.037173366'"
                )

cypher(db,query)

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'2hydog_c'})
                SET a.sid='3hydog_c_Recon3D', a.fullName='3-hydroxyglutaric acid', a.mid='3hydog', a.id='3hydog_c', a.inchi='1S/C5H8O5/c6-3(1-4(7)8)2-5(9)10/h3,6H,1-2H2,(H,7,8)(H,9,10)', a.avgmolweight='148.114',a.monoisotopicweight='148.037173366'"
)

cypher(db,query)



# delete 2HG_e
query <- paste0("
                MATCH (m:metBiGG:recon3D {id:'2hydog_e'})
                detach delete m")

cypher(db,query)

#################################################################################
# rename M00653 to 2hydog

query <- paste0("
                MATCH (a:metBiGG:recon3D {mid:'M00653'})
                SET a.sid='2hydog_c_Recon3D', a.mid='2hydog', a.id='2hydog_c'"
)

cypher(db,query)


####################
# Mitochondrion
####################


# 1
#####
# delete relationship to reaction that takes place in mitochondrion
query <- paste0("
                MATCH (:metBiGG:recon3D {id:'akg_c'})-[r2]-(:reBiGG:recon3D {id:'HMR_0718'})-[r1]-(:metBiGG:recon3D {id:'2hydog_c'})
                MATCH (h:metBiGG:recon3D {id:'L2hydog_m'}), (r:reBiGG:recon3D {id:'HMR_0718'}),(a:metBiGG:recon3D {id:'akg_m'})
                MERGE (h)-[:PARTICIPATES {source:r1.source, stoichiometry:r1.stoichiometry,taxID:r1.taxID}]->(r)-[:PARTICIPATES {source:r2.source, stoichiometry:r2.stoichiometry,taxID:r2.taxID}]->(a)
                delete r1,r2"
)
cypher(db,query)

# 2
#####
# add rections to D2HG_c from existing 2HG_c
query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'2hydog_c'})-[r]->(l) where 'reBiGG' in labels(l)
                with collect([type(r),l.sid,properties(r)]) as rels
                unwind rels as r
                MATCH (m:metBiGG:recon3D {compartment:'Cytoplasm'}) where m.enantiomere='D'
                MATCH (b:reBiGG:recon3D {sid:r[1]}) 
                CALL apoc.create.relationship(m, r[0], {source:r[2].source, stoichiometry:r[2].stoichiometry,taxID:r[2].taxID}, b) yield rel return count(*)"
)
cypher(db,query)

query <- paste0("
                MATCH (l)-[r]->(:metBiGG:recon3D {id:'2hydog_c'}) where 'reBiGG' in labels(l)
                with collect([type(r),l.sid,properties(r)]) as rels
                unwind rels as r
                MATCH (m:metBiGG:recon3D {compartment:'Cytoplasm'}) where m.enantiomere='D'
                MATCH (b:reBiGG:recon3D {sid:r[1]}) 
                CALL apoc.create.relationship(b, r[0], {source:r[2].source, stoichiometry:r[2].stoichiometry,taxID:r[2].taxID}, m) yield rel return count(*)"
)
cypher(db,query)

# delete D specific relationship to reaction 
query <- paste0("
                MATCH (:reBiGG:recon3D {id:'HMR_0719'})-[r1]-(:metBiGG:recon3D {id:'2hydog_c'})
                delete r1"
)
cypher(db,query)

# add non enantiomere specific rections to L2HG_c from existing 2HG_c
query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'2hydog_c'})-[r]->(l) where 'reBiGG' in labels(l)
                with collect([type(r),l.sid,properties(r)]) as rels
                unwind rels as r
                MATCH (m:metBiGG:recon3D {compartment:'Cytoplasm'}) where m.enantiomere='L'
                MATCH (b:reBiGG:recon3D {sid:r[1]}) 
                CALL apoc.create.relationship(m, r[0], {source:r[2].source, stoichiometry:r[2].stoichiometry,taxID:r[2].taxID}, b) yield rel return count(*)"
)
cypher(db,query)

query <- paste0("
                MATCH (l)-[r]->(:metBiGG:recon3D {id:'2hydog_c'}) where 'reBiGG' in labels(l)
                with collect([type(r),l.sid,properties(r)]) as rels
                unwind rels as r
                MATCH (m:metBiGG:recon3D {compartment:'Cytoplasm'}) where m.enantiomere='L'
                MATCH (b:reBiGG:recon3D {sid:r[1]}) 
                CALL apoc.create.relationship(b, r[0], {source:r[2].source, stoichiometry:r[2].stoichiometry,taxID:r[2].taxID}, m) yield rel return count(*)"
)
cypher(db,query)

# delete 2HG_c
query <- paste0("
                MATCH (m:metBiGG:recon3D {id:'2hydog_c'})
                detach delete m")

cypher(db,query)

# 3
#####
# query <- paste0("
#                 MATCH (h:metBiGG:recon3D {id:'L2hydog_m'}), (r:reBiGG:recon3D {id:'HMR_0718'}),(a:metBiGG:recon3D {id:'akg_m'})
#                 MERGE (h)-[:PARTICIPATES]->(r)-[:PARTICIPATES]->(a)"
# )
# cypher(db,query)

# 4
#####
query <- paste0("
                MATCH (h:metBiGG:recon3D {id:'L2hydog_m'}),(a:metBiGG:recon3D {id:'akg_m'})
                MERGE (r:reBiGG:recon3D {sid:'MDH2m_Recon3D',id:'MDH2m'})
                MERGE (a)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(h)"
)
cypher(db,query)

# 5
#####
query <- paste0("
                MATCH (:geneSymbol {sid:'MDH2'})-[:SYNONYM]-(g:gene {source:'ncbi'})
                MATCH (r:reBiGG:recon3D {id:'MDH2m'})
                MERGE (g)-[:CATALYZES {source:'manual fix', taxID:'9606'}]->(r)"
)
cypher(db,query)

# 6
#####
query <- paste0("
                MATCH (h:metBiGG:recon3D {id:'D2hydog_m'}),(a:metBiGG:recon3D {id:'akg_m'})
                MERGE (r:reBiGG:recon3D {sid:'D2HGDHm_Recon3D',id:'D2HGDHm'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)
# 7
#####
query <- paste0("
                MATCH (:geneSymbol {sid:'D2HGDH'})-[:SYNONYM]-(g:gene {source:'ncbi'})
                MATCH (r:reBiGG:recon3D {id:'D2HGDHm'})
                MERGE (g)-[:CATALYZES {source:'manual fix', taxID:'9606'}]->(r)"
)
cypher(db,query)
# 8
#####
query <- paste0("
                MATCH (h:metBiGG:recon3D {id:'fad_m'}),(a:metBiGG:recon3D {id:'fadh2_m'})
                MATCH (r:reBiGG:recon3D {id:'D2HGDHm'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)
# 9
#####
query <- paste0("
                MATCH (:metBiGG:recon3D {id:'fadh2_c'})-[r2]-(:reBiGG:recon3D {id:'HMR_0718'})-[r1]-(:metBiGG:recon3D {id:'fad_c'})
                MATCH (h:metBiGG:recon3D {id:'fad_m'}),(a:metBiGG:recon3D {id:'fadh2_m'})
                MATCH (r:reBiGG:recon3D {id:'HMR_0718'})
                MERGE (h)-[:PARTICIPATES {source:r1.source, stoichiometry:r1.stoichiometry,taxID:r1.taxID}]->(r)-[:PARTICIPATES {source:r2.source, stoichiometry:r2.stoichiometry,taxID:r2.taxID}]->(a)
                delete r1,r2"
)
cypher(db,query)
# 10
#####
# query <- paste0("
#                 MATCH (h:metBiGG:recon3D {id:'fad_m'}),(a:metBiGG:recon3D {id:'fadh2_m'})
#                 MATCH (r:reBiGG:recon3D {id:'HMR_0718'})
#                 MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
#                 )
# cypher(db,query)
# 11
#####
query <- paste0("
                MATCH (h:metBiGG:recon3D {id:'nadh_m'}),(a:metBiGG:recon3D {id:'nad_m'})
                MATCH (r:reBiGG:recon3D {id:'MDH2m'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)
# 12
#####
query <- paste0("
                MATCH (h:metBiGG:recon3D {id:'h_m'})
                MATCH (r:reBiGG:recon3D {id:'MDH2m'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)"
)
cypher(db,query)
# 13
#####
query <- paste0("
                MATCH (h:metBiGG:recon3D {id:'4ohbut_m'}),(a:metBiGG:recon3D {id:'sucsal_m'})
                MERGE (r:reBiGG:recon3D {sid:'ADHFE1m_Recon3D',id:'ADHFE1m'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)
# 14
#####
query <- paste0("
                MATCH (:geneSymbol {sid:'ADHFE1'})-[:SYNONYM]-(g:gene {source:'ncbi'})
                MATCH (r:reBiGG:recon3D {id:'ADHFE1m'})
                MERGE (g)-[:CATALYZES {source:'manual fix', taxID:'9606'}]->(r)"
)
cypher(db,query)

# 15
#####
query <- paste0("
                MATCH (h:metBiGG:recon3D {id:'D2hydog_m'}),(a:metBiGG:recon3D {id:'akg_m'})
                MATCH (r:reBiGG:recon3D {id:'ADHFE1m'})
                MERGE (a)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(h)"
)
cypher(db,query)


####################
# Cytosol
####################

# 1
#####
query <- paste0("
                MATCH (h:metBiGG:recon3D {id:'L2hydog_c'}),(a:metBiGG:recon3D {id:'akg_c'})
                MERGE (r:reBiGG:recon3D {sid:'MDH1c_Recon3D',id:'MDH1c'})
                MERGE (a)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(h)"
)
cypher(db,query)

# 2
#####
query <- paste0("
                MATCH (:geneSymbol {sid:'MDH1'})-[:SYNONYM]-(g:gene {source:'ncbi'})
                MATCH (r:reBiGG:recon3D {id:'MDH1c'})
                MERGE (g)-[:CATALYZES {source:'manual fix', taxID:'9606'}]->(r)"
)
cypher(db,query)

# 3
#####
query <- paste0("
                MATCH (h:metBiGG:recon3D {id:'nadh_c'}),(a:metBiGG:recon3D {id:'nad_c'})
                MATCH (r:reBiGG:recon3D {id:'MDH1c'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)

# 4
#####
query <- paste0("
                MATCH (h:metBiGG:recon3D {id:'h_c'})
                MATCH (r:reBiGG:recon3D {id:'MDH1c'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)"
)
cypher(db,query)

# 5
#####
query <- paste0("
                MATCH (h:metBiGG:recon3D {id:'L2hydog_c'}),(a:metBiGG:recon3D {id:'akg_c'})
                MERGE (r:reBiGG:recon3D {sid:'LDHAc_Recon3D',id:'LDHAc'})
                MERGE (a)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(h)"
)
cypher(db,query)

# 6
#####
query <- paste0("
                MATCH (:geneSymbol {sid:'LDHA'})-[:SYNONYM]-(g:gene {source:'ncbi'})
                MATCH (r:reBiGG:recon3D {id:'LDHAc'})
                MERGE (g)-[:CATALYZES {source:'manual fix', taxID:'9606'}]->(r)"
)
cypher(db,query)

# 7
#####
query <- paste0("
                MATCH (h:metBiGG:recon3D {id:'nadh_c'}),(a:metBiGG:recon3D {id:'nad_c'})
                MATCH (r:reBiGG:recon3D {id:'LDHAc'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)

# 8
#####
query <- paste0("
                MATCH (h:metBiGG:recon3D {id:'h_c'})
                MATCH (r:reBiGG:recon3D {id:'LDHAc'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)"
)
cypher(db,query)



####################
# Transport
####################

# 1a
#####

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'cl_e'}),(h:metBiGG:recon3D {id:'cl_c'})
                MERGE (r:reBiGG:recon3D {sid:'L2hydogOAT4cle_Recon3D',id:'L2hydogOAT4cle'})
                MERGE (a)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(h)"
)
cypher(db,query)

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'L2hydog_e'}),(h:metBiGG:recon3D {id:'L2hydog_c'})
                MERGE (r:reBiGG:recon3D {sid:'L2hydogOAT4cle_Recon3D',id:'L2hydogOAT4cle'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)

query <- paste0("
                MATCH (:geneSymbol {sid:'OAT4'})-[:SYNONYM]-(g:gene {source:'ncbi'})
                MATCH (r:reBiGG:recon3D {id:'L2hydogOAT4cle'})
                MERGE (g)-[:CATALYZES {source:'manual fix', taxID:'9606'}]->(r)"
)
cypher(db,query)
##

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'D2hydog_e'}),(h:metBiGG:recon3D {id:'D2hydog_c'})
                MERGE (r:reBiGG:recon3D {sid:'D2hydogOAT4cle_Recon3D',id:'D2hydogOAT4cle'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'cl_e'}),(h:metBiGG:recon3D {id:'cl_c'})
                MERGE (r:reBiGG:recon3D {sid:'D2hydogOAT4cle_Recon3D',id:'D2hydogOAT4cle'})
                MERGE (a)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(h)"
)
cypher(db,query)

query <- paste0("
                MATCH (:geneSymbol {sid:'OAT4'})-[:SYNONYM]-(g:gene {source:'ncbi'})
                MATCH (r:reBiGG:recon3D {id:'D2hydogOAT4cle'})
                MERGE (g)-[:CATALYZES {source:'manual fix', taxID:'9606'}]->(r)"
)
cypher(db,query)
##

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'3hydog_e'}),(h:metBiGG:recon3D {id:'3hydog_c'})
                MERGE (r:reBiGG:recon3D {sid:'3hydogOAT4cle_Recon3D',id:'3hydogOAT4cle'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'cl_e'}),(h:metBiGG:recon3D {id:'cl_c'})
                MERGE (r:reBiGG:recon3D {sid:'3hydogOAT4cle_Recon3D',id:'3hydogOAT4cle'})
                MERGE (a)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(h)"
)
cypher(db,query)

query <- paste0("
                MATCH (:geneSymbol {sid:'OAT4'})-[:SYNONYM]-(g:gene {source:'ncbi'})
                MATCH (r:reBiGG:recon3D {id:'3hydogOAT4cle'})
                MERGE (g)-[:CATALYZES {source:'manual fix', taxID:'9606'}]->(r)"
)
cypher(db,query)

# 1b
#####

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'estrones_e'}),(h:metBiGG:recon3D {id:'estrones_c'})
                MERGE (r:reBiGG:recon3D {sid:'L2hydogOAT4este_Recon3D',id:'L2hydogOAT4este'})
                MERGE (a)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(h)"
)
cypher(db,query)


query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'L2hydog_e'}),(h:metBiGG:recon3D {id:'L2hydog_c'})
                MERGE (r:reBiGG:recon3D {sid:'L2hydogOAT4este_Recon3D',id:'L2hydogOAT4este'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)

query <- paste0("
                MATCH (:geneSymbol {sid:'OAT4'})-[:SYNONYM]-(g:gene {source:'ncbi'})
                MATCH (r:reBiGG:recon3D {id:'L2hydogOAT4este'})
                MERGE (g)-[:CATALYZES {source:'manual fix', taxID:'9606'}]->(r)"
)
cypher(db,query)
##
query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'estrones_e'}),(h:metBiGG:recon3D {id:'estrones_c'})
                MERGE (r:reBiGG:recon3D {sid:'D2hydogOAT4este_Recon3D',id:'D2hydogOAT4este'})
                MERGE (a)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(h)"
)
cypher(db,query)

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'D2hydog_e'}),(h:metBiGG:recon3D {id:'D2hydog_c'})
                MERGE (r:reBiGG:recon3D {sid:'D2hydogOAT4este_Recon3D',id:'D2hydogOAT4este'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)

query <- paste0("
                MATCH (:geneSymbol {sid:'OAT4'})-[:SYNONYM]-(g:gene {source:'ncbi'})
                MATCH (r:reBiGG:recon3D {id:'D2hydogOAT4este'})
                MERGE (g)-[:CATALYZES {source:'manual fix', taxID:'9606'}]->(r)"
)
cypher(db,query)
##
query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'estrones_e'}),(h:metBiGG:recon3D {id:'estrones_c'})
                MERGE (r:reBiGG:recon3D {sid:'3hydogOAT4este_Recon3D',id:'3hydogOAT4este'})
                MERGE (a)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(h)"
)
cypher(db,query)

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'3hydog_e'}),(h:metBiGG:recon3D {id:'3hydog_c'})
                MERGE (r:reBiGG:recon3D {sid:'3hydogOAT4este_Recon3D',id:'3hydogOAT4este'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)

query <- paste0("
                MATCH (:geneSymbol {sid:'OAT4'})-[:SYNONYM]-(g:gene {source:'ncbi'})
                MATCH (r:reBiGG:recon3D {id:'3hydogOAT4este'})
                MERGE (g)-[:CATALYZES {source:'manual fix', taxID:'9606'}]->(r)"
)
cypher(db,query)

# 2
#####

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'akg_c'}),(h:metBiGG:recon3D {id:'akg_e'})
                MERGE (r:reBiGG:recon3D {sid:'L2hydogOAT1e_Recon3D',id:'L2hydogOAT1e'})
                MERGE (a)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(h)"
)
cypher(db,query)


query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'L2hydog_c'}),(h:metBiGG:recon3D {id:'L2hydog_e'})
                MERGE (r:reBiGG:recon3D {sid:'L2hydogOAT1e_Recon3D',id:'L2hydogOAT1e'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)

query <- paste0("
                MATCH (:geneSymbol {sid:'OAT1'})-[:SYNONYM]-(g:gene {source:'ncbi'})
                MATCH (r:reBiGG:recon3D {id:'L2hydogOAT1e'})
                MERGE (g)-[:CATALYZES {source:'manual fix', taxID:'9606'}]->(r)"
)
cypher(db,query)
##
query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'akg_c'}),(h:metBiGG:recon3D {id:'akg_e'})
                MERGE (r:reBiGG:recon3D {sid:'D2hydogOAT1e_Recon3D',id:'D2hydogOAT1e'})
                MERGE (a)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(h)"
)
cypher(db,query)

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'D2hydog_c'}),(h:metBiGG:recon3D {id:'D2hydog_e'})
                MERGE (r:reBiGG:recon3D {sid:'D2hydogOAT1e_Recon3D',id:'D2hydogOAT1e'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)

query <- paste0("
                MATCH (:geneSymbol {sid:'OAT1'})-[:SYNONYM]-(g:gene {source:'ncbi'})
                MATCH (r:reBiGG:recon3D {id:'D2hydogOAT1e'})
                MERGE (g)-[:CATALYZES {source:'manual fix', taxID:'9606'}]->(r)"
)
cypher(db,query)
##
query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'akg_c'}),(h:metBiGG:recon3D {id:'akg_e'})
                MERGE (r:reBiGG:recon3D {sid:'3hydogOAT1e_Recon3D',id:'3hydogOAT1e'})
                MERGE (a)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(h)"
)
cypher(db,query)

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'3hydog_c'}),(h:metBiGG:recon3D {id:'3hydog_e'})
                MERGE (r:reBiGG:recon3D {sid:'3hydogOAT1e_Recon3D',id:'3hydogOAT1e'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)

query <- paste0("
                MATCH (:geneSymbol {sid:'OAT1'})-[:SYNONYM]-(g:gene {source:'ncbi'})
                MATCH (r:reBiGG:recon3D {id:'3hydogOAT1e'})
                MERGE (g)-[:CATALYZES {source:'manual fix', taxID:'9606'}]->(r)"
)
cypher(db,query)

# 3
#####

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'na_c'}),(h:metBiGG:recon3D {id:'na_e'})
                MERGE (r:reBiGG:recon3D {sid:'3hydogNADC3e_Recon3D',id:'3hydogNADC3e'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-3',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'3',taxID:'9606'}]->(a)"
)
cypher(db,query)

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'3hydog_c'}),(h:metBiGG:recon3D {id:'3hydog_e'})
                MERGE (r:reBiGG:recon3D {sid:'3hydogNADC3e_Recon3D',id:'3hydogNADC3e'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)

query <- paste0("
                MATCH (:geneSymbol {sid:'NADC3'})-[:SYNONYM]-(g:gene {source:'ncbi'})
                MATCH (r:reBiGG:recon3D {id:'3hydogNADC3e'})
                MERGE (g)-[:CATALYZES {source:'manual fix', taxID:'9606'}]->(r)"
)
cypher(db,query)
##

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'na_c'}),(h:metBiGG:recon3D {id:'na_e'})
                MERGE (r:reBiGG:recon3D {sid:'L2hydogNADC3e_Recon3D',id:'L2hydogNADC3e'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-3',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'3',taxID:'9606'}]->(a)"
)
cypher(db,query)

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'L2hydog_c'}),(h:metBiGG:recon3D {id:'L2hydog_e'})
                MERGE (r:reBiGG:recon3D {sid:'L2hydogNADC3e_Recon3D',id:'L2hydogNADC3e'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)

query <- paste0("
                MATCH (:geneSymbol {sid:'NADC3'})-[:SYNONYM]-(g:gene {source:'ncbi'})
                MATCH (r:reBiGG:recon3D {id:'L2hydogNADC3e'})
                MERGE (g)-[:CATALYZES {source:'manual fix', taxID:'9606'}]->(r)"
)
cypher(db,query)
##
query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'na_c'}),(h:metBiGG:recon3D {id:'na_e'})
                MERGE (r:reBiGG:recon3D {sid:'D2hydogNADC3e_Recon3D',id:'D2hydogNADC3e'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-3',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'3',taxID:'9606'}]->(a)"
)
cypher(db,query)

query <- paste0("
                MATCH (a:metBiGG:recon3D {id:'D2hydog_c'}),(h:metBiGG:recon3D {id:'D2hydog_e'})
                MERGE (r:reBiGG:recon3D {sid:'D2hydogNADC3e_Recon3D',id:'D2hydogNADC3e'})
                MERGE (h)-[:PARTICIPATES {source:'manual fix', stoichiometry:'-1',taxID:'9606'}]->(r)-[:PARTICIPATES {source:'manual fix', stoichiometry:'1',taxID:'9606'}]->(a)"
)
cypher(db,query)


query <- paste0("
                MATCH (:geneSymbol {sid:'NADC3'})-[:SYNONYM]-(g:gene {source:'ncbi'})
                MATCH (r:reBiGG:recon3D {id:'D2hydogNADC3e'})
                MERGE (g)-[:CATALYZES {source:'manual fix', taxID:'9606'}]->(r)"
)
cypher(db,query)

