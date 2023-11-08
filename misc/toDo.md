## 20.09.17 - Meeting Ontology

- ontologie Definition
    - Knotenliste 
        - eine Zeile ist ein Metabolite
        - Felder: interner Name (unique id), Display Name, HMDB ID (wird benutzt um HMDB zu linken)
    - Kanten
        - VON (interner name), ZU (interner name), ANNOTATION STRING
- metMetabolon zu Ontologie
    - Kanten
        - VON (metMetabolon id), ZU (Ontologie interner Name)


-----

### Maria toDO

- ~~alle metMetabolons die gesehen wurden~~
    - ~~sid, Biochemical, HMDB~~
- ~~alle metHMDB~~
    - ~~parser umschreiben um Ontologien (5? Felder) zu parsen~~
        - ~~set Flag to parse Ontologies~~