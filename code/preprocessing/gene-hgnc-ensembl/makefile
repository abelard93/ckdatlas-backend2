.PHONY: all ensembl hgnc ncbi uniprot#corum

all:\
  hgnc\
  ensembl\
  ncbi \
	uniprot
  #corum

########################################################################
## HGNC
########################################################################
hgnc: \
  data/raw/hgnc_complete_set.txt\
  data/processed/hgnc_entrez.csv\
  data/processed/hgnc_ensembl.csv

data/raw/hgnc_complete_set.txt:
	mkdir -p data/raw
	wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt -O $@

data/processed/hgnc_entrez.csv data/processed/hgnc_ensembl.csv: scripts/genes_01_hgnc.R\
  data/raw/hgnc_complete_set.txt
	mkdir -p data/processed
	Rscript $^ data/processed/hgnc


########################################################################
## UNIPROT
########################################################################
uniprot: \
	data/raw/HUMAN_9606_idmapping_selected.tab\
	data/processed/uniprot.csv

data/raw/HUMAN_9606_idmapping_selected.tab:
		mkdir -p data/raw
		wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz -O $@.gz
		gunzip $@.gz

data/processed/uniprot.csv: data/raw/HUMAN_9606_idmapping_selected.tab
		mkdir -p data/processed
		echo "uniprot_id,ensembl_id" > $@
		cat $^ | cut -f 1,21 | tr "\\t" "," | tr ";" "|" | tr -d '[:blank:]' >> $@

########################################################################
## ENSEMBL
########################################################################
ensembl: \
	data/processed/ensembl_gene2transcript2protein.csv\
	data/processed/ensembl_ensembl2entrez.csv\
	data/processed/ensembl_gene.csv\
	data/processed/ensembl_transcript.csv\
	data/processed/ensembl_protein.csv


data/processed/ensembl_gene2transcript2protein.csv data/processed/ensembl_ensembl2entrez.csv data/processed/ensembl_gene.csv data/processed/ensembl_transcript.csv data/processed/ensembl_protein.csv: \
	scripts/genes_02_ensembl.R
	mkdir -p data/processed
	Rscript $^ data/processed/ensembl_



########################################################################
## NCBI
########################################################################
ncbi: \
	data/raw/Homo_sapiens.gene_info\
	data/processed/entrez.csv

data/raw/Homo_sapiens.gene_info:
	mkdir -p data/raw
	wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz -O $@.gz
	gunzip $@.gz

data/processed/entrez.csv: scripts/genes_03_entrez.R\
	data/raw/Homo_sapiens.gene_info
	mkdir -p data/processed
	Rscript $^ $@

########################################################################
## CORUM
########################################################################
#corum: \
#  data/raw/coreComplexes.txt\
#  data/processed/genes_corum.csv
#
#data/raw/coreComplexes.txt.zip:
#  mkdir -p data/raw
#  wget https://mips.helmholtz-muenchen.de/corum/download/coreComplexes.txt.zip -O $@
#
#data/raw/coreComplexes.txt: data/raw/coreComplexes.txt.zip
#  unzip $^ && rm $^
#
#data/processed/genes_corum.csv: scripts/genes_04_corum.R\
#  data/raw/coreComplexes.txt
#  mkdir -p data/processed
#  Rscript $^ $@
