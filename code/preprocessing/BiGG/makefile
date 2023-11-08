
RESULTDIST = data/processed/
MODELDIST = data/raw/
FILES = _genes.csv _met_re.csv _metabolites.csv _re_genes.csv _reactions.csv
RECONMODELS := $(wildcard $(MODELDIST)recon3D/*.xml)
RECONRDATA := $(patsubst $(MODELDIST)recon3D/%.xml, $(RESULTDIST)recon3D/%.Rdata,$(RECONMODELS))
RECONRESULT := $(foreach f,$(FILES),$(patsubst $(RESULTDIST)recon3D/%.Rdata, $(RESULTDIST)recon3D/%$(f),$(RECONRDATA)))
AGORAMODELS := $(wildcard $(MODELDIST)agora/*.xml)
AGORARDATA := $(patsubst $(MODELDIST)agora/%.xml, $(RESULTDIST)agora/%.Rdata,$(AGORAMODELS))
AGORARESULT := $(foreach f,$(FILES),$(patsubst $(RESULTDIST)agora/%.Rdata, $(RESULTDIST)agora/%$(f),$(AGORARDATA)))
AGORAINFO := $(patsubst $(MODELDIST)agora/%.xml, $(RESULTDIST)agora_info/%_info.csv,$(AGORAMODELS))

.PHONY: all

all: recon agora
########################################################################
## RECON
########################################################################

recon: \
	initiate_recon \
	genes \
  $(RECONRDATA) \
  $(RECONRESULT)

## INTITIATE
initiate_recon:
	mkdir -p data/raw/recon3D
	@echo "Deposit xml files in created folder and run 'make recon'"


## GENES
genes: \
  data/raw/hgnc_complete_set.txt\
  data/processed/genes_hgnc.rds

data/raw/hgnc_complete_set.txt:
	mkdir -p data/raw
	wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt -O $@

data/processed/genes_hgnc.rds: scripts/BiGG_01_hgnc.R\
  data/raw/hgnc_complete_set.txt
	mkdir -p data/processed
	Rscript $^ $@

## PARSE
$(RECONRDATA): scripts/BiGG_01_parse.R data/processed/genes_hgnc.rds
	@if [ ! -f $@ ]; then \
			mkdir -p data/processed/recon3D; \
			Rscript $^ $(patsubst $(RESULTDIST)recon3D/%.Rdata, $(MODELDIST)recon3D/%.xml,$@) $@; \
	fi;

## RESULTS
$(RECONRESULT): makeResultsR

makeResultsR: scripts/BiGG_02_neo4j.R $(RECONRDATA)
	mkdir -p data/processed/recon3D
	@for i in $(foreach f, $(RECONRDATA), $(shell basename -a -s .Rdata $f )) ; do \
		if [ ! -f $(RESULTDIST)recon3D/$${i}_met_re.csv -o ! -f $(RESULTDIST)recon3D/$${i}_genes.csv -o ! -f $(RESULTDIST)recon3D/$${i}_reactions.csv -o ! -f $(RESULTDIST)recon3D/$${i}_metabolites.csv -o ! -f $(RESULTDIST)recon3D/$${i}_re_gene.csv ]; then \
			Rscript scripts/BiGG_02_neo4j.R $(RESULTDIST)recon3D/$${i}.Rdata $${i} data/20190411_recon_genes.tsv data/20190411_recon_reactions.tsv data/20190411_recon_metabolites.tsv; \
		fi; \
	done;



########################################################################
## AGORA
########################################################################

agora: \
	initiate_agora \
	genes \
  $(AGORARDATA) \
  $(AGORARESULT) \
	$(AGORAINFO)

## INTITIATE
initiate_agora:
	mkdir -p data/raw/agora
	@echo "Deposit xml files in created folder and run 'make agora'"

## PARSE
$(AGORARDATA): scripts/BiGG_01_parse.R data/processed/genes_hgnc.rds
	@if [ ! -f $@ ]; then \
			mkdir -p data/processed/agora; \
			Rscript $^ $(patsubst $(RESULTDIST)agora/%.Rdata, $(MODELDIST)agora/%.xml,$@) $@; \
	fi;

## RESULTS
$(AGORARESULT): makeResults

makeResults: scripts/BiGG_02_neo4j.R $(AGORARDATA)
	mkdir -p data/processed/agora
	@for i in $(foreach f, $(AGORARDATA), $(shell basename -a -s .Rdata $f )) ; do \
		if [ ! -f $(RESULTDIST)agora/$${i}_met_re.csv -o ! -f $(RESULTDIST)agora/$${i}_genes.csv -o ! -f $(RESULTDIST)agora/$${i}_reactions.csv -o ! -f $(RESULTDIST)agora/$${i}_metabolites.csv -o ! -f $(RESULTDIST)agora/$${i}_re_gene.csv ]; then \
			Rscript scripts/BiGG_02_neo4j.R $(RESULTDIST)agora/$${i}.Rdata $${i}; \
		fi; \
	done;

$(AGORAINFO): makeInfo

makeInfo: scripts/BiGG_03_agora_info.R data/20190429_microbes_info.tsv
	mkdir -p data/processed/agora_info
	Rscript scripts/BiGG_03_agora_info.R data/20190429_microbes_info.tsv $(RESULTDIST)agora_info/;
