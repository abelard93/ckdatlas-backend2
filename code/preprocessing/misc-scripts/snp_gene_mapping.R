snp <- read.csv(file="SNP/snp_data_GRCh37.csv",header = T )
SNPids <- snp$sid

library("biomaRt")

# To show which marts are available
listMarts()

# You need the SNP mart
mart <- useMart("ENSEMBL_MART_SNP")

# Find homo sapiens
listDatasets(mart)

# This will be the dataset we want to use
dataset <- useDataset("hsapiens_snp", mart=mart)

# Show available filters
listFilters(dataset)

# Now list all available attributes
listAttributes(dataset)

# To get the ensembl gene id belonging to the SNPs
mapping <- getBM(attributes=c('refsnp_id', 'ensembl_gene_stable_id'), 
      filters = 'snp_filter', 
      values = SNPids, 
      mart = dataset)
