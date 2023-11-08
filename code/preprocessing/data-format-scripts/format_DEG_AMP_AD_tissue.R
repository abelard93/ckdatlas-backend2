library(jsonlite)

### ----------------------- LOAD DATA ----------------------------------
test <- fromJSON("DEG/rnaseq_differential_expression.json", simplifyVector = FALSE)

### ----------------------- FORMAT -------------------------------------
df <- rbindlist(test, fill = TRUE, idcol = F) 
# rename
colnames(df) <- c("gene_id","gene_symbol","logFC","FC","CI_l","CI_r","pvalue_fdr","tissue", "study","trait")
df$direction <- "unchanged"
df$direction[df$pvalue_fdr<=0.05&df$logFC>0] <- "up"
df$direction[df$pvalue_fdr<=0.05&df$logFC<0] <- "down"

### ----------------------- WRITE -------------------------------------
write.csv(df,"DEG/logsdon_deg_tissue_specific.csv",row.names = F,quote = T,na = "")

# GENERATE INFO FILE
data = c(
  "file_name","logsdon_deg_tissue_specific.csv",
  "sid", "logsdon_deg_tissue_specific",
  "organism", "Homo Sapiens",
  "node_type", "DEG",
  "cohort", "Mayo|MSBB|ROSMAP",
  "publication","logsdon_2019",
  "significance","FDR(p)<=0.05",
  "tax","hsa",
  "taxID","9606"
)
# make df
info <- data.frame(matrix(data = data, ncol = 2, byrow = T))
# OUTPUT
write.table(
  info,
  file = "DEG/logsdon_deg_tissue_specific_info.csv",
  sep = ",",
  col.names = F,
  row.names = F
)
