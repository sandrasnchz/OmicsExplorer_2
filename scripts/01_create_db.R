library(DBI)
library(duckdb)

con <- dbConnect(duckdb(), "../db/genomic.duckdb")

# =====================
# FILES METADATA
# =====================
dbExecute(con, "
CREATE TABLE IF NOT EXISTS files (
    id INTEGER,
    file_name TEXT,
    file_type TEXT,
    path TEXT,
    uploaded_at TIMESTAMP
);
")

# =====================
# GENES
# =====================
dbExecute(con, "
CREATE TABLE IF NOT EXISTS genes (
    gene_id TEXT,
    gene_name TEXT,
    hgnc_id TEXT,
    biotype TEXT,
    function_description TEXT,
    disease_description TEXT,
    omim_id TEXT
);
")

# =====================
# GENE HPO
# =====================
dbExecute(con, "
CREATE TABLE IF NOT EXISTS gene_hpo (
    gene_id TEXT,
    hpo_id TEXT,
    hpo_name TEXT
);
")

dbDisconnect(con)