library(DBI)

# =====================
# VARIANTS POR GEN (DINÁMICO)
# =====================
get_variants_by_gene <- function(con, gene){
  
  if(is.null(gene) || gene == ""){
    return(data.frame())
  }
  
  files <- list.files("../data/variants", full.names = TRUE)
  
  if(length(files) == 0){
    stop("No variant files found")
  }
  
  has_wes <- any(grepl("wes", files))
  has_wgs <- any(grepl("wgs", files))
  
  # =====================
  # WES + WGS
  # =====================
  if(has_wes & has_wgs){
    
    query <- paste0("
      
      WITH wes AS (
        SELECT *, 'WES' AS src FROM read_parquet('../data/variants/wes*.parquet')
      ),
      wgs AS (
        SELECT *, 'WGS' AS src FROM read_parquet('../data/variants/wgs*.parquet')
      ),
      combined AS (
        SELECT * FROM wes
        UNION ALL
        SELECT * FROM wgs
      ),
      
      normalized AS (
        SELECT *,
        
          REPLACE(CHILD_GT,'|','/') AS CHILD_GT_N,
          REPLACE(PARENT1_GT,'|','/') AS P1_GT_N,
          REPLACE(PARENT2_GT,'|','/') AS P2_GT_N
          
        FROM combined
      )
      
      SELECT DISTINCT
        ID, CHROM, POS, REF, ALT,
        UPPER(SYMBOL) AS \"Gene name\",
        Gene AS \"Gene ID\",
        OMIM_id
        
      FROM normalized
      
      WHERE LOWER(SYMBOL) LIKE LOWER('%", gene, "%')
      
    ")
    
  } else if(has_wes){
    
    query <- paste0("
      SELECT DISTINCT
        ID, CHROM, POS,
        UPPER(SYMBOL) AS \"Gene name\",
        Gene AS \"Gene ID\",
        OMIM_id
      FROM read_parquet('../data/variants/wes*.parquet')
      WHERE LOWER(SYMBOL) LIKE LOWER('%", gene, "%')
    ")
    
  } else if(has_wgs){
    
    query <- paste0("
      SELECT DISTINCT
        ID, CHROM, POS,
        UPPER(SYMBOL) AS \"Gene name\",
        Gene AS \"Gene ID\",
        OMIM_id
      FROM read_parquet('../data/variants/wgs*.parquet')
      WHERE LOWER(SYMBOL) LIKE LOWER('%", gene, "%')
    ")
  }
  
  dbGetQuery(con, query)
}


# =====================
# INFO DEL GEN
# =====================
get_gene_info_by_gene <- function(con, gene){
  
  if(is.null(gene) || gene == ""){
    return(data.frame())
  }
  
  query <- paste0("
    
    SELECT DISTINCT
      UPPER(SYMBOL) AS \"Gene name\",
      Gene AS \"Gene ID\",
      HGNC_ID,
      BIOTYPE,
      GENE_PHENO,
      Function_description,
      Disease_description,
      HPO_id,
      HPO_name,
      OMIM_id
      
    FROM read_parquet('../data/variants/*.parquet')
    
    WHERE LOWER(SYMBOL) LIKE LOWER('%", gene, "%')
    
  ")
  
  dbGetQuery(con, query)
}


# =====================
# LISTA DE GENES (para autocomplete)
# =====================
get_all_genes <- function(con){
  
  dbGetQuery(con, "
    SELECT DISTINCT UPPER(SYMBOL) AS gene
    FROM read_parquet('../data/variants/*.parquet')
    ORDER BY gene
  ")
}