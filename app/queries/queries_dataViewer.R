library(DBI)

# =====================
# VARIANTS BASE (SIN FILTROS)
# =====================
get_variants_with_inheritance <- function(con){
  
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
    
    query <- "
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
        
          CASE 
            WHEN REPLACE(CHILD_GT, '|', '/') IN ('1/0','0/1') THEN '0/1'
            ELSE REPLACE(CHILD_GT, '|', '/')
          END AS CHILD_GT_N,
          
          CASE 
            WHEN REPLACE(PARENT1_GT, '|', '/') IN ('1/0','0/1') THEN '0/1'
            ELSE REPLACE(PARENT1_GT, '|', '/')
          END AS P1_GT_N,
          
          CASE 
            WHEN REPLACE(PARENT2_GT, '|', '/') IN ('1/0','0/1') THEN '0/1'
            ELSE REPLACE(PARENT2_GT, '|', '/')
          END AS P2_GT_N
          
        FROM combined
      )
      
      SELECT 
        
        ID,
        CHROM,
        POS,
        REF,
        ALT,
        FILTER,
        
        CHILD_GT_N AS GT,
        CHILD_DP AS DP,
        CHILD_AD AS AD,
        CHILD_GQ AS GQ,
        
        UPPER(SYMBOL) AS \"Gene name\",
        Gene AS \"Gene ID\",
        
        Consequence,
        IMPACT,
        MAX_AF,
        VARIANT_CLASS,
        
        SIFT_pred,
        Polyphen2_HVAR_pred,
        
        -- SOURCE
        CASE
          WHEN COUNT(DISTINCT src) OVER (PARTITION BY ID) = 2 THEN 'BOTH'
          ELSE MAX(src) OVER (PARTITION BY ID)
        END AS source,
        
        -- INHERITANCE TYPE
        CASE
          WHEN CHILD_GT_N = '0/1' AND P1_GT_N = '0/0' AND P2_GT_N = '0/0' THEN 'de_novo'
          WHEN CHILD_GT_N = '1/1' AND P1_GT_N = '0/1' AND P2_GT_N = '0/1' THEN 'recessive'
          WHEN CHILD_GT_N = '0/1' AND (P1_GT_N = '0/1' OR P2_GT_N = '0/1') THEN 'dominant'
          ELSE 'other'
        END AS inheritance_type,
        
        -- INHERITANCE TEXT
        CONCAT(
          CASE
            WHEN CHILD_GT_N = '0/1' AND P1_GT_N = '0/0' AND P2_GT_N = '0/0' THEN 'de_novo'
            WHEN CHILD_GT_N = '1/1' AND P1_GT_N = '0/1' AND P2_GT_N = '0/1' THEN 'recessive'
            WHEN CHILD_GT_N = '0/1' AND (P1_GT_N = '0/1' OR P2_GT_N = '0/1') THEN 'dominant'
            ELSE 'other'
          END,
          ' (P1:', P1_GT_N, ' AD:', PARENT1_AD,
          ' | P2:', P2_GT_N, ' AD:', PARENT2_AD, ')'
        ) AS inheritance
        
      FROM (
        SELECT DISTINCT *
        FROM normalized
      )
    "
    
  } else if(has_wes){
    
    query <- "
      WITH normalized AS (
        SELECT *,
        REPLACE(CHILD_GT,'|','/') AS CHILD_GT_N,
        REPLACE(PARENT1_GT,'|','/') AS P1_GT_N,
        REPLACE(PARENT2_GT,'|','/') AS P2_GT_N
        FROM read_parquet('../data/variants/wes*.parquet')
      )
      
      SELECT *,
        'WES' AS source
      FROM (
        SELECT DISTINCT *
        FROM normalized
    "
    
  } else if(has_wgs){
    
    query <- "
      WITH normalized AS (
        SELECT *,
        REPLACE(CHILD_GT,'|','/') AS CHILD_GT_N,
        REPLACE(PARENT1_GT,'|','/') AS P1_GT_N,
        REPLACE(PARENT2_GT,'|','/') AS P2_GT_N
        FROM read_parquet('../data/variants/wgs*.parquet')
      )
      
      SELECT *,
        'WGS' AS source
      FROM (
        SELECT DISTINCT *
        FROM normalized
    "
  }
  
  dbGetQuery(con, query)
}


# =====================
# 🚀 VARIANTS FILTRADOS (NUEVO)
# =====================
get_variants_filtered <- function(con,
                                  gene = NULL,
                                  af = 0.01,
                                  impact = NULL,
                                  source = NULL,
                                  inheritance = NULL){
  
  base <- "SELECT * FROM ("
  
  query_main <- get_variants_with_inheritance_sql()
  
  query <- paste0(base, query_main, ") WHERE 1=1")
  
  if(!is.null(gene) && gene != ""){
    query <- paste0(query, " AND LOWER(\"Gene name\") LIKE LOWER('%", gene, "%')")
  }
  
  if(!is.null(af)){
    query <- paste0(query, " AND (MAX_AF IS NULL OR MAX_AF <= ", af, ")")
  }
  
  if(!is.null(impact)){
    query <- paste0(query, " AND IMPACT IN ('", paste(impact, collapse="','"), "')")
  }
  
  if(!is.null(source)){
    query <- paste0(query, " AND source IN ('", paste(source, collapse="','"), "')")
  }
  
  if(!is.null(inheritance)){
    query <- paste0(query, " AND inheritance_type IN ('", paste(inheritance, collapse="','"), "')")
  }
  
  dbGetQuery(con, query)
}


# =====================
# RNA
# =====================
get_rna <- function(con){
  
  files_sample <- list.files("../data/rnaseq", pattern="sample_", full.names=TRUE)
  files_ctrl   <- list.files("../data/rnaseq", pattern="controls", full.names=TRUE)
  
  if(length(files_sample) == 0){
    stop("No RNA sample data")
  }
  
  # SIN CONTROLES
  if(length(files_ctrl) == 0){
    
    return(dbGetQuery(con, "
      SELECT 
        gene_id,
        gene_name,
        tpm AS \"gene tpm\"
      FROM read_parquet('../data/rnaseq/sample_*.parquet')
    "))
  }
  
  # CON CONTROLES
  dbGetQuery(con, "
    SELECT 
      r.gene_id,
      r.gene_name,
      r.tpm AS \"gene tpm\",
      
      c.max_TPM,
      c.min_TPM,
      c.mean_TPM,
      c.median_TPM
      
    FROM read_parquet('../data/rnaseq/sample_*.parquet') r
    LEFT JOIN read_parquet('../data/rnaseq/controls*.parquet') c
    USING(gene_id)
  ")
}


# =====================
# DROP
# =====================
get_drop_expr <- function(con){
  dbGetQuery(con, "SELECT * FROM read_parquet('../data/drop/expression*.parquet')")
}

get_drop_splicing <- function(con){
  dbGetQuery(con, "SELECT * FROM read_parquet('../data/drop/splicing*.parquet')")
}

get_drop_mae <- function(con){
  dbGetQuery(con, "SELECT * FROM read_parquet('../data/drop/mae*.parquet')")
}