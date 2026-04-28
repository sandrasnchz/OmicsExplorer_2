library(DBI)

# =====================
# HELPER
# =====================
build_parquet_query <- function(files){
  files_sql <- paste0("'", normalizePath(files, winslash = "/"), "'", collapse = ", ")
  paste0("SELECT * FROM read_parquet([", files_sql, "])")
}


# =====================
# WES QC
# =====================
get_qc_wes <- function(con){
  
  files <- list.files("../data/qc", pattern="wes_qc", full.names=TRUE)
  
  if(length(files) == 0){
    stop("No WES QC data")
  }
  
  query <- paste0(
    "SELECT *, 'WES' AS type FROM (",
    build_parquet_query(files),
    ")"
  )
  
  dbGetQuery(con, query)
}


# =====================
# WGS QC
# =====================
get_qc_wgs <- function(con){
  
  files <- list.files("../data/qc", pattern="wgs_qc", full.names=TRUE)
  
  if(length(files) == 0){
    stop("No WGS QC data")
  }
  
  query <- paste0(
    "SELECT *, 'WGS' AS type FROM (",
    build_parquet_query(files),
    ")"
  )
  
  dbGetQuery(con, query)
}


# =====================
# COMBINED QC
# =====================
get_qc_combined <- function(con){
  
  wes_files <- list.files("../data/qc", pattern="wes_qc", full.names=TRUE)
  wgs_files <- list.files("../data/qc", pattern="wgs_qc", full.names=TRUE)
  
  queries <- c()
  
  if(length(wes_files) > 0){
    queries <- c(queries,
                 paste0("SELECT *, 'WES' AS type FROM (", build_parquet_query(wes_files), ")")
    )
  }
  
  if(length(wgs_files) > 0){
    queries <- c(queries,
                 paste0("SELECT *, 'WGS' AS type FROM (", build_parquet_query(wgs_files), ")")
    )
  }
  
  if(length(queries) == 0){
    stop("No QC data found")
  }
  
  final_query <- paste(queries, collapse = " UNION ALL ")
  
  dbGetQuery(con, final_query)
}