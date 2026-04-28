library(data.table)
library(dplyr)
library(stringr)

# =====================
# INFERENCIA HERENCIA
# =====================
infer_inheritance <- function(p1_gt, p2_gt, child_gt, p1_ad, p2_ad){
  
  # soporte alélico (ej: "10,5")
  get_alt_ratio <- function(ad){
    if(is.na(ad) | ad == "") return(NA)
    parts <- as.numeric(strsplit(ad, ",")[[1]])
    if(length(parts) != 2) return(NA)
    return(parts[2] / sum(parts))
  }
  
  r1 <- get_alt_ratio(p1_ad)
  r2 <- get_alt_ratio(p2_ad)
  
  # clasificación simple
  if(child_gt == "1/1" & p1_gt == "0/1" & p2_gt == "0/1"){
    inh <- "AR"
  } else if(child_gt == "0/1" & (p1_gt == "0/1" | p2_gt == "0/1")){
    inh <- "AD"
  } else {
    inh <- "UNKNOWN"
  }
  
  # añadir soporte
  support <- paste0(
    "P1(", p1_gt, "|", round(r1,2), ")_",
    "P2(", p2_gt, "|", round(r2,2), ")"
  )
  
  return(paste(inh, support, sep=" | "))
}

# =====================
# MAIN FUNCTION
# =====================
process_variants <- function(file, con, source){
  
  dt <- fread(file)
  
  # -------------------------
  # VARIANT ID
  # -------------------------
  dt[, variant_id := paste(CHROM, POS, REF, ALT, sep="_")]
  
  # -------------------------
  # DETECTAR COLUMNA GEN
  # -------------------------
  gene_col <- if("SYMBOL" %in% names(dt)) "SYMBOL"
  else if("gene_name" %in% names(dt)) "gene_name"
  else if("Gene" %in% names(dt)) "Gene"
  else stop("No gene column found")
  
  dt[, gene_name := toupper(get(gene_col))]
  
  # -------------------------
  # INHERITANCE
  # -------------------------
  dt[, inheritance := mapply(
    infer_inheritance,
    PARENT1_GT,
    PARENT2_GT,
    CHILD_GT,
    PARENT1_AD,
    PARENT2_AD
  )]
  
  # -------------------------
  # VARIANTS CORE
  # -------------------------
  variants_core <- dt[, .(
    variant_id,
    chrom = CHROM,
    pos = POS,
    ref = REF,
    alt = ALT,
    
    gene_id = Gene,
    gene_name,
    
    consequence = Consequence,
    impact = IMPACT,
    
    sift = SIFT_pred,
    polyphen = Polyphen2_HVAR_pred,
    
    child_gt = CHILD_GT,
    child_dp = as.numeric(CHILD_DP),
    
    inheritance,
    source
  )]
  
  variants_core <- unique(variants_core)
  
  DBI::dbWriteTable(con, "variants_core", variants_core, append=TRUE)
  
  # -------------------------
  # GENES
  # -------------------------
  genes <- dt[, .(
    gene_id = Gene,
    gene_name,
    hgnc_id = HGNC_ID,
    biotype = BIOTYPE,
    gene_pheno = GENE_PHENO,
    function_description = Function_description,
    disease_description = Disease_description,
    omim_id = OMIM_id
  )]
  
  genes <- unique(genes)
  
  DBI::dbWriteTable(con, "genes", genes, append=TRUE)
  
  # -------------------------
  # FREQUENCIES (LONG)
  # -------------------------
  af_cols <- grep("_AF$", names(dt), value=TRUE)
  
  if(length(af_cols) > 0){
    
    freq <- melt(
      dt,
      id.vars = "variant_id",
      measure.vars = af_cols,
      variable.name = "population",
      value.name = "af"
    )
    
    # -------------------------
    # LIMPIAR VALORES
    # -------------------------
    freq[, af := as.character(af)]
    freq[af == "-" | af == "", af := NA]
    
    # convertir a numérico
    freq[, af := as.numeric(af)]
    
    DBI::dbWriteTable(con, "variant_population", freq, append=TRUE)
  }
  
  message("✔ Variants procesadas e insertadas en DB")
}