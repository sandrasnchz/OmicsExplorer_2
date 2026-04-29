library(shiny)
library(DT)
library(dplyr)
library(data.table)
library(arrow)

# =====================
# UI
# =====================
loadUI <- function(id){
  ns <- NS(id)
  
  upload_item <- function(label, inputId, infoId){
    div(class="upload-item",
        
        div(style="display:flex; align-items:center; gap:8px;",
            tags$label(label),
            
            actionButton(
              ns(infoId),
              label = NULL,
              icon = icon("info-circle"),
              class = "btn btn-link",
              style = "padding:0;"
            )
        ),
        
        fileInput(ns(inputId), label = NULL)
    )
  }
  
  div(class="content",
      
      h2("📥 | LOAD DATA"),
      
      div(class="upload-box",
          
          div(class="section-title-big", "Upload data"),
          
          div(style="display:flex; gap:30px; justify-content:center; flex-wrap:wrap;",
              
              div(class="upload-col",
                  h4("Variants"),
                  upload_item("WES variants", "wes_file", "info_wes"),
                  upload_item("WGS variants", "wgs_file", "info_wgs")
              ),
              
              div(class="upload-col",
                  h4("RNA"),
                  upload_item("Gene TPM (sample)", "rna_tpm", "info_rna_tpm"),
                  upload_item("RNA controls TPM", "rna_controls", "info_rna_controls")
              ),
              
              div(class="upload-col",
                  h4("DROP"),
                  upload_item("Aberrant Expression", "drop_expr", "info_drop_expr"),
                  upload_item("Aberrant Splicing", "drop_splicing", "info_drop_splicing"),
                  upload_item("MAE", "drop_mae", "info_drop_mae")
              ),
              
              div(class="upload-col",
                  h4("Coverage"),
                  upload_item("WES coverage (.bw)", "cov_wes", "info_cov"),
                  upload_item("WGS coverage (.bw)", "cov_wgs", "info_cov"),
                  upload_item("RNA coverage (.bw)", "cov_rna", "info_cov"),
                  upload_item("Gene annotation (.gtf)", "cov_gtf", "info_cov")
              ),
              
              div(class="upload-col",
                  h4("Sample QC"),
                  upload_item("WES QC", "qc_wes", "info_qc"),
                  upload_item("WGS QC", "qc_wgs", "info_qc")
              )
          )
      ),
      
      div(class="section-center",
          div(class="section-title-big", "Loaded datasets")
      ),
      
      div(class="table-box",
          DTOutput(ns("table"))
      )
  )
}

# =====================
# SERVER
# =====================
loadServer <- function(id, con){
  moduleServer(id, function(input, output, session){
    
    # -------------------------
    # SAVE PARQUET
    # -------------------------
    save_parquet <- function(dt, folder, name){
      dir.create(folder, recursive = TRUE, showWarnings = FALSE)
      safe_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
      path <- file.path(folder, paste0(name, "_", safe_time, ".parquet"))
      write_parquet(dt, path)
      return(path)
    }
    
    # -------------------------
    # STORAGE
    # -------------------------
    rv <- reactiveValues(
      files_info = data.frame(
        id=character(),
        File=character(),
        Type=character(),
        Source=character(),
        path=character(),
        stringsAsFactors=FALSE
      )
    )
    
    add_file <- function(id, file, type, source, path){
      rv$files_info <- bind_rows(
        rv$files_info,
        data.frame(
          id=id,
          File=file,
          Type=type,
          Source=source,
          path=path,
          stringsAsFactors=FALSE
        )
      )
    }
    
    # -------------------------
    # INFO FILES
    # -------------------------
    get_file_metadata <- function(filepath){
      
      fname <- basename(filepath)
      fname_lower <- tolower(fname)
      
      # =========================
      # TYPE (orden IMPORTANTE)
      # =========================
      
      if (grepl("qc", fname_lower)) {
        type <- "QC"
        
      } else if (grepl("^sample_", fname_lower)) {
        type <- "RNA"
        
      } else if (grepl("controls", fname_lower)) {
        type <- "RNA"
        
      } else if (grepl("expression|splicing|mae", fname_lower)) {
        type <- "DROP"
        
      } else if (grepl("\\.bw|\\.gtf", fname_lower)) {
        type <- "Coverage"
        
      } else if (grepl("^wes|^wgs", fname_lower)) {
        type <- "Variants"
        
      } else {
        type <- "Unknown"
      }
      
      # =========================
      # SOURCE 
      # =========================
      
      if (grepl("qc", fname_lower)) {
        
        if (grepl("wes", fname_lower)) {
          source <- "WES"
        } else if (grepl("wgs", fname_lower)) {
          source <- "WGS"
        } else {
          source <- "QC"
        }
        
      } else if (grepl("^sample_", fname_lower)) {
        
        # extraer nombre de muestra 
        sample_name <- sub("^sample_([^_]+).*", "\\1", fname_lower)
        source <- sample_name
        
      } else if (startsWith(fname_lower, "controls")) {
        source <- "Controls"
        
      } else if (startsWith(fname_lower, "wes")) {
        source <- "WES"
        
      } else if (startsWith(fname_lower, "wgs")) {
        source <- "WGS"
        
      } else if (startsWith(fname_lower, "rna")) {
        source <- "RNA"
        
      } else if (startsWith(fname_lower, "expression")) {
        source <- "Expression"
        
      } else if (startsWith(fname_lower, "splicing")) {
        source <- "Splicing"
        
      } else if (startsWith(fname_lower, "mae")) {
        source <- "MAE"
        
      } else {
        source <- "Unknown"
      }
      
      return(list(type = type, source = source))
    }
    
    # ==========================
    # LOAD EXISTING FILES
    # ==========================
    load_existing_files <- function(){
      
      folders <- c(
        "../data/variants",
        "../data/rnaseq",
        "../data/drop",
        "../data/coverage",
        "../data/qc"
      )
      
      all_files <- list()
      
      for(folder in folders){
        
        if(!dir.exists(folder)) next
        
        files <- list.files(folder, full.names = TRUE)
        
        if(length(files) == 0) next
        
        df <- lapply(files, function(f){
          
          meta <- get_file_metadata(f)
          
          data.frame(
            id = paste0("existing_", digest::digest(f)),  # ID estable
            File = basename(f),
            Type = meta$type,
            Source = meta$source,
            path = f,
            stringsAsFactors = FALSE
          )
        })
        
        all_files[[length(all_files)+1]] <- bind_rows(df)
      }
      
      if(length(all_files) > 0){
        rv$files_info <- bind_rows(all_files)
      }
    }
    
    # Ejecutar al iniciar
    load_existing_files()
    
    # -------------------------
    # GENERIC UPLOAD (TABULAR)
    # -------------------------
    handle_upload <- function(file_input, folder, prefix, type, source){
      
      req(file_input)
      
      tryCatch({
        dt <- fread(file_input$datapath)
        
        path <- save_parquet(dt, folder, prefix)
        id <- paste0(prefix, "_", as.integer(Sys.time()))
        
        add_file(id, file_input$name, type, source, path)
        
        showNotification(
          paste("Loaded:", file_input$name),
          type = "message",
          duration = 3
        )
        
      }, error=function(e){
        showNotification(paste("Error:", e$message), type="error")
      })
    }
    
    # -------------------------
    # COVERAGE UPLOAD (BW / GTF)
    # -------------------------
    handle_bw_upload <- function(file_input, folder, prefix, type, source){
      
      req(file_input)
      
      tryCatch({

        folder <- normalizePath(folder, mustWork = FALSE)
        
        if(!dir.exists(folder)){
          dir.create(folder, recursive = TRUE)
        }
        
        safe_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
        ext <- tools::file_ext(file_input$name)
        
        new_name <- paste0(prefix, "_", safe_time, ".", ext)
        
        dest_path <- file.path(folder, new_name)

        success <- file.copy(file_input$datapath, dest_path, overwrite = TRUE)
        
        print(paste("Copy success:", success))
        print(paste("Saved at:", dest_path))
        
        if(!success){
          stop("File copy failed")
        }
        
        # guardar en tabla reactiva
        id <- paste0(prefix, "_", as.integer(Sys.time()))
        
        add_file(id, file_input$name, type, source, dest_path)
        
        showNotification("File uploaded correctly", type = "message")
        
      }, error=function(e){
        showNotification(paste("Error:", e$message), type="error")
        print(e)
      })
    }
    
    # =====================
    # INFO MODALS
    # =====================
    
    observeEvent(input$info_wes,{
      showModal(modalDialog(
        title="WES Variants",
        p("VEP annotated variants file."),
        tags$b("Required columns:"),
        tags$pre("CHROM POS REF ALT SYMBOL Consequence"),
        tags$b("Example:"),
        tags$pre("1 123 A G BRCA1 missense_variant"),
        easyClose=TRUE
      ))
    })
    
    observeEvent(input$info_wgs,{
      showModal(modalDialog(
        title="WGS Variants",
        p("VEP annotated variants file."),
        tags$b("Required columns:"),
        tags$pre("ID	CHROM	POS	REF	ALT	FILTER	PARENT1_GT	PARENT1_DP	PARENT1_AD	PARENT1_GQ	PARENT2_GT	PARENT2_DP	PARENT2_AD	PARENT2_GQ	CHILD_GT	CHILD_DP	CHILD_AD	CHILD_GQ	Gene	Location	Allele	Feature ... "),
        tags$b("Example:"),
        tags$pre("1_100148710_G_C	1	100148710	G	C	PASS	0/0	21	21,0	51	0/0	19	19,0	51	0/1	25	11,14	99	ENSG00000099260	1:100148711-100148721	-	ENST00000496843	Transcript ... "),
        easyClose=TRUE
      ))
    })
    
    
    observeEvent(input$info_rna_tpm,{
      showModal(modalDialog(
        title="RNA TPM (sample)",
        p("Gene expression values normalized as TPM."),
        tags$b("Required columns:"),
        tags$pre("gene_id gene_name gene_tpm"),
        tags$b("Example:"),
        tags$pre("ENSG000001 BRCA1 10.5"),
        easyClose=TRUE
      ))
    })
    
    observeEvent(input$info_rna_controls,{
      showModal(modalDialog(
        title="RNA Controls TPM",
        p("Summary statistics from control cohort."),
        tags$b("Required columns:"),
        tags$pre("gene_id gene_name max_TPM min_TPM mean_TPM  median_TPM"),
        tags$b("Example:"),
        tags$pre("ENSG000001 BRCA1 20.5 5.2 10.3 9.8"),
        easyClose=TRUE
      ))
    })
    
    observeEvent(input$info_drop_expr,{
      showModal(modalDialog(
        title="DROP - Aberrant Expression",
        p("Gene expression outliers."),
        tags$b("Required columns:"),
        tags$pre("hgncSymbol	geneID	sampleID	pValue	padjust	zScore	l2fc	rawcounts	meanRawcounts	normcounts	meanCorrected	theta	aberrant	AberrantBySample	AberrantByGene	padj_rank	foldChange"),
        tags$b("Example:"),
        tags$pre("RPL3P4	ENSG00000232573	ND0013	5.47792538302315E-12	8.59402812907001E-07	-8.5	-1.39	153	3809.48	1146.46	3015.01	92.13	TRUE	3	1	1	0.38"),
        easyClose=TRUE
      ))
    })
    
    observeEvent(input$info_drop_splicing,{
      showModal(modalDialog(
        title="DROP - Aberrant Splicing",
        p("Splicing outliers."),
        tags$b("Required columns:"),
        tags$pre("seqnames	start	end	width	strand	sampleID	hgncSymbol	type	pValue	padjust	psiValue	deltaPsi	counts	totalCounts	meanCounts	meanTotalCounts	nonsplitCounts	nonsplitProportion	nonsplitProportion_99quantile	annotatedJunction	pValueGene	padjustGene	potentialImpact	causesFrameshift	UTR_overlap	blacklist"),
        tags$b("Example:"),
        tags$pre("IGHG2	jaccard	1.2962E-15	3.2391E-11	0.6	-0.38	58	96	1485.97	1504.89	3	0.03	0.11	both	3.8886E-15	1.11E-11	annotatedIntron_reducedUsage	unlikely	no	TRUE"),
        easyClose=TRUE
      ))
    })
    
    observeEvent(input$info_drop_mae,{
      showModal(modalDialog(
        title="DROP - MAE",
        p("Monoallelic expression."),
        tags$b("Required columns:"),
        tags$pre("gene_name	ID	contig	position	variantID	refAllele	altAllele	refCount	altCount	totalCount	pvalue	padj	log2FC	altRatio	AF	AF_afr	AF_amr	AF_eas	AF_nfe	MAX_AF	rare	gene_type	other_names	N_var	cohort_freq	MAE	MAE_ALT"),
        tags$b("Example:"),
        tags$pre("HLA-A	ND0013_ND0013--ND0013	chr6	29912108	rs2231095	G	C	9	3994	4003	2.49631997819451E-41	6.74440536717421E-39	8.93119194898805	0.998	0.47	0.48	0.44	0.44	0.44	0.48	FALSE	protein_coding		1	0.25	TRUE	TRUE"),
        easyClose=TRUE
      ))
    })
    
    observeEvent(input$info_cov,{
      showModal(modalDialog(
        title = "Coverage",
        
        tags$h4("Coverage files (BigWig)"),
        p("Genome-wide coverage files used for visualization in the IGV viewer."),
        
        tags$b("Accepted format:"),
        tags$pre(".bw (BigWig)"),
        
        tags$b("Description:"),
        tags$ul(
          tags$li("Binary, indexed format for fast genomic coverage access"),
          tags$li("Generated from BAM/CRAM files"),
          tags$li("Contains coverage signal across the genome (not tabular data)")
        ),
        
        tags$b("Examples:"),
        tags$pre("WGS.bw\nWES.bw\nRNA.bw"),
        
        tags$hr(),
        
        tags$h4("Gene annotation (GTF)"),
        p("Annotation file used to display genes, transcripts and exons in the viewer."),
        
        tags$b("Accepted format:"),
        tags$pre(".gtf"),
        
        tags$b("Recommended:"),
        tags$ul(
          tags$li("Comprehensive gene annotation (CHR)"),
          tags$li("Matching genome build (e.g. GRCh37 / hg19)")
        ),
        
        tags$b("Description:"),
        tags$ul(
          tags$li("Text-based file describing gene structures"),
          tags$li("Includes gene coordinates, exons, transcripts, etc."),
          tags$li("Used for gene visualization and navigation")
        ),
        
        tags$b("Example:"),
        tags$pre("chr1\tensembl\tgene\t11869\t14409\t.\t+\t.\tgene_name \"DDX11L1\";"),
        
        tags$hr(),
        
        tags$h4("Important notes"),
        tags$ul(
          tags$li("Coverage files (.bw) are NOT tabular and cannot be read with fread"),
          tags$li("They are used directly in the genome browser (IGV.js)"),
          tags$li("Make sure all files use the same genome reference (e.g. GRCh37)")
        ),
        
        easyClose = TRUE
      ))
    })
    
    observeEvent(input$info_qc,{
      showModal(modalDialog(
        title="Sample QC",
        p("Sequencing quality metrics."),
        tags$b("Required columns:"),
        tags$pre("SAMPLE	MEAN_TARGET_COVERAGE	PCT_USABLE_BASES_ON_TARGET	FOLD_ENRICHMENT	PCT_TARGET_BASES_10X	PCT_TARGET_BASES_20X	PCT_TARGET_BASES_30X	PCT_TARGET_BASES_40X	PCT_TARGET_BASES_50X"),
        tags$b("Example:"),
        tags$pre("ND0013	47.588148	0.316725	61.775056	0.903266	0.790783	0.655728	0.521333	0.401723"),
        easyClose=TRUE
      ))
    })
    
    
    # =====================
    # VARIANTS
    # =====================
    observeEvent(input$wes_file,{
      handle_upload(input$wes_file, "../data/variants", "wes", "Variants", "WES")
    })
    
    observeEvent(input$wgs_file,{
      handle_upload(input$wgs_file, "../data/variants", "wgs", "Variants", "WGS")
    })
    
    # =====================
    # RNA TPM (ESPECIAL)
    # =====================
    observeEvent(input$rna_tpm,{
      req(input$rna_tpm)
      
      tryCatch({
        
        dt <- fread(input$rna_tpm$datapath)
        
        sample_name <- colnames(dt)[3]
        
        dt_clean <- dt %>%
          rename(tpm = all_of(sample_name)) %>%
          mutate(sample = sample_name)
        
        path <- save_parquet(dt_clean, "../data/rnaseq", paste0("sample_", sample_name))
        
        id <- paste0("rna_", as.integer(Sys.time()))
        
        add_file(id, input$rna_tpm$name, "RNA", sample_name, path)
        
        showNotification(
          paste("Loaded:", input$rna_tpm$name),
          type = "message",
          duration = 3
        )
        
      }, error=function(e){
        showNotification(paste("Error RNA:", e$message), type="error")
      })
    })
    
    # =====================
    # RESTO (GENÉRICO)
    # =====================
    observeEvent(input$rna_controls,{
      handle_upload(input$rna_controls, "../data/rnaseq", "controls", "RNA", "Controls")
    })
    
    observeEvent(input$drop_expr,{
      handle_upload(input$drop_expr, "../data/drop", "expression", "DROP", "Expression")
    })
    
    observeEvent(input$drop_splicing,{
      handle_upload(input$drop_splicing, "../data/drop", "splicing", "DROP", "Splicing")
    })
    
    observeEvent(input$drop_mae,{
      handle_upload(input$drop_mae, "../data/drop", "mae", "DROP", "MAE")
    })
    
    observeEvent(input$qc_wes,{
      handle_upload(input$qc_wes, "../data/qc", "wes_qc", "QC", "WES")
    })
    
    observeEvent(input$qc_wgs,{
      handle_upload(input$qc_wgs, "../data/qc", "wgs_qc", "QC", "WGS")
    })
    
    # =====================
    # COVERAGE EVENTS
    # =====================
    observeEvent(input$cov_wes,{
      handle_bw_upload(input$cov_wes, "../data/coverage", "WES", "Coverage", "WES")
    })
    
    observeEvent(input$cov_wgs,{
      handle_bw_upload(input$cov_wgs, "../data/coverage", "WGS", "Coverage", "WGS")
    })
    
    observeEvent(input$cov_rna,{
      handle_bw_upload(input$cov_rna, "../data/coverage", "RNA", "Coverage", "RNA")
    })
    
    observeEvent(input$cov_gtf,{
      handle_bw_upload(input$cov_gtf, "../data/coverage", "genes", "Coverage", "GTF")
    })
    
    # =====================
    # DELETE
    # =====================
    observeEvent(input$delete_row,{
      
      row <- rv$files_info %>% filter(id == input$delete_row)
      
      if(nrow(row)==1 && file.exists(row$path)){
        file.remove(row$path)
      }
      
      rv$files_info <- rv$files_info %>%
        filter(id != input$delete_row)
    })
    
    # =====================
    # TABLE
    # =====================
    output$table <- renderDT({
      
      if(nrow(rv$files_info)==0){
        return(datatable(data.frame(Message="No datasets loaded")))
      }
      
      df <- rv$files_info
      
      df$Action <- paste0(
        '<button class="btn-delete" data-id="', df$id, '">🗑️</button>'
      )
      
      datatable(
        df[, c("File","Type","Source","Action")],
        escape = FALSE,
        selection = "none",
        options = list(pageLength = 5),
        
        callback = JS(sprintf("
          table.on('click', '.btn-delete', function() {
            var id = $(this).data('id');
            Shiny.setInputValue('%s', id, {priority: 'event'});
          });
        ", session$ns("delete_row")))
      )
    })
  })
}