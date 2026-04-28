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
    # GENERIC UPLOAD (TABULAR)
    # -------------------------
    handle_upload <- function(file_input, folder, prefix, type, source){
      
      req(file_input)
      
      tryCatch({
        dt <- fread(file_input$datapath)
        
        path <- save_parquet(dt, folder, prefix)
        id <- paste0(prefix, "_", as.integer(Sys.time()))
        
        add_file(id, file_input$name, type, source, path)
        
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
        tags$pre("CHROM POS REF ALT SYMBOL Consequence"),
        tags$b("Example:"),
        tags$pre("1 123 A G BRCA1 missense_variant"),
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
        tags$pre("hgncSymbol geneID pValue padjust zScore l2fc"),
        tags$b("Example:"),
        tags$pre("BRCA1 ENSG000001 0.001 0.01 -3.5 -2.1"),
        easyClose=TRUE
      ))
    })
    
    observeEvent(input$info_drop_splicing,{
      showModal(modalDialog(
        title="DROP - Aberrant Splicing",
        p("Splicing outliers."),
        tags$b("Required columns:"),
        tags$pre("seqnames start end hgncSymbol pValue padjust deltaPsi"),
        tags$b("Example:"),
        tags$pre("1 12345 12400 BRCA1 0.002 0.02 -0.3"),
        easyClose=TRUE
      ))
    })
    
    observeEvent(input$info_drop_mae,{
      showModal(modalDialog(
        title="DROP - MAE",
        p("Monoallelic expression."),
        tags$b("Required columns:"),
        tags$pre("gene_name refCount altCount pvalue log2FC"),
        tags$b("Example:"),
        tags$pre("BRCA1 100 5 0.0001 -4.3"),
        easyClose=TRUE
      ))
    })
    
    observeEvent(input$info_cov,{
      showModal(modalDialog(
        title = "Coverage & Gene Annotation",
        
        tags$h4("Coverage files (BigWig)"),
        p("Genome-wide coverage files used for visualization in the IGV viewer."),
        
        tags$b("Accepted format:"),
        tags$pre(".bw (BigWig)"),
        
        tags$b("Description:"),
        tags$ul(
          tags$li("Binary, indexed format for fast genomic coverage access"),
          tags$li("Generated from BAM/CRAM files (e.g. using deepTools bamCoverage)"),
          tags$li("Contains coverage signal across the genome (not tabular data)")
        ),
        
        tags$b("Examples:"),
        tags$pre("WGS.bw\nWES.bw\nRNA.bw"),
        
        tags$hr(),
        
        tags$h4("🧬 Gene annotation (GTF)"),
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
        tags$pre("SAMPLE MEAN_TARGET_COVERAGE PCT_TARGET_BASES_10X ..."),
        tags$b("Example:"),
        tags$pre("sample1 120 0.98"),
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