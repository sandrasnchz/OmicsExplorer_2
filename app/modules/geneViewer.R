library(shiny)
library(shinycssloaders)
library(DT)
library(dplyr)

# =====================
# UI
# =====================
geneViewerUI <- function(id){
  ns <- NS(id)
  
  div(class="content",
      
      h2("🧬 | GENE VIEWER"),
      
      # ===== FILTER =====
      div(class="filter-box",
          textInput(
            ns("gene"),
            "Filter by gene:",
            placeholder = "e.g. DPM1"
          )
      ),
      
      # ===== SUMMARY =====
      div(class="table-box",
          withSpinner(
            uiOutput(ns("gene_summary")),
            type = 4,
            color = "#8b1e5b"
          )
      ),
      
      br(),
      
      # ===== INFO =====
      div(class="table-box",
          withSpinner(
            uiOutput(ns("gene_info")),
            type = 4,
            color = "#8b1e5b"
          )
      ),
      
      br(),
      
      # ===== LINKS + UCSC =====
      div(style="display:flex; gap:20px; align-items:stretch; flex-wrap:wrap;",
          
          div(class="table-box", style="flex:1; min-width:300px;",
              withSpinner(
                uiOutput(ns("external_links")),
                type = 4,
                color = "#8b1e5b"
              )
          ),
          
          div(class="table-box", style="flex:1; min-width:300px;",
              withSpinner(
                uiOutput(ns("ucsc_button")),
                type = 4,
                color = "#8b1e5b"
              )
          )
      )
  )
}

# =====================
# SERVER
# =====================
geneViewerServer <- function(id, pool, selected_gene){
  moduleServer(id, function(input, output, session){
    
    # =====================
    # SYNC INPUT
    # =====================
    observe({
      if(!is.null(selected_gene()) && selected_gene() != ""){
        updateTextInput(session, "gene", value = selected_gene())
      }
    })
    
    # =====================
    # CURRENT GENE
    # =====================
    current_gene <- reactive({
      if(!is.null(selected_gene()) && selected_gene() != ""){
        return(selected_gene())
      }
      if(!is.null(input$gene) && input$gene != ""){
        return(input$gene)
      }
      return(NULL)
    })
    
    # =====================
    # VARIANTS
    # =====================
    gene_variants <- reactive({
      gene <- current_gene()
      req(gene)
      
      tryCatch({
        get_variants_by_gene(pool, gene)
      }, error = function(e){
        print(e)
        return(NULL)
      })
    })
    
    # =====================
    # GENE INFO
    # =====================
    gene_info_filtered <- reactive({
      gene <- current_gene()
      req(gene)
      
      tryCatch({
        get_gene_info_by_gene(pool, gene)
      }, error = function(e){
        print(e)
        return(NULL)
      })
    })
    
    # =====================
    # SUMMARY
    # =====================
    output$gene_summary <- renderUI({
      
      df <- gene_variants()
      
      if(is.null(df) || nrow(df) == 0){
        return(tags$p("No gene found"))
      }
      
      gene <- unique(df$`Gene name`)[1]
      gene_id <- unique(df$`Gene ID`)[1]
      n_var <- nrow(df)
      
      tags$div(
        style="display:flex; justify-content:space-between; align-items:center;",
        
        tags$div(
          tags$h3(gene, style="color:#8b1e5b; margin-bottom:5px;"),
          tags$p(paste("Gene ID:", gene_id), style="margin:0; color:#555;")
        ),
        
        tags$div(
          style="text-align:right;",
          tags$div("Variants", style="font-size:12px; color:#777;"),
          tags$div(n_var, style="font-size:22px; font-weight:700; color:#8b1e5b;")
        )
      )
    })
    
    # =====================
    # GENE INFO (FIX OVERFLOW)
    # =====================
    output$gene_info <- renderUI({
      
      df <- gene_info_filtered()
      
      if(is.null(df) || nrow(df) == 0){
        return(NULL)
      }
      
      row <- df[1, ]
      
      tags$div(
        tags$h4("Gene information", style="color:#8b1e5b; margin-bottom:15px;"),
        
        # 🔥 CLAVE: evitar que se salga
        tags$div(
          style="overflow-x:auto;",
          
          tags$table(
            class = "table table-sm",
            
            tags$tr(tags$th("HGNC ID"), tags$td(row$HGNC_ID)),
            tags$tr(tags$th("Biotype"), tags$td(row$BIOTYPE)),
            tags$tr(tags$th("Gene phenotype"), tags$td(row$GENE_PHENO)),
            
            tags$tr(tags$th("Function"),
                    tags$td(style="max-width:600px; white-space:normal; word-break:break-word;",
                            row$Function_description)),
            
            tags$tr(tags$th("Disease"),
                    tags$td(style="max-width:600px; white-space:normal; word-break:break-word;",
                            row$Disease_description)),
            
            tags$tr(tags$th("HPO ID"), tags$td(row$HPO_id)),
            
            tags$tr(tags$th("HPO name"),
                    tags$td(style="max-width:600px; white-space:normal; word-break:break-word;",
                            row$HPO_name))
          )
        )
      )
    })
    
    # =====================
    # EXTERNAL LINKS
    # =====================
    output$external_links <- renderUI({
      
      df_info <- gene_info_filtered()
      df_var  <- gene_variants()
      
      if(is.null(df_var) || nrow(df_var) == 0) return(NULL)
      
      gene <- unique(df_var$`Gene name`)[1]
      
      genecards_url <- paste0(
        "https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
        gene
      )
      
      omim <- NULL
      if(!is.null(df_info) && "OMIM_id" %in% colnames(df_info)){
        omim <- unique(df_info$OMIM_id)[1]
      }
      
      tags$div(
        tags$h4("External resources", style="color:#8b1e5b; margin-bottom:10px;"),
        
        tags$div(style="display:flex; gap:10px; flex-wrap:wrap;",
                 
                 tags$a("GeneCards",
                        href = genecards_url,
                        target="_blank",
                        class="btn-download"),
                 
                 if(!is.null(omim) && !is.na(omim)){
                   tags$a("OMIM",
                          href = paste0("https://www.omim.org/entry/", omim),
                          target = "_blank",
                          class="btn-download")
                 }
        )
      )
    })
    
    # =====================
    # UCSC
    # =====================
    output$ucsc_button <- renderUI({
      
      df <- gene_variants()
      
      if(is.null(df) || nrow(df) == 0) return(NULL)
      
      chr <- unique(df$CHROM)[1]
      start <- min(df$POS, na.rm = TRUE)
      end <- max(df$POS, na.rm = TRUE)
      
      tags$div(
        tags$h4("Genome browser", style="color:#8b1e5b; margin-bottom:10px;"),
        
        tags$a(
          "Open in UCSC Genome Browser",
          href = paste0(
            "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=",
            chr, ":", start-1000, "-", end+1000
          ),
          target = "_blank",
          class = "btn-download"
        )
      )
    })
    
  })
}