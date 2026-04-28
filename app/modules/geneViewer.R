library(shiny)
library(DT)
library(dplyr)

# =====================
# UI
# =====================
geneViewerUI <- function(id){
  ns <- NS(id)
  
  div(class="content",
      
      h2("🧬 | GENE VIEWER"),
      
      textInput(
        ns("gene"),
        "Filter by gene:",
        placeholder = "e.g. DPM1"
      ),
      
      br(),
      
      uiOutput(ns("gene_summary")),
      
      br(),
      
      uiOutput(ns("gene_info")),
      
      br(),
      
      uiOutput(ns("external_links")),
      
      br(),
      
      uiOutput(ns("ucsc_button"))
  )
}

# =====================
# SERVER
# =====================
geneViewerServer <- function(id, con, selected_gene){
  moduleServer(id, function(input, output, session){
    
    observe({
      if(!is.null(selected_gene())){
        updateTextInput(session, "gene", value = selected_gene())
      }
    })
    
    # =====================
    # 🔁 SINCRONIZACIÓN INPUT
    # =====================
    observe({
      gene <- selected_gene()
      
      if(!is.null(gene) && gene != ""){
        updateTextInput(session, "gene", value = gene)
      }
    })
    
    
    # =====================
    # 🎯 GENE ACTUAL (CLAVE)
    # =====================
    current_gene <- reactive({
      
      # prioridad: navegación desde DataViewer
      if(!is.null(selected_gene()) && selected_gene() != ""){
        return(selected_gene())
      }
      
      # fallback: input manual
      if(!is.null(input$gene) && input$gene != ""){
        return(input$gene)
      }
      
      return(NULL)
    })
    
    
    # =====================
    # VARIANTS POR GEN (SQL)
    # =====================
    gene_variants <- reactive({
      
      gene <- current_gene()
      req(gene)
      
      tryCatch({
        get_variants_by_gene(con, gene)
      }, error = function(e){
        print(e)
        return(NULL)
      })
    })
    
    
    # =====================
    # INFO GEN (SQL)
    # =====================
    gene_info_filtered <- reactive({
      
      gene <- current_gene()
      req(gene)
      
      tryCatch({
        get_gene_info_by_gene(con, gene)
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
        tags$h3(gene),
        tags$p(paste("Gene ID:", gene_id)),
        tags$p(paste("Number of variants:", n_var))
      )
    })
    
    
    # =====================
    # GENE INFO
    # =====================
    output$gene_info <- renderUI({
      
      df <- gene_info_filtered()
      
      if(is.null(df) || nrow(df) == 0){
        return(NULL)
      }
      
      row <- df[1, ]
      
      tags$div(
        class = "box",
        tags$h4("Gene information"),
        
        tags$table(
          class = "table table-sm",
          
          tags$tr(tags$th("HGNC ID"), tags$td(row$HGNC_ID)),
          tags$tr(tags$th("Biotype"), tags$td(row$BIOTYPE)),
          tags$tr(tags$th("Gene phenotype"), tags$td(row$GENE_PHENO)),
          tags$tr(tags$th("Function"), tags$td(row$Function_description)),
          tags$tr(tags$th("Disease"), tags$td(row$Disease_description)),
          tags$tr(tags$th("HPO ID"), tags$td(row$HPO_id)),
          tags$tr(tags$th("HPO name"), tags$td(row$HPO_name))
        )
      )
    })
    
    
    # =====================
    # LINKS EXTERNOS
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
      
      tagList(
        tags$a("🧬 GeneCards",
               href = genecards_url,
               target="_blank",
               class="btn btn-outline-primary"),
        
        if(!is.null(omim) && !is.na(omim)){
          tags$a("📖 OMIM",
                 href = paste0("https://www.omim.org/entry/", omim),
                 target="_blank",
                 class="btn btn-outline-secondary")
        }
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
      
      tags$a(
        "🧬 View in UCSC Genome Browser",
        href = paste0(
          "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=",
          chr, ":", start-1000, "-", end+1000
        ),
        target = "_blank",
        class = "btn btn-primary"
      )
    })
    
  })
}