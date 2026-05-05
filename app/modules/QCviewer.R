library(shiny)
library(ggplot2)
library(dplyr)
library(DT)
library(plotly)
library(shinycssloaders)

# =====================
# UI
# =====================
QCviewerUI <- function(id){
  ns <- NS(id)
  
  div(class="content",
      
      h2("📊 | QC VIEWER"),
      
      # ===== PLOTS (EN PARALELO) =====
      div(style="display:flex; gap:20px; flex-wrap:wrap;",
          
          div(class="table-box", style="flex:1; min-width:400px;",
              h4("Mean Coverage", style="color:#8b1e5b;"),
              withSpinner(plotlyOutput(ns("mean_cov_plot")), type = 4, color="#8b1e5b")
          ),
          
          div(class="table-box", style="flex:1; min-width:400px;",
              h4("Coverage ≥20X", style="color:#8b1e5b;"),
              withSpinner(plotlyOutput(ns("cov20_plot")), type = 4, color="#8b1e5b")
          )
      ),
      
      br(),
      
      # ===== TABLA WES (FULL WIDTH) =====
      div(class="table-box",
          h4("WES QC", style="color:#8b1e5b;"),
          DTOutput(ns("wes_table"))
      ),
      
      br(),
      
      # ===== TABLA WGS (FULL WIDTH) =====
      div(class="table-box",
          h4("WGS QC", style="color:#8b1e5b;"),
          DTOutput(ns("wgs_table"))
      )
  )
}

# =====================
# SERVER
# =====================
QCviewerServer <- function(id, pool){
  moduleServer(id, function(input, output, session){
    
    # =====================
    # DATA
    # =====================
    wes <- reactive({
      tryCatch({ get_qc_wes(pool) }, error=function(e){ print(e); NULL })
    })
    
    wgs <- reactive({
      tryCatch({ get_qc_wgs(pool) }, error=function(e){ print(e); NULL })
    })
    
    combined <- reactive({
      tryCatch({ get_qc_combined(pool) }, error=function(e){ print(e); NULL })
    })
    
    
    # =====================
    # MEAN COVERAGE
    # =====================
    output$mean_cov_plot <- renderPlotly({
      
      df <- combined()
      validate(need(!is.null(df), "No QC data loaded"))
      
      p <- ggplot(df, aes(
        x = SAMPLE, 
        y = MEAN_TARGET_COVERAGE, 
        fill = type,
        text = paste0(
          "Sample: ", SAMPLE,
          "<br>Type: ", type,
          "<br>Mean coverage: ", round(MEAN_TARGET_COVERAGE, 2)
        )
      )) +
        geom_bar(stat="identity", position="dodge") +
        geom_hline(yintercept=30, linetype="dashed", color="#8b1e5b") +
        scale_fill_manual(values=c("#2c7fb8","#8b1e5b")) +
        theme_minimal() +
        labs(y="Mean Coverage", x="Sample")
      
      ggplotly(p, tooltip="text") %>%
        layout(legend=list(orientation="h"), margin=list(t=20))
    })
    
    
    # =====================
    # COVERAGE ≥20X
    # =====================
    output$cov20_plot <- renderPlotly({
      
      df <- combined()
      validate(need(!is.null(df), "No QC data loaded"))
      
      p <- ggplot(df, aes(
        x = SAMPLE, 
        y = PCT_TARGET_BASES_20X, 
        fill = type,
        text = paste0(
          "Sample: ", SAMPLE,
          "<br>Type: ", type,
          "<br>% ≥20X: ", round(PCT_TARGET_BASES_20X, 3)
        )
      )) +
        geom_bar(stat="identity", position="dodge") +
        geom_hline(yintercept=0.8, linetype="dashed", color="#8b1e5b") +
        scale_fill_manual(values=c("#2c7fb8","#8b1e5b")) +
        theme_minimal() +
        labs(y="% Bases ≥20X", x="Sample")
      
      ggplotly(p, tooltip="text") %>%
        layout(legend=list(orientation="h"), margin=list(t=20))
    })
    
    
    # =====================
    # TABLE WES
    # =====================
    output$wes_table <- renderDT({
      
      df <- wes()
      
      if(is.null(df)){
        return(datatable(data.frame(Message="No WES QC loaded")))
      }
      
      datatable(
        df,
        rownames = FALSE,
        options = list(
          scrollX = TRUE,
          pageLength = 5
        )
      )
    })
    
    
    # =====================
    # TABLE WGS
    # =====================
    output$wgs_table <- renderDT({
      
      df <- wgs()
      
      if(is.null(df)){
        return(datatable(data.frame(Message="No WGS QC loaded")))
      }
      
      datatable(
        df,
        rownames = FALSE,
        options = list(
          scrollX = TRUE,
          pageLength = 5
        )
      )
    })
    
  })
}