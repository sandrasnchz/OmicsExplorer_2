library(shiny)
library(ggplot2)
library(dplyr)
library(DT)
library(tidyr)

# =====================
# UI
# =====================
QCviewerUI <- function(id){
  ns <- NS(id)
  
  div(class="content",
      
      h2("📊 | QC VIEWER"),
      
      div(class="box",
          h3("Mean Coverage"),
          plotOutput(ns("mean_cov_plot"))
      ),
      
      div(class="box",
          h3("Coverage ≥20X"),
          plotOutput(ns("cov20_plot"))
      ),
      
      div(class="box",
          h3("WES QC"),
          DTOutput(ns("wes_table"))
      ),
      
      div(class="box",
          h3("WGS QC"),
          DTOutput(ns("wgs_table"))
      )
  )
}

# =====================
# SERVER
# =====================
QCviewerServer <- function(id, con){
  moduleServer(id, function(input, output, session){
    
    # =====================
    # DATA
    # =====================
    wes <- reactive({
      tryCatch({
        get_qc_wes(con)
      }, error=function(e){
        NULL
      })
    })
    
    wgs <- reactive({
      tryCatch({
        get_qc_wgs(con)
      }, error=function(e){
        NULL
      })
    })
    
    combined <- reactive({
      
      tryCatch({
        get_qc_combined(con)
      }, error=function(e){
        NULL
      })
    })
    
    
    # =====================
    # MEAN COVERAGE
    # =====================
    output$mean_cov_plot <- renderPlot({
      
      df <- combined()
      
      validate(
        need(!is.null(df), "No QC data loaded")
      )
      
      ggplot(df, aes(x = SAMPLE, y = MEAN_TARGET_COVERAGE, fill = type)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_hline(yintercept = 30, linetype = "dashed", color = "red") +
        theme_minimal() +
        labs(y = "Mean Coverage", x = "Sample")
    })
    
    
    # =====================
    # COVERAGE ≥20X
    # =====================
    output$cov20_plot <- renderPlot({
      
      df <- combined()
      
      validate(
        need(!is.null(df), "No QC data loaded")
      )
      
      ggplot(df, aes(x = SAMPLE, y = PCT_TARGET_BASES_20X, fill = type)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
        theme_minimal() +
        labs(y = "% Bases ≥20X", x = "Sample")
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
        options = list(scrollX=TRUE, pageLength=5)
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
        options = list(scrollX=TRUE, pageLength=5)
      )
    })
    
  })
}