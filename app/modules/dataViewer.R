library(shiny)
library(DT)
library(DBI)
library(dplyr)

# =====================
# UI
# =====================
dataViewerUI <- function(id){
  ns <- NS(id)
  
  div(class="content",
      
      h2("🔎 | DATA VIEWER"),
      
      # FILTRO GLOBAL
      textInput(
        ns("gene"),
        "Filter by gene:",
        placeholder = "e.g. DPM1"
      ),
      
      tabsetPanel(
        
        # =====================
        # VARIANTS
        # =====================
        tabPanel("Variants (WES + WGS)",
                 
                 br(),
                 
                 fluidRow(
                   
                   column(3,
                          sliderInput(ns("af"),
                                      "Max AF:",
                                      min = 0,
                                      max = 0.05,
                                      value = 0.01,
                                      step = 0.001)
                   ),
                   
                   column(3,
                          checkboxGroupInput(
                            ns("impact"),
                            "Impact:",
                            choices = c("HIGH","MODERATE","LOW","MODIFIER"),
                            selected = c("HIGH","MODERATE")
                          )
                   ),
                   
                   column(3,
                          checkboxGroupInput(
                            ns("source"),
                            "Source:",
                            choices = c("WES","WGS","BOTH"),
                            selected = c("WES","WGS","BOTH")
                          )
                   ),
                   
                   column(3,
                          checkboxGroupInput(
                            ns("inheritance"),
                            "Inheritance:",
                            choices = c("de_novo","recessive","dominant","other"),
                            selected = c("de_novo","recessive","dominant")
                          )
                   )
                 ),
                 
                 fluidRow(
                   
                   column(4,
                          selectInput(
                            ns("consequence"),
                            "Consequence:",
                            choices = c("ALL",
                                        "missense_variant",
                                        "stop_gained",
                                        "frameshift_variant",
                                        "splice_region_variant"),
                            selected = "ALL"
                          )
                   ),
                   
                   column(4,
                          selectInput(
                            ns("variant_class"),
                            "Variant class:",
                            choices = c("ALL","SNV","insertion","deletion"),
                            selected = "ALL"
                          )
                   ),
                   
                 ),
                 
                 br(),
                 
                 DTOutput(ns("variants"))
        ),
        
        # =====================
        # RNA
        # =====================
        tabPanel("RNA Data",
                 DTOutput(ns("rna"))
        ),
        
        # =====================
        # DROP
        # =====================
        tabPanel("DROP",
                 
                 selectInput(
                   ns("drop_type"),
                   "Select dataset:",
                   choices = c(
                     "Aberrant Expression"="expr",
                     "Aberrant Splicing"="splicing",
                     "MAE (Monoallelic Expression)"="mae"
                   )
                 ),
                 
                 DTOutput(ns("drop"))
        )
      )
  )
}


# =====================
# SERVER
# =====================
dataViewerServer <- function(id, con){
  moduleServer(id, function(input, output, session){
    
    # helper NULL-safe
    `%||%` <- function(a, b) if (is.null(a)) b else a
    
    # =====================
    # VARIANTS
    # =====================
    output$variants <- renderDT({
      
      req(con)
      
      tryCatch({
        
        df <- get_variants_with_inheritance(con)
        
        validate(
          need(nrow(df) > 0, "No variants loaded")
        )
        
        # =====================
        # FILTROS
        # =====================
        
        # filtro global por gen
        if(nzchar(input$gene)){
          df <- df %>%
            filter(grepl(input$gene, `Gene name`, ignore.case = TRUE))
        }
        
        # AF
        df <- df %>%
          filter(is.na(MAX_AF) | MAX_AF <= input$af)
        
        # impact
        df <- df %>%
          filter(IMPACT %in% (input$impact %||% unique(df$IMPACT)))
        
        # source
        df <- df %>%
          filter(source %in% (input$source %||% unique(df$source)))
        
        # inheritance
        df <- df %>%
          filter(inheritance_type %in% (input$inheritance %||% unique(df$inheritance_type)))
        
        # consequence
        if(!is.null(input$consequence) && input$consequence != "ALL"){
          df <- df %>%
            filter(grepl(input$consequence, Consequence, ignore.case = TRUE))
        }
        
        # variant class
        if(!is.null(input$variant_class) && input$variant_class != "ALL"){
          df <- df %>%
            filter(VARIANT_CLASS == input$variant_class)
        }
        
        validate(
          need(nrow(df) > 0, "No variants match filters")
        )
        
        # =====================
        # OUTPUT
        # =====================
        datatable(
          df,
          rownames = FALSE,
          options = list(
            scrollX = TRUE,
            pageLength = 10,
            autoWidth = TRUE,
            lengthMenu = c(10,25,50,100),
            dom = 'lftip'
          )
        )
        
      }, error=function(e){
        datatable(data.frame(Message="No variants loaded"))
      })
    })
    
    
    # =====================
    # RNA
    # =====================
    output$rna <- renderDT({
      
      req(con)
      
      tryCatch({
        
        df <- get_rna(con)
        
        validate(
          need(nrow(df) > 0, "No RNA data")
        )
        
        # filtro global por gen
        if(nzchar(input$gene)){
          df <- df %>%
            filter(grepl(input$gene, gene_name, ignore.case = TRUE))
        }
        
        validate(
          need(nrow(df) > 0, "No RNA matches gene filter")
        )
        
        datatable(
          df,
          rownames = FALSE,
          options = list(
            scrollX = TRUE,
            pageLength = 10
          )
        )
        
      }, error=function(e){
        datatable(data.frame(Message="No RNA data"))
      })
    })
    
    
    # =====================
    # DROP
    # =====================
    output$drop <- renderDT({
      
      req(con)
      
      tryCatch({
        
        df <- switch(input$drop_type,
                     "expr" = get_drop_expr(con),
                     "splicing" = get_drop_splicing(con),
                     "mae" = get_drop_mae(con)
        )
        
        validate(
          need(!is.null(df) && nrow(df) > 0, "No DROP data")
        )
        
        # filtro global por gen (adaptable)
        if(nzchar(input$gene)){
          
          if("gene_name" %in% colnames(df)){
            df <- df %>% filter(grepl(input$gene, gene_name, ignore.case = TRUE))
            
          } else if("hgncSymbol" %in% colnames(df)){
            df <- df %>% filter(grepl(input$gene, hgncSymbol, ignore.case = TRUE))
          }
        }
        
        validate(
          need(nrow(df) > 0, "No DROP matches gene filter")
        )
        
        datatable(
          df,
          rownames = FALSE,
          options = list(
            scrollX = TRUE,
            pageLength = 10
          )
        )
        
      }, error=function(e){
        datatable(data.frame(Message="No DROP data"))
      })
    })
    
  })
}

