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
      
      textInput(
        ns("gene"),
        "Filter by gene:",
        placeholder = "e.g. DPM1"
      ),
      
      tabsetPanel(id = ns("main_tabs"),
                  
                  tabPanel("Variants (WES + WGS)",
                           
                           br(),
                           
                           fluidRow(
                             column(3,
                                    sliderInput(ns("af"), "Max AF:",
                                                min = 0, max = 0.05,
                                                value = 0.01, step = 0.001)
                             ),
                             column(3,
                                    checkboxGroupInput(ns("impact"), "Impact:",
                                                       choices = c("HIGH","MODERATE","LOW","MODIFIER"),
                                                       selected = c("HIGH","MODERATE","LOW","MODIFIER"))
                             ),
                             column(3,
                                    checkboxGroupInput(ns("source"), "Source:",
                                                       choices = c("WES","WGS","BOTH"),
                                                       selected = c("WES","WGS","BOTH"))
                             ),
                             column(3,
                                    checkboxGroupInput(ns("inheritance"), "Inheritance:",
                                                       choices = c("de_novo","recessive","dominant","other"),
                                                       selected = c("de_novo","recessive","dominant", "other"))
                             )
                           ),
                           
                           fluidRow(
                             column(4,
                                    selectInput(ns("consequence"), "Consequence:",
                                                choices = c("ALL",
                                                            "missense_variant",
                                                            "stop_gained",
                                                            "frameshift_variant",
                                                            "splice_region_variant"),
                                                selected = "ALL")
                             ),
                             column(4,
                                    selectInput(ns("variant_class"), "Variant class:",
                                                choices = c("ALL","SNV","insertion","deletion"),
                                                selected = "ALL")
                             )
                           ),
                           
                           br(),
                           DTOutput(ns("variants"))
                  ),
                  
                  tabPanel("RNA Data",
                           DTOutput(ns("rna"))
                  ),
                  
                  tabPanel("DROP",
                           selectInput(ns("drop_type"),
                                       "Select dataset:",
                                       choices = c(
                                         "Aberrant Expression"="expr",
                                         "Aberrant Splicing"="splicing",
                                         "MAE (Monoallelic Expression)"="mae"
                                       )),
                           DTOutput(ns("drop"))
                  )
      )
  )
}

# =====================
# SERVER
# =====================
dataViewerServer <- function(id, con, selected_gene){
  moduleServer(id, function(input, output, session){
    
    `%||%` <- function(a, b) if (is.null(a)) b else a
    
    navigating <- reactiveVal(FALSE)
    
    # =====================
    # RESET CONTROLADO
    # =====================
    observeEvent(input$main_tabs, {
      if(!navigating()){
        selected_gene(NULL)
      }
      navigating(FALSE)
    })
    
    # =====================
    # VARIANTS
    # =====================
    output$variants <- renderDT({
      
      req(con)
      
      tryCatch({
        
        df <- get_variants_with_inheritance(con)
        validate(need(nrow(df) > 0, "No variants loaded"))
        
        # navegación
        if(!is.null(selected_gene())){
          df <- df %>% filter(grepl(selected_gene(), `Gene name`, ignore.case = TRUE))
        }
        
        # filtro manual
        if(nzchar(input$gene)){
          df <- df %>% filter(grepl(input$gene, `Gene name`, ignore.case = TRUE))
        }
        
        df <- df %>% filter(is.na(MAX_AF) | MAX_AF <= input$af)
        df <- df %>% filter(IMPACT %in% (input$impact %||% unique(df$IMPACT)))
        df <- df %>% filter(source %in% (input$source %||% unique(df$source)))
        df <- df %>% filter(inheritance_type %in% (input$inheritance %||% unique(df$inheritance_type)))
        
        if(input$consequence != "ALL"){
          df <- df %>% filter(grepl(input$consequence, Consequence, ignore.case = TRUE))
        }
        
        if(input$variant_class != "ALL"){
          df <- df %>% filter(VARIANT_CLASS == input$variant_class)
        }
        
        validate(need(nrow(df) > 0, "No variants match filters"))
        
        df$gene_hidden <- df$`Gene name`
        
        datatable(
          df,
          rownames = FALSE,
          escape = FALSE,
          selection = "none",
          options = list(
            scrollX = TRUE,
            pageLength = 10,
            dom = 'lftip',
            columnDefs = list(
              list(targets = ncol(df)-1, visible = FALSE)
            )
          ),
          
          callback = JS(sprintf("
            
            var format = function(gene) {
              return '<div style=\"padding:10px\">' +
                '<b>Explore:</b><br>' +
                '<button class=\"go-gene\" data-gene=\"'+gene+'\">🧬 Gene Viewer</button><br><br>' +
                '<button class=\"go-rna\" data-gene=\"'+gene+'\">RNA</button> ' +
                '<button class=\"go-drop\" data-gene=\"'+gene+'\" data-type=\"expr\">DROP Expr</button> ' +
                '<button class=\"go-drop\" data-gene=\"'+gene+'\" data-type=\"splicing\">Splicing</button> ' +
                '<button class=\"go-drop\" data-gene=\"'+gene+'\" data-type=\"mae\">MAE</button>' +
                '</div>';
            };

            table.on('click', 'tr', function () {
              var tr = $(this);
              var row = table.row(tr);
              var gene = row.data()[row.data().length - 1];
              
              if (row.child.isShown()) {
                row.child.hide();
                tr.removeClass('shown');
              } else {
                row.child(format(gene)).show();
                tr.addClass('shown');
              }
            });

            table.on('click', '.go-gene', function(e) {
              e.stopPropagation();
              var gene = $(this).data('gene');
              Shiny.setInputValue('%s', {gene: gene, tab: 'gene'}, {priority: 'event'});
              Shiny.setInputValue('menu', 'gene', {priority: 'event'});
            });

            table.on('click', '.go-rna', function(e) {
              e.stopPropagation();
              var gene = $(this).data('gene');
              Shiny.setInputValue('%s', {gene: gene, tab: 'rna'}, {priority: 'event'});
            });

            table.on('click', '.go-drop', function(e) {
              e.stopPropagation();
              var gene = $(this).data('gene');
              var type = $(this).data('type');
              Shiny.setInputValue('%s', {gene: gene, tab: 'drop', type: type}, {priority: 'event'});
            });

          ", 
                                session$ns("nav_click"),
                                session$ns("nav_click"),
                                session$ns("nav_click")))
        )
        
      }, error=function(e){
        datatable(data.frame(Message="No variants loaded"))
      })
    })
    
    
    # =====================
    # NAVEGACIÓN
    # =====================
    observeEvent(input$nav_click, {
      
      navigating(TRUE)
      selected_gene(input$nav_click$gene)
      
      if(input$nav_click$tab == "rna"){
        updateTabsetPanel(session, "main_tabs", selected = "RNA Data")
        
      } else if(input$nav_click$tab == "drop"){
        updateTabsetPanel(session, "main_tabs", selected = "DROP")
        updateSelectInput(session, "drop_type", selected = input$nav_click$type)
        
      } else if(input$nav_click$tab == "variants"){
        updateTabsetPanel(session, "main_tabs", selected = "Variants (WES + WGS)")
      }
    })
    
    
    # =====================
    # RNA
    # =====================
    output$rna <- renderDT({
      
      df <- get_rna(con)
      validate(need(nrow(df) > 0, "No RNA data"))
      
      if(!is.null(selected_gene())){
        df <- df %>% filter(grepl(selected_gene(), gene_name, ignore.case = TRUE))
      }
      
      if(nzchar(input$gene)){
        df <- df %>% filter(grepl(input$gene, gene_name, ignore.case = TRUE))
      }
      
      df$gene_hidden <- df$gene_name
      
      datatable(
        df,
        escape = FALSE,
        selection = "none",
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          columnDefs = list(list(targets = ncol(df)-1, visible = FALSE))
        ),
        
        callback = JS(sprintf("
          
          var format = function(gene){
            return '<div style=\"padding:10px\">' +
              '<b>Explore:</b><br>' +
              '<button class=\"go-var\" data-gene=\"'+gene+'\">🔬 Variants</button><br>' +
              '<button class=\"go-gene\" data-gene=\"'+gene+'\">🧬 Gene Viewer</button><br><br>' +
              '<button class=\"go-drop\" data-gene=\"'+gene+'\" data-type=\"expr\">DROP Expr</button> ' +
              '<button class=\"go-drop\" data-gene=\"'+gene+'\" data-type=\"splicing\">Splicing</button> ' +
              '<button class=\"go-drop\" data-gene=\"'+gene+'\" data-type=\"mae\">MAE</button>' +
              '</div>';
          };

          table.on('click','tr',function(){
            var tr=$(this); var row=table.row(tr);
            var gene=row.data()[row.data().length-1];
            
            if(row.child.isShown()){row.child.hide();tr.removeClass('shown');}
            else{row.child(format(gene)).show();tr.addClass('shown');}
          });

          table.on('click','.go-var',function(e){
            e.stopPropagation();
            var gene=$(this).data('gene');
            Shiny.setInputValue('%s',{gene:gene,tab:'variants'},{priority:'event'});
          });

          table.on('click','.go-gene',function(e){
            e.stopPropagation();
            var gene=$(this).data('gene');
            Shiny.setInputValue('%s',{gene:gene,tab:'gene'},{priority:'event'});
            Shiny.setInputValue('menu','gene',{priority:'event'});
          });
          
          table.on('click','.go-drop',function(e){
            e.stopPropagation();
            var gene=$(this).data('gene');
            var type=$(this).data('type');
            Shiny.setInputValue('%s',{gene:gene,tab:'drop',type:type},{priority:'event'});
          });

        ", 
                              session$ns("nav_click"),
                              session$ns("nav_click"),
                              session$ns("nav_click")
        ))
      )
    })
    
    
    # =====================
    # DROP
    # =====================
    output$drop <- renderDT({
      
      df <- switch(input$drop_type,
                   "expr" = get_drop_expr(con),
                   "splicing" = get_drop_splicing(con),
                   "mae" = get_drop_mae(con))
      
      validate(need(nrow(df) > 0, "No DROP data"))
      
      gene_col <- if("gene_name" %in% colnames(df)) "gene_name" else "hgncSymbol"
      
      if(!is.null(selected_gene())){
        df <- df %>% filter(grepl(selected_gene(), .data[[gene_col]], ignore.case = TRUE))
      }
      
      if(nzchar(input$gene)){
        df <- df %>% filter(grepl(input$gene, .data[[gene_col]], ignore.case = TRUE))
      }
      
      df$gene_hidden <- df[[gene_col]]
      
      datatable(
        df,
        escape = FALSE,
        selection = "none",
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          columnDefs = list(list(targets = ncol(df)-1, visible = FALSE))
        ),
        
        callback = JS(sprintf("
          
          var format = function(gene){
            return '<div style=\"padding:10px\">' +
              '<b>Explore:</b><br>' +
              '<button class=\"go-var\" data-gene=\"'+gene+'\">🔬 Variants</button><br>' +
              '<button class=\"go-gene\" data-gene=\"'+gene+'\">🧬 Gene Viewer</button><br><br>' +
              '<button class=\"go-rna\" data-gene=\"'+gene+'\">RNA</button>' +
              '</div>';
          };

          table.on('click','tr',function(){
            var tr=$(this); var row=table.row(tr);
            var gene=row.data()[row.data().length-1];
            
            if(row.child.isShown()){row.child.hide();tr.removeClass('shown');}
            else{row.child(format(gene)).show();tr.addClass('shown');}
          });

          table.on('click','.go-var',function(e){
            e.stopPropagation();
            var gene=$(this).data('gene');
            Shiny.setInputValue('%s',{gene:gene,tab:'variants'},{priority:'event'});
          });

          table.on('click','.go-gene',function(e){
            e.stopPropagation();
            var gene=$(this).data('gene');
            Shiny.setInputValue('%s',{gene:gene,tab:'gene'},{priority:'event'});
            Shiny.setInputValue('menu','gene',{priority:'event'});
          });
          
          table.on('click','.go-rna',function(e){
            e.stopPropagation();
            var gene=$(this).data('gene');
            Shiny.setInputValue('%s',{gene:gene,tab:'rna'},{priority:'event'});
          });

        ", 
                              session$ns("nav_click"),
                              session$ns("nav_click"),
                              session$ns("nav_click")
        ))
      )
    })
    
  })
}