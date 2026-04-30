library(shiny)
library(shinycssloaders)
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
      
      div(class="filter-box",
          textInput(
            ns("gene"),
            "Filter by gene:",
            placeholder = "e.g. DPM1"
          )
      ),
      
      tabsetPanel(id = ns("main_tabs"),
                  
                  tabPanel("Variants (WES + WGS)",
                           
                           div(class = "filter-box",
                               
                               fluidRow(
                                 column(3,
                                        sliderInput(ns("af"), "Max AF:",
                                                    min = 0, max = 1,
                                                    value = 0.05, step = 0.01)
                                 ),
                                 column(3,
                                        checkboxGroupInput(ns("impact"), "Impact:",
                                                           choices = c("HIGH","MODERATE","LOW","MODIFIER"),
                                                           selected = c("HIGH","MODERATE","LOW","MODIFIER"))
                                 ),
                                 column(3,
                                        checkboxGroupInput(ns("source"), "Source:",
                                                           choices = c("WES","WGS","BOTH"),
                                                           selected = c("BOTH"))
                                 ),
                                 column(3,
                                        checkboxGroupInput(ns("inheritance"), "Inheritance:",
                                                           choices = c("de_novo","recessive","dominant","other"),
                                                           selected = c("de_novo","recessive","dominant","other"))
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
                               )
                           ),
                           
                           br(),
                           withSpinner(DTOutput(ns("variants")), type = 4, color = "#8b1e5b")
                  ),
                  
                  tabPanel("RNA Data",
                           withSpinner(DTOutput(ns("rna")), type = 4, color = "#2c7fb8")
                  ),
                  
                  tabPanel("DROP",
                           selectInput(ns("drop_type"),
                                       "Select dataset:",
                                       choices = c(
                                         "Aberrant Expression"="expr",
                                         "Aberrant Splicing"="splicing",
                                         "MAE (Monoallelic Expression)"="mae"
                                       )),
                           withSpinner(DTOutput(ns("drop")), type = 4, color = "#2c7fb8")
                  )
      )
  )
}

# =====================
# SERVER
# =====================
dataViewerServer <- function(id, pool, selected_gene){
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
      
      req(pool)
      
      tryCatch({
        
        df <- get_variants_with_inheritance(pool)
        
        if(nrow(df) == 0){
          return(datatable(data.frame(Message="No variants loaded")))
        }
        
        if(!is.null(selected_gene()) && selected_gene() != ""){
          df <- df %>% filter(grepl(selected_gene(), SYMBOL, ignore.case = TRUE))
        }
        
        if(nzchar(input$gene)){
          df <- df %>% filter(grepl(input$gene, SYMBOL, ignore.case = TRUE))
        }
        
        df <- df %>% filter(is.na(MAX_AF) | MAX_AF <= input$af)
        
        if(length(input$impact) > 0){
          df <- df %>% filter(IMPACT %in% input$impact)
        }
        
        if(length(input$source) > 0){
          df <- df %>% filter(source %in% input$source)
        }
        
        if(length(input$inheritance) > 0){
          df <- df %>% filter(inheritance_type %in% input$inheritance)
        }
        
        if(input$consequence != "ALL"){
          df <- df %>% filter(grepl(input$consequence, Consequence, ignore.case = TRUE))
        }
        
        if(input$variant_class != "ALL"){
          df <- df %>% filter(VARIANT_CLASS == input$variant_class)
        }
        
        if(nrow(df) == 0){
          return(datatable(data.frame(Message="No variants match filters")))
        }
        
        cols_to_show <- c(
          "ID","CHROM","POS","REF","ALT","FILTER",
          "CHILD_GT_N","CHILD_DP","CHILD_AD","CHILD_GQ",
          "SYMBOL","Gene",
          "Consequence","IMPACT","MAX_AF","VARIANT_CLASS",
          "SIFT_pred","Polyphen2_HVAR_pred",
          "source","inheritance_type","inheritance",
          "AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF","AA_AF","EA_AF",
          "gnomAD_AF","gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF",
          "gnomAD_EAS_AF","gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF"
        )
        
        cols_to_show <- intersect(cols_to_show, colnames(df))
        df <- df[, cols_to_show]
        
        colnames(df)[colnames(df) == "CHILD_GT_N"] <- "GT"
        colnames(df)[colnames(df) == "CHILD_DP"]   <- "DP"
        colnames(df)[colnames(df) == "CHILD_AD"]   <- "AD"
        colnames(df)[colnames(df) == "CHILD_GQ"]   <- "GQ"
        colnames(df)[colnames(df) == "SYMBOL"]     <- "Gene name"
        colnames(df)[colnames(df) == "Gene"]       <- "Gene ID"
        
        # ORDEN CONTROLADO
        df <- df %>%
          mutate(
            AF = AF,
            AFR_AF = AFR_AF,
            AMR_AF = AMR_AF,
            EAS_AF = EAS_AF,
            EUR_AF = EUR_AF,
            SAS_AF = SAS_AF,
            AA_AF = AA_AF,
            EA_AF = EA_AF
          )
        
        df$gene_hidden <- df$`Gene name`
        
        datatable(
          df,
          rownames = FALSE,
          escape = FALSE,
          selection = "none",
          options = list(
            scrollX = TRUE,
            pageLength = 10,
            columnDefs = list(
              list(targets = (ncol(df)-17):(ncol(df)-1), visible = FALSE)
            )
          ),
          
          callback = JS(sprintf("
  
  var format = function(rowData) {

    var n = rowData.length;

    var gene = rowData[n - 1];

    // ===== ORIGINAL AF =====
    var AF      = rowData[n - 17];
    var AFR_AF  = rowData[n - 16];
    var AMR_AF  = rowData[n - 15];
    var EAS_AF  = rowData[n - 14];
    var EUR_AF  = rowData[n - 13];
    var SAS_AF  = rowData[n - 12];
    var AA_AF   = rowData[n - 11];
    var EA_AF   = rowData[n - 10];

    // ===== GNOMAD AF =====
    var g_AF      = rowData[n - 9];
    var g_AFR     = rowData[n - 8];
    var g_AMR     = rowData[n - 7];
    var g_ASJ     = rowData[n - 6];
    var g_EAS     = rowData[n - 5];
    var g_FIN     = rowData[n - 4];
    var g_NFE     = rowData[n - 3];
    var g_OTH     = rowData[n - 2];

    return '<div style=\"padding:10px\">' +

      '<b>Explore:</b><br>' +
      '<button class=\"go-gene\" data-gene=\"'+gene+'\">Gene Info</button><br><br>' +
      '<button class=\"go-rna\" data-gene=\"'+gene+'\">RNA Info</button><br><br>' +
      '<button class=\"go-drop\" data-gene=\"'+gene+'\" data-type=\"expr\">DROP Expression</button><br>' +
      '<button class=\"go-drop\" data-gene=\"'+gene+'\" data-type=\"splicing\">DROP Splicing</button><br>' +
      '<button class=\"go-drop\" data-gene=\"'+gene+'\" data-type=\"mae\">DROP MAE</button>' +

      '<hr>' +

      '<b>Additional info:</b><br>' +

      '<button class=\"show-freq\">Population frequencies</button> ' +
      '<button class=\"show-gnomad\">gnomAD frequencies</button>' +

      // ===== ORIGINAL =====
      '<div class=\"freq-box\" style=\"display:none;margin-top:10px\">' +
        '<b>Frequencies:</b><br>' +
        'AF: ' + AF + '<br>' +
        'AFR: ' + AFR_AF + '<br>' +
        'AMR: ' + AMR_AF + '<br>' +
        'EAS: ' + EAS_AF + '<br>' +
        'EUR: ' + EUR_AF + '<br>' +
        'SAS: ' + SAS_AF + '<br>' +
        'AA: ' + AA_AF + '<br>' +
        'EA: ' + EA_AF +
      '</div>' +

      // ===== GNOMAD =====
      '<div class=\"gnomad-box\" style=\"display:none;margin-top:10px\">' +
        '<b>gnomAD:</b><br>' +
        'AF: ' + g_AF + '<br>' +
        'AFR: ' + g_AFR + '<br>' +
        'AMR: ' + g_AMR + '<br>' +
        'ASJ: ' + g_ASJ + '<br>' +
        'EAS: ' + g_EAS + '<br>' +
        'FIN: ' + g_FIN + '<br>' +
        'NFE: ' + g_NFE + '<br>' +
        'OTH: ' + g_OTH +
      '</div>' +

    '</div>';
  };

  table.on('click','tr',function(){
    var tr=$(this); 
    var row=table.row(tr);
    
    if(row.child.isShown()){
      row.child.hide(); 
      tr.removeClass('shown');
    } else {
      row.child(format(row.data())).show(); 
      tr.addClass('shown');
    }
  });

  // ===== TOGGLES =====
  table.on('click','.show-freq',function(e){
    e.stopPropagation();
    $(this).parent().find('.freq-box').toggle();
  });

  table.on('click','.show-gnomad',function(e){
    e.stopPropagation();
    $(this).parent().find('.gnomad-box').toggle();
  });

  // ===== NAV =====
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
        
      }, error = function(e){
        print(e)
        datatable(data.frame(Message="Error loading variants"))
      })
    })
    
    # =====================
    # NAV
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
      
      df <- get_rna(pool)
      
      if(nrow(df)==0){
        return(datatable(data.frame(Message="No RNA data")))
      }
      
      if(!is.null(selected_gene()) && selected_gene()!=""){
        df <- df %>% filter(grepl(selected_gene(), gene_name, ignore.case=TRUE))
      }
      
      if(nzchar(input$gene)){
        df <- df %>% filter(grepl(input$gene, gene_name, ignore.case=TRUE))
      }
      
      df$gene_hidden <- df$gene_name
      
      datatable(
        df,
        escape = FALSE,
        selection = "none",
        options=list(
          scrollX=TRUE,pageLength=10,
          columnDefs=list(list(targets=ncol(df)-1,visible=TRUE))
        ),
        
        callback = JS(sprintf("
          
          var format=function(gene){
            return '<div style=\"padding:10px\">' +
              '<b>Explore:</b><br>' +
              '<button class=\"go-var\" data-gene=\"'+gene+'\">Variants Info</button><br>' +
              '<button class=\"go-gene\" data-gene=\"'+gene+'\">Gene Info</button><br><br>' +
              '<button class=\"go-drop\" data-gene=\"'+gene+'\" data-type=\"expr\">DROP Expression</button> ' +
              '<button class=\"go-drop\" data-gene=\"'+gene+'\" data-type=\"splicing\">DROP Splicing</button> ' +
              '<button class=\"go-drop\" data-gene=\"'+gene+'\" data-type=\"mae\">DROP MAE</button>' +
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
                   expr=get_drop_expr(pool),
                   splicing=get_drop_splicing(pool),
                   mae=get_drop_mae(pool))
      
      if(nrow(df)==0){
        return(datatable(data.frame(Message="No DROP data")))
      }
      
      gene_col <- if("gene_name" %in% colnames(df)) "gene_name" else "hgncSymbol"
      
      if(!is.null(selected_gene()) && selected_gene()!=""){
        df <- df %>% filter(grepl(selected_gene(), .data[[gene_col]], ignore.case=TRUE))
      }
      
      if(nzchar(input$gene)){
        df <- df %>% filter(grepl(input$gene, .data[[gene_col]], ignore.case=TRUE))
      }
      
      df$gene_hidden <- df[[gene_col]]
      
      datatable(
        df,
        escape=FALSE,
        selection="none",
        options=list(
          scrollX=TRUE,pageLength=10,
          columnDefs=list(list(targets=ncol(df)-1,visible=FALSE))
        ),
        
        callback = JS(sprintf("
          
          var format=function(gene){
            return '<div style=\"padding:10px\">' +
              '<b>Explore:</b><br>' +
              '<button class=\"go-var\" data-gene=\"'+gene+'\">Variants Info</button><br>' +
              '<button class=\"go-gene\" data-gene=\"'+gene+'\">Gene Info</button><br><br>' +
              '<button class=\"go-rna\" data-gene=\"'+gene+'\">RNA Info</button>' +
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