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
                           withSpinner(DTOutput(ns("rna")), type = 4, color = "#8b1e5b")
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
    # FUNCIÓN PARA OBTENER FLAGS DE DROP (PARA RESUMEN EN RNA)
    # =====================
    get_drop_flags <- function(pool){
      
      expr <- get_drop_expr(pool) %>% mutate(drop_expr = TRUE)
      spl  <- get_drop_splicing(pool) %>% mutate(drop_splicing = TRUE)
      mae  <- get_drop_mae(pool) %>% mutate(drop_mae = TRUE)
      
      gene_col_expr <- if("gene_name" %in% colnames(expr)) "gene_name" else "hgncSymbol"
      gene_col_spl  <- if("gene_name" %in% colnames(spl)) "gene_name" else "hgncSymbol"
      gene_col_mae  <- if("gene_name" %in% colnames(mae)) "gene_name" else "hgncSymbol"
      
      expr <- expr %>% select(gene = all_of(gene_col_expr), drop_expr)
      spl  <- spl  %>% select(gene = all_of(gene_col_spl), drop_splicing)
      mae  <- mae  %>% select(gene = all_of(gene_col_mae), drop_mae)
      
      df <- full_join(expr, spl, by="gene") %>%
        full_join(mae, by="gene") %>%
        mutate(across(starts_with("drop"), ~replace_na(., FALSE)))
      
      return(df)
    }
    
    
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
        
        df$gene_hidden <- df$`Gene name`
        
        datatable(
          df,
          rownames = FALSE,
          escape = FALSE,
          selection = "none",
          extensions = 'FixedHeader',
          options = list(
            scrollX = TRUE,
            pageLength = 10,
            fixedHeader = TRUE,
            columnDefs = list(
              list(targets = (ncol(df)-17):(ncol(df)-1), visible = FALSE)
            )
          ),
          
          callback = JS(sprintf("

var format = function(rowData){

  var n = rowData.length;
  var gene = rowData[n - 1];
  var uid = Math.random().toString(36).substring(2,9);

  // ===== ORIGINAL =====
  var AF      = rowData[n - 17];
  var AFR_AF  = rowData[n - 16];
  var AMR_AF  = rowData[n - 15];
  var EAS_AF  = rowData[n - 14];
  var EUR_AF  = rowData[n - 13];
  var SAS_AF  = rowData[n - 12];
  var AA_AF   = rowData[n - 11];
  var EA_AF   = rowData[n - 10];

  // ===== GNOMAD =====
  var g_AF  = rowData[n - 9];
  var g_AFR = rowData[n - 8];
  var g_AMR = rowData[n - 7];
  var g_ASJ = rowData[n - 6];
  var g_EAS = rowData[n - 5];
  var g_FIN = rowData[n - 4];
  var g_NFE = rowData[n - 3];
  var g_OTH = rowData[n - 2];

  return '<div style=\"padding:10px\">' +

    '<ul class=\"nav nav-tabs\">' +

      '<li class=\"nav-item\">' +
        '<a class=\"nav-link active freq-tab\" data-target=\"#pop-'+uid+'\">Population</a>' +
      '</li>' +

      '<li class=\"nav-item\">' +
        '<a class=\"nav-link freq-tab\" data-target=\"#gnom-'+uid+'\">gnomAD</a>' +
      '</li>' +

    '</ul>' +

    '<div class=\"tab-content\" style=\"margin-top:10px\">' +

      '<div class=\"tab-pane active\" id=\"pop-'+uid+'\">' +
        'AF: '+AF+'<br>' +
        'AFR: '+AFR_AF+'<br>' +
        'AMR: '+AMR_AF+'<br>' +
        'EAS: '+EAS_AF+'<br>' +
        'EUR: '+EUR_AF+'<br>' +
        'SAS: '+SAS_AF+'<br>' +
        'AA: '+AA_AF+'<br>' +
        'EA: '+EA_AF +
      '</div>' +

      '<div class=\"tab-pane\" id=\"gnom-'+uid+'\">' +
        'AF: '+g_AF+'<br>' +
        'AFR: '+g_AFR+'<br>' +
        'AMR: '+g_AMR+'<br>' +
        'ASJ: '+g_ASJ+'<br>' +
        'EAS: '+g_EAS+'<br>' +
        'FIN: '+g_FIN+'<br>' +
        'NFE: '+g_NFE+'<br>' +
        'OTH: '+g_OTH +
      '</div>' +

    '</div>' +

    '<hr>' +

    '<div class=\"explore-box\">' +
      '<div class=\"explore-title\">Explore</div>' +
      '<div class=\"explore-desc\">Navigate to gene viewer or RNA data table</div>' +
      '<div style=\"display:flex; gap:10px; margin-top:8px;\">' +
      '<button class=\"go-gene\" data-gene=\"'+gene+'\">View gene info</button><br>' +
      '<button class=\"go-rna\" data-gene=\"'+gene+'\">Go to RNA Data table</button>' +
      '</div>' +
      
    '</div>' +

  '</div>';
};

// abrir/cerrar
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

// tabs
table.on('click','.freq-tab',function(e){
  e.stopPropagation();
  
  var container = $(this).closest('div');
  
  container.find('.nav-link').removeClass('active');
  $(this).addClass('active');
  
  var target = $(this).data('target');
  container.find('.tab-pane').removeClass('active');
  container.find(target).addClass('active');
});

// navegación
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
                                session$ns("nav_click")
          ))
        )
        
      }, error = function(e){
        print(e)
        datatable(data.frame(Message="Error loading variants"))
      })
    })
    
    # =====================
    # NAVEGACIÓN CONTROLADA DESDE BOTONES EN VARIANTES Y RNA
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
    # RNA DATA
    # =====================
    output$rna <- renderDT({
      
      df <- get_rna(pool)
      
      if(nrow(df)==0){
        return(datatable(data.frame(Message="No RNA data")))
      }
      
      # ===== JOIN DROP =====
      drop_flags <- get_drop_flags(pool)
      
      df <- df %>%
        left_join(drop_flags, by = c("gene_name" = "gene"))
      
      # ===== CREAR COLUMNA RESUMEN =====
      df <- df %>%
        mutate(
          drop_expr = replace_na(drop_expr, FALSE),
          drop_splicing = replace_na(drop_splicing, FALSE),
          drop_mae = replace_na(drop_mae, FALSE),
          
          DROP_status = paste0(
            "expr: ", drop_expr, "<br>",
            "splicing: ", drop_splicing, "<br>",
            "mae: ", drop_mae
          )
        )
      
      # ===== FILTROS =====
      if(!is.null(selected_gene()) && selected_gene()!=""){
        df <- df %>% filter(grepl(selected_gene(), gene_name, ignore.case=TRUE))
      }
      
      if(nzchar(input$gene)){
        df <- df %>% filter(grepl(input$gene, gene_name, ignore.case=TRUE))
      }
      
      # ===== HIDDEN COLUMN =====
      df$gene_hidden <- df$gene_name
      
      df <- df %>% select(-drop_expr, -drop_splicing, -drop_mae)
      
      datatable(
        df,
        escape = FALSE,
        selection = "none",
        options=list(
          scrollX = TRUE,
          autoWidth = TRUE,
          pageLength = 10,
          columnDefs=list(
            list(targets = which(colnames(df) == "gene_hidden") - 1, visible = FALSE)
          )
        ),
        
        callback = JS(sprintf("
      
      var format=function(rowData){

        var gene=rowData[rowData.length-1];
        var uid = Math.random().toString(36).substring(2,9);

        return '<div style=\"padding:10px\">' +

          '<b>DROP summary</b><br><br>' +

          '<ul class=\"nav nav-tabs\" role=\"tablist\">' +

            '<li class=\"nav-item\">' +
              '<a class=\"nav-link active drop-tab\" data-gene=\"'+gene+'\" data-type=\"expr\" data-uid=\"'+uid+'\" href=\"#\">Expression</a>' +
            '</li>' +

            '<li class=\"nav-item\">' +
              '<a class=\"nav-link drop-tab\" data-gene=\"'+gene+'\" data-type=\"splicing\" data-uid=\"'+uid+'\" href=\"#\">Splicing</a>' +
            '</li>' +

            '<li class=\"nav-item\">' +
              '<a class=\"nav-link drop-tab\" data-gene=\"'+gene+'\" data-type=\"mae\" data-uid=\"'+uid+'\" href=\"#\">MAE</a>' +
            '</li>' +

          '</ul>' +

          '<div class=\"tab-content\" style=\"margin-top:10px\">' +
            '<div class=\"tab-pane active\" id=\"expr-'+uid+'\">Click tab</div>' +
            '<div class=\"tab-pane\" id=\"splicing-'+uid+'\"></div>' +
            '<div class=\"tab-pane\" id=\"mae-'+uid+'\"></div>' +
          '</div>' +

          '<hr>' +

          '<div class=\"explore-box\">' +
            '<div class=\"explore-title\">Explore variants</div>' +
            '<div class=\"explore-desc\">Go to variants table filtered by this gene</div>' +
            '<button class=\"go-var\" data-gene=\"'+gene+'\">View variants</button>' +
          '</div>' +

        '</div>';
      };

      table.on('click','tr',function(){
        var tr=$(this); var row=table.row(tr);
        if(row.child.isShown()){
          row.child.hide();tr.removeClass('shown');
        } else {
          row.child(format(row.data())).show();tr.addClass('shown');
        }
      });

      table.on('click','.drop-tab',function(e){
        e.preventDefault();
        e.stopPropagation();
        
        var gene=$(this).data('gene');
        var type=$(this).data('type');
        var uid=$(this).data('uid');
        
        var target = '#'+type+'-'+uid;
        var container = $(target);
        
        $(this).closest('ul').find('.nav-link').removeClass('active');
        $(this).addClass('active');
        
        container.siblings().removeClass('active');
        container.addClass('active');
        
        container.html('Loading...');
        
        Shiny.setInputValue('%s',{
          gene:gene,
          type:type,
          uid:uid
        },{priority:'event'});
      });

      table.on('click','.go-var',function(e){
        e.stopPropagation();
        var gene=$(this).data('gene');
        Shiny.setInputValue('%s',{gene:gene,tab:'variants'},{priority:'event'});
      });

      Shiny.addCustomMessageHandler('drop_render', function(msg) {
        $(msg.target).html(msg.html);
      });

      setTimeout(function(){
        table.columns.adjust();
      }, 300);

    ",
                              session$ns("drop_request"),
                              session$ns("nav_click")
        ))
      )
    })
    
    # =====================
    # DROP LOADER (FUERA)
    # =====================
    observeEvent(input$drop_request, {
      
      gene <- input$drop_request$gene
      type <- input$drop_request$type
      uid  <- input$drop_request$uid
      
      target <- paste0("#", type, "-", uid)
      
      df <- switch(type,
                   expr = get_drop_expr(pool),
                   splicing = get_drop_splicing(pool),
                   mae = get_drop_mae(pool))
      
      gene_col <- if("gene_name" %in% colnames(df)) "gene_name" else "hgncSymbol"
      
      df <- df %>% filter(.data[[gene_col]] == gene)
      
      if(nrow(df) == 0){
        html <- "<i>No data</i>"
      } else {
        html <- paste0(
          "<table style='font-size:12px;width:100%'>",
          "<tr>",
          paste0("<th>", colnames(df), "</th>", collapse=""),
          "</tr>",
          paste0(
            apply(df, 1, function(row){
              paste0("<tr>",
                     paste0("<td>", row, "</td>", collapse=""),
                     "</tr>")
            }), collapse=""
          ),
          "</table>"
        )
      }
      
      session$sendCustomMessage(
        "drop_render",
        list(html = html, target = target)
      )
    })
    
  })
}