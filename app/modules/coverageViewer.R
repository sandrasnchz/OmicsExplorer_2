library(shiny)

# =====================
# UI
# =====================
coverageViewerUI <- function(id){
  ns <- NS(id)
  
  div(class="content",
      
      h2("📈 | COVERAGE VIEWER"),
      
      # ===== INPUT =====
      div(class="filter-box",
          textInput(
            ns("region"),
            "Search by region or gene name",
            placeholder = "e.g. 1:1-1000000 or DPM1"
          )
      ),
      
      # ===== ACTION BOX =====
      div(class="table-box",
          
          tags$h4("Genome visualization", style="color:#8b1e5b;"),
          
          tags$p(
            "Open IGV with WES, WGS and RNA coverage tracks.",
            style="color:#555;"
          ),
          
          actionButton(
            ns("open"),
            "Open IGV on the explorer",
            class="btn-download"
          )
      )
  )
}

# =====================
# SERVER
# =====================
coverageViewerServer <- function(id){
  moduleServer(id, function(input, output, session){
    
    observeEvent(input$open, {
      
      # -------------------------
      # BUSCAR ARCHIVOS
      # -------------------------
      coverage_dir <- normalizePath("../data/coverage")
      
      files <- list.files(coverage_dir, full.names = TRUE)
      
      wes_path <- tail(files[grepl("^WES.*\\.bw$", basename(files))], 1)
      wgs_path <- tail(files[grepl("^WGS.*\\.bw$", basename(files))], 1)
      rna_path <- tail(files[grepl("^RNA.*\\.bw$", basename(files))], 1)
      
      # -------------------------
      # VALIDACIÓN 
      # -------------------------
      if (length(wes_path) == 0 || length(wgs_path) == 0 || length(rna_path) == 0) {
        showNotification("Missing coverage files", type = "error")
        return()
      }
      
      # -------------------------
      # URLs 
      # -------------------------
      wes_url <- paste0("coverage/", basename(wes_path))
      wgs_url <- paste0("coverage/", basename(wgs_path))
      rna_url <- paste0("coverage/", basename(rna_path))
      
      # -------------------------
      # HTML IGV
      # -------------------------
      html <- paste0(
        "<!DOCTYPE html>
<html>
<head>
<script src='https://cdn.jsdelivr.net/npm/igv@2.15.4/dist/igv.min.js'></script>
</head>
<body>

<div id='igvDiv' style='height:600px;'></div>

<script>

var options = {
  genome: 'hg19',
  locus: '", input$region, "',
  tracks: [

    {
      name:'WGS',
      type:'wig',
      format:'bigwig',
      url:'", wgs_url, "',
      color:'black'
    },

    {
      name:'WES',
      type:'wig',
      format:'bigwig',
      url:'", wes_url, "',
      color:'blue'
    },

    {
      name:'RNA',
      type:'wig',
      format:'bigwig',
      url:'", rna_url, "',
      color:'red'
    }

  ]
};

igv.createBrowser(document.getElementById('igvDiv'), options);

</script>

</body>
</html>"
      )
      
      writeLines(html, "www/tmp_igv.html")
      
      url <- paste0(
        session$clientData$url_protocol, "//",
        session$clientData$url_hostname, ":",
        session$clientData$url_port,
        "/tmp_igv.html"
      )
      
      browseURL(url)
      
    })
    
  })
}