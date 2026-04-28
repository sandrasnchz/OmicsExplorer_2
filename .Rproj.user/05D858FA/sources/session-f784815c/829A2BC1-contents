library(shiny)

# =====================
# UI
# =====================
coverageViewerUI <- function(id){
  ns <- NS(id)
  
  tagList(
    
    h2("📈 | Coverage Viewer"),
    
    textInput(ns("region"), "Find by region or gene name", placeholder = "e.g. 1:1-1000000 or DPM1"),
    
    actionButton(ns("open"), "Open IGV on the explorer")
  )
}

# =====================
# SERVER
# =====================
coverageViewerServer <- function(id){
  moduleServer(id, function(input, output, session){
    
    observeEvent(input$open, {
      
      # -------------------------
      # BUSCAR ARCHIVOS AUTOMÁTICAMENTE
      # -------------------------
      coverage_dir <- normalizePath("../data/coverage")
      
      files <- list.files(coverage_dir, full.names = TRUE)
      
      wes_path <- tail(files[grepl("^WES.*\\.bw$", basename(files))], 1)
      wgs_path <- tail(files[grepl("^WGS.*\\.bw$", basename(files))], 1)
      rna_path <- tail(files[grepl("^RNA.*\\.bw$", basename(files))], 1)
      gtf_path <- tail(files[grepl("^genes.*\\.gtf$", basename(files))], 1)
      
      # -------------------------
      # VALIDACIÓN
      # -------------------------
      if (length(wes_path) == 0 || length(wgs_path) == 0 || length(rna_path) == 0) {
        showNotification("Missing coverage files", type = "error")
        return()
      }
      
      # -------------------------
      # URLs (mismo puerto → sin CORS)
      # -------------------------
      wes_url <- paste0("coverage/", basename(wes_path))
      wgs_url <- paste0("coverage/", basename(wgs_path))
      rna_url <- paste0("coverage/", basename(rna_path))
      
      gtf_url <- if(length(gtf_path) > 0) {
        paste0("coverage/", basename(gtf_path))
      } else {
        NULL
      }
      
      # -------------------------
      # TRACK GTF (opcional)
      # -------------------------
      gtf_track <- ""
      
      if (!is.null(gtf_url)) {
        gtf_track <- paste0(
          "{
            name:'Genes',
            type:'annotation',
            format:'gtf',
            url:'", gtf_url, "'
          },"
        )
      }
      
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

", gtf_track, "

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