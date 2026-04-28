library(shiny)
options(shiny.maxRequestSize = 20000 * 1024^2)
library(shinyjs)
library(DBI)
library(duckdb)

APP_DIR <- normalizePath(getwd())

# =========================
# CONNECTION TO DATABASE
# =========================
con <- dbConnect(
  duckdb(),
  normalizePath("../db/genomic.duckdb")
)

# =========================
# SCRIPTS
# =========================
source("modules/utils_processing.R")
source("queries/queries_dataViewer.R")
source("queries/queries_geneViewer.R")
source("queries/queries_QC.R")

# =========================
# MODULES
# =========================
source("modules/home.R")
source("modules/introduction.R")
source("modules/loadData.R")
source("modules/dataViewer.R")
source("modules/geneViewer.R")
source("modules/QCviewer.R")
source("modules/coverageViewer.R")

# =========================
# SERVIR COVERAGE DESDE SHINY (SIN CORS)
# =========================
addResourcePath(
  "coverage",
  normalizePath("../data/coverage")
)

# =========================
# UI
# =========================
ui <- fluidPage(
  
  useShinyjs(),
  
  tags$head(
    tags$link(rel="stylesheet", type="text/css", href="styles.css"),
    
    tags$script(HTML("
      function changePage(page){
        Shiny.setInputValue('menu', page, {priority: 'event'});

        document.querySelectorAll('.menu-item').forEach(el => {
          el.classList.remove('active');
        });

        document.getElementById('btn_' + page).classList.add('active');
      }
    "))
  ),
  
  div(class="app-container",
      
      # SIDEBAR
      div(class="sidebar",
          
          div(id="btn_home", class="menu-item active",
              onclick="changePage('home')", "🏠 Home"),
          
          div(id="btn_intro", class="menu-item",
              onclick="changePage('intro')", "📄 Introduction"),
          
          div(id="btn_load", class="menu-item",
              onclick="changePage('load')", "📥 Load Data"),
          
          div(id="btn_viewer", class="menu-item",
              onclick="changePage('viewer')", "🔎 Data Viewer"),
          
          div(id="btn_gene", class="menu-item",
              onclick="changePage('gene')", "🧬 Gene Viewer"),
          
          div(id="btn_qc", class="menu-item",
              onclick="changePage('qc')", "📊 QC Viewer"),
          
          div(id="btn_coverage", class="menu-item",
              onclick="changePage('coverage')", "📈 Coverage Viewer"),
          
          div(id="btn_plots", class="menu-item",
              onclick="changePage('plots')", "📉 Plots")
      ),
      
      # MAIN
      div(class="main",
          uiOutput("page")
      )
  )
)

# =========================
# SERVER
# =========================
server <- function(input, output, session){
  
  current_page <- reactiveVal("home")
  
  observeEvent(input$menu, {
    current_page(input$menu)
  })
  
  # =========================
  # PAGE RENDER
  # =========================
  output$page <- renderUI({
    
    switch(current_page(),
           
           "home" = homeUI("home"),
           "intro" = introUI("intro"),
           "load" = loadUI("loadData"),
           "viewer" = dataViewerUI("viewer"),
           "gene" = geneViewerUI("gene"),
           "qc" = QCviewerUI("qc"),
           "coverage" = coverageViewerUI("coverage"),
           
           "plots" = div("Plots")
    )
  })

  # =========================
  # ACTIVE MODULES
  # =========================
  selected_gene <- reactiveVal(NULL)
  
  
  homeServer("home")
  introServer("intro")
  dataViewerServer("viewer", con, selected_gene)
  geneViewerServer("gene", con, selected_gene)
  QCviewerServer("qc", con)
  coverageViewerServer("coverage")
  
  # =========================
  # CLOSE DB PROPERLY
  # =========================
  session$onSessionEnded(function(){
    dbDisconnect(con, shutdown = TRUE)
  })
}

# =========================
# RUN APP
# =========================
shinyApp(ui, server)