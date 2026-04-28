homeUI <- function(id){
  ns <- NS(id)
  
  div(class="content",
      
      h2("🏠 | HOME"),
      
      # TÍTULO CENTRAL
      div(class="home-center",
          h1("OmicsExplorer"),
          p(class="subtitle", "Multi-omics integration for clinical variant interpretation"),
          
          div(class="description-box",
              "OmicsExplorer is an interactive Shiny application designed to support the analysis and prioritization of genetic variants through the integration of genomic (WES/WGS) and transcriptomic (RNA-seq) data."
          )
      ),
      
      # FEATURES + USE CASE
      div(class="home-grid",
          
          # LEFT
          div(
            div(class="section-title", "🧪 Key Features"),
            div(class="card",
                tags$ul(
                  tags$li("Integration of genomic and transcriptomic data"),
                  tags$li("Variant filtering and functional annotation"),
                  tags$li("Gene-level exploration and visualization"),
                  tags$li("Interactive plots and data export")
                )
            )
          ),
          
          # RIGHT
          div(
            div(class="section-title", "🧬 Use Case"),
            div(class="card",
                "Designed for researchers and clinical geneticists working on rare diseases, OmicsExplorer facilitates the identification and prioritization of candidate variants through an intuitive and unified analysis environment."
            )
          )
      ),
      
      # FOOTER
      div(class="footer",
          "Developed as part of a Master’s Thesis in Bioinformatics and Computational Biology (UAM)"
      )
  )
}

homeServer <- function(id){
  moduleServer(id, function(input, output, session){})
}