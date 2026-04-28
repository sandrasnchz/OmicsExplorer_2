introUI <- function(id){
  ns <- NS(id)
  
  div(class="content",
      
      h2("📄 | INTRODUCTION"),
      
      # TÍTULO CENTRAL
      div(class="home-center",
          h1("OmicsExplorer"),
          p(class="subtitle", "Multi-omics integration for clinical variant interpretation")
      ),
      
      # CAJA TEXTO
      div(class="intro-box",
          
          p(strong("OmicsExplorer"), 
            " is an interactive application developed for the integration and clinical interpretation of genomic and transcriptomic data, with a particular focus on patients with rare diseases without a definitive diagnosis."
          ),
          
          p("This tool has been developed as part of a Master’s Thesis within the Master’s Degree in Bioinformatics and Computational Biology at the Universidad Autónoma de Madrid, in collaboration with the Instituto de Salud Carlos III, and is framed within the SpainUDP (Spain Undiagnosed Rare Diseases Program)."),
          
          p("Advances in Next-Generation Sequencing technologies (NGS), such as Whole Exome Sequencing (WES), Whole Genome Sequencing (WGS), and RNA-seq, have enabled the identification of genetic variants and transcriptomic alterations associated with disease. However, the interpretation of these data remains a major challenge in clinical genomics."),
          
          p("OmicsExplorer addresses this challenge by providing a unified environment for the exploration and prioritization of candidate variants through the integration of genomic annotations and transcriptomic evidence. The platform incorporates filtering strategies, functional annotations, and interactive visualizations to support the identification of clinically relevant variants."),
          
          p("By combining multi-omics data with intuitive analytical tools, OmicsExplorer aims to facilitate data interpretation and support clinical decision-making in the context of precision medicine.")
      )
  )
}

introServer <- function(id){
  moduleServer(id, function(input, output, session){})
}