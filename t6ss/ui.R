library('shiny')
library('genoPlotR')

# actionButton with dark color
nx.actionButton = function (inputId, label, icon = NULL) {
  if (!is.null(icon))
    buttonContent <- list(icon, label)
  else
    buttonContent <- label
  tags$button(id = inputId,
              type = "button",
              class = "btn btn-primary action-button",
              buttonContent)
}

shinyUI(
  fluidPage(
    theme = 'index.min.css',
    tags$head(includeScript("www/ga.js")),
    includeHTML('header.html'),
    headerPanel(title = "T6SS Predictor", windowTitle = 'Predict T6SS Proteins'),
    sidebarPanel(
      radioButtons(
        "seqType",
        "1. Are you providing a genome sequence (FASTA nucleotide) or predicted proteins (FASTA amino acid)?",
        c("Proteins (faa)" = "protein",
          "Genome (fna)" = "genomic")
        
      ),
      strong("2."),
      #br(),
      fileInput('fastaFile', 'Upload FASTA file'),
      a('Example genome sequence file', href='bgt8_genomic.fna'),
      br(),
      a('Example protein sequence file', href='bgt8_cds.faa'),
      br(),
      conditionalPanel(
        condition = "input.seqType == 'protein'",
        strong("3. Click the button and be patient. Predictions take ~5 mins"),
        hr()
        
      ),
      conditionalPanel(
        radioButtons(
          "predict",
          "3. Are you providing your own annotation (GFF)?",
          c("No" = "pred",
            "Yes" = "nopred")
        ),
        condition = "input.seqType == 'genomic'",
      conditionalPanel(
        condition = "input.predict == 'pred'",
        strong("4. Click the button and be patient. Predictions take ~5 mins"),
        hr()

      ),
      conditionalPanel(
        condition = "input.predict == 'nopred'",
        strong('4. Upload your GFF file'),
        fileInput('gffFile', 'Upload GFF file'),
        a('Example GFF file', href='bgt8.gff'),
        br(),
        strong("5. Click the button and be patient. Predictions take ~5 mins"),
        hr()

      )),
      
      nx.actionButton('submit', 'Predict T6SS'),
      tags$hr()
      
      
    ),
    
    mainPanel(# css hack to move the progress bar to a lower place
      # from https://gist.github.com/johndharrison/9578241
      # tags$link(rel = 'stylesheet', type = 'text/css', href = 'progbar.css'),
      tabsetPanel(
        id = 'vchot6ss',
        tabPanel(
          "Introduction",
          p(
            "Type VI Secretion Systems (T6SS) are bacterial weaponry designed to poison, and potentially kill, neighboring bacteria and eukaryotic cells. T6SS relies on a set of structural proteins, encoded in what's called the Large Cluster, and a set of effector/immunity pairs encoded on Auxiliary clusters. This tool attempts to predict Auxiliary Cluster loci and their VgrG and Effector proteins"
          ),
          img(src = 't6loci.png', align = "center")
        ),
        
        tabPanel(
          "T6SS Predictions",
          h3("Predicted Type VI Secretion System Loci and Proteins"),
          tags$hr(),
          plotOutput("plot"),
          tags$hr(),
          uiOutput("dlProts"),
          uiOutput("dlPreds")
        )
        
      ),
      
      includeHTML('footer.html')
      
    )
  )
)
