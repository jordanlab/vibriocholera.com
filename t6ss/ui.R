library('shiny')
library('shinyjs')
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
#
shinyUI(
  fluidPage(
    useShinyjs(),
    theme = 'cerulean.css',
    tags$head(includeScript("www/ga.js")),
    tags$head(tags$script(HTML(
     "var _paq = window._paq || [];
       _paq.push(['setCookieDomain', '*.vibriocholera.com']);
       _paq.push(['trackPageView']);
       _paq.push(['enableLinkTracking']);
       (function() {
         var u='//matomo.chande.science/';
         _paq.push(['setTrackerUrl', u+'matomo.php']);
         _paq.push(['setSiteId', '2']);
         var d=document, g=d.createElement('script'), s=d.getElementsByTagName('script')[0];
         g.type='text/javascript'; g.async=true; g.defer=true; g.src=u+'matomo.js'; s.parentNode.insertBefore(g,s);
       })();
     "
    ))),
    includeHTML('header.html'),
    headerPanel(title = "T6SS Predictor", windowTitle = 'Predict T6SS Proteins'),
    sidebarPanel(
      strong("1. Upload your genome assembly "),
      #br(),
      fileInput('fastaFile', 'Upload FASTA file', accept = c(".fasta", ".fas", ".fna")),
      a('Example genome sequence file', href='BGT49_PacBio_Assembly.fasta'),
      p("You may upload files up to 8MB"),
      p("Do not upload proteins or predicted CDS sequences"),
      br(),
      strong("2. Click the button and be patient. Predictions take ~5 mins"),
       hr(),
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
          "Prediction results", value="results",
          h3("Predicted Type VI Secretion System Loci and Proteins"),
          # tags$hr(),
          htmlOutput("link")
      )),


      includeHTML('footer.html')

    )
  )
)

