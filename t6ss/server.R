library('shiny')
library('genoPlotR')
library(uuid)
library(shinyjs)

options(shiny.maxRequestSize=8*1024^2, shiny.reactlog=FALSE)
shinyServer(function(input, output, session) {
  hideTab("vchot6ss", "results")
  observeEvent(input$submit, {
      shiny::validate(
        need(input$fastaFile, "Please select a file to upload")
      )
      id = UUIDgenerate(use.time=F)
      dataFolder = file.path('/home/blast/prediction_server/server/www/results/',id)
      dir.create(dataFolder)
      file.copy('/home/blast/prediction_server/server/index.html', dataFolder)

      showTab("vchot6ss", "results")
      updateTabsetPanel(session, "vchot6ss", selected = "results")
      # print(input$tabBox_id)
      # sleep(4)
      sub_path = file.path('/home/blast/prediction_server/server/submissions',paste0(id,".toProcess"))
      # out_path = file.path('/home/blast/prediction_server/server/submissions',id,"output")
      dir.create(sub_path)
      filename = paste0(input$fastaFile$name,".toProcess")
      fasta = file.path(sub_path,filename)
      file.copy(input$fastaFile$datapath, fasta)
      queueSize <- length(list.files('/home/blast/prediction_server/server/submissions', "*.toProces"))
      if (queueSize-1 == 1){
        outText <- paste0("There is <b>1</b> submission ahead of you")
      }
      else{
        outText <- paste0("There are <b>", queueSize-1, "</b> submissions ahead of you")
      }
       output$link <- renderUI({
        HTML(paste("","","<h3>Your submission has been added to the queue!</h3><br/>",
          "<h4>" ,outText,
          "</h4><br/>","<h3><a href='https://t6ss.vibriocholera.com/results/",id,"/index.html' target=_blank>Click here for your results</a></h3><br/><p>Results will be available for at least 7 days</p>", sep = ''))
        })
    })

})
