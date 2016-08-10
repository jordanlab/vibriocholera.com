library('shiny')
library('genoPlotR')
t6pred <- function(outdir){
  annotate <- paste('/home/blast/prediction_server/server/predict_t6.pl -fasta ',outdir,'/prots.faa', sep="")
  system(annotate)
}
plotGff <- function(outdir){
  genplotdata <-
    paste(
      '/home/blast/prediction_server/server/gff.pl -fasta '
      ,outdir,'/prots.faa', sep="")
  system(genplotdata)
  filelist = dir(outdir, pattern = "*.ptt")
  dna <- list()
  annot <- list()
  for (i in 1:length(filelist)) {
    file = paste(outdir, filelist[i], sep = "")
    dna[[i]] <- read_dna_seg_from_ptt(file)
    mid_pos <- middle(dna[[i]])
    annot[[i]] <-
      annotation(x1 = mid_pos,
                 text = dna[[i]]$name,
                 rot = "45")
  }
  
  plot_gene_map(dna_segs = dna, annotations = annot)
}
plotSkel <- function(){
  filelist = dir('/home/blast/prediction_server/server/skel/', pattern = "*.ptt")
  dna <- list()
  annot <- list()
  for (i in 1:length(filelist)) {
    file = paste('/home/blast/prediction_server/server/skel/', filelist[i], sep = "")
    dna[[i]] <- read_dna_seg_from_ptt(file)
    mid_pos <- middle(dna[[i]])
    annot[[i]] <-
      annotation(x1 = mid_pos,
                 text = dna[[i]]$name,
                 rot = "45")
  }
  
  plot_gene_map(dna_segs = dna, annotations = annot)
}


shinyServer(function(input, output, session) {
  observe({
    # switch tab
    if (!is.null(input$fastaFile)  & input$submit != 0L) {
      updateTabsetPanel(session, "vchot6ss", selected = "T6SS Predictions")
    }
  })
  
  
  annotPlot = reactive({
    seqCond = input$seqType
    predCond = input$predict
    if (seqCond == 'protein' &
        predCond == 'pred' &
        input$submit != 0L & !is.null(input$fastaFile)) {

      outdir = substr(input$fastaFile$datapath, 1, nchar(input$fastaFile$datapath) - 1)
      cpProt <- paste('cp ', input$fastaFile$datapath, ' ', outdir, '/prots.faa', sep="")
      system(cpProt)
      sed <- paste("sed -i 's/ #.*.//g' ",outdir,"/prots.faa", sep="")
      system(sed)
      t6pred(outdir)
      plotSkel()
    }
    else if (predCond == 'pred' &
             seqCond == 'genomic' &
             input$submit != 0L & !is.null(input$fastaFile)) {
      outdir = substr(input$fastaFile$datapath, 1, nchar(input$fastaFile$datapath) - 1)
      prediction <- paste('prodigal -q -c -i ', input$fastaFile$datapath, ' -a ', outdir, '/prots.faa -f gff -o ',outdir,'/prots.gff 2> /dev/null', sep="")
      system(prediction)
      sed <- paste("sed -i 's/ #.*.//g' ",outdir,"/prots.faa", sep="")
      system(sed)
      t6pred(outdir)
      plotGff(outdir)
    }
    else if (predCond == 'nopred' &
             seqCond == 'genomic' &
             input$submit != 0L & !is.null(input$fastaFile) & !is.null(input$gffFile)) {
      outdir = substr(input$fastaFile$datapath, 1, nchar(input$fastaFile$datapath) - 1)
      prediction <- paste('prodigal -q -c -i ', input$fastaFile$datapath, ' -a ', outdir, '/prots.faa -f gff -o ',outdir,'/prots.gff 2> /dev/null', sep="")
    }
    
    # else if (predCond == 'nopred' &
    #          input$submit != 0L &
    #          !is.null(input$fastaFile) & !is.null(input$gffFile)) {
    #   outdir = substr(input$fastaFile$datapath,
    #                   1,
    #                   nchar(input$fastaFile$datapath) - 1)
    #   annotate <-
    #     paste(
    #       '/home/blast/prediction_server/server/predict_t6.pl -predict no -fasta ',
    #       input$fastaFile$datapath,
    #       ' -gff ',
    #       input$gffFile$datapath
    #     )
    #   #system('echo "',input$fastaFile$datapath,'\n$(date) >> ~/vibrio_project/webtest.txt"', intern = TRUE)
    #   system(annotate)
    #   genplotdata <-
    #     paste(
    #       '/home/blast/prediction_server/server/gff.pl -fasta ',
    #       input$fastaFile$datapath
    #     )
    #   system(genplotdata)
    #   filelist = dir(outdir, pattern = "*.ptt")
    #   dna <- list()
    #   annot <- list()
    #   for (i in 1:length(filelist)) {
    #     file = paste(outdir, filelist[i], sep = "")
    #     dna[[i]] <- read_dna_seg_from_ptt(file)
    #     mid_pos <- middle(dna[[i]])
    #     annot[[i]] <-
    #       annotation(x1 = mid_pos,
    #                  text = dna[[i]]$name,
    #                  rot = "45")
    #   }
    #   
    #   plot_gene_map(dna_segs = dna, annotations = annot)
    # }
    else{
      return(NULL)
    }
  })
  
  output$plot = renderPlot({
    annotPlot()
  })
  output$buttons <- renderUI({
    downloadButton('dlFasta', 'Download Fasta Amino Acid',
                   class = 'btn btn-primary btn-large')
    downloadButton('dlPreds', 'Download predicted T6SS proteins',
                   class = 'btn btn-primary btn-large')
  })
  output$dlFasta = downloadHandler(
    filename = function() {
      paste("proteins",
            paste(collapse = '-'),
            '-',
            gsub(' ', '-', gsub(':', '-', Sys.time())),
            '.faa',
            sep = '')
    },
    content = function(file) {
      outdir = substr(input$fastaFile$datapath, 1, nchar(input$fastaFile$datapath) -
                        1)
      copyFile = paste(outdir, "prots.faa", sep = "")
      file.copy(copyFile, file)
    }
  )
  output$dlPreds = downloadHandler(
    filename = function() {
      paste("predictions",
            paste(collapse = '-'),
            '-',
            gsub(' ', '-', gsub(':', '-', Sys.time())),
            '.faa',
            sep = '')
    },
    content = function(file) {
      outdir = substr(input$fastaFile$datapath, 1, nchar(input$fastaFile$datapath) -
                        1)
      copyFile = paste(outdir, "predictions.faa", sep = "")
      file.copy(copyFile, file)
    }
  )
  
})
