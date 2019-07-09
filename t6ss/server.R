library('shiny')
library('genoPlotR')
library(uuid)
library(shinyjs)

options(shiny.maxRequestSize=8*1024^2, shiny.reactlog=FALSE)
# t6pred <- function(outdir){
#   annotate <- paste('/home/blast/prediction_server/server/predict_t6.pl -fasta ',outdir,'/prots.faa', sep="")
#   system(annotate)
# }
# plotGff <- function(outdir){
#   genplotdata <-
#     paste(
#       '/home/blast/prediction_server/server/gff.pl -fasta '
#       ,outdir,'/prots.faa', sep="")
#   system(genplotdata)
#   filelist = dir(outdir, pattern = "*.ptt")
#   dna <- list()
#   annot <- list()
#   for (i in 1:length(filelist)) {
#     file = paste(outdir,"/", filelist[i], sep = "")
#     dna[[i]] <- read_dna_seg_from_ptt(file)
#     mid_pos <- middle(dna[[i]])
#     annot[[i]] <-
#       annotation(x1 = mid_pos,
#                  text = dna[[i]]$name,
#                  rot = "45")
#   }

#   plot_gene_map(dna_segs = dna, annotations = annot)
# }
# plotSkel <- function(outdir){
#   genplotdata <-
#     paste(
#       '/home/blast/prediction_server/server/gff.pl -fasta '
#       ,outdir,'/prots.faa', sep="")
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

#   plot_gene_map(dna_segs = dna, annotations = annot)
# }


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
      # onclick("submit", runjs(paste0('"window.open("',"'", "https://t6ss.vibriocholera.com/results/",id,"/index.html')")))
      # onclick("submit", runjs(window.open('test.txt')))

      output$link <- renderUI({
        HTML(paste("","","<h3>Your submission has been added to the queue!</h3><br/>",
          "<h4>" ,outText,
          "</h4><br/>","<h3><a href='https://t6ss.vibriocholera.com/results/",id,"/index.html' target=_blank>Click here for your results</a></h3><br/><p>Results will be available for at least 7 days</p>", sep = ''))
        # a("", href=, target="_blank"),
        })
      # prodigal = paste0("prodigal -f gff -i ", fasta, " -a ",
      #   paste0(fasta,".faa"), " -d ", paste0(fasta,".fna"), ' -o ', paste0(fasta,".gff"))
      # cat(prodigal)
      # system(prodigal)

      # for (hmm in c('PAAR','aux3','lipase','vgrg','hcp','unknown','hydrolase','lysm','ntpase','transferase')){
      #   run_hmm = paste0(
      #   'hmmscan --cpu 4 -o /dev/null --tblout ', paste0(out_path,".",hmm), ' /home/blast/prediction_server/server/hmm_profiles/',hmm,'.hmm ', paste0(fasta,".faa"))
      #   cat(run_hmm)
      #   cat("\n")
      #   system(run_hmm)
      # }
      # combine = paste0('cat ', out_path,'.* >', fasta,".output")
      # cat(paste(combine, "\n"))
      # system(combine)
      # parse = paste0("/home/blast/prediction_server/server/parse_predictions.py -i ", fasta,".output -a ", fasta,".faa -n ", fasta,".fna -g ", fasta,".gff -f",fasta, ' -o ', sub_path)
      # cat(paste(parse, "\n"))
      # system(parse)
      # output_results = paste0("cat /home/blast/prediction_server/server/done.html ", sub_path, "/*.html > ",dataFolder,"/index.html" )
      # cat(paste(output_results, "\n"))
      # system(output_results)
      # file.rename(file.path(sub_path,"nucleotides.fna"), file.path(dataFolder,"nucleotides.fna") )
      # file.rename(file.path(sub_path,"proteins.faa"), file.path(dataFolder,"proteins.faa"))

    })
# observe({
#     print(input$vchot6ss)
#   })
#     observe({
#       if(req(input$vchot6ss == "results")){
#               sub_path = file.path('/home/blast/prediction_server/server/submissions',id)
#       out_path = file.path('/home/blast/prediction_server/server/submissions',id,"output")
#       dir.create(sub_path)
#       filename = input$fastaFile$name
#       fasta = file.path(sub_path,filename)
#       file.copy(input$fastaFile$datapath, fasta)
#       prodigal = paste0("prodigal -f gff -i ", fasta, " -a ",
#         paste0(fasta,".faa"), " -d ", paste0(fasta,".fna"), ' -o ', paste0(fasta,".gff"))
#       cat(prodigal)
#       system(prodigal)

#       for (hmm in c('PAAR','aux3','lipase','vgrg','hcp','unknown','hydrolase','lysm','ntpase','transferase')){
#         run_hmm = paste0(
#         'hmmscan --cpu 4 -o /dev/null --tblout ', paste0(out_path,".",hmm), ' /home/blast/prediction_server/server/hmm_profiles/',hmm,'.hmm ', paste0(fasta,".faa"))
#         cat(run_hmm)
#         cat("\n")
#         system(run_hmm)
#       }
#       combine = paste0('cat ', out_path,'.* >', fasta,".output")
#       cat(paste(combine, "\n"))
#       system(combine)
#       parse = paste0("/home/blast/prediction_server/server/parse_predictions.py -i ", fasta,".output -a ", fasta,".faa -n ", fasta,".fna -g ", fasta,".gff -f",fasta, ' -o ', sub_path)
#       cat(paste(parse, "\n"))
#       system(parse)
#       output_results = paste0("cat /home/blast/prediction_server/server/done.html ", sub_path, "/*.html > ",dataFolder,"/index.html" )
#       cat(paste(output_results, "\n"))
#       system(output_results)
#       file.rename(file.path(sub_path,"nucleotides.fna"), file.path(dataFolder,"nucleotides.fna") )
#       file.rename(file.path(sub_path,"proteins.faa"), file.path(dataFolder,"proteins.faa"))
#     }
#    })


})

 #  annotPlot = reactive({
 # #    seqCond = input$seqType
 # #    predCond = input$predict
 # #    if (seqCond == 'protein' &
 # #        predCond == 'pred' &
 # #        input$submit != 0L & !is.null(input$fastaFile)) {
 # #      outdir = dirname(input$fastaFile$datapath)
 # #      cpProt <- paste('cp ', input$fastaFile$datapath, ' ', outdir, '/prots.faa', sep="")
	# # cat(cpProt)
 # #      system(cpProt)
 # #      sed <- paste("sed -i 's/ #.*.//g' ",outdir,"/prots.faa", sep="")
 # #      system(sed)
 # #      t6pred(outdir)
 # #    }
 # #    else if (predCond == 'pred' &
 # #             seqCond == 'genomic' &
 # #             input$submit != 0L & !is.null(input$fastaFile)) {
 # #      outdir = dirname(input$fastaFile$datapath)
 # #      prediction <- paste('prodigal -q -c -i ', input$fastaFile$datapath, ' -a ', outdir, '/prots.faa -f gff -o ',outdir,'/prots.gff 2> /dev/null', sep="")
 # #      system(prediction)
 # #      sed <- paste("sed -i 's/ #.*.//g' ",outdir,"/prots.faa", sep="")
 # #      system(sed)
 # #      t6pred(outdir)
 # #      plotGff(outdir)
 # #    }
 # #    else if (predCond == 'nopred' &
 # #             seqCond == 'genomic' &
 # #             input$submit != 0L & !is.null(input$fastaFile) & !is.null(input$gffFile)) {
 # #      outdir = dirname(input$fastaFile$datapath)
 # #      prediction <- paste('prodigal -q -c -i ', input$fastaFile$datapath, ' -a ', outdir, '/prots.faa -f gff -o ',outdir,'/prots.gff 2> /dev/null', sep="")
 # #    }

 #    # else if (predCond == 'nopred' &
 #    #          input$submit != 0L &
 #    #          !is.null(input$fastaFile) & !is.null(input$gffFile)) {
 #    #   outdir = substr(input$fastaFile$datapath,
 #    #                   1,
 #    #                   nchar(input$fastaFile$datapath) - 1)
 #    #   annotate <-
 #    #     paste(
 #    #       '/home/blast/prediction_server/server/predict_t6.pl -predict no -fasta ',
 #    #       input$fastaFile$datapath,
 #    #       ' -gff ',
 #    #       input$gffFile$datapath
 #    #     )
 #    #   #system('echo "',input$fastaFile$datapath,'\n$(date) >> ~/vibrio_project/webtest.txt"', intern = TRUE)
 #    #   system(annotate)
 #    #   genplotdata <-
 #    #     paste(
 #    #       '/home/blast/prediction_server/server/gff.pl -fasta ',
 #    #       input$fastaFile$datapath
 #    #     )
 #    #   system(genplotdata)
 #    #   filelist = dir(outdir, pattern = "*.ptt")
 #    #   dna <- list()
 #    #   annot <- list()
 #    #   for (i in 1:length(filelist)) {
 #    #     file = paste(outdir, filelist[i], sep = "")
 #    #     dna[[i]] <- read_dna_seg_from_ptt(file)
 #    #     mid_pos <- middle(dna[[i]])
 #    #     annot[[i]] <-
 #    #       annotation(x1 = mid_pos,
 #    #                  text = dna[[i]]$name,
 #    #                  rot = "45")
 #    #   }
 #    #
 #    #   plot_gene_map(dna_segs = dna, annotations = annot)
 #    # }
 #    else{
 #      return(NULL)
 #    }
 #  })


# })
