footer <- function(){
  tags$div(
    class = "footer",
    style = "text-align:center",
    tags$div(
      class = "foot-inner",
      list(
        hr(),
        "pcaExplorer is a project developed by Federico Marini in the Bioinformatics division of the ",
        tags$a(href="http://www.unimedizin-mainz.de/imbei","IMBEI"),
        ". ",br(),
        "Development of the pcaExplorer package is on ",
        tags$a(href="https://github.com/federicomarini/pcaExplorer", "GitHub")
      )
    )
  )
}



if ( !requireNamespace('shiny',quietly = TRUE) ) {
  stop("pcaExplorer requires 'shiny'. Please install it using
       install.packages('shiny')")
}

library("shinyBS") # to correctly display the bsTooltips

# get modes and themes for the ace editor
modes <- shinyAce::getAceModes()
themes <- shinyAce::getAceThemes()

## upload max 300mb files - can be changed if necessary
options(shiny.maxRequestSize=300*1024^2)













shinyServer(function(input, output, session) {

  ## for server setup:
  dds=NULL
  rlt=NULL
  countmatrix=NULL
  coldata=NULL
  pca2go=NULL
  annotation=NULL








  ## placeholder for the figures to export
  exportPlots <- reactiveValues(
    samplessamples_heatmap=NULL,
    reads_barplot=NULL,

    samplesPca=NULL,
    samplesZoom=NULL,
    samplesScree=NULL,
    samplesHiload=NULL,
    samplesOutlier=NULL,

    genesPca=NULL,
    genesZoom=NULL,
    genesProfile=NULL,
    genesBoxplot=NULL,
    genesHeatmap=NULL,

    genefinder_countsplot=NULL
  )

  if(!is.null(dds)){
    if(is.null(sizeFactors(dds)))
      dds <- estimateSizeFactors(dds)
  }

  ## reactive values to use in the app
  values <- reactiveValues()
  values$mydds <- dds
  values$myrlt <- rlt
  values$mycountmatrix <- countmatrix
  values$mymetadata <- coldata
  values$mypca2go <- pca2go
  values$myannotation <- annotation

  user_settings <- reactiveValues(save_width = 15, save_height = 11)


  if(!is.null(dds)){
    if(!is(dds,"DESeqDataSet"))
      stop("dds must be a DESeqDataSet object. If it is a simple counts matrix, provide it to the countmatrix parameter!")

    if(is.null(sizeFactors(dds)))
      dds <- estimateSizeFactors(dds)
  }
  if(!is.null(rlt)){
    if(!is(rlt,"DESeqTransform"))
      stop("dds must be a DESeqTransform object")
  }

  # compute only rlt if dds is provided but not cm&coldata
  if(!is.null(dds) & (is.null(countmatrix) & is.null(coldata)) & is.null(rlt)){
    withProgress(message = "computing rlog transformed values...",
                 value = 0,
                 {
                   values$myrlt <- rlogTransformation(dds)
                 })
  }

  output$color_by <- renderUI({
    if(is.null(values$mydds))
      return(NULL)
    poss_covars <- names(colData(values$mydds))
    selectInput('color_by', label = 'Group/color by: ',
                choices = c(NULL, poss_covars), selected = NULL,multiple = TRUE)
  })



  ## Render the UI element to upload the count matrix
  output$upload_count_matrix <- renderUI({
    if (!is.null(dds) | !is.null(countmatrix)) {
      NULL
    } else {
      return(fileInput(inputId = "uploadcmfile",
                       label = "Upload a count matrix file",
                       accept = c("text/csv", "text/comma-separated-values",
                                  "text/tab-separated-values", "text/plain",
                                  ".csv", ".tsv"), multiple = FALSE))
    }
  })

  readCountmatrix <- reactive({
    if (is.null(input$uploadcmfile))
      return(NULL)
    cm <- utils::read.delim(input$uploadcmfile$datapath, header = TRUE,
                            as.is = TRUE, sep = "\t", quote = "",
                            row.names = 1, # https://github.com/federicomarini/pcaExplorer/issues/1
                            ## TODO: tell the user to use tsv, or use heuristics
                            ## to check what is most frequently occurring separation character? -> see sepGuesser.R
                            check.names = FALSE)

    return(cm)
  })


  output$upload_metadata <- renderUI({
    if (!is.null(dds) | !is.null(coldata)) {
      NULL
    } else {
      return(fileInput(inputId = "uploadmetadatafile",
                       label = "Upload a sample metadata matrix file",
                       accept = c("text/csv", "text/comma-separated-values",
                                  "text/tab-separated-values", "text/plain",
                                  ".csv", ".tsv"), multiple = FALSE))
    }
  })

  readMetadata <- reactive({
    if (is.null(input$uploadmetadatafile))
      return(NULL)
    coldata <- utils::read.delim(input$uploadmetadatafile$datapath, header = TRUE,
                                 as.is = TRUE, sep = "\t", quote = "",
                                 check.names = FALSE)

    return(coldata)
  })


  output$upload_annotation <- renderUI({
    if (!is.null(annotation)) {
      NULL
    } else {
      return(fileInput(inputId = "uploadannotationfile",
                       label = "Upload an annotation file",
                       accept = c("text/csv", "text/comma-separated-values",
                                  "text/tab-separated-values", "text/plain",
                                  ".csv", ".tsv"), multiple = FALSE))
    }
  })

  readAnnotation <- reactive({
    if (is.null(input$uploadannotationfile))
      return(NULL)
    annodata <- utils::read.delim(input$uploadannotationfile$datapath, header = TRUE,
                                  as.is = TRUE, sep = "\t", quote = "",
                                  check.names = FALSE)

    return(annodata)
  })


  output$printdds <- renderPrint({

    shiny::validate(
      need(!is.null(values$mydds),
           "Upload your dataset, as a count matrix or passing it as a parameter, as well as the design information"
      )
    )

    values$mydds

  })


  output$printanno <- DT::renderDataTable({
    shiny::validate(
      need(!is.null(values$myannotation),
           "Upload your annotation table as a matrix/data frame or passing it as a parameter"
      )
    )
    DT::datatable(values$myannotation,options = list(pageLength=10))

  })


  createDDS <- reactive({
    if(is.null(countmatrix) | is.null(coldata))
      return(NULL)

    dds <- DESeqDataSetFromMatrix(countData = countmatrix,
                                  colData = coldata,
                                  design=~1)
    dds <- estimateSizeFactors(dds)

    return(dds)


  })

  createRLT <- reactive({
    if(is.null(countmatrix) | is.null(coldata))
      return(NULL)

    rlt <- rlogTransformation(values$mydds)

    return(rlt)
  })


  observeEvent(createDDS,
               {
                 if(!is.null(values$mycountmatrix) & !is.null(values$mymetadata))
                   values$mydds <- createDDS()
               })

  observeEvent(createRLT,
               {
                 if(!is.null(values$mycountmatrix) & !is.null(values$mymetadata))
                   values$myrlt <- createRLT()
               })

  # useful when count matrix is uploaded by hand
  sneakpeek <- reactiveValues()
  observeEvent(input$uploadcmfile,
               {

                 sneakpeek$cm <- readCountmatrix()
               })

  output$sneakpeekcm <- DT::renderDataTable({
    head(sneakpeek$cm,10)
  })




  # as in http://stackoverflow.com/questions/29716868/r-shiny-how-to-get-an-reactive-data-frame-updated-each-time-pressing-an-actionb
  observeEvent(input$uploadcmfile,
               {
                 values$mycountmatrix <- readCountmatrix()
                 if(!is.null(values$mymetadata)){
                   withProgress(message="Computing the objects...",value = 0,{

                     values$mydds <- DESeqDataSetFromMatrix(countData = values$mycountmatrix,
                                                            colData = values$mymetadata,
                                                            design=~1)
                     values$myrlt <- rlogTransformation(values$mydds)})
                 }
               })

  observeEvent(input$uploadmetadatafile,
               {
                 values$mymetadata <- readMetadata()
                 if(!is.null(values$mycountmatrix)){
                   withProgress(message="Computing the objects...",value = 0,{

                     values$mydds <- DESeqDataSetFromMatrix(countData = values$mycountmatrix,
                                                            colData = values$mymetadata,
                                                            design=~1)
                     values$myrlt <- rlogTransformation(values$mydds)})
                 }
               })

  observeEvent(input$uploadannotationfile,
               {
                 values$myannotation <- readAnnotation()
               })



  output$showuploaded1 <- renderPrint({
    head(values$mycountmatrix)
  })
  output$showuploaded2 <- renderPrint({
    values$mymetadata
  })
  output$showuploaded3 <- renderPrint({
    values$mydds
  })
  output$showuploaded4 <- renderPrint({
    values$myrlt
  })


  colSel <- reactive({
    # find out how many colors to generate: if no factor is selected, either
    # return all say steelblue or all different

    if(!is.null(input$color_by)) {
      expgroups <- as.data.frame(colData(values$myrlt)[,input$color_by])
      expgroups <- interaction(expgroups)
    } else {
      expgroups <- factor(colnames(values$myrlt))
      # return(rep("steelblue",ncol(values$myrlt))) # to return all same
    }

    nrgroups <- length(levels(expgroups))

    if(input$col_palette=="hue"){
      return(hue_pal()(nrgroups))
    }
    # hue_pal()(ncol(values$myrlt)/2) # or somewhat other way
    if(input$col_palette=="set1"){
      if(nrgroups <= 9) { # max color nr allowed for set1
        return(brewer_pal(palette = "Set1")(nrgroups))
      } else {
        return(hue_pal()(nrgroups)) # plus print message?
      }
    }
    # (ncol(values$myrlt)/2) # or somewhat other way
    if(input$col_palette=="rainbow"){
      return(rainbow(nrgroups))
    }
  })


  output$sessioninfo <- renderPrint({
    sessionInfo()
  })



  output$showdata <- renderPrint({
    values$mydds
  })

  output$showcoldata <- DT::renderDataTable({
    totreads <- (colSums(counts(values$mydds)))
    df <- data.frame(
      colData(values$mydds),
      "Total number of reads"=totreads
    )

    datatable(df)
  })



  current_countmat <- reactive({
    if(input$countstable_unit=="raw_counts")
      return(counts(values$mydds,normalized=FALSE))
    if(input$countstable_unit=="normalized_counts")
      return(counts(values$mydds,normalized=TRUE))
    if(input$countstable_unit=="rlog_counts")
      return(assay(values$myrlt))
    if(input$countstable_unit=="log10_counts")
      return(log10(1 + counts(values$mydds,normalized=TRUE)))
    if(input$countstable_unit=="tpm_counts")
      return(NULL) ## TODO!: assumes length of genes/exons as known, and is currently not required in the dds

  })

  output$showcountmat <- DT::renderDataTable({
    datatable(current_countmat())
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(input$countstable_unit,"table.csv")
    },
    content = function(file) {
      write.csv(current_countmat(), file)
    }
  )

  output$download_genefinder_countstable <- downloadHandler(
    filename = function() {
      paste0("genefinder_pcaE_","table.csv")
    },
    content = function(file) {

      anno_id <- rownames(values$myrlt)
      anno_gene <- values$myannotation$gene_name

      if(is.null(input$color_by) & input$genefinder!="")
        return(NULL)
      if(is.null(input$color_by) & input$genefinder=="")
        return(NULL)
      if(input$genefinder=="")
        return(NULL)
      if(!input$genefinder %in% anno_gene & !input$genefinder %in% anno_id)
        return(NULL)

      if (input$genefinder %in% anno_id) {
        selectedGene <- rownames(values$myrlt)[match(input$genefinder,rownames(values$myrlt))]
        selectedGeneSymbol <- values$myannotation$gene_name[match(selectedGene,rownames(values$myannotation))]
      }
      if (input$genefinder %in% anno_gene) {
        selectedGeneSymbol <- values$myannotation$gene_name[which(values$myannotation$gene_name==input$genefinder)]
        if (length(selectedGeneSymbol) > 1) return(ggplot() + annotate("text",label=paste0("Type in a gene name/id of the following:\n",paste(selectedGene,collapse=", ")),0,0) + theme_bw())
        selectedGene <- rownames(values$myannotation)[which(values$myannotation$gene_name==input$genefinder)]
      }
      genedata <- plotCounts(values$mydds,gene=selectedGene,intgroup = input$color_by,returnData = TRUE)
      genedata

      write.csv(genedata, file)
    }
  )






  output$corrplot <- renderPlot({
    if(input$compute_pairwisecorr)
      pair_corr(current_countmat(),method=input$corr_method)
  })

  output$heatcorr <- renderPlot({
    if(input$compute_pairwisecorr)
      pheatmap(cor(current_countmat()))
  })


  output$pairwise_plotUI <- renderUI({
    if(!input$compute_pairwisecorr) return()

    plotOutput("corrplot", height = "1000px")
    # )
  })


  output$heatcorr_plotUI <- renderUI({
    if(!input$compute_pairwisecorr) return()

    plotOutput("heatcorr")
  })


  # overview on number of detected genes on different threshold types
  output$detected_genes <- renderPrint({
    t1 <- rowSums(counts(values$mydds))
    t2 <- rowMeans(counts(values$mydds,normalized=TRUE))

    thresh_rowsums <- input$threshold_rowsums
    thresh_rowmeans <- input$threshold_rowmeans
    abs_t1 <- sum(t1 > thresh_rowsums)
    rel_t1 <- 100 * mean(t1 > thresh_rowsums)
    abs_t2 <- sum(t2 > thresh_rowmeans)
    rel_t2 <- 100 * mean(t2 > thresh_rowmeans)

    cat("Number of detected genes:\n")
    # TODO: parametrize the thresholds
    cat(abs_t1,"genes have at least a sample with more than",thresh_rowsums,"counts\n")
    cat(paste0(round(rel_t1,3),"%"), "of the",nrow(values$mydds),
        "genes have at least a sample with more than",thresh_rowsums,"counts\n")
    cat(abs_t2,"genes have more than",thresh_rowmeans,"counts (normalized) on average\n")
    cat(paste0(round(rel_t2,3),"%"), "of the",nrow(values$mydds),
        "genes have more than",thresh_rowsums,"counts (normalized) on average\n")
    cat("Counts are ranging from", min(counts(values$mydds)),"to",max(counts(values$mydds)))
  })



  output$heatmapsampledist <- renderPlot({
    if (!is.null(input$color_by)){
      expgroups <- as.data.frame(colData(values$myrlt)[,input$color_by])
      # expgroups <- interaction(expgroups)
      rownames(expgroups) <- colnames(values$myrlt)
      colnames(expgroups) <- input$color_by

      pheatmap(as.matrix(dist(t(assay(values$myrlt)))),annotation_col = expgroups)
    } else {
      pheatmap(as.matrix(dist(t(assay(values$myrlt)))))
    }
  })



  output$reads_barplot <- renderPlot({
    rr <- colSums(counts(values$mydds))/1e6
    if(is.null(names(rr)))
      names(rr) <- paste0("sample_",1:length(rr))
    rrdf <- data.frame(Reads=rr,Sample=names(rr),stringsAsFactors = FALSE)
    if (!is.null(input$color_by)) {
      selGroups <- as.data.frame(colData(values$mydds)[input$color_by])
      rrdf$Group <- interaction(selGroups)
      p <- ggplot(rrdf,aes_string("Sample",weight="Reads")) + geom_bar(aes_string(fill="Group")) + theme_bw()
      p
    } else {
      p <- ggplot(rrdf,aes_string("Sample",weight="Reads")) + geom_bar() + theme_bw()

      exportPlots$reads_barplot <- p
      p
    }
  })

  output$reads_summary <- renderPrint({
    print(colSums(counts(values$mydds)))
    summary(colSums(counts(values$mydds))/1e6)
  })








  ## SAMPLES VIEW
  output$samples_pca <- renderPlot({
    res <- pcaplot(values$myrlt,intgroup = input$color_by,ntop = input$pca_nrgenes,
                   pcX = as.integer(input$pc_x),pcY = as.integer(input$pc_y),
                   text_labels = input$sample_labels,
                   point_size = input$pca_point_size, title="Samples PCA",
                   ellipse = input$pca_ellipse, ellipse.prob = input$pca_cislider)


    res <- res + theme_bw()
    exportPlots$samplesPca <- res
    res
  })

  output$samples_pca_zoom <- renderPlot({

    shiny::validate(
      need(!is.null(input$pca_brush),
           "Zoom in by brushing in the main plot panel above"
      )
    )
    # if(is.null(input$pca_brush))
    # return(ggplot() + annotate("text",label="zoom in by brushing",0,0) + theme_bw())

    res <- pcaplot(values$myrlt,intgroup = input$color_by,ntop = input$pca_nrgenes,
                   pcX = as.integer(input$pc_x),pcY = as.integer(input$pc_y),
                   text_labels = input$sample_labels,
                   point_size = input$pca_point_size, title="Samples PCA - zoom in",
                   ellipse = input$pca_ellipse, ellipse.prob = input$pca_cislider
    )
    res <- res + xlim(input$pca_brush$xmin,input$pca_brush$xmax) + ylim(input$pca_brush$ymin,input$pca_brush$ymax)
    res <- res + theme_bw()
    exportPlots$samplesZoom <- res
    res
  })

  output$samples_scree <- renderPlot({
    rv <- rowVars(assay(values$myrlt))
    select <- order(rv, decreasing = TRUE)[seq_len(min(input$pca_nrgenes,length(rv)))]
    pca <- prcomp(t(assay(values$myrlt)[select, ]))

    res <- pcascree(pca,type = input$scree_type, pc_nr = input$scree_pcnr, title="Scree plot for the samples PCA")
    res <- res + theme_bw()
    exportPlots$samplesScree <- res
    res
  })


  output$geneshiload <- renderPlot({
    rv <- rowVars(assay(values$myrlt))
    select <- order(rv, decreasing = TRUE)[seq_len(min(input$pca_nrgenes,length(rv)))]
    pca <- prcomp(t(assay(values$myrlt)[select, ]))

    par(mfrow=c(2,1))
    hi_loadings(pca,whichpc = as.integer(input$pc_x),topN = input$ntophiload,annotation = values$myannotation)
    hi_loadings(pca,whichpc = as.integer(input$pc_y),topN = input$ntophiload,annotation = values$myannotation)

  })


  output$ui_outliersamples <- renderUI({
    available_samples <- c("",colnames(values$myrlt))

    selectInput("outlierselection",label = "Select which sample(s) to remove - suspected outliers",choices = available_samples,multiple = TRUE)

  })

  output$samples_outliersremoved <- renderPlot({

    shiny::validate(
      need(input$outlierselection!="",
           message = "Select at least one sample to plot the new PCA where the selection is removed")
    )

    currentrlt <- values$myrlt
    allsamples <- colnames(currentrlt)

    outliersamples <- input$outlierselection
    currentrlt <- currentrlt[,setdiff(allsamples,outliersamples)]

    res <- pcaplot(currentrlt,intgroup = input$color_by,ntop = input$pca_nrgenes,
                   pcX = as.integer(input$pc_x),pcY = as.integer(input$pc_y),
                   text_labels = input$sample_labels,
                   point_size = input$pca_point_size, title="Samples PCA",
                   ellipse = input$pca_ellipse, ellipse.prob = input$pca_cislider
    )
    res <- res + theme_bw()
    # exportPlots$samplesPca <- res
    exportPlots$samplesOutlier <- res
    res

  })


  output$pca3d <- renderScatterplotThree({
    pcaplot3d(values$myrlt,intgroup = input$color_by,ntop = input$pca_nrgenes,
              pcX = as.integer(input$pc_x),pcY = as.integer(input$pc_y),pcZ = as.integer(input$pc_z))

  })




  ## GENES VIEW
  output$genes_biplot <- renderPlot({
    if(!is.null(input$color_by)) {
      expgroups <- as.data.frame(colData(values$myrlt)[,input$color_by])
      expgroups <- interaction(expgroups)
      expgroups <- factor(expgroups,levels=unique(expgroups))

    } else {
      expgroups <- colnames(values$myrlt)
    }
    colGroups <- colSel()[factor(expgroups)]

    res <- genespca(values$myrlt,
                    ntop = input$pca_nrgenes,
                    choices = c(as.integer(input$pc_x),as.integer(input$pc_y)),
                    biplot = TRUE,
                    arrowColors = factor(colGroups,levels=unique(colGroups)),
                    groupNames = expgroups,
                    alpha=input$pca_point_alpha,coordEqual=FALSE,useRownamesAsLabels=FALSE,labels.size=input$pca_label_size,
                    point_size=input$pca_point_size,varname.size=input$pca_varname_size, scaleArrow = input$pca_scale_arrow,annotation=values$myannotation)
    exportPlots$genesPca <- res
    res
  })


  output$genes_biplot_zoom <- renderPlot({
    # if(is.null(input$pcagenes_brush)) return(ggplot() + annotate("text",label="zoom in by brushing",0,0) + theme_bw())

    shiny::validate(
      need(
        !is.null(input$pcagenes_brush),
        "Zoom in by brushing in the main panel - this will also allow displaying the gene names"
      )
    )

    if(!is.null(input$color_by)) {
      expgroups <- as.data.frame(colData(values$myrlt)[,input$color_by])
      expgroups <- interaction(expgroups)
      expgroups <- factor(expgroups,levels=unique(expgroups))
    } else {
      expgroups <- colnames(values$myrlt)
    }
    colGroups <- colSel()[factor(expgroups)]

    res <- genespca(values$myrlt,
                    ntop = input$pca_nrgenes,
                    choices = c(as.integer(input$pc_x),as.integer(input$pc_y)),
                    biplot = TRUE,
                    arrowColors = factor(colGroups,levels=unique(colGroups)),
                    groupNames = expgroups,
                    alpha=input$pca_point_alpha,coordEqual=FALSE,
                    var.axes=input$variable_labels, # workaround for a ggplot2 bug/missing thing: here details: https://github.com/hadley/ggplot2/issues/905
                    labels.size=input$pca_label_size,varname.size=input$pca_varname_size,
                    scaleArrow = input$pca_scale_arrow,point_size=input$pca_point_size,annotation=values$myannotation)

    res <- res +
      xlim(input$pcagenes_brush$xmin,input$pcagenes_brush$xmax) +
      ylim(input$pcagenes_brush$ymin,input$pcagenes_brush$ymax)
    exportPlots$genesZoom <- res
    res
  })


  output$genes_profileexplorer <- renderPlot({
    shiny::validate(
      need(
        !is.null(input$pcagenes_brush),
        "Zoom in by brushing in the main panel"
      )
    )
    shiny::validate(
      need(
        length(input$color_by)>0,
        "Select an experimental factor in the Group/color by element in the sidebar"
      )
    )

    geneprofiler(values$myrlt,
                 genelist = curData_brush()$ids,
                 intgroup = input$color_by,
                 plotZ = input$zprofile)

  })


  output$genes_biplot_boxplot <- renderPlot({
    # if(length(input$color_by)==0) return(ggplot() + annotate("text",label="select an experimental factor",0,0) + theme_bw())
    # if(is.null(input$pcagenes_zoom_click)) return(ggplot() + annotate("text",label="click to generate the boxplot\nfor the selected gene",0,0) + theme_bw())

    shiny::validate(
      need(
        length(input$color_by)>0,
        "Select an experimental factor in the Group/color by element in the sidebar"
      )
    )
    shiny::validate(
      need(
        !is.null(input$pcagenes_zoom_click),
        "Click the plot above to generate the boxplot for the selected gene"
      )
    )

    selectedGene <- curData_zoomClick()$ids
    selectedGeneSymbol <- values$myannotation$gene_name[match(selectedGene,rownames(values$myannotation))]
    # plotCounts(dds_cleaner,)

    shiny::validate(
      need(nrow(curData_zoomClick()) >0,message = "Click closer to a gene to get the boxplot")

    )

    genedata <- plotCounts(values$mydds,gene=selectedGene,intgroup = input$color_by,returnData = TRUE)

    onlyfactors <- genedata[,match(input$color_by,colnames(genedata))]
    genedata$plotby <- interaction(onlyfactors)

    genedata$sampleID <- rownames(genedata)

    # input$plot_style chooses the style of plotting
    if(input$plot_style=="boxplot"){
      res <- ggplot(genedata,aes_string(x="plotby",y="count",fill="plotby")) +
        geom_boxplot(outlier.shape = NA,alpha=0.7) + theme_bw()
      if(input$ylimZero_genes){
        res <- res + scale_y_log10(name="Normalized counts - log10 scale",limits=c(0.4,NA))
      } else {
        res <- res + scale_y_log10(name="Normalized counts - log10 scale")
      }



      res <- res +
        labs(title=paste0("Normalized counts for ",selectedGeneSymbol," - ",selectedGene)) +
        scale_x_discrete(name="") +
        geom_jitter(aes_string(x="plotby",y="count"),position = position_jitter(width = 0.1)) +
        scale_fill_discrete(name="Experimental\nconditions")
      # if(input$addsamplelabels){
      #   res <- res + geom_text(aes(label=sampleID),hjust=-.1,vjust=0)
      # }

      exportPlots$genesBoxplot <- res
      res
    } else if(input$plot_style=="violin plot"){
      res <- ggplot(genedata,aes_string(x="plotby",y="count",fill="plotby")) +
        geom_violin(aes_string(col="plotby"),alpha = 0.6) + theme_bw()
      if(input$ylimZero_genes){
        res <- res + scale_y_log10(name="Normalized counts - log10 scale",limits=c(0.4,NA))
      } else {
        res <- res + scale_y_log10(name="Normalized counts - log10 scale")
      }


      res <- res +
        labs(title=paste0("Normalized counts for ",selectedGeneSymbol," - ",selectedGene)) +
        scale_x_discrete(name="") +
        geom_jitter(aes_string(x="plotby",y="count"),alpha = 0.8,position = position_jitter(width = 0.1)) +
        scale_fill_discrete(name="Experimental\nconditions") + scale_color_discrete(guide="none")
      # if(input$addsamplelabels){
      #   res <- res + geom_text(aes(label=sampleID),hjust=-.1,vjust=0)
      # }

      exportPlots$genesBoxplot <- res
      res
    }
  })


  # for reading in the brushed/clicked points
  curData_brush <- reactive({
    df2 <- genespca(values$myrlt,
                    ntop = input$pca_nrgenes,
                    choices = c(as.integer(input$pc_x),as.integer(input$pc_y)),
                    biplot = TRUE,
                    # arrowColors = colGroups,
                    alpha=input$pca_point_alpha,
                    returnData=TRUE,annotation=values$myannotation)
    df2$geneName <- values$myannotation$gene_name[match(rownames(df2),rownames(values$myannotation))]
    res <- brushedPoints(df2, input$pcagenes_brush,xvar="xvar",yvar="yvar",)
    res
  })


  curData_click <- reactive({
    df2 <- genespca(values$myrlt,
                    ntop = input$pca_nrgenes,
                    choices = c(as.integer(input$pc_x),as.integer(input$pc_y)),
                    biplot = TRUE,
                    # arrowColors = colGroups,
                    alpha=input$pca_point_alpha,
                    returnData=TRUE,annotation=values$myannotation)
    df2$geneName <- values$myannotation$gene_name[match(rownames(df2),rownames(values$myannotation))]
    res <- nearPoints(df2, input$pcagenes_click,
                      threshold = 20, maxpoints = 3,
                      addDist = TRUE)
    # res <- brushedPoints(df2, input$pcagenes_brush,xvar="xvar",yvar="yvar",)
    res
  })


  # data to be used for plotting the picked gene from the zoomed panel
  curData_zoomClick <- reactive({
    df2 <- genespca(values$myrlt,
                    ntop = input$pca_nrgenes,
                    choices = c(as.integer(input$pc_x),as.integer(input$pc_y)),
                    biplot = TRUE,
                    # arrowColors = colGroups,
                    alpha=input$pca_point_alpha,
                    returnData=TRUE,annotation=values$myannotation)
    df2$geneName <- values$myannotation$gene_name[match(rownames(df2),rownames(values$myannotation))]
    res <- nearPoints(df2, input$pcagenes_zoom_click,
                      threshold = 20, maxpoints = 1,
                      addDist = TRUE)
    # res <- brushedPoints(df2, input$pcagenes_brush,xvar="xvar",yvar="yvar",)
    res
  })



  output$pca_brush_out <- DT::renderDataTable({
    datatable(curData_brush(),options = list(pageLength = 50))
  })

  output$pca_click_out <- DT::renderDataTable({
    datatable(curData_click(),options = list(pageLength = 50))
  })




  output$heatzoomd3 <- renderD3heatmap({
    shiny::validate(
      need(
        !is.null(input$pcagenes_brush),
        "Brush the main panel above to generate a heatmap"
      )
    )

    # if(is.null(input$pcagenes_brush)) return(NULL)

    brushedObject <- curData_brush()
    shiny::validate(
      need(
        nrow(brushedObject) > 1,
        "Brush to include at least two genes"
      )
    )

    selectedGenes <- brushedObject$ids
    toplot <- assay(values$myrlt)[selectedGenes,]
    rownames(toplot) <- values$myannotation$gene_name[match(rownames(toplot),rownames(values$myannotation))]

    mycolss <- c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#fee090","#fdae61","#f46d43","#d73027","#a50026") # to be consistent with red/blue usual coding

    d3heatmap(toplot,Colv = as.logical(input$heatmap_colv),colors = mycolss)
  })


  output$heatzoom <- renderPlot({
    # if(is.null(input$pcagenes_brush)) return(NULL)
    shiny::validate(
      need(
        !is.null(input$pcagenes_brush),
        "Brush the main panel above to generate a heatmap"
      )
    )

    brushedObject <- curData_brush()
    shiny::validate(
      need(
        nrow(brushedObject) > 1,
        "Brush to include at least two genes"
      )
    )
    selectedGenes <- brushedObject$ids
    toplot <- assay(values$myrlt)[selectedGenes,]
    rownames(toplot) <- values$myannotation$gene_name[match(rownames(toplot),rownames(values$myannotation))]
    # pheatmap(toplot,cluster_cols = as.logical(input$heatmap_colv))
    NMF::aheatmap(toplot,Colv = as.logical(input$heatmap_colv))
    ## aheatmap is actually consistent in displaying the clusters with most of other heatmap packages
    ## keep in mind: pheatmap does somehow a better job if scaling/centering
  })


  #     displayed by default, with possibility to select from gene id provided as row names of the objects
  #     output$ui_selectID <- renderUI({
  #       allIDs <- withProgress(message = "loading the names in the UI",value = 0,
  #                              {
  #                                rownames(values$myrlt)
  #                              }
  #
  #
  #       )
  #       selectInput("selectID",label = "Select ID",choices = c("",allIDs),selected=NULL)
  #     })
  #     # additionally displayed if an annotation is provided
  #     output$ui_selectName <- renderUI({
  #       shiny::validate(
  #         need(
  #           !is.null(values$myannotation),
  #           "If you provide an annotation table, you could search by the corresponding name/ID"
  #         )
  #       )
  #
  #       selectInput("selectName",label = "Select gene name",choices = c("",values$myannotation$gene_name),selected=NULL)
  #       # selectInput("selectName")
  #     })
  #
  #     output$debugf <- renderPrint({
  #       input$selectID
  #     })


  #     output$newgenefinder_plot <- renderPlot({
  #       if (input$selectID == "")
  #         return(ggplot() + annotate("text",label="Type in a gene name/id",0,0) + theme_bw())
  #
  #
  #       if(!is.null(values$myannotation)){
  #         if(input$selectName != "") {
  #           selectedGeneName <- input$selectName
  #           selectedGene <- rownames(values$myannotation)[which(values$myannotation$gene_name==input$selectName)]
  #         } else {
  #           if (input$selectID == "") {
  #             return(ggplot() + annotate("text",label="Type in a gene name/id",0,0) + theme_bw())
  #           } else {
  #             selectedGene <- input$selectID
  #             selectedGeneName <- ifelse(!is.null(values$myannotation),
  #                                        values$myannotation$gene_name[match(selectedGene,rownames(values$myannotation))],
  #                                        "")
  #           }
  #         }
  #       } else if (input$selectID != ""){
  #         selectedGene <- input$selectID
  #         selectedGeneName <- ifelse(!is.null(values$myannotation),
  #                                    values$myannotation$gene_name[match(selectedGene,rownames(values$myannotation))],
  #                                    "")
  #       } else {
  #         return(ggplot() + annotate("text",label="Type in a gene name/id",0,0) + theme_bw())
  #       }
  #
  #
  #       anno_id <- rownames(values$myannotation)
  #       anno_gene <- values$myannotation$gene_name
  #
  #       if(is.null(input$color_by))
  #         return(ggplot() + annotate("text",label="Select a factor to plot your gene",0,0) + theme_bw())
  # #       if(is.null(input$color_by) & (input$selectName=="" | input$selectID ==""))
  # #         return(ggplot() + annotate("text",label="Select a gene and a factor to plot gene",0,0) + theme_bw())
  # #       if((input$selectName=="" | input$selectID ==""))
  # #         return(ggplot() + annotate("text",label="Type in a gene name/id",0,0) + theme_bw())
  #       # if(!input$genefinder %in% anno_gene & !input$genefinder %in% anno_id)
  #         # return(ggplot() + annotate("text",label="Gene not found...",0,0) + theme_bw())
  #
  # #       if (input$genefinder %in% anno_id) {
  # #         selectedGene <- rownames(values$myrlt)[match(input$genefinder,rownames(values$myrlt))]
  # #         selectedGeneSymbol <- values$myannotation$gene_name[match(selectedGene,rownames(values$myannotation))]
  # #       }
  # #       if (input$genefinder %in% anno_gene) {
  # #         selectedGeneSymbol <- values$myannotation$gene_name[which(values$myannotation$gene_name==input$genefinder)]
  # #         if (length(selectedGeneSymbol) > 1) return(ggplot() + annotate("text",label=paste0("Type in a gene name/id of the following:\n",paste(selectedGene,collapse=", ")),0,0) + theme_bw())
  # #         selectedGene <- rownames(values$myannotation)[which(values$myannotation$gene_name==input$genefinder)]
  # #       }
  #
  #
  #
  #
  #
  #       genedata <- plotCounts(values$mydds,gene=selectedGene,intgroup = input$color_by,returnData = TRUE)
  #
  #       onlyfactors <- genedata[,match(input$color_by,colnames(genedata))]
  #       genedata$plotby <- interaction(onlyfactors)
  #
  #       p <- ggplot(genedata,aes(x=plotby,y=count,fill=plotby)) + geom_boxplot() + labs(title=paste0("Normalized counts for ",selectedGeneName," - ",selectedGene)) +  scale_x_discrete(name="") + geom_jitter(aes(x=plotby,y=count),position = position_jitter(width = 0.1)) + scale_fill_discrete(name="Experimental\nconditions")
  #
  #       if(input$ylimZero)
  #       {
  #         p <- p + scale_y_log10(name="Normalized counts - log10 scale",limits=c(1,max(genedata$count)))
  #       } else {
  #         p <- p + scale_y_log10(name="Normalized counts - log10 scale")
  #       }
  #       exportPlots$genefinder <- p
  #
  #       p
  #     })
  #


  ## GENE FINDER
  output$searchresult <- renderPrint({

    if(is.null(input$color_by)) return("Select a factor to plot your gene")
    if(input$genefinder=="")
      return("Type in the gene name/id you want to plot")

    foundGeneID <- input$genefinder %in% rownames(values$myrlt)
    foundGeneName <- input$genefinder %in% values$myannotation$gene_name
    if(!foundGeneID){
      foundGeneID <- toupper(input$genefinder) %in% toupper(rownames(values$myrlt))
      if(foundGeneID){
        return(paste0("Maybe you mis-spelled the name of your gene. Did you mean ",
                      unique(rownames(values$myannotation)[which(toupper(input$genefinder)==toupper(rownames(values$myannotation)))]),"?"))
      } else {
        foundGeneNAME <- input$genefinder %in% values$myannotation$gene_name
        if(!foundGeneNAME){
          foundGeneNAME <- toupper(input$genefinder) %in% toupper(values$myannotation$gene_name)
          if(foundGeneNAME){
            return(paste0("Maybe you mis-spelled the name of your gene. Did you mean ",
                          unique(values$myannotation$gene_name[which(toupper(input$genefinder)==toupper(values$myannotation$gene_name))]),"?"))
          } else {return("Could not find the gene you typed!")}
        } else {
          fgn <- values$myannotation$gene_name[which(values$myannotation$gene_name==input$genefinder)]
          if (length(fgn) > 1) return(paste0("Found more than one gene with the selected gene name. Select one of the following: ",paste(selectedGene,collapse=", ")))
          selectedGene <- rownames(values$myannotation)[which(values$myannotation$gene_name==input$genefinder)]

          fg <- rownames(values$myannotation)[match(fgn,values$myannotation$gene_name)]
          return(paste0("I found the gene! Plotting ", fg, " - ", values$myannotation$gene_name[match(fg,rownames(values$myannotation))],"..."))

        }}
    } else {
      fg <- rownames(values$myannotation)[match(input$genefinder,rownames(values$myrlt))]
      return(paste0("I found the gene! Plotting ", fg, " - ", values$myannotation$gene_name[match(fg,rownames(values$myannotation))],"..."))

    }
  })



  output$genefinder_plot <- renderPlot({
    anno_id <- rownames(values$myrlt)
    anno_gene <- values$myannotation$gene_name

    if(is.null(input$color_by) & input$genefinder!="")
      return(ggplot() + annotate("text",label="Select a factor to plot your gene",0,0) + theme_bw())
    if(is.null(input$color_by) & input$genefinder=="")
      return(ggplot() + annotate("text",label="Select a gene and a factor to plot gene",0,0) + theme_bw())
    if(input$genefinder=="")
      return(ggplot() + annotate("text",label="Type in a gene name/id",0,0) + theme_bw())
    if(!input$genefinder %in% anno_gene & !input$genefinder %in% anno_id)
      return(ggplot() + annotate("text",label="Gene not found...",0,0) + theme_bw())

    if (input$genefinder %in% anno_id) {
      selectedGene <- rownames(values$myrlt)[match(input$genefinder,rownames(values$myrlt))]
      selectedGeneSymbol <- values$myannotation$gene_name[match(selectedGene,rownames(values$myannotation))]
    }
    if (input$genefinder %in% anno_gene) {
      selectedGeneSymbol <- values$myannotation$gene_name[which(values$myannotation$gene_name==input$genefinder)]
      if (length(selectedGeneSymbol) > 1) return(ggplot() + annotate("text",label=paste0("Type in a gene name/id of the following:\n",paste(selectedGene,collapse=", ")),0,0) + theme_bw())
      selectedGene <- rownames(values$myannotation)[which(values$myannotation$gene_name==input$genefinder)]
    }
    genedata <- plotCounts(values$mydds,gene=selectedGene,intgroup = input$color_by,returnData = TRUE)
    onlyfactors <- genedata[,match(input$color_by,colnames(genedata))]
    genedata$plotby <- interaction(onlyfactors)
    genedata$sampleID <- rownames(genedata)


    # input$plot_style chooses the style of plotting
    if(input$plot_style=="boxplot"){
      res <- ggplot(genedata,aes_string(x="plotby",y="count",fill="plotby")) +
        geom_boxplot(outlier.shape = NA,alpha=0.7) + theme_bw()
      if(input$ylimZero){
        res <- res + scale_y_log10(name="Normalized counts - log10 scale",limits=c(0.4,NA))
      } else {
        res <- res + scale_y_log10(name="Normalized counts - log10 scale")
      }

      res <- res +
        labs(title=paste0("Normalized counts for ",selectedGeneSymbol," - ",selectedGene)) +
        scale_x_discrete(name="") +
        geom_jitter(aes_string(x="plotby",y="count"),position = position_jitter(width = 0.1)) +
        scale_fill_discrete(name="Experimental\nconditions")
      if(input$addsamplelabels){
        res <- res + geom_text(aes_string(label="sampleID"),hjust=-.1,vjust=0)
      }
      exportPlots$genefinder_countsplot <- res
      res
    } else if(input$plot_style=="violin plot"){
      res <- ggplot(genedata,aes_string(x="plotby",y="count",fill="plotby")) +
        geom_violin(aes_string(col="plotby"),alpha = 0.6) + theme_bw()
      if(input$ylimZero){
        res <- res + scale_y_log10(name="Normalized counts - log10 scale",limits=c(0.4,NA))
      } else {
        res <- res + scale_y_log10(name="Normalized counts - log10 scale")
      }

      res <- res +
        labs(title=paste0("Normalized counts for ",selectedGeneSymbol," - ",selectedGene)) +
        scale_x_discrete(name="") +
        geom_jitter(aes_string(x="plotby",y="count"),alpha = 0.8,position = position_jitter(width = 0.1)) +
        scale_fill_discrete(name="Experimental\nconditions") + scale_color_discrete(guide="none")
      if(input$addsamplelabels){
        res <- res + geom_text(aes_string(label="sampleID"),hjust=-.1,vjust=0)
      }
      exportPlots$genefinder_countsplot <- res
      res
    }
  })


  output$genefinder_table <- DT::renderDataTable({
    anno_id <- rownames(values$myrlt)
    anno_gene <- values$myannotation$gene_name

    if(is.null(input$color_by) & input$genefinder!="")
      return(NULL)
    if(is.null(input$color_by) & input$genefinder=="")
      return(NULL)
    if(input$genefinder=="")
      return(NULL)
    if(!input$genefinder %in% anno_gene & !input$genefinder %in% anno_id)
      return(NULL)

    if (input$genefinder %in% anno_id) {
      selectedGene <- rownames(values$myrlt)[match(input$genefinder,rownames(values$myrlt))]
      selectedGeneSymbol <- values$myannotation$gene_name[match(selectedGene,rownames(values$myannotation))]
    }
    if (input$genefinder %in% anno_gene) {
      selectedGeneSymbol <- values$myannotation$gene_name[which(values$myannotation$gene_name==input$genefinder)]
      if (length(selectedGeneSymbol) > 1) return(ggplot() + annotate("text",label=paste0("Type in a gene name/id of the following:\n",paste(selectedGene,collapse=", ")),0,0) + theme_bw())
      selectedGene <- rownames(values$myannotation)[which(values$myannotation$gene_name==input$genefinder)]
    }
    genedata <- plotCounts(values$mydds,gene=selectedGene,intgroup = input$color_by,returnData = TRUE)
    genedata
  })


  ## PCA2GO
  output$ui_computePCA2GO <- renderUI({
    if(is.null(pca2go))
      actionButton("computepca2go","Compute the PCA2GO object",icon=icon("spinner"),class = "btn btn-primary")
  })

  annoSpecies_df <- data.frame(species=c("","Anopheles","Arabidopsis","Bovine","Worm",
                                         "Canine","Fly","Zebrafish","E coli strain K12",
                                         "E coli strain Sakai","Chicken","Human","Mouse",
                                         "Rhesus","Malaria","Chimp","Rat",
                                         "Yeast","Streptomyces coelicolor", "Pig","Toxoplasma gondii",
                                         "Xenopus"),
                               pkg=c("","org.Ag.eg.db",	"org.At.tair.db", "org.Bt.eg.db",	"org.Ce.eg.db",
                                     "org.Cf.eg.db",	"org.Dm.eg.db", "org.Dr.eg.db",	"org.EcK12.eg.db",
                                     "org.EcSakai.eg.db","org.Gg.eg.db","org.Hs.eg.db",	"org.Mm.eg.db",
                                     "org.Mmu.eg.db","org.Pf.plasmo.db","org.Pt.eg.db","org.Rn.eg.db",
                                     "org.Sc.sgd.db","org.Sco.eg.db",	"org.Ss.eg.db","org.Tgondii.eg.db",
                                     "org.Xl.eg.db"),
                               stringsAsFactors = FALSE)
  annoSpecies_df <- annoSpecies_df[order(annoSpecies_df$species),]
  annoSpecies_df <- annoSpecies_df[annoSpecies_df$species %in% c("","Human", "Mouse", "Rat", "Fly", "Chimp"),]

  output$ui_selectspecies <- renderUI({
    if(is.null(values$mypca2go)) {
      selectInput("speciesSelect",label = "Select the species of your samples",choices = annoSpecies_df$species,selected="")
    }
  })
  output$ui_inputtype <- renderUI({
    if(is.null(values$mypca2go)) {
      selectInput("idtype",label = "Select the input type of your identifiers",
                  choices = c("ENSEMBL","SYMBOL","REFSEQ","ENTREZID"), selected = "ENSEMBL")
    }
  })


  output$speciespkg <- renderText({

    if(!is.null(pca2go))
      return("pca2go object provided")

    if(!is.null(values$mypca2go))
      return("pca2go object computed or provided")

    shiny::validate(
      need(input$speciesSelect!="",
           "Select a species - requires the corresponding annotation package"
      )
    )

    annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==input$speciesSelect]

    shiny::validate(
      need(require(annopkg,character.only=TRUE),
           paste0("The package ",annopkg, " is not installed/available. Try installing it with biocLite('",annopkg,"')"))
    )

    retmsg <- paste0(annopkg," - package available and loaded")
    # if (!require(annopkg,character.only=TRUE)) {
    # stop("The package",annopkg, "is not installed/available. Try installing it with biocLite() ?")
    # }
    retmsg <- paste0(retmsg," - ",gsub(".eg.db","",gsub("org.","",annopkg)))
    retmsg

  })




  computedPCA2GO <- eventReactive( input$computepca2go, {
    annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==input$speciesSelect]
    withProgress(message = "Computing the PCA2GO object...",
                 value = 0,
                 {
                   pcpc <- limmaquickpca2go(values$myrlt,background_genes = rownames(values$mydds),
                                            inputType = input$idtype,
                                            organism = gsub(".eg.db","",gsub("org.","",annopkg)))
                 })
    pcpc
  })


  observeEvent(input$computepca2go,
               {
                 values$mypca2go <- computedPCA2GO()
               })



  output$pca2go <- renderPlot({
    shiny::validate(
      need(
        !is.null(values$mypca2go),
        "Please provide a pca2go object to the app or alternatively click on the action button - could take some time to compute live!"
      )
    )
    # if(is.null(pca2go))
    # return(ggplot() + annotate("text",label="Provide a pca2go object to the app",0,0) + theme_bw())
    res <- pcaplot(values$myrlt,intgroup = input$color_by,
                   ntop = attr(values$mypca2go,"n_genesforpca"),
                   pcX = as.integer(input$pc_x),pcY = as.integer(input$pc_y),text_labels = input$sample_labels,
                   point_size = input$pca_point_size, title=paste0("PCA on the samples - ",attr(values$mypca2go,"n_genesforpca"), " genes used")

    )
    res
  })


  output$dt_pchor_pos <- DT::renderDataTable({
    if(is.null(values$mypca2go)) return(datatable(NULL))
    goe <- values$mypca2go[[paste0("PC",input$pc_x)]][["posLoad"]]
    if(input$compact_pca2go)
      return(datatable(goe[,c("GO.ID","Term","Significant","p.value_elim")],options = list(pageLength = 5)))
    datatable(goe)
  })

  output$dt_pchor_neg <- DT::renderDataTable({
    if(is.null(values$mypca2go)) return(datatable(NULL))
    goe <- values$mypca2go[[paste0("PC",input$pc_x)]][["negLoad"]]
    if(input$compact_pca2go)
      return(datatable(goe[,c("GO.ID","Term","Significant","p.value_elim")],options = list(pageLength = 5)))
    datatable(goe)
  })

  output$dt_pcver_pos <- DT::renderDataTable({
    if(is.null(values$mypca2go)) return(datatable(NULL))
    goe <- values$mypca2go[[paste0("PC",input$pc_y)]][["posLoad"]]
    if(input$compact_pca2go)
      return(datatable(goe[,c("GO.ID","Term","Significant","p.value_elim")],options = list(pageLength = 5)))
    datatable(goe)
  })

  output$dt_pcver_neg <- DT::renderDataTable({
    if(is.null(values$mypca2go)) return(datatable(NULL))
    goe <- values$mypca2go[[paste0("PC",input$pc_y)]][["negLoad"]]
    if(input$compact_pca2go)
      return(datatable(goe[,c("GO.ID","Term","Significant","p.value_elim")],options = list(pageLength = 5)))
    datatable(goe)
  })

  output$enrichinfo <- renderPrint({
    cat("enrich info:\n")
    # str(goEnrichs)
    class(input$pc_x)
    head(values$mypca2go[[paste0("PC",input$pc_x)]][["posLoad"]])
    class(datatable(values$mypca2go[[paste0("PC",input$pc_x)]][["posLoad"]]))
  })






  ## MULTIFACTOR EXPLORATION
  output$intro_multifac <- renderText({
    if(!is.null(values$mydds))
      shiny::validate(
        need(ncol(colData(values$mydds)) > 1,
             message = "To use this section, you need a dataset where more than one experimental factor is available.")
      )
    return("Refer to the Instructions section if you need help on using this section")
  })


  output$covar1 <- renderUI({
    # if(is.null(values$myrlt))
    # return(NULL)
    poss_covars <- names(colData(values$mydds))
    selectInput('covar1', label = 'Select factor 1: ',
                choices = c(NULL, poss_covars), selected = NULL,multiple = FALSE)
  })

  output$covar2 <- renderUI({
    # if(is.null(values$myrlt))
    # return(NULL)
    poss_covars <- names(colData(values$mydds))
    selectInput('covar2', label = 'Select factor 2: ',
                choices = c(NULL, poss_covars), selected = NULL,multiple = FALSE)
  })


  output$c1levels <- renderUI({
    if(is.null(input$covar1))
      return(NULL)
    fac1lev <- levels(colData(values$myrlt)[[input$covar1]])
    selectInput('covar1levels', label = 'Factor 1 available levels: ',
                choices = c(NULL, fac1lev), selected = NULL,multiple = TRUE) # actually 2
  })

  output$c2levels <- renderUI({
    if(is.null(input$covar2))
      return(NULL)
    fac2lev <- levels(colData(values$myrlt)[[input$covar2]])
    selectInput('covar2levels', label = 'Factor 2 available levels: ',
                choices = c(NULL, fac2lev), selected = NULL,multiple = TRUE) # 2 or more are allowed!
  })

  output$colnames1 <- renderUI({
    if(is.null(values$myrlt))
      return(NULL)
    if(is.null(input$covar1))
      return(NULL)
    if(is.null(input$covar2))
      return(NULL)

    fac1 <- input$covar1
    fac2 <- input$covar2

    fac1_touse <- input$covar1levels
    fac2_touse <- input$covar2levels

    preselected_fac1 <- colnames(values$myrlt)[colData(values$myrlt)[[fac1]] %in% fac1_touse]
    preselected_fac2 <- colnames(values$myrlt)[colData(values$myrlt)[[fac2]] %in% fac2_touse]
    presel <- intersect(preselected_fac1,preselected_fac2)
    mysamples <- colData(values$myrlt)[presel,] # check that the repl are balanced

    presel1 <- colnames(values$myrlt)[(colData(values$myrlt)[[fac1]] %in% fac1_touse[1]) & colData(values$myrlt)[[fac2]] %in% fac2_touse]

    selectInput('picksamples1', label = 'Combine samples from Factor1-Level1 in the selected order: ',
                choices = c(NULL, presel1), selected = NULL,multiple = TRUE)
  })


  output$colnames2 <- renderUI({
    if(is.null(values$myrlt))
      return(NULL)
    if(is.null(input$covar1))
      return(NULL)
    if(is.null(input$covar2))
      return(NULL)

    fac1 <- input$covar1
    fac2 <- input$covar2

    fac1_touse <- input$covar1levels
    fac2_touse <- input$covar2levels

    preselected_fac1 <- colnames(values$myrlt)[colData(values$myrlt)[[fac1]] %in% fac1_touse]
    preselected_fac2 <- colnames(values$myrlt)[colData(values$myrlt)[[fac2]] %in% fac2_touse]
    presel <- intersect(preselected_fac1,preselected_fac2)
    mysamples <- colData(values$myrlt)[presel,] # check that the repl are balanced

    presel2 <- colnames(values$myrlt)[(colData(values$myrlt)[[fac1]] %in% fac1_touse[2]) & colData(values$myrlt)[[fac2]] %in% fac2_touse]
    selectInput('picksamples2', label = 'Combine samples from Factor1-Level2 in the selected order: ',
                choices = c(NULL, presel2), selected = NULL,multiple = TRUE)


  })



  composedMat <- eventReactive( input$composemat, {
    exprmat <- t(assay(values$myrlt))
    exprmat <- exprmat[,rowSums(counts(values$mydds) > 5)>2]

    withProgress(message = "Composing the matrix...",
                 value = 0,
                 {
                   pcmat <- cbind(exprmat[input$picksamples1,],
                                  exprmat[input$picksamples2,])
                 })
    pcmat
  })


  obj3 <- reactive({

    pcmat <- composedMat()
    aval <- 0.3
    fac2pal <- alpha(c("green","red","blue","orange","violet"),aval) # 5 are enough

    # colData(values$myrlt)[input$covar2][rownames(pcmat),]
    max.type <- apply(pcmat[,1:(ncol(pcmat)/2)],2,which.max)

    fac2_col <- factor(colData(values$myrlt)[input$covar2][rownames(pcmat),],
                       levels=unique(as.character(colData(values$myrlt)[input$covar2][rownames(pcmat),])))
    tcol.justMax <- fac2pal[fac2_col][max.type]
    # tcol.justMax <- ifelse(max.type <= 4,"green",ifelse(max.type <= 8,"red",ifelse(max.type <= 12,"blue","orange")))

    max.type2 <- apply(pcmat[,((ncol(pcmat)/2)+1):ncol(pcmat)],2,which.max)
    # tcol2.justMax <- ifelse(max.type2 <= 4,alpha("green",aval),ifelse(max.type2 <= 8,alpha("red",aval),ifelse(max.type2 <= 12,alpha("blue",aval),alpha("orange",aval))))

    tcol2.justMax <- fac2pal[fac2_col][max.type2]

    # using the median across replicates
    celltypes <- gsub("_R.","",rownames(pcmat))

    tcol <- tcol.justMax
    tcol2 <- tcol2.justMax
    # pcmat
    return(list(pcmat,tcol,tcol2))
  })


  output$pcamultifac <- renderPlot({
    pcmat <- obj3()[[1]]
    tcol <- obj3()[[2]]
    tcol2 <- obj3()[[3]]
    pres <- prcomp(t(pcmat),scale=FALSE)

    plot.index <- c(as.integer(input$pc_x_multifac),as.integer(input$pc_y_multifac))
    offset <- ncol(pcmat)/2
    gene.no <- offset
    pcx <- pres$x
    # set.seed(11)
    # for (i in 1:ncol(pcx)) {
    #   pcx[,i] <- pcx[,i] + rnorm(nrow(pcx),sd=diff(range(pcx[,i]))/100)
    # }
    plot(pcx[(offset+1):ncol(pcmat),plot.index[1]][1:gene.no],pcx[(offset+1):ncol(pcmat),plot.index[2]][1:gene.no],xlim=range(pcx[,plot.index[1]]),ylim=range(pcx[,plot.index[2]]),pch=20,col=tcol,cex=0.3)#,type="n")
    #plot(0,type="n",xlim=range(pres$x[,plot.index]),ylim=range(pres$x[,plot.index]))
    lcol <- ifelse(tcol != tcol2,"black","grey")
    for (i in 1:gene.no) {
      lines(pcx[c(i,offset+i),plot.index[1]],pcx[c(i,offset+i),plot.index[2]],col=lcol[i])
    }
    points(pcx[1:offset,plot.index[1]][1:gene.no],pcx[1:offset,plot.index[2]][1:gene.no],pch=20,col=tcol,cex=0.3)
    points(pcx[(offset+1):ncol(pcmat),plot.index[1]][1:gene.no],pcx[(offset+1):ncol(pcmat),plot.index[2]][1:gene.no],pch=20,col=tcol2,cex=0.3)

  })


  output$multifaczoom <- renderPlot({
    if(is.null(input$pcamultifac_brush)) return(NULL)
    pcmat <- obj3()[[1]]
    tcol <- obj3()[[2]]
    tcol2 <- obj3()[[3]]
    pres <- prcomp(t(pcmat),scale=FALSE)

    plot.index <- c(as.integer(input$pc_x_multifac),as.integer(input$pc_y_multifac))
    offset <- ncol(pcmat)/2
    gene.no <- offset
    pcx <- pres$x

    plot(pcx[(offset+1):ncol(pcmat),plot.index[1]][1:gene.no],
         pcx[(offset+1):ncol(pcmat),plot.index[2]][1:gene.no],
         xlim=c(input$pcamultifac_brush$xmin,input$pcamultifac_brush$xmax),
         ylim=c(input$pcamultifac_brush$ymin,input$pcamultifac_brush$ymax),
         pch=20,col=tcol,cex=0.3)#,type="n")
    #plot(0,type="n",xlim=range(pres$x[,plot.index]),ylim=range(pres$x[,plot.index]))
    lcol <- ifelse(tcol != tcol2,"black","grey")
    for (i in 1:gene.no) {
      lines(pcx[c(i,offset+i),plot.index[1]],pcx[c(i,offset+i),plot.index[2]],col=lcol[i])
    }
    points(pcx[1:offset,plot.index[1]][1:gene.no],pcx[1:offset,plot.index[2]][1:gene.no],pch=20,col=tcol,cex=0.3)
    points(pcx[(offset+1):ncol(pcmat),plot.index[1]][1:gene.no],pcx[(offset+1):ncol(pcmat),plot.index[2]][1:gene.no],pch=20,col=tcol2,cex=0.3)

  })


  #       output$plot_brushinfo <- renderPrint({
  #         cat("input$pcagenes_brush:\n")
  #         str(input$pcagenes_brush)
  #       })

  curData_brush_multifac <- reactive({
    pcmat <- obj3()[[1]]
    tcol <- obj3()[[2]]
    tcol2 <- obj3()[[3]]

    pres <- prcomp(t(pcmat),scale=FALSE)

    plot.index <- c(as.integer(input$pc_x_multifac),as.integer(input$pc_y_multifac))
    offset <- ncol(pcmat)/2
    gene.no <- offset
    pcx <- pres$x


    firstPCselected <- c(
      pcx[1:offset,plot.index[1]][1:gene.no],
      pcx[(offset+1):ncol(pcmat),plot.index[1]][1:gene.no])

    secondPCselected <- c(
      pcx[1:offset,plot.index[2]][1:gene.no],
      pcx[(offset+1):ncol(pcmat),plot.index[2]][1:gene.no]
    )

    pcspcs <- data.frame(firstPC=firstPCselected,secondPC=secondPCselected,geneID=colnames(pcmat))
    rownames(pcspcs) <- c(paste0(colnames(pcmat)[1:gene.no],"_WT"),
                          paste0(colnames(pcmat)[(gene.no+1):(2*gene.no)],"_G37"))

    if(!is.null(values$myannotation))
      pcspcs$geneName <- values$myannotation$gene_name[match(pcspcs$geneID,rownames(values$myannotation))]


    res <- brushedPoints(pcspcs, input$pcamultifac_brush,xvar="firstPC",yvar="secondPC",)
    res
  })



  output$pcamultifac_out <- DT::renderDataTable({
    datatable(curData_brush_multifac())
  })


  output$downloadData_brush_multifac <- downloadHandler(
    filename = function() { paste(input$brushedPoints_filename_multifac, '.csv', sep='') },
    content = function(file) {
      if(length(input$pcamultifac_out_rows_selected)){
        data <- curData_brush_multifac()[input$pcamultifac_out_rows_selected,]
      } else {
        data <- curData_brush_multifac()
      }
      write.csv(data, file, quote=FALSE)
    }
  )







  ## REPORT EDITOR
  ### yaml generation
  rmd_yaml <- reactive({
    paste0("---",
           "\ntitle: '", input$report_title,
           "'\nauthor: '", input$report_author,
           "'\ndate: '", Sys.Date(),
           "'\noutput:\n  html_document:\n    toc: ", input$report_toc, "\n    number_sections: ", input$report_ns, "\n    theme: ", input$report_theme, "\n---\n\n",collapse = "\n")
  })


  # rmd_full <- reactive({
  #   paste0(rmd_yaml(),"\n",
  #          readLines("reportTemplate.Rmd"))
  # })
  # output$loadedRmd <- renderPrint({
  #   # rmd_yaml() # or rmd_full()
  #   paste0(
  #     # rmd_yaml(),
  #     paste0(readLines("reportTemplate.Rmd"),collapse = "\n"))
  #   # head(paste0(rmd_yaml(),
  #   # readLines("reportTemplate.Rmd")),collapse="\n")
  # })

  ### loading report template
  # update aceEditor module
  observe({
    # loading rmd report from disk
    inFile <- system.file("extdata", "reportTemplate.Rmd",package = "pcaExplorer")

    isolate({
      if(!is.null(inFile) && !is.na(inFile)) {

        rmdfilecontent <- paste0(readLines(inFile),collapse="\n")

        shinyAce::updateAceEditor(session, "acereport_rmd", value = rmdfilecontent)
      }
    })
  })


  ### ace editor options
  observe({
    autoComplete <- if(input$enableAutocomplete) {
      if(input$enableLiveCompletion) "live" else "enabled"
    } else {
      "disabled"
    }

    updateAceEditor(session, "acereport_rmd", autoComplete = autoComplete,theme=input$theme, mode=input$mode)
    # updateAceEditor(session, "plot", autoComplete = autoComplete)
  })

  #Enable/Disable R code completion
  rmdOb <- aceAutocomplete("acereport_rmd")
  observe({
    if(input$enableRCompletion) {
      rmdOb$resume()
    } else {
      rmdOb$suspend()
    }
  })

  ## currently not working as I want with rmarkdown::render, but can leave it like this - the yaml will be taken in the final version only
  output$knitDoc <- renderUI({
    input$updatepreview_button
    return(isolate(HTML(knit2html(text = input$acereport_rmd, fragment.only = TRUE, quiet = TRUE))))
  })





  ## STATE SAVING
  ### to environment
  observe({
    if(is.null(input$exit_and_save) || input$exit_and_save ==0 ) return()

    # quit R, unless you are running an interactive session
    if(interactive()) {
      # flush input and values to the environment in two distinct objects (to be reused later?)
      isolate({
        assign(paste0("pcaExplorer_inputs_",
                      gsub(" ","_",gsub("-","",gsub(":","-",as.character(Sys.time()))))),
               reactiveValuesToList(input), envir = .GlobalEnv)
        assign(paste0("pcaExplorer_values_",
                      gsub(" ","_",gsub("-","",gsub(":","-",as.character(Sys.time()))))),
               reactiveValuesToList(values), envir = .GlobalEnv)
        stopApp("pcaExplorer closed, state successfully saved to global R environment.")
      })
    } else {
      stopApp("pcaExplorer closed")
      q("no")
    }
  })

  ### to binary data
  saveState <- function(filename) {
    isolate({
      LiveInputs <- reactiveValuesToList(input)
      # values[names(LiveInputs)] <- LiveInputs
      r_data <- reactiveValuesToList(values)
      save(LiveInputs, r_data , file = filename)
    })
  }

  output$state_save_sc <- downloadHandler(
    filename = function() {
      paste0("pcaExplorerState_",gsub(" ","_",gsub("-","",gsub(":","-",as.character(Sys.time())))),".RData")
    },
    content = function(file) {
      saveState(file)
    }
  )








  ## all download handlers
  output$downloadData_brush <- downloadHandler(
    filename = function() { paste(input$brushedPoints_filename, '.csv', sep='') },
    content = function(file) {
      if(length(input$pca_brush_out_rows_selected)){
        data <- curData_brush()[input$pca_brush_out_rows_selected,]
      } else {
        data <- curData_brush()
      }
      write.csv(data, file, quote=FALSE)
    }
  )

  output$downloadData_click <- downloadHandler(
    filename = function() { paste(input$clickedPoints_filename, '.csv', sep='') },
    content = function(file) {
      write.csv(curData_click(), file, quote=FALSE)
    }
  )



  output$download_samplessamplesheat <- downloadHandler(filename=function(){
    input$filename_samplessamplesheat
  },
  content = function(file){
    pdf(file)

    if (!is.null(input$color_by)){
      expgroups <- as.data.frame(colData(values$myrlt)[,input$color_by])
      # expgroups <- interaction(expgroups)
      rownames(expgroups) <- colnames(values$myrlt)
      colnames(expgroups) <- input$color_by

      pheatmap(as.matrix(dist(t(assay(values$myrlt)))),annotation_col = expgroups)
    } else {
      pheatmap(as.matrix(dist(t(assay(values$myrlt)))))
    }

    dev.off()
  })



  output$download_readsbarplot <- downloadHandler(
    filename = function() { input$filename_readsbarplot },
    content = function(file) {
      ggsave(file, exportPlots$reads_barplot, width = input$export_width, height = input$export_height, units = "cm")
    })

  output$download_samplesPca <- downloadHandler(
    filename = function() { input$filename_samplesPca },
    content = function(file) {
      ggsave(file, exportPlots$samplesPca, width = input$export_width, height = input$export_height, units = "cm")
    })

  output$download_samplesScree <- downloadHandler(
    filename = function() { input$filename_samplesScree },
    content = function(file) {
      ggsave(file, exportPlots$samplesScree, width = input$export_width, height = input$export_height, units = "cm")
    })

  output$download_samplesPcazoom <- downloadHandler(
    filename = function() { input$filename_samplesPcazoom },
    content = function(file) {
      ggsave(file, exportPlots$samplesZoom, width = input$export_width, height = input$export_height, units = "cm")
    })

  output$download_samplesPca_hiload <- downloadHandler(filename=function(){
    input$filename_samplesPca_hiload
  },
  content = function(file){
    pdf(file)

    rv <- rowVars(assay(values$myrlt))
    select <- order(rv, decreasing = TRUE)[seq_len(min(input$pca_nrgenes,length(rv)))]
    pca <- prcomp(t(assay(values$myrlt)[select, ]))

    par(mfrow=c(2,1))
    hi_loadings(pca,whichpc = as.integer(input$pc_x),topN = input$ntophiload,annotation = values$myannotation)
    hi_loadings(pca,whichpc = as.integer(input$pc_y),topN = input$ntophiload,annotation = values$myannotation)

    dev.off()
  })

  output$download_samplesPca_sampleout <- downloadHandler(
    filename = function() { input$filename_samplesPca_sampleout },
    content = function(file) {
      ggsave(file, exportPlots$samplesOutlier, width = input$export_width, height = input$export_height, units = "cm")
    })




  output$download_genesPca <- downloadHandler(
    filename = function() { input$filename_genesPca },
    content = function(file) {
      ggsave(file, exportPlots$genesPca, width = input$export_width, height = input$export_height, units = "cm")
    })

  output$download_genesZoom <- downloadHandler(
    filename = function() { input$filename_genesZoom },
    content = function(file) {
      ggsave(file, exportPlots$genesZoom, width = input$export_width, height = input$export_height, units = "cm")
    })


  output$download_genesPca_profile <- downloadHandler(
    filename=function(){
      input$filename_genesPca_profile
    },
    content = function(file){
      pdf(file)
      geneprofiler(values$myrlt,
                   genelist = curData_brush()$ids,
                   intgroup = input$color_by,
                   plotZ = input$zprofile)
      dev.off()
    })

  output$download_genesPca_countsplot <- downloadHandler(
    filename = function() { input$filename_genesPca_countsplot },
    content = function(file) {
      ggsave(file, exportPlots$genesBoxplot, width = input$export_width, height = input$export_height, units = "cm")
    })




  output$download_genesHeatmap <- downloadHandler(
    filename=function(){
      input$filename_genesHeatmap
    },
    content = function(file){
      pdf(file)
      brushedObject <- curData_brush()

      selectedGenes <- brushedObject$ids
      toplot <- assay(values$myrlt)[selectedGenes,]
      rownames(toplot) <- values$myannotation$gene_name[match(rownames(toplot),rownames(values$myannotation))]
      aheatmap(toplot,Colv = as.logical(input$heatmap_colv))
      dev.off()
    })


  output$download_genefinder_countsplot <- downloadHandler(
    filename = function() { input$filename_genefinder_countsplot },
    content = function(file) {
      ggsave(file, exportPlots$genefinder_countsplot, width = input$export_width, height = input$export_height, units = "cm")
    })




  # Generate and Download module
  output$saveRmd <- downloadHandler(
    filename = function() {
      if(input$rmd_dl_format == "rmd") {
        "report.Rmd"
      } else {
        "report.html"
      }
    },
    content = function(file) {

      # knit2html(text = input$rmd, fragment.only = TRUE, quiet = TRUE))

      tmp_content <-
        paste0(rmd_yaml(),
               input$acereport_rmd,collapse = "\n")
      # input$acereport_rmd
      if(input$rmd_dl_format == "rmd") {
        cat(tmp_content,file=file,sep="\n")
      } else {
        # write it somewhere too keeping the source
        # tmpfile <- tempfile()
        # file.create(tmpfile)
        # fileConn<- file(tempfile())
        # writeLines(tmp_content, fileConn)
        # close(fileConn)
        if(input$rmd_dl_format == "html") {
          cat(tmp_content,file="tempreport.Rmd",sep="\n")
          rmarkdown::render(input = "tempreport.Rmd",
                            output_file = file,
                            # fragment.only = TRUE,
                            quiet = TRUE)
        }
      }
    })




}) # end of pcaExplorer(dds,rlt,countmatrix,coldata,pca2go,annotation)
