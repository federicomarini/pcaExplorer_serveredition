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


library(shinydashboard)
library(shiny)
library(d3heatmap)
library(pcaExplorer)
library(DESeq2)
library(ggplot2)
library(shinyAce)
library(DT)
library(knitr)
library(rmarkdown)
library(pheatmap)
library(threejs)
library(scales)
library(genefilter)

# DESeq2, SummarizedExperiment, GenomicRanges, IRanges,
# S4Vectors, genefilter, ggplot2 (>= 2.0.0), d3heatmap, scales,
# NMF, plyr, topGO, limma, GOstats, GO.db, AnnotationDbi, shiny
# (>= 0.12.0), shinydashboard, shinyBS, ggrepel, DT, shinyAce,
# threejs, biomaRt, pheatmap, knitr, rmarkdown, tidyr, grDevices,
# methods
# Suggests: testthat, BiocStyle, airway, org.Hs.eg.db
#






shinydashboard::dashboardPage(

  dashboardHeader(
    title = paste0("pcaExplorer - Interactive exploration of Principal Components ",
                   "of Samples and Genes in RNA-seq data - version ",
                   packageVersion("pcaExplorer")),
    titleWidth = 900,

    # task menu for saving state to environment or binary data
    shinydashboard::dropdownMenu(type = "tasks",icon = icon("cog"),badgeStatus = "success",
                                 notificationItem(
                                   text = actionButton("exit_and_save","Exit pcaExplorer & save",
                                                       class = "btn_no_border",
                                                       onclick = "setTimeout(function(){window.close();}, 100); "),
                                   icon = icon("sign-out"),status = "primary"),
                                 menuItem(
                                   text = downloadButton("state_save_sc","Save State as .RData"))
    )
  ),

  dashboardSidebar(
    width = 280,

    menuItem("App settings",icon = icon("cogs"),
             selectInput('pc_x', label = 'x-axis PC: ', choices = 1:8, selected = 1),
             shinyBS::bsTooltip(
               "pc_x", paste0("Select the principal component to display on the x axis"),
               "right", options = list(container = "body")),
             selectInput('pc_y', label = 'y-axis PC: ', choices = 1:8, selected = 2),
             shinyBS::bsTooltip(
               "pc_y", paste0("Select the principal component to display on the y axis"),
               "right", options = list(container = "body")),
             uiOutput("color_by"),
             shinyBS::bsTooltip(
               "color_by", paste0("Select the group of samples to stratify the analysis. Can also assume multiple values"),
               "right", options = list(container = "body")),
             numericInput('pca_nrgenes', label = 'Nr of (most variant) genes:', value = 300,min = 50,max = 20000),
             shinyBS::bsTooltip(
               "pca_nrgenes", paste0("Number of genes to select for computing the principal components. The top n genes are",
                                     " selected ranked by their variance inter-samples"),
               "right", options = list(container = "body")),
             numericInput('pca_point_alpha', label = 'Alpha: ', value = 1,min = 0,max = 1,step = 0.01),
             shinyBS::bsTooltip(
               "pca_point_alpha", paste0("Color transparency for the plots. Can assume values from 0 (transparent) ",
                                         "to 1 (opaque)"),
               "right", options = list(container = "body")),
             numericInput('pca_label_size', label = 'Labels size: ', value = 2,min = 1,max = 8),
             shinyBS::bsTooltip(
               "pca_label_size", paste0("Size of the labels for the samples in the principal components plots"),
               "right", options = list(container = "body")),
             numericInput('pca_point_size', label = 'Points size: ', value = 2,min = 1,max = 8),
             shinyBS::bsTooltip(
               "pca_point_size", paste0("Size of the points to be plotted in the principal components plots"),
               "right", options = list(container = "body")),
             numericInput('pca_varname_size', label = 'Variable name size: ', value = 4,min = 1,max = 8),
             shinyBS::bsTooltip(
               "pca_varname_size", paste0("Size of the labels for the genes PCA - correspond to the samples names"),
               "right", options = list(container = "body")),
             numericInput('pca_scale_arrow', label = 'Scaling factor : ', value = 1,min = 0.01,max = 10),
             shinyBS::bsTooltip(
               "pca_scale_arrow", paste0("Scale value for resizing the arrow corresponding to the variables in the ",
                                         "PCA for the genes. It should be used for mere visualization purposes"),
               "right", options = list(container = "body")),
             selectInput("col_palette","Color palette",choices = list("hue","set1","rainbow")),
             selectInput("plot_style","Plot style for gene counts",choices = list("boxplot","violin plot")),
             shinyBS::bsTooltip(
               "col_palette", paste0("Select the color palette to be used in the principal components plots. The number of ",
                                     "colors is selected automatically according to the number of samples and to the levels ",
                                     "of the factors of interest and their interactions"),
               "right", options = list(container = "body"))
    ),
    menuItem("Plot export settings", icon = icon("paint-brush"),

             numericInput("export_width",label = "Width of exported figures (cm)",value = 10,min = 2),
             shinyBS::bsTooltip(
               "export_width", paste0("Width of the figures to export, expressed in cm"),
               "right", options = list(container = "body")),
             numericInput("export_height",label = "Height of exported figures (cm)",value = 10,min = 2),
             shinyBS::bsTooltip(
               "export_height", paste0("Height of the figures to export, expressed in cm"),
               "right", options = list(container = "body"))

    )
  ),

  dashboardBody(

    ## Define output size and style of error messages
    tags$head(
      tags$style(HTML("
                      .shiny-output-error-validation {
                      font-size: 15px;
                      color: forestgreen;
                      text-align: center;
                      }
                      "))
      ),

    ## main structure of the body for the dashboard
    tabBox(
      width=12,

      tabPanel(
        "Data Upload", icon = icon("upload"),
        uiOutput("upload_count_matrix"),
        shinyBS::bsTooltip(
          "upload_count_matrix", paste0("Select file containing the count matrix"),
          "right", options = list(container = "body")),
        uiOutput("upload_metadata"),
        shinyBS::bsTooltip(
          "upload_metadata", paste0("Select file containing the samples metadata"),
          "right", options = list(container = "body")),
        uiOutput("upload_annotation"),
        shinyBS::bsTooltip(
          "upload_annotation", paste0("Select file containing the annotation data"),
          "right", options = list(container = "body")),
        # wellPanel(
        #
        #   checkboxInput("header_ct", "Header", TRUE),
        #   radioButtons("sep_ct", "Separator", c( Tab="\t",Comma=",", Semicolon=";", Space = " "), selected = "\t", inline = TRUE)
        #
        # ),
        # ## TODO: replace with heuristics for detecting separator with a guess
        p("Preview on the uploaded data"),

        verbatimTextOutput("printdds"),
        DT::dataTableOutput("printanno"),
        # this last one will just be visible if the user uploads the count matrix separately via button
        DT::dataTableOutput("sneakpeekcm")
      ),

      tabPanel(
        "Instructions",  icon = icon("info-circle"),
        includeMarkdown(system.file("extdata", "instructions.md",package = "pcaExplorer")),
        footer()
      ),

      tabPanel(
        "Counts Table",
        icon = icon("table"),
        h3("Counts table"),

        selectInput("countstable_unit", label = "Data scale in the table",
                    choices = list("Counts (raw)" = "raw_counts",
                                   "Counts (normalized)" = "normalized_counts",
                                   "Regularized logarithm transformed" = "rlog_counts",
                                   "Log10 (pseudocount of 1 added)" = "log10_counts",
                                   "TPM (Transcripts Per Million)" = "tpm_counts")),

        DT::dataTableOutput("showcountmat"),

        downloadButton("downloadData","Download", class = "btn btn-success"),
        hr(),
        h3("Sample to sample scatter plots"),
        selectInput("corr_method","Correlation method palette",choices = list("pearson","spearman")),
        p("Compute sample to sample correlations on the normalized counts - warning, it can take a while to plot all points (depending mostly on the number of samples you provided)."),
        actionButton("compute_pairwisecorr", "Run", class = "btn btn-primary"),
        uiOutput("pairwise_plotUI"),
        uiOutput("heatcorr_plotUI")

      ),

      tabPanel(
        "Data Overview", icon = icon("eye"),
        h1("Sneak peek in the data"),
        h3("Design metadata"),
        DT::dataTableOutput("showcoldata"),

        h3("Sample to sample distance heatmap"),
        fluidRow(
          column(
            width=8,
            plotOutput("heatmapsampledist"),
            div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                downloadButton("download_samplessamplesheat", "Download Plot"),
                textInput("filename_samplessamplesheat",label = "Save as...",value = "pcae_sampletosample.pdf")))
        ),
        hr(),
        h3("General information on the provided SummarizedExperiment/DESeqDataSet"),
        shiny::verbatimTextOutput("showdata"),
        h3("Number of million of reads per sample"),
        fluidRow(
          column(
            width=8,
            plotOutput("reads_barplot"),
            div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                downloadButton("download_readsbarplot", "Download Plot"),
                textInput("filename_readsbarplot",label = "Save as...",value = "pcae_readsbarplot.pdf")))),
        h3("Basic summary for the counts"),
        p("Number of uniquely aligned reads assigned to each sample"),
        verbatimTextOutput("reads_summary"),
        wellPanel(
          fluidRow(
            column(
              width = 4,
              numericInput("threshold_rowsums","Threshold on the row sums of the counts",value = 0, min = 0)),
            column(
              width = 4,
              numericInput("threshold_rowmeans","Threshold on the row means of the normalized counts",value = 0, min = 0))
          )),
        p("According to the selected filtering criteria, this is an overview on the provided count data"),
        verbatimTextOutput("detected_genes")
        # DT::dataTableOutput("reads_samples"),


      ),

      tabPanel(
        "Samples View",
        icon = icon("share-alt"),
        p(h1('Principal Component Analysis on the samples'),
          "PCA projections of sample expression profiles onto any pair of components."),
        fluidRow(
          column(
            width = 4,
            wellPanel(checkboxInput("sample_labels","Display sample labels",value = TRUE),
                      checkboxInput("pca_ellipse","draw a confidence ellipse for each group",value = FALSE),
                      sliderInput("pca_cislider", "select the confidence interval level", min=0,max=1,value=0.95)))),
        fluidRow(
          column(
            width = 6,
            plotOutput('samples_pca',brush = "pca_brush"),
            div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                downloadButton("download_samplesPca", "Download Plot"),
                textInput("filename_samplesPca",label = "Save as...",value = "samplesPca.pdf"))
          ),
          column(
            width= 6,
            plotOutput("samples_scree"),
            div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                downloadButton("download_samplesScree", "Download Plot"),
                textInput("filename_samplesScree",label = "Save as...",value = "samplesScree.pdf")),
            wellPanel(fluidRow(
              column(
                width = 6,
                radioButtons("scree_type","Scree plot type:",
                             choices=list("Proportion of explained variance"="pev",
                                          "Cumulative proportion of explained variance"="cev"),"pev")
              ),
              column(
                width = 6,
                numericInput("scree_pcnr","Number of PCs to display",value=8,min=2)
              )
            ))
          )
        ),
        hr(),
        fluidRow(
          column(
            width = 6,
            plotOutput("samples_pca_zoom"),
            div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                downloadButton("download_samplesPcazoom", "Download Plot"),
                textInput("filename_samplesPcazoom",label = "Save as...",value = "samplesPcazoom.pdf"))
          ),
          column(
            width = 6,
            numericInput("ntophiload", "Nr of genes to display (top & bottom)",value = 10, min = 1, max=40),
            plotOutput("geneshiload"),
            div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                downloadButton("download_samplesPca_hiload", "Download Plot"),
                textInput("filename_samplesPca_hiload",label = "Save as...",value = "pcae_hiload.pdf"))
          )
        ),
        hr(),
        fluidRow(
          column(
            width = 6,
            p(h4('Outlier Identification'), "Toggle which samples to remove - suspected to be considered as outliers"),
            uiOutput("ui_outliersamples"),
            plotOutput("samples_outliersremoved"),
            div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                downloadButton("download_samplesPca_sampleout", "Download Plot"),
                textInput("filename_samplesPca_sampleout",label = "Save as...",value = "samplesPca_sampleout.pdf"))
          )

        ),

        fluidRow(
          column(
            width = 8,
            selectInput("pc_z","Select the principal component to display on the z axis",choices = 1:8,selected = 3),
            scatterplotThreeOutput("pca3d")
          )
        )

      ),

      tabPanel(
        "Genes View",
        icon = icon("yelp"),
        p(h1('Principal Component Analysis on the genes'), "PCA projections of genes abundances onto any pair of components."),

        fluidRow(checkboxInput("variable_labels","Display variable labels",value = TRUE)),
        fluidRow(
          checkboxInput("ylimZero_genes","Set y axis limit to 0",value=TRUE)),

        fluidRow(
          column(
            width = 6,
            h4("Main Plot - interact!"),
            plotOutput('genes_biplot',brush = 'pcagenes_brush',click="pcagenes_click"),
            div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                downloadButton("download_genesPca", "Download Plot"),
                textInput("filename_genesPca",label = "Save as...",value = "genesPca.pdf"))),
          column(
            width = 6,
            h4("Zoomed window"),
            plotOutput("genes_biplot_zoom",click="pcagenes_zoom_click"),
            div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                downloadButton("download_genesZoom", "Download Plot"),
                textInput("filename_genesZoom",label = "Save as...",value = "genesPca_zoomed.pdf")))
        ),

        fluidRow(
          column(
            width = 6,
            h4("Profile explorer"),

            checkboxInput("zprofile","Display scaled expression values",value=TRUE),
            plotOutput("genes_profileexplorer"),
            div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                downloadButton("download_genesPca_profile", "Download Plot"),
                textInput("filename_genesPca_profile",label = "Save as...",value = "genesPca_profile.pdf")))
          ,
          column(
            width = 6,
            h4("Boxplot of selected gene"),

            plotOutput("genes_biplot_boxplot"),
            div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                downloadButton("download_genesPca_countsplot", "Download Plot"),
                textInput("filename_genesPca_countsplot",label = "Save as...",value = "genesPca_countsplot.pdf")))
        ),

        fluidRow(
          column(
            width = 6,
            h4("Zoomed heatmap"),
            plotOutput("heatzoom"),
            div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                downloadButton("download_genesHeatmap","Download Plot"),
                textInput("filename_genesHeatmap",label = "Save as...",value = "genesHeatmap.pdf"))),
          column(
            width = 6,
            h4("Zoomed interactive heatmap"),
            fluidRow(radioButtons("heatmap_colv","Cluster samples",choices = list("Yes"=TRUE,"No"=FALSE),selected = TRUE)),
            fluidRow(d3heatmapOutput("heatzoomd3")))),

        hr(),
        box(
          title = "Table export options", status = "primary", solidHeader = TRUE,
          collapsible = TRUE, collapsed = TRUE, width = 12,
          fluidRow(
            column(
              width = 6,
              h4("Points selected by brushing - clicking and dragging:"),
              DT::dataTableOutput("pca_brush_out"),
              downloadButton('downloadData_brush', 'Download brushed points'),
              textInput("brushedPoints_filename","File name...")),
            column(
              width = 6,
              h4("Points selected by clicking:"),
              DT::dataTableOutput("pca_click_out"),
              downloadButton('downloadData_click', 'Download clicked (or nearby) points')),
            textInput("clickedPoints_filename","File name...")
          )
        )
      ),


      tabPanel(
        "Gene Finder",
        icon = icon("crosshairs"),
        fluidRow(
          h1("GeneFinder"),
          wellPanel(width=5,
                    textInput("genefinder",label = "Type in the name of the gene to search",value = NULL),
                    shinyBS::bsTooltip(
                      "genefinder", paste0("Type in the name of the gene to search. If no annotation is ",
                                           "provided, you need to use IDs that are the row names of the ",
                                           "objects you are using - count matrix, SummarizedExperiments ",
                                           "or similar. If an annotation is provided, that also contains ",
                                           "gene symbols or similar, the gene finder tries to find the ",
                                           "name and the ID, and it suggests if some characters are in a ",
                                           "different case"),
                      "right", options = list(container = "body")),
                    checkboxInput("ylimZero","Set y axis limit to 0",value=TRUE),
                    checkboxInput("addsamplelabels","Annotate sample labels to the dots in the plot",value=TRUE)),

          #               fluidRow(
          #                 column(
          #                   width = 6,
          #                   uiOutput("ui_selectID")
          #                 ),
          #                 column(
          #                   width = 6,
          #                   uiOutput("ui_selectName")
          #                 )
          #               ),
          # verbatimTextOutput("debugf"),

          verbatimTextOutput("searchresult"),
          verbatimTextOutput("debuggene"),

          # plotOutput("newgenefinder_plot"),
          column(
            width = 8,
            plotOutput("genefinder_plot"),
            div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                downloadButton("download_genefinder_countsplot", "Download Plot"),
                textInput("filename_genefinder_countsplot",label = "Save as...",value = "pcae_genefinder.pdf"))),
          column(
            width = 4,
            DT::dataTableOutput("genefinder_table"),
            downloadButton("download_genefinder_countstable", "Download Table")

          )
        )
      ),

      tabPanel(
        "PCA2GO",
        icon = icon("magic"),
        h1("pca2go - Functional annotation of Principal Components"),
        h4("Functions enriched in the genes with high loadings on the selected principal components"),
        # verbatimTextOutput("enrichinfo"),
        wellPanel(column(
          width = 6,
          uiOutput("ui_selectspecies")
        ),
        column(
          width = 6,
          uiOutput("ui_inputtype")
        ),

        shinyBS::bsTooltip(
          "ui_selectspecies", paste0("Select the species for the functional enrichment analysis, ",
                                     "choosing among the ones currently supported by limma::goana. ",
                                     "Alternatively, for other species, it can be possible to use one ",
                                     "of the available annotation packages in Bioconductor, and pre-",
                                     "computing the pca2go object in advance"),
          "bottom", options = list(container = "body")),
        verbatimTextOutput("speciespkg"),
        checkboxInput("compact_pca2go","Display compact tables",value=FALSE),
        shinyBS::bsTooltip(
          "compact_pca2go", paste0("Should I display all the columns? If the information content of the ",
                                   "tables is somehow too much for the screen width, as it can be for ",
                                   "objects generated by pca2go with the topGO routines, the app can ",
                                   "display just an essential subset of the columns"),
          "bottom", options = list(container = "body")),

        uiOutput("ui_computePCA2GO"),
        shinyBS::bsTooltip(
          "ui_computePCA2GO", paste0("Compute a pca2go object, using the limma::goana function, ",
                                     "after selecting the species of the experiment under investigation"),
          "bottom", options = list(container = "body"))),

        fluidRow(
          column(width = 3),
          column(
            width = 6,
            DT::dataTableOutput("dt_pcver_pos")),
          column(width = 3)
        ),

        fluidRow(
          column(4,
                 DT::dataTableOutput("dt_pchor_neg")),
          column(4,
                 plotOutput("pca2go")),
          column(4,
                 DT::dataTableOutput("dt_pchor_pos"))
        ),
        fluidRow(
          column(width = 3),
          column(
            width = 6,
            DT::dataTableOutput("dt_pcver_neg")),
          column(width = 3)
        )
      ),



      tabPanel(
        "Multifactor Exploration",
        icon = icon("th-large"),
        h1("Multifactor exploration of datasets with 2 or more experimental factors"),

        verbatimTextOutput("intro_multifac"),

        wellPanel(fluidRow(
          column(
            width = 6,
            uiOutput("covar1")
          ),
          column(
            width = 6,
            uiOutput("covar2")
          )
        ),
        fluidRow(
          column(
            width = 6,
            uiOutput("c1levels")
          ),
          column(
            width = 6,
            uiOutput("c2levels")
          )
        ),
        fluidRow(
          column(
            width = 6,
            uiOutput("colnames1"),
            uiOutput("colnames2")
          )
        ),

        shinyBS::bsTooltip(
          "covar1", paste0("Select the first experimental factor"),
          "bottom", options = list(container = "body")),
        shinyBS::bsTooltip(
          "covar2", paste0("Select the second experimental factor"),
          "bottom", options = list(container = "body")),
        shinyBS::bsTooltip(
          "c1levels", paste0("For factor 1, select two levels to contrast"),
          "bottom", options = list(container = "body")),
        shinyBS::bsTooltip(
          "c2levels", paste0("For factor 2, select two or more levels to contrast"),
          "bottom", options = list(container = "body")),
        shinyBS::bsTooltip(
          "colnames1", paste0("Combine samples belonging to Factor1-Level1 samples for each level in Factor 2"),
          "bottom", options = list(container = "body")),
        shinyBS::bsTooltip(
          "colnames2", paste0("Combine samples belonging to Factor1-Level2 samples for each level in Factor 2"),
          "bottom", options = list(container = "body"))),

        actionButton("composemat","Compose the matrix",icon=icon("spinner"),class = "btn btn-primary"),
        shinyBS::bsTooltip(
          "composemat", paste0("Select first two different experimental factors, for example ",
                               "condition and tissue. For each factor, select two or more ",
                               "levels. The corresponding samples which can be used are then displayed ",
                               "in the select boxes. Select an equal number of samples for each of ",
                               "the levels in factor 1, and then click the button to compute the ",
                               "new matrix which will be used for the visualizations below"),
          "bottom", options = list(container = "body")),

        wellPanel(fluidRow(
          column(4,
                 selectInput('pc_x_multifac', label = 'x-axis PC: ', choices = 1:8,
                             selected = 1)
          ),
          column(4,
                 selectInput('pc_y_multifac', label = 'y-axis PC: ', choices = 1:8,
                             selected = 2)
          ))),

        # fluidRow(verbatimTextOutput("multifacdebug")),
        fluidRow(
          column(6,
                 plotOutput('pcamultifac',brush = 'pcamultifac_brush')),
          column(6,
                 plotOutput("multifaczoom"))
        ),
        fluidRow(downloadButton('downloadData_brush_multifac', 'Download brushed points'),
                 textInput("brushedPoints_filename_multifac","File name..."),
                 DT::dataTableOutput('pcamultifac_out'))

      ),

      tabPanel(
        "Report Editor",
        icon = icon("pencil"),


        fluidRow(
          column(
            width = 6,
            box(
              title = "markdown options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 9,
              radioButtons("rmd_dl_format", label = "Choose Format:", c("HTML" = "html", "R Markdown" = "rmd"), inline = TRUE),
              textInput("report_title", "Title: "),
              textInput("report_author", "Author: "),
              radioButtons("report_toc", "Table of Contents", choices = list("Yes" = "true", "No" = "false")),
              radioButtons("report_ns", "Number sections", choices = list("Yes" = "true", "No" = "false")),
              selectInput("report_theme", "Theme", choices = list("Default" = "default", "Cerulean" = "cerulean",
                                                                  "Journal" = "journal", "Flatly" = "flatly",
                                                                  "Readable" = "readable", "Spacelab" = "spacelab",
                                                                  "United" = "united", "Cosmo" = "cosmo")),
              radioButtons("report_echo", "Echo the commands in the output", choices = list("Yes" = "TRUE", "No" = "FALSE")))),
          column(
            width = 6,
            box(
              title = "editor options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 9,
              checkboxInput("enableAutocomplete", "Enable AutoComplete", TRUE),
              conditionalPanel(
                "input.enableAutocomplete",
                wellPanel(
                  checkboxInput("enableLiveCompletion", "Live auto completion", TRUE),
                  checkboxInput("enableRCompletion", "R code completion", TRUE)
                )
              ),

              selectInput("mode", "Mode: ", choices=modes, selected="markdown"),
              selectInput("theme", "Theme: ", choices=themes, selected="solarized_light"))
          )
          # ,
          # column( # kept for debugging purposes!
          #   width = 6,
          #   verbatimTextOutput("loadedRmd")
          # )
        ),
        fluidRow(
          column(3,
                 actionButton("updatepreview_button", "Update report",class = "btn btn-primary"),p()
          ),
          column(3, downloadButton("saveRmd", "Generate & Save",class = "btn btn-success"))
        ),

        tabBox(
          width = NULL,
          id="report_tabbox",
          tabPanel("Report preview",
                   icon = icon("file-text"),
                   htmlOutput("knitDoc")
          ),

          tabPanel("Edit report",
                   icon = icon("pencil-square-o"),
                   aceEditor("acereport_rmd", mode="markdown",theme = "solarized_light",autoComplete = "live",
                             value=readLines(system.file("extdata", "reportTemplate.Rmd",package = "pcaExplorer")),
                             height="800px"))
        )
      ),

      tabPanel(
        "About", icon = icon("institution"),
        includeMarkdown(system.file("extdata", "about.md",package = "pcaExplorer")),
        hr(),
        #             shiny::verbatimTextOutput("showuploaded1"),
        #             shiny::verbatimTextOutput("showuploaded2"),
        #             shiny::verbatimTextOutput("showuploaded3"),
        #             shiny::verbatimTextOutput("showuploaded4"),

        h4("Session Info"),
        verbatimTextOutput("sessioninfo"),
        footer()
      )

      #           tabPanel(
      #             "Session manager",
      #             ## will put here the things to save/restore the sessions
      #             p("something"),
      #           )
    )
      ),
  skin="blue"
      )
