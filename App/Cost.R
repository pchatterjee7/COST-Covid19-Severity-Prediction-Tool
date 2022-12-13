 # Tutorial with helpful Shiny information: https://shiny.rstudio.com/tutorial/written-tutorial/lesson1/

# Load packages

library(shiny)
library(shinythemes)
#library(flexdashboard)
library(plotly)
library(Seurat)
library(tidyverse)
library(rsconnect)
library(hdf5r)

# Define user interface

# Shiny uses fluidPage to create a display that automatically adjusts to the dimensions 
# of your user’s browser window. 

# To add a widget to your app, place a widget function in sidebarPanel or mainPanel 
# in your ui object.Each widget function requires several arguments. 

# The first two arguments for each widget are 
# 1) a name for the widget: The user will not see this name, but you can use it to access 
# the widget’s value. The name should be a character string.
# 2) a label: This label will appear with the widget in your app. It should be a character string, 
# but it can be an empty string "".

# Valid Shiny themes are: cerulean, cosmo, cyborg, darkly, flatly, journal, lumen, paper, 
# readable, sandstone, simplex, slate, spacelab, superhero, united, yeti

ui <- fluidPage(
  
  theme = shinytheme("cerulean"),
  titlePanel(
    tags$div(class = "jumbotron text-center", style = "margin-bottom:0px;margin-top:0px",
             tags$h2(class = 'jumbotron-heading', style = 'margin-bottom:0px;margin-top:0px'),
             h1(HTML("<font color=#000000>CO</font>VID-19 <font color=#000000>S</font>everity Prediction <font color=#000000>T</font>ool <font color=#000000>(COST)</font>")),
             #'COVID-19 Severity Prediction Tool (COST)'),
             p("An app for triage."),
             p('Created by Mariam Nawaz, Neha Bhatia, and Paramita Chatterjee')
    )),
    sidebarPanel(
      h3("What This Tool Does"),
      p("This tool takes in single-cell RNA-seq data and outputs a COVID-19 severity prediction based on the expression of gene markers identified by Random Forest classifier."),
      h3("How to Use This Tool"),
      p('Upload a .h5 file of a single-cell PBMC sample from a COVID-19 patient and indicate whether the sample has been subsetted for NK cells. If not, no worries! Our tool will do the rest. Press the "predict!" button to determine the expected COVID-19 severity for the sample. Once the pipeline is complete, a button will pop up below your result to download information on the genes we analyzed in your sample.'),
      strong("Note: Your input MUST be a .h5 file of PBMCs or pure NK cells (preferred) from a COVID-19 patient or the tool will not function correctly."),
      br(),
      radioButtons("subset_option", label = h5("Are NK cells subsetted in your input data?"),
                   choices = list("Yes" = 1, "No" = 2), 
                   selected = 1),
      br(),
      fileInput(inputId = "file", "Upload .h5 file", multiple = FALSE, accept = c(".h5")),
      h5("Click the button below to predict COVID-19 severity"),
      actionButton("button", "Predict!")),
    
    
    mainPanel(
      tabsetPanel(type = "tabs",

                  tabPanel("Severity Prediction",
                             plotlyOutput("gauge"), tags$head(tags$style("#result{color: black;
                                 font-size: 20px;
            font-style: bold;
            }")),
                           uiOutput("download"),
                           br(),
                           uiOutput("text1"),
                           uiOutput("result")

                           ),
                  tabPanel("How It Works",
                           br(),
                           p("Acute SARS-CoV-2 infection affects the number of circulating NK cells and their phenotype. Severe COVID-19 results in Natural Killer (NK) cells with high levels of cytotoxic proteins. NK cell dysfunction has been observed in COVID-19 patients."),
                           p("Here, we have created a workflow to predict COVID-19 disease severity based on NK cell transcriptome expression profiles of Single-cell RNA-seq data from pure NK cells or Peripheral Blood Mononuclear Cells (PBMCs). We use the SEURAT pipeline for the analysis. The input file is a H5 matrix, and the output is a COVID-19 severity prediction and .csv file with detailed marker expression information."),
                           img(src="chart-modified.jpeg")),
                  tabPanel("About Us",
                           br(),
                           p("Hi! Welcome to COST (COVID-19 Severity Prediction Tool). We are a team of three students at Georgia Tech who developed this app to help people predict the severity of their COVID-19"),
                           h3("About the Team:"),
                           p(HTML("<strong>Mariam Nawaz</strong> is a M.S. Bioinformatics student at 
                                  the Georgia Institute of Technology. Her current research involves 
                                  using supervised machine learning to help identify biomarkers at the 
                                  transcriptomic level in single-cell data. Her hobbies include 
                                  photography, playing outdoor games, and cooking.")),
                           p(a("LinkedIn", href = "https://www.linkedin.com/in/mariamnawaz1/")),
                           p(a("Bhasin Lab", href = "https://www.bhasinlab.org/")),
                           p(HTML("<strong>Neha Bhatia</strong> is a M.S Bioinformatics 
                           student at the Georgia Institute of Technology. She has a B.S in Biology from Georgia Tech 
                           and has worked in the pharmaceutical and health insurance 
                           industries. She is broadly interested in predictive 
                           healthcare, and her current project integrates ‘omics data 
                           from the Center of Health Discovery and Well Being to profile
                           immune health over time. Following her graduation in spring 2023, she will be 
                           working as a data analyst in Government and Public Services 
                           at Deloitte.")),
                           p(a("LinkedIn", href = "https://www.linkedin.com/in/nehabhatia8/")),
                           p(a("Gibson Lab", href = "https://ggibsongt.wixsite.com/gibsongatech")),
                           p(HTML("<strong>Paramita Chatterjee</strong> has 
                           9 years of genomics research and quality control experience 
                           at Georgia Tech and in the food microbiology industry. 
                           In addition to working in academia, Paramita has worked in 
                           the manufacturing and food manufacturing industries.
                           She has a B.S. in Microbiology, an M.S. in Microbiology, and 
                           an MBA with HR and marketing specialization from Madurai 
                           Kamaraj University, India. Currently, she is a full time 
                           Research Scientist at GT and a PhD student at the 
                           GT-Bioinformatics department. She is interested in 
                           bioinformatics, predictive health genomics and 
                           biomedical health informatics applicable to the cell therapy 
                           research and therapeutic cell manufacturing, as well as 
                           being integrated in the product design and development, 
                           which involves experimental and analytics. She wants to 
                           contribute to research leading to efficient and affordable 
                           next generation personalized cell therapy and precision 
                           medicine.")),
                           p(a("LinkedIn", href = "https://www.linkedin.com/in/paramitachatterjee2022/")),
                           p(a("MC3M Lab", href = "https://cellmanufacturing.gatech.edu/")),
                           p(a("Gibson Lab", href = "https://ggibsongt.wixsite.com/gibsongatech"))
                           )
      ),
      
      
    )
  )
  

# Define server logic

server <- function(input, output) {
  options(shiny.maxRequestSize = 3000*1024^2)
  base_plot <- plot_ly(
    type = "pie",
    rotation = 108,
    direction = "clockwise",
    hole = 0.1,
    textinfo = "label",
    textposition = "outside",
    hoverinfo = "none",
    domain = list(x = c(0, 0.48), y = c(0, 1)),
    marker = list(colors = c('rgb(255, 255, 255)', 'rgb(255, 255, 255)', 'rgb(255, 255, 255)', 'rgb(255, 255, 255)', 'rgb(255, 255, 255)', 'rgb(255, 255, 255)', 'rgb(255, 255, 255)')),
    showlegend = FALSE
  )
  
  reference = readRDS("data/berlin_bonn_merged_mildvssevere_dataset.rds")
  outputs <- eventReactive(input$button, {
    withProgress(message = "Initializing...", value = 0.125, {
    inFile <- input$file
    if (is.null(inFile)){
      return(NULL)
    }
    #1) Read .h5 files
    filepaths <- as.character(input$file$datapath[1])
    incProgress(amount = .125, message = "Reading input file")
    #print(filepaths)
    obj.gse = vector(mode = "list", length = 2)
    obj.gse[[1]] = Read10X_h5(filepaths)
    # Function to check if the input file has Antibody capture data along with the Expression data
    
    check <- function(x) {
      tryCatch(
        expr = {
          obj.gse[[1]] = obj.gse[[1]]$'Gene Expression' # select only gene expression data
        },
        error = function(x) {
          return(obj.gse[[1]])
        } 
      )
    }
    obj.gse[[1]] = check(obj.gse[[1]]) # Calling the check()
    obj.gse[[1]] = CreateSeuratObject(counts = obj.gse[[1]], project = "input_sample", min.cells = 3, min.features = 200)
    #obj.gse[[1]] = obj.gse[[1]]$'Gene Expression' # select only gene expression data
    #obj.gse[[1]] = CreateSeuratObject(counts = obj.gse[[1]], project = "input_sample", min.cells = 3, min.features = 200)
    obj.gse[[1]][["percent.mt"]] <- PercentageFeatureSet(obj.gse[[1]], pattern = "^MT-") # calculate percentage mitochondrial count for each cell in each file simultaneously
    obj.gse[[2]] = reference
    

    #1.1) Subset NK cells
    if (input$subset_option == "2"){
      obj.gse[[1]] = subset(obj.gse[[1]], (NCAM1 > 0) & (GNLY > 0) & (NKG7 > 0) & (KLF2 > 0) , slot = "data")
    }

    #2) Quality Control
    

    #2.2) Filter based on nFeature_RNA,and percent.mt for each dataset
    obj.gse[[1]]= subset(obj.gse[[1]], subset = nFeature_RNA > 200 & nFeature_RNA < 2400 & percent.mt < 20)
    

    #3) Anchor-based integration
    incProgress(amount = .5, message = "Performing anchor-based integration")

    #3.1) Normalize and identify variable features for input dataset

    obj.gse[[1]] <- NormalizeData(obj.gse[[1]])
    obj.gse[[1]] <- FindVariableFeatures(obj.gse[[1]], selection.method = "vst", nfeatures = 2000)
    
    

    # 3.2) Select features that are repeatedly variable across datasets for integration

    features <- SelectIntegrationFeatures(object.list = obj.gse)
  
    #3.3) Perform Integration
    
    tryCatch(FindIntegrationAnchors(object.list = obj.gse, anchor.features = features),
             error = function(e)
               stop("Please make sure your input file contains single-cell RNA-seq data of PBMCs or an NK cell subset for a COVID-19 patient")
)

    anchors = FindIntegrationAnchors(object.list = obj.gse, anchor.features = features)

    #3.4) Create an 'integrated' data assay

    gse <- IntegrateData(anchorset = anchors)

    #3.5) Perform integrated analysis

    #3.5.1) Specify that we will perform downstream analysis on the corrected data. (Note that the original unmodified data still resides in the 'RNA' assay)

    DefaultAssay(gse) <- "integrated"

    #3.5.2) Visualization and clustering
    
    incProgress(amount= 0.70, message = "Performing clustering")
    gse <- ScaleData(gse, features = rownames(gse), verbose = FALSE)
    gse <- RunPCA(gse, verbose = FALSE)

    #3.5.2.2) Clustering

    gse <- RunUMAP(gse, reduction = "pca", dims = 1:30)
    gse <- FindNeighbors(gse, reduction = "pca", dims = 1:30)
    gse <- FindClusters(gse, resolution = 0.5)
    gse$all = "all"
    gse$dataType = "reference"
    gse$dataType[gse$orig.ident == "input_sample"] = "input"


    #4.1) Scale RNA Assay Data for Heatmaps

    gse = ScaleData(gse, features = rownames(gse), assay="RNA", verbose = FALSE)

    #4.3) Calculating expression of Markers in the input dataset
    incProgress(amount = 0.90, message = "Calculating expression of severe and mild COVID-19 markers")})
    Idents(gse) = "dataType"
    #DefaultAssay(gse) = "RNA"

    severe_markers = c("IGKC","S100A8","S100A9","MT-RNR2","RPL10","MT-ND2","EEF1A1","RPLP1","MT-ND4","MT-ND2","AREG","MT-RNR2")
    mild_markers = c("MTRNR2L8","MTRNR2L12","HBB","IFITM1","IFI44L","XAF1","SLFN5","ATP5F1E","MX1","IFITM3","ISG15")
    avg_exp_severe = data.frame(AverageExpression(gse, features = severe_markers, slot = "data", assay = "RNA"))
    avg_exp_mild = data.frame(AverageExpression(gse, features = mild_markers, slot = "data", assay = "RNA"))
    # Adding Gene names as a column
    avg_exp_severe$gene_name = rownames(avg_exp_severe)
    avg_exp_mild$gene_name = rownames(avg_exp_mild)
    
   

    # Create an empty dataframe to store markers and their expressions in the input data that cross the threshold 
    
    sig_severe_markers_in_inputData = data.frame(matrix(ncol = 2, nrow = 0))
    colnames(sig_severe_markers_in_inputData) = c("Gene_Name", "Expression")
    
    sig_mild_markers_in_inputData = data.frame(matrix(ncol = 2, nrow = 0))
    colnames(sig_mild_markers_in_inputData) = c("Gene_Name", "Expression")
    
    #sig_severe_markers_in_inputData$Severity <- "Severe"
    #sig_mild_markers_in_inputData$Severity <- "Mild"
    
    
    # Select Severe Markers present in the input data that cross the threshold
    
    for (row in 1:nrow(avg_exp_severe)) {
      expression = avg_exp_severe[row, "RNA.input"]
      if (expression > 0.5) {
        sig_severe_markers_in_inputData[nrow(sig_severe_markers_in_inputData) + 1,] = 
          c(avg_exp_severe$gene_name[row], avg_exp_severe$RNA.input[row])
      }   
    }
    
    for (row in 1:nrow(avg_exp_mild)) {
      expression = avg_exp_mild[row, "RNA.input"]
      if (expression > 0.5) {
        sig_mild_markers_in_inputData[nrow(sig_mild_markers_in_inputData) + 1,] = 
          c(avg_exp_mild$gene_name[row], avg_exp_mild$RNA.input[row])
      }   
    }
    
    sig_mild_markers_in_inputData$Severity <- "Mild"
    sig_severe_markers_in_inputData$Severity <- "Severe"
    df_download <- rbind(sig_severe_markers_in_inputData, sig_mild_markers_in_inputData)
    
    #4.4) Calculating expression of Markers in the input dataset

    if (nrow(sig_severe_markers_in_inputData) > nrow(sig_mild_markers_in_inputData)) {
      prediction = "Severe"
    } else if (nrow(sig_severe_markers_in_inputData) < nrow(sig_mild_markers_in_inputData)) {
      prediction = "Mild"
    } else {
      prediction = "Unsure"
    }
    
    final_outputs <- list(df_download, prediction, sig_mild_markers_in_inputData$Gene_Name, sep=" ", sig_severe_markers_in_inputData$Gene_Name)
    return(final_outputs)
  })
  
  output$gauge <- renderPlotly({
    if (outputs()[[2]] == "Mild"){
        base_plot <- add_trace(
          base_plot,
          type = "pie",
          values = c(50, 16.67,16.67,16.67),
          labels = c("Predicted COVID-19 Severity", "Mild", "Moderate","Severe"),
          rotation = 90,
          direction = "clockwise",
          hole = 0.3,
          textinfo = "label",
          textposition = "inside",
          hoverinfo = "none",
          domain = list(x = c(0, 0.48), y = c(0, 1)),
          marker = list(colors = c('rgb(255, 255, 255)', 'rgb(68,214,44)', 'rgb(255,222,0)', 'rgb(255,6,0)')),
          showlegend= FALSE
        ) %>% add_annotations(ax=-0.5,
                              ay=0.5,
                              axref='x',
                              ayref='y',
                              x=-0.9,
                              y=0.7,
                              xref='x',
                              yref='y',
                              showarrow=TRUE,
                              arrowhead=3,
                              arrowsize=1,
                              arrowwidth=4,
                              text="") %>% layout(xaxis=list(showgrid = F, range = list(-1,1), visible = F),
                                                  yaxis=list(showgrid = F, range = list(0,1), visible = F))
      
    }
    if (outputs()[2] == "Severe"){
    base_plot <- add_trace(
      base_plot,
      type = "pie",
      values = c(50, 16.67,16.67,16.67),
      labels = c("Predicted COVID-19 Severity", "Mild", "Unknown","Severe"),
      rotation = 90,
      direction = "clockwise",
      hole = 0.3,
      textinfo = "label",
      textposition = "inside",
      hoverinfo = "none",
      domain = list(x = c(0, 0.48), y = c(0, 1)),
      marker = list(colors = c('rgb(255, 255, 255)', 'rgb(68,214,44)', 'rgb(255,222,0)', 'rgb(255,6,0)')),
      showlegend= FALSE
    ) %>% add_annotations(ax=-0.5,
                          ay=0.5,
                          axref='x',
                          ayref='y',
                          x=-0.9,
                          y=0.7,
                          xref='x',
                          yref='y',
                          showarrow=TRUE,
                          arrowhead=3,
                          arrowsize=1,
                          arrowwidth=4,
                          text="") %>% layout(xaxis=list(showgrid = F, range = list(-1,1), visible = F),
                                              yaxis=list(showgrid = F, range = list(0,1), visible = F))
    
    }
    if(outputs()[[2]] == "Unsure"){
      base_plot <- add_trace(
        base_plot,
        type = "pie",
        values = c(50, 16.67,16.67,16.67),
        labels = c("Predicted COVID-19 Severity", "Mild", "Unknown","Severe"),
        rotation = 90,
        direction = "clockwise",
        hole = 0.3,
        textinfo = "label",
        textposition = "inside",
        hoverinfo = "none",
        domain = list(x = c(0, 0.48), y = c(0, 1)),
        marker = list(colors = c('rgb(255, 255, 255)', 'rgb(68,214,44)', 'rgb(255,222,0)', 'rgb(255,6,0)')),
        showlegend= FALSE
      ) %>% add_annotations(ax=-0.5,
                            ay=0.5,
                            axref='x',
                            ayref='y',
                            x=-0.5,
                            y=0.8,
                            xref='x',
                            yref='y',
                            showarrow=TRUE,
                            arrowhead=3,
                            arrowsize=1,
                            arrowwidth=4,
                            text="") %>% layout(xaxis=list(showgrid = F, range = list(-1,1), visible = F),
                                                yaxis=list(showgrid = F, range = list(0,1), visible = F))

    }
    base_plot
  })
  
  output$text1 <- renderUI({
    if(outputs()[[2]] == "Mild"){
      str1 <- "You have above threshold expression of the following genes that are indicative of mild COVID-19:"
    }
    if(outputs()[[2]] == "Severe"){
      str1 <- "You have above threshold expression of the following genes that are indicative of severe COVID-19:"
    }
    if(outputs()[[2]] == "Unsure"){
      str1 <- "You do not have above threshold expression of gene expressing either severe or mild COVID-19"
    }
    
    HTML(paste(str1))
  })

  output$result <- renderUI({
    if(outputs()[[2]] == "Mild"){
      
      str2 <- paste0(outputs()[[3]], sep="<br/>")
    }
    if(outputs()[[2]] == "Severe"){
      str2 <- paste0(outputs()[[4]], sep="<br/>")
    }
    
    HTML(paste(str2))
})

  
  output$download <- renderUI({
    req(input$file, input$button, outputs())
    downloadButton('download_item', label = 'Download gene expression information (.csv)') })

  output$download_item <- downloadHandler(
    filename = function(file) {
      paste0("COVID_severity_expression", Sys.Date(), ".csv")
    },
    content = function(con) {
      write.csv(outputs()[[1]], con, row.names = TRUE)
    }
  )
  
}

# Run the app

shinyApp(ui = ui, server = server)
