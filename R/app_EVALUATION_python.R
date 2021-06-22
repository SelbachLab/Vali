library(shiny)
library(RSQLite)
library(pracma)
library(data.table)
library(parallel)
library(gtools)
library(h2o)# Machine Learning module, for Scoring
library(rawDiag)# RawFile Reader https://github.com/fgcz/rawDiag
library(Peptides)# for MS1, currently not implemented
library(enviPat)# for isotope pattern distribution calculation used in MS1 assessment
library(shinyWidgets)
require(rawDiag)
require(compiler)
library(tcltk)

# library(compiler)# 
# library(gplots)
# library(drc)
# library(shinydashboard)


options(warn=-1)
# library(gridExtra)
#humpe
args1 = commandArgs(trailingOnly=TRUE)
use_ggplot <- F
plotheight=400

inactivity <- "function idleTimer() {
  var t = setTimeout(logout, 500000);
window.onmousemove = resetTimer; // catches mouse movements
window.onmousedown = resetTimer; // catches mouse movements
window.onclick = resetTimer;     // catches mouse clicks
window.onscroll = resetTimer;    // catches scrolling
window.onkeypress = resetTimer;  //catches keyboard actions

function logout() {
window.close();  //close the window 
}

function resetTimer() {
clearTimeout(t);
t = setTimeout(logout, 500000);  // time is in milliseconds (1000 is 1 second)
}
}
idleTimer();"

version = "v0.1.002"
print("Version")
#install.packages("lsa")
minthresh = 5 # minimal amount of datapoints required over whole run per species per rawfile
if(length(args1) >0){
  print(args1)
  #SystemPath <- dirname(dirname(args[1]))
  # print(SystemPath)
  
  
}
print("Evaluating Path")
if(!exists("SystemPath")){

  try(SystemPath <- dirname(dirname(rstudioapi::documentPath())))
  SystemPath <- "/Users/henno/Dropbox/RPackages/Vali_git"
  
}
print("Evaluating Path 2")

SystemPath <- gsub("\\\\","/",SystemPath)
SystemPath <- gsub("/$","",SystemPath)

SystemPath <<- SystemPath


pythonpath <- c(paste(SystemPath,"PortablePython2.7.6.1/App/python.exe",sep = "/"),paste(SystemPath,"python-3.6.6.amd64/python.exe",sep = "/"))



# StandardSettings: (Needs an Update) 
inputListStdSet <- list()
inputListStdSet$PeakElution <- "NA"
inputListStdSet$RetentionTimeWindow <- 10
inputListStdSet$Match.Count <- 3
inputListStdSet$p.value <- 0.01
inputListStdSet$PeakWidth <- 1
inputListStdSet$IncludeUnsignificantPeaks <-F
inputListStdSet$ValuesAcrossSampleThreshold <- 2
enableBookmarking("server")

sy <- Sys.info()
PresetPaths <- F
if(PresetPaths){
  if(any(grep("windows",Sys.info(),ignore.case = T))){
    mainPath  = "E:/Projects/Picky_PRM_Test/PRManalyzerTest/"
    maxquant = mainPath
    inclusionList = "E:/Projects/Picky_PRM_Test/PRMUPS_TRUE/"
  }else{
    mainPath <- ("/Users/henno/Documents/Skr?ipte/R-Functions/Selbach-Functions/DataAnalysis/Picky_Analyzer/TestRun")
    maxquant <- "/Users/henno/Documents/Skripte/R-Functions/Selbach-Functions/DataAnalysis/Miffy/PRM_TEST_YEAST_2/PRM_T"
    inclusionList <- "/Users/henno/Documents/Skripte/R-Functions/Selbach-Functions/DataAnalysis/Miffy/PRM_TEST_YEAST_2/PRMUPS_TRUE"
  }
}else{
  mainPath <- "/~"
  maxquant <- "/~"
  inclusionList <- "/~"
  
}

print("SystemPath:")
print(SystemPath)


sourcescript = paste(SystemPath,"/R/EvaluationScript_PRM_sqlite.R",sep = "")
Py_MSFR_Module = c(paste(SystemPath,"/Py/MSFileReader_PRM_2.py",sep = ""),paste(SystemPath,"/Py/MSFileReader_PRM_3.py",sep = "")) # Newer with random matches but not so well working

# Setting Paths, reading functions
print(getwd())
print(paste("Loading SourceScript",sourcescript))
# if(!file.exists(sourcescript)){
#   setwd("../")
#   if(!file.exists(sourcescript)){
#     print("Cannot find EvaluationScript_PRM_sqlite. ")
#     stop()
#   }
# }
# mainPath <- "./"
#setwd(mainPath)
source(sourcescript)
print("Loaded SourceScript")

# Start of Shiny Script 

# ui ------ 

ui <- fluidPage(
  tags$head(tags$style(
    HTML('#title {
           color: black;
           font-size: 40px;
           font-style: bold;
          }'))),
  #tags$head(tags$style(".shiny-progress {top: 50% !important;left: 50% !important;margin-top: -100px !important;margin-left: -250px !important; color: blue;font-size: 20px;font-style: italic;}")),
  tags$script(inactivity),   
  
  fluidRow(column(1,icon("jedi","fa-4x")),column(9, titlePanel(paste("Vali - PRM Validation Tool"),version) ),column(2, actionButton("Export", "Export Tables!",style="margin-top:25px;margin-right:5px",icon = icon("file-export"))) ),
  
  fluidRow(column(1, actionButton("View", icon = icon("eye"),label = ""),style = "margin-top: 25px"),
           column(7,uiOutput("Sessions")),
           column(2,textInput("SearchSessions",label = "Filter",value = "")),
           column(2, actionButton("RemoveEntry",label = "Remove Analysis",icon=icon("trash"),width = "100%"),style = "margin-top: 25px")
  ),
  # Paths
  mainPanel(
    tabsetPanel(
      tabPanel("Validation Tool",icon=icon("binoculars") ,
               
               wellPanel(
                 #ROW 2------
                 #fluidRow(column(7,uiOutput("precursors")),         
                 fluidRow(column(7,uiOutput("pc1")),
                          column(3,selectInput("sortType","Sorting",c("FDR","DL_Scores","RF_Scores","Count","SCA","Rating","mz"),selected = "DL_Scores")),
                          column(2,prettyRadioButtons("Decreasing","+/-",c("+","-"),selected = "-",inline = T))
                 ),

                 
                 # chromatogram --------- plot
                 wellPanel(style = "background-color: #ffffff;",
                           fluidRow(
                             column(3,
                                    # mainPanel(h5("Precursor Rating"),width = "100%"),
                                    wellPanel(
                                    fluidRow(column(12,textOutput("PrecursorQuality"),style="margin-bottom:5px")),
                                    fluidRow(column(4,actionButton("BadPeak","",icon = icon("thumbs-down"),style="color: #ff2600; border-color: #ff2600")),
                                             column(4,actionButton("notdefinedPeak","",icon = icon("stethoscope"),style="color: #2e6da4; border-color: #2e6da4")),
                                             column(4,actionButton("GoodPeak","",icon = icon("thumbs-up"),style="color: #00b19c; border-color: #00b19c"))))     ,                       

                                    # checkboxInput("PRMonly","PRM",value = T),
                                    conditionalPanel("false",
                                                     checkboxInput("SignificantOnly","SignificantOnly",value = F)
                                    ),
                                    wellPanel(
                                      fluidRow(column(4,mainPanel(h5("Peak Assignment"),width = "100%"))),#column(8,checkboxInput("allshifts","All m/z",T))),
                                      
                                      fluidRow(column(4,actionButton("ManualSetPeak","Set",style="background-color: #00b19c;color: #fff")),column(8,actionButton("RemovePeak","Remove",width = "100%",style="background-color: #ff2600;color: #fff"))),
                                      conditionalPanel("false",
                                                       fluidRow(column(12,actionButton("SetSearchArea","Area Search",width = "100%",style="background-color: #2e6da4;color: #fff"),style = "margin-top: 5px"))
                                                       
                                                       ),
                                      fluidRow(column(6,actionButton("Requantify","Requantify"),style = "margin-top: 10px")),
                                      fluidRow(column(12,selectInput("Requantify_Priority","Peak Priority",c("DL_Scores","Intensity","Light","Heavy","none")),style = "margin-top: 10px")),
                                      fluidRow(column(6,actionButton("reset", "Reset Assignments"),style = "margin-top: 10px")),
                                      fluidRow(column(12,switchInput("supersmooth_I_set","Smoother",value = F),style="margin-top:20px")),
                                      fluidRow(column(12,switchInput("CenterPeak",label = "Center Peak",width="auto")))
                                    ),
                                    wellPanel(
                                      fluidRow(conditionalPanel("false",column(12,switchInput("Align",label = "RT Alignment",FALSE)))
                                               ),
                                      
                                      fluidRow(column(12,selectInput("secPlotType","Added Info",    c("FDR",
                                                                                                      "RF_Scores",
                                                                                                      "DL_Scores",
                                                                                                      # "GLM_Scores",
                                                                                                      # "GBM_Scores",
                                                                                                      "SCAall",
                                                                                                      "SCAcut",
                                                                                                      "MatchCount"#,
                                                                                                      #"CorrelationScoreSumSpecial2"#,
                                                                                                      # "DiffSum",
                                                                                                      # "Peaks1",
                                                                                                      # "Peaks2"
                                                                                                      
                                      ),selected = "DL_Scores")))
                                    ,
                                      fluidRow(column(12,prettyRadioButtons("MStype","MS Level",c("MS1","MS2"),"MS2",inline=T)),
                                              column(12,uiOutput("FDRpep"))),
                                    fluidRow(column(12,actionButton("ExcludePrecursor","Exclude Precursor")))
                                    
                                    
                                    )
                                  
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    # sliderInput(inputId = "CorrelationCutoff",# id for identifiying the input later on
                                    #             label = "Minimal accepted Pearson's R", # Title or Explanation
                                    #             value = 0.5,min = 0.3,max = 1 ,step = 0.05# input specific arguments
                                    # ),
                                    # sliderInput(inputId = "fdrset",# id for identifiying the input later on
                                    #             label = "FDR Correlation", # Title or Explanation
                                    #             value = 0.01,min = 0,max = 0.4 ,step = 0.001# input specific arguments
                                    # )
                             ),
                             column(9,
                                    wellPanel(
                                    fluidRow(
                                      column(10,uiOutput("rawfile")),
                                      column(2,switchInput("allRawFiles","all",FALSE),style="margin-top:25px")),
                                  
                                    
                                    fluidRow(                              
                                      column(3,uiOutput("pc2")),
                                      column(2,switchInput("SimultaneousMassShift","all",value = T),style = "margin-top: 25px"),
                                      column(1,actionButton("Minus","", icon = icon("minus"),style="color: #fff; background-color: #D3D3D3; border-color: #fff"),style = "margin-top: 25px"),
                                      column(1,actionButton("Plus","", icon= icon("plus"),style="color: #fff; background-color: #D3D3D3; border-color: #fff"),style = "margin-top: 25px"),
                                      column(1,actionButton("left","", icon= icon("arrow-left"),style="color: #fff; background-color: #888800; border-color: #fff"),style = "margin-top: 25px"),
                                      column(1,actionButton("right","", icon= icon("arrow-right"),style="color: #fff; background-color: #888800; border-color: #fff"),style = "margin-top: 25px"),
                                      column(1,actionButton("vertical","", icon= icon("arrows-v"),style="color: #fff; background-color: #008800; border-color: #fff"),style = "margin-top: 25px"),
                                      column(1,actionButton("horiz","", icon= icon("arrows-h"),style="color: #fff; background-color: #008800; border-color: #fff"),style = "margin-top: 25px"),
                                      column(1,actionButton("height","", icon= icon("arrow-up"),style="color: #fff; background-color: #000088; border-color: #fff"),style = "margin-top: 25px"),
                                      
                                      
                                      
                                    )
                                    
                                    
                                    ,
                                    
                                    plotOutput("PlotOutput", height = plotheight,
                                               dblclick = "plot1_dblclick",
                                               brush = brushOpts(
                                                 id = "plot1_brush",
                                                 clip = F,
                                                 resetOnNew = TRUE
                                               )
                                    )),
                                    wellPanel(
                                      
                                      fluidRow(
                                        column(8,uiOutput("selected_transitions")),
                                        column(1,actionButton("all_Transitions",label = "all"),style = "margin-top: 25px"),
                                        column(1,checkboxInput("nomodifications","no mods",value = FALSE),style = "margin-top: 25px"),
                                        column(1,checkboxInput("top5","top10",value = FALSE),style = "margin-top: 25px"),
                                        tags$style(type='text/css', ".selectize-input { font-size: 12px; line-height: 12px;} .selectize-dropdown { font-size: 12px; line-height: 12px; }")
                                        
                                        
                                      )
                                      
                                    )
                                    
                                    
                                    ,
                                    )
                             #,
                             # column(width = 4, class = "well",
                             #        h4("Brush and double-click to zoom"),
                             #        plotOutput("plot1", height = 300,
                             #                   dblclick = "plot1_dblclick",
                             #                   brush = brushOpts(
                             #                     id = "plot1_brush",
                             #                     resetOnNew = TRUE
                             #                   )
                             #        ))
                             
                             
                           ))),
               #               wellPanel(
               # # PeptidePlot ------
               # 
               #               ),
               
               wellPanel(
                 #------ ROW1
                 fluidRow(
                   # column(6,  uiOutput("precursors"),uiOutput("rawfile"),
                   #        uiOutput("PeakElutionSet")
                   #        
                   #       
                   #               
                   #               
                   # ),#,
                   # column(6,  )
                   
                   
                   # column(2,actionButton("SaveSettings", "Save Sets")),
                   column(3, 
                          # actionButton("PeaksReset","Reset Peaks"),
                          conditionalPanel("false",
                                           selectInput("plottype",NULL,c("Transition Plot","Volcano"),selected = "Transition Plot")
                                           
                          )
                          
                   )
                   
                   
                   
                   
                 )
               )  
               
      ),
      # tabPanel("Quantitative ReadOut",
      #          wellPanel(style = "background-color: #ffffff;",
      #                    fluidRow(
      #                      column(5,
      #                             # checkboxInput("IncludeUnsignificantPeaks","Include Unsignificant Peaks",value = F),
      #                             # sliderInput(inputId = "ValuesAcrossSampleThreshold",label = "ValuesAcrossSampleThreshold",min = 0,max = 40,step = 1,value = inputListStdSet$ValuesAcrossSampleThreshold),
      #                             # # sliderInput(inputId = "Match.Count",label = "Match Count",min = 0,max = 40,step = 1,value = inputListStdSet$Match.Count),
      #                             # checkboxInput("Normalize","Normalize (FOT)",value = F),
      #                             # checkboxInput("log10","log10",value = T),
      #                             
      #                             selectInput("quantitationType","Quantitation Type",c("XIC","Intensities"),selected = c("XIC"))
      #                             
      #                             
      #                      ),
      #                      column(7,plotOutput("PeptidePlot",height = 600))
      #                      
      #                    )
      #                    
      #                    
      #                    
      #                    
      #                    
      #                    # RetentionTimeWindow = input$RetentionTimeWindow,alpha = input$p.value,MatchCount = input$Match.Count,PeakWidth = input$PeakWidth
      #                    
      #          )
      # ),
      # TAB PANEl ------
      tabPanel("Rawfile Scanner",icon=icon("folder-open"),
               wellPanel(
              wellPanel(
                 fluidRow(column(2,actionButton("mainPathButton","Browse",icon=icon("folder"),style="margin-top:25px")),column(10,textInput(inputId = "mainPath","Main Folder",value = mainPath,width = '100%'))) ,# Protein Selecter
                 conditionalPanel("false",
                 fluidRow(column(2,actionButton("MaxQuantButton","Browse",icon=icon("folder"),style="margin-top:25px")),column(10,textInput(inputId = "MaxQuant","MaxQuant txt",value = maxquant,width = '100%')))
                 ),# Protein Selecter
                 fluidRow(column(2,actionButton("inclusionListButton","Browse",icon=icon("folder"),style="margin-top:25px")),column(10,textInput(inputId = "inclusionList","Picky Export Folder",value = inclusionList,width = '100%'))) ,# Protein Selecter
                 
                 actionButton("saveSessionPaths",icon = icon("bed"),label = "Save Paths",style="margin-bottom=25px")),
                 # fluidRow(column(1)),
                 wellPanel(
                   fluidRow(
                     column(2,checkboxInput("DIA","DIA",FALSE)),
                     
                     column(3,conditionalPanel("input.DIA==false",
                                               sliderInput("ppm1","MS1 ppm",0,50,10)
                                               
                     )
                     
                     ),
                     column(3,
                            sliderInput("ppm2","MS2 ppm",0,50,10)
                            
                     ),
                     column(2,
                            numericInput(inputId = "Threads",label = "Threads",min = 1,max = 10,value = 1)),
                     column(2,                 checkboxInput("Recompile","Recompile Database",value = FALSE))
                     
                   )
                 )
                ,
                 actionButton("goButton", "Go!",width = "100%",icon = icon("bullhorn")),
                 #,
                 #actionButton("savepaths", "Only Save Path",width = "100%"),
                 
                 
                 textOutput("analyzedListPath")
               )
               
               
      ),
      tabPanel("Advanced Settings",icon=icon("skull-crossbones"),
               # sliderInput(inputId = "PeakWidth",label = "Peak Detection Width",min = 0,max = 40,step = 1,value = inputListStdSet$PeakWidth),
               conditionalPanel("false",
                                selectInput("quantitationType","Quantitation Type",c("XIC","Intensities"),selected = c("XIC"))
                                ),
               wellPanel(
                 switchInput("ApplyMaximumWidth","ApplyMaximumWidth"),
                 numericInput("MinPeakWidth","MinPeakWidth [min]",min=0,value=0.16),
                 numericInput("MaxPeakWidth","MaxPeakWidth [min]",min=0,value=1/60*30)
               )
               ,
               switchInput("EnableNeighborZeroImputation","NeighbourZeroImputation",value = T),
               sliderInput(inputId = "RetentionTimeWindow",label = "Peak Window",min = 0,max = 10,step = 0.1,value = inputListStdSet$RetentionTimeWindow),
               wellPanel(
                 numericInput("supersmooth_bw_set","Smoother span",min=0.001,value = 0.1)
                 
                 
               )
               
               
      )
    ),  style='width: 100%; height: 100%'
  )
  
  
)


# Title
# tableOutput("FilteredTable")
# verbatimTextOutput("ret")

#server ------
server <- function(input, output, session){
  ## Removes Peaks from a m/z
  observeEvent(c(input$reset,input$Requantify_Priority),{
    print("Reseting")
    DBL <<- dblistT(dbpath())
    db <- dbConnect(SQLite(),dbpath())
    dbname <- dbtaName(AnalyzedTransitions(),dbpath())
    sapply(dbname,function(na){
      gr <- any(DBL==na)
      if(gr){
          print(paste("Removing",na))
          dbRemoveTable(db,na)
      
      }else{
        print(na,"not found.")
      }
    })
    
    dbDisconnect(db)
    rm(db)
  })
  # Vali Database Path
  dbpath <- reactive({
    validate(need(file.exists(input$mainPath),"No directory available."))
    setwd(input$mainPath);paste(TransitionListsDir(),"PickyAnalyzer.sqlite",sep = "")
    })
  
  # Sessions Settings #######
  # library for reanalysis:
  ST <- reactive({
    
    try({fread(paste(input$inclusionList,"SpectraTable.txt",sep="/"),sep = "\t")})
  })
  
  TT <- reactive({try({fread(paste(input$inclusionList,"TransitionTable.txt",sep="/"),sep="\t")})})
  Models <- reactive({
    Models.h2o <- NULL
    if(require(h2o)){
      try({
        list.files(paste(SystemPath,"Model/MojoBackup",sep = "/"),pattern = "model_",full.names =T)
        
        h2o.init(nthreads = threads, max_mem_size = '2g', ip = "127.0.0.1", port = 4321)
        # 
        # # localH2O <- h2o.init(nthreads = threads)
        
        # Models.h2o <- lapply(mdpath,function(x){
        #   try(h2o.loadModel(x))
        # })
        Models.h2o <- lapply(mdpath,function(x){
          try(h2o.import_mojo(x))
        })
        names(Models.h2o) <- strsplitslot(basename(mdpath),1,"_")
        Models.h2o <- Models.h2o[sapply(Models.h2o,class)!="try-error"]
      })
      
    }
    Models.h2o
  })
  
  ## Opens Folder from current Session
  observeEvent(input$View,{
    path <- input$mainPath
    if(file.exists(path)){
      system(paste("open",path))
      
    }
  })
  ## Reads Sessions
  Paths <- reactive({
    pa <- paste(SystemPath,"Sessions",sep = "/")
    Paths <- list(mainPath = "no sessions detected")
    temp <- Paths
    if(file.exists(pa)){
      temp2 <- list.files(pa,full.names = T)
      validate(need(length(temp2)>0,"No sessions found. You can initiate a new analysis in 'the RawFile Scanner' Tab."))
      if(length(temp2) > 0){
        temp <- lapply(temp2,function(x){
          load(x)
          Li$SessionPath <- x
          Li <- Li
        })
        names(temp) <- temp2
      }
      
      
    }
    
    return(temp)
    
  })
  ## UI output
  SessionChoices <- reactive({
    input$saveSessionPaths
    Paths <- sapply(Paths(),function(x){x$mainPath})
    Paths <- unique(Paths)
    Paths <- Paths[file.exists(Paths)]
    validate(need(length(Paths) > 0,"No Sessions available"))
    Paths <- Paths[order(file.info(Paths)$mtime,decreasing = T)]
    if(input$SearchSessions!=""){
      PathsPeu <- grep(input$SearchSessions,Paths,value = T,ignore.case = T)
      if(length(PathsPeu)>0){
        Paths <- PathsPeu
      }
    }
    
    Paths
  })
  output$Sessions <- renderUI({
    
    
    selectInput(inputId = "Sessions",label = "Loaded Session",choices = SessionChoices(), selected = SessionChoices()[1],width = "100%")
    
  })
  ## Load Session and Settings 
  observe({
    pa <- input$Sessions
    sessions <- Paths()
    
    if(length(pa) > 0){
      temp <- lapply(sessions,function(x){
        x$mainPath
      })
      tempFileInfo <- sapply(temp,function(x){file.info(x)$mtime})
      ul <- names(temp[unlist(duplicated(temp))])
      
      # sapply(ul,unlink)
      #temp <- temp[!unlist(duplicated(temp))]
      
      selsessions <-sessions[temp == pa]
      selsessions <- selsessions[[1]]
      updateTextInput(session,"mainPath",value = selsessions$mainPath)
      updateTextInput(session,"MaxQuant",value = selsessions$MaxQuant)
      updateTextInput(session,"inclusionList",value = selsessions$inclusionList)
      updateSliderInput(session,"ppm",value = selsessions$ppm)
      updateSliderInput(session,"ppm2",value = selsessions$ppm2)
      updateNumericInput(session,"Threads",value = selsessions$threads)
    }
    
    
    
  })
  ## Remove Session
  observeEvent(input$RemoveEntry,{
    isfun <- input$Sessions
    PAS <- Paths()
    # 
    sapply(PAS,function(x){
      if(x$mainPath==isfun){
        x <<- x
        unlink(x$SessionPath)
      }
    })
    updateSelectInput(session,"Sessions",choices = c(SessionChoices()))
    
  })
  # Unknown Settings ##########
  observeEvent(input$ExcludePrecursor,{
    
    if(any(dblistT(db = dbpath())== "ExcludedPrecursors")){
      append = T
    }else{
      append = F
    }
    precu <- precursors.i()
    precu <- data.frame(precursorID = precu)
    # print(precu)
    # print(dbpath())
    
    dbwrite(x = precu, db =dbpath(),name = "ExcludedPrecursors",append = append)
    print("Excluded Precursor")
    precursors()
    precursorRemove()
  })
  precursorRemove <- reactive({
    if(any(dblistT(db=dbpath())=="ExcludedPrecursors")){
      cl <- class(try(EP <- dbread("ExcludedPrecursors",dbpath())))
      return(EP)
    }else{return(NULL)}
    
    
  })
  
  #### UI ------
  
  # Precursors for Precursor List
  precursors <- reactive({
    # pf<- dbReadTable(db,"PrecursorSelectionTable")
    er <- try({
      # dbps <<- dbpath()
      # wds <<- getwd()
      # dbReadTable(db,"PrecursorSelectionTable")
      PrecursorSelectionTable <- dbread("PrecursorSelectionTable",dbpath())
      PrecursorSelectionTable <- PrecursorSelectionTable[PrecursorSelectionTable$Dim >= minthresh,]
      PrecursorSelectionTableS <- PrecursorSelectionTable
      Precursors <- unique(PrecursorSelectionTableS$Species)
      
      PrecursorsT <- strsplit(Precursors,"_",fixed = T)
      PrecursorsT <- lapply(PrecursorsT,function(x){
        x <<- x
        if(length(x)> 6){
          x[6] <- paste(x[6:length(x)],collapse = "_")
          x <- x[1:6]
        }
        if(length(x)< 6){
          x <- c(x,rep(NA,6-length(x)))
        }
        return(x)
      })
      
      
      PrecursorsT <- t(as.data.frame(PrecursorsT,header = F,row.names = NULL,check.names = F))
      colnames(PrecursorsT) <- c("mz","Charge","Sequence","Accession","GeneSymbol","misc")
      
      rownames(PrecursorsT) <- Precursors
      PrecursorSelectionTableS <- unique(PrecursorSelectionTableS)
      # PrecursorSelectionTableS <- data.table(PrecursorSelectionTableS)
      # PrecursorSelectionTableS[,{
      #   temp <<- .SD
      #   list(max(FDR))
      # },.(rf,Species)]
      try({PrecursorSelectionTableS <-aggregate(PrecursorSelectionTableS[,-match(c("rf","Species"),colnames(PrecursorSelectionTableS))],list(Species = PrecursorSelectionTableS$Species,rf = PrecursorSelectionTableS$rf),function(x){list(min =min(as.numeric(x),na.rm= T),max=max(as.numeric(x),na.rm = T),Stuff=paste(sort(unique(x))),collapse = "|")})},silent = T)
      
      PrecursorsT <- PrecursorsT[match(PrecursorSelectionTableS$Species,rownames(PrecursorsT)),]
      Species <- rownames(PrecursorsT)
      rownames(PrecursorsT) <- NULL
      # Precursors <<- Precursors 
      PrecursorsT <- as.data.frame(PrecursorsT,make.names = F)
      PrecursorsT$Species <- Species
      # PrecursorsT <- unique(PrecursorsT)
      
      # PrecursorsT <- PrecursorsT,PrecursorSelectionTableS
      
      seSplit <- strsplit(as.character(PrecursorsT$Sequence),"#")
      Label <- sapply(seSplit,function(x){
        if(length(x)==1){
          return("")
        } 
        if(length(x) == 2){
          return(x[2])
        }
      })
      Label[lengths(Label)==0]<-""
      se <- sapply(seSplit,function(x){x[1]})
      PrecursorsT$Sequence <- se
      PrecursorsT$Label <- sapply(Label,function(x){return(x[1])})
      
      
      #######
      # PrecursorsT$Species <- NULL
      # rownames(PrecursorsT) <- NULL
      
      
      
      
      PrecursorSelectionTableS <- PrecursorSelectionTableS[order(as.character(PrecursorsT$Sequence)),]
      PrecursorsT <- PrecursorsT[order(as.character(PrecursorsT$Sequence)),]
      
      PrecursorSelectionTableS <- PrecursorSelectionTableS[order(as.character(PrecursorsT$Accession)),]
      PrecursorsT <- PrecursorsT[order(as.character(PrecursorsT$Accession)),]
      
      PrecursorSelectionTableS <- PrecursorSelectionTableS[order(as.character(PrecursorsT$GeneSymbol)),]
      PrecursorsT <- PrecursorsT[order(as.character(PrecursorsT$GeneSymbol)),]
      
      
      ordervec <- NULL
      decreasingType <- input$Decreasing == "-"
      sortType <- input$sortType
      if(sortType == "FDR"){
        ordervec <- as.numeric(PrecursorSelectionTableS$FDR[,1])
      }
      if(sortType == "DL_Scores"){
        ordervec <- as.numeric(PrecursorSelectionTableS$DL_Scores[,2])
      }
      if(sortType == "Count"){
        ordervec <- as.numeric(PrecursorSelectionTableS$Count[,2])
      }
      if(sortType == "SCA"){
        ordervec <- as.numeric(PrecursorSelectionTableS$SCA[,2])
      }
      if(sortType == "Rating"){
        ordervec <- as.character(PrecursorSelectionTableS$Rating[,3])
      }
      if(sortType == "mz"){
        ordervec <- as.numeric(gsub("mz","",PrecursorsT$mz[,1]))
      }
      ordervec <<- ordervec
      PrecursorsT <<- PrecursorsT
      if(length(ordervec)>0){
        PrecursorsT <- PrecursorsT[order(ordervec,decreasing = decreasingType),]
      }
      
      
      # print("precursorRemove()")
      # print(precursorRemove())
      if(length(precursorRemove())> 0){
        preR <<- precursorRemove()
        # PrecursorsT1 <<- PrecursorsT
        MEX <- is.na(match(rownames(PrecursorsT),preR$precursorID))
        PrecursorsT <- PrecursorsT[MEX,]
      }
      PrecursorsT <- PrecursorsT[!is.na(PrecursorsT$mz),]
    },silent = T)
    if(class(er)=="try-error"){
      validate(
        need(class(er)!="try-error","Error in PrecursorSelectionTable")
      )
      PrecursorsT <- NULL
    }
    return(PrecursorsT)
    
  })
  ##  Sequences for Ui from precursors()
  output$pc1 <- renderUI({

    precursorTable <- precursors()
    validate(
      need(is.data.frame(precursorTable),"no Sequences")
    )
    
    precursorTable <- precursorTable[!is.na(precursorTable$mz),]
    if(length(precursorTable$misc)==0){
      precursorTable$misc <- 1:dim(precursorTable)[1]
    }
    ##mancage Order
    sevec <- c((paste(as.character(precursorTable$Sequence),
                      as.character(precursorTable$GeneSymbol)#,
                      # as.character(precursorTable$GeneSymbol)
    )))
    sevec <- (unlist(sevec[sevec!=""]))
    sevec <- sevec[!duplicated(sevec)]
    
    selectInput("sequence", "Sequence", 
                choices  = sevec,selected = sevec[1])
  })
  ##  mz for Ui from precursors()
  mz <- reactive({
    # STi <<- ST()
    # TTi <<- TT()
    # mdi <<- Models()
    
    validate(
      need(is.character(input$sequence),"sequence is not a character")
    )
    se <- input$sequence
    if(length(se)==0){
      return(NULL)
    }
    se <- sapply(strsplit(input$sequence," "),function(x){x[1]})
    mz <- precursors()$mz[as.character(precursors()$Seq)==se]
    return(unique(mz))
  })
  ##  mz for Ui from precursors()
  output$pc2 <- renderUI({
    choicesvec <- gsub("mz","",unique(mz()))
    
    selectInput("mz", "m/z", 
                choices  = sort(choicesvec),multiple = input$SimultaneousMassShift,selected = choicesvec)
  })
  ## precursor Selection
  precursors.i <- reactive({
    
    prec <<- precursors()
    mztest<<- mz()
    
    validate(
      need(is.character(input$sequence),"sequence is not a character"),
      need(is.data.frame(precursors()),"Precursors is not a character"),
      # need(all(is.factor(mz())),"Factor Problem in mz()"),
      need(all(is.character(input$mz)),"Character Problem in mz()")
      
    )
    
    se <- sapply(strsplit(input$sequence," "),function(x){x[1]})
    if(length(se) == 0){
      return(se)
    }
    se.i <- se
    
    if(input$SimultaneousMassShift){
      se <<- se
      pr <<- precursors()
      mzi <<- input$mz
      # mza <- mz()
      # mzii <<- mzi
      # mzai <<- mza
      # pr.i <<- pr
      mza <<- paste("mz",mzi,sep = "")
      selectedInfo <- apply(pr,1,function(x){
        x <- x
        any(x[names(x) == "mz"]==mza)&x[names(x) == "Sequence"]==se
      })
      validate(need(any(selectedInfo),"no matching mz"))
      # selectedInfo <- sapply(pr$mz,function(x){
      #   any(mza==x) &
      #   which(mza == x & pr$Sequence == se)
      # })
      # selectedInfo <- pr[unlist(selectedInfo[lengths(selectedInfo) > 0]),]
      selectedInfo <- pr[selectedInfo,]
      # selectedInfotemp <<- selectedInfo
      
    }else{
      pre <- precursors()
      mzi <- as.character(input$mz)
      selectedInfo <- pre[pre$mz == paste("mz",mzi,sep = "") &pre$Sequence==se,]
      
    }
    if(dim(selectedInfo)[1]==0){
      showNotification(session,paste("No match found:",paste("mz = ",mzi),";",paste("Sequence =",se)))
      
    }
    selectedInfo <- unique(selectedInfo)
    return(as.character((selectedInfo$Species)))
    
  })
  ## Precursor selction Ui
  output$precursors <- renderUI({
    Precursors <- grep("^mz",dblistT(db = dbpath()),invert = F,value = T,fixed = F)
    
    
    
    selectInput("precursors.i", "Select Precursors", 
                choices  = c(setdiff(Precursors,c("msmsScans","msmsScansRev"))))
    
  })
  
  # TransitionsSelection ########
  ## Control Variable for transition selection
  LoadValues <- reactiveValues(x = NULL)
  ## selectizeInput vor transitions
  output$selected_transitions = renderUI({
    print("Waiting for Transitions")
    validate(need(length(colnames(TransitionsSelected()))>0,""),
             need(length(colnames(TransitionsSelected()))>0,"")
    )
    
    sel <- colnames(TransitionsSelected())
    if(length(LoadValues$x) > 0){
      Hm <- inputList()
      sel <- Hm$selected_transitions
      LoadValues$x <- NULL
    }
    
    
    
    
    trans <- colnames(Transitions())
    if(input$nomodifications){
      sel <- grep(".",sel,invert = T,value = T,fixed = T)
      sel <- grep("-",sel,invert = T,value = T,fixed = T)
      
    }
    if(input$top5){
      if(length(sel) > 10){
        sel <- sel[1:10]
      }
    }
  
    trans <- trans[trans!="rawfile"]
    print("Transitions Done")
    selectizeInput("selected_transitions", "fragments",
                   trans,selected=sel,multiple = T)
  })
  
  # Other Variables
  ## FDR
  output$FDRpep <- renderUI({
    
    # setwd(input$mainPath)
    numericInput("FDRpep", "FDR", 
                 min = 0,max = 0.5 ,value = 0.01,step = 0.001)
  })
  ## Rawfiles
  rawfiles <- reactive({
    print("Rawfiles")
    PrecursorSelectionTable <<- dbread("PrecursorSelectionTable",dbpath())
    RF <- unique(PrecursorSelectionTable$rf)
    # try({
    #   anatemp <- AnalyzedTransitions()
    #   RF <- lapply(anatemp,function(x){unique(x$rawfile)})
    #   
    #   # RF <- lapply(AnalyzedTransitions(),function(x){unique(x$rawfile)})
    #   RF <- unique(unlist(RF))
    #   RFiles <<- RF
    #   RF <- RF[RF!= ""]
    # })
    
    return(RF)
  })
  ## Rawfiles UI
  output$rawfile <- renderUI({
    # SEL <- unique(unlist(sapply(precursors.i.selected(),function(x){which(x == PrecursorSelectionTable$Species)})))
    # rf <- PrecursorSelectionTable[SEL,]
    # RF <- rf$rf[rf$Count > minthresh]
    # RF <- RF[RF != ""]
    # 
    selectedRF <- input$rawfile.i
    # # RF <<- RF
    # selectInput("rawfile.i", "Select Raw-File", 
    #             choices  = sort(RF),
    #             selected = selectedRF)
    selectInput("rawfile.i", "Raw-File",
                choices  = sort(rawfiles()),
                selected = selectedRF)
    
  })
  
  # Plotting Stuff ############
  # Transitions Plotting Function
  TransPlotReactive <- reactive({
    print("Traceplot start")
    #Obtain Data
    validate(
      need(is.list(AnalyzedTransitions()),"No Fragments"),
      need(is.list(PeaksFun()),""),
      need(file.exists(dbpath()),""),
      need(is.character(input$rawfile.i),"")
      
      
      
    )
    
    {
      print("Transplot Reactive")
      validate(need(file.exists(input$mainPath),"No valid directory."))
      setwd(input$mainPath)
      test <- T
      if(test){
        
        ana <<- AnalyzedTransitions()
        cat("\r Loading IL")
        # RA <<- intialRanges()
        IL <- IL()#read.csv(st<<<-paste(input$inclusionList,"SpectraTable.txt",sep = "/"),sep = "\t")#
        ILs <<- IL
        ilt <<- ILthermo()# Is this really necessary? Would lead to errors, if no QExactive is used.
        rfall.i <<- input$rawfile.i#rawfiles()
        rfall <<- rawfiles()
        ms1List <<- ms1Scan()
        PeaksFunFun <<- PeaksFun()
        # xlPreset <<- xlrange()
        TransCole <<- TransCol()
        precursors.i.selectedvec <<- ""#precursors.i.selected()
        dbpathvec <<- dbpath()
        SimultaneousMassShiftvec <<- input$SimultaneousMassShift
        rfi <<- input$rawfile.i
        setallrf <<- input$allRawFiles
        # TRANS <<- input$selected_transitions
        ppm1 <<- input$ppm1
        MStype <<- input$MStype
        # ylim <<- input$ylim
        # xlim <<- input$xlim
        # showNotification(paste(xl,collapse = " "))
        PeaksServer <<- input$PeaksServer
        CenterPeak <<- input$CenterPeak
        Align <<- input$Align
        RetentionTimeWindow <<- input$RetentionTimeWindow
        plottype <<- input$plottype
        FDRpep <<- input$FDRpep
        p.value <<- input$p.value
        secPlotType<<- input$secPlotType
        SignificantOnly <<- input$SignificantOnly
        # rangessetx <<- ranges$x
        # rangessety <<- ranges$y
      }else{
        ana <- AnalyzedTransitions()
        cat("\r Loading IL")
        IL <- IL()#read.csv(st<<-paste(input$inclusionList,"SpectraTable.txt",sep = "/"),sep = "\t")#
        ILs <- IL
        ilt <- ILthermo()# Is this really necessary? Would lead to errors, if no QExactive is used.
        rfall.i <- input$rawfile.i#rawfiles()
        rfall <- rawfiles()
        ms1List <- ms1Scan()
        PeaksFunFun <- PeaksFun()
        # xlPreset <- xlrange()
        TransCole <- TransCol()
        precursors.i.selectedvec <- ""#precursors.i.selected()
        dbpathvec <- dbpath()
        SimultaneousMassShiftvec <- input$SimultaneousMassShift
        rfi <- input$rawfile.i
        setallrf <- input$allRawFiles
        # TRANS <- input$selected_transitions
        ppm1 <- input$ppm1
        MStype <- input$MStype
        # ylim <- input$ylim
        # xlim <- input$xlim
        # showNotification(paste(xl,collapse = " "))
        PeaksServer <- input$PeaksServer
        CenterPeak <- input$CenterPeak
        Align <- input$Align
        RetentionTimeWindow <- input$RetentionTimeWindow
        plottype <- input$plottype
        FDRpep <- input$FDRpep
        p.value <- input$p.value
        secPlotType<- input$secPlotType
        SignificantOnly <- input$SignificantOnly
        rangessetx <- ranges$x
        rangessety <- ranges$y
      }
      if(input$allRawFiles){
        rfall.i <-  rfall
      }
      
      # 
      xl <- range(unlist(lapply(ana,function(x){range(x$RT_Used,na.rm = T)})),na.rm = T)
      yl <- range(unlist(lapply(ana,function(x){
        temp <- x[,1:(which(colnames(x)=="charge")-1)]
        temp <- max(as.numeric(unlist(temp)),na.rm = T)
        yl <- c(0,temp)
      })))
      
      
      anatemp <- ana
      
      # print("Updated ana")
      if(length(ana) == 0){
        return(NULL)
      }
      print("Addidtional settings")
      SpecIDS <- sapply(strsplit(names(ana),"_"),function(x){gsub(".txt$","",x[length(x)])})
      mz      <- gsub("^mz","",sapply(strsplit(names(ana),"_"),function(x){x[1]}))
      
      
      
      
      cat("\r Loading IL2")
      
      ilt2 <- ILs
      ilt2$Mass..m.z. <- ILs$mz
      ilt2$Species <- ILs$Charge
      ilt2$CS..z. <- ILs$Charge
      ilt2$Polarity <- "Positive"
      ilt2$Start..min. <- NA
      ilt2$End.min. <- NA
      ilt <- ilt2
      
      cat("\r Selecting Mass")
      if(any(mz==ILs$mz)){
        SelectedMass <- ILs[match(mz,ILs$mz),]
      }else{
        ppmwindow <- as.numeric(mz)[1]/(1/ppm1*1000000)
        ppmwindow <- c(as.numeric(mz)[1]-ppmwindow/2,as.numeric(mz)[1]+ppmwindow/2)
        wh <- which(ILs$mz >=min(ppmwindow) & ILs$mz <= max(ppmwindow))
        
        
        SelectedMass <- ILs[wh,]
        if(length(wh)==0){
          
          print("Problem in matching SelectedMass with SPectratable")
        }
      }
      # print(SelectedMass)
      
      TransitionColors <- T
      # if(length(ana) > 1){
      #   TransitionColors <- F
      # }
      # if(SimultaneousMassShiftvec ){
      #   par(mai = c(0.1,0.8,0.1,0.9),mfrow = c(length(ana),1))
      #   
      # }else{
      #   par(mai = c(0.8,0.8,0.1,0.9),mfrow = c(length(ana),1))
      #   
      # }
      # if(setallrf){
      #   rfall.i.initialized <- rfall.i
      #   rfall.i <- rfall
      #   requiredPlots <- length(rfall.i) * length(ana)
      #   
      #   
      #   
      #   par(mai= c(0.1,0.1,0.1,0),mfcol = c(ceiling(requiredPlots^0.5),ceiling(requiredPlots^0.5)))
      #   CenterPeakAlternative <- T
      #   
      # }else{
      #   CenterPeakAlternative <-F
      # }
      cat("\r started loop")
      TransplotListALL <- list()
      rf <- rfall.i
      ITplotting <- 0
      if(length(session) > 0){
        progress <- Progress$new(session,min=0, max=length(rfall.i))
        on.exit(progress$close())
        progress$set(message = 'Starting Export',
                     detail = "",value = 0)
      }
      
      it <- 0
      for(rf in rfall.i){
        it <- it+1
        if(MStype == "MS1"){
          print("Preparing MS1 scan")
          for(ms1 in ms1List){
            ms1scan_table <- data.table(ms1)
            if(length(ms1scan_table)==0){
              plot(1,main = "no data",axes = F,xlab = "",ylab="",type = "n")
              return(NULL)
            }else{
              
              
              #rfall <- rf
              ms1scan_table <- ms1scan_table[rawfile==rf,]
              # yl <- ranges$y
              # yl[1] <-0
              # 
              # xl <- ranges$x
              # ms1scan_table <<- ms1scan_table
              # xl <- xl
              # yl <- yl
              # if(length(xl) != 2){
              # }
              # if(length(yl)!=2){
              # 
              # }
              
              maxval <- unique(unlist(ms1scan_table[,max(.SD,na.rm = T),.SDcols=grep("mz.",colnames(ms1scan_table))]))
              yl <- c(0,maxval)
              ITplotting <-ITplotting+1
              xl <- range(ms1scan_table$RT,na.rm = T)
              
              
              tempList <-list(ms1scan_table=ms1scan_table,xl=xl,yl=yl)
              class(tempList) <- c("ms1scanlist","list")
              TransplotListALL[[ITplotting]] <- tempList
              
              
              # detecting MS2 ANalysis range:
              # try({
              #   anafun <- sapply(ana,function(x){c(x$mz[1],range(x$RT_Used,na.rm = T))})
              #   mzms1 <- ms1scan_table$mz[1]
              #   limits <- anafun[round(anafun[1,],1)==round(mzms1,1),]
              #   
              #   if(length(limits)>0){
              #     RA <- range(anatemp$RT_Used[anatemp$rawfile==rf],na.rm = T)
              #     abline(v=RA,col = "blue",lty = "dotted")
              #     mtext( paste(unique(ms1scan_table$mz),unique(rf),sep= "\n"),cex = 0.7,adj = 0,col = "blue",line = -1,xpd = NA)
              #   }
              #   
              # })
              
              
            }
          }
          
          
          
          
          
        }
        if(MStype == "MS2"){
          
          # print(rf)
          for(i in 1:length(ana)){
            ITplotting <-ITplotting+1
            print(paste("ITplotting",ITplotting))
            
            ILtemp  <- SelectedMass[i,]
            
            add = F
            # rawfile.i
            subDat <- ana[[i]]
            subDat <- subDat[subDat$rawfile == rf,]
            
            if(dim(subDat)[1]==0){
              # print("NO DATA")
              return("NO DATA")
            }
            # yl <- ylim
            # xl <- xlim
            # showNotification(paste(xl,collapse = " "))
            Pe <- PeaksServer
            Ce <- CenterPeak
            
            # yl <-
            PeaksSel <- NA
            
            # Defining Peaks:
            PeaksTable <- PeaksFunFun[[i]]
            Peaks <- PeaksTable$Peaks
            names(Peaks) <- PeaksTable$rawfile
            PeaksSel <- Peaks[names(Peaks) == rf]
            # print("subset2")
            if(!Align){
              Ppos <- subset(PeaksTable,select = c("Q1","Q2"))
              
            }else{
              if(any(colnames(PeaksTable) == "Q1align")){
                Ppos <- subset(PeaksTable,select = c("Q1align","Q2align"))
              }else{
                print("WARNING no alignment information available")
                Ppos <- subset(PeaksTable,select = c("Q1","Q2"))
                
              }
              
            }

            Pposi <- Ppos
            Ppos <- Ppos[PeaksTable$rawfile == rf,]
            Ppos <- unlist(Ppos)
            Ppos <- range(unlist(Ppos),na.rm = T)
            
            
            # Defining Ranges: XL
            
            # if(length(rangessetx) > 0){
            #   xl <- rangessetx
            # }
            # if(length(xl)==0&SimultaneousMassShiftvec){
            #   xl <- xlPreset
            # }
            
            # XL Center Peak:
            print("Center Peak")
            PeaksSel <- PeaksSel
            
            if(Ce){
              print(Ppos)
              if(length(Ppos)==2&all(!is.infinite(Ppos))&all(!is.na(Ppos))){
                xl <- Ppos+diff(Ppos)*c(-0.1,0.1)
              }
            }
            
            
            cat("\rStarting TracePlot")
            colAll = 2
            #LoadSettings(type = plottype,frame = allrf,FDRCutOff = FDRpep,Ppos = as.numeric(Ppos),colmap = TransCole,xl =xl,yl = yl,TransitionColors = TransitionColors,col = colAll,add = add,secPlotType = secPlotType,p.value = p.value,RetTime = c(ILtemp$Start,ILtemp$End),PRMonly = T,AddSecondPlot = precursors.i.selected() != "all",onlySignificant = SignificantOnly,xlab = "Retention time [min]",ylab = "Intensity",blankplot = allrf)
            # subDat <- xhui
            
            
            ylset <<- yl
            xlset <<- xl
            # if(length(xl)==2&all(!is.na(xl))){
            #   sel <-subDat$RT_Used>=min(xl)&subDat$RT_Used<=max(xl)
            #   subDat<-subDat[sel,]
            #
            # }
            cat("\rStart Traceplot")
            
            # p1<- TransitionGGplotG(subDat,ILtemp,Ppos,secPlotType)
            # p1
            
            
            cat("\rEnded TracePlot")
            
            
            
            TransplotList <- list(x=subDat,
                                  type = plottype,
                                  frame = setallrf,
                                  FDRCutOff = FDRpep,
                                  Ppos = as.numeric(Ppos),
                                  colmap = TransCole,
                                  xl =xl,
                                  yl = yl,
                                  TransitionColors = TransitionColors,
                                  col = colAll,
                                  add = add,
                                  secPlotType = secPlotType,
                                  p.value = p.value,
                                  RetTime = c(ILtemp$Start,ILtemp$End),
                                  PRMonly = T,
                                  AddSecondPlot = precursors.i.selectedvec != "all",
                                  onlySignificant = SignificantOnly,
                                  xlab = "Retention time [min]",
                                  ylab = "Intensity",
                                  blankplot = setallrf,
                                  namesAna = names(ana)[i],
                                  SelectedPrecursor = precursors.i.selectedvec[1],
                                  Peaks=Peaks,
                                  PeaksSel= PeaksSel,
                                  SelectedMass=SelectedMass,
                                  rf = rf)
            class(TransplotList) <- c("TransplotList","list")
            
            TransplotListALL[[ITplotting]] <- TransplotList
            
            
          }
        }
        
        progress$set(message = 'Preparing Data',
                     detail = "rf",value =it )
        
      }
      
      
      # try({
      #   # print("Updating Trans Selection")
      #   
      #   try(rm(transname),silent = T)
      #   try(transname <- dbtaName(ana,dbpathvec))
      #   if(exists("transname")){
      #     dbp <- dbpathvec
      #     dbTranssuccess <- dbwrite(x = data.frame(TRANS),name = gsub("^PEAKS","TRANS",transname),db =dbp,overwrite = T)
      #     
      #   }
      #   
      # })
    }
    
    
    TransplotListALL
  })
  # Transitions Plotting
  ## Precursor Colors
  TransCol <- reactive({
    return(cbind(colnames(Transitions()),rainbow(length(colnames(Transitions())))))
  })
  
  output$PlotOutput <- renderPlot({
    CheckTransPlot <<- TransPlotReactive()
    print("TransPlotReactive Result:")
    print(length(CheckTransPlot))
    validate(need(is.list(CheckTransPlot),"No Data"))
    # add validate here 
    
    # validate(
    #   need(is.character(input$sequence),"sequence is not a character")
    # )
    # 
    # validate(
    #   need(length(precursors())>0, "No sequence selected"), # display custom message in need
    #   need(length(input$sequence)>0, "No sequence selected"), # display custom message in need
    #   need(req(precursors()),"no Sequences"),
    #   need(is.list(CheckTransPlot),"")
    #   
    #   
    # )
    
    # fix xl
    maxXL <- max(unlist(lapply(CheckTransPlot,function(x){x$xl})))
    # ra <<- ranges()
    CheckTransPlot <- lapply(CheckTransPlot,function(x){
      x <- x
      # x$xl <- c(0,maxXL)
      if(length(ranges$x) ==2 ){
        x$xl <- ranges$x
      }
      if(length(ranges$y) ==2 ){
        x$yl <- ranges$y
      }
      x
    })
    
    
    # x <- CheckTransPlot[[1]]
    if(input$SimultaneousMassShift){
      par(mai = c(0.5,0.8,0.1,0.9),mfrow = c(length(CheckTransPlot),1),mgp=c(1.4,0.5,0))
    }else{
      par(mai = c(0.8,0.8,0.1,0.9),mfrow = c(length(CheckTransPlot),1))
    }
    
    if(input$allRawFiles){
      
      requiredPlots <- length(CheckTransPlot)
      
      
      
      par(mai= c(0.5,0.1,0.1,0),mfcol = c(ceiling(requiredPlots^0.5),ceiling(requiredPlots^0.5)),mgp=c(0,0.2,0.5))
      CenterPeakAlternative <- F
      
    }else{
      CenterPeakAlternative <-F
    }
    selected_transitions <<- input$selected_transitions
    lapply(CheckTransPlot,plot,maximalnumber = 500)
    if(length(CheckTransPlot)>0&!input$allRawFiles){
      par(new=T,mfrow = c(1,1))
      plot(1, type  = "n",xlab = "",ylab = "",main = "",frames = "",xlim =CheckTransPlot[[1]]$xl,frame = F,axes = F)
    }
    
    
    
    
  })
  # Peptide
  output$PeptidePlot <- renderPlot({
    # print("Starting Protein Plot")
    PeaksTable <- PeaksFun()
    hum <- ""#precursors.i.selected()
    if(!is.data.frame(PeaksTable)){
      PeaksTable <- PeaksTable[[which(names(PeaksTable) == hum)]]
      
    }
    slt <- input$selected_transitions
    Ma <- match(colnames(PeaksTable),slt)
    Ma <- !is.na(Ma)
    Ma[which(colnames(PeaksTable) == "Q1"):dim(PeaksTable)[2]] <- T
    PeaksTable <- PeaksTable[,Ma]
    # print("Starting Protein Quan")
    hu <- NULL
    if(sum(!is.na(PeaksTable$Q1))<input$ValuesAcrossSampleThreshold){
      plot(1,xlab = "",ylab = "",main = "Nope, no plot for you: Not enough data...",type = "n",frame =F,axes = F)
      
    }else{
      
      
      pl  <- ProteinQuan(PeaksTable,ValuesAcrossSampleThreshold = input$ValuesAcrossSampleThreshold,name= input$rawfile.i,log10 = input$log10)
      # print("Starting Protein Quan Plot")
      layout(matrix(c(2,1),ncol = 1,nrow = 2),heights = c(0.2,1))
      
      RT <- pl$RT
      usr <- par()$mai
      par(mai = c(2,usr[2],0,usr[4]))
      hu<-plotPeptideList(list(pl),colmap = TransCol(),alpha = input$FDRpep)
      usrfun <- par()$usr[1:2]
      par(mai = c(0,usr[2],0,usr[4]))
      try({
        
        
        plot(hu[,1],RT,frame =F,pch = 3,axes = F,xlim = usrfun,type = "n")
        abline(h=pretty(RT),col = "grey",lty = "dotted")
        QU <- quantile(RT,probs = c(0.10,0.9),na.rm = T)
        try(polygon(c(0,length(RT)+2,length(RT)+2,0),c(QU[1],QU[1],QU[2],QU[2]),col = "lightgrey",border = NA))
        points(hu[,1],RT,pch = 3)
        
        axis(2,las = 2)
      })
      
      
    }
    return(hu)
    
    
  })
  # total scanrange of current run
  intialRanges <- reactive({
    print("Initialranges")
    xr <- range(unlist(lapply(AnalyzedTransitions(),function(x){range(x$RT_Used,na.rm = T)})),na.rm = T)
    
    if(input$MStype == "MS2"){
      yr <- range(unlist(lapply(AnalyzedTransitions(),
                                function(x){x <- x[,1:(which(colnames(x)=="charge")-1)];range(as.numeric(unlist(x)),na.rm = T)})),na.rm = T)
      yr[1] <- 0
    }
    if(input$MStype == "MS1"){
      ms1datL <- ms1Scan()
      yr <- lapply(ms1datL,function(ms1dat){
        yr <- range(unlist(ms1dat[,{.SD},.SDcols=grep("^mz.",colnames(ms1dat))]),na.rm = T)
        return(yr)
      })
      yr <- range(unlist(yr),na.rm = T)
      
      yr[1] <- 0
    }
    # yr <- range(unlist(lapply(AnalyzedTransitions(),
    #                           function(x){x <- x[,1:(which(colnames(x)=="charge")-1)];range(as.numeric(unlist(x)),na.rm = T)})),na.rm = T)
    # yr[1] <- 0
    return(list(x=xr,y = yr))
  })
  # Ranges from manual selection
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  
  ## navigation arrows ----------
  # Zoom
  observeEvent(input$Plus,{
    
    if(length(ranges$x) == 0){
      ranges$x <-  intialRanges()$x
    }
    d <- diff(ranges$x)
    ranges$x <- ranges$x+c(+d/20,-d/20)
    
  })
  # Unzoom
  observeEvent(input$Minus,{
    if(length(ranges$x) == 0){
      ranges$x <-  intialRanges()$x
    }
    d <- diff(ranges$x)
    ranges$x <- ranges$x+c(-d/5,+d/5)
    
  })
  # Move to the right
  observeEvent(input$right,{
    if(length(ranges$x) == 0){
      ranges$x <-  intialRanges()$x
    }
    
    d <- diff(ranges$x)
    ranges$x <- ranges$x+c(+d/20,+d/20)
    
  })
  # Move to the left
  
  observeEvent(input$left,{
    if(length(ranges$x) == 0){
      ranges$x <-  intialRanges()$x
    }
    d <- diff(ranges$x)
    ranges$x <- ranges$x+c(-d/20,-d/20)
    
  })
  # Full height
  observeEvent(input$vertical,{
    
    ranges$y <-  intialRanges()$y
    
  })
  # Full width
  observeEvent(input$horiz,{
    ranges$x <- intialRanges()$x
  })
  # Increase height
  observeEvent(input$height,{
    if(length(ranges$y) == 0){
      ranges$y <-  intialRanges()$y
    }
    d <- diff(ranges$y)
    
    ranges$y <- ranges$y+c(0,-d/20)
    
  })
  ## interactive Zoom with mouse------------
  observeEvent(input$plot1_dblclick, {
    # print("START reavtiveValue MOD")
    
    brush <- input$plot1_brush
    brvec <- brush
    # ravec <<- ranges$x
    if(length(ranges$x)==0){
      ravecs <- intialRanges()
    }else{
      ravecs<- ranges
    }
    
    ConvertBrushToRT <- function(ravec,brush,type = "x"){
      ravec <<- ravec
      brush <<- brush
      
      ravecdiff <- diff(ravec)
      RAVEC <- NULL
      if(type=="x"){
        RAVEC <- min(ravec)+ravecdiff*c(brush$xmin, brush$xmax)
      }
      if(type=="y"){
        RAVEC <- min(ravec)+ravecdiff*c(brush$ymin, brush$ymax)      
      }
      print(RAVEC)
      if(type!="x"&type!="y"){
        warning("Type does not match x or y in ConvertBrushToRT. No Change applied")
      }
      return(RAVEC)
    }
    
    
    if(use_ggplot){
      if (!is.null(brush)) {
        if(input$SimultaneousMassShift){
          # usr<- par()
          # ana <- AnalyzedTransitions()
          # Diff <- length(ana)
          # bre <- lapply(ana,function(x){
          #   x <- x
          #   yr <- max(x[,1:(which(colnames(x) == "charge")-1)],na.rm = T)
          # })
          # bre <- max(unlist(bre))
          # # ym <- brush$ymin +brush$ymax/Diff
          # # brmax <<- brush$ymax
          # # brmin <<- brush$ymin
          # # ranges$x <- c(brush$xmin, brush$xmax)
          # # ranges$y <- c(ym, brush$ymax)
          
          ranges$x <- ConvertBrushToRT(ravecs$x,brush)
          # ranges$y <- c(brush$ymin, brush$ymax)
        }else{
          ranges$x <- ConvertBrushToRT(ravecs$x,brush)
          ranges$y <- ConvertBrushToRT(ravecs$y,brush,"y")
        }
        
        
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    }else{
      if (!is.null(brush)) {
        if(input$SimultaneousMassShift){
          # usr<- par()
          # ana <- AnalyzedTransitions()
          # Diff <- length(ana)
          # bre <- lapply(ana,function(x){
          #   x <- x
          #   yr <- max(x[,1:(which(colnames(x) == "charge")-1)],na.rm = T)
          # })
          # bre <- max(unlist(bre))
          # # ym <- brush$ymin +brush$ymax/Diff
          # # brmax <<- brush$ymax
          # # brmin <<- brush$ymin
          # # ranges$x <- c(brush$xmin, brush$xmax)
          # # ranges$y <- c(ym, brush$ymax)
          ranges$x <- c(brush$xmin, brush$xmax)
          # ranges$y <- c(brush$ymin, brush$ymax)
        }else{
          ranges$x <- c(brush$xmin, brush$xmax)
          ranges$y <- c(brush$ymin, brush$ymax)
        }
        
        
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    }
    
  })
  #### Saving Sessions ---------
  observeEvent(input$saveSessionPaths,{
    fo <- paste(SystemPath,"Sessions",sep = "/")
    if(!file.exists(fo)){
      dir.create(fo)
    }
    Fi <- tempfile(tmpdir = fo)
    Fi <- paste(Fi,basename(input$mainPath),basename(input$inclusionList),input$ppm1,input$ppm2,sep = "#")
    Li <- list(mainPath = input$mainPath,MaxQuant = input$MaxQuant,inclusionList=input$inclusionList,ppm = input$ppm1,ppm2 = input$ppm2,threads = input$Threads)
    save(Li,file = Fi)
    session$reload()
    
  })
  # Go of Rawfile Scanner Tool
  observeEvent(input$goButton ,{
    fo <- paste(SystemPath,"Sessions",sep = "/")
    if(!file.exists(fo)){
      dir.create(fo)
    }
    Fi <- tempfile(tmpdir = fo)
    Fi <- paste(Fi,basename(input$mainPath),basename(input$inclusionList),input$ppm1,input$ppm2,sep = "#")
    
    Li <- list(mainPath = input$mainPath,MaxQuant = input$MaxQuant,inclusionList=input$inclusionList,ppm = input$ppm1,ppm2 = input$ppm2,threads = input$Threads)
    save(Li,file = Fi)
    
    try({
      
      Paths <- sapply(Paths(),function(x){x$mainPath})
      Paths <- unique(Paths)
      Paths <- Paths[file.exists(Paths)]
      
      updateSelectInput(inputId = "Sessions",choices = Paths,selected = Fi$mainPath)
    })
    
    
    wd <- getwd()
    try(PrepareTransitionList(input$mainPath,input$MaxQuant,input$inclusionList,ppm = input$ppm1,ppm2 = input$ppm2,
                              session = session,threads = input$Threads,
                              #pythonpath =pythonpath[as.numeric(input$Python3)+1],
                              test = input$Recompile,
                              useDIA=input$DIA))
    
    session$reload()
    validate(need(file.exists(wd),"No valid directory."))
    
    setwd(wd)
    # CheckTransPlot()
  })
  
  
  # PREPARING DATA #############
  
  # Prepare Transitionlist  (checking Path)
  TransitionListsDir <- reactive({
    Path <- "Transition matrix is not compiled."
    validate(need(file.exists(input$mainPath),"No valid directory."))
    
    setwd(input$mainPath)
    Path <- "./PRM_Analyzer_Matches/"
    return(Path)
  })
  ## show directory Path
  output$analyzedListPath <- renderText({TransitionListsDir()})
  # Read Inclusion List
  IL <- reactive({
    IL <- read.csv(paste(input$inclusionList,"SpectraTable.txt",sep = "/"),sep = "\t")
    IL <- IL
    return(IL)
  })
  # Read InclusionsList Thermo (should become obsolete)
  ILthermo <- reactive({
    P <- paste(input$inclusionList,"InclusionList.txt",sep = "/")
    if(!file.exists(P)){
      P <- paste(input$inclusionList,"InclusionList.csv",sep = "/")
    }
    if(!file.exists(P)){
      showNotification("Couldn't find InclusionList.txt")
      P <- list.files(input$inclusionList,pattern="InclusionList",full.names = T)[1]
      showNotification(paste("Couldn't find InclusionList.txt. Trying",P))
      
      
    }
    IL <- read.csv(P,sep = "\t")
    # ILthermo <<- ILthermo
    return(IL)
  })
  # Noramlization based on msms.txt (requires MQ)
  MSMS <- reactive({
    print("MSMS Stuff?")
    # Creates RawFile Specific Normalization Factors
    evipath <- paste(input$MaxQuant,"evidence.txt",sep = "/")
    if(file.exists(evipath)){
      msms <- fread(evipath,sep = "\t",stringsAsFactors = F)
      FunProduct <- function(y,...){median(as.numeric(y),...)}
      FunPrecursor <- function(y,...){sum(as.numeric(y),...)}
      msmsReferenz <- aggregate(msms$Intensity,list(msms$`Raw file`),function(x){
        tempProd <- sapply(strsplit(as.character(x),";"),FunProduct,na.rm = T)
        tempProd <- FunPrecursor(tempProd,na.rm = T)
      })
      msmsReferenz[,2] <- mean(msmsReferenz[,2] )/msmsReferenz[,2]
      return(msmsReferenz)
    }else{
      
      ana <- AnalyzedTransitions()
      
      SUMI <- lapply(ana,function(x){
        
        x <- data.table(x)
        x[,.(x={
          temp <- .SD[,1:(which(colnames(.SD)=="charge")-1)]
          temp <- unlist(temp)
          temp <- temp[temp!= -1]
          median(unlist(temp),na.rm = T) 
        }),rawfile]
      })
      SUMI <- SUMI[[1]]
      SUMI$x <- mean(SUMI$x,na.rm = T)/SUMI$x
      SUMI <- data.frame(SUMI)
      SUMI$Group.1 <- SUMI$rawfile
      return(SUMI)
      
    }
    
    
    
  })
  # Reads TableList  from database-----
  analyzedList <- reactive({
    print("Reading analyzedList")
    validate(need(file.exists(dbpath()),"No Database available."),
             need(file.exists(input$mainPath),"Please set an valid analysing folder path."))
    
    setwd(input$mainPath)
    dbp <- dbpath()
    analyzedList <-  grep("^mz",dblistT(db = dbp),invert = F,value = T,fixed = F)
    
    
    return(analyzedList)    
    
    
  })
  
  ## MS1
  # Reads MS1 Scan Table 
  ms1Scan <- reactive({
    # added 20190513
    validate(
      need(is.character(input$sequence),"sequence is not a character")
    )
    precur <<- precursors.i()
    
    ms1scanList <- lapply(precur,function(ms1_mz){
      # se <<- sapply(strsplit(ms1_se," "),function(x){x[1]})
      ms1_mz <<- ms1_mz
      db <- dbConnect(SQLite(),dbpath())
      dbl <- dbListTables(db)
      
      ms1_mz <- paste(unlist(strsplit(ms1_mz,"_"))[1:3],collapse = "_")
      
      ms1scan_Name<- grep(gsub("^mz","ms1scans",ms1_mz),dbl,value = T)[1]# 
      
      # ms1scan_Name <- ms1scan_Name
      # validate(need(!is.na(ms1scan_Name),"No MS1 tables available."))
      
      if(is.na(ms1scan_Name)){
        print(paste("No MS1 Information Available.",ms1_mz,sep= " "))
        return(NULL)
      }
      
      ms1scan <- dbReadTable(db,ms1scan_Name)
      dbDisconnect(db)
      ms1scan <- data.table(ms1scan)
    })
    
    return(ms1scanList)
  })
  
  ## MS2
  # Reads Transitions 
  AnalyzedTransitions1 <- reactive({
    print("Starting AnalyzedTransitions1")
    validate(need(is.character(precursors.i()),""),
             need(file.exists(dbpath()),""),
             need(is.character(input$rawfile.i),""))
    # print("Pausing AnalyzedTransitions1")
    # Sys.sleep(2)
    if(precursors.i() == "all"){
      # not available anymore
      tl <- analyzedList()
    }else{
      Pre <- precursors.i()
      aL <- analyzedList()
      tln <- aL[match(Pre,aL)]
      if(all(is.na(tln))){
        showNotification("No entry match in the database.")
        return(NULL)
      }
      
      tl <- lapply(tln,function(tlnx){
        
        Tempx <- dbread(tlnx,db = dbpath())
        
        if(length(Tempx$RT2) ==0){
          Tempx$RT2 <- as.numeric(Tempx$RT_min)
        }
        Tempx$Peaks <- NA
        Tempx$Peaks1 <- NA
        Tempx$Peaks2 <- NA
        if(length(Tempx$RF_Scores)==0){
          Tempx$RF_Scores <- 0
        }
        if(length(Tempx$DL_Scores)==0){
          Tempx$DL_Scores <- -5
        }
        
        if(input$Align){
          Tempx$RT_Used <- as.numeric(Tempx$RT2)
          Tempx$Peaks  <- as.numeric(Tempx$Peaks2) 
        }else{
          Tempx$RT_Used <- as.numeric(Tempx$RT_min)
          Tempx$Peaks  <- as.numeric(Tempx$Peaks1)
        }
        fifu <- paste(precursors.i(),input$rawfile.i,sep = "_")
        
        if(!file.exists(PATHI <-paste(input$mainPath,"settings",gsub(".txt$","",fifu),sep = "/"))){
          tr <- colnames(Tempx)[1:(which("charge" == colnames(Tempx))-1)]
        }
        
        if(!file.exists(PATHI <-paste(input$mainPath,"settings",gsub(".txt$","",fifu),sep = "/"))){
          tr <-  colnames(Tempx)[1:which("charge" == colnames(Tempx))-1]
          selec <- tr
          try({
            try(transname <- gsub("^mz","TRANS",tlnx))
            dbl <- dblistT(dbpath())
            if(any(dbl==transname)){
              dbp <- dbpath()
              TransSaved <- NULL
              try(selec <- unlist(dbread(gsub("^PEAKS","TRANS",transname),db =dbp)))
            }
          })
          print("Update SelectizeInput")
          # tri <<- tr
          selec.i <- selec
          selec.i <- selec.i[selec.i!="rawfile"]
          tr <- tr[tr!="rawfile"]
          updateSelectizeInput(session,"selected_transitions",choices =tr,selected =selec )

        }
        return(Tempx)
      })
      
      names(tl) <- tln
      
    }
    tl <<- tl
    # decide for FDR scoring column:
    tl <- lapply(tl,function(x){
      x$FDR <- 1
      try({x$FDR <- x$FDR_DL_all})# Decide here for Selected FDR Scoring model
      x[x$rawfile != "",]
      x$FDR <- p.adjust(x$FDR,"none")
      x
    })
    print("AnalyzedTransitions1 Done")
    return(tl)
  })
  # Modify Transitions# should become obsolete
  AnalyzedTransitions <- reactive({
    A1 <- AnalyzedTransitions1()
    validate(
      need(is.list(A1),"No Fragments available")
    )
    if(length(AnalyzedTransitions1()) == 0){return(NULL)}
    if(length(AnalyzedTransitions1()) > 0){
      print("Starting AnalyzedTransitions")
      
      A  <<- AnalyzedTransitions1()
      st <- NULL
      if(length(st) > 0){
        A <- lapply(A,function(x){
          tr <- colnames(x)[1:(which("charge" == colnames(x))-1)]
          rm <- which(is.na(match(tr,st)))
          if(length(rm) > 0){
            x <- x[,-rm]
          }
          return(x)
        })
      }
      
      # Apply ZeroNeighbour Imputation
      if(input$EnableNeighborZeroImputation){
        progress <- Progress$new(session,max=length(A))
        progress$set(message= "Applying ZeroNeighbour Imputation",value=0)
        on.exit(progress$close())
        itx <<- 0
        A <- lapply(A,function(tempi){
          itx <<- 1
          
          tempi <<- tempi
          # stop()
          
          tra <- tempi[,1:(which(colnames(tempi)=="charge")-1)]
          cat("\rFinding Zeros")
          tra2 <- apply(tra,2,function(x){
            x <<- x
            xfu <- which(x<=0)
            whichxfu <- which(diff(xfu)>1)
            if(length(whichxfu)>0){
              zero <- sapply(whichxfu,function(i){
                i <- i
                # cat("\r",i)
                be <- F
                af <- F
                
                
                if(x[i]<=0){
                  try(be <- x[i-1]>1)
                  try(af <- x[i+1]>1)
                  if(length(be)==0|is.na(af)){
                    be <- F
                    af <- F
                  }
                  
                  
                  if(be&af){
                    T
                  }else{
                    F
                  }
                  
                }else{
                  F
                }
                
              })
              
              x[whichxfu][zero] <- NA
            }
            
            x
            
          })
          cat("\rPerforming Imputation")
          tra3 <- Neighbour_Mean_imputation(tra2)
          tempi[,1:(which(colnames(tempi)=="charge")-1)] <- tra3
          progress$set(value=)
          
          tempi
        })
        
        # A <- lapply(A,function(tempi){
        #   tempi <<- tempi
        #   # stop()
        #   tra <- tempi[,1:(which(colnames(tempi)=="charge")-1)]
        #   tra2 <- apply(tra,2,function(x){
        #     x <- x
        #     
        #     zero <- sapply(1:length(x),function(i){
        #       i <- i
        #       # cat("\r",i)
        #       be <- F
        #       af <- F
        #       if(x[i]<=0){
        #         try(be <- x[i-1]>1)
        #         try(af <- x[i+1]>1)
        #         if(length(be)==0|is.na(af)){
        #           be <- F
        #           af <- F
        #         }
        #         
        #         
        #         if(be&af){
        #           T
        #         }else{
        #           F
        #         }
        #         
        #       }else{
        #         F
        #       }
        #       
        #     })
        #     # plot(x,type="b")
        #     # x[zero] <- NA
        #     # points(x,col=2)
        #     x[zero] <- NA
        #     
        #     x
        #   })
        #   
        #   tra3 <- Neighbour_Mean_imputation(tra2)
        #   tempi[,1:(which(colnames(tempi)=="charge")-1)] <- tra3
        #   tempi
        # })
        
      }
      if(input$supersmooth_I_set){
        it <<- 0
        A <- lapply(A,function(tempi){
          it <<- it+1
          tempi <<- tempi
          co <- colnames(tempi)
          tempi <- data.table(tempi)
          tempi2 <- tempi[,{
            # te <<- .SD
            tra <- .SD[,.SD,.SDcols=1:(which(colnames(tempi)=="charge")-1)]
            tri <- .SD[,.SD,.SDcols=(which(colnames(tempi)=="charge"):dim(.SD)[2])]
            
            tra <- apply(tra,2,function(int){
              int <- int
              rtf <- .SD$RT_Used
              sp <- input$supersmooth_bw_set
              
              intv <- supsmu(rtf,int,span=sp)$y
              intv <- intv/max(intv,na.rm = T)*max(int,na.rm = T)
              intv
            })
            
            tra <- as.data.table(tra)
            cbind(tra,tri)
            
          },rawfile]
          
          setcolorder(tempi2,co)
          data.frame(tempi2)
          # recalculate stuff:
          recal <- F
          if(recal){
            # ta_init <- dbread(mzname <- gsub("^PEAKS","mz",dbname),dbpath())
            mzname <- names(A)[it]
            SpecIDi <- as.numeric(sapply(strsplit(dbname,"_"),"[[",6))
            
            if(!is.na(SpecIDi)){
              spectrum <- TT()[SpecID==SpecIDi]
              peptidelength <- nchar(ST()[SpecID==SpecIDi]$Sequence)
              ta <- data.table(tempi2)
              tacols <- colnames(ta)
              
              ta[,c("SCAall","SCAcut"):={
                temp <<- .SD
                transipos <- temp[,.SD,.SDcols=1:(which("charge"==tacols)-1)]
                M <- match(colnames(transipos),make.names(spectrum$Matches))
                SCAall <- apply(transipos,1,function(x){NormSpecAngle(x,spectrum$Intensities_Raw[M])})
                SCAcut <- apply(transipos,1,function(x){NormSpecAngle(x[1:6],spectrum$Intensities_Raw[M][1:6])})
                
                list(SCAall,SCAcut)
                
              },.(scan,rawfile)]
              # scoring: 
            }
            
            # if(length(ta)>0&is.data.frame(ta)){
            #   ta$
            # }
            
            mztable <- ta
            MZFUN <- MakeFeatureTable(ta,peptidelength,mzname,mdpar = Models()$DeepLearning@parameters$x)
            mztable.h2o <- as.h2o(MZFUN)
            s2 <-system.time(
              PredictionScores <- sapply(Models(),function(x){
                x <<- x
                PREDI <- as.data.frame(h2o.predict(x,mztable.h2o))
                PREDI$predict
              })
            )
            PredictionScores <- as.data.frame(PredictionScores)
            
            
            ta$RF_Scores <- PredictionScores$DRF
            ta$DL_Scores <- PredictionScores$DeepLearning
            tempi2 <<- tempi2
            ta <<- ta
            setcolorder(ta,colnames(tempi2))
            setdiff(colnames(ta),colnames(tempi2))
            tempi2 <- ta
          }
          data.frame(tempi2)
  

        })
        
      }
      return(A)
    }else{return(AnalyzedTransitions1())}
    cat("AnalyzedTransitione Done")
  })
  # Transitions Selecter -------  
  Transitions <- reactive({
    print("Transitions")
    SelectedTrans <- AnalyzedTransitions1()
    if(length(SelectedTrans)> 0){
      coln <- lapply(SelectedTrans,colnames)
      coln <- sapply(coln,function(x){any(x!=coln[[1]])})
      if(any(coln)){
        showNotification("Lists with different transition order. Library might be corrupt. Chosing first one.", duration = 1, closeButton =  TRUE, type = "warning")
        
      }
      SelectedTrans <- SelectedTrans[[length(SelectedTrans)]]
      SelectedTrans <- SelectedTrans[,1:(which("charge" == colnames(SelectedTrans))-1)]
    }
    
    return(SelectedTrans)
  })
  # Check function, updating available transition checkbox
  observeEvent(input$all_Transitions,{
    print("OBSERVING ALL TRANSITIONS")
    validate(need(length(colnames(Transitions()))>0),"No Transitions selected")
    
    updateSelectizeInput(session,"selected_transitions",selected = colnames(Transitions()))
  })
  
  # available Transitions: 
  TransitionsSelected <- reactive({
    print("Transitions Check")
    SelectedTrans <- AnalyzedTransitions()
    if(length(SelectedTrans)> 0){
      SelectedTrans <- SelectedTrans[[1]]
      SelectedTrans <- SelectedTrans[,1:(which("charge" == colnames(SelectedTrans))-1)]
    }
    return(SelectedTrans)
  })
  
  ## Automated PEAKS Finding MODUL ########
  
  # Preselecting Candidates:
  RTcandidates <- reactive({
    print("RT RTcandidates")
    RTcandidatesFun(AnalyzedTransitions(),input$FDRpep)
  })
  # MAIN Peak detection/quantification Modul 
  PeaksFun <- reactive({
    print("Detecting Peaks PeaksFun")
    # reactivity
    try(MSP <- ManualSetPeak(),silent = T)
    try(REQ <- Requantify(),silent = T)
    try(SSA <- SetSearchArea(),silent = T)
    try(RMP <- ManualRemovePeak(),silent = T)
    try(input$reset)
    # Evaluation
    validate(need(is.list(AnalyzedTransitions()),""),
             need(file.exists(dbpath()),""))
    
    ### Start
    ana <<- AnalyzedTransitions()
    print("Extracting RT with best Score")
    dbp <- dbpath()
    # Looping trough Transitiontables
    CANDIDATE_RT <<- RTcandidates()
    RTwin <<- input$RetentionTimeWindow
    qType <<- input$quantitationType
    print("QuantitationLoop")
    
    Reanalysis_Check <- T
    try({
      dbl <- dblistT(db =dbp)
      dbt <- sapply(1:length(ana),function(a){dbtaName(ana[a],dbp)})
      CheckCall <- sapply(dbt,function(x){any(x==dbl)})
      if(all(CheckCall)){
        Reanalysis_Check <- F
      }
    })
    # Reanalysis_Check <- T
    print("Reanalysis_Check:")
    # print(Reanalysis_Check)
    DPW  <- try({
      dbp <<- dbp
      qType <<- qType
      RTwin <<- RTwin
      Reanalysis_Check <<- Reanalysis_Check
      # input <- list(MinPeakWidth=10,MaxPeakWidth=20,supersmooth_bw_set=F,ApplyMaximumWidth=10)
      DPlist <- DetectPeakWrapper(ana = ana,CANDIDATE_RT = CANDIDATE_RT,dbp = dbp,
                                  RetentionTimeWindow = RTwin,QType = qType,Reanalysis = Reanalysis_Check,
                                  RT_BASED_onbestScore = T,
                                  MinPeakWidth=input$MinPeakWidth,
                                  MaxPeakWidth=input$MaxPeakWidth,
                                  supersmooth_I_set = F,#input$supersmooth_I_set,
                                  supersmooth_bw_set = input$supersmooth_bw_set,
                                  ApplyMaximumWidth = input$ApplyMaximumWidth,
                                  Requantify_Priority= input$Requantify_Priority,
                                  session=session
                                  
      )
    })
    
    if(class(DPW)=="try-error"){
      validate(need(F,"Error in DetectPeakWrapper function."))
    }
    
    
    return(DPlist)
    
  })
  
  
  ## Automated Ares PEAKS Finding MODUL ########
  RTcandidatesAreaSearch <- reactive({
    print("RT candidates AreaSearch")
    RTcandidatesFun(AnalyzedTransitions(),input$FDRpep,ranges$x)
  })
  # Set Search in this Ares
  SetSearchArea <- eventReactive(input$SetSearchArea,{
    showNotification("This function is not supported in this version of Vali.")
    validate(need(F,"This function is broken at the moment"))
    progress <- Progress$new(session)
    progress$set(message= "Rescanning Window")
    on.exit(progress$close())
    print("Detecting Peaks with SetSearchArea")
    
    try(MSP <- ManualSetPeak(),silent = T)
    try(REQ <- Requantify(),silent = T)
    try(SSA <- SetSearchArea(),silent = T)
    try(RMP <- ManualRemovePeak(),silent = T)
    try(input$reset)
    ana <- AnalyzedTransitions()
    
    print("Extracting RT with best Score")
    dbp <- dbpath()
    # Looping trough Transitiontables
    CANDIDATE_RT <- RTcandidatesAreaSearch()
    #
    print("QuantitationLoop")
    
    DPlist <- DetectPeakWrapper(ana = ana,CANDIDATE_RT = CANDIDATE_RT,dbp = dbp,
                                RetentionTimeWindow = RTwin,QType = qType,Reanalysis = Reanalysis_Check,
                                RT_BASED_onbestScore = T,
                                MinPeakWidth=input$MinPeakWidth,
                                MaxPeakWidth=input$MaxPeakWidth,
                                supersmooth_I_set = F,#input$supersmooth_I_set,
                                supersmooth_bw_set = input$supersmooth_bw_set,
                                ApplyMaximumWidth = input$ApplyMaximumWidth,
                                Requantify_Priority= input$Requantify_Priority,
                                session=session
                                
                                
    )
    
    
    
    return(DPlist)
    
    
  })
  ## Requantifiy and Manual Set Peak Modul
  # Requnatify
  Requantify <- eventReactive(c(input$Requantify,input$reset),{
    print("Requantify")
    ana <- AnalyzedTransitions()
    
    if(!input$SimultaneousMassShift){#was allshifts before
      ana <- ana[names(ana) == precursors.i.selected()]
    }
    analength <- length(ana)
    anaFun <<- ana
    dbp <- dbpath()
    # Required Input:
    # dbp # Database name
    # Ana Selected # Output
    # list with detected peptides
    DPALL <<- lapply(1:analength,function(x){
      dbt <- dbtaName(anaFun[x],dbp)
      DP <- dbread(x = tempcheck <<- dbtaName(ana[x],dbpath()),db =dbp)
    })
    
    CHECK <- sapply(DPALL,function(x){dim(x)[1]})
    validate(need(length(unique(CHECK))==1&length(CHECK)==analength,"Not All Species processed"))
    # Checking Boundaries 
    
    sapply(DPALL,function(x){x$Q1})
    sapply(DPALL,function(x){x$Q2})
    
    Q1vec <- sapply(DPALL,function(x){x$Q1[DPALL[[1]]$QuantitationType=="XIC"]})
    Q2vec <- sapply(DPALL,function(x){x$Q2[DPALL[[1]]$QuantitationType=="XIC"]})
    Q2vec[Q2vec==1] <- NA
    rf <- DPALL[[1]]$rawfile[DPALL[[1]]$QuantitationType=="XIC"]
    Q1vec <- apply(Q1vec,1,min,na.rm = T)
    Q1vec[is.infinite(Q1vec)] <- NA
    Q2vec <- apply(Q2vec,1,max,na.rm = T)
    Q2vec[is.infinite(Q2vec)] <- NA
    
    for(i in 1:length(DPALL)){
      tempdata <- DPALL[[i]] 
      NAME <- dbtaName(anaFun[i],dbp)
      tempdatadt <- data.table(tempdata)
      # Finding requantify candidates
      tempdatadt[,c("Q1","Q2"):={
        temp <<- .SD
        rafi <<- .BY$rawfile
        if(all(is.na(Q1))){
          # Qs <- 
          sel <- rf==rafi
          Q1s <- Q1vec[sel]
          Q2s <- Q2vec[sel]
          
        }else{
          Q1s <- Q1
          Q2s <- Q2
        }
        
        # print(Q1)
        list(Q1s,Q2s)
      },rawfile]
      Requantified <- tempdatadt[,{
        temp <<- .SD
        if(any(temp$Requantify=="+")&any(!is.na(temp$Q1))){
          tempa <- anaFun[[i]]
          SplitList <- SplitTransitionInfo(tempa)
          X_limit <- c(temp$Q1[1],temp$Q2[1])
          PeakDetected <- DetectPeak(pe <- min(X_limit,na.rm = T)+diff(X_limit)/2,diff(X_limit)/2,SplitList$Transitions,
                                     SplitList$Info$RT_Used,
                                     presetQuantiles = X_limit,scores = SplitList$Info$DL_Scores)#$quantile
          INFO <- SplitTransitionInfo(temp,"Q1")
          
          XIC <- c(PeakDetected$XIC,INFO$Info[INFO$Info$QuantitationType=="XIC",])
          Intensities <- c(PeakDetected$intensity,INFO$Info[INFO$Info$QuantitationType=="Intensities",])
          QuantifiedTable <- rbind(XIC,Intensities)
          # temp <- as.data.table(QuantifiedTable)
        }
        
        temp
      },rawfile]
      
      # Readjusting Order
      setcolorder(tempdatadt,colnames(tempdata))
      tempdatadt$Requantify=""
      tempdatadt$Requantify[is.na(tempdata$Q1[match(tempdatadt$rawfile,tempdata$rawfile)])&!is.na(Q1)] <- "+"
      # Requantify
      
      ###
      
      # ###
      
      # setcolorder(tempdatadt,colnames(tempdata))
      dbwrite(tempdatadt,NAME,dbp,overwrite=T)
      
    }
    
  })
  # Manual
  ManualSetPeak <- eventReactive(    input$ManualSetPeak  ,{
    print("ManualSetPeak")
    # Creating CUrrent Name
    
    X_limit <<- ranges$x# Set RangeControl
    
    
    if(length(X_limit) == 0){
      progress <- Progress$new(session)
      progress$set(message= "Please define a RT window.")
      return(NULL)
    }
    ana <- AnalyzedTransitions()
    
    if(!input$SimultaneousMassShift){#was allshifts before
      ana <- ana[names(ana) == precursors.i.selected()]
      
    }
    anaFun <<- ana
    analength <- length(ana)
    # ana <<- ana
    dbp <- dbpath()
    for(a in 1:analength){
      cat(a,"Number of Peaks set")
      # print(dbtaName(ana[a],dbpath(),input))
      dbl <- dblistT(db =dbp)
      dbt <- dbtaName(anaFun[a],dbp)
      # grep("PEAKS",dbListTables(db),value= T)s
      # any(dbListTables(db)==dbt)
      # X_limit <<- X_limit
      if(any(dbl == dbt)&length(X_limit) ==2){
        print(paste("Setting Peak",names(ana)[a]))
        
        try(rm("DP"),silent = T)
        # print(dbtaName(ana[a],dbpath(),input))
        
        DP <- dbread(x = tempcheck <<- dbtaName(ana[a],dbp),db =dbset <<- dbp)
        DPNames <- colnames(DP)
        # print(colnames(DP))
        if(input$CenterPeak){
          progress <- Progress$new(session)
          progress$set(message= "Set Peak")
        }
        rfaall <- input$rawfile.i
        # DP <<- DP
        if(input$allRawFiles){
          seldp <- T
          rfaall <- unique(DP$rawfile)
        }else{
          seldp <- DP$rawfile == rfaall
          
        }
        if(any(seldp)|1){
          # print(rfaall)
          for(rfa in rfaall){
            # Quantify Peak:
            tempa <- ana[[a]]
            # tempa <- tempa[[1]]
            
            # DP <<- DP
            dbtaNameVec <- dbtaName(ana[a],dbpath())
            # stop()
            
            tempa <- tempa[tempa$rawfile == rfa,]
            transManual <- tempa[,1:(which("charge" == colnames(tempa))-1)]
            info <- tempa[,which("charge" == colnames(tempa)):dim(tempa)[2]]
            # tempa <- tempa
            transManual <- transManual
            print("\rDETECT PEAK MANUAL Start")
            X_limit <- X_limit
            PeakDetected <- DetectPeak(pe <- min(X_limit,na.rm = T)+diff(X_limit)/2,diff(X_limit)/2,transManual,info$RT_Used,presetQuantiles = X_limit)#$quantile
            # print(PeakDetected)
            print("\rDETECT PEAK MANUAL Finish")
            
            infoquantile <- info[info$RT_min>= min(PeakDetected$quantile,na.rm = T)&info$RT_min<=max(PeakDetected$quantile,na.rm= T),]
            qe2 <-  range(infoquantile$RT2)
            infoquantile <- infoquantileParser(infoquantile = infoquantile)
            infoquantile <<- unique(infoquantile)[1,]
            # dbtaNameVec <<- dbtaNameVec
            # infoquantile <<- infoquantile
            # qe2 <<- qe2
            DP <- dbread(x = dbtaNameVec,db =dbpath())
            DPinit <- DP
            # rfa <<- rfa
            seldp <- DP$rawfile == rfa
            # if(length(seldp)==0){
            #   DP <- rbind(DP,NA)
            #   DP$rawfile[length(DP$rawfile)]<- rfa
            # }
            print("Adding Quantitation Values")
            PeakDetected <<- PeakDetected
            for(qit in c("XIC","Intensities")){
              # print(qit)
              if(qit == "XIC"){
                INT <- PeakDetected$XIC
              }
              if(qit == "Intensities"){
                INT <- PeakDetected$intensity
              }
              seli <- seldp & DP$QuantitationType == qit
              Ma <- match(colnames(DP)[1:which(colnames(DP) == "Q1")-1],names(INT))
              INT <- INT[Ma]
              # INT <<- INT
              # PeakDetected <<- PeakDetected
              # pe <<- pe
              # qit <<- qit
              # rfa <<- rfa
              infoquantile <<- infoquantile
              
              INSERTVEC <- data.frame(lapply(list(as.numeric(INT), # Quantitative values of fragments
                                                  as.numeric(PeakDetected$quantile), # Peak Quantiles
                                                  as.numeric(qe2),#Q1Allign
                                                  as.numeric(pe), # no clue
                                                  qit,# Quantitation Type
                                                  as.numeric(infoquantile[1,]),
                                                  rfa # rawfile,
              ),function(x) t(data.frame(x))),stringsAsFactors = F)
              if(length(DP$Requantify)>0){
                INSERTVEC$Requantify <- "Manual"
              }
              # seli <<- seli
              # DP <<- DP
              # DPtest <- DP
              # DP <<- DP
              # seli <<- seli
              if(any(seli)&length(INSERTVEC) == dim(DP)[2]){
                print("Successful Update of PeaksTable")
                DP[seli,] <- INSERTVEC
              }else{
                # transManual <<- transManual
                # tempa <<- tempa
                # DPProblem <<- DP
                # INSERTVEC <<- INSERTVEC
                # tm <<- transManual
                # ir <<- info$RT_min
                # X_limit <<- X_limit
                # stop()
                print("Problem in Manual Set Peaks Option")
                
              }
            }
            DPout <<- DP
            print("Writing to db")
            dbwrite(x = DP,name = dbtaNameVec,db =dbpath(),overwrite = T)
            print("Changed")
            DPinit <<- DPinit
            # print(paste(DPinit$Q1,DPinit$Q2))
            # print(paste(DPout$Q1,DPout$Q2))
            # print(dbtaNameVec)
            print("Finished writing to db")
            # DPtest <- dbread(x = dbtaNameVec,db =dbpath())
            
            # progress <- Progress$new(session)
            # on.exit(progress$close())
            # progress$set(message= "Peak Updated")
            # print(DP[sel,])
            # return(c(min(X_limit)+diff(X_limit)/2,X_limit))
          }
        }else{
          print("No Manual Set")
          
        }
        
        
      }else{
        print("NoEntry in the database");showNotification("No entry match",duration = 3)
      }
      
    }
    return(NULL)
  })
  # Manual RemovePeak
  ManualRemovePeak <- eventReactive(
    input$RemovePeak
    ,{
      print("Initiating removal of enries.")
      ana <- AnalyzedTransitions()
      # ana <<- ana[names(ana) == precursors.i.selected()]
      p <- dbpath()
      ObjectExists <- sapply(dbtaName(ana,p),function(x){
        any(dblistT(db =p) == x)
      })
      
      if(any(ObjectExists)){
        print("Removing PEak")
        
        try(rm("DP"))
        dbtaNames <- dbtaName(ana,p)
        sapply(dbtaNames,function(x){
          try({
            
            
            DP <- dbread(x = x,db =p)
            
            "Animal_20200206_HZ_SC_PRM_brUbiLP_mutated_B_Lys8_2.raw"
            
            DP[DP$rawfile == input$rawfile.i,1:(which(colnames(DP)== "Q2align"))] <- NA
            
            dbwrite(x = DP,name = x,db =dbpath(),overwrite = T)
            print("Removed Peak")
          })
          NULL
          
        })
        
        
        
        
        return(NULL)
        
        
      }else{
        print("No Table to remove.")
      }
      print("RemovePeak finished.")
      
      
    })
  
  # RAWFILEREADER Module ########
  # setting Paths:
  observeEvent(input$mainPathButton,{
    BrowseFun("mainPath",session = session)
  })
  observeEvent(input$MaxQuantButton,{
    BrowseFun("MaxQuant",session = session)
  })
  
  observeEvent(input$inclusionListButton,{
    BrowseFun("inclusionList",session = session)
    
  })
  ### Save Settings for later #-----------
  
  # input List for automated Export
  inputList <- reactive({
    validate(need(is.character(precursors.i()),"No Precursors."),
             need(is.character(input$rawfile.i),"No Raw Files selected."))
    print("List for automated Export")
    fifu <- paste(precursors.i(),input$rawfile.i,sep = "_")
    if(length(fifu) > 0){
      wd <- getwd()
      validate(need(file.exists(input$mainPath),"No valid directory."))
      
      setwd(input$mainPath)
      tempdir <- "settings"
      
      dir.create(tempdir)
      setwd("settings")
      if(file.exists(gsub(".txt$","",fifu))){
        load(gsub(".txt$","",fifu))
        inputList <- inputList
        if(inputList$PeakElution == ""){
          inputList$PeakElution <- "NA"
        }
        
        # Peak
        updateSliderInput(session,"PeakWidth",value = inputList$PeakWidth)
        updateCheckboxInput(session,"CenterPeak",value = inputList$CenterPeak)
        updateSliderInput(session,"RetentionTimeWindow",value = inputList$RetentionTimeWindow)
        updateSliderInput(session,"Match.Count",value = inputList$Match.Count)
        
        # FDR
        updateSliderInput(session,"p.value",value = inputList$p.value)
        # Transitions
        updateSelectInput(session,"selected_transitions",selected = inputList$selected_transitions)
        updateSliderInput(session,"ValuesAcrossSampleThreshold",value = inputList$ValuesAcrossSampleThreshol)
        updateCheckboxInput(session,"SignificantOnly",value = inputList$SignificantOnly)
        
        # updateSliderInput(session,"xlim",value = inputList$xlim)
        # updateSliderInput(session,"ylim",value = inputList$ylim)
        # General Settings: quantitationType, Align
        
        
      }
      setwd(wd)
      LoadValues$x <- 1
      return(inputList)
    }else{return(NULL)}
  })
  # Peak Assignments
  observeEvent(input$GoodPeak,{
    print("GoodPeak")
    dbn <- dbtaName(AnalyzedTransitions(),dbpath())
    for(dbni in dbn){
      SavePrecursorRating(rating = "good",db = dbpath(),dbtaName = dbni )
    }
  })
  observeEvent(input$BadPeak,{
    print("BAdPeak")
    dbn <- dbtaName(AnalyzedTransitions(),dbpath())
    
    for(dbni in dbn){
      SavePrecursorRating(rating = "bad",db = dbpath(),dbtaName = dbni )
    }  })
  observeEvent(input$notdefinedPeak,{
    print("SavePrecursor")
    dbn <- dbtaName(AnalyzedTransitions(),dbpath())
    
    for(dbni in dbn){
      SavePrecursorRating(rating = "not defined",db = dbpath(),dbtaName = dbni )
    }  })
  output$PrecursorQuality <- renderText({
    print("Quality")
    PrecursorQuality <- list(rating = "nodefined")
    input$notdefinedPeak
    input$BadPeak
    input$GoodPeak
    validate(need(file.exists(dbpath()),"No Path Set"),
             need(is.list(AnalyzedTransitions()),"No Transitions available.")
    )
    print("Waiting for names")
    dbf <- paste("RATING",dbtaName(AnalyzedTransitions(),dbpath()),sep = "_")
    print("Received Names")
    
    PrecursorQuality <- ""
    if(any(!is.na(match(dbf,dblistT(dbpath()))))){
      try(PrecursorQuality <- sapply(dbf,function(dbfi){    try(PrecursorQuality<- dbread(dbfi,dbpath())$rating)}))
      PrecursorQuality <- paste(unique(PrecursorQuality),collapse = " & ")
    }
    return(paste("Rating:",unique(PrecursorQuality)))
  })
  
  # Export ------
  observeEvent(input$Export,{
    validate(need(file.exists(input$mainPath),"No valid directory."))
    
    setwd(input$mainPath)
    dir.create("settings")
    dir.create("export")
    # setwd("export")
    # unlink("TransitionTable.txt")
    # unlink("PrecursorTable.txt")
    # setwd("../")
    if(length(session) > 0){
      progress <- Progress$new(session,min=0, max=length(analyzedList()))
      on.exit(progress$close())
      progress$set(message = 'Starting Export',
                   detail = "",value = 0)
    }
    
    fipath = "./Export/TransitionList.txt"        
    try(unlink(fipath),silent = T)
    # pdf("export/OutputResults.pdf")
    huha <- analyzedList()
    print("Export Setting Paths")
    
    print("Export: iterating entries")
    anaexport <<- analyzedList()
    dbp <<- dbpath()
    FDRpep <- input$FDRpep
    MinPeakWidth_vec <<- input$MinPeakWidth
    MaxPeakWidth_vec <<- input$MaxPeakWidth
    supersmooth_bw_set_vec <<- input$supersmooth_bw_set
    ApplyMaximumWidth_vec <<-   input$ApplyMaximumWidth
    Requantify_Priority <<- input$Requantify_Priority
    for(i in 1:length(anaexport)){
      fifu <- anaexport[i]
      cat("\rworking on ",fifu)
      validate(need(file.exists("settings"),"No valid directory."))
      
      setwd("settings")
      finame <- gsub(".txt$","",fifu)
      # print("Loading Stuff")
      if(file.exists(finame)){
        load(finame)
      }else{
        inputList <- inputListStdSet
      }
      
      # print("Export: Creating DP Table")
      # Creating DP Table
      inputList$FDR <- FDRpep
      if(length(session) > 0){
        progress$set(message = 'Export, working on:',
                     detail = finame,value = i)
      }
      
      setwd("../")
      if(i==1){
        unlink("./export/Peaks.txt")
        
      }
      
      PeaksName <- dbtaNameExport(dbp,ana = fifu)
      PeaksNameRatin <- paste("RATING_",PeaksName,sep = "")
      
      
      if(any(PeaksNameRatin== dblistT(dbp))){
        # print("check Rating")
        PrecursorRating <- dbread(PeaksNameRatin,db = dbp)$rating
        # dbread(PeaksName,dbpath)
        # print("Found Rating")
        if(length(PrecursorRating) == 0){
          PrecursorRating <- "not set"
        }
      }else{
        # print("not set Rating")
        PrecursorRating <- "not set"
      }
      
      
      # check if Peptide (including or excluding PW) was already ana in the database  
      TYPIES <- which(gsub(".$","",PeaksName)==gsub(".$","",dblistT(dbp)))
      Lo <- dblistT(dbp)[TYPIES]
      
      LoopCount <- 1
      if(length(TYPIES) >1){
        LoopCount <- length(TYPIES)
      }
      # Loop is to account several Peak tables caused by different checked PeakWidths in the database, need to manually be selected afterwards
      # print("Export: PeakCount")
      try(rm("tempa"))
      # dbp <- dbpath()
      
      for(PeakCount in 1:LoopCount){
        tempa <<- dbread(x = fifu,dbp)
        if(length(tempa$RF_Scores)==0){
          tempa$RF_Scores <- 0
        }
        if(length(tempa$DL_Scores)==0){
          tempa$DL_Scores <- 0
        }
        
        if(dim(tempa)[1]>=minthresh){
          
          if(length(TYPIES) > 0){
            # print("FOUND")
            # tempa$SCA_mz_fdr <- NULL # compatibility with datasets analyzed with different versions
            # tempa$SCAfdr <- NULL# compatibility with datasets analyzed with differents
            # tempa$FDR_DL_all <- NULL# compatibility with datasets analyzed with differents
            # tempa$FDR_RF_all <- NULL# compatibility with datasets analyzed with differents
            # 
            # tempa$DL_Scores <- -10
            cl <- class(try(DP <- dbread(Lo[PeakCount],dbp)))
            PeakType <- "Manual"
            pw <- inputList$PeakWidth
            # pw <- as.numeric(strsplitslot(Lo[PeakCount],8,"_"))
            # tempa <- dbread(x = fifu,dbpath())#
            
          }else{
            # dbp <- dbpath()
            tempa$RT_Used <- tempa$RT_min
            tempa$FDR <- tempa$FDR_RF_all
            FDRcut <- FDRpep
            CANDIDATE_RT <- RTcandidatesFun(list(tempa),FDRcut)
            LI <- list(tempa)
            names(LI) <- fifu
            # DPlist_I <- DetectPeakWrapper(ana = ana,CANDIDATE_RT = CANDIDATE_RT,dbp = dbp,
            #                             RetentionTimeWindow = RTwin,QType = "Intensities",Reanalysis = T,
            #                             RT_BASED_onbestScore = T,
            #                             MinPeakWidth=input$MinPeakWidth,
            #                             MaxPeakWidth=input$MaxPeakWidth,
            #                             supersmooth_I_set = input$supersmooth_I_set,
            #                             supersmooth_bw_set = input$supersmooth_bw_set,
            #                             ApplyMaximumWidth = input$ApplyMaximumWidth
            # )
            tempList <- list(tempa)
            names(tempList) <- fifu
            tempList <- tempList
            DPlist_XIC <- DetectPeakWrapper(ana = tempList,CANDIDATE_RT = CANDIDATE_RT,dbp = dbp,
                                            RetentionTimeWindow = RTwin,QType = "Intensities",
                                            Reanalysis = F,
                                            RT_BASED_onbestScore = T,
                                            MinPeakWidth=MinPeakWidth_vec,
                                            MaxPeakWidth=MaxPeakWidth_vec,
                                            supersmooth_I_set = F,#input$supersmooth_I_set,
                                            supersmooth_bw_set = supersmooth_bw_set_vec,
                                            ApplyMaximumWidth =  ApplyMaximumWidth_vec,
                                            Requantify_Priority= Requantify_Priority
                                            
            )

            # DPlist_I <- DetectPeakWrapper(LI,CANDIDATE_RT,dbp,RetentionTimeWindow,"Intensities",Reanalysis = T)
            # DPlist_XIC <- DetectPeakWrapper(LI,CANDIDATE_RT,dbp,RetentionTimeWindow,"XIC",Reanalysis = T)
            # DP <- rbind(DPlist_I[[1]],DPlist_XIC[[1]])
            DP <- data.table(DPlist_XIC[[1]],stringsAsFactors = F)
            
            tempa$PEP[is.na(tempa$PEP)] <- 1
            # cl <- class(try(DP <- dbread(Lo[PeakCount],dbp)))
            cl <- "everythingok"
            # cl <- class(try(DP <- PeakFun(tempa,RetentionTimeWindow = inputList$RetentionTimeWindow,alpha = input$FDRpep,PeakWidth = inputList$PeakWidth,SimilarRT = F,session = session)))
            PeakType <- "Automated"
            pw <- inputList$PeakWidth
          }
          try({
            # tlnx <<- tlnx
            try(transname <- gsub("^mz","TRANS",fifu))
            # print(transname)
            dbl <- dblistT(dbp)
            # print(dbl)
            if(any(dbl==transname)){
              # print("FoundTrans")
              dbp <- dbpath()
              TransSaved <- NULL
              try(selec <- unlist(dbread(gsub("^PEAKS","TRANS",transname),db =dbp)))
              transis <- colnames(tempa)[1:(which(colnames(tempa) == "charge")-1)]
              notused <- which(is.na(match(transis,selec)))
              if(length(notused)> 0){
                showNotification("Excluding Fragments", duration = 1, closeButton =  TRUE, type = "warning")
                tempa <- tempa[,-notused]
              }
              
            }
          })
          # Writing Out FragmentTable
          try({
            tempa_dt <- data.table(tempa)
            DPdt  <- data.table(DP)
            # tempa <- cbind(tempa,DPdt[match(tempa$rawfile,DP$rawfile),.(Q1,Q2,Q1align,Q2align)])
            tempa_dt[,Peak :=RT_min%between%DPdt[unlist(.BY)==rawfile,.(Q1,Q2)][1,],rawfile]
            tempa_dt <- tempa_dt[Peak == TRUE,]
            tempa_dt[,precursor_ID :=fifu]
            tempa_dt$Precursor_Rating <- PrecursorRating
            rf <- tempa_dt$rawfile
            tempa_dt$rawfile <- NULL
            tempa_dt$rawfile <- rf
            id <- 1:(grep("^charge$",names(tempa_dt))-1)
            
            # check order of COlumns:
            # tempa_dt <- tempa_dt
            # id <- id
            NAMES <- c(colnames(tempa_dt)[id],"charge","mz","index","RT_min","MatchCount","R","R_p","IntensitySum","SCAall","SCAcut","ZeroScore","DiffSum","FDR","PEP","SCA","CorrelationScoreMedian","CorrelationScoreMean","CorrelationScoreSD","CorrelationScoreSum","CorrelationScoreSumSpecial","CorrelationScoreSumSpecial2","Score","Peaks1","Peaks2","ScaBefore","ScaBefore2","ScaAfter","ScaAfter2","RF_Scores","DL_Scores","GBM_Scores","GLM_Scores","FDR_RFspl_peptide","FDR_RF_peptide","FDR_RF_all","FDR_DL_all","FDR_GLM_all","FDR_GBM_all","SCA_mz_fdr","SCAfdr","Peak","precursor_ID","Precursor_Rating","rawfile")
            NAMES <- unique(NAMES)
            MORDER <- match(NAMES,colnames(tempa_dt))
            if(any(is.na(MORDER))){
              try({showNotification(session,paste("missing columns in table",fifu))},silent = T)
              for(i in which((is.na(MORDER)))){
                tempa_dt <- cbind(tempa_dt,as.double(NA))
              }
              colnames(tempa_dt)[(dim(tempa_dt)[2]-(sum(is.na(MORDER))-1)):dim(tempa_dt)[2]] <- NAMES[is.na(MORDER)]
            }
            tempa_dt <- tempa_dt[,.SD,.SDcols=NAMES]
            tempa_dt <- setcolorder(tempa_dt,NAMES)
            
            HUI <- melt(tempa_dt,measure.vars = id,id.vars = setdiff(1:dim(tempa_dt)[2],id))
            HUI <- HUI[value >-1,]
            # if(dim(HUI)[2]==46){
            #   NAMES46 <<- colnames(HUI)
            # }
            # if(dim(HUI)[2]==47){
            #   NAMES47 <<- colnames(HUI)
            # }
            
            
            
            if(dim(HUI)[1]>0){
              if(file.exists("./export/Peaks.txt")){
                header_b <- F
              }else{
                header_b <- T
                
              }
              write.table(HUI,"./export/Peaks.txt",sep = "\t",row.names = F,quote = F,append = T,col.names  = header_b)
              
            }
          })
          
          
          
          
          # print(head(DP))
          if(cl != "try-error"){
            DPA <<- DP
            DP <- DP
            if(!all(is.na(DP$DL_Scores))){
              
              DP$QuantitationType[DP$QuantitationType=="intensity"] <-  "Intensities"
              fifu <- fifu
              fifu <- gsub("__","-",fifu)
              DP <- data.frame(DP)
              # trans <- DP[,1:(which(colnames(DP) == "Q1")-1)]
              # info <- DP[,(which(colnames(DP) == "Q1")):dim(DP)[2]]
              
              
              DPmelt <- melt(data.table(DP)[!is.na(DP$DL_Scores),],id.vars=(which(colnames(DP) == "Q1")):dim(DP)[2])
              data.table::setnames(DPmelt,"value","Intensity")
              data.table::setnames(DPmelt,"variable","Fragment")
              
              # xic <- trans[info$QuantitationType == "XIC",]
              # int <- trans[info$QuantitationType == "Intensities",]
              # 
              AddInf <- unlist(strsplit(fifu,"_"))
              AddInf[1] <- gsub("^mz","",AddInf[1])
              AddInf[length(AddInf)] <- gsub(".txt","",AddInf[length(AddInf)])
              if(length(AddInf)< 6){
                AddInf <- c(AddInf,rep(NA,6-length(AddInf)))
                
              }
              if(length(AddInf)> 6){
                AddInf[6] <- paste(AddInf[6:AddInf[length(AddInf)]],collapse = "_")
                AddInf <- AddInf[1:6]
              }
              # longtab <- c()
              # for(itr in 1:dim(trans)[2]){
              #   longtab <- rbind(longtab,cbind(xic[,itr],int[,itr],colnames(trans)[itr],info))
              # }
              # melt(trans)
              # print("longtab naming")
              # colnames(longtab)[1:3] <- c("XIC","Intensities","Match") 
              # rownames(longtab) <- NULL
              # longtab$Intensities[is.na(longtab$Intensities)] <- 0
              # longtab$XIC[is.na(longtab$XIC)] <- 0
              # Remove <-longtab$XIC == 0&longtab$Intensities == 0
              # Remove[is.na(Remove)] <- F
              # longtab <- longtab[!Remove,]
              # longtab$QuantitationType <- NULL
              longtab <- DPmelt
              MaIn <- matrix(AddInf,dim(longtab)[1],length(AddInf),byrow = T)
              MaIn <- cbind(MaIn,PeakType,pw,PrecursorRating)
              longtab$PeakType <- PeakType
              longtab$PrecursorRating <- PrecursorRating
              
              InfInsert <- data.table(t(matrix((AddInf))))
              longtab[,c("m/z","z","Sequence","Accession","Gene Symbol","Library_ScanID"):={
                InfInsert
              },]
              if(length(MaIn)>0){
                
                # colnames(MaIn) <- c("m/z","z","Sequence","Accession","Gene Symbol","Library_ScanID","PeakType","Peak Detection Width","Rating")
                # longtab <- cbind(longtab,MaIn)
                # longtab <- unique(longtab)
                # print("Writing")
                # AddInf <- unlist(strsplit(fifu,"_"))
                
                if(!file.exists(fipath)){
                  try(write.table(longtab,file = fipath,quote = F,row.names = F,append = F,sep = "\t"))
                }else{
                  try(write.table(longtab,file = fipath,quote = F,row.names = F,col.names = F,append = T,sep = "\t"))
                }
              }
            }
            
            
          }
        }
      }
      
    }    
    print("Export: longtab")
    try({
      # longtab <- fread(fipath,sep = "\t",stringsAsFactors = F)
      # # se <- c("XIC","Intensities","Peaks")
      # # print("subset1")
      # longtab$Gene.Symbol[is.na(longtab$Gene.Symbol)] <- "unknown"
      # longtab$Accession[is.na(longtab$Accession)] <- "unknown"
      # 
      # longtab[,max(Sequence),.(Sequence,`Gene Symbol`,Accession)]
      # try({
      #   CondenseRaw <- CondenseNames(uniRaw)
      #   CondenseRaw <- CondenseRaw[match(LoTab$rawfile,uniRaw)]
      # },silent = T)
      # 
      # 
      # try(write.table(LoTab,file = "./export/ProteinList.txt",quote = F,row.names = F,col.names = F,append = T,sep = "\t"))
      # # Calculate Ratios:
      try({
        peakspath <- "./export/Peaks.txt"
        if(file.exists(peakspath)){
          peaks <- fread(peakspath,sep ="\t",stringsAsFactors = F)
          peaks[,Precursor := gsub("_.*$","",precursor_ID)]
          peaks[,PrecursorInfo := substr(precursor_ID,gregexpr("_",precursor_ID)[[1]][1],nchar(precursor_ID))]
          peaks[,Sequence:= strsplit(precursor_ID,"_")[[1]][3],precursor_ID]
          peaks[,Label:= strsplit(Sequence,"#")[[1]][2],precursor_ID]
          peaks[,Sequence:= strsplit(Sequence,"#")[[1]][1],precursor_ID]
          
          peaks[,Charge:= strsplit(precursor_ID,"_")[[1]][[2]],precursor_ID]
          peaks[,Gene:= strsplit(precursor_ID,"_")[[1]][[5]],precursor_ID]
          
          if(any(is.na(peaks$mz))){
            peaks$mz <- as.numeric(gsub("^mz","",peaks$Precursor))
          }
          peaks <- peaks[!is.na(value),]
          Ratios <- PeaksProcessingLM(peaks)
          fwrite(Ratios,"./export/Ratios.txt",sep = "\t")
        }
        
      })
    })
    showNotification("Finished Export. You can find the tables in the export folder in your experiment location.",duration = NULL,type = "message")
  })
  
  
  
}



shinyApp(ui = ui,server = server)

