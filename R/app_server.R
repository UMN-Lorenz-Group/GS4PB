#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @importFrom utils read.csv read.table write.table write.csv as.roman
#' @importFrom doParallel registerDoParallel
#' @importFrom reticulate import py_module_available py_run_string
#' @importFrom stats cor model.matrix as.formula dist qchisq quantile setNames
#' @importFrom NAM snpQC
#' @importFrom rrBLUP mixed.solve A.mat
#' @importFrom bWGR emRR emBB emBL emCV
#' @importFrom BGLR BGLR Multitrait
#' @importFrom grDevices dev.off png
#' @importFrom graphics par plot.default
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom plotly plot_ly plotlyOutput
#' @importFrom shiny isolate
#' @importFrom shinyjs enable
#' @importFrom fs path_home
#' @importFrom callr r_bg
#' @importFrom DT DTOutput renderDT
#' @noRd


app_server <- function(input, output, session) {

  # Allow the user to select directories from the home directory or root
  #
  
  volumes <- c(
	  Home = if (dir_exists(fs::path_home())) fs::path_home() else "/home",
	  Root = "/",
	  Cyverse_Data_Store = if (dir_exists("/data-store")) "/data-store" else NULL
  )
  
  
  # volumes <- c(Home = fs::path_home(), "Root" = "/","Cyverse_Data_Store"="/data-store/")
  
  # Set up shinyFiles to allow directory selection
  shinyDirChoose(input, "directory", roots = volumes, session = session)
  
  # Reactive variable to store the selected directory path
  selected_dir <- reactive({
    req(input$directory)
    dir_path <- parseDirPath(volumes, input$directory)
    dir_path
  })
  
  # Display selected directory path
  output$dir_path <- renderText({
    selected_dir()
  })
  
 
 
  
### GenoData  
  
  GenoTas <- reactive({
    genoFile <- input$infileVCF
    ext <- tools::file_ext(genoFile$datapath)
    req(genoFile)
    validate(need(ext == "vcf", "Please upload a vcf file"))
    withProgress(message = 'Reading Data', value = 0, {
      getTasObj(genoFile$datapath)})
  })
  
  Geno <- eventReactive(input$infileVCF,{
    withProgress(message = 'Converting VCF to Dataframe', value = 0, {
      gt2d <- getGenoTas_to_DF(GenoTas())})
    
    gt2d

  })
  
  
  ## Tas to DF
  #Geno_DF <- reavtiveVal(NULL)
  Geno_DF <- reactive({
    # browser()
    if(!is.null(GenoTas())){getGenoTas_to_DF(GenoTas())}
  })
 
  
  temp_file0a <- reactiveVal('none')
  
 # Start the process in a separate R process
  
  observe({
    
    temp_file0a(tempfile())
    if(!is.null(Geno_DF())){
      sink(temp_file0a())
      cat(getGenoQCStats(Geno_DF()))
      sink()
    }
    
  })
  
  ### Test output 
  # Periodically read the file and update the UI
  mssgGenoStatsOut <- reactiveVal(NULL)
  
  output$messageGenoStats <- renderText({
    
    invalidateLater(1000, session) # Update every second
    if(file.exists(temp_file0a())){
      lines <- readLines(temp_file0a(), warn = FALSE)
      mssgGenoStatsOut(paste(lines, collapse = "\n"))
      return(paste(lines, collapse = "\n"))
    
    }else {
      return("Waiting for output...")
    }
  })
  
  
  
#### Target Table 

  TargetTab <- reactive({
    TargetFile <- input$infileTargetTable
    
    ext <- tools::file_ext(TargetFile$datapath)
    req(TargetFile)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    read.csv(TargetFile$datapath, header = input$header)
  })
  
  
  TargetIDs <- reactiveVal(NULL)
  
  observeEvent(input$infileTargetTable,{updateSelectInput(inputId="TargetIDCol",choices=colnames(TargetTab()))})
  
  
### TargetID Col 
  
  TargetIDCol <- reactive({
    req(input$TargetIDCol)
    input$TargetIDCol
  })
  
  observe({
    req(TargetTab())
    req(TargetIDCol)
    TargetIDs(as.character(TargetTab()[,TargetIDCol()]))
  })
  
  
  buildLib <- reactive({
    libFile <- input$inHapLib
    ext <- tools::file_ext(libFile$datapath)
    
    req(libFile)
    validate(need(ext == "phase", "Please upload a .phase file"))
    
    read.table(libFile$datapath, header = FALSE)
  })
  
  hapLibLoaded <- reactiveVal(FALSE)
  hapLibStats <- reactiveVal(NULL)
  
  observeEvent(input$inHapLib,{ 
    #browser()
     shinyjs::enable("impute_APIdata")
     hapLibStats(getLibStats(buildLib()))
     hapLibLoaded(TRUE)
     write.table(buildLib(),"lib.phase",quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)
    
  })
 
    

#### Target Table 
  
  temp_file0d <- reactiveVal('none')
  
  # Start the process in a separate R process
  
  observe({
    
    temp_file0d(tempfile())
    if(!is.null(TargetTab())){
      sink(temp_file0d())
      cat(paste("Table with information on ",nrow(TargetTab())," Target lines",sep=""))
      sink()
    }
    
  })
  
  ### Test output 
  # Periodically read the file and update the UI
  mssgTargetStatsOut <- reactiveVal(NULL)
  
  output$messageTargetStats <- renderText({
    
    invalidateLater(1000, session) # Update every second
    if(file.exists(temp_file0d())){
      lines <- readLines(temp_file0d(), warn = FALSE)
      mssgTargetStatsOut(paste(lines, collapse = "\n"))
      
      return(paste(lines, collapse = "\n"))
      
    }else {
      return("Waiting for output...")
    }
  })
  
  
  
  
  
 
### Filter Data
  
  ### Filter 1  
  observeEvent(input$infileVCF, {
    updateNumericInput(inputId = "siteMinCnt",value= 0.80,min =0, max=1)}) 
  
  
  siteMinCnt <- reactive(input$siteMinCnt)
  MAF <- reactive(input$MAF)
  
 
  
  GenoFilt1 <- eventReactive(input$FilterSites,{
    withProgress(message = 'Filtering Sites', value = 0, {
      getFilteredSitesGenoData(GenoTas(),round(siteMinCnt()*(ncol(Geno())-5),digits = 0),MAF())}) 
  })
  
  
 
  
######## 
  
  temp_file0b <- reactiveVal('none')
  
  
# Start the process in a separate R process
  
  observeEvent(input$FilterSites, {
    
    temp_file0b(tempfile())
    if(!is.null(GenoFilt1_DF())){
      sink(temp_file0b())
      cat(getGenoQCStatsFilt1(Geno_DF(),GenoFilt1_DF(),siteMinCnt(),MAF()))
      sink()
    }

  })
  
### Test output 
# Periodically read the file and update the UI

  mssgGenoFilt1StatsOut <- reactiveVal(NULL)
  
output$messageGenoFilt1 <- renderText({
  
  invalidateLater(1000, session) # Update every second
  if(file.exists(temp_file0b())){
    lines <- readLines(temp_file0b(), warn = FALSE)
    mssgGenoFilt1StatsOut(paste(lines, collapse = "\n"))
    return(paste(lines, collapse = "\n"))
  }else {
    return("Waiting for output...")
  }
})

  
  
### Filter 2
  
  minNotMissing <- reactive(input$minNotMissing)
  
  setTasGenoFilt1 <- reactive(input$setGenoFilt1Tas)
  
  
  GenoFilt2 <- eventReactive(input$FilterTaxa,{
    
    if(setTasGenoFilt1()== FALSE){        
      withProgress(message = 'Filtering Taxa', value = 0, {getFilteredTaxaGenoData(GenoTas(),minNotMissing()) })
    }else if(setTasGenoFilt1()== TRUE){
      withProgress(message = 'Filtering Taxa', value = 0, {getFilteredTaxaGenoData(GenoFilt1(),minNotMissing())})
    }
  })
  
  setTasGenoFilt2 <- reactive(input$setGenoFilt2Tas)
  
######## 
  
  temp_file0c <- reactiveVal('none')
  
# Start the process in a separate R process
  
  observeEvent(input$FilterTaxa,{
    temp_file0c(tempfile())
    if(!is.null(GenoFilt2_DF())){
      sink(temp_file0c())
      cat(getGenoQCStatsFilt2(GenoFilt1_DF(),GenoFilt2_DF(),minNotMissing()))
      sink()
    }
  })
  
### Test output 
# Periodically read the file and update the UI
  
  mssgGenoFilt2StatsOut <- reactiveVal(NULL)
  
  output$messageGenoFilt2 <- renderText({
    
    invalidateLater(1000, session) # Update every second
    if(file.exists(temp_file0c())){
      lines <- readLines(temp_file0c(), warn = FALSE)
      mssgGenoFilt2StatsOut(paste(lines, collapse = "\n"))
      return(paste(lines, collapse = "\n"))
    }else {
      return("Waiting for output...")
    }
  })
  
  
 
  
  GenoFilt1_DF <- reactive({
    #browser()
    if(!is.null(GenoFilt1())){ getGenoTas_to_DF(GenoFilt1())}
  })
  
  GenoFilt2_DF <- reactive({
    if(!is.null(GenoFilt2())){ getGenoTas_to_DF(GenoFilt2())}
  })
  
 
  
  
### Imputation
  
  FiltGeno <- reactive({ 
    if(setTasGenoFilt1()== FALSE  & setTasGenoFilt2()== FALSE){        
      GenoTas()
    }else if(setTasGenoFilt1()== TRUE  & setTasGenoFilt2()== FALSE){        
      GenoFilt1()
    }else if(setTasGenoFilt2()== TRUE){
      GenoFilt2()
    }
    
  })
  
  l<- reactive(input$l)
  k <- reactive(input$k)
  
  nN <- reactive(input$nN)
  Dist <- reactive(input$Dist)
  impMethod <- reactive(input$imputeMet)
  
  
  GenoImp <-  eventReactive(input$Impute,{
    
    withProgress(message = 'Imputing Genotypic Scores', value = 0, {
      
      print(impMethod())
      if(impMethod()=='Numeric'){ 
        
        df<- getImputedData_Num(FiltGeno(),nN(),Dist())
        
      }else if(impMethod()=='LDKNNI'){ 
        df <-  getImputedData_LDKNNI(FiltGeno(),l(),k())
        
      }
      
      df
    })
    
  },ignoreNULL = TRUE)
  
  
  GenoImp_DF1 <- reactiveVal(NULL)
  
  # GenoImp_DF1 <- eventReactive(input$Impute,{
  #   getGenoTas_to_DF(GenoImp())
  #})
  
  
   observeEvent(input$Impute,{
     GenoImp_DF1(getGenoTas_to_DF(GenoImp()))
   })
  
  ######## 
  
  temp_file0e <- reactiveVal('none')
  
  
  # Start the process in a separate R process
  
  observeEvent(input$Impute, {
    
    temp_file0e(tempfile())
    if(!is.null(GenoImp_DF1())){
      sink(temp_file0e())
      cat(getGenoImp1Stats(Geno_DF(),GenoImp_DF1()))
      sink()
    }
    
  })
  
  ### Test output 
  ## Periodically read the file and update the UI
  
  mssgImpGeno1StatsOut <- reactiveVal(NULL) 
  
  output$messageImpGeno1 <- renderText({
    
    invalidateLater(1000, session) # Update every second
    if(file.exists(temp_file0e())){
      lines <- readLines(temp_file0e(), warn = FALSE)
      mssgImpGeno1StatsOut(paste(lines, collapse = "\n"))
      return(paste(lines, collapse = "\n"))
    }else {
      return("Waiting for output...")
    }
  })
  
  
  
  ### if you use only reactive, this will throw an error ncol(GenoImp_DF())-5
  
  #genoImpHead <- eventReactive(input$Impute,{paste("Genotype Table with ",ncol(GenoImp_DF())-5," lines and ",nrow((GenoImp_DF()))," markers",sep="")})
  
  # genoImpHead <- reactive({paste("Genotype Table with ",ncol(GenoImp_DF())-5," lines and ",nrow((GenoImp_DF()))," markers",sep="")})
  
 
  
  observeEvent(input$imputeMet,{
    
    if(impMethod()=="AlphaPlantImpute"){ 
      
      reticulate::py_run_string("import sys")
      
      API_exists <- reticulate::py_module_available("alphaplantimpute2")
      print(paste("API Available",API_exists,sep=" "))
      
      output$APIUI <- renderUI({
       
        # Check if AlphaPlantImpute is selected
  
        fluidPage(
          
          fluidRow(
           box(width = 12,title="AlphaPlantImpute",solidHeader=TRUE,status="primary",
              tags$strong(tags$h4("AlphaPlantImpute (Gonen et al., 2018) implements imputation in two steps:"))
            )
          ),
          
          tags$br(),
          
          fluidRow(
            box(
              width = 12,
              title= tags$strong("Step 1: Building the haplotype library"),
              solidHeader=TRUE,
              status="primary",
              collapsible= TRUE,
              tags$ul(
                tags$li("This step constructs a haplotype library using the uploaded and filtered genotype file (available in the left panel)."),
                tags$li("Set the parameters and build the library."),
                tags$li("Once the library construction is complete, the 'Impute' button will become active."),
                tags$li("If you already have a prebuilt library in `.phase` format, you can upload it to skip the time-consuming library construction step."),
                tags$li("Uploading a `.phase` file also activates the 'Impute' button.")
              )
            )
          ),
          
          tags$br(),
          
          fluidRow(
            box(
              width = 12,
              title=tags$strong("Step 2: Genotype file imputation"),
              solidHeader=TRUE,
              status="primary",
              collapsible=TRUE,
              tags$ul(
                tags$li("The imputation process uses the constructed or uploaded haplotype library."),
                tags$li("If a founders file is provided, a faster pedigree-based imputation will be performed."),
                tags$li("If no founders file is uploaded, a population-based imputation is applied.")
              )
            )
          ),
          
        tags$br(),
         
        fluidRow(
            box(width=6,title=tags$strong(tags$h3("Build Haplotype Library ")),solidHeader = TRUE,status="primary",
                tags$br(),
                numericInput(inputId="HDthresh",label = "Enter non-missing threshold for high density genotypes",value=0.9,min=0.5,max=1),
                numericInput(inputId="nHap",label = "Enter number of haplotypes",value=20,min=2,max=100),
                numericInput(inputId="nSampRnds",label = "Enter number of sample rounds",value=5,min=2,max=100),
                tags$br(),
                actionButton("build_library", "Build Library"),
              
            ),
            box(width=6,title=tags$strong(tags$h3("Impute Genotype Table ")),solidHeader = TRUE,status="primary",
                tags$br(),
                fileInput("founder_file", "Select Founder Files", multiple = FALSE,accept = ".txt"),
                tags$br(),
                tags$br(),
                shinyBS::bsButton("impute_APIdata",label="Impute Data",style="primary",disabled=TRUE),
            ),
        ),
        tags$br(),
        
        
       fluidRow(
            
        box(width=6,title="API Build Status",solidHeader=TRUE,status="info",
                   verbatimTextOutput("message")
        ),
        box(width=6,title="API Impute Status",solidHeader=TRUE,status="info",  
                   verbatimTextOutput("message2")
        )
       ),
      tags$br(),
      
      fluidRow(
        
        box(width=12,title="APImputed Geno Stats",solidHeader=TRUE,status="info",
             verbatimTextOutput("messageAPIImpStats")
        )
      ), 
     )
     })
      
    } 
    
  })
  
###
  
  observeEvent(input$imputeMet,{
    if(impMethod()=="AlphaPlantImpute" && (!is.null(FiltGeno()))){
      getGenoData_API(FiltGeno())
    }
  })
  
  ####  
  
  founder <- reactive({
    
    # Attempt to retrieve file information
    founderFile <- input$founder_file
    
    # Proceed only if a file is uploaded
    if (!is.null(founderFile)) {
      # Extract file extension
      ext <- tools::file_ext(founderFile$datapath)
      
      # Validate file type
      validate(need(ext == "txt", "Please upload a txt file"))
      
      # Return the file path if the file is valid
      return(founderFile$datapath)
    }
    
    # Return NULL if no file is uploaded
    return(NULL)
  })
  
  #### 
  temp_file <- reactiveVal()
  nHap <- reactive(input$nHap)
  nSampRnds <- reactive(input$nSampRnds)
  HDthresh <- reactive(input$HDthresh)
  processComplete <- reactiveVal(FALSE)
   
  rProcess <- eventReactive(input$build_library,{
    # Path for the temporary file
    temp_file(tempfile())
    
    # Start the process in a separate R process
    callr::r_bg(function(nHap, nSampRnds, HDthresh, temp_file) {
      library(reticulate)
      library(GS4PB)
      # reticulate::use_virtualenv("./renv/python/virtualenvs/renv-python-3.12", required = TRUE)
      
      reticulate::use_condaenv("GS4PB_CondaEnv", required = TRUE)  
      
      sys <- reticulate::import("sys")
      api2 <- reticulate::import("alphaplantimpute2.alphaplantimpute2")
      
      sys$argv <- c('alphaplantimpute2', '-createlib', '-out', 'lib',
                    '-genotypes', 'current_GenoTable.genotypes',
                    '-n_haplotypes', as.character(nHap),
                    '-n_sample_rounds', as.character(nSampRnds),
                    '-hd_threshold', as.character(HDthresh),
                    '-seed', '42')
      
      sink(temp_file)
      tryCatch({
        api2$main()
      }, finally = {
        sink(NULL)
      })
    }, args = list(nHap(), nSampRnds(), HDthresh(), temp_file()), stdout = temp_file(), stderr = temp_file())
    
  })
   
###
  mssgAPIBuildOut <- reactiveVal(NULL)
  
# Reactive expression to check process status and read temp file

  processStatus <- reactive({

    invalidateLater(1000, session)  # Check every second
    if(!is.null(rProcess()) && !rProcess()$is_alive()) {
      rProcess()$kill()  # Ensures the background process is terminated
    }
    # Check if the process is running

    if(!is.null(rProcess()) && rProcess()$is_alive()) {
      if(file.exists(temp_file())){
        lines <- readLines(temp_file(), warn = FALSE)
        return(paste(lines, collapse = "\n"))
      }
      return("Waiting for output...")
    }
    # Check if the process has completed
    if(!is.null(rProcess()) && !rProcess()$is_alive()){
      if(file.exists(temp_file())){
        lines <- readLines(temp_file(), warn = FALSE)
        txt <- paste(lines, collapse = "\n")
        # Mark process as complete and update output
        processComplete(TRUE)
        shinyjs::enable("impute_APIdata")
      	# Ensure the background process is terminated
        mssgAPIBuildOut(paste(txt, "Build Completed", collapse = "\n"))
        return(paste(txt, "Build Completed", collapse = "\n"))
      }
      return("Waiting for output...")
    }
    return("Error determining process status")
 })


  ##
   output$message <- renderText({
     invalidateLater(1000, session)
     if(hapLibLoaded()) {
       mssgLibStats <- hapLibStats()
      return(mssgLibStats)
     }else{
       processStatus()      	
     }	     
   })
  
######   
  
  
  temp_file2 <- reactiveVal()
  process2Complete <- reactiveVal(FALSE)
  outputRead <- reactiveVal(FALSE)
  
  
# Start the background process and return rProcess

 rProcess2 <- eventReactive(input$impute_APIdata,{
    
    #browser()
    ## Extract necessary input values
    
    founder_file <- founder()
    
    # Path for the temporary file
    
    temp_file2(tempfile())
    
    # Start the process in a separate R process
    callr::r_bg(function(founder_file,temp_file){
      library(reticulate)
      library(GS4PB)
      reticulate::use_condaenv("GS4PB_CondaEnv")
      #reticulate::use_virtualenv("./renv/python/virtualenvs/renv-python-3.12", required = TRUE)
      sys <- import("sys")
      api2 <- import("alphaplantimpute2.alphaplantimpute2")
      
      if(is.null(founder_file)){
        sys$argv <- c('alphaplantimpute2', '-impute', '-out', 'imputed_out',
                      '-genotypes', 'current_GenoTable.genotypes',
                      '-libphase','lib.phase'
        )
      }else {
        sys$argv <- c('alphaplantimpute2', '-impute', '-out', 'imputed_out',
                      '-genotypes', 'current_GenoTable.genotypes',
                      '-founders', as.character(founder_file),
                      '-libphase','lib.phase'
        )
      }
      
      sink(temp_file)
      tryCatch({
        api2$main()
      }, finally = {
        sink(NULL)
      })
    }, args = list(founder_file,temp_file2()), stdout = temp_file2(), stderr = temp_file2())
    
  })
  
  mssgAPIImputeOut <- reactiveVal(NULL)
  
  # Reactive expression to check process status and read temp file

  processStatus2 <- reactive({
    
    invalidateLater(1000, session)  # Check every second
    
    
    if(!is.null(rProcess2()) && !rProcess2()$is_alive()) {
      rProcess2()$kill()  # Ensures the background process is terminated
    }
    
    # Check if the process is running
    
    if(!is.null(rProcess2()) && rProcess2()$is_alive()) {
      if(file.exists(temp_file2())){
        lines <- readLines(temp_file2(), warn = FALSE)
        return(paste(lines, collapse = "\n"))
      }
      return("Waiting for output...")
    }
    
    # Check if the process has completed
    
    if(!is.null(rProcess2()) && !rProcess2()$is_alive()){
      if(file.exists(temp_file2())){
        lines <- readLines(temp_file2(), warn = FALSE)
        txt <- paste(lines, collapse = "\n")
        # Mark process as complete and update output
        process2Complete(TRUE)
        # Ensure the background process is terminated
        mssgAPIImputeOut(paste(txt, "Imputation Completed", collapse = "\n"))
        return(paste(txt, "Imputation Completed", collapse = "\n"))
      }
      return("Waiting for output...")
    }
    
    return("Error determining process status")
  })
  
  output$message2 <- renderText({
    processStatus2()
  })
  
  # Observe and finalize output processing
  observe({
    req(process2Complete())
    outputRead(TRUE)  # Mark output as read
    
  })
  
###### STOutPlots 

  output$STPlotsUI <- renderUI({
    
    
      if(nSelTraits()==1){
        div(
          plotlyOutput("plots", height = "400px", width = "100%"),
          style = "display: flex; flex-wrap: wrap; justify-content: center;" # Prevents horizontal scrolling issues
        )
          
      }else if(nSelTraits()>1){
       div(
        fluidRow(
          tagList(lapply(1:nSelTraits(), function(i) {
            plotlyOutput(outputId = paste0("plot", i), height = "400px", width = "100%")
          }))
        ),
        style = "display: flex; flex-wrap: wrap; justify-content: center;" # Prevents horizontal scrolling issues
       )
      }
  })
  
 
  output$MTPlotsUI <- renderUI({
    
    div(
      fluidRow(
        tagList(
          lapply(1:nSelTraits(), function(i) {
            plotlyOutput(outputId = paste0("plotMT", i), height = "400px", width = "100%")
          })
         )
      ),
      style = "display: flex; flex-wrap: wrap; justify-content: center;"  # Prevents horizontal scrolling issues
    )
    
  })
  
   
  
#####  
  GenoImp_DF2 <- reactiveVal(NULL)
  
  observe({
    req(outputRead()) # Ensure output has been marked as read
    ImpGeno <- getImpGenoData_API()  # Retrieve imputed genotype data
    GenoImp_DF2(ImpGeno)
    
  })

### Get Imputed Geno Stats 
  
  temp_file0e2 <- reactiveVal('none')

# Start the process in a separate R process
  
 
  observe({ 
    
    temp_file0e2(tempfile())
    if(!is.null(GenoImp_DF2())){
      sink(temp_file0e2())
      cat(getGenoImp1Stats(Geno_DF(),GenoImp_DF2()))
      sink()
    }
  })
  
### Test output 
## Periodically read the file and update the UI
  
  mssgImpGeno2StatsOut <- reactiveVal(NULL) 
  
  output$messageAPIImpStats <- renderText({
    
    invalidateLater(1000, session) # Update every second
    if(file.exists(temp_file0e2())){
      lines <- readLines(temp_file0e2(), warn = FALSE)
      mssgImpGeno2StatsOut(paste(lines, collapse = "\n"))
      return(paste(lines, collapse = "\n"))
    }else {
      return("Waiting for output...")
    }
  })
  

### Set Imputed Geno Data to geno data imputed using TASSEL / API methods  
  
  GenoImp_DF <- reactive({
    
    if(!is.null(GenoImp_DF1()) && nrow(GenoImp_DF1()) > 0){
        # GenoImp_DF1 is not NULL and has data
        return(GenoImp_DF1())
       
    }else if(!is.null(GenoImp_DF2()) && nrow(GenoImp_DF2()) > 0) {
        # GenoImp_DF2 is not NULL and has data, and GenoImp_DF1 was NULL or had no data
        return(GenoImp_DF2())
       
    }else {
      # Both are NULL or have no data, return NULL or a default value/data frame
      return(NULL)  # Or any default value you deem appropriate
    }
  })
  
  
  # genoImpHead <- reactive(paste("Genotype Table with ",ncol(GenoImp_DF())-5," lines and ",nrow((GenoImp_DF()))," markers",sep=""))
  # GenoImp_DF <- eventReactive(input$Impute,{getGenoTas_to_DF(GenoImp(),Geno())})
  # GenoImp_DF <- reactive(getGenoTas_to_DF(GenoImp(),Geno()))
  # GenoImp_DF <- reactive(getGenoTas_to_DF(GenoImp()))
  # GenoImp_DF <- reactive(as.data.frame(as.matrix(GenoImp())))
  
  
  
### Merge and Process Data
  
  setTasImpGeno <- reactive(input$setGenoImpTas) 
  
### Check for which genotype matrix to set   
  
  GenoPre <- reactive({
    
    if( setTasGenoFilt1() == FALSE  && setTasGenoFilt2() == FALSE && setTasImpGeno()== FALSE ){        
      Geno_DF()
    }else if(setTasGenoFilt1() == TRUE  && setTasGenoFilt2() == FALSE && setTasImpGeno()== FALSE){        
      GenoFilt1_DF()
    }else if(setTasGenoFilt2() == TRUE && setTasImpGeno()== FALSE){
      GenoFilt2_DF()
    }else if(setTasImpGeno()== TRUE){
      GenoImp_DF()
    }
    
  })
 

  
  
###### 
  
  #Pheno 
  
  ## PhenoSE
  
  Pheno <- reactive({
    
    phenoFile <- input$infileBLUEsSE
    
    ext <- tools::file_ext(phenoFile$datapath)
    req(phenoFile)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    read.csv(phenoFile$datapath, header = input$header)
    
  })
  
###
  
  observeEvent(input$infileBLUEsSE, {
    updateSelectInput(inputId="trait",choices = colnames(Pheno())[2:ncol(Pheno())])})
  
  observeEvent(input$infileBLUEsSE, {
    updateSelectInput(inputId="fixed",choices = c("NULL",colnames(Pheno())[2:ncol(Pheno())]),selected="NULL")})
  
###  


  observeEvent(input$infileBLUEsSE, {
    shinyWidgets::updateNumericRangeInput(inputId ="traitColsNum",value= c(2,ncol(Pheno())))})
  
  observeEvent(input$infileBLUEsSE, {
   updateSelectInput(inputId="IDColSE",choices=colnames(Pheno()))})
  
  observeEvent(input$infileBLUEsSE, {
   updateSelectInput(inputId="strainSE",choices=colnames(Pheno()))})
 
##
  observeEvent(input$iPh,{
    showModal(modalDialog(
      title = "Information",
      "Format:",
      tags$ul(
        tags$li("The pheno file should have the following format: 'GermplasmId','Trait 1',Trait 2',..."),
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  
  observeEvent(input$iPhME,{
    showModal(modalDialog(
      title = "Information",
      "Format:",
      tags$ul(
        tags$li("The pheno file from METs should have the following format: 'uniqID','Strain','Loc','Year','Test/Other factors','Trait 1','Trait 2',..."),
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  #
  
 
###
  
 
  IDColSE <- reactive(input$IDColSE)
  StrainSE <- reactive(input$strainSE)
  
 
  TraitColsSE <- reactiveVal(NULL) #(c(input$traitColsNum[1]:input$traitColsNum[2]))
  TraitsPhSE <- reactiveVal(NULL) # (colnames(Pheno())[TraitColsSE()])
  TraitsPhSEVec <- reactiveVal(NULL) # (paste(TraitsPhSE()," ",sep="",collapse=","))
 
  
  # TraitsPhSEVec(paste(TraitsPhSE()," ",sep="",collapse=","))
  
####  
  
  temp_file3d <- reactiveVal('none')
  
  observeEvent(input$traitColsNum,{
    
     
    TraitColsSE(c(input$traitColsNum[1]:input$traitColsNum[2]))
    TraitsPhSE(colnames(Pheno())[TraitColsSE()])
    TraitsPhSEVec(paste(TraitsPhSE()," ",sep="",collapse=","))
    
    # Path for the temporary file
    temp_file3d(tempfile())
    
    if(length(IDColSE())>0 & !is.null(Pheno())){
      sink(temp_file3d())
      lenIDColSE <- reactive(length(unique(Pheno()[,IDColSE()])))
      if(lenIDColSE()== nrow(Pheno())){
        cat(paste("The Unique ID column has ",lenIDColSE()," unique number of entries and matches the total number of entries in the data table",sep=""))
        cat("\n")
        cat(paste("The phenotypes data table has data for ",length(TraitColsSE())," including ",TraitsPhSEVec(),"\n",sep=""))
        cat("\n")
      }else{
        cat(paste("The Unique ID column has ",lenIDColSE()," unique number of entries, which doesn't match the total number of entries in the data table",sep=""))
        cat("\n")
        cat(paste("The phenotypes data table has data for ",length(TraitColsSE())," including ",TraitsPhSEVec(),"\n",sep=""))
        cat("\n")
      }
      sink()
    }
  })  
  
  ### Test output 
  # Periodically read the file and update the UI
  
  mssgPhenoStatsOut <- reactiveVal(NULL)
  
  output$messagePhSE <- renderText({
    
    invalidateLater(1000, session) # Update every second
    if(file.exists(temp_file3d())){
      lines <- readLines(temp_file3d(), warn = FALSE)
      mssgPhenoStatsOut(paste(lines, collapse = "\n"))
      return(paste(lines, collapse = "\n"))
    }else {
      return("Waiting for output...")
    }
  })
  
  # 
  # phenoHead <- eventReactive(input$infileBLUEsSE,{ paste("Table with ",nrow(Pheno())," lines and ",ncol(Pheno())-1," traits",sep="")})
  # output$PhenoHeader <- renderText({phenoHead()})
  # output$PhenoTable <- renderTable({as.data.frame((Pheno())[1:5,1:3])})

### PhenoME  
  
  PhenoME <- reactive({
    
    phenoFileME <- input$infileBLUEsME
    
    ext <- tools::file_ext(phenoFileME$datapath)
    req(phenoFileME)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    read.csv(phenoFileME$datapath, header = input$headerME)
    
  })
  
  
  observeEvent(input$infileBLUEsME, {
   updateSelectInput(inputId = "traitCols",choices = colnames(PhenoME()))
  })
  
  observeEvent(input$infileBLUEsME, {
    updateSelectInput(inputId="IDColME",choices = colnames(PhenoME()))
  })
  
  observeEvent(input$infileBLUEsME, {
    updateSelectInput(inputId="strainME",choices = colnames(PhenoME()))
  })
  
  StrainME <- reactive(input$strainME)
  
  TraitCols <- reactive(input$traitCols)

  observeEvent(input$traitCols, {
     updateSelectInput(inputId ="traitME",choices = colnames(PhenoME())[which(colnames(PhenoME()) %in% TraitCols())])
  })
  
  
  
  
#### 
   temp_file3 <- reactiveVal()
   IDColME <- reactive(input$IDColME)
   
   observeEvent(input$IDColME,{
     
     # Extract necessary input values
     # Path for the temporary file
     temp_file3(tempfile())
     
     # Start the process in a separate R process
    
     if(length(IDColME())>0 & !is.null(PhenoME())){
       sink(temp_file3())
       lenIDColME <- reactive(length(unique(PhenoME()[,IDColME()])))
       if(lenIDColME() == nrow(PhenoME())){ 
         cat(paste("The Unique ID column has ",lenIDColME()," unique number of entries and matches the total number of entries in the data table",sep=""))
       }else{
         cat(paste("The Unique ID column has ",lenIDColME()," unique number of entries, which doesn't match the total number of entries in the data table",sep=""))
       }
       sink()
      }
    })
   
### Test output 
# Periodically read the file and update the UI
     
     output$messagePhME <- renderText({
       
       invalidateLater(1000, session) # Update every second
       if(file.exists(temp_file3())){
         lines <- readLines(temp_file3(), warn = FALSE)
         mssgPhenoStatsOut(paste(lines, collapse = "\n"))
         return(paste(lines, collapse = "\n"))
       }else {
         return("Waiting for output...")
       }
     })
     
   
   
  
## Trait 
  
  nTraits <- eventReactive(input$infileBLUEsSE,{(ncol(Pheno())-1)})
  Trait <- reactive(input$trait)
  
### 
 
#####  
 
  temp_file3a <- reactiveVal('none')
  observeEvent(input$trait,{
    # Path for the temporary file
    temp_file3a(tempfile())
  
    if(length(Trait())>0 & !is.null(Pheno())){
      sink(temp_file3a())
      cat("Summary of selected trait values : \n")
      cat(paste("Trait",names(round(summary(Pheno()[,Trait()[1]]),digits=2)),collapse="\t"))
      cat("\n")
      for(nT in 1:length(Trait())){
        cat(paste(Trait()[nT],paste(round(summary(Pheno()[,Trait()[nT]]),digits=2),collapse="\t"),sep="\t"))
        cat("\n")
      }
      sink()
    }
  })
  
  ### Test output 
  # Periodically read the file and update the UI

  
  mssgTraitsStatsOut <- reactiveVal(NULL)
  
  output$messageTrtSE <- renderText({
    
    invalidateLater(1000, session) # Update every second
    if(file.exists(temp_file3a())){
      lines <- readLines(temp_file3a(), warn = FALSE)
      mssgTraitsStatsOut(paste(lines, collapse = "\n"))
      return(paste(lines, collapse = "\n"))
    }else {
      return("Waiting for output...")
    }
  })
  
  nSelTraits <- reactive(length(Trait()))
  
  # output$summaryHeader <- renderText({summaryHead()})
  # SummaryTxt <- eventReactive(input$trait,{ 
  #   do.call(rbind,lapply(Trait(),function(x) round(summary(Pheno()[,x]),digits=2)))
  # })
  
  # observeEvent(input$trait,{ 
  #   output$Summary <- renderTable({
  #     trTable <- cbind(unlist(Trait()),SummaryTxt())
  #     print(trTable)
  #   })
  # })
  
  
## Trait ME
  
  TraitMEIn <- reactive({input$traitME})
  TraitME <- reactiveVal(NULL)
  nSelTraitsME <- reactiveVal(NULL)
  # 
  observeEvent(input$traitME,{
   
    TraitME(TraitMEIn())
    if(!is.null(TraitMEIn())) {
      nSelTraitsME(length(TraitMEIn()))
    }
  })
 
  #nSelTraitsME <- reactive(length(TraitME()))
 
  IDColsME <- reactive({
     colnames(PhenoME())[which(!as.character(colnames(PhenoME())) %in% as.character(TraitCols()))]
  })

##### 
#####  
 
  
  phenoMEData <- reactive({
    print("phME In")
    #browser()
    if(!is.null(TraitME())){
      getPhenoMEData(PhenoME(),TraitME(),nSelTraitsME(),IDColsME(),StrainME())
    }else if(is.null(TraitME())){ 
      checkPhDataME("Select trait and merge data")
    }
  })
  
  
  
  
  temp_file3b <- reactiveVal('none')
  IDColME <- reactive(input$IDColME)
  
  observeEvent(input$traitME,{
   # browser()
    # Extract necessary input values
    # Path for the temporary file
    temp_file3b(tempfile())
    
    # Start the process in a separate R process
    
    if(length(TraitME())>0 & !is.null(PhenoME())){
        sink(temp_file3b())
        cat("Summary of selected trait values across all locations : \n")
        
        if(nSelTraitsME()==1){
         cat(paste(names(summary((PhenoME()[,TraitME()]))),"\t",sep=""))
         cat("\n")
         cat(paste(round(summary(as.numeric(PhenoME()[,TraitME()])),digits = 3),"\t",sep=""))
        }else if(nSelTraitsME()>1){ 
          cat(paste(names(summary((PhenoME()[,TraitME()[1]]))),"\t",sep=""))
          cat("\n")
          
          for(nSel in 1:nSelTraitsME()){
            
            cat(paste(round(summary(as.numeric(PhenoME()[,TraitME()[nSel]])),digits = 3),"\t",sep=""))
            cat("\n")
          }
        }
        sink()
    }
  })
  
  ### Test output 
  # Periodically read the file and update the UI
  
  
  output$messageTrtME <- renderText({
    
    invalidateLater(1000, session) # Update every second
    if(file.exists(temp_file3b())){
      lines <- readLines(temp_file3b(), warn = FALSE)
      mssgTraitsStatsOut(paste(lines, collapse = "\n"))
      return(paste(lines, collapse = "\n"))
    }else {
      return("Waiting for output...")
    }
  })
#   
#   
# ### Merge geno and pheno data    

 # MergeOR <- reactive({input$MergeGP | input$MergeGPME})
  
  isPhenoME <- reactive(input$chkPhenoME)

  mergedData <-  reactiveVal(NULL)
  
  observeEvent(input$MergeGP,{ 
       
    print("mergeIn")
    if(!is.null(dim(GenoPre()))){
     t_Geno <- reactive(GenoPre())
    }else{t_Geno <-reactive(NULL)}

    
     withProgress(message = 'Merging Pheno and Geno Data', value = 0,{

       if(!is.null(t_Geno()) & !is.null(Pheno())){
           
        mergeOut <-  tryCatch({
    		    withCallingHandlers({
      			  getMergedData(t_Geno(),Pheno(),TargetIDs())
      			 }, warning = function(w){
        			# This function will execute if a warning is issued
        			# Capture the warning message
        			warningMessage <- conditionMessage(w)
        			print(paste("Warning captured:", warningMessage))

        			# Invoke the default warning handler (optional)
        			invokeRestart("muffleWarning")
      		   }, error=function(e) {
               print(paste("Error encountered:", e$message))
            })
          })
        } else if(is.null(t_Geno()) | is.null(Pheno())){
         return(NULL)
        } 
         
         mergedData(mergeOut)
    })
  })
      
#####   
  checkPhDataME <- reactiveVal(NULL)
  
  observeEvent(input$MergeGPME,{
        
       # browser()
        print("mergeIn")
        if(!is.null(dim(GenoPre()))){
          t_Geno <- reactive(GenoPre())
        }else{t_Geno <-reactive(NULL)}
        if(!is.null(TraitME())){
         withProgress(message = 'Merging Pheno and Geno Data', value = 0,{
    
          
          if(isPhenoME()){
            
            mergeOutME <- tryCatch({
              withCallingHandlers({
                getMergedDataME(phenoMEData(),t_Geno(),TargetIDs())
              }, warning = function(w){
                # This function will execute if a warning is issued
                # Capture the warning message
                warningMessage <- conditionMessage(w)
                print(paste("Warning captured:", warningMessage))
                
                # Invoke the default warning handler (optional)
                invokeRestart("muffleWarning")
              }, error=function(e) {
                print(paste("Error encountered:", e$message))
              })
            })
            
            mergedData(mergeOutME)
            
          }      
        })
        }else if(is.null(TraitME())){ 
          checkPhDataME("Select trait and merge data")
        }
  })
  
  output$MssgTrait <- renderText({
    checkPhDataME()
  })
###
# 
# 

####

  output$ReadResult <- renderPrint({
    ReadReslt <- mergedData()
    print(paste("Class of ReadReslt:", class(ReadReslt)))  # Debugging line
    if (!is.null(ReadReslt)) {
      if (is.character(ReadReslt)) {
        return(paste("Error:", ReadReslt))
      } else {
        return(paste("Merged Data Dim:", paste(dim(ReadReslt[[1]]), collapse = " x ")))
      }
    } else {
      return("No data available.")
    }
  })

  
  
###   
 
  processedData <- reactive({
    print("prIn")
    if(!isPhenoME()){
      getProcessedData(mergedData(),Trait())
    }else if(isPhenoME()){
      mergedData()
    }
  })
  
### Plot distribution of trait values across Env  
  
  DT_Filt_List <- reactive({
    print("prIn2")
    if(isPhenoME() && class(processedData())=="list"){
      processedData()[[1]]
    }
  })
  
  genoDat_List <- reactive({ 
    if(isPhenoME() && class(processedData())=="list"){
      processedData()[[2]]
    }
  })

  
  
###  
  
  plotData <- reactiveVal(NULL)
  
 # observeEvent(input$traitME,{
  
  observe({ 
     
    req(nSelTraitsME())
    req(DT_Filt_List())
    req(TraitME())
    
     
    
    print("in")
    if(nSelTraitsME()==1 && !is.null(TraitME())){
      print("in1")
      tryCatch({
        withCallingHandlers({
          plotData(getPhenoDistbnPlots(DT_Filt_List(),TraitME(),1,outResults_Dir()))
         }, warning = function(w){
          # This function will execute if a warning is issued
          # Capture the warning message
          warningMessage <- conditionMessage(w)
          print(paste("Warning captured:", warningMessage))
          
          # Invoke the default warning handler (optional)
          invokeRestart("muffleWarning")
         }, error=function(e) {
          print(paste("Error encountered:", e$message))
        })
      })
      
    }else if(nSelTraitsME()>1 & !is.null(DT_Filt_List()) & !is.null(TraitME())){
       
         plotList <- lapply(1:nSelTraitsME(),function(loc_i){
           
         tryCatch({
             withCallingHandlers({
               getPhenoDistbnPlots(DT_Filt_List(),TraitME()[loc_i],loc_i,outResults_Dir())
             }, warning = function(w){
               # This function will execute if a warning is issued
               # Capture the warning message
               warningMessage <- conditionMessage(w)
               print(paste("Warning captured:", warningMessage))
               
               # Invoke the default warning handler (optional)
               invokeRestart("muffleWarning")
             }, error=function(e) {
               print(paste("Error encountered:", e$message))
             })
          })
        })
        plotData(plotList)
    }else if(is.null(TraitME())){ 
        plotData(NULL)
    }
  })
    
   
# Render the main plot

 output$distPlots <- renderPlot({ 
    if(nSelTraitsME()==1){
     print("plLD")
    return(plotData())
   }
 })
 
 # Define a reactive expression to generate renderPlot functions for each trait
 plotRenderList <- reactive({
   if (nSelTraitsME() > 1){
     print("plLD_Mult")
     lapply(1:nSelTraitsME(), function(i) {
       renderPlot({
         plotData()[[i]]
       })
     })
   }
  })
 
 # Render individual plots for each trait using renderPlot functions generated in plotRenderList
 
 outputListPhME <- reactiveVal(list())
 
 dummyPlt <- reactive({ 
  
   if(nSelTraitsME()>1){
     outputListPhME(lapply(1:nSelTraitsME(), function(i) {
       plotOutput(noquote(paste("'","distPlot",i,"'",sep="")))}))
   }
 })
 

# Link plotRenderList and outputList
 
 #observeEvent(input$traitME,{
   
 observe({
   req(nSelTraitsME())
   req(TraitME())
   req(dummyPlt())
   
   if(!is.null(TraitME()) &  nSelTraitsME()>1){
     if(nSelTraitsME()>1 & length(plotRenderList()) >1 & length(outputListPhME()) >1){
           tempList <-  plotRenderList()
           outputListPhME(tempList) 
     }else{print(paste("The Render plot list with length", length(plotRenderList())," or outputPlot list with length",length(outputListPhME()),
          "do not meet the requirements",sep=""))
     }
   }
 })
     
 

#### 

output$LocDistributionME <- renderUI({
   print("outLD")
   if(!is.null(TraitME())){
     if (nSelTraitsME() == 1) {
       plotOutput("distPlots")
     }else if (nSelTraitsME() > 1) { 
       plotOutList <- outputListPhME()
       do.call(tagList,plotOutList)
     }
   }
 })

# 

#### 3

temp_file3e <- reactiveVal('none')
observeEvent(input$trait,{
  # Path for the temporary file
  temp_file3e(tempfile())
  
  if(length(Trait())>0 & !is.null(Pheno())){
    sink(temp_file3e())
    cat("Summary of selected trait values : \n")
    cat(paste("Trait",paste(names(summary(Pheno()[,Trait()[1]])),collapse="\t"),sep="\t"))
    cat("\n")
    for(nT in 1:length(Trait())){
      cat(paste(Trait()[nT],paste(round(summary(Pheno()[,Trait()[nT]]),digits=2),collapse="\t"),sep="\t"))
      cat("\n")
    }
    sink()
  }
})

### Test output 
# Periodically read the file and update the UI

mssgMergedPhenoStatsOut <- reactiveVal(NULL)

output$messageMergeTrtSE <- renderText({
  
  invalidateLater(1000, session) # Update every second
  if(file.exists(temp_file3e())){
    lines <- readLines(temp_file3e(), warn = FALSE)
    mssgMergedPhenoStatsOut(paste(lines, collapse = "\n"))
    return(paste(lines, collapse = "\n"))
  }else {
    return("Waiting for output...")
  }
})

####

temp_file3g <- reactiveVal('none')
observeEvent(input$traitME,{
  # Path for the temporary file
  temp_file3g(tempfile())
  
  # Start the process in a separate R process
  
  if(length(TraitME())>0 & !is.null(PhenoME())){
    sink(temp_file3g())
    cat("Summary of selected trait values across all locations : \n")
    
    if(nSelTraitsME()==1){
      cat(paste("Trait",paste(names(summary((PhenoME()[,TraitME()]))),"\t",sep=""),sep="\t"))
      cat("\n")
      cat(paste(TraitME(),paste(round(summary(as.numeric(PhenoME()[,TraitME()])),digits = 3),"\t",sep=""),sep="\t"))
    }else if(nSelTraitsME()>1){ 
      cat(paste("Trait",paste(names(summary((PhenoME()[,TraitME()[1]]))),"\t",sep=""),sep="\t"))
      cat("\n")
      
      for(nSel in 1:nSelTraitsME()){
        
        cat(paste(TraitME()[nSel],paste(round(summary(as.numeric(PhenoME()[,TraitME()[nSel]])),digits = 3),"\t",sep=""),sep="\t"))
        cat("\n")
      }
    }
    sink()
  }
})

### Test output 
# Periodically read the file and update the UI

output$messageMergeTrtME <- renderText({
  
  invalidateLater(1000, session) # Update every second
  if(file.exists(temp_file3g())){
    lines <- readLines(temp_file3g(), warn = FALSE)
    mssgMergedPhenoStatsOut(paste(lines, collapse = "\n"))
    return(paste(lines, collapse = "\n"))
  }else {
    return("Waiting for output...")
  }
})



#### Enviromics 
 
LocCoords <- reactive({
   
   locCoordsTab <- input$inFileLocCoord
   ext <- tools::file_ext(locCoordsTab$datapath)
   req(locCoordsTab)
   validate(need(ext == "csv", "Please upload a csv file"))
   
   read.csv(locCoordsTab$datapath, header = input$header)
   
})

startDate <- reactive(input$startDate)
endDate <- reactive(input$endDate) 


temp_file4 <- reactiveVal('none')
envData <- reactiveVal(NULL)


observeEvent(input$getEnvData,{
  # Path for the temporary file
# browser()
  temp_file4(tempfile())
  if(!is.null(LocCoords())){
    sink(temp_file4())
    
     withProgress(message = 'Collecting Weather Data', value = 0, {
        envDat <- getEnvData(LocCoords(),startDate(),endDate())
     })
   
    sink()
    
    envData(envDat)
  }
})

### Test output 
# Periodically read the file and update the UI
envDatStatus <- reactive({
  invalidateLater(1000, session)
  if(file.exists(temp_file4())){
    lines <- readLines(temp_file4(), warn = FALSE)
    return(paste(lines, collapse = "\n"))
  }else {
    return("...")
  }
})

output$messageEnvDat <- renderText({
 envDatStatus()  

})

output$EnvTab <- renderDataTable({
  if(!is.null(envData())){
   envData()[1:250,]
  }
})

processWthDat <- reactive(input$processWth)



### EnvK 

 
 gaussVar <- reactive(input$GaussKE)

 
 EnvK <- reactiveVal(NULL)
 observeEvent(input$getEnvK,{
   #browser()
    withProgress(message='Estimating Environmental Relationships',value=0,{
      EnvK(getEnvKernel(envData(),processWthDat(),gaussVar())) 
      
    })
 })



 EnvK_Mod <- reactiveVal(NULL)

 observeEvent(input$OtherLoc,{
          EnvK_Mod(syncEnvPhenoDat(EnvK(),LocCoords(),OtherLoc()))
 })


OtherLoc <- reactive(input$OtherLoc)

# Render the main plot

output$envPlot <- renderPlot({ 
  print("envPl")
  if(!is.null(EnvK())){ 
   if(OtherLoc()==FALSE){
     tryCatch({
       withCallingHandlers({
         plotEnvRel(EnvK(),outResults_Dir())
       }, warning = function(w){
         # This function will execute if a warning is issued
         # Capture the warning message
         warningMessage <- conditionMessage(w)
         print(paste("Warning captured:", warningMessage))
         
         # Invoke the default warning handler (optional)
         invokeRestart("muffleWarning")
       }, error=function(e) {
         print(paste("Error encountered:", e$message))
       })
     })
   }else if(OtherLoc()==TRUE){
     tryCatch({
       withCallingHandlers({
         plotEnvRel(EnvK_Mod(),outResults_Dir())
       }, warning = function(w){
         # This function will execute if a warning is issued
         # Capture the warning message
         warningMessage <- conditionMessage(w)
         print(paste("Warning captured:", warningMessage))
         
         # Invoke the default warning handler (optional)
         invokeRestart("muffleWarning")
       }, error=function(e) {
         print(paste("Error encountered:", e$message))
       })
     })
     EnvK(EnvK_Mod())
   }
  }
})


#####  Merge Geno Pheno Data

######## 

temp_file3f <- reactiveVal('none')


# Start the process in a separate R process
# 
   

observe({
    
  if(!isPhenoME()){ 
   
   req(mergedData())
   mergeStatus <- reactive(mergedData()[[4]])
 
   #missingStatus <- reactive(mergedData()[[3]])
  
   missingStatus <- reactive({
    # Get the third element from mergedData() which is a list with two elements
    data <- mergedData()[[3]]
    return(data)
   })
 
####
  
  temp_file3f(tempfile())

    sink(temp_file3f())

    if(mergeStatus()==TRUE){
        cat(paste("Data merge succeeded.\n"))
        cat("\n")
        #cat(paste("The merged data table has data for ",nrow(processedData()[[1]])," genotypes and ", ,ncol(processedData()[[1]])," and markers \n",sep=""))
        cat(paste("The merged data table has training data for ",nrow(mergedData()[[1]])," genotypes and ",ncol(mergedData()[[1]])-(nTraits()+2)," markers \n",sep=""))
  	    cat("\n")
	      if(!is.null(mergedData()[[2]])){
	       cat(paste("The merged data table has target data for ",nrow(mergedData()[[2]])," genotypes and ",ncol(mergedData()[[2]])," markers \n",sep=""))
	       cat("\n")

	      }else{
	       cat(paste("Target data is not defined. Please upload target data or target line IDs.\n"))
	        cat("\n")
	      }

  	
      cat(paste("The following IDs were missing in geno table, but were present in pheno table: ",paste(missingStatus()[[1]],collapse=","),"\n"))
      cat(paste("The following IDs were missing in pheno table, but were present in geno table: ",paste(missingStatus()[[2]],collapse=","),"\n"))
    }else if(mergeStatus()==FALSE){
      cat(paste("Data merge failed.Examine the input files for id match. \n"))
    }

    sink()
    
  }

})
  
  
### Test output 
## Periodically read the file and update the UI

mssgMergedGenoStatsOut <- reactiveVal(NULL)

output$messageMergedGeno <- renderText({

  invalidateLater(1000, session) # Update every second
  if(file.exists(temp_file3f())){
    lines <- readLines(temp_file3f(), warn = FALSE)
    mssgMergedGenoStatsOut(paste(lines, collapse = "\n"))
    return(paste(lines, collapse = "\n"))
  }else {
    return("Waiting for output...")
  }
})

#####


temp_file3h <- reactiveVal('none')


# Start the process in a separate R process
# 
observe({
  if(isPhenoME()){
   req(mergedData())
   mergeStatus <- reactive(mergedData()[[4]])
   missingStatus <- reactive(mergedData()[[3]])
  
   temp_file3h(tempfile())
   sink(temp_file3h())
  
   if(mergeStatus()==TRUE){
      cat(paste("Data merge succeeded.\n"))
      cat("\n")
      #cat(paste("The merged data table has data for ",nrow(processedData()[[1]])," genotypes and ", ,ncol(processedData()[[1]])," and markers \n",sep=""))
      #cat(paste("The merged data table has training data for ",nrow(mergedData()[[1]])," genotypes and ",ncol(mergedData()[[1]])," markers \n",sep=""))
      #cat("\n")
      if(!is.null(mergedData()[[2]])){
        cat(paste("The merged data table has target data for ",nrow(mergedData()[[2]][[1]])," genotypes and ",ncol(mergedData()[[2]][[1]])," markers \n",sep=""))
        cat("\n")
        
      }else{cat(paste("Target data is not defined. Please upload target data or target line IDs.\n"))}
      
      cat(paste("The following IDs were missing in geno table, but were present in pheno table: ",paste(unlist(missingStatus()[[1]]),collapse=", "),"\n",sep=""))
      cat(paste("The following IDs were missing in pheno table, but were present in geno table: ",paste(unlist(missingStatus()[[2]]),collapse=", "),"\n"),sep="")
    }else if(mergeStatus()==FALSE){
      cat(paste("Data merge failed.Examine the input files for id match. \n"))
    }
  
    sink()
  }
})


output$messageMergedGenoME <- renderText({
  invalidateLater(1000, session) # Update every second
  if(file.exists(temp_file3h())){
    lines <- readLines(temp_file3h(), warn = FALSE)
    mssgMergedGenoStatsOut(paste(lines, collapse = "\n"))
    return(paste(lines, collapse = "\n"))
  }else {
    return("Waiting for output...")
  }
})


#########


  
 observeEvent(input$MergeGP, {
    updateNumericInput(inputId ="noCandidates",value= nrow(processedData()[[1]]),min =2, max=nrow(processedData()[[1]]))}) 
  
 observeEvent(input$MergeGP,{
    updateNumericInput(inputId ="noToSelect",value= 100,min =2, max=nrow(processedData()[[1]]))})
  
  
# TS Optimization
  
  noCandidates <- reactive(input$noCandidates)
  nTrainToSelect <- reactive(input$noToSelect)
  optimCriteria <- reactive(input$optCriteria)
  
  predictionData <- reactive({getPredictionData(processedData(),noCandidates())}) 
  
  GAParameters <- reactiveValues(npop=100,nelite=10,mutprob=0.5,mutintensity=1,niterations=100,minitbefstop=50,tabu=TRUE,tabumemsize=1,plotiters=FALSE,errorstat="PEVMEAN2",lambda=1e-6,mc.cores=1)
  
  
  
  checkSEData <- reactiveVal(NULL) 
  
  observe({ 
    
    if(!is.null(checkSTData()) && !is.null(checkMTData())){
      checkSEData(checkSTData() || checkMTData())
     }
    
  })
  
 
  temp_file4b <- reactiveVal(NULL)
  
  rProcess4b <- eventReactive(input$Optimize,{
    
   # browser()
    temp_file4b(tempfile())
    
    
    if(checkSEData()==TRUE){
      
    # Start the process in a separate R process
    b <- callr::r_bg(function(getOptimalTS,getRandomTS,getTSComparisons,getTSComparisonsMT,predictionData,Trait,nTraits,noCandidates,nTrainToSelect,GAParameters,nSelTraits,TargetIDs,temp_file4b_path){
        
        library(bWGR)
        library(NAM)
        library(STPGA)
        library(BGLR)
        
        sink(temp_file4b_path)
        
        cat(paste("Running optimizations for selection of training sets"))    
        
        #GAParams <- isolate(reactiveValuesToList(GAParameters))
        ## Removed isolate to avoid path.expand invalid error
        GAParams <- (GAParameters)
       
        Train_STPGA <- tryCatch({
          withCallingHandlers({
           getOptimalTS(predictionData,unlist(Trait),nTraits,noCandidates,nTrainToSelect,GAParams)
          }, warning = function(w){
            # This function will execute if a warning is issued
            # Capture the warning message
            warningMessage <- conditionMessage(w)
            print(paste("Warning captured:", warningMessage))
            # Invoke the default warning handler (optional)
            invokeRestart("muffleWarning")
          }, error=function(e) {
            print(paste("Error encountered:", e$message))
          })
        })
        
        
        cat(paste("Running random selection of training sets"))
        
        Train_Random <- tryCatch({
          withCallingHandlers({
            getRandomTS(predictionData,unlist(Trait),nTraits,noCandidates,nTrainToSelect)
            
          }, warning = function(w){
            # This function will execute if a warning is issued
            # Capture the warning message
            warningMessage <- conditionMessage(w)
            print(paste("Warning captured:", warningMessage))
            
            # Invoke the default warning handler (optional)
            invokeRestart("muffleWarning")
          }, error=function(e) {
            print(paste("Error encountered:", e$message))
          })
        })
          
       ####   
        cat(paste("Running crossvalidations for TS set comparisons"))
        if (!is.null(Train_STPGA) && !is.null(Train_Random) &&
            length(Train_STPGA) > 0 && length(Train_Random) > 0) {

          TSOptOutputList <- if (nSelTraits == 1) {
            getTSComparisons(
              predictionData, Train_STPGA, Train_Random,
              unlist(Trait), nTraits, TargetIDs
            )
          } else if (nSelTraits > 1){
            getTSComparisonsMT(
              predictionData, Train_STPGA, Train_Random,
              unlist(Trait), nTraits, TargetIDs
            )
          }

         # Optional debugging or status confirmation
          cat("TSOptOutputList updated successfully.\n")
         
         # Return the final TSOptOutputList
          return(TSOptOutputList)
          cat("\n")
        }else{return(NULL)}
        
       sink()
        
      }, args = list(
        getOptimalTS = getOptimalTS,
        getRandomTS = getRandomTS,
        getTSComparisons = getTSComparisons,
        getTSComparisonsMT=getTSComparisonsMT,
        predictionData = predictionData(),
        Trait=unlist(Trait()),
        nTraits= nTraits(),
        noCandidates= noCandidates(),
        nTrainToSelect = nTrainToSelect(),
        GAParameters = reactiveValuesToList(GAParameters),
        nSelTraits=nSelTraits(),
        TargetIDs=TargetIDs(),
        temp_file4b_path = temp_file4b()
       ), stdout = temp_file4b(), stderr = temp_file4b())
 
    }else if(checkSEData()==FALSE){
      
      sink(temp_file4b())      
      cat("Current data is not suitable for optimization of training sets.\n")
      sink()
      
    }
    
  })
  
  
  SEOptTSMssg <- reactiveVal(NULL)
  TSOptOutputTab <- reactiveVal(NULL)
  processComplete4b <- reactiveVal(FALSE)
  
  # Reactive expression to check process status and read temp file
  processStatus4b <- reactive({
    invalidateLater(1000, session) # Update every second
    
    if (!is.null(rProcess4b()) && rProcess4b()$is_alive()) {
      if (file.exists(temp_file4b())) {
        lines <- readLines(temp_file4b(), warn = FALSE)
        return(paste(lines, collapse = "\n"))
      } else {
        return("Waiting for output...")
      }
    } else if (!is.null(rProcess4b()) && !rProcess4b()$is_alive()) {
      if (file.exists(temp_file4b())) {
        lines <- readLines(temp_file4b(), warn = FALSE)
        txt <- paste(lines, collapse = "\n")
        processComplete4b(TRUE)
        TSOptOutputTab(rProcess4b()$get_result()) 
        SEOptTSMssg(paste(txt, "\n Computations completed", collapse = "\n"))
        return(paste(txt, "\n Computations completed", collapse = "\n"))
      } else {
        return("Waiting for output...")
      }
    } else {
      return("Error")
    }
  })
 
  
  output$messageOptTS <- renderText({ 
    processStatus4b()
  })
  
  
  tsOptHead <- eventReactive(input$Optimize,{paste("Cross Validation Accuracies")})
  output$tsOptHeader <- renderText({tsOptHead()}) 
  
  output$PredAccuracyTable <- renderTable({ 
    if(!is.null(TSOptOutputTab())){
      TSTable <- TSOptOutputTab()
      colnames(TSTable) <- c("Ridge Regression","BayesB","Bayes LASSO")
      rownames(TSTable) <- c("Optimal Set","Random Set")
      TSTable
    }
  }, colnames=TRUE,rownames=TRUE,digits=2)
  
  
  Mssg1optTS <- reactiveVal(NULL)
  Mssg1optTSOK <- reactiveVal(NULL)
  
  #observeEvent(input$MergeGP,{
    
    observe({ 
    
    okMssgTS <- "Data OK for optimization of training set.\n"
    MssgTS <- "Data not suitable for optimizaiton of training set using the current method"
    
    if(!is.null(checkSEData())){
      if(checkSEData()){ 
        Mssg1optTSOK(okMssgTS)
        Mssg1optTS(NULL)
        
      }else{
        Mssg1optTSOK(NULL)
        Mssg1optTS(MssgTS)
        
      }
    }
  })

########  
  
  output$MssgOptTS <- renderText({
    Mssg1optTS()
  })
  
  output$MssgOptTSOK <- renderText({
    Mssg1optTSOK()
  })
  
######  MssgOptTS
  
## CVR
  
  k <- reactive(input$k)
  nIter <- reactive(input$nIter)
  
### Conflicting statement messages for ST CV and GP pages 
  
  
  Mssg1CVR <- reactiveVal(NULL)
  Mssg1GP <- reactiveVal(NULL)
  
  Mssg1CVROK <- reactiveVal(NULL)
  Mssg1GPOK <- reactiveVal(NULL)
  
  Mssg1MergeOK <- reactiveVal(NULL)
  
  checkSTData <- reactiveVal(NULL)
  
  observeEvent(input$MergeGP,{
  
    okMTraitMssg <- "Data OK for cross validation of single trait GP models for multiple traits.\n"
    okMssg <- "Data OK for single trait CV"
   
    okMTraitMssgGP <- "Data OK for single trait genomic prediction of multiple traits.\n"
    okMssgGP <- "Data OK for single trait GP"
    
  
    # Update the reactive values based on conditions
    if(input$MergeGP){
      if(nSelTraits()==1 && !isPhenoME()){
        Mssg1CVROK(okMssg)
        Mssg1GPOK(okMssgGP)
        Mssg1CVR(NULL)
        Mssg1GP(NULL)
        checkSTData(TRUE)
        Mssg1MergeOK("")
      }else if(nSelTraits() > 1 && !isPhenoME()){
        Mssg1CVROK(okMTraitMssg)
        Mssg1GPOK(okMTraitMssgGP)
        Mssg1CVR(NULL)
        Mssg1GP(NULL)
        checkSTData(TRUE)
        Mssg1MergeOK("")
      }else if(nSelTraits() == 0 && !isPhenoME()){
        Mssg1MergeOK("You have not selected any trait. Please select a trait and merge.")
      }
    }
   
    print(Mssg1CVR())
  })
 
  
  
  observeEvent(input$MergeGPME,{
    
    meTraitMssg <- "Warning! You've loaded pheno data from multiple environments.\n Try CV for multi-environmental models. \n"
    meTraitMssgGP <- "Warning! You've loaded pheno data from multiple environments.\n Try genomic prediction for multi-environmental data. \n"
    
    
    if(nSelTraitsME()>=1 && isPhenoME()){
      Mssg1CVROK(NULL)
      Mssg1GPOK(NULL)
      Mssg1CVR(meTraitMssg)
      Mssg1GP(meTraitMssgGP)
      checkSTData(FALSE)
      Mssg1MergeOK("")
    }else if(nSelTraitsME()==0 && isPhenoME()){
      Mssg1MergeOK("You have not selected any trait. Please select a trait and merge.")
    }
    
  })
  
### Render messages 
   
  
  output$MssgST <- renderText({
      Mssg1CVR()
  })
  
  output$MssgSTGP <- renderText({
    Mssg1GP()
  })
  
  output$MssgSTOK <- renderText({
    Mssg1CVROK()
  })
  
  output$MssgSTGPOK <- renderText({
    Mssg1GPOK()
  })
  
  output$MssgMergedData <- renderText({ 
    Mssg1MergeOK()
  })
  
 ##  Messages for MT CV and GP pages 
  
   
   Mssg2CVR <- reactiveVal(NULL)
   Mssg2GP <- reactiveVal(NULL)
   Mssg2CVROK <- reactiveVal(NULL)
   Mssg2GPOK <- reactiveVal(NULL)
   
   checkMTData <- reactiveVal(NULL)
   
  observeEvent(input$MergeGP,{
  
     # Define the messages
     
     mTraitMssg <- "Warning! Select more than one trait for cross validation of multi-trait models.\n"
     okMssg <- "Data OK for Multi-trait CV"
     

     mTraitMssgGP <- "Warning! Select more than one trait for genomic prediction of multi-trait GP models.\n"
     okMssgGP <- "Data OK for Multi-trait GP"
     
     if(nSelTraits()==1 && !isPhenoME()){
       
       Mssg2CVR(mTraitMssg)
       Mssg2GP(mTraitMssgGP)
       Mssg2CVROK(NULL)
       Mssg2GPOK(NULL)
       
       checkMTData(FALSE)
     }else if(nSelTraits() > 1 && !isPhenoME()){
       
       Mssg2CVR(NULL)
       Mssg2GP(NULL)
       Mssg2CVROK(okMssg)
       Mssg2GPOK(okMssgGP)
       
       checkMTData(TRUE)
     }
     
     print(Mssg2CVR())
   })
   
    
  observeEvent(input$MergeGPME,{ 
    mTraitMssg <- "Warning! Select more than one trait for cross validation of multi-trait models.\n"
    mTraitMssgGP <- "Warning! Select more than one trait for genomic prediction of multi-trait GP models.\n"
  
    meTraitMssg <- "Warning! You've loaded pheno data from multiple environments. Try CV for multi-environmental models. \n"
    meTraitMssgGP <- "Warning! You've loaded pheno data from multiple environments. Try GP for multi-environmental models. \n"
    
   
    # Update the reactive values based on conditions
      if(nSelTraitsME()==1 && isPhenoME()){
        
        Mssg2CVR(paste(mTraitMssg, meTraitMssg, sep = " "))
        Mssg2GP(paste(mTraitMssgGP, meTraitMssgGP, sep = " "))
        Mssg2CVROK(NULL)
        Mssg2GPOK(NULL)
        
        checkMTData(FALSE)
      }else if(nSelTraitsME()>1 && isPhenoME()){
        
        Mssg2CVR(meTraitMssg)
        Mssg2GP(meTraitMssgGP)
        Mssg2CVROK(NULL)
        Mssg2GPOK(NULL)
        
        checkMTData(FALSE)
      }
  })
  
## Render messages to MT pages 
  
  output$MssgMT <- renderText({
      Mssg2CVR()
  })
 
   output$MssgMTGP <- renderText({
     Mssg2GP()
  })
   
   output$MssgMTOK <- renderText({
   
     Mssg2CVROK()
   })
   
   output$MssgMTGPOK <- renderText({
     
     Mssg2GPOK()
   })
  
#####
###
   
   temp_file5a <- reactiveVal('none')
   processComplete5a <- reactiveVal(FALSE)
   cvrOutputListST <- reactiveVal(NULL)
   
   # Start the background process and return rProcess
   
   rProcess5a <- eventReactive(input$CrossValidationST,{
     
     # Path for the temporary file
     temp_file5a(tempfile())
    
     if(checkSTData()==TRUE){
       # Start the process in a separate R process
       callr::r_bg(function(getemCVR, predictionData, Trait, nTraits, k, nIter,nSelTraits,temp_file5a_path){
         library(bWGR)
         library(rrBLUP)
         library(NAM)
         sink(temp_file5a_path)
         
         if(nSelTraits==1){
           cat(paste("Running Cross Validation for ",unlist(Trait),"\n",sep=""))    
           PATab <- getemCVR(predictionData,unlist(Trait),nTraits,k,nIter)
           PATable <- round(PATab[c("emRR","emBB","emBL")],digits=2) 
           names(PATable)<- c("Ridge Regression","BayesB","Bayes LASSO")
           PATable <- rbind.data.frame(names(PATable),PATable)
           colnames(PATable) <- c("Ridge Regression","BayesB","Bayes LASSO")
           rownames(PATable) <- c("Model","Accuracy")
           return(PATable[2,])
         }
         
         if(nSelTraits >1){
           PATableComb <- c()
           for(nSelT in 1:nSelTraits){
             cat(paste("Running Cross Validation for ",unlist(Trait)[nSelT],"\n",sep=""))    
             PATab <- getemCVR(predictionData,unlist(Trait)[nSelT],nTraits,k,nIter)
             PATable <- round(PATab[c("emRR","emBB","emBL")],digits=2)
             PATableComb <- rbind.data.frame(PATableComb,PATable)
             cat("\n")
           }
           colnames(PATableComb) <- c("Ridge Regression","BayesB","Bayes LASSO")
           rownames(PATableComb) <- unlist(Trait)
           return(PATableComb)
         }
         
          sink()
         
       }, args = list(
         getemCVR = getemCVR,
         predictionData = predictionData(),
         Trait = unlist(Trait()),
         nTraits = nTraits(),
         k = k(),
         nIter = nIter(),
         nSelTraits=nSelTraits(),
         temp_file5a_path = temp_file5a()
       ), stdout = temp_file5a(), stderr = temp_file5a())
       
     }else if(checkSTData()==FALSE){
       
       sink(temp_file5a())
       cat("Current data is not suitable for single trait in single environment CV.\n")
       sink()
       
     }
     
   })
   
   STCVMssg <- reactiveVal(NULL)
   
   # Reactive expression to check process status and read temp file
   processStatus5a <- reactive({
     invalidateLater(1000, session) # Update every second
     
     if (!is.null(rProcess5a()) && rProcess5a()$is_alive()) {
       if (file.exists(temp_file5a())) {
         lines <- readLines(temp_file5a(), warn = FALSE)
         return(paste(lines, collapse = "\n"))
       } else {
         return("Waiting for output...")
       }
     } else if (!is.null(rProcess5a()) && !rProcess5a()$is_alive()) {
       if (file.exists(temp_file5a())) {
         lines <- readLines(temp_file5a(), warn = FALSE)
         txt <- paste(lines, collapse = "\n")
         
         processComplete5a(TRUE)
         cvrOutputListST(rProcess5a()$get_result()) /
           STCVMssg(paste(txt, "\n Cross Validation Completed", collapse = "\n"))  
           return(paste(txt, "\n Cross Validation Completed", collapse = "\n"))
         
       }else {
         return("Waiting for output...")
       }
     }else {
       return("Error")
     }
   })
   
   # Update the UI with process status
   output$messageST5 <- renderText({
     processStatus5a()
   })
   
   
######
  
  temp_file5b <- reactiveVal('none')
  processComplete5b <- reactiveVal(FALSE)
  cvrOutputListMT <- reactiveVal(NULL)
  
  # Start the background process and return rProcess
  rProcess5b <- eventReactive(input$CrossValidationMT, {
    
    # Path for the temporary file
    temp_file5b(tempfile())
  
    if(checkMTData()==TRUE){
    # Start the process in a separate R process
    callr::r_bg(function(getMTCVR, predictionData, Trait, nTraits, k, nIter, temp_file5b_path){
      library(BGLR)
      library(rrBLUP)
      library(NAM)
     
      sink(temp_file5b_path)
      result <- getMTCVR(predictionData, Trait, nTraits, k, nIter)
      sink()
      result
    }, args = list(
      getMTCVR = getMTCVR,
      predictionData = predictionData(),
      Trait = unlist(Trait()),
      nTraits = nTraits(),
      k = k(),
      nIter = nIter(),
      temp_file5b_path = temp_file5b()
    ), stdout = temp_file5b(), stderr = temp_file5b())
    
    }else if(checkMTData()==FALSE){
      
      sink(temp_file5b())
      cat("Current data is not suitable for multi-trait CV.\n")
      sink()
      
    }
    
 })
  
  MTCVMssg <- reactiveVal(NULL)
  # Reactive expression to check process status and read temp file
  processStatus5b <- reactive({
    invalidateLater(1000, session) # Update every second
    
    if (!is.null(rProcess5b()) && rProcess5b()$is_alive()) {
      if (file.exists(temp_file5b())) {
        lines <- readLines(temp_file5b(), warn = FALSE)
        return(paste(lines, collapse = "\n"))
      } else {
        return("Waiting for output...")
      }
    } else if (!is.null(rProcess5b()) && !rProcess5b()$is_alive()) {
      if (file.exists(temp_file5b())) {
        lines <- readLines(temp_file5b(), warn = FALSE)
        txt <- paste(lines, collapse = "\n")
        processComplete5b(TRUE)
        cvrOutputListMT(rProcess5b()$get_result()) /
        MTCVMssg(paste(txt,"\n","Cross Validation Completed", collapse = "\n"))
        return(paste(txt,"\n","Cross Validation Completed", collapse = "\n"))
      } else {
        return("Waiting for output...")
      }
    } else {
      return("Error")
    }
  })
  
  # Update the UI with process status
  output$messageMT5 <- renderText({
    processStatus5b()
  })

###
  
 

##### CVR ME Version 
  
  
  observeEvent(input$infileBLUEsME, {
    updateSelectInput(inputId ="YearMECV",choices = c("All",YearsME()) )
  })
  
  observeEvent(input$infileBLUEsME, {
    updateSelectInput(inputId ="LocationMECV",choices = c("All",LocationsME()))
    
  })
  
  
  observeEvent(input$infileBLUEsME, {
    updateSelectInput(inputId ="fixedMECV",choices = IDColsME())
  })
  
  observeEvent(input$infileBLUEsME, {
    updateSelectInput(inputId ="EnvVarIDCV",choices = IDColsME())
  })
  
  kCV <- reactive(input$kME)
  nIterCV <- reactive(input$nIterME)

  LocationMECV <- reactiveVal(NULL)
  YearMECV <- reactiveVal(NULL)
  fixMECV <- reactiveVal(NULL)
  varEnvCV <- reactiveVal(NULL)
  
###  
  processComplete6 <- reactiveVal(FALSE)
  MECV_Out <- reactiveVal(NULL)
  
  MECV_Out_SwitchTab <- reactiveVal(NULL)
  
###  
  
  CV_Switch <- reactiveVal(0)
  CV_run_Tab <- reactiveVal(0)
  CV_out_Tab <- reactiveVal(0)
   
  temp_file6 <- reactiveVal('none')
  CVMet <- reactive(input$CVMet)
  factVar <- reactive(input$CVFactor)
 
### 
   
  observeEvent(input$CrossValidationME,{
    
    processComplete6(FALSE)
    MECV_Out(NULL)
    temp_file6('none')
   
    if(CVMet()!= "CV_LOFO"){
      LocationMECV("All")
      YearMECV("All")
      fixMECV("Loc")
      varEnvCV("Loc")
    }else if(CVMet() == "CV_LOFO"){
      YearMECV(input$YearMECV)
      LocationMECV(input$LocationMECV)
      fixMECV(input$fixedMECV)
      varEnvCV(input$EnvVarIDCV)
    }
    
     CV_run_Tab(CV_run_Tab()+1)
  })
  
 
#### 
 observeEvent(input$CVMet,{
#   
    if(!is.null(MECV_Out()) & CV_Switch()>0){
        MECV_Out_SwitchTab(MECV_Out())
        MECV_Out(NULL)
        CV_Switch(CV_Switch()+1)
      }else if(is.null(MECV_Out()) & CV_Switch()>0){
        MECV_Out_SwitchTab(MECV_Out())
      }
     
 })
  
#### 
 
  observeEvent(input$infileBLUEsME, {
    updateSelectInput(inputId ="CVFactor",choices = IDColsME())
  })
  
### 
  
  
  checkMEData <- reactiveVal(FALSE)
  Mssg3CVR <- reactiveVal(NULL)
  Mssg3GP <- reactiveVal(NULL)
  
  Mssg3CVROK <- reactiveVal(NULL)
  Mssg3GPOK <- reactiveVal(NULL)
  
  observeEvent(input$MergeGP,{
    
    
    noMETraitMssg <- "Warning! You've not loaded pheno data from multiple environments. Try CV for single environmental models. \n "
    
    noMETraitMssgGP <- "Warning! You've not loaded pheno data from multiple environments. Try GP for single environmental models. \n "
 
    
    
    if(!isPhenoME() && input$MergeGP){
      Mssg3CVROK(NULL)
      Mssg3GPOK(NULL)
      
      Mssg3CVR(noMETraitMssg)
      Mssg3GP(noMETraitMssgGP)
      
      checkMEData(FALSE)
    }
    print(Mssg3CVR())
    
  })
  
  observeEvent(input$MergeGPME,{
  
    okMssgME <- "Data OK for Multi-environmental CV"
    okMssgMEGP <- "Data OK for Multi-environmental GP"
  
    if(isPhenoME()){
      checkMEData(TRUE)
      
      Mssg3CVROK(okMssgME)
      Mssg3GPOK(okMssgMEGP) 
      
      Mssg3CVR(NULL)
      Mssg3GP(NULL)
    }
  
  })
  
  
  
  output$MssgMEOK <- renderText({
    Mssg3CVROK()
  })
  
  output$MssgME <- renderText({
    Mssg3CVR()
  })
  
  output$MssgMEGP <- renderText({
    Mssg3GP()
  })
  
  output$MssgMEGPOK <- renderText({
    Mssg3GPOK()
  })
  
 
# Start the background process and return rProcess
 
 rProcess6 <- eventReactive(input$CrossValidationME,{
   
   
   #browser()
   
    # Path for the temporary file
    temp_file6(tempfile())
    
    #ME_argList <- list(DT_1_Filt_List,genoDat_List,traits,KG,KE,CVMet,factVar,KMethod,FitEnvModels,fixedME,envVar,IDColsME,LocME,YrME)
    
    if(checkMEData()==TRUE){
    # Start the process in a separate R process ## fitMEModels_LOF_CV,fitMEModels_CV
    callr::r_bg(function(getME_CV,DT_1_Filt_List,genoDat_List,traits,KG,KE,CVMet,factVar,k,niter,KMethod,FitEnvModels,fixedME,envVar,IDColsME,IDColME,LocME,YrME,temp_file6_path,FN) {
      library(BGLR)
      library(BGGE)
      library(EnvRtype)
      library(foreach)
      library(doParallel)
      registerDoParallel(5)
      source(FN)
      sink(temp_file6_path)
      result <- getME_CV(DT_1_Filt_List,genoDat_List,traits,KG,KE,CVMet,factVar,k,niter,KMethod,FitEnvModels,fixedME,envVar,IDColsME,IDColME,LocME,YrME)
      sink()
      result
    }, args = list(
      getME_CV = getME_CV,
      DT_1_Filt_List = DT_Filt_List(),
      genoDat_List = genoDat_List(),
      traits = TraitME(),
      KG = NULL,
      KE= NULL, 
      CVMet= CVMet(),
      factVar= factVar(),  
      k= kCV(),
      niter=nIterCV(),
      KMethod= "GK",
      FitEnvModels= fitEnvCovs(),
      fixedME= fixMECV(),
      envVar= varEnvCV(),
      IDColsME= IDColsME(),
      IDColME= IDColME(),
      LocME= LocationMECV(),
      YrME=YearMECV(),
      temp_file6_path = temp_file6(),
      FN=FN
    ), stdout = temp_file6(), stderr = temp_file6())
      
      
      # result <- getME_CV(DT_Filt_List(),genoDat_List(),TraitME(),NULL,NULL,CVMet(),factVar(),kCV(),nIterCV(),"GK",fitEnvCovs(),fixMECV(),varEnvCV(),IDColsME(),IDColME(),LocationMECV(),YearMECV())
      
      
    } else if(checkMEData()==FALSE){ 
      
       sink(temp_file6())
       cat("Current data is not suitable for multi-environmental CV")
       sink()
      
    }
  })

#### 
 
 
 
 MECVMssg <- reactiveVal(NULL)
  
# Reactive expression to check process status and read temp file

  processStatus6 <- reactive({
    invalidateLater(1000, session) # Update every second
    
    if (!is.null(rProcess6()) && rProcess6()$is_alive()) {
      if (file.exists(temp_file6())) {
        lines <- readLines(temp_file6(), warn = FALSE)
        return(paste(lines, collapse = "\n"))
      } else {
        return("Waiting for output...")
      }
    }else if (!is.null(rProcess6()) && !rProcess6()$is_alive()) {
      if (file.exists(temp_file6())) {
        lines <- readLines(temp_file6(), warn = FALSE)
        txt <- paste(lines, collapse = "\n")
        processComplete6(TRUE)
        MECV_Out(rProcess6()$get_result())
        # CV_out_Tab(CV_out_Tab()+1)
        MECVMssg(paste(txt, "Cross Validation Completed", collapse = "\n"))
        return(paste(txt, "Cross Validation Completed", collapse = "\n"))
      } else {
        return("Waiting for output...")
      }
    } else {
      return("Error")
    }
  })
  
  # Update the UI with process status
  output$messageME6 <- renderText({
    processStatus6()
  })
  
  

  MECV_Out_List_Init <- reactive({ 
    a <- list() 
    nCVMet <- length(c("CV1","CV2","CV0","CV00","CV_LOFO"))
    for(i in 1:nCVMet){a[[i]] <- 0}
    names(a) <- c("CV1","CV2","CV0","CV00","CV_LOFO")
    a
  })
  
  MECV_Out_List <- reactiveVal()
  
  # Initialize the reactive value within an observe or observeEvent
  observeEvent(input$infileBLUEsME,{
    MECV_Out_List(MECV_Out_List_Init())
  })
  
######
  
  output$emCVRST <- renderTable({
     cvrOutputListST()
  },colnames=TRUE,rownames=TRUE)
  
  
  observe({ 
    req(cvrOutputListST())
    #browser()
    if(!is.null(cvrOutputListST())){
      writeCVOutTable(cvrOutputListST(),"ST",outResults_Dir())
      
    }
  })
  
  observe({ 
    req(cvrOutputListMT())
    if(!is.null(cvrOutputListMT())){
      writeCVOutTable(cvrOutputListMT(),"MT",outResults_Dir())
    
    }
  })
  
 
  
  output$emCVRMT <- renderTable({
    if(nSelTraits()>1){
      PATable <- cvrOutputListMT()
      #colnames(PATable) <- rep("",ncol(PATable))
      print.data.frame(as.data.frame(PATable))
    }
   },colnames=TRUE,rownames=TRUE) 
  
  
  
  observe({ 
    req(gpOutputListST())
    
    if(!is.null(gpOutputListST())){
      writeGPOutTable(gpOutputListST(),"ST",outResults_Dir())
      
    }
  })
    

  observe({ 
    req(gpOutputListMT())
    
    if(!is.null(gpOutputListMT())){
      writeGPOutTable(gpOutputListMT(),"MT",outResults_Dir())
      
    }
  }) 
  
  
  observe({ 
    req(outputListMETab())
    
    if(!is.null(outputListMETab())){
      
      oFNME <-  paste(outResults_Dir(),"/MultiEnv_GP_Ouput.txt",sep="") 
      write.table(as.data.frame(outputListMETab()),oFNME, row.names = FALSE,sep="\t")
      
    }
  }) 
  

## & CV_run_Tab()==CV_out_Tab() & CV_run_Tab()==CV_out_Tab()
  
  
  cvrOutputListME <- reactive({
    if(!is.null(MECV_Out())){
      getOutTab_ME_CV(MECV_Out(),CVMet(),TraitME())
    }else{NULL}
  })
  
  
 
dummyOut <- reactive({ 
   if(!is.null(cvrOutputListME())){
      b <- MECV_Out_List()
      cvmetVar <- CVMet()
      cvmetVarInd <- which(names(b) %in% cvmetVar)
      b[[cvmetVarInd]] <-  cvrOutputListME()
      MECV_Out_List(b) 
      write.csv(as.data.frame(cvrOutputListME()),paste(outResults_Dir(),"/","ME_",cvmetVar,"_OutTable.csv",sep=""))
   }
})
 
observe({
  MECV_Out_List()
  dummyOut()
})

####

output$emCVRME <- renderTable({
    
    current_cvmet <- CVMet()
    mec_list <- MECV_Out_List()
    print(paste("MECList",mec_list))
    print(current_cvmet)
    current_cvmet_ind <- which(c("CV1","CV2","CV0","CV00","CV_LOFO") %in% current_cvmet)
    
    if(isPhenoME() & !is.null(mec_list[[current_cvmet_ind]]) & length(mec_list[[current_cvmet_ind]]) > 1){
      print("CVMEcheckin")
      PATable <- mec_list[[current_cvmet_ind]]
      print.data.frame(as.data.frame(PATable))
      
    }
   
  },colnames=TRUE,rownames=TRUE)
  
###  
  
  observeEvent(input$CrossValidationST, {
    updateTextInput(session,"CVEventST", value = "ST")
  })
  
  observeEvent(input$CrossValidationMT,{
    updateTextInput(session,"CVEventMT", value = "MT")
  })
  
  observeEvent(input$CrossValidationME,{
    updateTextInput(session,"CVEventME", value = "ME")
  })
  

  
#### Optimized Training Sets
  
  TS <- reactive(input$TrainSet) 

  optTS <-  reactive({ 
    if(TS() == "Complete Input Genotype Set"){
      optimTS <- NULL
    }else if(TS() == "Optimal Train Set from Step 3"){
      
      optimTS <- Train_STPGA()[[2]]
    }else if(TS() == "Random Train Set from Step 3"){
    
        optimTS <- Train_Random()[[1]]
    }
    (optimTS)
  })
  
  ## GP  Models   
  
  ### /// & Fixed()!= "NULL" ///| Fixed()== "NULL"
  
  GPModelST <- reactive(input$GPModelST)
  Fixed <- reactive(input$fixed)
  
  fixedData_List <- reactive({ 
    
    if(!is.null(Fixed()) & Fixed()!="NULL"){
      getFixedData_List(predictionData(),Trait(),Fixed(),TargetIDs())
    }else if(is.null(Fixed()) | Fixed() == "NULL"){ 
      "NULL" 
    } 
  })
  
  
  # 
  # outputList <- eventReactive(input$RunPredictionsST,{
  #   if(nSelTraits()==1){
  #     #browser()
  #     outputDF <-  withProgress(message = 'Running Computations', value = 0, {
  #       getRankedPredictedValues(predictionData(),nTraits(),unlist(Trait()),GPModelST(),fixedX=Fixed(),fixedData=fixedData_List(),optTS())})
  #     return(outputDF)
  #   }else if(nSelTraits()>1){
  #     outputDFComb <- c()
  #     for(nSelT in 1:nSelTraits()){
  #       i <- reactive(nSelT)
  #       outputDF <-  withProgress(message = 'Running Computations', value = 0, {
  #       getRankedPredictedValues(predictionData(),nTraits(),unlist(Trait())[i()],GPModelST(),fixedX=Fixed(),fixedData=fixedData_List(),optTS())})
  #       outputDFComb <- cbind(outputDFComb,outputDF[,2])
  #     } 
  #     outputDFComb <- cbind(outputDF[,1],outputDFComb,outputDF[,3])
  #     colnames(outputDFComb) <- c("GermplasmID",unlist(Trait()),"UpperBound of Reliability")
  #     return(outputDFComb)
  #   }
  #   
  # })
  
  # 
  # 
  # 
  # outputList <- eventReactive(input$RunPredictionsST, {
  #   # Case for single selected trait
  #   if (nSelTraits() == 1) {
  #     
  #     outputDF <- withProgress(message = 'Running Computations', value = 0, {
  #       getRankedPredictedValues(predictionData(), nTraits(), unlist(Trait()), GPModelST(), fixedX = Fixed(), fixedData = fixedData_List(), optTS())
  #     
  #     })
  #     return(outputDF)
  #     
  #     # Case for multiple selected traits
  #   } else if (nSelTraits() > 1) {
  #     outputDFComb <- NULL
  #     
  #     for (nSelT in 1:nSelTraits()) {
  #       # Use the index directly in the loop
  #       outputDF <- withProgress(message = paste('Running Computations for trait', nSelT), value = 0, {
  #         getRankedPredictedValues(predictionData(), nTraits(), unlist(Trait())[nSelT], GPModelST(), fixedX = Fixed(), fixedData = fixedData_List(), optTS())
  #       })
  #       
  #       # Combine columns - the second column is trait-specific values
  #       if (is.null(outputDFComb)) {
  #         outputDFComb <- outputDF  # Initialize with the first set (including UB of reliability)
  #       }else {
  #         outputDFComb <- merge(outputDFComb, outputDF,by="LineID")  # Append the trait-specific predicted values and UB reliability
  #       }
  #       
  #     }
  #     
  #     # Add proper column names for the combined data frame
  #     colnames(outputDFComb) <- c("GermplasmID", unlist(Trait()), "UpperBound of Reliability")
  #     
  #     return(outputDFComb)
  #   }
  # })
  # 
  
  
  temp_file7a <- reactiveVal('none')
  processComplete7a <- reactiveVal(FALSE)
  gpOutputListST <- reactiveVal(NULL)
  
  # Start the background process and return rProcess
  
  rProcess7a <- eventReactive(input$RunPredictionsST,{
    
    # Path for the temporary file
    temp_file7a(tempfile())

    if(checkSTData()==TRUE){
      # Start the process in a separate R process
      callr::r_bg(function(getRankedPredictedValues, predictionData, nTraits, Trait, GPModelST,Fixed,fixedData_List,optTS,temp_file7a_path,nSelTraits,cleanREPV2){
        library(bWGR)
        library(rrBLUP)
        library(NAM)
        sink(temp_file7a_path)
        
        if(nSelTraits==1){
          cat(paste("Running Computations for ",unlist(Trait),"\n",sep=""))    
          outputDFList <-   getRankedPredictedValues(predictionData, nTraits, unlist(Trait), GPModelST, fixedX = Fixed, fixedData = fixedData_List, optTS)
          return(outputDFList)
        }
        
        if(nSelTraits >1){
        
          trOutputDFComb <- NULL
          tstOutputDFComb <- NULL
          
          for(nSelT in 1:nSelTraits){
            cat(paste("Running Computations for ",unlist(Trait)[nSelT],"\n",sep=""))    
            outDF_List <- getRankedPredictedValues(predictionData, nTraits, unlist(Trait)[nSelT], GPModelST, fixedX = Fixed, fixedData = fixedData_List, optTS)
           
            if (is.null(trOutputDFComb)) {
              trOutputDFComb <- outDF_List[[1]]
              tstOutputDFComb <- outDF_List[[2]]
            }else {
              trOutputDFComb <- merge(trOutputDFComb,outDF_List[[1]],by="LineID")  # Append the trait-specific predicted values and UB reliability
             
              if(class(outDF_List[[2]])=="data.frame"){
                tstOutputDFComb <- merge(tstOutputDFComb,outDF_List[[2]],by="LineID")
              }
              
            }
          }
          
           cNamestrOut <- colnames(trOutputDFComb)
           colIndTr <- c(grep("LineID",cNamestrOut),grep("Predicted",cNamestrOut))
          
           cNameststOut <- colnames(tstOutputDFComb)
           colIndTst <- c(grep("LineID",cNameststOut),grep("Predicted",cNameststOut),grep("Reliability",cNameststOut)[1])
           
           outputDFList <- list(trOutputDFComb[,colIndTr],tstOutputDFComb[,colIndTst])
           cat("\n")
           return(outputDFList)
        }
        
        sink()
        
      }, args = list(
        getRankedPredictedValues = getRankedPredictedValues,
        predictionData = predictionData(),
        nTraits = nTraits(),
        Trait = unlist(Trait()),
        GPModelST=GPModelST(),
        Fixed=Fixed(),
        fixedData_List=fixedData_List(),
        optTS=optTS(),
        temp_file7a_path = temp_file7a(),
        nSelTraits=nSelTraits(),
        cleanREPV2=cleanREPV2
      ), stdout = temp_file7a(), stderr = temp_file7a())
      
    }else if(checkSTData()==FALSE){
      
      sink(temp_file7a())
      cat("Current data is not suitable for single trait in single environment GP.\n")
      sink()
      
    }
    
  })
  
####  
  
  STGPMssg <- reactiveVal(NULL)
  
  # Reactive expression to check process status and read temp file
  processStatus7a <- reactive({
    invalidateLater(1000, session) # Update every second
    
    if (!is.null(rProcess7a()) && rProcess7a()$is_alive()) {
      if (file.exists(temp_file7a())) {
        lines <- readLines(temp_file7a(), warn = FALSE)
        return(paste(lines, collapse = "\n"))
      } else {
        return("Waiting for output...")
      }
    } else if (!is.null(rProcess7a()) && !rProcess7a()$is_alive()) {
      if (file.exists(temp_file7a())) {
        lines <- readLines(temp_file7a(), warn = FALSE)
        txt <- paste(lines, collapse = "\n")
        processComplete7a(TRUE)
        gpOutputListST(rProcess7a()$get_result()) 
        STGPMssg(paste(txt, "\n Computations Completed", collapse = "\n"))
          return(paste(txt, "\n Computations Completed", collapse = "\n"))
      } else {
        return("Waiting for output...")
      }
    } else {
      return("Error")
    }
  })
  
  # Update the UI with process status
  output$messageSTGPRun <- renderText({
    processStatus7a()
  })
  
  
  

  ### GPModel MT
  
  GPModelMT <- reactive(input$GPModelMT)
  
  ## 
  
  
  # outputListMT <- eventReactive(input$RunPredictionsMT,{ 
  #   
  #   withProgress(message = 'Running Computations', value = 0, {
  #     getRankedPredictedValuesMT(predictionData(),nTraits(),unlist(Trait()),GPModelMT(),optTS())
  #   })
  # })
  
  
  temp_file7b <- reactiveVal('none')
  processComplete7b <- reactiveVal(FALSE)
  gpOutputListMT <- reactiveVal(NULL)
  
  # Start the background process and return rProcess
  
  rProcess7b <- eventReactive(input$RunPredictionsMT,{
    
    # browser()
    # Path for the temporary file
    temp_file7b(tempfile())
    
    if(checkMTData()==TRUE){
      
      # Start the process in a separate R process
      callr::r_bg(function(getRankedPredictedValuesMT, predictionData, nTraits, Trait, GPModelMT,optTS,temp_file7b_path,nSelTraits,cleanREPV2){
        
        library(BGLR)
        library(rrBLUP)
        library(NAM)
        library(sommer)
        sink(temp_file7b_path)
        
       #a <-  getRankedPredictedValuesMT(predictionData(),nTraits(),unlist(Trait()),GPModelMT(),optTS())
        
        if(nSelTraits>1){
          cat(paste("Running computations for multi-trait genomic prediction of ",paste(unlist(Trait),collapse=", "),"\n",sep=""))    
          outputDFList <- getRankedPredictedValuesMT(predictionData,nTraits,unlist(Trait),GPModelMT,optTS)
          cat("\n")
          return(outputDFList)
        }
     
        sink()
        
      }, args = list(
        getRankedPredictedValuesMT = getRankedPredictedValuesMT,
        predictionData = predictionData(),
        nTraits = nTraits(),
        Trait = unlist(Trait()),
        GPModelMT=GPModelMT(),
        optTS=optTS(),
        temp_file7b_path = temp_file7b(),
        nSelTraits=nSelTraits(),
        cleanREPV2=cleanREPV2
      ), stdout = temp_file7b(), stderr = temp_file7b())
      
    }else if(checkSTData()==FALSE){
      
      sink(temp_file7b())
      cat("Current data is not suitable for multi-trait genomic prediction.\n")
      sink()
      
    }
    
  })
  
  
  MTGPMssg <- reactiveVal(NULL)
  # Reactive expression to check process status and read temp file
  processStatus7b <- reactive({
    invalidateLater(1000, session) # Update every second
    
    if (!is.null(rProcess7b()) && rProcess7b()$is_alive()) {
      if (file.exists(temp_file7b())) {
        lines <- readLines(temp_file7b(), warn = FALSE)
        return(paste(lines, collapse = "\n"))
      } else {
        return("Waiting for output...")
      }
    } else if (!is.null(rProcess7b()) && !rProcess7b()$is_alive()) {
      if (file.exists(temp_file7b())) {
        lines <- readLines(temp_file7b(), warn = FALSE)
        txt <- paste(lines, collapse = "\n")
        processComplete7b(TRUE)
        gpOutputListMT(rProcess7b()$get_result()) 
        MTGPMssg(paste(txt, "\n Computations completed", collapse = "\n"))
        return(paste(txt, "\n Computations completed", collapse = "\n"))
      } else {
        return("Waiting for output...")
      }
    } else {
      return("Error")
    }
  })
  
  # Update the UI with process status
  output$messageMTGPRun <- renderText({
    processStatus7b()
  })
  

###### GPModel ME #YearME,LocationME,fixedME
  
  yrCol <- reactive(which(colnames(PhenoME()) %in% "Year"))
  locCol <- reactive(which(colnames(PhenoME()) %in% "Loc"))
 
  YearsME <- reactive({
    if(length(yrCol()) >0){
      levels(factor(PhenoME()[,yrCol()]))
    }
  })
  LocationsME <- reactive({
    if(length(locCol()) >0){
     levels(factor(PhenoME()[,locCol()]))
    }
  })
 
  
  observeEvent(input$infileBLUEsME, {
    updateSelectInput(inputId ="YearME",choices = c("All",YearsME()) )
  })
  
  observeEvent(input$infileBLUEsME, {
    updateSelectInput(inputId ="LocationME",choices = c("All",LocationsME()))
    
  })
  
  YearME <- reactive(input$YearME)
  LocationME <- reactive(input$LocationME)
  
  observeEvent(input$infileBLUEsME, {
    updateSelectInput(inputId ="fixedME",choices = IDColsME())
  })
  
  observeEvent(input$infileBLUEsME, {
    updateSelectInput(inputId ="EnvVarID",choices = IDColsME())
  })
  
  
  fixME <- reactive(input$fixedME)
  varEnv <- reactive(input$EnvVarID)
  GPModelME <- reactive(input$MEModelEnvR)
  
  GenK <- reactiveVal(NULL) 
  
  KGMethod <- reactive({
    if(input$GKernelMet =="Linear"){"GB"}else if(input$GKernelMet =="Gaussian"){"GK"}
  })
  
  fitEnvCovs <- reactive(input$fitEnvCov)
  
   
####   
  
  temp_file7c <- reactiveVal('none')
  processComplete7c <- reactiveVal(FALSE)
  gpOutputListSTME <- reactiveVal(NULL)
 
  
  # Start the background process and return rProcess
  
  rProcess7c <- eventReactive(input$RunPredictionsME,{
    
    #browser()
    # Path for the temporary file
    temp_file7c(tempfile())
    
    if(checkMEData()==TRUE){
      
      # Start the process in a separate R process
      callr::r_bg(function(getMEPred,DT_Filt_List,genoDat_List,TraitME,KG,EnvK,KGMethod,fitEnvCovs,fixME,varEnv,IDColsME,LocationME,YearME,temp_file7c_path,fitMEModels_Predictions){
        
        library(BGGE)
        library(EnvRtype)
        #library(bWGR)
        library(NAM)
        library(sommer)
        #source("GS_Pipeline_Jan_2022_FnsApp.R")
        sink(temp_file7c_path)
        
         cat(paste("Running computations for multi-environment genomic prediction of ",paste(unlist(TraitME),collapse=", "),"\n",sep=""))    
        
          if(fitEnvCovs == FALSE){
            outputDFList <- getMEPred(DT_Filt_List,genoDat_List,TraitME,KG=NULL,KE=NULL,KMethod= KGMethod,FitEnvModels = fitEnvCovs,fixedME=fixME,envVar=varEnv,IDColsME =IDColsME,LocME=LocationME,YrME=YearME)
          }else{ 
            outputDFList <- getMEPred(DT_Filt_List,genoDat_List,TraitME,KG=NULL,KE=EnvK,KMethod= KGMethod,FitEnvModels = fitEnvCovs,fixedME=fixME,envVar=varEnv,IDColsME =IDColsME,LocME=LocationME,YrME=YearME)
          }
         cat("\n")
         return(outputDFList)
        sink()
        
      }, args = list(
        getMEPred = getMEPred,
        DT_Filt_List = DT_Filt_List(),
        genoDat_List = genoDat_List(),
        TraitME = unlist(TraitME()),
        KG=GenK(),
        EnvK=EnvK(),
        KGMethod=KGMethod(),
        fitEnvCovs= fitEnvCovs(),
        fixME =fixME(),
        varEnv=varEnv(),
        IDColsME=IDColsME(),
        LocationME=LocationME(),
        YearME= YearME(),
        temp_file7c_path = temp_file7c(),
        fitMEModels_Predictions=fitMEModels_Predictions
      ), stdout = temp_file7c(), stderr = temp_file7c())
      
   }else if(checkMEData()==FALSE){
      
      sink(temp_file7c())      
      cat("Current data is not suitable for multi-environment genomic prediction.\n")
      sink()
      
    }
    
  })
  
  
  STMEGPMssg <- reactiveVal(NULL)
  # Reactive expression to check process status and read temp file
  processStatus7c <- reactive({
    invalidateLater(1000, session) # Update every second
    
    if (!is.null(rProcess7c()) && rProcess7c()$is_alive()) {
      if (file.exists(temp_file7c())) {
        lines <- readLines(temp_file7c(), warn = FALSE)
        return(paste(lines, collapse = "\n"))
      } else {
        return("Waiting for output...")
      }
    } else if (!is.null(rProcess7c()) && !rProcess7c()$is_alive()) {
      if (file.exists(temp_file7c())) {
        lines <- readLines(temp_file7c(), warn = FALSE)
        txt <- paste(lines, collapse = "\n")
        processComplete7c(TRUE)
        gpOutputListSTME(rProcess7c()$get_result()) 
        STMEGPMssg(paste(txt, "\n Computations completed", collapse = "\n"))
        return(paste(txt, "\n Computations completed", collapse = "\n"))
      } else {
        return("Waiting for output...")
      }
    } else {
      return("Error")
    }
  })
  
  # Update the UI with process status
  output$messageMEGPRun <- renderText({
    processStatus7c()
  })
  
  outputListMETab <- reactive({
    #browser()
    if(!is.null(gpOutputListSTME())){
     getCombinedTab(gpOutputListSTME(),TraitME(),IDColsME(),IDColME(),fitEnvCovs())
    }
  })

  
  
  

## Render Table for ST  
  
  
  STGPOut <- reactive(input$STGPOut)
  
  outputListST <- reactiveVal(NULL)
  
  output$Ranked_Lines_for_SelectionST <- renderDataTable({ 
    if(STGPOut()=="Target"){
      if(class(gpOutputListST()[[2]])=="data.frame"){ 
        outputListST(gpOutputListST()[[2]])
        return(gpOutputListST()[[2]][1:10,])
        
      }else if(class(gpOutputListST()[[2]])=="character"){
        gpOutputListST()[[2]]
      }
    }else if(STGPOut()=="Training"){
      outputListST(gpOutputListST()[[1]])
        return(gpOutputListST()[[1]][1:15,])
       
    }
  })
  
  
  
## Render Table for MT
  
  MTGPOut <- reactive(input$MTGPOut)
  outputListMT <- reactiveVal(NULL)
  
  output$Ranked_Lines_for_SelectionMT <- renderDataTable({ 
    
  if(nSelTraits()>1){
    
    if(MTGPOut()=="Target"){
      if(class(gpOutputListMT()[[2]])=="data.frame"){ 
        outputListMT(gpOutputListMT()[[2]])
        return(gpOutputListMT()[[2]][1:20,])
      }else if(class(gpOutputListMT()[[2]])=="character"){
        gpOutputListMT()[[2]]
      }
    }else if(MTGPOut()=="Training"){
      outputListMT(gpOutputListMT()[[1]])
      return(gpOutputListMT()[[1]][1:20,])
      
    }
  }else{NULL}
})
  

### 
  
 ## Render Table for ME
 # "Homogeneous Variance (G+E+GxE)","Heterogeneous Variance (G+E+GxEi)")
  
  output$Ranked_Lines_for_SelectionSTME <- renderTable({ 
    if(GPModelME()=="Main Effect (G+E)" || GPModelME()=="Main Effect (G+E+W)" ){
        (gpOutputListSTME()[[1]][[3]][[1]][1:15,])
    }else if(GPModelME()=="Homogeneous Variance (G+E+GxE)" || GPModelME()=="Homogeneous Variance (G+E+GxE+W)"){
        (gpOutputListSTME()[[1]][[3]][[2]][1:15,])
      
    }else if(GPModelME()=="Heterogeneous Variance (G+E+GxEi)" || GPModelME()== "Heterogeneous Variance (G+E+GxEi+W)"){
        (gpOutputListSTME()[[1]][[3]][[3]][1:15,])
    }
  })

  
  
  
  library(plotly)
  library(shiny)
  
  # Inside the server function of your Shiny app
  observe({
    
    req(gpOutputListST())  # Ensure gpOutputListST() exists before proceeding
   
    # browser()
    
    output$plots <- renderPlotly({  # Use renderPlotly for plotly plots
      if(nSelTraits() == 1 && STGPOut() == "Target" && !is.null(gpOutputListST()[[2]])){
        indX <- grep("Predicted Value", colnames(gpOutputListST()[[2]]))
        indY <- grep("Reliability", colnames(gpOutputListST()[[2]]))

        # Extract values from the data
        predicted_values <- gpOutputListST()[[2]][,indX]
        upper_bound_reliability <- gpOutputListST()[[2]][, indY]
        id <- gpOutputListST()[[2]][, "LineID"]  # Assuming "LineID" exists in the data

        # Create a data frame for plotly
        df <- data.frame(Predicted_Value = predicted_values,
                         Upper_Bound_of_Reliability = upper_bound_reliability,
                         ID = id)

        # Create the plot using plotly
        p <- plot_ly(data = df,
                     x = ~Predicted_Value,
                     y = ~Upper_Bound_of_Reliability,
                     type = 'scatter',
                     mode = 'markers',
                     text = ~paste("ID:", ID),  # Text to display on hover
                     hoverinfo = 'text')

        # Customize layout
        p <- p %>%
          layout(
                  #title = paste("Upper Bound of Reliability vs Predicted Values for", Trait()),
            title = list(
              text = paste("Upper Bound of Reliability vs Predicted Values for", Trait()), 
              y = 0.95,  # Adjust the vertical placement of the title (closer to top = 1)
              x = 0.5,  # Center the title horizontally
              xanchor = 'center',
              yanchor = 'top'
            ),
            margin = list(t = 100), 
            xaxis = list(title = "Predicted Value"),
            yaxis = list(title = "Upper Bound of Reliability"),
            font = list(size = 16))

       # Return the plotly object
        
        outFile_path <- paste(outResults_Dir(),"/","SingleTrait_GP_Plot_for_",Trait(),".html",sep="")
        
        # Save the plot as an HTML file
        saveWidget(p, file = outFile_path)
        
        p  
        
      }else {
        return(NULL)  # Return NULL if the condition is not met
      }
    })
    
  })

###  
  
 observe({
    req(gpOutputListST()) 
    req(STGPOut())
    req(nSelTraits())
    print(Trait())
    nPlots <- nSelTraits()  # Get the number of selected traits
    req(nPlots)  # Ensure nPlots is available
    
    # Hide all plotly outputs initially
    lapply(1:10, function(i) {
      shinyjs::hide(paste0("plot", i))  # Hides plots that are not needed
    })
    
    # Render only the required number of plots based on nSelTraits()

    if(nSelTraits()>1 && STGPOut()=="Target" && class(gpOutputListST()[[2]])=="data.frame"){
      
      indX <- grep("Predicted", colnames(gpOutputListST()[[2]]))
      indY <- grep("Reliability", colnames(gpOutputListST()[[2]]))
      
      lapply(1:nPlots, function(loc_i) { 
       
      
        shinyjs::show(paste0("plot", loc_i))  # Show only necessary plots
     
        #output[[noquote(paste("'","plot",loc_i,"'",sep=""))]]
        
        output[[paste0("plot",loc_i)]] <- renderPlotly({
       
        # Extract values from the data
        predicted_values <- gpOutputListST()[[2]][,indX[loc_i]]
        upper_bound_reliability <- gpOutputListST()[[2]][,indY[1]]
        id <- gpOutputListST()[[2]][, "LineID"]  # Assuming "LineID" exists in the data
        
        # Create a data frame for plotly
        df <- data.frame(Predicted_Value = predicted_values,
                       Upper_Bound_of_Reliability = upper_bound_reliability,
                       ID = id)
      
      # Create the plot using plotly
        p <- plot_ly(data = df,
                   x = ~Predicted_Value,
                   y = ~Upper_Bound_of_Reliability,
                   type = 'scatter',
                   mode = 'markers',
                   text = ~paste("ID:", ID),  # Text to display on hover
                   hoverinfo = 'text')
      
      # Customize layout
        p <- p %>% 
             layout(
               title = paste("Upper Bound of Reliability vs Predicted Values for ", unlist(Trait())[loc_i],sep=""),
               margin = list(t = 100),
               xaxis = list(title = "Predicted Value"),
               yaxis = list(title = "Upper Bound of Reliability"),
               font = list(size = 14))
    
      outFile_path <- paste(outResults_Dir(),"/","SingleTrait_GP_Plot_for_",unlist(Trait())[loc_i],".html",sep="")
      
      # Save the plot as an HTML file
      saveWidget(p, file = outFile_path)
      
      p
      
      
     })
    
    })
   } 
})
 
 
 
 # 
 # indX <- grep("Predicted",colnames(gpOutputListMT()[[2]]))
 # indY <- grep("Reliability",colnames(gpOutputListMT()[[2]]))
 
 # lapply(1:nSelTraits(),function(loc_i){
 #  
 #  output[[noquote(paste("'","plot",loc_i,"'",sep=""))]] <- renderPlot({
 #    plot.default(gpOutputListMT()[[2]][,indX[loc_i]],gpOutputListMT()[[2]][,indY[1]],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main=paste("Upper bound of Reliability vs Predicted Values for ",Trait()[loc_i()],sep=""),font.lab=2,cex.lab=1.25)
 #  })
 #  
 #     })
  
### 
 
    
  observe({
    req(gpOutputListMT()) 
    req(MTGPOut())
    req(nSelTraits())
    print(Trait())
    nPlots <- nSelTraits()  # Get the number of selected traits
    req(nPlots)  # Ensure nPlots is available
    
    # Hide all plotly outputs initially
    lapply(1:10, function(i) {
      shinyjs::hide(paste0("plotMT", i))  # Hides plots that are not needed
    })
    
    if(nSelTraits()>1 && MTGPOut()=="Target" && !is.null(gpOutputListMT()[[2]])){
     
     indX <- grep("Predicted", colnames(gpOutputListMT()[[2]]))
     indY <- grep("Reliability", colnames(gpOutputListMT()[[2]]))
    
     lapply(1:nPlots,function(loc_i){
    
       shinyjs::show(paste0("plotMT", loc_i))  # Show only necessary plots
       output[[paste0("plotMT",loc_i)]] <- renderPlotly({
        
       # Extract values from the data
       predicted_values <- gpOutputListMT()[[2]][,indX[loc_i]]
       upper_bound_reliability <- gpOutputListMT()[[2]][,indY[1]]
       id <- gpOutputListMT()[[2]][,"LineID"]  # Assuming "LineID" exists in the data
        
       # Create a data frame for plotly
       df <- data.frame(Predicted_Value = predicted_values,
                         Upper_Bound_of_Reliability = upper_bound_reliability,
                         ID = id)
        
       # Create the plot using plotly
        p <- plot_ly(data = df,
                     x = ~Predicted_Value,
                     y = ~Upper_Bound_of_Reliability,
                     type = 'scatter',
                     mode = 'markers',
                     text = ~paste("ID:", ID),  # Text to display on hover
                     hoverinfo = 'text')
        
        # Customize layout
        p <- p %>% 
          layout(
            title = paste("Upper Bound of Reliability vs Predicted Values for ", unlist(Trait())[loc_i],sep=""),
            margin = list(t = 100),
            xaxis = list(title = "Predicted Value"),
            yaxis = list(title = "Upper Bound of Reliability"),
            font = list(size = 14))
        
        outFile_path <- paste(outResults_Dir(),"/","Multitrait_GP_Plot_for_",unlist(Trait())[loc_i],".html",sep="")
        
        # Save the plot as an HTML file
        saveWidget(p, file = outFile_path)
        
        p
        
        
      })
      
     })
  
    }
    
  })
  
### 
  
######## Table Export Handlers
  
  output$ExportImpGenoDF <- downloadHandler(
    filename = function() {
      "ImputedGenotypesDF.txt"
    },
    content = function(file) {
      
      write.table(as.data.frame(GenoImp_DF()), file, row.names = FALSE)
    }
  )
  
  
  output$ExportOptTS <- downloadHandler(
    filename = function() {
      "OptimalTS.csv"
    },
    content = function(file) {
      
      write.csv(as.data.frame((Train_STPGA()[[2]])), file, row.names = FALSE)
    }
  )
  
  
  output$ExportOut <- downloadHandler(
    filename = function() {
      "PredictionsST.csv"
    },
    content = function(file) {
      
      write.csv(as.data.frame(outputListST()), file, row.names = FALSE)
    }
  )
  
  output$ExportOutMT <- downloadHandler(
    filename = function() {
      "PredictionsMT.csv"
    },
    content = function(file) {
      
      write.csv(as.data.frame(outputListMT()), file, row.names = FALSE)
    }
  )
  
### Export ME output
  
  output$ExportOutME <- downloadHandler(
    
    filename = function() {
      "PredictionsME.csv"
    },
    content = function(file) {
      write.csv(as.data.frame(outputListMETab()), file, row.names = FALSE)
    }
  )
  
 
##### Documentation  
  
  outResults_Dir <- reactiveVal(NULL)
  
  observe({
    
   
    sessionID <- substring(session$token,1,10)
    
    log_file <- paste("GS4PB_App_log-", Sys.Date(),"-",sessionID,".txt",sep="")
    
    req(selected_dir())
    # Use the selected output directory
    output_dir <- selected_dir()
    
    # Create the "Results" folder inside the selected directory if it doesn't exist
    results_dir <- file.path(output_dir, "Results")
    if (!dir.exists(results_dir)) {
      dir.create(results_dir, recursive = TRUE)
    }
    req(results_dir)
    
    outResults_Dir(results_dir)
    
    print(outResults_Dir())
    
    # Define the path to save the file in the Results folder
    output_path <- file.path(results_dir,log_file)
    
    print(output_path)
    
####  Doc parameters for geno data
    
    if(impMethod() == "LDKNNI"){ 
      
      impDocOut <-  c( 
        paste("Imputed geno data using default method : ",impMethod()," with parameters including number of high LD sites : ", l()," and # neighboring samples : ",k(),sep=""),  
        paste(mssgImpGeno1StatsOut(),"\n",sep="")
      )
    }else if (impMethod() == "Numeric"){
      impDocOut <-  c( 
        paste("Imputed geno data using method : ",impMethod()," with parameters including number of nearest neighbors : ",nN()," and distance metric : ", Dist(),sep=""),  
        paste(mssgImpGeno1StatsOut(),"\n",sep="")
      )
      
    }else if (impMethod() == "AlphaPlantImpute") {
      
      impDocOut <- c(
        paste("Imputed geno data using method : ",impMethod(),"\n",sep=""),
        paste("Built library using the following parameters including nHaplotypes : ",nHap()," nSampleRounds : ", nSampRnds(), " hd Threshold : ",HDthresh(),sep=""),
        #paste(mssgAPIBuildOut(),"\n",sep=""),
        paste("Imputed using founder file ",as.character(founder()),sep=""),       
        paste(mssgAPIImputeOut(),"\n",sep="") 
      )
    }
    
    
    ##*(ncol(Geno())-5)
    
    docmnGeno <- c(
      
      paste("GS4PB App Log on : ", date(), "\n",sep=""),
      "###### Impute Geno Data ########" ,
      paste("Loaded genodata from file : ", input$infileVCF$name,"\n",sep=""),
      paste(mssgGenoStatsOut(), "\n",sep=""),
      "###### Filter Geno Data ########",
      paste("Filtered sites using 'Min Site Count': ",siteMinCnt()," and 'Min Allele Freq' : ", MAF() ,sep=""),  
      paste(mssgGenoFilt1StatsOut(),"\n",sep=""),    
      "########",
      paste("Filtered taxa using 'Min Not Missing' : ",minNotMissing() ,sep=""),
      paste(mssgGenoFilt2StatsOut(),"\n",sep=""),
      "###### Impute Geno Data ########",
      impDocOut,
      "###### Target Genotypes ########",
      paste("Loaded target lines from file : ", input$infileTargetTable$name,"\n",sep=""),
      paste(mssgTargetStatsOut(), "\n",sep="")
      
    )
    
    
    
    ### Set doc parameters for pheno data 
    if(isPhenoME()){ 
      PhFile <- input$infileBLUEsME$name
      PhData <- paste("Pheno data for multi-environmental trials loaded from file ",PhFile,"\n",sep="")
      IDCols <- IDColME() 
      strainCol <- StrainME()
      traitCols <- TraitCols()
      traits <- TraitME()
      
    }else if(!isPhenoME()){ 
      
      PhFile <- input$infileBLUEsSE$name
      PhData <- paste("Pheno data from single environment loaded from file ",PhFile,"\n",sep="")
      IDCols <- IDColSE() 
      strainCol <- StrainSE()
      traitCols <- TraitsPhSEVec()
      traits <- Trait()
    }
    
    
    docmnPh <- c(
      "###### Pheno Data ########" ,
      as.character(PhData),
      paste("Unique ID column was set to ",IDCols,sep=""),
      paste("Strain column was set to ",strainCol,sep=""),
      paste("The following trait(s) were selected : ",traitCols,sep=""),
      "\n",
      paste(mssgPhenoStatsOut(),"\n",sep=""),
      "\n",
      "##### Traits Data ######",
      paste("The following trait(s) were selected : ",traits,sep=""),
            paste(mssgTraitsStatsOut(),"\n",sep="")
            
    )
      
    if(gaussVar()){
        kern <- "Gaussian"
    }else if( !gaussVar()){
        kern <- "Linear"
    }
      
      docmnEnv <-  c(   
        paste("Environmental coordinates were uploaded from file : " ,input$inFileLocCoord$name,"\n",sep=""),
        "######\n", 
        paste("Enviromics variable collection was set for the period starting from ", input$startDate, "to ", input$endDate," ","\n",sep=""),
        paste("Use processed weather variables is set to : ",input$processWth,"\n",sep=""),
        paste(kern," kernel was used for estimating environmental relationship matrix",sep=""),"\n"
      ) 
      
      
      
      docmnGPMerge <- c( 
        "##### Geno Data in Merged Table ######",
        paste(mssgMergedGenoStatsOut(),"\n",sep=""),
        "##### Pheno Data in Merged Table ######",
        paste(mssgMergedPhenoStatsOut(),"\n",sep="")
      )
      
      ####
      
      docmnOTP <- c("")
      
      ####
      docmnCV <- c(
        paste("Single Trait CV was performed for ",STCVMssg(),sep=""),"\n",
        paste("Multi trait CV was performed for ",MTCVMssg(),sep=""),"\n",
        paste("Multi-environmental CV was performed for ",MECVMssg(),sep=""),"\n")
      
      ####
      docmnGP <- c(
        paste("Single Trait GP was performed for ",STGPMssg(),sep=""),"\n",
        paste("Multi trait GP was performed for ",MTGPMssg(),sep=""),"\n")
        #paste("Multi-environmental GP was performed for ",MEGPMssg(),sep=""),"\n")
      
      #### 
      
      documentation <- c(docmnGeno,docmnPh,docmnEnv,docmnGPMerge,docmnOTP,docmnCV,docmnGP)    
      
     
      
      writeLines(paste(documentation,collapse="\n"),output_path)
      
      # Copy the processed file to the "Results" directory
      # file.copy(log_file, output_path, overwrite = TRUE)
      
      
  }) 
    

}






