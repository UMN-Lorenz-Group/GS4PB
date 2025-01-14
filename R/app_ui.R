#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinyjs
#' @import shinyFiles
#' @import renv
#' @import fs
#' @import DT
#' @import htmlwidgets
#' @import plotly
#' @noRd


app_ui <- function(request){
  
 fluidPage( 
    
    tags$style(HTML(".message-text { color: red; }")),
    useShinyjs(),
    tags$br(),
  
    tags$br(),
    tags$head(
      tags$style(HTML("
             #title-row {
                      color: white;
                      background-color: black;
                      padding: 10px;
                      font-size: 20px;
                      text-align: center;
      }
      "))
    ),
    
    fluidRow(
      id = "title-row", 
      column(3),
      column(6, div(
        titlePanel(tags$strong("GS4PB App")),
        tags$p("Genomic Selection For Plant Breeding Applications")
      ))
    ),
    
    
    
    tabsetPanel(id="inData",
                tabPanel("Home",
                         fluidRow(
                           column(5),column(width=4,tags$h3(tags$strong("GS4PB App")))
                           # column(6,offset=6),actionButton(inputId ="Home_Data", "next")
                         ),
                         fluidRow(
                           column(2),column(width=8,tags$h5("The GS4PB App implements a genomic selection pipeline. 
                        The GS pipeline involves steps starting from data management to making selection decisions based on genome estimated breeding values.The steps 
                        in the pipeline are depicted in the schematic below.") 
                           )),
                         tags$br(),
                         tags$br(),
                         fluidRow(
                           column(2),column(width=10,tags$a(tags$img(src="www/GSPipeline.png",title="Genomic Selection Pipeline",width="850",height="500")))), 
                         
                         fluidRow(
                           column(2),column(width=9,tags$h6(("Contributors: Vishnu Ramasubramanian, Cleiton Wartha, Paolo Vitale, Sushan Ru,and Aaron Lorenz")))),
                         tags$br(),        
                         fluidRow(
                           column(2),column(width=9,tags$h6("Contact: Vishnu Ramasubramanian - vramasub@umn.edu")))
                ),
                
                tabPanel("Start GS Pipeline",
                         sidebarLayout(
                           sidebarPanel(
                             tags$br(),
                             
                             # Text input for setting the output directory, where log file and results folder will be saved
                             # textInput("outdir", "Set Output Directory", value = "C:/myoutputdir
                             
                             tags$h5(tags$strong("Set output dir to save log and result files ")),
                             tags$br(),
                             # Button to trigger directory selection
                             shinyDirButton("directory", "Select Output Directory", "Please select a folder"),
                             tags$br(),
                             # Display selected directory path
                             verbatimTextOutput("dir_path")
                           ),
                           mainPanel(
                             tags$br(),
                             tags$h4(tags$strong("GS Pipeline Instructions")),
                             
                             tags$br(),
                             fluidRow(
                               column(width=12,tags$h5("The GS pipeline involves steps ranging from genotypic and phenotypic data processing, merging geno and pheno data, opimizing training populations, 
                            cross validation of genomic prediction models and genomic prediction (figure below). Some of the steps are mandatory, whereas other steps are optional.  
                            Once a step is done, the user can navigate to the next step in the tab panel. 
                            The user also needs to select a directory for storing log file and results. The log file records all the input parameter values 
                            and all of the results generated are stored in the 'Results' folder. It is also possible to use the app for genotypic/phenotypic data processing alone.") 
                               )),
                             tags$br(),
                             tags$br(),
                             fluidRow(
                               column(width=12,tags$a(tags$img(src="www/GSPipelineInstructions.png",title="GS Pipeline Steps",width="850",height="500")))), 
                           ),
                         ),      
                ),    
                ## 
                tabPanel("Genotypic Data Processing",
                         tabsetPanel(id="GenDatProcess",
                                     
                                     ## Tab for loading data
                                     
                                     tabPanel("Load Genotypic Data",
                                              sidebarLayout(
                                                sidebarPanel(
                                                  # fluidRow(
                                                  #   column(1),column(width=8,
                                                  #                    div(style="display: inline-block; vertical-align: top;",
                                                  #                        checkboxInput("InGenoFormat1","",TRUE)
                                                  #                    ),
                                                  #                    div(style="display: inline-block; vertical-align: top; padding-left: 10px;",
                                                  #                        tags$h6(tags$strong("1) Training & Target genotypic data in a combined VCF & target line IDs in a csv file"))
                                                  #                    )
                                                  #   )), 
                                                  
                                                  tags$br(),
                                                  selectInput(inputId="TargetIDCol","Choose Target ID Col",choices=NULL,multiple=FALSE),
                                                  fluidRow(
                                                    tags$h6(tags$strong(HTML("&nbsp;&nbsp;&nbsp; Select the column with targetIDs, if available.")))), 
                                                  fluidRow(
                                                    tags$h6(tags$strong(HTML("&nbsp;&nbsp;&nbsp; If not, leave it blank.")))),
                                                  tags$br()
                                                  
                                                ),
                                                
                                                mainPanel(
                                                  
                                                  fluidRow(
                                                    column(3),column(width=4,tags$h3(tags$strong("Load Genotypic Data"))),   
                                                  ), 
                                                  tags$br(),
                                                  fluidRow(
                                                    column(1),column(width=10,tags$h5(tags$strong("Upload genotypic data in VCF file format. Genotypic data could include the training and target data in the same file 
                                                                        or just the training data alone.If your data includes a target population, please upload the line IDs of target population in the right panel. 
                                                                        If not leave it blank. Once you upload the target line IDs in csv format, select tehe column containing the target IDs")))
                                                  ),
                                                  tags$br(),
                                                  fluidRow(
                                                    column(1),column(width=4,tags$h4(tags$strong("Genotypic Data"))),
                                                    column(6,offset=0,tags$h4(tags$strong("Target Population IDs"))),
                                                  ),
                                                  
                                                  
                                                  tags$br(),
                                                  # conditionalPanel(condition="input.InGenoFormat1 == true",
                                                  tags$br(),
                                                  fluidRow(
                                                    column(width=5,tags$h5(tags$strong("Upload Complete Genotypic Data (Training +Target) (VCF)"))),
                                                    column(6),column(width=5,tags$h5(tags$strong("Upload Target Line IDs (CSV) file"))),
                                                  ),
                                                  fluidRow(
                                                    column(width=5,fileInput("infileVCF", "", accept = ".vcf")),
                                                    column(6),column(width=5,fileInput("infileTargetTable", "", accept = ".csv")),
                                                  ),
                                                  
                                                  fluidRow(
                                                    column(1),column(width=5,checkboxInput("header", "Header", TRUE)),
                                                    column(6),column(width=5,checkboxInput("header", "Header", TRUE)) 
                                                  ),
                                                  #),
                                                  
                                                  tags$head(
                                                    tags$style(HTML("
                                          #messageGenoStats {
                                          
                                          /*max-height: 1200px; Set maximum height */
                                          overflow-y: scroll; /* Enable vertical scrolling */
                                          overflow-x: scroll;  /*Hide horizontal scrolling */
                                          overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                          width: 300px; 
                                          /*max-width: 100%; */
                                          padding: 6px 12px;
                                          height: 250px;
                                         }
                                        "))
                                                  ),
                                                  
                                                  tags$head(
                                                    tags$style(HTML("
                                #messageTargetStats {
                                          
                                          /*max-height: 1200px; Set maximum height */
                                          overflow-y: scroll; /* Enable vertical scrolling */
                                          overflow-x: scroll;  /*Hide horizontal scrolling */
                                          overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                          width: 300px; 
                                          /*max-width: 100%; */
                                          padding: 6px 12px;
                                          height: 250px;
                                         }
                                        "))
                                                  ),
                                                  
                                                  fluidRow(column(width=5,tags$h6(verbatimTextOutput("messageGenoStats"))),
                                                           column(7),column(width=4,tags$h6(verbatimTextOutput("messageTargetStats")))),
                                                  tags$br(),
                                                  tags$br(),
                                                )
                                              )
                                     ),  
                                     
                                     
                                     ###  
                                     tabPanel("Filter Genotypic Data",
                                              
                                              tags$br(),
                                              tags$br(),
                                              fluidRow(column(1),
                                                       column(width=10,"Filter sites (markers) and taxa (lines) in genotype table using rTASSEL (Monier et al. 2021), if you uploaded a raw genotype file or a QC genotype file that requires filtering. For example, it is common to 
                    set minimum number of sites that are not missing to 80% of the lines in the input genotype data and set the MAF threshold to 0.02/0.05. 
                    If you uploaded a QC genotype file, you can choose to skip this step")),
                                              
                                              tags$br(),
                                              tags$br(),
                                              
                                              fluidRow(
                                                column(1),column(width=3,tags$strong(tags$h4("Filter Genotype Table Sites"))),
                                                #Filt2
                                                column(5,offset=2,tags$strong(tags$h4("Filter Genotype Table Taxa")))),
                                              tags$br(),
                                              fluidRow(
                                                #Filt1
                                                column(1),column(width=3,
                                                                 numericInput(inputId="siteMinCnt","Minimum Site Count (TASSEL)",value = 0,min =0, max=0)),
                                                #Filt2
                                                column(5,offset=2,numericInput(inputId="minNotMissing","Minimum Not Missing (TASSEL)",value = 0.9,min =0.5, max=1))),
                                              tags$br(),
                                              
                                              tags$br(),
                                              fluidRow(
                                                #filt1
                                                column(1),column(width=3,
                                                                 numericInput(inputId="MAF","Minimum Allele Frequency",value =0.02,min =0, max=0.5)),
                                                #Filt2
                                                column(5,offset=2,
                                                       actionButton(inputId="FilterTaxa","Filter Genotype Table Taxa"))),
                                              
                                              tags$br(),
                                              fluidRow(
                                                #Filt1
                                                column(1),column(width=3,
                                                                 actionButton(inputId="FilterSites","Filter Genotype Table Sites"))),
                                              
                                              tags$br(),
                                              tags$br(),
                                              fluidRow(
                                                column(1),column(width=3,
                                                                 checkboxInput("setGenoFilt1Tas", "Use Filtered Genotypes From Filter Sites For Next Steps", TRUE)),
                                                column(5,offset=2,
                                                       checkboxInput("setGenoFilt2Tas", "Use Filtered Genotypes From Filter Taxa For Next Steps", TRUE))),  
                                              tags$br(),
                                              tags$br(),
                                              
                                              tags$br(),
                                              tags$br(),
                                              # 
                                              
                                              tags$head(
                                                tags$style(HTML("
                                          #messageGenoFilt1 {
                                          
                                          /*max-height: 1200px; Set maximum height */
                                          overflow-y: scroll; /* Enable vertical scrolling */
                                          overflow-x: scroll;  /*Hide horizontal scrolling */
                                          overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                          width: 350px; 
                                          /*max-width: 100%; */
                                          padding: 6px 12px;
                                          height: 400px;
                                         }
                                        "))
                                              ),
                                              
                                              tags$head(
                                                tags$style(HTML("
                                          #messageGenoFilt2 {
                                          
                                          /*max-height: 1200px; Set maximum height */
                                          overflow-y: scroll; /* Enable vertical scrolling */
                                          overflow-x: scroll;  /*Hide horizontal scrolling */
                                          overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                          width: 350px; 
                                          /*max-width: 100%; */
                                          padding: 6px 12px;
                                          height: 400px;
                                         }
                                        "))
                                              ),
                                              
                                              
                                              
                                              fluidRow(column(1),column(width=5,tags$h6(verbatimTextOutput("messageGenoFilt1"))),
                                                       column(6),column(width=5,tags$h6(verbatimTextOutput("messageGenoFilt2")))),
                                              
                                              
                                              tags$br(),
                                              tags$br(),
                                              tags$br()
                                     ),
                                     
                                     tabPanel("Impute Genotypic Data",  
                                              sidebarLayout(
                                                sidebarPanel(
                                                  tags$strong(tags$h4("Impute Genotype Table ")),
                                                  tags$br(),
                                                  selectInput(inputId="imputeMet","Select Imputation Method",choices=c("LDKNNI","Numeric","AlphaPlantImpute"),multiple=FALSE,selected="LDKNNI"),
                                                  tags$br(),
                                                  
                                                  # 
                                                  conditionalPanel(condition = "input.imputeMet == 'LDKNNI'",
                                                                   numericInput(inputId="l","Number of High LD Sites",value = 30,min =30,max=1000),
                                                                   numericInput(inputId="k","Neighboring Samples",value =10,min =10, max=100),
                                                  ),
                                                  
                                                  
                                                  conditionalPanel(condition = "input.imputeMet == 'Numeric'",
                                                                   numericInput(inputId="nN","Number of Nearest Neighbors",value = 5,min =2,max=100),
                                                                   selectInput(inputId="Dist","Distance",choices=c("Euclidean", "Manhattan", "Cosine"),multiple=FALSE,selected="Euclidean"),
                                                  ),
                                                  
                                                  conditionalPanel(condition = "input.imputeMet == 'Numeric' || input.imputeMet == 'LDKNNI'",
                                                                   actionButton(inputId="Impute","Impute Genotype Scores"),
                                                                   tags$br(),
                                                                   tags$br(),
                                                                   tags$br(),
                                                  ),
                                                  tags$br(),
                                                  tags$br(),
                                                  conditionalPanel(condition = "input.imputeMet == 'AlphaPlantImpute'",
                                                                   fileInput("inHapLib", "Load Prebuilt Haplotype Library (.phase)", accept = ".phase"),
                                                                   tags$br(),
                                                                   tags$br(),
                                                                   tags$br(),
                                                  ),
                                                  
                                                  checkboxInput("setGenoImpTas", "Use Imputed Genotypes For Next Steps", TRUE), 
                                                  tags$br(),
                                                  downloadButton("ExportImpGenoDF", "Export Imputed Genotypes"),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$br()
                                                ),
                                                mainPanel(
                                                  
                                                  conditionalPanel(condition = "input.imputeMet == 'LDKNNI' || input.imputeMet == 'Numeric'",
                                                                   tags$br(),
                                                                   tags$br(),
                                                                   fluidRow(
                                                                     column(width=10,tags$h5(
                                                                       "Impute missing genotype scores using rTASSEL (Monier et al. 2021), if you uploaded a raw genotype file or a QC genotype file with missing scores. 
                     Available options include numeric imputation and imputation using LD-K-Nearest Neighbors method are availble. 
                     For LD-KNN imputation set parameters l and K (Money et al. 2015). l corresponds to the number of high LD sites  
                     and k corresponds to the number of neighboring samples that are used to impute scores."))
                                                                   ),
                                                                   tags$br(),
                                                                   
                                                                   tags$head(
                                                                     tags$style(HTML("
                                           
                                            #messageImpGeno1 {
                                             /*max-height: 1200px; Set maximum height */
                                              overflow-y: scroll; /* Enable vertical scrolling */
                                              overflow-x: scroll;  /*Hide horizontal scrolling */
                                              overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                              width: 350px; 
                                              /*max-width: 100%; */
                                              padding: 6px 12px;
                                              height: 400px;
                                             }"))
                                                                   ),
                                                                   
                                                                   tags$br(),
                                                                   tags$br(),
                                                                   fluidRow(column(1),column(width=4,tags$h6(verbatimTextOutput("messageImpGeno1")))),
                                                                   tags$br(),
                                                                   tags$br(),
                                                  ),
                                                  conditionalPanel(condition = "input.imputeMet == 'AlphaPlantImpute'",
                                                                   uiOutput("APIUI"),
                                                  )
                                                )
                                              )),
                         )),
                
                ## Tab to load phenotypic data and select trait  
                tabPanel("Phenotypic Data Processing ",                  
                         
                         tabsetPanel(id="Phenotypic Data",
                                     
                                     tabPanel("Load Phenotypic Data and Traits",             
                                              ####
                                              sidebarLayout(
                                                sidebarPanel(
                                                  checkboxInput("chkPhenoME",tags$h5(tags$strong("Load Pheno Data from Multi-Environs")),FALSE),
                                                  tags$br(),
                                                  tags$br(),
                                                  conditionalPanel(condition="input.chkPhenoME == false",
                                                                   selectInput(inputId="IDColSE","Choose Uniq ID Col",choices=NULL,multiple=FALSE),
                                                                   selectInput(inputId="strainSE","Choose Strain Col",choices=NULL,multiple=FALSE),
                                                                   shinyWidgets::numericRangeInput(inputId="traitColsNum","Choose all the trait Cols",value=NULL,min=1,max=1),
                                                                   tags$br(),
                                                  ),
                                                  conditionalPanel(condition="input.chkPhenoME == true",
                                                                   selectInput(inputId="traitCols","Choose all the trait Cols",choices=NULL,multiple=TRUE),
                                                                   selectInput(inputId="IDColME","Choose Uniq ID Col",choices=NULL,multiple=FALSE),
                                                                   selectInput(inputId="strainME","Choose Strain Col",choices=NULL,multiple=FALSE),
                                                                   tags$br(),
                                                                   
                                                  ),
                                                  
                                                ),
                                                mainPanel(
                                                  tags$br(),
                                                  fluidRow(
                                                    column(width=10,tags$h5(
                                                      "Load phenotypic data and select traits for cross validation and genomic prediction. Check the box in the side panel to load phenotypic data from multi-environmental trials. 
                      For pheno file format specific to SE or ME trials, press the 'i' info button. Once the pheno file is loaded, select one or more traits."))
                                                  ),
                                                  tags$br(),
                                                  
                                                  conditionalPanel(condition="input.chkPhenoME == false",
                                                                   
                                                                   tags$br(),
                                                                   fluidRow(
                                                                     column(width=5,tags$h4(tags$strong("Phenotypic Data"))),
                                                                     column(6),column(7),column(width=5,tags$h4(tags$strong("Select trait(s)"))),
                                                                   ),
                                                                   tags$br(),
                                                                   
                                                                   fluidRow(
                                                                     column(width=4,fileInput("infileBLUEsSE", "Choose Phenotype File (CSV)", accept = ".csv")),
                                                                     column(width=1,actionButton("iPh", "i")),column(7),
                                                                     column(width=5,selectInput(inputId="trait","Choose One or More Traits",choices=NULL,multiple=TRUE)),
                                                                   ),
                                                                   fluidRow(
                                                                     column(width=5,checkboxInput("header", "Header", TRUE)),
                                                                   ),
                                                                   tags$br(),
                                                                   
                                                                   
                                                                   
                                                                   tags$br(),
                                                                   tags$head(
                                                                     tags$style(HTML("
                                                #messagePhSE {
                                                  
                                                  /*max-height: 1200px; Set maximum height */
                                                  overflow-y: scroll; /* Enable vertical scrolling */
                                                  overflow-x: scroll;  /*Hide horizontal scrolling */
                                                  overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                                  width: 350px; 
                                                  /*max-width: 100%; */
                                                  padding: 6px 12px;
                                                  height: 150px;
                                                 }
                                                "))
                                                                   ),
                                                                   
                                                                   tags$head(
                                                                     tags$style(HTML("
                                              #messageTrtSE {
                                              
                                              /*max-height: 1200px; Set maximum height */
                                              overflow-y: scroll; /* Enable vertical scrolling */
                                              overflow-x: scroll;  /*Hide horizontal scrolling */
                                              overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                              width: 350px; 
                                              /*max-width: 100%; */
                                              padding: 6px 12px;
                                              height: 150px;
                                             }
                                            "))
                                                                   ),
                                                                   
                                                                   fluidRow(
                                                                     column(width = 5,
                                                                            verbatimTextOutput("messagePhSE")
                                                                     ), column(7),column(width=5, verbatimTextOutput("messageTrtSE")),
                                                                     
                                                                   ),
                                                                   
                                                                   
                                                  ),
                                                  
                                                  conditionalPanel(condition="input.chkPhenoME == true",
                                                                   
                                                                   tags$br(),
                                                                   fluidRow(
                                                                     column(width=5,tags$h4(tags$strong("Phenotypic Data"))),
                                                                     column(6),column(width=3,tags$h4(tags$strong("Select trait(s)"))),
                                                                   ),
                                                                   fluidRow(
                                                                     column(width=3,fileInput("infileBLUEsME", "Choose Phenotype File (CSV)", accept = ".csv")),
                                                                     column(width=1,actionButton("iPhME", "i")),
                                                                     column(7),column(width=5,selectInput(inputId="traitME","Choose One or More Traits",choices=NULL,multiple=TRUE)),
                                                                   ),
                                                                   
                                                                   fluidRow(
                                                                     column(width=5,checkboxInput("headerME", "Header", TRUE)),
                                                                   ),
                                                                   
                                                                   
                                                                   tags$br(),
                                                                   tags$br(),
                                                                   tags$head(
                                                                     tags$style(HTML("
                                          #messagePhME {
                                          
                                          /*max-height: 1200px; Set maximum height */
                                          overflow-y: scroll; /* Enable vertical scrolling */
                                          overflow-x: scroll;  /*Hide horizontal scrolling */
                                          overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                          width: 350px; 
                                          /*max-width: 100%; */
                                          padding: 6px 12px;
                                          height: 150px;
                                         }
                                        "))
                                                                   ),
                                                                   
                                                                   tags$head(
                                                                     tags$style(HTML("
                                          #messageTrtME {
                                          
                                          /*max-height: 1200px; Set maximum height */
                                          overflow-y: scroll; /* Enable vertical scrolling */
                                          overflow-x: scroll;  /*Hide horizontal scrolling */
                                          overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                          width: 350px; 
                                          /*max-width: 100%; */
                                          padding: 6px 12px;
                                          height: 150px;
                                         }
                                        "))
                                                                   ),
                                                                   
                                                                   
                                                                   
                                                                   fluidRow(
                                                                     
                                                                     column(width = 5,
                                                                            verbatimTextOutput("messagePhME")
                                                                     ), column(7),column(width=5, verbatimTextOutput("messageTrtME")),
                                                                     
                                                                   ),
                                                                   
                                                                   tags$br()
                                                                   # uiOutput("LocDistribution")
                                                                   
                                                  )
                                                )
                                              )
                                     ))),
                
                
                ### GenoPheno Tab
                tabPanel("Geno-Pheno Data Merge",
                         
                         sidebarLayout(
                           sidebarPanel(
                             
                             conditionalPanel(condition="input.chkPhenoME == false",
                                              tags$br(),
                                              tags$h5(tags$strong(" Merge Geno and Pheno Data")),
                                              tags$br(),
                                              actionButton("MergeGP", "Merge"),
                                              tags$br(),
                             ),
                             
                             conditionalPanel(condition="input.chkPhenoME == true",
                                              tags$br(),
                                              tags$h5(tags$strong(" Merge Geno and Pheno Data")),
                                              tags$br(),
                                              actionButton("MergeGPME", "Merge"),
                                              tags$br(),
                             )
                             
                           ),
                           mainPanel(
                             
                             conditionalPanel(condition="input.chkPhenoME == false",
                                              
                                              tags$br(),
                                              fluidRow(
                                                column(width=10, tags$p("Merge genotypic and phenotypic data tables and check if the merge is successful. A successful merge is contingent on proper match of Line IDs. 
                              LineIDs with no matching IDs in either of the tables are printed in the message box.") 
                                                )),
                                              
                                              fluidRow(
                                                column(1),column(width = 4,tags$h5(tags$strong("Geno Data in Merged Table"))),
                                                column(6),column(width = 4, tags$h5(tags$strong("Pheno Data in Merged Table")))
                                              ),
                                              
                                              tags$br(),
                                              tags$head(
                                                tags$style(HTML("
                                        #messageMergedGeno{

                                                  /*max-height: 1200px; Set maximum height */
                                                  overflow-y: scroll; /* Enable vertical scrolling */
                                                  overflow-x: scroll;  /*Hide horizontal scrolling */
                                                  overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                                  width: 350px;
                                                  /*max-width: 100%; */
                                                  padding: 6px 12px;
                                                  height: 150px;
                                                 }
                                          "))
                                              ),
                                              
                                              tags$head(
                                                tags$style(HTML("
                                          #messageMergeTrtSE{

                                              /*max-height: 1200px; Set maximum height */
                                              overflow-y: scroll; /* Enable vertical scrolling */
                                              overflow-x: scroll;  /*Hide horizontal scrolling */
                                              overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                              width: 350px;
                                              /*max-width: 100%; */
                                              padding: 6px 12px;
                                              height: 150px;
                                             }
                                            "))
                                              ),
                                              
                                              fluidRow(
                                                column(width = 5,verbatimTextOutput("messageMergedGeno")),column(6),column(7),column(width=5, verbatimTextOutput("messageMergeTrtSE")),
                                                
                                              ),
                                              
                                              tags$br(),
                                              fluidRow(
                                                column(width = 4), column(6),column(width=4, verbatimTextOutput("ReadResult"))),
                                              
                                              tags$br(),
                                              
                                              tags$head(
                                                tags$style(HTML("
                                        #MssgMergedData {
                                         color: red;
                                        }")
                                                )
                                              ),
                                              fluidRow(
                                                column(2),column(width=8,tags$h5(tags$strong(textOutput("MssgMergedData"))))),
                                              
                                              tags$br(),
                                              tags$br(),
                             ),
                             conditionalPanel(condition="input.chkPhenoME == true",
                                              
                                              tags$br(),
                                              fluidRow(
                                                column(width=10, tags$p("Merge genotypic and phenotypic data tables and check if the merge is successful. A successful merge is contingent on proper match of Line IDs. 
                                        LineIDs with no matching IDs in either of the tables are printed in the message box.") 
                                                )),
                                              
                                              fluidRow(
                                                column(1),column(width = 4,tags$h5(tags$strong("Geno Data in Merged Table"))),
                                                column(6),column(width = 4, tags$h5(tags$strong("Pheno Data in Merged Table")))
                                              ),
                                              
                                              tags$br(),
                                              tags$head(
                                                tags$style(HTML("
                                          #messageMergedGenoME{

                                                    /*max-height: 1200px; Set maximum height */
                                                    overflow-y: scroll; /* Enable vertical scrolling */
                                                    overflow-x: scroll;  /*Hide horizontal scrolling */
                                                    overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                                    width: 350px;
                                                    /*max-width: 100%; */
                                                    padding: 6px 12px;
                                                    height: 150px;
                                                  }
                                                "))
                                              ),
                                              
                                              tags$head(
                                                tags$style(HTML("
                                          #messageMergeTrtME{

                                              /*max-height: 1200px; Set maximum height */
                                              overflow-y: scroll; /* Enable vertical scrolling */
                                              overflow-x: scroll;  /*Hide horizontal scrolling */
                                              overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                              width: 350px;
                                              /*max-width: 100%; */
                                              padding: 6px 12px;
                                              height: 150px;
                                             }
                                            "))
                                              ),
                                              tags$br(),
                                              
                                              fluidRow(
                                                column(width = 5,verbatimTextOutput("messageMergedGenoME")),column(6),column(7),column(width=5, verbatimTextOutput("messageMergeTrtME")),
                                                
                                              ),
                                              tags$br(),
                                              
                                              tags$head(
                                                tags$style(HTML("
                                              #MssgTrait {
                                               color: red;
                                              }")
                                                )
                                              ), 
                                              
                                              conditionalPanel(condition="input.traitME == 'NULL'",
                                                               fluidRow(
                                                                 column(1),column(width=8,tags$h5(tags$strong(textOutput("MssgTrait"))))),
                                              ),
                                              uiOutput("LocDistributionME")
                                              
                             )
                           )
                         )
                ),
                
                tabPanel("Enviromics Data Processing ",                  
                         
                         tabsetPanel(id="Enviromics Data",
                                     
                                     tabPanel("Get Enviromics Data",
                                              
                                              sidebarLayout(
                                                sidebarPanel(
                                                  
                                                  tags$h5(tags$strong("Upload Location Co-ordinates File (csv)")),
                                                  fileInput("inFileLocCoord","Location Coordinate File",accept = ".csv"),
                                                  tags$br(),
                                                  dateInput("startDate","Start Date for Daily Weather"),
                                                  tags$br(),
                                                  dateInput("endDate","End Date for Daily Weather"),
                                                  tags$br(),
                                                  tags$br(),
                                                  actionButton(inputId="getEnvData","Get Enviromics Data"),
                                                  tags$br(),
                                                  tags$br(),
                                                  
                                                ),
                                                mainPanel(
                                                  tags$strong(tags$h4("Enviromics Data Collection and Processing")),
                                                  fluidRow(
                                                    
                                                    column(width=10, tags$p("Extract weather data from NASA POWER database to explore weather patterns and estimate environmental relationship.
                                            The user needs to upload the co-ordinates file with the following columns 'Location','Country','Lat','Long','OtherLocName','LocName'. 
                                            The 'OtherLocName' and 'LocName' columns could contain abbreviations of location names.") 
                                                    )),
                                                  
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$br(),
                                                  
                                                  tags$head(
                                                    tags$style(HTML("
                                     #messageEnvDat{
                                          
                                          /*max-height: 1200px; Set maximum height */
                                          overflow-y: scroll; /* Enable vertical scrolling */
                                          overflow-x: scroll;  /*Hide horizontal scrolling */
                                          overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                          width: 600px; 
                                          /*max-width: 100%; */
                                          padding: 6px 12px;
                                          height: 300px;
                                         }
                                   "))
                                                  ),
                                                  
                                                  fluidRow(column(width = 4,verbatimTextOutput("messageEnvDat"))),
                                                  DTOutput("EnvTab")
                                                  #plotOutput("envPlot",height = "500px", width = "700px")
                                                  
                                                )
                                              )
                                     ),
                                     
                                     tabPanel("Get Enviromental Relationship",
                                              
                                              sidebarLayout(
                                                sidebarPanel( 
                                                  
                                                  tags$h5(tags$strong("Parameters for Environmental Relationship")),
                                                  checkboxInput("processWth","Process Raw Weather Data",value=FALSE),
                                                  checkboxInput("GaussKE","Use Gaussian Kernel",value=FALSE),
                                                  tags$br(),
                                                  actionButton(inputId="getEnvK","Get Env K"),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$h5(tags$strong(" Match Env Relationship and Pheno Data")),
                                                  checkboxInput("OtherLoc","Use Other Location Name",value=FALSE),
                                                ),
                                                mainPanel(
                                                  tags$strong(tags$h4("Estimate Enviromental Relationship")),
                                                  fluidRow(
                                                    
                                                    column(width=10, tags$p("Using the extracted weather data from NASA POWER database, estimate environmental relationship matrix using the kernel 
                                            method. The user needs to select the weather covariates and decide whether to use processed weather variables or not. The user also needs to set the kernel 
                                            for estimating environmental relationship.") 
                                                    )),
                                                  
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$br(),
                                                  
                                                  tags$head(
                                                    tags$style(HTML("
                                     #messageEnvK{
                                          
                                          /*max-height: 1200px; Set maximum height */
                                          overflow-y: scroll; /* Enable vertical scrolling */
                                          overflow-x: scroll;  /*Hide horizontal scrolling */
                                          overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                          width: 600px; 
                                          /*max-width: 100%; */
                                          padding: 6px 12px;
                                          height: 300px;
                                         }
                                   "))
                                                  ),
                                                  
                                                  fluidRow(column(width = 4,verbatimTextOutput("messageEnvK"))),
                                                  plotOutput("envPlot",height = "500px", width = "700px")
                                                  
                                                  
                                                )
                                              )
                                     )
                         )
                ),
                
                
                ### TS optimization tab
                tabPanel("Optimize Training Population",
                         sidebarLayout(
                           sidebarPanel(
                             numericInput(inputId="noCandidates","Candidate Set Size",value = 0,min =0, max=0),
                             tags$br(),
                             numericInput(inputId="noToSelect","Training Set Size",value =0,min =0, max=0),
                             tags$br(),
                             selectInput(inputId="optCriteria","Select Optimization Criteria",choices=c("PEVMEAN2","PEVMAX2","CDMEAN2","CDMAX2","DOpt","AOpt","EOpt"),multiple=TRUE,selected="PEVMEAN2"),
                             tags$br(),
                             actionButton(inputId="Optimize","Optimize Training Sets"),
                             tags$br(),
                             tags$br(),
                             tags$br(),
                             downloadButton("ExportOptTS", "Export Optimal TrainSet"),
                             tags$br(),
                             tags$br(),
                             tags$br(),
                             checkboxInput("setGA", "Use Default Genetic Algorithm Parameters", TRUE),
                             conditionalPanel(condition="input.setGA == false",
                                              numericInput(inputId="noPop","GA Population Size",value = 100,min =1,max=1000),
                                              
                                              numericInput(inputId="noElite","Elite Population Size",value =10,min =1, max=50),
                                              
                                              numericInput(inputId="mutProb","Mutation Probability",value =0.5,min =0.1, max=1),
                                              
                                              numericInput(inputId="mutInt","Mutation Intensity",value =1,min =0, max=1),
                                              
                                              numericInput(inputId="nIterGA","# Iterations",value =100,min =50, max=1000),
                                              
                                              selectInput(inputId="tabu","Tabu Search",choices=c("TRUE","FALSE"),selected="TRUE"),
                                              
                                              numericInput(inputId="tabumemsize","Tabu Memory Size",value =1,min =1, max=10),
                                              
                                              numericInput(inputId="lambda","Shrinkage Parameter (lambda)",value =1e-6,min =1e-9, max=1),
                                              
                                              numericInput(inputId="mcCores","#Cores",value =10,min =1, max=20)
                                              
                             )                          
                             
                           ), mainPanel(
                             fluidRow(
                               column(width=10,tags$h3(tags$strong("Select Optimal Training Population using 'STPGA' "))),
                               #actionButton("TS_Trait", "prev"),
                               #actionButton("TS_CVR", "next")
                             ),
                             fluidRow(
                               
                               column(width=10, tags$p("A training population for genomic prediction is selected by optimizing an objective function using a genetic algorithm.
                              The optimization method uses the 'GenAlgForSubsetSelection' function implemented in the 'STPGA' R package (Akdemir 2017).
                              In this implementation, the parameters for the genetic algorithm are set to default values.") 
                               )),
                             tags$br(),
                             tags$br(),
                             fluidRow(
                               
                               column(width=10, tags$p("The user needs to set the candidate and training population size. Candidate population refers to the set of genotypes with phenotypic data. The default size is set to the number of genotypes in the input genotypic data that are not in the target population set. 
                             Training population size refers to the number of genotypes that are to be selected in the training population. The default value is set to 80% of the candidate set.") 
                               )),
                             tags$br(),
                             tags$br(),
                             tags$head(
                               tags$style(HTML("
                                              #MssgOptTS {
                                                  color: red;
                                              }")
                               )
                             ), 
                             tags$head(
                               tags$style(HTML("
                                #MssgOptTSOK {
                                       color: green;
                                }")
                               )
                             ), 
                             
                             fluidRow(
                               column(1),column(width=8,tags$h5(tags$strong(textOutput("MssgOptTSOK"))))), 
                             
                             fluidRow(
                               column(1),column(width=8,tags$h5(tags$strong(textOutput("MssgOptTS"))))), 
                             
                             
                             tags$head(
                               tags$style(HTML("
                                #messageOptTS {
                                          
                                          /*max-height: 1200px; Set maximum height */
                                          overflow-y: scroll; /* Enable vertical scrolling */
                                          overflow-x: scroll;  /*Hide horizontal scrolling */
                                          overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                          width: 300px; 
                                          /*max-width: 100%; */
                                          padding: 6px 12px;
                                          height: 250px;
                                         }
                                        "))
                             ),
                             
                             fluidRow(column(2),column(width=6,tags$h6(verbatimTextOutput("messageOptTS")))),
                             tags$br(),
                             tags$br(),
                             fluidRow(
                               column(2),
                               column(6, div(
                                 tags$h5(tags$strong(textOutput("tsOptHeader"))),
                                 tags$h6(tableOutput("PredAccuracyTable")), 
                                 
                               )))
                           )
                         )), 
                
                
                ### CV tab
                tabPanel("Cross Validations",
                         tabsetPanel(id="CV",
                                     tabPanel("Single Trait",
                                              sidebarLayout(
                                                sidebarPanel(
                                                  numericInput(inputId="k",label = "Enter k for k-fold cross validation",value=2,min=2,max=10),
                                                  numericInput(inputId="nIter",label = "Enter n for n-iterations of each cycle of cross validation",value=2,min=2,max=10),
                                                  tags$br(),
                                                  tags$br(),
                                                  actionButton("CrossValidationST", "Run Cross Validation for Single Trait Models"),
                                                  tags$br(),
                                                  tags$br(),
                                                ), 
                                                mainPanel(
                                                  fluidRow(
                                                    column(1),column(width=7,tags$h3(tags$strong("Cross Validation of GP Models"))), 
                                                  ),
                                                  fluidRow(
                                                    column(width=10, tags$p("Perform an optional k-fold cross validation with training dataset to identify the model with the best prediction accuracy.
                             K-fold cross validation is performed using 'emCV' function implemented in the 'bWGR' package (Xavier et al. 2020)."))),
                                                  tags$br(),
                                                  fluidRow(
                                                    column(1),column(width=8,tags$h5(tags$strong(textOutput("cvrHeader"))))),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$head(
                                                    tags$style(HTML("
                                              #MssgST {
                                               color: red;
                                              }")
                                                    )
                                                  ),
                                                  tags$head(
                                                    tags$style(HTML("
                                              #MssgSTOK {
                                               color: green;
                                              }")
                                                    )
                                                  ),
                                                  fluidRow(
                                                    column(1),column(width=8,tags$h5(tags$strong(textOutput("MssgSTOK"))))),
                                                  fluidRow(
                                                    column(1),column(width=8,tags$h5(tags$strong(textOutput("MssgST"))))),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$head(
                                                    tags$style(HTML("
                                     #messageST5{
                                        
                                          /*max-height: 1200px; Set maximum height */
                                          overflow-y: scroll; /* Enable vertical scrolling */
                                          overflow-x: scroll;  /*Hide horizontal scrolling */
                                          overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                          width: 600px; 
                                          /*max-width: 100%; */
                                          padding: 6px 12px;
                                          height: 150px;
                                      }
                                   "))
                                                  ),
                                                  
                                                  fluidRow(column(width = 6,verbatimTextOutput("messageST5"))),
                                                  tags$br(),
                                                  fluidRow(
                                                    column(2),column(width=8,tags$h5(tableOutput("emCVRST")))),
                                                ),
                                              )),
                                     
                                     tabPanel("Multi-trait",
                                              sidebarLayout(
                                                sidebarPanel(
                                                  numericInput(inputId="k",label = "Enter k for k-fold cross validation",value=2,min=2,max=10),
                                                  numericInput(inputId="nIter",label = "Enter n for n-iterations of each cycle of cross validation",value=2,min=2,max=10),
                                                  tags$br(),
                                                  tags$br(),
                                                  actionButton("CrossValidationMT", "Run CV for MT Models"),
                                                  tags$br()
                                                ),
                                                mainPanel(
                                                  fluidRow(
                                                    column(1),column(width=7,tags$h3(tags$strong("Cross Validation of GP Models"))), 
                                                  ),
                                                  fluidRow(
                                                    column(width=10, tags$p("Perform an optional k-fold cross validation with training dataset to identify the model with the best prediction accuracy.
                                 For multi-trait models, 'Multitrait' and 'mmer' functions implemented in BGLR and Sommer packages are used for cross validations."))),
                                                  tags$br(),
                                                  fluidRow(
                                                    column(1),column(width=8,tags$h5(tags$strong(textOutput("cvrHeaderMT"))))),
                                                  tags$br(),
                                                  tags$head(
                                                    tags$style(HTML("
                                              #MssgMTOK {
                                               color: green;
                                              }")
                                                    )
                                                  ),
                                                  tags$head(
                                                    tags$style(HTML("
                                              #MssgMT {
                                               color: red;
                                              }")
                                                    )
                                                  ),
                                                  fluidRow(
                                                    column(1),column(width=8,tags$h5(tags$strong(textOutput("MssgMTOK"))))),
                                                  fluidRow(
                                                    column(1),column(width=8,tags$h5(tags$strong(textOutput("MssgMT"))))),
                                                  tags$br(),
                                                  
                                                  tags$head(
                                                    tags$style(HTML("
                                     #messageMT5{
                                        
                                          /*max-height: 1200px; Set maximum height */
                                          overflow-y: scroll; /* Enable vertical scrolling */
                                          overflow-x: scroll;  /*Hide horizontal scrolling */
                                          overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                          width: 600px; 
                                          /*max-width: 100%; */
                                          padding: 6px 12px;
                                          height: 150px;
                                      }
                                   "))
                                                  ),
                                                  
                                                  fluidRow(column(width = 6,verbatimTextOutput("messageMT5"))),
                                                  
                                                  tags$br(),
                                                  fluidRow(
                                                    column(3),column(width=8,tags$h5(tableOutput("emCVRMT")))),
                                                  tags$br()
                                                )
                                              )
                                     ),
                                     tabPanel("Multi Environment",
                                              sidebarLayout(
                                                
                                                sidebarPanel(
                                                  
                                                  tags$h4(tags$strong("Crossvalidation (CV) Method")),
                                                  tags$br(),
                                                  selectInput(inputId = "CVMet","Select CV Method",choices=c("CV1","CV2","CV0","CV00","CV_LOFO"),selected="CV1",multiple=FALSE),
                                                  
                                                  conditionalPanel(condition = "input.CVMet=='CV1' || input.CVMet=='CV2' || input.CVMet== 'CV0' || input.CVMet=='CV00'", 
                                                                   numericInput(inputId="kME",label = "Enter k for k-fold cross validation",value=2,min=2,max=10),
                                                                   numericInput(inputId="nIterME",label = "Enter n for n-iterations of each cycle of cross validation",value=2,min=2,max=10),
                                                                   tags$br(),
                                                                   tags$br(),
                                                                   checkboxInput("fitEnvCov", "Include Enviromics Kernel from Step 2", FALSE),
                                                                   
                                                  ),
                                                  conditionalPanel(condition = "input.CVMet=='CV_LOFO'",
                                                                   tags$h5(tags$strong("Subset Data")),
                                                                   selectInput(inputId="YearMECV","Select Year/Years",choices="All",selected="All",multiple=TRUE),
                                                                   tags$br(),
                                                                   selectInput(inputId="LocationMECV","Select Location/Locations",choices="All",selected="All",multiple=TRUE),
                                                                   tags$br(),
                                                                   tags$h5(tags$strong("Choose Covariates")),
                                                                   selectInput(inputId = "EnvVarIDCV", "Environmental Factor", choices = NULL, multiple = FALSE),
                                                                   tags$br(),
                                                                   selectInput(inputId = "fixedMECV", "Fixed Effect", choices = NULL, multiple = FALSE),
                                                                   tags$br(),
                                                                   tags$h5(tags$strong("Leave One Factor Level Out")),
                                                                   selectInput(inputId = "CVFactor","Select Factor for LOFO CV",choices=NULL,selected=NULL,multiple=FALSE),
                                                                   tags$br(),  
                                                                   tags$br(),
                                                                   checkboxInput("fitEnvCov", "Include Enviromics Kernel from Step 2", FALSE),
                                                  ),
                                                  tags$br(),
                                                  actionButton("CrossValidationME", "Run CV for ME GP Models"),
                                                  tags$br()
                                                ),
                                                mainPanel(
                                                  fluidRow(
                                                    column(1),column(width=7,tags$h3(tags$strong("Cross Validation of GP Models"))), 
                                                  ),
                                                  
                                                  fluidRow(
                                                    column(width=9, tags$p("Perform crossvalidation for multi-environmental models to identify the model 
                                         with the best prediction accuracy. Crossvalidation for multi-environmental data is performed using the BGGE/EnvRtype packages.")),
                                                    column(width=10,
                                                           tags$ul(
                                                             tags$li("CV1, where novel genotypes in tested environments are predicted."),
                                                             tags$li("CV2, where tested genotypes in tested environments are predicted."),
                                                             tags$li("CV0, where tested genotypes in untested novel environments are predicted."),
                                                             tags$li("CV00, where novel genotypes in novel environments are predicted."),
                                                             tags$li("CV LOFO (Leave One Factor Out), eg: Leave One Test Out/ Leave One Line Out cross validation."),
                                                           )
                                                    )
                                                  ),
                                                  tags$br(),
                                                  tags$head(
                                                    tags$style(HTML("
                                              #MssgMEOK {
                                               color: green;
                                              }")
                                                    )
                                                  ),
                                                  
                                                  tags$head(
                                                    tags$style(HTML("
                                              #MssgME {
                                               color: red;
                                              }")
                                                    )
                                                  ),
                                                  fluidRow(
                                                    column(1),column(width=8,tags$h5(tags$strong(textOutput("MssgMEOK"))))),
                                                  fluidRow(
                                                    column(1),column(width=8,tags$h5(tags$strong(textOutput("MssgME"))))),
                                                  
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$head(
                                                    tags$style(HTML("
                                     #messageME6{
                                        
                                          /*max-height: 1200px; Set maximum height */
                                          overflow-y: scroll; /* Enable vertical scrolling */
                                          overflow-x: scroll;  /*Hide horizontal scrolling */
                                          overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                          width: 600px; 
                                          /*max-width: 100%; */
                                          padding: 6px 12px;
                                          height: 150px;
                                      }
                                   "))
                                                  ),
                                                  
                                                  fluidRow(column(width = 6,verbatimTextOutput("messageME6"))),
                                                  # 
                                                  
                                                  tags$br(),
                                                  fluidRow(
                                                    column(3),column(width=8,tags$h5(tableOutput("emCVRME")))),
                                                  tags$br()
                                                )
                                              )
                                     )
                         )
                ), 
                ### GP Tab 
                tabPanel("Genomic Prediction",
                         tabsetPanel(id="GPModels",
                                     tabPanel("Single Trait",
                                              
                                              sidebarLayout(
                                                sidebarPanel(
                                                  selectInput(inputId="TrainSet","Choose Training Set",c("Complete Input Genotype Set","Optimal Train Set from Step 3","Random Train Set from Step 3"),selected="Complete Input Genotype Set"),
                                                  tags$br(),
                                                  tags$br(),
                                                  selectInput(inputId="fixed","Choose Fixed Effect",choices=NULL),
                                                  tags$br(),
                                                  tags$br(),
                                                  selectInput(inputId="GPModelST","Choose Prediction Model for Single Trait",c("rrBLUP (rrBLUP)","rrBLUP (bWGR)","BayesB (bWGR)","BayesLASSO (bWGR)")),
                                                  tags$br(),
                                                  tags$br(),
                                                  actionButton("RunPredictionsST", "Predict Single Trait!"),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$h5(tags$strong("View Predictions Table for Target/Training Set")),
                                                  selectInput(inputId="STGPOut","Choose Data Set ",c("Target","Training"),selected = "Target",multiple = FALSE),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$br(),
                                                  downloadButton("ExportOut", "Export Output Table")
                                                  
                                                ), mainPanel(
                                                  fluidRow(
                                                    column(width=8,tags$h3(tags$strong("Train Genomic Prediction Model"))), 
                                                    #actionButton("GP_CVR", "prev"),
                                                    # actionButton("GP_VP", "next")
                                                  ),
                                                  fluidRow(
                                                    column(width=10, tags$p(" Select the statistical method to train the genomic prediction model and predict the values of target lines. rrBLUP method is implemented using the
                                                    rrBLUP (Endelman 2011) package. Expectation maximization based RR-BLUP, BayesB and BayesLASSO methods are implemented using the bWGR package (Xavier et al. 2020).
                                                "))),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$head(
                                                    tags$style(HTML("
                                                   #MssgSTGP{
                                                    color: red;
                                                   }")
                                                    )
                                                  ),
                                                  tags$head(
                                                    tags$style(HTML("
                                                   #MssgSTGPOK{
                                                    color: green;
                                                   }")
                                                    )
                                                  ),
                                                  fluidRow(
                                                    column(2),column(width=8,tags$h5(tags$strong(textOutput("MssgSTGPOK"))))),
                                                  
                                                  fluidRow(
                                                    column(2),column(width=8,tags$h5(tags$strong(textOutput("MssgSTGP"))))),
                                                  tags$br(),
                                                  tags$br(),
                                                  
                                                  tags$head(
                                                    tags$style(HTML("
                                                      #messageSTGPRun{
                                        
                                                      /*max-height: 1200px; Set maximum height */
                                                      overflow-y: scroll; /* Enable vertical scrolling */
                                                      overflow-x: scroll;  /*Hide horizontal scrolling */
                                                      overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                                      width: 600px; 
                                                      /*max-width: 100%; */
                                                      padding: 6px 12px;
                                                      height: 150px;
                                                    }
                                                 "))
                                                  ),
                                                  
                                                  fluidRow(column(width = 6,verbatimTextOutput("messageSTGPRun"))),
                                                  tags$br(),
                                                  
                                                  tags$h4(textOutput("RankedLinesHeader")),
                                                  tableOutput("Ranked_Lines_for_SelectionST"),
                                                  tags$br(),
                                                  tags$br(),
                                                  
                                                  plotlyOutput("plots", height = "400px", width = "100%"),
                                                  #uiOutput("scatter")
                                                  
                                                  fluidRow(
                                                    lapply(1:10, function(i) {
                                                      plotlyOutput(outputId = paste0("plot", i), height = "400px", width = "100%")
                                                    })
                                                  )
                                                )
                                              )
                                     ),
                                     
                                     tabPanel("Multi-trait",
                                              
                                              sidebarLayout(
                                                sidebarPanel(
                                                  selectInput(inputId="TrainSetMT","Choose Training Set",c("Complete Input Genotype Set","Optimal Train Set from Step 3","Random Train Set from Step 3")),
                                                  tags$br(),
                                                  tags$br(),
                                                  selectInput(inputId="GPModelMT","Choose Prediction Model for Multiple Traits",c("BRR (BGLR)","RKHS (BGLR)","Spike-Slab(BGLR)")),
                                                  actionButton("RunPredictionsMT", "Predict Multiple Traits!"),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$h5(tags$strong("View Predictions Table for Target/Training Set")),
                                                  selectInput(inputId="MTGPOut","Choose Data Set ",c("Target","Training"),selected = "Target",multiple = FALSE),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$br(),
                                                  downloadButton("ExportOutMT", "Export Output Table")
                                                  
                                                ), mainPanel(
                                                  fluidRow(
                                                    column(width=8,tags$h3(tags$strong("Train Multi-trait Genomic Prediction Model"))), 
                                                    #actionButton("GP_ST", "prev"),
                                                    #actionButton("GP_MTME", "next")
                                                  ),
                                                  fluidRow(
                                                    column(width=10, tags$p("Select a statistical method to train the genomic prediction model and predict the values of target lines. Multi-trait predictions are implemented using the BGLR and Sommer packages. The Multitrait function in BGLR implements Bayesian Ridge Regression, RKHS, and Spike-Slab methods (Perez-Rodriguez P and de Los Campos G, 2022). Multi-trait GBLUP is implemented using the 'mmer' function in 'sommer' package (Covarrubias-Pazaran G 2016)."))),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$head(
                                                    tags$style(HTML("
                                                     #MssgMTGPOK {
                                                        color: green;
                                                     }")
                                                    )
                                                  ),
                                                  tags$head(
                                                    tags$style(HTML("
                                                     #MssgMTGP {
                                                        color: red;
                                                     }")
                                                    )
                                                  ),
                                                  fluidRow(
                                                    column(2),column(width=8,tags$h5(tags$strong(textOutput("MssgMTGPOK"))))),
                                                  fluidRow(
                                                    column(2),column(width=8,tags$h5(tags$strong(textOutput("MssgMTGP"))))),
                                                  tags$br(),
                                                  
                                                  tags$br(),
                                                  
                                                  tags$head(
                                                    tags$style(HTML("
                                                      #messageMTGPRun{
                                        
                                                      /*max-height: 1200px; Set maximum height */
                                                      overflow-y: scroll; /* Enable vertical scrolling */
                                                      overflow-x: scroll;  /*Hide horizontal scrolling */
                                                      overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                                      width: 600px; 
                                                      /*max-width: 100%; */
                                                      padding: 6px 12px;
                                                      height: 150px;
                                                    }
                                                 "))
                                                  ),
                                                  
                                                  fluidRow(column(width = 6,verbatimTextOutput("messageMTGPRun"))),
                                                  tags$br(),
                                                  tags$br(),
                                                  tableOutput("Ranked_Lines_for_SelectionMT"),
                                                  tags$br(),
                                                  
                                                  fluidRow(
                                                    lapply(1:10, function(i) {
                                                      plotlyOutput(outputId = paste0("plotMT", i), height = "400px", width = "100%")
                                                    })
                                                  )
                                                )
                                              )
                                     ),
                                     
                                     tabPanel("Multi-environment",
                                              sidebarLayout(
                                                sidebarPanel(
                                                  tags$h5(tags$strong("Subset Data")),
                                                  selectInput(inputId="YearME","Select Year/Years",choices="All",selected="All",multiple=TRUE),
                                                  tags$br(),
                                                  selectInput(inputId="LocationME","Select Location/Locations",choices="All",selected="All",multiple=TRUE),
                                                  tags$br(),
                                                  tags$h5(tags$strong("Choose Covariates")),
                                                  selectInput(inputId = "EnvVarID", "Environmental Factor", choices = NULL, multiple = FALSE),
                                                  tags$br(),
                                                  selectInput(inputId = "fixedME", "Fixed Effect", choices = NULL, multiple = FALSE),
                                                  tags$br(),
                                                  
                                                  # bsCollapse(
                                                  #   id = "collapseCov",
                                                  #   bsCollapsePanel(tags$h5(tags$strong("Choose Covariates")), 
                                                  #   #tags$h5(tags$strong("Choose Covariates")),
                                                  #     selectInput(inputId="EnvVarID","Environmental Factor",choices=NULL,multiple=FALSE),
                                                  #     tags$br(),
                                                  #     selectInput(inputId="fixedME","Fixed Effect",choices=NULL,multiple=FALSE),
                                                  #     tags$br()
                                                  #  ,style = "primary")),
                                                  # tags$br(),
                                                  
                                                  
                                                  # bsCollapse(
                                                  #   id = "collapseCov",
                                                  #   bsCollapsePanel(
                                                  #     
                                                  #     style = "primary"
                                                  #   )
                                                  # ),
                                                  
                                                  tags$br(),
                                                  selectInput(inputId="Package","Choose Package for MultiEnvironment GP Modeling ",c("BGGE-EnvRType")),
                                                  
                                                  conditionalPanel(condition="input.Package == 'BGGE-EnvRType'",
                                                                   
                                                                   selectInput(inputId="GKernelMet","Choose Genotype Kernel Method",c("Linear","Gaussian"),selected="Linear",multiple=FALSE),
                                                                   tags$br(),
                                                                   checkboxInput("fitEnvCov", "Include Enviromics Kernel from Step 2", FALSE),
                                                                   tags$br(),
                                                                   tags$br(),
                                                  ),
                                                  # conditionalPanel(condition="input.Package == 'SOMMER'",
                                                  #     selectInput(inputId="MEModelSom","Choose Model ",c("Main Effect","Homogeneous Variance (CS)","Heterogeneous Variance (CS+DG)","Unstructured (US)")),
                                                  #     tags$br(),
                                                  # ),
                                                  tags$br(),
                                                  tags$br(),
                                                  actionButton("RunPredictionsME", "Fit Multi-environmental Model!"),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$h5(tags$strong("View Output Table from Model")),
                                                  conditionalPanel(condition="input.fitEnvCov == false",
                                                                   selectInput(inputId="MEModelEnvR","Choose ME Model ",c("Main Effect (G+E)","Homogeneous Variance (G+E+GxE)","Heterogeneous Variance (G+E+GxEi)"),selected = "Main Effect (G+E)",multiple = FALSE)
                                                  ),
                                                  conditionalPanel(condition="input.fitEnvCov == true",
                                                                   selectInput(inputId="MEModelEnvR","Choose ME Model ",c("Main Effect (G+E+W)","Homogeneous Variance (G+E+GxE+W)","Heterogeneous Variance (G+E+GxEi+W)"),selected = "Main Effect (G+E+W)",multiple = FALSE)
                                                  ),
                                                  tags$br(),
                                                  tags$br(),
                                                  downloadButton("ExportOutME", "Export Predictions Table")
                                                ),
                                                mainPanel(
                                                  
                                                  fluidRow(
                                                    column(width=10, tags$p("Environments are defined as combinations of year and locations in which phenotypic data are collected.
                                                  Multi-environmental models are implemented using the EnvRType/BGGE pipeline as well as 'mmer' function in 'sommer' package.
                                                  Fit genomic prediction models taking into account only the main effects or the main effects + GxE effects.
                                                  The user needs to select one or many years and locations and the type of variance-covariance structure for fitting the model"))),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$head(
                                                    tags$style(HTML("
                                                     #MssgMEGPOK {
                                                        color: green;
                                                     }")
                                                    )
                                                  ),
                                                  tags$head(
                                                    tags$style(HTML("
                                                     #MssgMEGP {
                                                        color: red;
                                                     }")
                                                    )
                                                  ),
                                                  fluidRow(
                                                    column(2),column(width=8,tags$h5(tags$strong(textOutput("MssgMEGPOK"))))),
                                                  fluidRow(
                                                    column(2),column(width=8,tags$h5(tags$strong(textOutput("MssgMEGP"))))),
                                                  tags$br(),
                                                  tags$br(),
                                                  tags$head(
                                                    tags$style(HTML("
                                                      #messageMEGPRun{
                                        
                                                      /*max-height: 1200px; Set maximum height */
                                                      overflow-y: scroll; /* Enable vertical scrolling */
                                                      overflow-x: scroll;  /*Hide horizontal scrolling */
                                                      overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
                                                      width: 600px; 
                                                      /*max-width: 100%; */
                                                      padding: 6px 12px;
                                                      height: 150px;
                                                    }
                                                 "))
                                                  ),
                                                  
                                                  fluidRow(column(width = 6,verbatimTextOutput("messageMEGPRun"))),
                                                  tags$br(),
                                                  tags$br(),
                                                  tableOutput("Ranked_Lines_for_SelectionSTME"),
                                                  
                                                )
                                              )
                                     )
                                     
                         )
                )
    )
    
  )
 
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "GS4PB"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
