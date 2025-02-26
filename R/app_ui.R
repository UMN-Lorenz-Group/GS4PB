
#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinyFiles
#' @import shinydashboard
#' @import fs
#' @import htmlwidgets
#' @importFrom renv config
#' @importFrom shinyjs useShinyjs
#' @importFrom DT renderDT DTOutput
#' @importFrom plotly plotlyOutput
#' @importFrom shinyWidgets numericRangeInput
#' @importFrom shinyFiles shinyDirButton
#' @importFrom fs path
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd

###

app_ui <- function(request){


# Custom CSS to make primary box black

	custom_css <- tags$style(HTML("
	   
	   /* Center-align box titles */
	   .box.box-primary .box-header h3 {
			text-align: center;
			font-weight: bold;
			width: 100%;
	   }
	  
	   .box.box-primary {
		border-top-color: black !important;
	   }
	   
	  .box.box-primary > .box-header {
		background-color: black !important;
		color: white !important;
		border-bottom: 2px solid #4CAF50;
	  }
	  
	")) 

##
 


  dashboardPage(skin="black",
  dashboardHeader(title = "GS4PB App"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Start GS Pipeline", tabName = "gs_pipeline", icon = icon("play-circle")),
      
      menuItem("Genotypic Data Processing",tabName="genotypic_data",icon = icon("database"),
        menuSubItem("Load Genotypic Data",tabName="loadGeno",icon =icon("database")),
        menuSubItem("Filter Genotypic Data",tabName="filterGeno",icon=icon("filter")),
        menuSubItem("Impute Genotypic Data",tabName="imputeGeno",icon=icon("pen"))
      ),
      menuItem("Phenotypic Data Processing", tabName = "phenotypic_data", icon = icon("chart-bar")),
      menuItem("Merge Geno-Pheno Data", tabName = "geno_pheno_merge", icon = icon("project-diagram")),
     
      menuItem("Enviromics Data Processing", tabName = "enviromics", icon = icon("cloud"),
         menuSubItem("Get Enviromics Data",tabName="EnvrData",icon=icon("cloud")),
         menuSubItem("Get Environmental Relationship",tabName="EnvK",icon=icon("cloud"))
      ),
      
      menuItem("Optimize Training Population", tabName = "optimize_training", icon = icon("users")),
      
      menuItem("Cross Validations", tabName = "cross_validations", icon = icon("check-circle"),
        menuSubItem("Single Trait CV",tabName="STCV",icon=icon("check-circle")),
        menuSubItem("Multi-trait CV",tabName="MTCV",icon=icon("check-circle")), 
        menuSubItem("Multi-environment CV",tabName="MECV",icon=icon("check-circle"))
      ),
      
      menuItem("Genomic Prediction", tabName = "genomic_prediction", icon = icon("chart-line"), 
        menuSubItem("Single Trait GP",tabName="STGP",icon=icon("chart-line")),
        menuSubItem("Multi-trait GP",tabName="MTGP",icon=icon("chart-line")), 
        menuSubItem("Multi-environment GP",tabName="MEGP",icon=icon("chart-line"))
      )
    )
  ),
  
  dashboardBody(
    useShinyjs(),	
    custom_css,  # Include custom CSS
    tabItems(
      tabItem(
        tabName = "home",
        fluidRow(
          column(12, align = "center", tags$h2(tags$strong("GS4PB App")))),
        
  	    fluidRow(
		      box(width = 12, title = tags$strong("Genomic Selection For Plant Breeding Applications"), solidHeader = TRUE, status = "primary",
          tags$h4("The GS4PB App implements a genomic selection pipeline. 
                        The GS pipeline involves steps starting from data management to making selection decisions based on genome estimated breeding values.The steps 
                        in the pipeline are depicted in the schematic below.")
        )),
        tags$br(),
        fluidRow(
					   box(width = 12, title = "", solidHeader = TRUE, status = "primary",
                         tags$a(tags$img(src="www/GSPipeline.png",
                                         title="Genomic Selection Pipeline",
                                         class = "img-fluid",  # Bootstrap class for responsive images
                                         style = "max-width: 100%; height: auto;"
                         ))                              
                )
        ),
		    tags$br(),
		    
        fluidRow(
                column(12),align="left",tags$h5(("Contributors: Vishnu Ramasubramanian, Cleiton Wartha, Paolo Vitale, Sushan Ru,and Aaron Lorenz"))),
        tags$br(),

		    fluidRow(
            column(12),align="left",tags$h5("Contact: Vishnu Ramasubramanian - vramasub@umn.edu"))
		    
      ),
      
      tabItem(
        tabName = "gs_pipeline",
        fluidRow(
          column(4,
            box(width = 12, title = "Set output dir to save log and result files", solidHeader = TRUE, status = "primary",
                tags$br(),
                # Button to trigger directory selection
                shinyDirButton("directory", "Select Output Directory", "Please select a folder"),
                tags$br(),
                # Display selected directory path
                verbatimTextOutput("dir_path")
            )
           ),
          column(8,
         
           box(width = 12, title = "GS Pipeline Instructions", solidHeader = TRUE, status = "primary",
            tags$h4("The GS pipeline involves steps ranging from genotypic and phenotypic data processing, merging geno and pheno data, optimizing training populations, 
                    cross-validation of genomic prediction models, and genomic prediction."),
            tags$img(src="www/GSPipelineInstructions.png", 
                     class = "img-fluid",  # Bootstrap class for responsive images
                     style = "max-width: 100%; height: auto;")
                     
           )
          )
        )
       ),
      
      tabItem(
        tabName = "loadGeno",
        fluidRow(
          box(width=12,title="Load Training and Target Data",solidHeader = TRUE,status="primary",
          tags$h4("Upload genotypic data in VCF file format. Genotypic data could include the training and target data in the same file 
                              or just the training data alone.If your data includes a target population, please upload the line IDs of target population in the right panel. 
                              If not leave it blank. Once you upload the target line IDs in csv format, 
                              select tehe column containing the target IDs")
          )
        ),
        tags$br(),
        fluidRow(
          box(width = 6, title = "Load Genotypic Data", solidHeader = TRUE, status = "primary",
              fileInput("infileVCF", "Upload Genotypic Data (VCF Format)", accept = ".vcf"),
              tags$br(),
              verbatimTextOutput("messageGenoStats"),
              tags$br(),
          ),
          box(width = 6, title = "Load Target Info", solidHeader = TRUE, status = "primary",
              fileInput("infileTargetTable", "Upload Target Line IDs (CSV)", accept = ".csv"),
              selectInput(inputId="TargetIDCol","Select Target ID Col",choices=NULL,multiple=FALSE),
              verbatimTextOutput("messageTargetStats")
          ),
        )    
      ),
    
      tabItem(
        tabName = "filterGeno",
        fluidRow(
          box(width=12,title="Filter Genotypic Data",solidHeader = TRUE,status="primary",
               tags$h4("Filter sites (markers) and taxa (lines) in genotype table using rTASSEL (Monier et al. 2021), if you uploaded a raw genotype file or a QC genotype file that requires filtering. For example, it is common to 
                    set minimum number of sites that are not missing to 80% of the lines in the input genotype data and set the MAF threshold to 0.02/0.05. 
                    If you uploaded a QC genotype file, you can choose to skip this step")
             )
        ),
        tags$br(),
       
        fluidRow(
          
              box(width=6,title= tags$strong(tags$h4("Filter Genotype Table Sites")),solidHeader=TRUE,status="primary",
                  numericInput(inputId="siteMinCnt","Minimum Site Count (TASSEL)",value = 0,min =0, max=0),
                  numericInput(inputId="MAF","Minimum Allele Frequency",value =0.02,min =0, max=0.5),
                  actionButton(inputId="FilterSites","Filter Genotype Table Sites"),
                  checkboxInput("setGenoFilt1Tas", "Use Filtered Genotypes From Filter Sites For Next Steps", TRUE)
              ),
          
              box(width=6,title= tags$strong(tags$h4("Filter Genotype Table Taxa")),solidHeader=TRUE,status="primary",
                 numericInput(inputId="minNotMissing","Minimum Not Missing (TASSEL)",value = 0.9,min =0.5, max=1),
                 tags$br(),
                 tags$br(),
                 tags$br(),
                 actionButton(inputId="FilterTaxa","Filter Genotype Table Taxa"),
                 checkboxInput("setGenoFilt2Tas", "Use Filtered Genotypes From Filter Taxa For Next Steps", TRUE)
            )
        ),
     
        box(width=6,title="Filtered Geno Sites",solidHeader = TRUE,status = "info",
         tags$h6(verbatimTextOutput("messageGenoFilt1"))
        ),
        box(width=6,title="Filtered Geno Taxa",solidHeader = TRUE,status = "info",
         tags$h6(verbatimTextOutput("messageGenoFilt2"))
        ),
       
        tags$br(),
        tags$br(),
        
      ),
      
      tabItem(
        tabName = "imputeGeno",
        
        fluidRow(
          column(4, 
                 box(title = "Select Imputation Method", width = 12,solidHeader = TRUE,status="primary",
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
                )
          ),
          column(8, 
                     conditionalPanel(condition = "input.imputeMet == 'LDKNNI' || input.imputeMet == 'Numeric'",
                                      
                                    box(width = 12,title = "Impute Genotype Table",solidHeader = TRUE,status="primary",      
                                     tags$h4(
                                          "Impute missing genotype scores using rTASSEL (Monier et al. 2021), if you uploaded a raw genotype file or a QC genotype file with missing scores. 
                                           Available options include numeric imputation and imputation using LD-K-Nearest Neighbors method are availble. 
                                           For LD-KNN imputation set parameters l and K (Money et al. 2015). l corresponds to the number of high LD sites  
                                           and k corresponds to the number of neighboring samples that are used to impute scores.")
                                    ),
                                      
                                   
                                    tags$br(),
                                    box(width=12,title="",solidHeader=TRUE,status="info",
                                           tags$h6(verbatimTextOutput("messageImpGeno1"))
                                    ),
                                    tags$br(),
                                    
                     ),
                     conditionalPanel(condition = "input.imputeMet == 'AlphaPlantImpute'",
                                      uiOutput("APIUI"),
                     )
            )
        ) 
     ),
  
		tabItem(
		  tabName = "phenotypic_data",
		        fluidRow(
		          column(4,
		          box(
		            title = "Data Input Options", width = 12, solidHeader = TRUE, status = "primary",
		            checkboxInput("chkPhenoME", "Load Pheno Data from Multi-Environs", FALSE),
		            br(),
		            conditionalPanel(condition = "input.chkPhenoME == false",
		                             selectInput("IDColSE", "Choose Uniq ID Col", choices = NULL, multiple = FALSE),
		                             selectInput("strainSE", "Choose Strain Col", choices = NULL, multiple = FALSE),
		                             numericRangeInput("traitColsNum", "Choose all the trait Cols", value = NULL, min = 1, max = 1),
		                             br()
		            ),
		            conditionalPanel(condition = "input.chkPhenoME == true",
		                             selectInput("traitCols", "Choose all the trait Cols", choices = NULL, multiple = TRUE),
		                             selectInput("IDColME", "Choose Uniq ID Col", choices = NULL, multiple = FALSE),
		                             selectInput("strainME", "Choose Strain Col", choices = NULL, multiple = FALSE),
		                             br(),
		                             tags$br(),
		                    )
		           )
		          ),
              column(8,
		           box(
		            title = "Phenotypic Data Processing", width = 12, solidHeader = TRUE, status = "primary",
		            tags$h4("Load phenotypic data and select traits for cross-validation and genomic prediction. Check the box to load data from multi-environmental trials."),
		            tags$br(),
		            conditionalPanel(condition = "input.chkPhenoME == false",
		                tags$head(
		                  tags$style(HTML("
                             .info-button {
                             font-size: 10px;         /* Smaller font size for superscript effect */
                             margin-left: 5px;        /* Space from the text */
                             padding: 2px 5px;        /* Adjust button size */
                            vertical-align: super;   /* Align as superscript */
                            }
                    "))),
		                
		                fluidRow(
		                    box(
		                      fluidRow(
		                        tags$label("Choose Phenotype File (CSV)"),
		                        actionButton("iPh", "i", class = "info-button"),  # Small "i" as superscript
		                        fileInput("infileBLUEsSE", label = NULL, accept = ".csv")
		                      ),
		                      checkboxInput("header", "Header", TRUE),
		                      tags$br(),
		                      verbatimTextOutput("messagePhSE")
		                    ), 
		                        # fluidRow(fileInput("infileBLUEsSE", "Choose Phenotype File (CSV)", accept = ".csv"),
		                        #  actionButton("iPh", "i")),
		                    box(
		                       selectInput("trait", "Choose One or More Traits", choices = NULL, multiple = TRUE),
		                       tags$br(),
		                       tags$br(),
		                       verbatimTextOutput("messageTrtSE")
		                    )
		                 )
		            ),
		            conditionalPanel(condition="input.chkPhenoME == true",
		                   tags$br(),
		                   
		                   tags$head(
		                     tags$style(HTML("
                             .info-button {
                              font-size: 10px;         /* Smaller font size for superscript effect */
                              margin-left: 5px;        /* Space from the text */
                              padding: 2px 5px;        /* Adjust button size */
                              vertical-align: super;   /* Align as superscript */
                             }
                        "))
		                   ),
		                   fluidRow(
		                     box(width = 6,title = "Phenotypic Data", 
		                     fluidRow(
		                           tags$label("Choose Phenotype File (CSV)"),
		                           actionButton("iPhME", "i", class = "info-button"),  # Small "i" as superscript
		                           fileInput("infileBLUEsME", label = NULL, accept = ".csv")
		                      ),
		                      checkboxInput("headerME", "Header", TRUE),
		                      tags$br(),
		                      verbatimTextOutput("messagePhME")
		                     ),
		                      
		                      # fileInput("infileBLUEsME", "Choose Phenotype File (CSV)", accept = ".csv"),
		                      # actionButton("iPhME", "i"),
		                     box(width=6,title="Select trait(s)",
		                      selectInput(inputId="traitME","Choose One or More Traits",choices=NULL,multiple=TRUE),
		                      tags$br(),
		                      tags$br(),
		                      tags$br(),
		                      tags$br(),
		                      verbatimTextOutput("messageTrtME")
		                     )
		                    )
		                )
		              )        
		            )
          )
     ),
 
    tabItem(
      tabName = "geno_pheno_merge",
        fluidRow(
          column(3, 
                 box(title = "Merge Genotypic and Phenotypic Data", width = 12,solidHeader = TRUE, status = "primary",
                     tags$br(),
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
                 )     
        
           ),
          
           column(8, 
                conditionalPanel(condition="input.chkPhenoME == false",
                  box(title = "Merge Genotypic and Phenotypic Data", width = 12,solidHeader = TRUE, status = "primary",
                     
                     tags$h4("Merge genotypic and phenotypic data tables and check if the merge is successful. A successful merge is contingent on proper match of Line IDs. 
                              LineIDs with no matching IDs in either of the tables are printed in the message box."),
                      
                     tags$br(), 
                     fluidRow(
                      box(width=6,title=tags$h5(tags$strong("Geno Data in Merged Table")),solidHeader=T,status="primary",
                                verbatimTextOutput("messageMergedGeno")
                          ),
                                         
                      box(width=6,title=tags$h5(tags$strong("Pheno Data in Merged Table")),solidHeader=T,status="primary",
                                          verbatimTextOutput("messageMergeTrtSE")
                      )
                     ),                   
                  
                  tags$br(),
                              
                  tags$head(
                            tags$style(HTML("
                                #MssgMergedData {
                                color: red;
                                }")
                             )
                  ),
                        
                  fluidRow(
                       tags$h5(tags$strong(textOutput("MssgMergedData")))
                      )
                  )
                 ),
                  
                 conditionalPanel(condition="input.chkPhenoME == true",
                                       
                          box(width=12,title="Merge Geno and Multi-environmental Pheno Data",solidHeader = T,status="primary",
                                           tags$h4("Merge genotypic and phenotypic data tables and check if the merge is successful. A successful merge is contingent on proper match of Line IDs. 
                                            LineIDs with no matching IDs in either of the tables are printed in the message box."),
                                   
                            fluidRow(
                                       box(width=6, title=tags$h5(tags$strong("Geno Data in Merged Table")),solidHeader=T,status="info",
                                            verbatimTextOutput("messageMergedGenoME") 
                                       ),
                                       box(width=6,title= tags$h5(tags$strong("Pheno Data in Merged Table")),solidHeader=T,status="info",
                                            verbatimTextOutput("messageMergeTrtME")
                                       )
                            ),
                                       
                            tags$br(),
                                       
                            tags$head(
                                         tags$style(HTML("
                                              #MssgTrait {
                                               color: red;
                                              }")
                                         )
                                      ), 
                           ),
                          
                          conditionalPanel(condition="input.traitME == 'NULL'",
                                                        fluidRow(
                                                          column(1),column(width=6,tags$h5(tags$strong(textOutput("MssgTrait"))))),
                          ),
                                      
                          box(width=12,title="Distribution of Phenotypic Values",solidHeader=TRUE, status = "primary",
                                         uiOutput("LocDistributionME")
                          )
                     )
              )
       )
      ),
      
      tabItem(
        tabName = "EnvrData",
        
        fluidRow(
          column(4,
            box(width = 12, title = "Enviromics Data Collection", solidHeader = TRUE, status = "primary",
              fileInput("inFileLocCoord", "Upload Location Coordinates File (CSV)", accept = ".csv"),
              dateInput("startDate", "Start Date for Weather Data"),
              dateInput("endDate", "End Date for Weather Data"),
              actionButton("getEnvData", "Get Enviromics Data"),
              verbatimTextOutput("messageEnvDat")
            )
          ),
      
       column(8,
           box(width=12,title= tags$strong(tags$h4("Enviromics Data Collection and Processing")),solidHeader=TRUE,status="primary",
            tags$h4("Extract weather data from NASA POWER database to explore weather patterns and estimate environmental relationship.
                                            The user needs to upload the co-ordinates file with the following columns 'Location','Country','Lat','Long','OtherLocName','LocName'. 
                                            The 'OtherLocName' and 'LocName' columns could contain abbreviations of location names.") 
           ),
           tags$br(),
           box(width=12,title="",solidHeader = T,status = "primary",collapsible=TRUE,
             verbatimTextOutput("messageEnvDat")
           ),
           
           tags$br(),
             tabBox(
               title = "Env Data",
               side = "left",
               width = 12,  # Adjust as needed
               tabPanel(
                 "Table View",
                 div(
                   DTOutput("EnvTab"),
                   style = "overflow-x: auto; width: 100%;"  # Prevents horizontal scrolling issues
                 )
               )
             )
         )
       )
      ),
		
		  tabItem(
		    tabName="EnvK",title = "Enviromental Relationship", solidHeader = TRUE, status = "primary",
		    
		    fluidRow( 
		      
		      column(4, 
		          box(width=12,title = "Enviromental Relationship", solidHeader = TRUE, status = "primary",
		            tags$h5(tags$strong("Parameters for Environmental Relationship")),
		            checkboxInput("processWth","Process Raw Weather Data",value=FALSE),
		            checkboxInput("GaussKE","Use Gaussian Kernel",value=FALSE),
		            tags$br(),
		            actionButton(inputId="getEnvK","Get Env K"),
		            tags$br(),
		            tags$br(),
		            tags$h5(tags$strong(" Match Env Relationship and Pheno Data")),
		            checkboxInput("OtherLoc","Use Other Location Name",value=FALSE),
		           )
		      ),
		      
		      column(8, 
		            box(width=12,title="Estimate Environmental Relationship",solidHeader=TRUE,status="primary",
		                
		               tags$h4("Using the extracted weather data from NASA POWER database, estimate environmental relationship matrix using the kernel 
                                            method. The user needs to select the weather covariates and decide whether to use processed weather variables or not. The user also needs to set the kernel 
                                            for estimating environmental relationship.") 
		              ),
		                
		                tags$br(),
		                tags$br(),
		                
		             box(width = 12,title="",solidHeader=T,status="info",verbatimTextOutput("messageEnvK")),
		             
		             tags$br(),
		             box(width=12,title="Environmental Relationship",solidHeader=TRUE,status="primary",collapsible=TRUE,
		               div(  
		                  plotOutput("envPlot",height = "500px", width = "100%"), 
		                  style= "display: flex; justify-content: center;"
		               )
		             )
		                
		        )
		      )
		  ),
      
      tabItem(
        tabName = "optimize_training",
         fluidRow(
              column(4,
                  box(width = 12, title = "Optimize TP Parameters", solidHeader = TRUE, status = "primary",
                  
                  numericInput(inputId="noCandidates","Candidate Population Size",value = 0,min =0, max=0),
                  numericInput(inputId="noToSelect","Training Population Size",value =0,min =0, max=0),
                  selectInput(inputId="optCriteria","Select Optimization Criteria",choices=c("PEVMEAN2","PEVMAX2","CDMEAN2","CDMAX2","DOpt","AOpt","EOpt"),multiple=TRUE,selected="PEVMEAN2"),
                  tags$br(),
                  actionButton(inputId="Optimize","Optimize Training Population"),
                  tags$br(),
                 
                  downloadButton("ExportOptTS", "Export Optimal Training Population"),
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
            )
        ),
        column(8,
                 
            box(width = 12, title = "Optimize Training Population using 'STPGA'", solidHeader = TRUE, status = "primary",
      
                   tags$h4("A training population for genomic prediction is selected by optimizing an objective function using a genetic algorithm.
                              The optimization method uses the 'GenAlgForSubsetSelection' function implemented in the 'STPGA' R package (Akdemir 2017).
                              In this implementation, the parameters for the genetic algorithm are set to default values."
                          ), 
                    
                      tags$br(),
                      tags$br(),
                      tags$h4("The user needs to set the candidate and training population size. Candidate population refers to the set of genotypes with phenotypic data. The default size is set to the number of genotypes in the input genotypic data that are not in the target population set. 
                             Training population size refers to the number of genotypes that are to be selected in the training population. The default value is set to 80% of the candidate set."
                              ), 
                      
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
                      
                      
                      tags$h5(tags$strong(textOutput("MssgOptTSOK"))), 
                                            
                      tags$h5(tags$strong(textOutput("MssgOptTS"))), 
                      
                      fluidRow(box(width=6,solidHeader = T,status="info",
                                   tags$h6(verbatimTextOutput("messageOptTS")))
                      ),
            ),
                      tags$br(),
                      tags$br(),
            
            fluidRow(
                    tabBox(
                      title = "Opt TS Out",
                      side = "left",
                      width = 12,  # Adjust as needed
                      tabPanel(
                        "Optimal TS Table",
                        div(
                          tableOutput("PredAccuracyTable"),
                          style = "overflow-x: auto; width: 100%;"  # Prevents horizontal scrolling issues
                        )
                      )
                    )
              )
            )
          )
     ),
	
		tabItem(
		    tabName="STCV",
		  
		    fluidRow(
		      column(4, 
		        box(width=12,title="Single Trait CV parameters",solidHeader = T,status="primary",
		            numericInput(inputId="k",label = "Enter k for k-fold cross validation",value=2,min=2,max=10),
		            numericInput(inputId="nIter",label = "Enter n for n-iterations of each cycle of cross validation",value=2,min=2,max=10),
		            tags$br(),
		            tags$br(),
		            actionButton("CrossValidationST", "Run Cross Validation for Single Trait Models"),
		            tags$br(),
		        )  
		     ),
		     column(8,
		        box(width=12,title="Single Trait Cross Validation",solidHeader = T,status="primary",
		            
		            tags$h4("Perform a k-fold cross validation with training dataset to identify the model with the best prediction accuracy.
                             K-fold cross validation is performed using 'emCV' function implemented in the 'bWGR' package (Xavier et al. 2020)."),
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
		            
		            fluidRow(
		              box(width=10,title="",solidHeader=T,status="info",
		               verbatimTextOutput("messageST5")
		              )
		            ),
		        ),
		            tags$br(),
		            tabBox(
		              title = "Single Trait CV ",
		              side = "left",
		              width = 12,  # Adjust as needed
		              tabPanel(
		                "GP Output Table",
		                div(
		                  tags$h5(tableOutput("emCVRST")),
		                  style = "overflow-x: auto; width: 100%;"  # Prevents horizontal scrolling issues
		                )
		              )
		            )
		        )
		     )
		),
		
### MTCV		
		  tabItem(
		   tabName="MTCV",
		   fluidRow(
		     column(4, 
		            box(width=12,title="Multitrait CV Parameters",solidHeader=TRUE,status="primary",
  		            numericInput(inputId="k",label = "Enter k for k-fold cross validation",value=2,min=2,max=10),
  		            numericInput(inputId="nIter",label = "Enter n for n-iterations of each cycle of cross validation",value=2,min=2,max=10),
  		            tags$br(),
  		            tags$br(),
  		            actionButton("CrossValidationMT", "Run CV for MT Models"),
  		            tags$br()  
		            )
		     ),
		     column(8, 
		            box(width=12,title="Cross Validation of Multitrait GP Models",solidHeader=TRUE,status="primary",
		                fluidRow(
		                  column(width=10, tags$h4("Perform a k-fold cross validation with training dataset to identify the model with the best prediction accuracy.
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

		                fluidRow(
		                  box(width=10,title="",solidHeader=T,status="info",collapsible=T,
		                        verbatimTextOutput("messageMT5")
		                  )
		                ),
		            ),   
		                tags$br(),
		                
		                tabBox(
		                  title = "Multitrait CV ",
		                  side = "left",
		                  width = 12,  # Adjust as needed
		                  tabPanel(
		                    "Table View",
		                    div(
		                      tags$h5(tableOutput("emCVRMT")),
		                      style = "overflow-x: auto; width: 100%;"  # Prevents horizontal scrolling issues
		                    )
		                  )
		                 )
		            
		     )
		    )
		 ),
		 
## ME-CV
		 tabItem(
		    tabName="MECV",
		    fluidRow(
		      column(4, 
		             box(width=12,title="Multi-environment CV Parameters",solidHeader=TRUE,status="primary",
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
		                 
		             )     
		    
		  ),
      column(8, 
             box(width=12,title="Cross Validation of Multi-environment GP Models",solidHeader=TRUE,status="primary",
             
               tags$h4("Perform crossvalidation for multi-environmental models to identify the model 
                         with the best prediction accuracy. Crossvalidation for multi-environmental data is performed using the BGGE/EnvRtype packages."
               ),
                
               column(width=10,
                          tags$ul(
                            tags$li("CV1, where novel genotypes in tested environments are predicted."),
                            tags$li("CV2, where tested genotypes in tested environments are predicted."),
                            tags$li("CV0, where tested genotypes in untested novel environments are predicted."),
                            tags$li("CV00, where novel genotypes in novel environments are predicted."),
                            tags$li("CV LOFO (Leave One Factor Out), eg: Leave One Test Out/ Leave One Line Out cross validation."),
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
                 
                 fluidRow(
                   box(width=10,title="",solidHeader = T,status="info",collapsible=T,
                       verbatimTextOutput("messageME6"))
                 ),
              ),
                 tags$br(),
               
               tabBox(
                   title = "Multi-environmental CV ",
                   side = "left",
                   width = 12,  # Adjust as needed
                   tabPanel(
                     "Table View",
                     div(
                       tags$h5(tableOutput("emCVRME")),
                       style = "overflow-x: auto; width: 100%;"  # Prevents horizontal scrolling issues
                     )
                   )
                 )
              )
         )   
    ),
		  
#### ST GP
     
		 tabItem(
		  tabName="STGP",
		   
		   fluidRow(
		     column(4,
		      box(width=12,title="Single Trait GP Parameters",solidHeader=T,status="primary", 
		          selectInput(inputId="TrainSet","Choose Training Set",c("Complete Input Genotype Set","Optimal Train Set from Step 3","Random Train Set from Step 3"),selected="Complete Input Genotype Set"),
		          selectInput(inputId="fixed","Choose Fixed Effect",choices=NULL),
		          selectInput(inputId="GPModelST","Choose Prediction Model for Single Trait",c("rrBLUP (rrBLUP)","rrBLUP (bWGR)","BayesB (bWGR)","BayesLASSO (bWGR)")),
		          tags$br(),
		          tags$br(),
		          actionButton("RunPredictionsST", "Predict Single Trait!"),
		          tags$br(),
		          tags$br(),
		          tags$h5(tags$strong("View Predictions Table for Target/Training Set")),
		          selectInput(inputId="STGPOut","Choose Data Set ",c("Target","Training"),selected = "Target",multiple = FALSE),
		          tags$br(),
		          tags$br(),
		          downloadButton("ExportOut", "Export Output Table")  
		      )
		     ),
		     column(8,
		       box(width=12,title="Single Trait Genomic Prediction",solidHeader=T,status="primary",
		           tags$h4(" Select the statistical method to train the genomic prediction model and predict the values of target lines. rrBLUP method is implemented using the
                                                    rrBLUP (Endelman 2011) package. Expectation maximization based RR-BLUP, BayesB and BayesLASSO methods are implemented using the bWGR package (Xavier et al. 2020).
                                                "),
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
		         
		           
		           fluidRow(
		         
		             box(width=12,title="",solidHeader = TRUE,status="info",collapsible=TRUE,
		               verbatimTextOutput("messageSTGPRun")
		           )),
		           
		        ),
		           tags$br(),
		           tags$h4(textOutput("RankedLinesHeader")),
		           
		           tags$br(),
		           tabBox(
		             title = "Single Trait GP Results",
		             side = "left",
		             width = 12,  # Adjust as needed
		             tabPanel(
		               "Output Table",
		               div(
		                 DTOutput("Ranked_Lines_for_SelectionST"),
		                 style = "overflow-x: auto; width: 100%;"  # Prevents horizontal scrolling issues
		               )
		             ),
		             tabPanel(
		               "Output Plot",
		                uiOutput("STPlotsUI")
		               
		             )
		           ),
		           tags$br(),
		           
#### 
		       )  
		    )
		),
		
### MT-GP Panel
		tabItem(
		 tabName="MTGP",
		 fluidRow(
		  column(4,
		        box(width=12,title="Multitrait GP Parameters",solidHeader=T,status="primary",
  		        selectInput(inputId="TrainSetMT","Choose Training Set",c("Complete Input Genotype Set","Optimal Train Set from Step 3","Random Train Set from Step 3")),
  		        tags$br(),
  		        tags$br(),
  		        selectInput(inputId="GPModelMT","Choose Prediction Model for Multiple Traits",c("BRR (BGLR)","RKHS (BGLR)","Spike-Slab(BGLR)","Mmer (Sommer)")),
  		        actionButton("RunPredictionsMT", "Predict Multiple Traits!"),
  		        tags$br(),
  		        tags$br(),
  		        tags$h5(tags$strong("View Predictions Table for Target/Training Set")),
  		        selectInput(inputId="MTGPOut","Choose Data Set ",c("Target","Training"),selected = "Target",multiple = FALSE),
  		        tags$br(),
  		        tags$br(),
  		        downloadButton("ExportOutMT", "Export Output Table")
		        )
		  ),
		  column(8,
		       box(width=12,title="Multitrait Genomic Prediction",solidHeader=T,status="primary",
		           tags$h4("Select a statistical method to train the genomic prediction model 
		                      and predict the values of target lines. Multi-trait predictions are 
		                      implemented using the BGLR and Sommer packages. The Multitrait function in BGLR 
		                      implements Bayesian Ridge Regression, RKHS, and Spike-slab methods (Perez-Rodriguez P and de Los Campos G, 2022). Multi-trait GBLUP is implemented using the 'mmer' 
		                      function in 'sommer' package (Covarrubias-Pazaran G 2016)."),
		      
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
		         
		       fluidRow(
		               box(width=12,title="",solidHeader = FALSE,status="info",collapsible=TRUE,
		                div(
  		                tags$h5(tags$strong(textOutput("MssgMTGPOK"))),
  		                tags$h5(tags$strong(textOutput("MssgMTGP")))
		                  ,style="overflow-wrap: anywhere;width: 100%;"
		                ),
		                tags$br(),
		                verbatimTextOutput("messageMTGPRun")
		              )
		       ),
		      ), 
		       tags$br(),
      	
		       tabBox(
      		         title = "Multitrait GP Results",
      		         side = "left",
      		         width = 12,  # Adjust as needed
      		         tabPanel(
      		           "GP Output Table",
      		           div(
      		             DTOutput("Ranked_Lines_for_SelectionMT"),
      		             style = "overflow-x: auto; width: 100%;"  # Prevents horizontal scrolling issues
      		           )
      		         ),
      		         tabPanel(
      		           "GP Output Plots",
      		           uiOutput("MTPlotsUI")
      		         )
      		 )
		    )
		  )
		),
		
### ME GP Panel

	 tabItem(
		  tabName="MEGP",
		  fluidRow(
		    column(4,
		           box(width=12,title="Multi-environmental GP Parameters",solidHeader=T,status="primary",
		               
		            selectInput(inputId="Package","Choose Package for MultiEnvironment GP Modeling ",c("BGGE-EnvRType")),
		            conditionalPanel(condition="input.Package == 'BGGE-EnvRType'",
		                  selectInput(inputId="GKernelMet","Choose Genotype Kernel Method",c("Linear","Gaussian"),selected="Linear",multiple=FALSE),
		                  tags$br(),
		                  checkboxInput("fitEnvCov", "Include Enviromics Kernel from Step 2", FALSE),
		                  
		            ),
		            tags$br(),
		           ),
		           box(width=12,title="Subset Data",solidHeader=T,status="primary",collapsible = TRUE,collapsed = TRUE,
		             selectInput(inputId="YearME","Select Year/Years",choices="All",selected="All",multiple=TRUE),
		             selectInput(inputId="LocationME","Select Location/Locations",choices="All",selected="All",multiple=TRUE),
		           ),
		           tags$br(),
		           box(width=12,title="Choose Covariates",solidHeader=T,status="primary",collapsible = TRUE,collapsed = TRUE,
		               selectInput(inputId = "EnvVarID", "Environmental Factor", choices = NULL, multiple = FALSE),
		               selectInput(inputId = "fixedME", "Fixed Effect", choices = NULL, multiple = FALSE),
		           ),
		           tags$br(),
		           box(width=12,title="Fit ME Model",solidHeader=T,status="primary",collapsible = FALSE,collapsed = FALSE,
		              actionButton("RunPredictionsME", "Fit Multi-environmental Model!"),
		           tags$br(),
		           ),
		           box(width=12,title="View Output Table",solidHeader=T,status="primary",collapsible = TRUE,collapsed = TRUE,
		               #tags$h5(tags$strong("View Output Table from Model")),
		               conditionalPanel(condition="input.fitEnvCov == false",
		                                selectInput(inputId="MEModelEnvR","Choose ME Model ",c("Main Effect (G+E)","Homogeneous Variance (G+E+GxE)","Heterogeneous Variance (G+E+GxEi)"),selected = "Main Effect (G+E)",multiple = FALSE)
		               ),
		               conditionalPanel(condition="input.fitEnvCov == true",
		                                selectInput(inputId="MEModelEnvR","Choose ME Model ",c("Main Effect (G+E+W)","Homogeneous Variance (G+E+GxE+W)","Heterogeneous Variance (G+E+GxEi+W)"),selected = "Main Effect (G+E+W)",multiple = FALSE)
		               ),
		             tags$br(),
		             downloadButton("ExportOutME", "Export Predictions Table")
		           )
		    ),
		    column(8,
		           box(width=12,title="Multi-environment Genomic Prediction",solidHeader=T,status="primary",
		               tags$h4("Environments are defined as combinations of year and locations in which phenotypic data are collected.
                              Multi-environmental models are implemented using the EnvRType/BGGE pipeline as well as 'mmer' function in 'sommer' package.
                              Fit genomic prediction models taking into account only the main effects or the main effects + GxE effects.
                              The user needs to select one or many years and locations and the type of variance-covariance structure for fitting the model"
		               ),
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

		                fluidRow(
		                 box(width=12,title="",solidHeader = T, status="info",collapsible=TRUE,
		                  verbatimTextOutput("messageMEGPRun")
		                 )
		                )
		           ),
		               tags$br(),
		               tabBox(
		                 title = "Multi-environmental GP Table ",
		                 side = "left",
		                 width = 12,  # Adjust as needed
		                 tabPanel(
		                   "Table View",
		                   div(
		                     tableOutput("Ranked_Lines_for_SelectionSTME"),
		                     style = "overflow-x: auto; width: 100%;"  # Prevents horizontal scrolling issues
		                   )
		                 )
		               )
		        )
		     )
		  )
	 		   #FR            
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






