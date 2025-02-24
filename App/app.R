if(!require(shiny)){
  install.packages("shiny")
}

if(!require(shinyBS)){
  install.packages("shinyBS")
}


if(!require(shinyjs)){
  install.packages("shinyjs")
}

if(!require(shinyWidgets)){
  install.packages("shinyWidgets")
}


if(!require(reticulate)){
  install.packages("reticulate")
}


if(!require(renv)){
  install.packages("renv")
}


library(shiny) 
library(shinyjs)
library(shinyBS)
library(shinyWidgets)
library(reticulate) 
library(renv)

#renv::activate()

options(java.parameters = "-Xmx30g")
installed_packages <- renv::dependencies()$Package
lapply(installed_packages, library, character.only = TRUE)


# reticulate::use_virtualenv("./pyEnv", required = TRUE)
# reticulate::py_config() 
#setwd("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2App/App/")

### Change source file path to working directory

FN <- paste(getwd(),"/GS_Pipeline_Jan_2022_FnsApp.R",sep="")
source(FN)

PN <-  paste(getwd(),"/GSPipeline.png",sep="")

####


reticulate::use_virtualenv("./renv/python/virtualenvs/renv-python-3.12", required = TRUE)

## Set the file upload limit to X MB
options(shiny.maxRequestSize=1000*1024^2)

ui <- fluidPage( 
  #theme = bslib::bs_theme(bootswatch = "solar"),
  tags$style(HTML(".message-text { color: red; }")),
  useShinyjs(),
  
  fluidRow(
    column(4),
    column(6, div(
      id = "app-title",
      titlePanel(tags$strong("SOYGEN2 Genomic Selection App")),
      tags$p("An app to perform genomic predictions given genotypic and phenotypic data")
    ))),
  
  tabsetPanel(id="inData",
                tabPanel("Home",
                       fluidRow(
                         column(5),column(width=4,tags$h3(tags$strong("SOYGEN2 App")))
                         # column(6,offset=6),actionButton(inputId ="Home_Data", "next")
                       ),
                       fluidRow(
                         column(2),column(width=9,tags$h5("The SOYGEN2 App implements a genomic selection pipeline designed by public soybean breeders in the US MidWest. 
                        The GS pipeline involves steps starting from data management to making selection decisions based on genome estimated breeding values.The steps 
                        in the pipeline are depicted in the schematic below. The current version starts with the upload of filtered and imputed data exported from databases.") 
                         )),
                       tags$br(),
                       tags$br(),
                       fluidRow(
                         column(2),column(width=9,tags$a(tags$img(src="GSPipeline.png",title="Genomic Selection Pipeline",width="750",height="400")))), 
                       
                       fluidRow(
                         column(2),column(width=9,tags$h6(("Contributors: Vishnu Ramasubramanian, Cleiton Wartha, Paolo Vitale, Sushan Ru,and Aaron Lorenz")))),
                       tags$br(),        
                       fluidRow(
                         column(2),column(width=9,tags$h6("Contact: Vishnu Ramasubramanian - vramasub@umn.edu")))
               ),
          ## 
               tabPanel("Genotypic Data Processing",
                    tabsetPanel(id="GenDatProcess",
              
                    ## Tab for loading data
              
                    tabPanel("Load Genotypic Data",
                      sidebarLayout(
                        sidebarPanel(
                          fluidRow(
                            column(1),column(width=8,
                                             div(style="display: inline-block; vertical-align: top;",
                                                 checkboxInput("InGenoFormat1","",TRUE)
                                             ),
                                             div(style="display: inline-block; vertical-align: top; padding-left: 10px;",
                                                 tags$h6(tags$strong("1) Training & Target genotypic data in a combined VCF & target line IDs in a csv file"))
                                             )
                            )), 
                          
                          tags$br(),
                          
                          selectInput(inputId="TargetIDCol","Choose Target ID Col",choices=NULL,multiple=FALSE),
                          
                          # fluidRow(
                          #   column(1),column(width=8,
                          #                    div(style="display: inline-block; vertical-align: top;",
                          #                        checkboxInput("InGenoFormat2", "", FALSE)
                          #                    ),
                          #                    div(style="display: inline-block; vertical-align: top; padding-left: 10px;",              
                          #                        tags$h6(tags$strong("2) Training & Target genotypic data in separate VCF files")))
                          #   )),
                          # 
                          # fluidRow(
                          #   column(1),column(width=8,
                          #                    div(style="display: inline-block; vertical-align: top;",
                          #                        checkboxInput("InGenoFormat3","",FALSE),
                          #                    ),
                          #                    div(style="display: inline-block; vertical-align: top;padding-left: 10px;",         
                          #                        tags$h6(tags$strong("3) Training genotypic data only in VCF format for Crossvalidation"))), 
                          #   )), 
                          tags$br(),
                          tags$br()
                        ),
                        
                       mainPanel(
                       
                         fluidRow(
                           column(3),column(width=4,tags$h3(tags$strong("Load Genotypic Data"))),   
                         ), 
                         tags$br(),
                         fluidRow(
                           column(2),column(width=8,tags$h4(tags$strong("You can upload genotypic data in the format in the side panel. Select your input format and load data.")))
                         ),
                        tags$br(),
					              fluidRow(
                         column(1),column(width=3,tags$h4(tags$strong("Genotypic Data"))),
                         column(5,offset=2,tags$h4(tags$strong("Target Population IDs"))),
                        ),
                        
					                 
                        tags$br(),
                        conditionalPanel(condition="input.InGenoFormat1 == true",
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
                                        
                                       
                                        
                        ),
         
					              
# 					             conditionalPanel(condition="input.InGenoFormat2 == true",
#                                         fluidRow(
#                                           column(1),column(width=4,tags$h5(tags$strong("Upload Training Genotypic Data (VCF)"))),
#                                           column(6),column(width=4,tags$h5(tags$strong("Upload Target Genotypic Data (VCF)"))),
#                                         ),
#                                         fluidRow(
#                                           # column(1),column(width=5,fileInput("infileVCF", "Choose Training Genotype File (VCF)", accept = ".vcf")),
#                                            column(1),column(width=5,fileInput("infileVCF", "", accept = ".vcf")),
#                                            column(6),column(width=5,fileInput("infileTargetVCF","",accept = ".vcf")),
#                                         ),
#                                         tags$br(),
#                         ),
#                        
#                         conditionalPanel(condition="input.InGenoFormat3 == true",
#                                         fluidRow(
#                                           column(1),column(width=4,tags$h5(tags$strong("Upload Training Genotypic Data (VCF)"))),
#                                           #column(6),column(width=4,tags$h5(tags$strong("Upload Target Genotypic Data (.vcf)"))),
#                                         ),
#                                         fluidRow(
#                                           column(1),column(width=5,fileInput("infileVCF", "", accept = ".vcf")),
#                                           #column(6),column(width=5,fileInput("infileTargetVCF", "Choose Target Genotype  File (VCF)", accept = ".vcf")),
#                                         ),
#                                         
#                                         tags$br(),
#                        ),
                       
					            
                       
                       # fluidRow(
                       #   column(1),column(width=5,tags$h4(tags$strong(textOutput("GenoHeader")))),
                       #   column(6),column(width=5,tags$h4(tags$strong(textOutput("TargetHeader")))),
                       # ),
                       # fluidRow(
                       #   column(1),column(width=5,tags$h6(tableOutput("GenoTable"))),
                       #   column(6),column(width=5,tags$h6(tableOutput("TargetTable")))
                       # ),

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
                       # fluidRow(column(1),column(width=4,tags$h4(tags$strong(textOutput("GenoFiltHeader")))),
                       #          column(7),column(width=4,tags$h4(tags$strong(textOutput("GenoFiltHeader2"))))),      
                       # tags$br(),
                       # fluidRow(column(1),column(width=4,tags$h6(tableOutput("FilteredGenoTable"))),
                       #          column(7),column(width=4,tags$h6(tableOutput("FilteredGenoTable2")))),
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
                     
                     
                     
                        fluidRow(column(1),column(width=4,tags$h6(verbatimTextOutput("messageGenoFilt1"))),
                                 column(7),column(width=4,tags$h6(verbatimTextOutput("messageGenoFilt2")))),
                       
                     
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
                                            fileInput("inHapLib", "Load Haplotype Library (.phase)", accept = ".phase"),
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
                                            # 
                                            # fluidRow(column(2),column(width=8,tags$h4(tags$strong(textOutput("GenoImpHeader"))))),
                                            # tags$br(),
                                            # fluidRow(column(2),column(width=8,tags$h6(tags$strong(tableOutput("ImputedGenoTable"))))),
                                            # tags$br(),
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
                                  numericRangeInput(inputId="traitColsNum","Choose all the trait Cols",value=NULL,min=1,max=1),
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
                  
                 conditionalPanel(condition="input.chkPhenoME == false",
                                        
                 tags$br(),
                 fluidRow(
                           column(width=5,tags$h4(tags$strong("Phenotypic Data"))),
                           column(6),column(width=3,tags$h4(tags$strong("Select trait(s)"))),
                         ),
                         tags$br(),
                                
                         fluidRow(
                           column(width=3,fileInput("infileBLUEsSE", "Choose Phenotype File (CSV)", accept = ".csv")),
                           column(width=1,actionButton("iPh", "i")),column(7),
                           column(width=5,selectInput(inputId="trait","Choose One or More Traits",choices=NULL,multiple=TRUE)),
                         ),
                         fluidRow(
                           column(width=5,checkboxInput("header", "Header", TRUE)),
                         ),
                         tags$br(),
                         
                         # fluidRow(
                         #  column(width=6,tags$h4(tags$strong(textOutput("PhenoHeader")))),
                         #   column(7),column(width=5,div( 
                         #     tags$h4(tags$strong(textOutput("summaryHeader")))))
                         # ),
                         # fluidRow(
                         #   column(width=5,tags$h6(tableOutput("PhenoTable"))),
                         #   tags$h6(tableOutput("Summary"))
                         # ),
                 
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
                   
                   tags$br(),
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
                                        
                                        # fluidRow(
                                        #   column(width=5,selectInput(inputId="traitCols","Choose all the trait Cols",choices=NULL,multiple=TRUE)),
                                        # ),
                                        # 
                                        # fluidRow(
                                        #   column(width=5,selectInput(inputId="IDColME","Choose Uniq ID Col",choices=NULL,multiple=FALSE)),
                                        # ),
                                        # 
                                        # fluidRow(
                                        #   column(width=5,selectInput(inputId="strainME","Choose Strain Col",choices=NULL,multiple=FALSE)),
                                        # ),
                                        
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
                                        
                                      tags$br(),
                                      uiOutput("LocDistribution")
                                        
                     )
                 )
               )
            ),
            tabPanel("Enviromics Data",
                
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
                  fluidRow(
                    
                    column(width=10, tags$p("Extract weather data from NASA POWER database and estimate environmental relationship matrix using the envRType pipeline.
                                            The user needs to upload the co-ordinates file with the following columns 'Location','Country','Lat','Long','OtherLocName','LocName'. 
                                            The 'OtherLocName' and 'LocName' columns could contain abbreviations of location names.") 
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
                  
                    fluidRow(column(width = 5,verbatimTextOutput("messageEnvK"))),
                   
                   plotOutput("envPlot",height = "500px", width = "700px")
                  # uiOutput("envKPanel")
                  
                )
              )
             )
            )
          )
        ,
            
          
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
                           fluidRow(
                             column(2),column(width=8,tags$h5(tags$strong(textOutput("Mssg"))))),
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
                                              #MssgMT {
                                               color: red;
                                              }")
                                 )
                                ),
                                fluidRow(
                                column(2),column(width=8,tags$h5(tags$strong(textOutput("MssgMT"))))),
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
                                   # fluidRow(
                                   #   column(1),column(width=8,tags$h5(tags$strong(textOutput("cvrHeaderME"))))),
                                   # tags$br(),
                                   # fluidRow(
                                   #   column(2),column(width=8,tags$h5(tags$strong(textOutput("MssgME"))))),
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
                                                tags$h4(textOutput("RankedLinesHeader")),
                                                tableOutput("Ranked_Lines_for_SelectionST"),
                                                uiOutput("scatter")
                                              )
                                            )
                                   ),
                                   
                                   tabPanel("Multi-trait",
                                            
                                            sidebarLayout(
                                              sidebarPanel(
                                                selectInput(inputId="TrainSetMT","Choose Training Set",c("Complete Input Genotype Set","Optimal Train Set from Step 3","Random Train Set from Step 3")),
                                                tags$br(),
                                                tags$br(),
                                                selectInput(inputId="GPModelMT","Choose Prediction Model for Multiple Traits",c("BRR (BGLR)","RKHS (BGLR)","Spike-Slab(BGLR)","Mmer (Sommer)")),
                                                actionButton("RunPredictionsMT", "Predict Multiple Traits!"),
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
                                                  column(width=10, tags$p("Select the statistical method to train the genomic prediction model and predict the values of target lines.Multi-trait predictions are implemented using the BGLR and Sommer packages. The Multitrait function in BGLR implements Bayesian Ridge Regression, RKHS, and Spike-Slab methods (Perez-Rodriguez P and de Los Campos G, 2022).The multi-trait GBLUP is implemented using the 'mmer' function in 'sommer' package (Covarrubias-Pazaran G 2016)."))),
                                                tags$br(),
                                                tags$br(),
                                                #tags$h4(textOutput("RankedLinesHeader")),
                                                tableOutput("Ranked_Lines_for_SelectionMT"),
                                                uiOutput("scatterMT")
                                                
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
                                                

                                                selectInput(inputId="Package","Choose Package for MultiEnvironment GP Modeling ",c("BGGE-EnvRType","SOMMER")),

                                                conditionalPanel(condition="input.Package == 'BGGE-EnvRType'",

                                                     selectInput(inputId="GKernelMet","Choose Genotype Kernel Method",c("Linear","Gaussian"),selected="Linear",multiple=FALSE),
                                                     tags$br(),
                                                     checkboxInput("fitEnvCov", "Include Enviromics Kernel from Step 2", FALSE),

                                                     tags$br(),
                                                     tags$h5(tags$strong("View Output Table from Model")),
                                                     selectInput(inputId="MEModelEnvR","Choose ME Model ",c("Main Effect (G+E)","Homogeneous Variance (G+E+GxE)","Heterogeneous Variance (G+E+GxEi)"),selected = "Main Effect (G+E)",multiple = FALSE),
                                                     tags$br(),
                                                     tags$br(),
                                                ),
                                                conditionalPanel(condition="input.Package == 'SOMMER'",
                                                    selectInput(inputId="MEModelSom","Choose Model ",c("Main Effect","Homogeneous Variance (CS)","Heterogeneous Variance (CS+DG)","Unstructured (US)")),
                                                    tags$br(),
                                                ),
                                                tags$br(),
                                                tags$br(),
                                                actionButton("RunPredictionsME", "Fit Multi-environmental Model!"),
                                                tags$br(),
                                                tags$br(),
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
                                                #tags$h4(textOutput("RankedLinesHeader")),
                                                tableOutput("Ranked_Lines_for_SelectionSTME"),
                                                #tableOutput("Ranked_Lines_for_SelectionMTME")

                                               )
                                        )
                                   )

                       )
              )
  )
  
)           

###

server <- function(input,output,session){
  
  # Geno
  
  # Geno <- reactive({
  #   
  #   genoFile <- input$infileVCF
  #   ext <- tools::file_ext(genoFile$datapath)
  #   req(genoFile)
  #   validate(need(ext == "vcf", "Please upload a vcf file"))
  #   
  #  # withProgress(message = 'Reading Data', value = 0, {
  #  #   NUST_Genotypes_VCF <- read.table(genoFile$datapath)})
  #   
  #   # withProgress(message = 'Converting VCF to Dataframe', value = 0, {
  #   #   gt2d <- VCFtoDF(genoFile$datapath) }) 
  #   
  #   
  #    gt2d
  # })
  
  # InGenoFormat1 <- reactive(input$InGenoFormat1) 
  # InGenoFormat2 <- reactive(input$InGenoFormat2) 
  # InGenoFormat3 <- reactive(input$InGenoFormat3) 
  
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
  # genoHead <- eventReactive(input$infileVCF,{ paste("Genotype Table with ",ncol(Geno())-5," lines and ",nrow((Geno()))," markers",sep="")})
  # output$GenoHeader <- renderText({genoHead()})
  # output$GenoTable <- renderTable({as.data.frame((Geno())[1:5,1:5])})
  # 
  
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
  
  output$messageGenoStats <- renderText({
    
    invalidateLater(1000, session) # Update every second
    if(file.exists(temp_file0a())){
      lines <- readLines(temp_file0a(), warn = FALSE)
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
  
  
  TargetIDs <- reactive(NULL)
  
  
  observeEvent(input$infileTargetTable,{updateSelectInput(inputId="TargetIDCol",choices=colnames(TargetTab()))})
  
  TargtIDCol <- reactive(input$TargetIDCol)
  TargetIDs <- eventReactive(input$infileTargetTable,{
                  as.character(TargetTab()[,TargtIDCol()])
               })
  
  buildLib <- reactive({
    libFile <- input$inHapLib
    ext <- tools::file_ext(libFile$datapath)
    
    req(libFile)
    validate(need(ext == "phase", "Please upload a .phase file"))
    
    read.table(libFile$datapath, header = FALSE)
  })
  
  observeEvent(input$inHapLib,{ 
     shinyjs::enable("impute_APIdata")
   
     write.table(buildLib(),"lib.phase",quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)
     
     
     # createAlert(session, "alert", "alert1", title = "Process Completed",
     #            content = "The build file has loaded. You can now impute", style = "success", append = FALSE)
     # 
    
  })
  
  
  
#### Target Table 

  
  # TargetHead <- eventReactive(input$infileTargetTable,{ paste("Table with information on ",nrow(TargetTab())," Target lines",sep="")})
  # output$TargetHeader <- renderText({TargetHead()})
  # output$TargetTable <-  renderTable({
  #      TargetTableOut <- as.data.frame(TargetTab()[1:5,]) 
  #      #colnames(TargetTableOut) <- ""
  #      TargetTableOut
  # })
  
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
  
  output$messageTargetStats <- renderText({
    
    invalidateLater(1000, session) # Update every second
    if(file.exists(temp_file0d())){
      lines <- readLines(temp_file0d(), warn = FALSE)
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

output$messageGenoFilt1 <- renderText({
  
  invalidateLater(1000, session) # Update every second
  if(file.exists(temp_file0b())){
    lines <- readLines(temp_file0b(), warn = FALSE)
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
  
  output$messageGenoFilt2 <- renderText({
    
    invalidateLater(1000, session) # Update every second
    if(file.exists(temp_file0c())){
      lines <- readLines(temp_file0c(), warn = FALSE)
      return(paste(lines, collapse = "\n"))
    }else {
      return("Waiting for output...")
    }
  })
  
  
  
  ## Tas to DF
  
  Geno_DF <- reactive({
   # browser()
    if(!is.null(GenoTas())){getGenoTas_to_DF(GenoTas())}
  })
  
  GenoFilt1_DF <- reactive({
    #browser()
    if(!is.null(GenoFilt1())){ getGenoTas_to_DF(GenoFilt1())}
  })
  
  GenoFilt2_DF <- reactive({
    if(!is.null(GenoFilt2())){ getGenoTas_to_DF(GenoFilt2())}
  })
  
  # Geno_DF <- reactiveVal(NULL)
  # observeEvent(input$infileVCF, {Geno_DF(getGenoTas_to_DF(GenoTas()))})
  # # 
  # GenoFilt1_DF <- eventReactive(input$FilterSites,{getGenoTas_to_DF(GenoFilt1())})
  # GenoFilt2_DF <- eventReactive(input$FilterTaxa,{getGenoTas_to_DF(GenoFilt2())})
  # 
  
  # GenoFilt1_DF <- reactive(getGenoTas_to_DF(GenoFilt1(),Geno()))
  # GenoFilt2_DF <- reactive(getGenoTas_to_DF(GenoFilt2(),Geno()))
  
  # Geno_DF <- reactive(getGenoTas_to_DF(GenoTas()))
  # GenoFilt1_DF <- reactive(getGenoTas_to_DF(GenoFilt1()))
  # GenoFilt2_DF <- reactive(getGenoTas_to_DF(GenoFilt2()))
  
 
  
  # Geno_DF <- reactive(as.data.frame(as.matrix(GenoTas())))
  # GenoFilt1_DF <- reactive(as.data.frame(as.matrix(GenoFilt1())))
  # GenoFilt2_DF <- reactive(as.data.frame(as.matrix(GenoFilt2())))
  # 
  
  # ####  
  # genoFiltHead <- eventReactive(input$FilterSites,{ paste("Filter Sites: Genotype Table with ",ncol(GenoFilt1_DF())-5," lines and ",nrow((GenoFilt1_DF()))," markers",sep="")})
  # output$GenoFiltHeader <- renderText({genoFiltHead()})
  # output$FilteredGenoTable <- renderTable({as.data.frame((GenoFilt1_DF())[1:5,1:6])})
  # 
  # ###
  # genoFiltHead2 <- eventReactive(input$FilterTaxa,{ paste("Filter Taxa: Genotype Table with ",ncol(GenoFilt2_DF())-5," lines and ",nrow((GenoFilt2_DF()))," markers",sep="")})
  # output$GenoFiltHeader2 <- renderText({genoFiltHead2()})
  # output$FilteredGenoTable2 <- renderTable({as.data.frame((GenoFilt2_DF())[1:5,1:6])})
  # 
  #####
  
  # Debug output to show the value of input.InGenoFormat3
  # output$debugOutput <- renderText({
  #   paste("input.InGenoFormat3:", input$InGenoFormat3)
  # })
  # 
  
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
  
  
  GenoImp_DF1 <- eventReactive(input$Impute,{
    getGenoTas_to_DF(GenoImp())
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
  
  output$messageImpGeno1 <- renderText({
    
    invalidateLater(1000, session) # Update every second
    if(file.exists(temp_file0e())){
      lines <- readLines(temp_file0e(), warn = FALSE)
      return(paste(lines, collapse = "\n"))
    }else {
      return("Waiting for output...")
    }
  })
  
  
  
  ### if you use only reactive, this will throw an error ncol(GenoImp_DF())-5
  
  #genoImpHead <- eventReactive(input$Impute,{paste("Genotype Table with ",ncol(GenoImp_DF())-5," lines and ",nrow((GenoImp_DF()))," markers",sep="")})
  
  # genoImpHead <- reactive({paste("Genotype Table with ",ncol(GenoImp_DF())-5," lines and ",nrow((GenoImp_DF()))," markers",sep="")})
  
  # 
  # genoImpHead <- reactive({
  #   # Ensure GenoImp_DF is only called when it has updated data
  #   df <- as.data.frame(GenoImp_DF1()) # This call ensures dependency
  #   
  #   if(is.data.frame(df) && ncol(df) > 5) {
  #     #     # Successful case
  #     paste("Genotype Table with ", ncol(df) - 5, " lines and ", nrow(df), " markers", sep = "")
  #   } else {
  #     #     # Fallback or error message
  #     "Waiting for data"
  #   }
  # })
  # # 
  # # UI and server logic to render the header text and imputed genotype table
  # output$GenoImpHeader <- renderText({
  #   genoImpHead()  # This will now wait for GenoImp_DF to complete
  # })
  # 
  # 
  # # output$GenoImpHeader <- renderText({genoImpHead()})
  # output$ImputedGenoTable <- renderTable({as.data.frame((GenoImp_DF1())[1:5,1:10])})
  # # 
  ### 
  
  observeEvent(input$imputeMet,{ 
    
    if(impMethod()=="AlphaPlantImpute"){ 
      # 
      
      reticulate::py_run_string("import sys")
      
      API_exists <- py_module_available("alphaplantimpute2")
      print(paste("API Available",API_exists,sep=" "))
      
      output$APIUI <- renderUI({
       
        # Check if AlphaPlantImpute is selected
  
        fluidPage(
          
          fluidRow(
            #column(1),
            column(width=4,tags$strong(tags$h3("Build Haplotype Library "))),
            column(7),column(width=4,tags$strong(tags$h3("Impute Genotype Table ")))
          ),
          
          tags$br(),
          
          fluidRow(
            #column(1),
            column(width=4,
                   numericInput(inputId="nHap",label = "Enter number of haplotypes",value=20,min=2,max=100),
            ),
            column(7),column(width=5,
                             tags$strong(("Input founders file for pedigree based imputation. No need to load any file for population based imputation.")))
          ),
          tags$br(),
          
          fluidRow(
            #column(1)
            column(width=4,
                   numericInput(inputId="nSampRnds",label = "Enter number of sample rounds",value=5,min=2,max=100),
            ),
            column(7),column(width=4, 
                             fileInput("founder_file", "Select Founder Files", multiple = FALSE,accept = ".txt")),
          ),
          tags$br(),
          
          fluidRow( 
            #column(1),
            column(width=3,
                   numericInput(inputId="HDthresh",label = "Enter non-missing threshold for high density genotypes",value=0.9,min=0.5,max=1),
            )),
          tags$br(),
          
          # Build library button
          fluidRow( 
            column(2),column(width=3,
                             actionButton("build_library", "Build Library"),
            ),
            column(7),column(width=3,
                             bsButton("impute_APIdata",label="Impute Data",style="primary",disabled=TRUE),
            ),
          ),   
          tags$br(),
          tags$br(),
          
          tags$head(
            tags$style(HTML("
             #message {
        
              /*max-height: 1200px; Set maximum height */
              overflow-y: scroll; /* Enable vertical scrolling */
              overflow-x: scroll;  /*Hide horizontal scrolling */
              overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
              width: 400px; 
              /*max-width: 100%; */
              padding: 6px 12px;
              height: 300px;
            }
          "))
         ),
          
          tags$head(
            tags$style(HTML("
        #message2 {
        
        overflow-y: scroll; /* Enable vertical scrolling */
        overflow-x: scroll;  /*Hide horizontal scrolling */
        overflow-wrap: anywhere; /* Ensure long words do not cause horizontal scrolling */
        width: 400px; 
        /*max-width: 100%; */
        padding: 6px 12px;
        height: 300px;
      }
    "))
          ),
          fluidRow(
            
            column(width = 5,
                   verbatimTextOutput("message")
            ),
            column(width = 1),  # Adjust the width or remove if not needed
            column(width = 5,
                   verbatimTextOutput("message2")
            )
          ),
          
          
          
        )
      })
      
    } 
    
  })
  
  
  
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
  
  observeEvent(input$build_library,{
    #browser()
    # Extract necessary input values
    # Path for the temporary file
    temp_file(tempfile())
    
    # Start the process in a separate R process
    rProcess <- callr::r_bg(function(nHap, nSampRnds, HDthresh, temp_file) {
      library(reticulate)
      #reticulate::use_virtualenv("./pyEnv", required = TRUE)
      reticulate::use_virtualenv("./renv/python/virtualenvs/renv-python-3.12", required = TRUE)
      
      sys <- reticulate::import("sys")
      api2 <- reticulate::import("alphaplantimpute2.alphaplantimpute2")
      
      sys$argv <- c('alphaplantimpute2', '-createlib', '-out', 'lib',
                    '-genotypes', 'current_GenoTable.genotypes',
                    '-n_haplotypes', as.character(nHap),
                    '-n_sample_rounds', as.character(nSampRnds),
                    '-hd_threshold', as.character(HDthresh),
                    '-seed', '42')
      
      sink(temp_file)
      api2$main()
      sink()
    }, args = list(nHap(), nSampRnds(), HDthresh(), temp_file()), stdout = temp_file(), stderr = temp_file())
    
### Test output 
    
    on.exit({
      if (exists("rProcess") && rProcess$is_alive()) {
        rProcess$kill()
      }
    }, add = TRUE)
    
    
    ###
    
    processAlive <- reactive(rProcess$is_alive())
    # Periodically read the file and update the UI
    output$message <- renderText({
      
      invalidateLater(1000, session) # Update every second
      
      if(file.exists(temp_file())){
        lines <- readLines(temp_file(), warn = FALSE)
        return(paste(lines, collapse = "\n"))
      }else {
        return("Waiting for output...")
      }
      
      if(!processAlive()){
        
        if (file.exists(temp_file())) {
          lines <- readLines(temp_file(), warn = FALSE)
          txt <- paste(lines,collapse="\n")
          shinyjs::enable("impute_APIdata")
          
          # createAlert(session, "alert", "alert1", title = "Process Completed",
          #             content = "The build process is complete. You can now impute", style = "success", append = FALSE)
          # 
          return(paste0(txt,"\n","Build Completed",collapse=" "))
          
          
        }else {
          return("Waiting for output...")
        }
      }
    })
    
  })
  
  ##
  
  temp_file2 <- reactiveVal()
  processComplete2 <- reactiveVal(FALSE)
  ouputRead <- reactiveVal(FALSE)
  
  
  # Start the background process and return rProcess
  rProcess2 <- eventReactive(input$impute_APIdata,{
    
    library(reticulate)
    #reticulate::use_virtualenv("./pyEnv", required = TRUE)
    reticulate::use_virtualenv("./renv/python/virtualenvs/renv-python-3.12", required = TRUE)
    sys <- import("sys")
    api2 <- import("alphaplantimpute2.alphaplantimpute2")
    
    ## Extract necessary input values
    
    founder_file <- founder()
    
    # Path for the temporary file
    
    temp_file2(tempfile())
    
    # Start the process in a separate R process
    callr::r_bg(function(founder_file,temp_file){
      library(reticulate)
      reticulate::use_virtualenv("./pyEnv", required = TRUE)
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
      api2$main()
      # print("Imputing")
      # Sys.sleep(2)
      sink()
    }, args = list(founder_file,temp_file2()), stdout = temp_file2(), stderr = temp_file2())
    
  })
  
  
  # Reactive expression to check process status and read temp file
  processStatus2 <- reactive({
    invalidateLater(1000,session) # Update every second
    
    if (!is.null(rProcess2()) && rProcess2()$is_alive() && !ouputRead()){
      if (file.exists(temp_file2())) {
        lines <- readLines(temp_file2(), warn = FALSE)
        return(paste(lines,collapse = "\n"))
      }else {
        return("Waiting for output...")
      }
    }else if(!is.null(rProcess2()) && !rProcess2()$is_alive() && !ouputRead()){
      
      if (file.exists(temp_file2())) {
        lines <- readLines(temp_file2(), warn = FALSE)
        txt <- paste(lines, collapse = "\n")
        processComplete2(TRUE)
        return(paste(txt,"Imputation Completed",collapse="\n"))
        
      }else {
        return("Waiting for output...")
      }
      
    }else if(!is.null(rProcess2()) && !rProcess2()$is_alive() && ouputRead()){
      
      if (file.exists(temp_file2())) {
        lines <- readLines(temp_file2(), warn = FALSE)
        txt <- paste(lines, collapse = "\n") 
        
        return(paste0(txt,"\n","Imputation Completed","\n", "Imputed Geno Table has Dim:",paste(dim(GenoImp_DF()),collapse=" ")))
      }else{
        return("Output not read")
      }
    }else{return("Var Error")}
  })
  
  # Update the UI with process status
  output$message2 <- renderText({
    processStatus2()
  })
  
  
  # Process and return the output once the background process is complete
  observe({
    # Ensure the background process is completed
    
    req(processComplete2())
    ouputRead(TRUE)
  })
  
  ### Read Imputed genotypic data and prepare output 
  
  # #impGeno <- reactive({
  GenoImp_DF2 <- reactive({
    # Ensure outputRead is set to TRUE
    req(ouputRead())
    # ... Your existing code to process and return the output ...
    
    ImpGeno <- getImpGenoData_API()
    ImpGeno
  })
  
  GenoImp_DF <- reactive({
    if (!is.null(GenoImp_DF1()) && nrow(GenoImp_DF1()) > 0) {
      # GenoImp_DF1 is not NULL and has data
      return(GenoImp_DF1())
    } else if (!is.null(GenoImp_DF2()) && nrow(GenoImp_DF2()) > 0) {
      # GenoImp_DF2 is not NULL and has data, and GenoImp_DF1 was NULL or had no data
      return(GenoImp_DF2())
    } else {
      # Both are NULL or have no data, return NULL or a default value/data frame
      return(NULL)  # Or any default value you deem appropriate
    }
  })
  
  
  # genoImpHead <- reactive(paste("Genotype Table with ",ncol(GenoImp_DF())-5," lines and ",nrow((GenoImp_DF()))," markers",sep=""))
  # GenoImp_DF <- eventReactive(input$Impute,{getGenoTas_to_DF(GenoImp(),Geno())})
  #
  #GenoImp_DF <- reactive(getGenoTas_to_DF(GenoImp(),Geno()))
  #GenoImp_DF <- reactive(getGenoTas_to_DF(GenoImp()))
  #GenoImp_DF <- reactive(as.data.frame(as.matrix(GenoImp())))
  
  
  
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
    updateNumericRangeInput(inputId ="traitColsNum",value= c(2,ncol(Pheno())))})
  
  observeEvent(input$infileBLUEsSE, {
   updateSelectInput(inputId="IDColSE",choices=colnames(Pheno()))})
  
  observeEvent(input$infileBLUEsSE, {
   updateSelectInput(inputId="strainSE",choices=colnames(Pheno()))})
  
 
###
  
 
  IDColSE <- reactive(input$IDColSE)
  StrainSE <- reactive(input$strainSE)
  
  TraitColsSE <- reactive(c(input$traitColsNum[1]:input$traitColsNum[2]))
  TraitsPhSE <- reactive(colnames(Pheno())[TraitColsSE()])
  TraitsPhSEVec <- reactive(paste(TraitsPhSE()," ",sep="",collapse=","))
  
####  
  
  temp_file3d <- reactiveVal('none')
  
  observeEvent(input$traitColsNum,{
   
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

  output$messagePhSE <- renderText({
    
    invalidateLater(1000, session) # Update every second
    if(file.exists(temp_file3d())){
      lines <- readLines(temp_file3d(), warn = FALSE)
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
  
    # observeEvent(input$infileBLUEsSE,{
    #   if(nSelTraits()>1){
    #     lapply(1:nSelTraits(),function(i){
    #       loc_i <- reactive(i)
    #       distOutput[[noquote(paste("'","plot",loc_i(),"'",sep=""))]] <- renderPlot({
    #         getDistributionPlot() 
    #        })
    #       
    #     })
    #   }
    # })  
    # 
  
 
  
  
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
         return(paste(lines, collapse = "\n"))
       }else {
         return("Waiting for output...")
       }
     })
     
   
   
  
  ## Trait 
  
  nTraits <- eventReactive(input$infileBLUEsSE,{(ncol(Pheno())-1)})
  Trait <- reactive(input$trait)
  
  ### 
  
  MssgTargetSet <- eventReactive(input$trait,{
    if(nrow(processedData()[[2]])==0 | is.null(processedData()[[2]])){
      response <- "The current genotypic file doesn't have target line genotypes.
        Load genotypic data file with Target IDs"
      return(response)
    }
  })

  # output$MssgTarget <- renderUI({
  #   message <- MssgTargetSet()
  #   HTML(paste0("<div class='message-text'>", message, "</div>"))
  # })
  
  #summaryHead <- eventReactive(input$trait,{paste("Summary of trait values")})
 
 
  temp_file3a <- reactiveVal('none')
  observeEvent(input$trait,{
    # Path for the temporary file
    temp_file3a(tempfile())
    
    if(length(Trait())>0 & !is.null(Pheno())){
      sink(temp_file3a())
      cat("Summary of selected trait values : \n")
      cat(paste(names(round(summary(Pheno()[,Trait()[1]]),digits=2)),collapse="\t"))
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
  
  output$messageTrtSE <- renderText({
    
    invalidateLater(1000, session) # Update every second
    if(file.exists(temp_file3a())){
      lines <- readLines(temp_file3a(), warn = FALSE)
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
  
  TraitME <- reactive({input$traitME})
  nSelTraitsME <- reactiveVal(NULL)
  # 
  observeEvent(input$traitME,{
    if(!is.null(TraitME())){
      nSelTraitsME(length(TraitME()))
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
    getPhenoMEData(PhenoME(),TraitME(),nSelTraitsME(),IDColsME(),StrainME())
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
      return(paste(lines, collapse = "\n"))
    }else {
      return("Waiting for output...")
    }
  })
  
  
### Merge geno and pheno data    
 
  isPhenoME <- reactive(input$chkPhenoME)
  
  mergedData <-  reactive({
    print("mergeIn")
    #browser()
    t_Geno <- reactive(GenoPre())
    
    withProgress(message = 'Merging Pheno and Geno Data', value = 0, {
    
      if(isPhenoME()==FALSE){
        
        getMergedData(t_Geno(),Pheno(),TargetIDs())
      }else if(isPhenoME() == TRUE){ 
        getMergedDataME(phenoMEData(),t_Geno(),TargetIDs())
      }
    })
  
  })
  
###   
  

  processedData <- reactive({ 
    print("prIn")
    if(isPhenoME() == FALSE){
      getProcessedData(mergedData(),Trait())
    }else if(isPhenoME() == TRUE){ 
       mergedData()
    }
  })
  
  
  DT_Filt_List <- reactive({
    #browser()
    print("prIn2")
    if(isPhenoME()){
      processedData()[[1]]
    }
  })
  
  genoDat_List <- reactive({ 
    if(isPhenoME()){
      processedData()[[2]]
    }
  })
  
  
  plotData <- reactiveVal(NULL)
  
  observeEvent(input$traitME,{
    print("in")
    if(nSelTraitsME()==1){
      print("in1")
      plotData(getPhenoDistbnPlots(DT_Filt_List(),TraitME(),1))
      
    }else if(nSelTraitsME()>1 & !is.null(DT_Filt_List())){
       
         plotList <- lapply(1:nSelTraitsME(),function(loc_i){
             getPhenoDistbnPlots(DT_Filt_List(),TraitME()[loc_i],loc_i)
         })
         plotData(plotList)
     }
  })
    
    # else if(nSelTraitsME()>1){
    #     lapply(1:nSelTraitsME(),function(i){
     #       loc_i <- reactive(i)
    #       output[[noquote(paste("'","distplot",loc_i(),"'",sep=""))]] <- renderPlot({
    #         getPhenoDistbnPlots(DT_Filt_List(),TraitME()[loc_i],loc_i)
    #       })
    #     })
    #   }
    
  #})

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
 
 
 
 # 
 # outputListPhME <- reactive({
 #   if(nSelTraitsME()>1){
 #     lapply(1:nSelTraitsME(), function(i){
 #       plotOutput(noquote(paste("'","distPlot",i,"'",sep="")))
 #     })
 #   }
 # })
 
 # Link plotRenderList and outputList
 
 observeEvent(input$traitME,{
   #browser()
   dummyPlt()
   if(!is.null(TraitME())){
     if(nSelTraitsME()>1 & length(plotRenderList()) >1 & length(outputListPhME()) >1){
           tempList <-  plotRenderList()
           outputListPhME(tempList) 
     }else{print(paste("The Render plot list with length", length(plotRenderList())," or outputPlot list with length",length(outputListPhME()),
          "do not meet the requirements",sep=""))}
     
   }
 })
     
       # tempList2 <-  outputListPhME()
       # # 
       # # for(i in 1:nSelTraitsME()){
       # #   tempList2[[noquote(paste("'","distPlot",i,"'",sep=""))]] <-  tempList[[i]]
       # # }
 


#### 

output$LocDistribution <- renderUI({
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
#  
#  output$LocDistribution <- renderUI({
#   print("outLD")
#    if(nSelTraitsME()==1){
#      plotOutput("distPlots")
#    }
#    
#  })
 # 
 # output$LocDistribution <- renderUI({
 #   print("outLD")
 #   if(nSelTraitsME()==1){
 #     plotOutput("distPlots")
 #   }else if(nSelTraitsME()>1){ 
 #     distPlotOut_List <- lapply(1:nSelTraitsME(),function(x) {
 #     i <- reactive(x)
 #     plotOutput(noquote(paste("'","distPlot",i(),"'",sep="")))})
 #     do.call(tagList,distPlotOut_List)
 #   }
 # })
 # 
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


observeEvent(input$getEnvK,{
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

output$messageEnvK <- renderText({
 envDatStatus()  

})


processWthDat <- reactive(input$processWth)
gaussVar <- reactive(input$GaussKE)

# 
# EnvK <- eventReactive(input$getEnvK,{
#   print(dim(envData()))
#  getEnvKernel(envData(),processWthDat,gaussVar)
# })

EnvK <- reactive({
  #browser()
  withProgress(message='Estimating Environmental Relationships',value=0,{
     getEnvKernel(envData(),processWthDat(),gaussVar())})
 })




OtherLoc <- reactive(input$OtherLoc)
EnvK_Mod <- reactive({syncEnvPhenoDat(EnvK(),LocCoords(),OtherLoc())})



# Render the main plot
output$envPlot <- renderPlot({ 
  
  print("envPl")
  if(!is.null(envData())){ 
   if(OtherLoc()==FALSE){
      plotEnvRel(EnvK())
   }else if(OtherLoc()==TRUE){
     plotEnvRel(EnvK_Mod())
   }
  }
})
# 
#  
# output$envKPanel <- renderUI({
#     plotOutput("envPlot",height = "500px", width = "700px")
# }) 

##### 
  

  
  observeEvent(input$trait, {
    updateNumericInput(inputId ="noCandidates",value= nrow(processedData()[[1]]),min =2, max=nrow(processedData()[[1]]))}) 
  
  observeEvent(input$trait,{
    updateNumericInput(inputId ="noToSelect",value= 100,min =2, max=nrow(processedData()[[1]]))})
  
  
  
  # TS Optimization
  
  noCandidates <- reactive(input$noCandidates)
  nTrainToSelect <- reactive(input$noToSelect)
  optimCriteria <- reactive(input$optCriteria)
  
  predictionData <- reactive({getPredictionData(processedData(),noCandidates())}) 
  
  GAParameters <- reactiveValues(npop=100,nelite=10,mutprob=0.5,mutintensity=1,niterations=100,minitbefstop=50,tabu="TRUE",tabumemsize=1,plotiters="FALSE",errorstat="PEVMEAN2",lambda=1e-6,mc.cores=10)
  
  
  # 
  Train_STPGA <-   eventReactive(input$Optimize,{
    withProgress(message = 'Running Optimizations', value = 0, {
      getOptimalTS(predictionData(),unlist(Trait()),nTraits(),noCandidates(),nTrainToSelect(),isolate(reactiveValuesToList(GAParameters)))})
  })
  

  Train_Random <-   eventReactive(input$Optimize,{
    withProgress(message = 'Running Random Set', value = 0, { 
      getRandomTS(predictionData(),unlist(Trait()),nTraits(),noCandidates(),nTrainToSelect())})
  })
  
####
  
  TSOptOutputList <- eventReactive(input$Optimize,{
    if(nSelTraits()==1){
      withProgress(message = 'Performing CrossValidations', value = 0, {
        getTSComparisons(predictionData(),Train_STPGA(),Train_Random(),unlist(Trait()),nTraits(),TargetIDs())
      })
    }else if(nSelTraits()>1){
      withProgress(message = 'Performing CrossValidations', value = 0, {
        getTSComparisonsMT(predictionData(),Train_STPGA(),Train_Random(),unlist(Trait()),nTraits(),TargetIDs())
      })
    }
    
  })
  
  tsOptHead <- eventReactive(input$Optimize,{paste("Correlation between observed and predicted values")})
  output$tsOptHeader <- renderText({tsOptHead()}) 
  output$PredAccuracyTable <- renderTable({
    TSTable <- TSOptOutputList()
    
    colnames(TSTable) <- c("Ridge Regression","BayesB","Bayes LASSO")
    rownames(TSTable) <- c("Optimal Set","Random Set")
    TSTable},
    colnames=TRUE,rownames=TRUE,digits=2)
  
  
  
  
  ## CVR
  
  k <- reactive(input$k)
  nIter <- reactive(input$nIter)
  
  ### Conflicting statement messages
  
  
  Mssg2CVR <- eventReactive(input$CrossValidationMT,{
    
    mTraitMssg <- "Select more than one trait for cross validation of multi-trait models."
    meTraitMssg <- "You've loaded pheno data from multiple environments. Try CV for multi-environmental models. "
    if(nSelTraits()==1 && isPhenoME()){
      paste(mTraitMssg,meTraitMssg,sep=" ")
    }else if(nSelTraits()==1 && !isPhenoME()){
      mTraitMssg
    }else if(nSelTraits()>1 && isPhenoME()){
      meTraitMssg
    }else{""}  
  })
  
  output$MssgMT <- renderText({
    Mssg2CVR()
    
    #HTML(paste0("<div class='message-text'>", message, "</div>"))
  })
  
  
  
  # cvrOutputListMT <- eventReactive(input$CrossValidationMT,{ withProgress(message = 'Running CrossValidations', value = 0, {
  #         getMTCVR(predictionData(),unlist(Trait()),nTraits(),k(),nIter()) })
  # })
  # 
  # temp_file5 <- reactiveVal('none')
  # processComplete5 <- reactiveVal(FALSE)
  # cvrOutputListMT <- reactiveVal(NULL)
  # 
  # 
  # # Start the background process and return rProcess
  # rProcess5 <- eventReactive(input$CrossValidationMT,{
  #    
  #   # Path for the temporary file
  #   
  #   temp_file5(tempfile())
  #   
  #   # Start the process in a separate R process
  #   callr::r_bg(function(getMTCVR,predictionData,Trait,nTraits,k,nIter,temp_file5,cvrOutputListMT){
  # 
  #     library(BGLR)
  #     #library(sommer)
  #     sink(temp_file5)
  #     cvrOutputListMT(getMTCVR(predictionData,Trait,nTraits,k,nIter))
  #     sink()
  #   }, args = list(getMTCVR,predictionData(),unlist(Trait()),nTraits(),k(),nIter(),temp_file5(),cvrOutputListMT()), stdout = temp_file5(), stderr = temp_file5())
  #   
  # })
  # 
  # 
  # # Reactive expression to check process status and read temp file
  # processStatus5 <- reactive({
  #   invalidateLater(1000,session) # Update every second
  #   
  #   if (!is.null(rProcess5()) && rProcess5()$is_alive()){
  #     if (file.exists(temp_file5())) {
  #       lines <- readLines(temp_file5(), warn = FALSE)
  #       return(paste(lines,collapse = "\n"))
  #     }else {
  #       return("Waiting for output...")
  #     }
  #   }else if(!is.null(rProcess5()) && !rProcess5()$is_alive()){
  #     
  #     if (file.exists(temp_file5())) {
  #       lines <- readLines(temp_file5(), warn = FALSE)
  #       txt <- paste(lines, collapse = "\n")
  #       processComplete5(TRUE)
  #       return(paste(txt,"Crossvalidation Completed",collapse="\n"))
  #       
  #     }else {
  #       return("Waiting for output...")
  #     }
  #   }else{return("Error")}
  # })
  # 
  # # Update the UI with process status
  # output$messageMT5 <- renderText({
  #   processStatus5()
  # })
  # 
###
  
  temp_file5 <- reactiveVal('none')
  processComplete5 <- reactiveVal(FALSE)
  cvrOutputListMT <- reactiveVal(NULL)
  
  # Start the background process and return rProcess
  rProcess5 <- eventReactive(input$CrossValidationMT, {
    # Path for the temporary file
    temp_file5(tempfile())
    
    # Start the process in a separate R process
    callr::r_bg(function(getMTCVR, predictionData, Trait, nTraits, k, nIter, temp_file5_path) {
      library(BGLR)
      library(rrBLUP)
      # library(sommer)
      sink(temp_file5_path)
      # Simulate a long-running process
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
      temp_file5_path = temp_file5()
    ), stdout = temp_file5(), stderr = temp_file5())
  })
  
  # Reactive expression to check process status and read temp file
  processStatus5 <- reactive({
    invalidateLater(1000, session) # Update every second
    
    if (!is.null(rProcess5()) && rProcess5()$is_alive()) {
      if (file.exists(temp_file5())) {
        lines <- readLines(temp_file5(), warn = FALSE)
        return(paste(lines, collapse = "\n"))
      } else {
        return("Waiting for output...")
      }
    } else if (!is.null(rProcess5()) && !rProcess5()$is_alive()) {
      if (file.exists(temp_file5())) {
        lines <- readLines(temp_file5(), warn = FALSE)
        txt <- paste(lines, collapse = "\n")
        processComplete5(TRUE)
        cvrOutputListMT(rProcess5()$get_result())
        return(paste(txt, "Cross-validation Completed", collapse = "\n"))
      } else {
        return("Waiting for output...")
      }
    } else {
      return("Error")
    }
  })
  
  # Update the UI with process status
  output$messageMT5 <- renderText({
    processStatus5()
  })

  
  
  
###
  
  cvrHead <- eventReactive(input$CVEvent,{
    paste("Correlation between observed and predicted values for ",paste(Trait(),collapse=" and "),sep="")})
  
  output$cvrHeader <- renderText({cvrHead()}) 
  
  cvrOutputListST <- eventReactive(input$CrossValidationST,{
    
    if(nSelTraits()==1){
      PATab <- withProgress(message = 'Running CrossValidations', value = 0, {
        getemCVR(predictionData(),unlist(Trait()),nTraits(),k(),nIter())})
      PATable <- round(PATab[c("emRR","emBB","emBL")],digits=2) 
      names(PATable)<- c("Ridge Regression","BayesB","Bayes LASSO")
      PATable <- rbind.data.frame(names(PATable),PATable)
      colnames(PATable) <- c("Ridge Regression","BayesB","Bayes LASSO")
      rownames(PATable) <- c("Model","Accuracy")
      return(PATable[2,])
    }
    if(nSelTraits()>1){
      PATableComb <- c()
      for(nSelT in 1:nSelTraits()){
        i <- reactive(nSelT)
        PATab <- withProgress(message = 'Running CrossValidations', value = 0, {
          getemCVR(predictionData(),unlist(Trait())[i()],nTraits(),k(),nIter())})
        PATable <- round(PATab[c("emRR","emBB","emBL")],digits=2)
        #PATable2 <- c(unlist(Trait())[i()],PATable)
        PATableComb <- rbind(PATableComb,PATable)
        
      }
     
      #PATableComb <- rbind(c("Trait","Ridge Regression","BayesB","Bayes LASSO"),PATableComb)
      #colnames(PATableComb) <- rep("",ncol(PATableComb))
      
      colnames(PATableComb) <- c("Ridge Regression","BayesB","Bayes LASSO")
      rownames(PATableComb) <- unlist(Trait())
      return(PATableComb)
    }
    
  })

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
  
 
# 
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
  
 
# Start the background process and return rProcess
 
 rProcess6 <- eventReactive(input$CrossValidationME,{
   
    # Path for the temporary file
    temp_file6(tempfile())
    
    #ME_argList <- list(DT_1_Filt_List,genoDat_List,traits,KG,KE,CVMet,factVar,KMethod,FitEnvModels,fixedME,envVar,IDColsME,LocME,YrME)
    
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
  })


  
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
        return(paste(txt, "Cross-validation Completed", collapse = "\n"))
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
  
  
  output$emCVRMT <- renderTable({
    if(nSelTraits()>1){
      PATable <- cvrOutputListMT()
      #colnames(PATable) <- rep("",ncol(PATable))
      print.data.frame(as.data.frame(PATable))
    }
   },colnames=TRUE,rownames=TRUE) 
  
  
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
      write.csv(as.data.frame(cvrOutputListME()),"ME_CV_OutTable.csv")
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
    # if(isPhenoME() && (MECV_Out_List()[[CVMet()]]!=0)){
    # PATable <- MECV_Out_List()[[CVMet()]]
  },colnames=TRUE,rownames=TRUE)
  

# 
# output$emCVRME <- renderTable({
# 
  #   if(isPhenoME() & !is.null(MECV_Out())){
  #     PATable <- cvrOutputListME()
  #     print.data.frame(as.data.frame(PATable))
  #   }
  # },colnames=TRUE,rownames=TRUE)
  
  
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
  
  
  
  outputList <- eventReactive(input$RunPredictionsST,{ 
    if(nSelTraits()==1){
      outputDF <-  withProgress(message = 'Running Computations', value = 0, {
        getRankedPredictedValues(predictionData(),nTraits(),unlist(Trait()),GPModelST(),fixedX=Fixed(),fixedData=fixedData_List(),optTS())})
      return(outputDF)
    }else if(nSelTraits()>1){
      outputDFComb <- c()
      for(nSelT in 1:nSelTraits()){
        i <- reactive(nSelT)
        outputDF <-  withProgress(message = 'Running Computations', value = 0, {
           getRankedPredictedValues(predictionData(),nTraits(),unlist(Trait())[i()],GPModelST(),fixedX=Fixed(),fixedData=fixedData_List(),optTS())})
        outputDFComb <- cbind(outputDFComb,outputDF[,2])
      } 
      outputDFComb <- cbind(outputDF[,1],outputDFComb,outputDF[,3])
      colnames(outputDFComb) <- c("GermplasmID",unlist(Trait()),"UpperBound of Reliability")
      return(outputDFComb)
    }
    
  })
  
  
  ### GPModel MT
  
  GPModelMT <- reactive(input$GPModelMT)
  
  ## 
  
  
  outputListMT <- eventReactive(input$RunPredictionsMT,{ 
    
    withProgress(message = 'Running Computations', value = 0, {
      getRankedPredictedValuesMT(predictionData(),nTraits(),unlist(Trait()),GPModelMT(),optTS())
    })
  })
  
  ### GPModel ME
  
  #YearME,LocationME,fixedME
  
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
  
  KGMethod <- reactive({
    if(input$GKernelMet =="Linear"){"GB"}else if(input$GKernelMet =="Gaussian"){"GK"}
  })
  
  fitEnvCovs <- reactive(input$fitEnvCov)
   
  outputListSTME <- eventReactive(input$RunPredictionsME,{ 
    print("runMEIn")
    #browser()
    withProgress(message = 'Running Computations', value = 0, {
      if(fitEnvCovs() == FALSE){
        getMEPred(DT_Filt_List(),genoDat_List(),TraitME(),KG=NULL,KE=NULL,KMethod= KGMethod(),FitEnvModels = fitEnvCovs(),fixedME=fixME(),envVar=varEnv(),IDColsME = IDColsME(),LocME=LocationME(),YrME=YearME())
      }else{ 
        getMEPred(DT_Filt_List(),genoDat_List(),TraitME(),KG=NULL,KE=EnvK_Mod(),KMethod= KGMethod(),FitEnvModels = fitEnvCovs(),fixedME=fixME(),envVar=varEnv(),IDColsME = IDColsME(),LocME=LocationME(),YrME=YearME())
      }
    })
  })
  
  #outputListSTME <- reactive({ getMEPred(DT_Filt_List(),genoDat_List(),TraitME(),KG=NULL,KE=NULL,KMethod= KGMethod(),FitEnvModels = fitEnvCovs(),fixedME=fixME(),envVar=varEnv())})

  outputListMETab <- reactive({
    #browser()
    getCombinedTab(outputListSTME(),TraitME(),IDColsME(),IDColME(),fitEnvCovs())
  })


## Render Table for ST  
  
  output$Ranked_Lines_for_SelectionST <- renderTable({
    if(nSelTraits()>=1){
     as.data.frame(outputList()[1:20,])
    }  
  })
  
  ## Render Table for MT
  
  output$Ranked_Lines_for_SelectionMT <- renderTable({ if(nSelTraits()>1){
    
    (outputListMT()[1:20,])
  }
  }) 
  
  #%>% bindCache(outputListMT(),processedData(),nTraits(),unlist(Trait()),GPModelMT()) %>% bindEvent()
  
### 
  
 ## Render Table for ME
 # "Homogeneous Variance (G+E+GxE)","Heterogeneous Variance (G+E+GxEi)")
  
  output$Ranked_Lines_for_SelectionSTME <- renderTable({ 
    if(GPModelME()=="Main Effect (G+E)"){
        (outputListSTME()[[1]][[3]][[1]][1:20,])
    }else if(GPModelME()=="Homogeneous Variance (G+E+GxE)"){
        (outputListSTME()[[1]][[3]][[2]][1:20,])
      
    }else if(GPModelME()=="Heterogeneous Variance (G+E+GxEi)"){
        (outputListSTME()[[1]][[3]][[3]][1:20,])
    }
  })

  
  
  output$scatter <- renderUI({ 
    if(nSelTraits()==1){
      plotOutput("plots")
    }
  }) 
  output$scatterMT <- renderUI({  
    if(nSelTraits()>1){
      plotOut_List <- lapply(1:nSelTraits(),function(x) { 
        i <- reactive(x)
        plotOutput(noquote(paste("'","plot",i(),"'",sep="")))})
      do.call(tagList,plotOut_List)
    }
  })
  
  output$plots <- renderPlot({
    if(nSelTraits()==1){
      plot.default(outputList()[,2],outputList()[,3],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main=paste("Upper bound of Reliability vs Predicted Values for ",Trait(),sep=""),font.lab=2,cex.lab=1.25)
    } 
  })
 
  observeEvent(input$RunPredictionsST,{
    if(nSelTraits()>1){
      lapply(1:nSelTraits(),function(i){
        loc_i <- reactive(i)
        output[[noquote(paste("'","plot",loc_i(),"'",sep=""))]] <- renderPlot({
          plot.default(outputList()[,(loc_i()+1)],outputList()[,(nSelTraits()+2)],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main=paste("Upper Bound of Reliability vs Predicted Values for ",Trait()[loc_i()],sep=""),font.lab=2,cex.lab=1.25)
        })
        
      })
    }
  }) 
  
  
  #  
  # local({
  
  #printPlots2 <- reactive({  
  observeEvent(input$RunPredictionsMT,{ 
    
    lapply(1:nSelTraits(),function(i){
      loc_i <- reactive(i)
      output[[noquote(paste("'","plot",loc_i(),"'",sep=""))]] <- renderPlot({
        plot.default(outputListMT()[,(loc_i()+1)],outputListMT()[,(nSelTraits()+2)],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main=paste("Upper bound of Reliability vs Predicted Values for ",Trait()[loc_i()],sep=""),font.lab=2,cex.lab=1.25)
      })
      
    })
    
  })
  
  # XLim <- c(min(c(unlist(outputListMT()[,2:3])))-2,max(c(unlist(outputListMT()[,2:3])))+2)
  # YLim <- c(min(c(unlist(outputListMT()[,4]))),max(c(unlist(outputListMT()[,4]))))
  # plot.default(outputListMT()[,2],outputListMT()[,4],type="p",xlim=XLim,ylim=YLim,xlab=paste("Predicted Value of ",paste(unlist(Trait()),collapse=" and "),sep=""),ylab="Upper Bound of Reliability",main="Upper bound of Reliability vs Predicted Values")
  # points(outputListMT()[,3],outputListMT()[,4],type="p",col="red")
  # legend(XLim[2]-2,YLim[2],legend=c(unlist(Trait())),xpd=TRUE,lty=c(1,1),lwd=3,col=c("black","red"),cex=0.7)
  #          }
  # })
  
  switch_page <- function(i) {
    sel <- reactive(i)
    PageID <- c("Home","Load Data","Set Target & Trait","Optimize Training Population",
                "Cross Validations","Genomic Prediction","Visualize Predictions")
    updateTabsetPanel(inputId = "inData", selected = PageID[sel()])
  }
  
  
  ########
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
      
      write.csv(as.data.frame(outputList()), file, row.names = FALSE)
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
  
  
  
  # session$allowReconnect("force")
  # options(shiny.launch.browser=FALSE)
  
}

shinyApp(server=server,ui=ui)





## 
# observeEvent(input$infileBLUEs, {
#   updateTabsetPanel(session,"inData",
#                     output$PhenoTable <- renderDataTable({as.data.frame(Pheno()[1:5,])})) }
# )
# 



# observeEvent(input$infileBLUEs,{cat(paste("Phenotype Table with ",nrow(Pheno())," lines and ",ncol(Pheno())-1," traits\n",sep=""))})



#Meta 
# MetaTab <- reactive({
# 

#metaFile <- input$infileMetaTable
#   
#   ext <- tools::file_ext(metaFile$datapath)
#   req(metaFile)
#   validate(need(ext == "csv", "Please upload a csv file"))
#   #NUST_Meta_Table <- 
#   read.csv(metaFile$datapath)
#   
# })


#fluidRow(
#column(3),column(width=9,tags$h4(tags$strong(textOutput("boxHeader"))))),
#tags$br(),
#plotOutput("plotBox")

#  
# withProgress(message = 'Processing Data', value = 0, {
#  NUST_Data_Table_Num_Filt <- getProcessedData(Geno(),Pheno(),MetaTab())
#     
#     
#   })

# fileInput("infileMetaTable", "Choose Meta Information File (CSV)", accept = ".csv"),
# checkboxInput("header", "Header", TRUE)

# metaHead <- eventReactive(input$infileMetaTable,{ paste("Table with meta information on ",nrow(MetaTab())," lines",sep="")})
# output$MetaHeader <- renderText({metaHead()})
# output$MetaTable <-  renderTable({as.data.frame(MetaTab()[1:5,1:5])})
# 
# tags$h3(textOutput("MetaHeader")),
# tableOutput("MetaTable"),

# output$Ranked_Lines_for_Selection <- renderDataTable({cbind.data.frame(rownames(outputList())[1:50],outputList()[1:50,])})
# output$scatter <- renderPlot({plot.default(outputList()[,1],outputList()[,2],type="p",xlab="Observed Value",ylab='Predicted Value',main="Predicted vs Observed Values")})


#### Convert VCF to Dataframe
######################################################################
###                         VCFToSimple                            ###
###   a function to generate simple format from vcf                ###
######################################################################
# function to convert vcf into simple data format
# need vcfR library and dplyr library
# input: vcf file path and name



# fluidRow(
#    column(3),column(width=5,tags$li(" A '.csv' file containing information on the target data, whose genetic values are to be predicted using the genomic prediction model. 
#                        The first column in the target data information file should have the line IDs."))),
#      
#  fluidRow(
#    column(3),column(width=5,fileInput("infileTargetTable", "Choose Target Data Information File (CSV)", accept = ".csv"),
#      checkboxInput("header", "Header", TRUE)
#    )),
#   fluidRow(
#        column(3),column(width=5,div(
#          tags$h4(tags$strong(textOutput("TargetHeader"))),
#          tags$h6(tableOutput("TargetTable")))))


#  fluidRow(
#   column(6),column(width=5,tags$h6(tableOutput("PhenoTable")))),
#  tags$br(),




# 
# fluidRow(column(2),column(width=5,checkboxInput("header", "Header", TRUE)),
#   column(6),column(width=5,checkboxInput("header", "Header", TRUE))),
# 
# fluidRow(
#   column(6),column(width=5,div(
#     tags$h4(tags$strong(textOutput("GenoHeader"))),
#     tags$h6(tableOutput("GenoTable"))))),
#  tags$br(),
#  tags$br(),
#  tags$br(),
#  fluidRow(
#  column(6),column(width=5,tags$ul(tags$li("A '.csv' file containing phenotypic information. Phenotype file should have the line ID in the first column and phenotypic trait values after that with the trait ID as column name.")))),
#  tags$br(),
#  fluidRow(
#  column(6),column(width=5,fileInput("infileBLUEs", "Choose Phenotype File (CSV)", accept = ".csv"),
#    checkboxInput("header", "Header", TRUE))),
#      
#  tags$br(),
#  fluidRow(
#   column(6),column(width=5,tags$h4(tags$strong(textOutput("PhenoHeader"))))),
#  fluidRow(
#   column(6),column(width=5,tags$h6(tableOutput("PhenoTable")))),
#  tags$br(),
#  
#  fluidRow(
#    column(3),column(width=5,tags$li(" A '.csv' file containing information on the target data, whose genetic values are to be predicted using the genomic prediction model. 
#                        The first column in the target data information file should have the line IDs."))),
#      
#  fluidRow(
#    column(3),column(width=5,fileInput("infileTargetTable", "Choose Target Data Information File (CSV)", accept = ".csv"),
#      checkboxInput("header", "Header", TRUE)
#    )),
#   fluidRow(
#        column(3),column(width=5,div(
#          tags$h4(tags$strong(textOutput("TargetHeader"))),
#          tags$h6(tableOutput("TargetTable")))))


# if(nSelTraits()>1){
#   PATable <- cvrOutputListMT()
#   colnames(PATable) <- rep("",ncol(PATable))
#   print.data.frame(as.data.frame(PATable))
# } 
# 
#})
# if(nSelTraits()==1){
#   cvrOutputListST <- reactive({
#     withProgress(message = 'Running CrossValidations', value = 0, {
#       getemCVR(processedData(),unlist(Trait()),nTraits(),k(),nIter())})
#     
#   }) 
#   PATable <- rbind(c("RR","BB","BL"),round(cvrOutputListST()[c("emRR","emBB","emBL")],digits=2)) 
#   rownames(PATable)<- c("Prediction Model","Prediction Accuracy")
#   colnames(PATable) <- rep("",ncol(PATable))
#   PATable
#   
#   
# }
# }) # %>% bindEvent(input$CrossValidationST)

#%>% bindCache(processedData(),unlist(Trait()),nTraits(),k(),nIter()) 


#%>% bindCache(cvrOutputListMT(),processedData(),Trait(),nTraits(),k(),nIter(),nSelTraits())

## One plot in renderUI
# if(nSelTraits()==1){
#   plotOutput("plots")
# } 
## Multiple plots in renderUI
# if(nSelTraits()>1){
#   plotOutput("plots1")
#   plotOutput("plots2")
#   # plotOut_List <- lapply(1:nSelTraits(),function(x) {plotname <- paste("plot",x,sep="") 
#   # plotOutput(plotname)})
#   # do.call(tagList,plotOut_List)
# }



## Multiple plots in renderPlot      

# for(i in 1:nSelTraits()){ 
#   local({
#     loc_i <- i
#     plotName <- paste("plot",loc_i,sep="")
#     output[[plotName]] <- renderPlot({
#         plot.default(outputListMT()[,2],outputListMT()[,(nSelTraits()+2)],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main="Upper bound of Reliability vs Predicted Values")
#     })
#   })
# }

#    reactive(
# ### One plot in renderPlot
#      if(nSelTraits()==1){
#         output$plots <- renderPlot({
#          plot.default(outputList()[,2],outputList()[,3],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main="Upper bound of Reliability vs Predicted Values")
#         })
#      }
#    )



# output$Ranked_Lines_for_Selection <- renderDataTable({as.data.frame(outputListMT()[1:50,])})
# output$scatter <- renderPlot({plot.default(outputListMT()[,2],outputListMT()[,(nTraits()+2)],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main="Upper bound of Reliability vs Predicted Values")})
# 


#   outputListMT <- reactive({
# 
#     withProgress(message = 'Running Computations', value = 0, {
#       getRankedPredictedValuesMT(processedData(),nTraits(),unlist(Trait()),GPModelMT(),optTS=NULL)
#     })
#    })  %>% bindCache(processedData(),Trait(),nTraits(),GPModelMT()) %>% bindEvent(input$RunPredictionsMT)



# output$Ranked_Lines_for_Selection <- renderDataTable({as.data.frame(outputList()[1:50,])})
# output$scatter <- renderPlot({ plot.default(outputList()[,2],outputList()[,3],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main="Upper bound of Reliability vs Predicted Values")})




#%>% bindCache(TSOptOutputList(),processedData(),unlist(Trait()),nTraits(),noCandidates(),nTrainToSelect()) %>% bindEvent(input$noToSelect,input$Optimize,input$trait,input$infileVCF,input$infileTargetTable,input$infileBLUEs)

# TSTable_Mod <- cbind(c("STPGA","Random"),apply(TSTable,2,function(x) round(as.numeric(x),digits=2)))
#colnames(TSTable_Mod) <- c("TS","RR","BB","BL")
# setGAPar <- reactive(input$setGA)


# output$TSGA <- renderUI({ 
#   if(setGAPar() ==TRUE){} 
#   
#   if(setGAPar() ==FALSE){
#   
#     numericInput(inputId="noPop","GA Population Size",value = 100,min =1,max=1000)
#   
#     numericInput(inputId="noElite","Elite Population Size",value =10,min =1, max=50)
#     numericInput(inputId="mutProb","Mutation Probability",value =0.5,min =0.1, max=1)
#     
#     numericInput(inputId="mutInt","Mutation Intensity",value =1,min =0, max=1)
#    
#     numericInput(inputId="nIter","# Iterations",value =100,min =50, max=1000)
#  
#     selectInput(inputId="tabu","Tabu Search",choices=c("TRUE","FALSE"),selected="TRUE")
#     
#     numericInput(inputId="tabumemsize","Tabu Memory Size",value =1,min =1, max=10)
#    
#     numericInput(inputId="lambda","Shrinkage Parameter (lambda)",value =1e-6,min =1e-9, max=1)
#     
#     numericInput(inputId="mcCores","#Cores",value =10,min =1, max=20)
#    
#     
#   }
#   
# })
# 
# Mssg1CVR <- eventReactive(input$CrossValidationST,{
# if(nSelTraits>1){
#    paste("Performing CrossValidation for the first trait")
#   
#  }
# }) 

# cvrOutputListST <- eventReactive(input$CrossValidationST,{
# 
#       withProgress(message = 'Running CrossValidations', value = 0, {
#       getemCVR(processedData(),unlist(Trait()),nTraits(),k(),nIter())})
# 
# })
# 

# %>% bindEvent(input$CrossValidationMT)

#nSelTraits <- eventReactive(input$CrossValidationST,{(1)})


# observeEvent(input$Home_Data,{ switch_page(2)})
# observeEvent(input$Data_Home, switch_page(1))
# 
# observeEvent(input$Data_Trait, switch_page(3))
# observeEvent(input$Trait_Data, switch_page(2))
# 
# observeEvent(input$Trait_TS, switch_page(4))
# observeEvent(input$TS_Trait, switch_page(3))
# 
# observeEvent(input$TS_CVR, switch_page(5))
# observeEvent(input$CVR_TS, switch_page(4))
# 
# observeEvent(input$CVR_GP, switch_page(6))
# observeEvent(input$GP_CVR,switch_page(5))
# 
# observeEvent(input$GP_VP, switch_page(7))
# observeEvent(input$VP_GP, switch_page(6))

# output$ExportOut <- downloadHandler(
#   filename = function() {
#     "Predictions.csv"
#   },
#   content = function(file) {
#     if(nSelTraits()==1){
#       write.csv(as.data.frame(outputList()), file, row.names = FALSE)
#     }else if(nSelTraits()>1){
#       write.csv(as.data.frame(outputListMT()), file, row.names = FALSE)
#     }
#   }
# )
