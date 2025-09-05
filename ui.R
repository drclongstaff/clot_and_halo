library(shiny)
library(plotly)
library(tibble)
library(purrr)
library(dplyr)
library(bslib)
library(shinyWidgets)

#locations of the functions
source("./Functions/Load_file.R")
source("./Functions/Splin_calc.R")
source("./Functions/Data_proc.R")
#source("./Functions/Plots_all.R")
#source("./Functions/Plot_one.R")
source("./Functions/Making_plots.R")

#√ersion
vers <- "version 0.3.4,"
ThisApp <- "Clot or Halo"

fluidPage(
  includeCSS("./www/styles3.css"), #make a few changes to the colours and fonts
  #theme = bs_theme(brand = "./www/_brand.yml"),
  
  titlePanel (h2(ThisApp,  align = "center")),
  
  sidebarLayout(#position="right",
    
      #sidebarPanel(style="background: #DCF0DC", 
      sidebarPanel(style="background: #C7E3F0",              
            #tags$style(".well {background-color:#F6FFF6;}"),#colour for sidebar
           
      conditionalPanel(condition="input.tabselected==1", #conditional panels are different in each tab
                #wellPanel(style="background: #F6FFF6",
                tags$h4("Load your data"),        
            fluidRow(
                      #Will call a function in the server to detect and load the user's data file
                      column(6,fileInput("data", label = "Select data", accept = c(".csv", ".txt", ".xlsx"))),
        
                      column(5,numericInput("sheet", "Excel sheet", value = 1, min = 1, step = 1) )
               # )
                    ),
      
            
            wellPanel(style="background: #EFF5F8",
            #a section to adjust analysis settings for % lysis and threshold
            tags$h4("Analysis settings"),
     
            fluidRow(

                      column(4,  offset= 0, numericInput("ini", 
                                           label = h5("%change"),step = 5,
                                           value = 50)),
        
                      column(8, offset = 0, numericInput("thresh", label = h5("Interpolation threshold"), value = 0.05, min = 0, step = 0.05))
                      )
                    ),
            
            #Baseline settings
            #wellPanel(style="background: #F6FFF6",
            tags$h4("Baseline options"),
      
            fluidRow(
             
            radioGroupButtons(
                inputId = "abini",
                label = NULL,
                choices = c("global zero", 
                            "nth absorbance", "min+offset"),
                selected="global zero",
                status = "danger",
                individual = TRUE,
                checkIcon = list(
                  yes = tags$i(class = "fa fa-check-square", 
                               style = "color: steelblue"),
                  no = tags$i(class = "fa fa-square-o", 
                              style = "color: steelblue"))
              )
              
            ),
            
            
            fluidRow(
        
                      column(4, numericInput("back",
                               label = h5("global zero"), step = 0.005, value = 0.042)),
        
                      column(4, numericInput("arow",
                               label = h5("nth point"), value = 1)),
        
                      column(4, numericInput("off", label = h5("offset"), value = 0, min = 0, step = 0.005) ),
                     
                            ),
                      
                      tags$h4("Maximum absorbance options"),
                      
                      fluidRow(
                        
                        radioGroupButtons(
                          inputId="abfin", label = NULL, 
                          choices = c(  "each curve max", "global max" ), 
                          selected = "each curve max", 
                          status = "warning",
                          individual = TRUE,
                          checkIcon = list(
                            yes = tags$i(class = "fa fa-check-square", 
                                         style = "color: steelblue"),
                            no = tags$i(class = "fa fa-square-o", 
                                        style = "color: steelblue"))
                          )
                          
                         
                        ),
                        fluidRow(column(6, offset=4, numericInput("finabs",
                                               label = h5("global max value"), step = 0.05, value = 0.4))
                                 ),
                       
      
              
      
              #select raw data or spline fitting
            
            #tags$br(),
            wellPanel(style="background: #EFF5F8",  
            tags$h4("Use raw data or spline fit additional points"),
      
      
              fluidRow(
                
                      column(6, radioButtons("spln", label = NULL, inline = FALSE, 
                               choices = c("raw", "spline"))),
        
                      ),
            
           
            
      
              tags$h5("Values for spline fit"),
      
              fluidRow(

                      column(4, offset = 0, numericInput("npoints", "n points", min=2, step = 5, value = 100)),
        
                      column(4, numericInput("zero", "start", min=0, step = 5, value = 0)),
        
                      column(4, numericInput("trunc", "truncate", min = 0, step = 10, value = 500))
                      
                      )  
            )
          ),  #end of conditional panel 1 
    
      
      
      conditionalPanel(condition="input.tabselected==2",
                       wellPanel(style="background: #EFF5F8",
                                
                       
                          #This is the output to select the well which comes from the server
                          uiOutput("what"),
                                  
                         
                       )
                     
                      ),  #end of panel 2  
      
      
        conditionalPanel(condition="input.tabselected==3", 
                        # wellPanel(style="background: #F6FFF6",
                                   
                           helpText(h4("All results for")),
                          #displays the name of the file being analysed 
                          uiOutput("which2")
                          
                                #   )
                         ),  #end of panel 3 
    
    
      conditionalPanel(condition="input.tabselected==4", 
                       #wellPanel(style="background: #F6FFF6",  
                     
                        
                                  
                         helpText(h4("Raw data")),
                        
                          #displays the name of the raw data file
                          uiOutput("whichraw")
                       
                      
                        
                         ),  #end of panel 4 
    
    
        conditionalPanel(condition="input.tabselected==5", 
                        
                         wellPanel(style="background: #EFF5F8",
                        #This is the explore panel 
                         helpText(h3("Heat map or scatter plot")),
                     
                         #which plot to show
                         radioButtons(inputId = "heat", label = "Choose plot", 
                                      choices = c("heat", "scatter"), selected = "heat")
                         ),
                         
                        wellPanel(style="background: #EFF5F8",
                        helpText(h4("Select heat plot output")),
                        #selects what set of results to be displayed in the heat ma[] 
                        uiOutput("whatheat")
                        ),
                         helpText(h4("")),
                         
                        wellPanel(style="background: #EFF5F8", 
                        
                        helpText(h4("Select scatter plot output")),
                        #select the variables to go in the scatter plot 
                          uiOutput("whatx"), 
                         
                          uiOutput("whaty")
                          )
                        ),    #end of panel 5 
      
      
        conditionalPanel(condition="input.tabselected==6", 
                       #  wellPanel(style="background: #F6FFF6", 
                      #the matrix of settings to help reproducibility
                        helpText(h4("Settings")),
                        helpText(h5("This table records analysis settings to aid quality control and reproducibility "))
                       #  )
                        ),  #end of panel 6
    
    
        conditionalPanel(condition="input.tabselected==7",
                      
                         #some blurb and promotional stuff
                         wellPanel(style="background: #EFF5F8", 
                        helpText(h5("Please cite this reference if you find the app useful:")),
                        #tags$br(),
                        helpText(h5("Longstaff C, Development of a Shiny app tool to simplify and standardize the analysis 
                                   of hemostasis assay data:", 
                                   tags$a(href="https://doi.org/10.1111/jth.13656","J Thromb Haemost, 15: 1044-6, 2017")
                                   #communication from the SSC of the ISTH, J Thromb Haemost, 15: 1044-6, 2017))
                         ))),
                        
                        tags$br(),
                        wellPanel(style="background: #EFF5F8", 
                       tags$i("Please contact me with issues relating to:"),
                       helpText(h5(ThisApp, vers,
                                   " last accessed", Sys.Date()),
                                tags$a(href="mailto: drclongstaff@gmail.com", "drclongstaff@gmail.com")),
                       
                       helpText(h5("You may be interested in a simple version of this clotting app and other apps that 
                                   are useful in hemostasis assays. "),
                       
                       tags$a(href="https://drclongstaff.github.io/shiny-clots/", "Links to other apps and help notes")),
                              #),
                       tags$br(),
                       tags$a(href="https://www.youtube.com/@colinlongstaff7270", "Youtube channel of help videos"))
                      
                        ) #end of conditional panel 7
     
    
          ), #end of sidebar panels
    
    
  mainPanel( 
      tabsetPanel(#type="tab",
       # tags$br(),    
            tabPanel("Plots", value = 1,
                      tags$br(), 
                        
  
                              fluidRow(
                          
                                      column(5, helpText(h4(uiOutput("which")))),      
                                      column(4, card_title(h4("Plot all curves"))),
                                      column(3, card_header(numericInput("numrows",label = "Plot number of rows", value = 8, min = 1, step = 1),
                                      class = "bg-light")),
                                      
                                      ),
                          
                     card(
                       #can make plots expandable 
                       #height = 550,
                       full_screen = TRUE, #is expandable         
                     
                     card_body(
                            
                                      plotOutput(outputId = "plot")

                                       ) 
                            ),     #end of card
                         
                        #select which results to show in the table
                        radioButtons(inputId="tabRes", label = h4("Select a parameter", align = "center"), 
                                     choices = c("Column names"=1, "Chosen zero"=2, "Time to % clotting   "=3, "Reading at % clotting"=4, 
                                                          "Reading at peak      "=5, "Reading peak-zero"=6,
                                                          "Time to peak from zero"=7, "Time of max rate"=8
                                                                                                                        
                                                ), selected = 3, inline = TRUE, width = "100%"),
                         #text3 displayes the results selected
                                     h5(textOutput("text3")),

                        #this is the results table generated in the server
                         tableOutput("resultsTable"), align = "center"
                      
                 
                    ),  #the end of the main page in tab 1
                   
                
                tabPanel("Curve", value = 2,
                         
                         card(
                           #can make plots expandable 
                           #height = 550,
                           full_screen = TRUE, #is expandable         
                           
                           card_body(
                         
                         #the single plot of the selected well
                         plotOutput(outputId = "myplot")
                           )
                         ),
                         #the results table
                            fluidRow(
                                  column(5, offset=4,
                                  h4("Results table"), tableOutput("curveTable"))
                                    )
                      
                           ), #end of main panel in tab 2
                 
                  #tabPanel("All Results", value = 3, dataTableOutput("contents"), align = "centre"),
                  tabPanel("All Results", value = 3, DT::DTOutput("contents"), align = "centre"),
                
                  tabPanel("Raw Data", value = 4, tableOutput("raw")),
                
                  tabPanel("Explore", value = 5,
                         
                         #Using plotly
                         plotlyOutput(outputId = "exploreplot")
  
                           ),  #end of panel 5
                
                  tabPanel("Settings", value = 6,  
                         
                           #headings for the table showing data analysed
                           fluidRow(
                                      column(4,helpText(h5("file name, size, Excel sheet"))),
                                      column(6, h5(uiOutput("which1"))),
                                      column(2, h5(uiOutput("whichsheet")))
                                  ),
                            
                           tableOutput("settings")
                           
                           ),  #end of panel 6
                
                    tabPanel("Help", value = 7,
                         tags$h4("Load data"),
                         tags$blockquote(h5("►An example set of data is provided and loads when the app is started. Different settings can be explored with these data ",
                                            tags$br(),
                                            "►Load your own data and the file type will be detected (csv, txt, xlsx are accepted and an Excel sheet number can be specified)",
                                            tags$br(),
                                            "►Data should be supplied as one column of time and any number of columns of absorbances as shown in the example data provided",
                                            tags$br(),
                                            "►Data files should not contain any gaps or spaces between characters and avoid unusual characters",
                                            tags$br(),
                                            "►The default setting for clotting is 50% but can be adjusted",
                                            tags$br(),
                                            "►The threshold value determines a lower value below which interpolation does not take place and can be increased if the analysis fails ",
                                            tags$br(),
                                            "►Noisy data or empty cells can cause the program to crash but more robust (though less accurate) curve analysis can be achieved by using a high threshold absorbance (e.g. 1)",
                                            tags$br(),
                                            "►Accuracy may be improved, if interpolation is not used, by spline fitting and adding extra points (see below) "
                                            
                                        )),
                         
                         tags$h4("Baseline options"),                 
                         tags$blockquote(h5 ("►These selections adjust the baseline for the start of clotting",
                                             tags$br(),
                                             "►'Global zero' sets a baseline absorbance for all the curves to the same selected value",
                                             tags$br(),
                                             "►'nth point' subtracts the absorbance value from a selected nth point from each curve (e.g first point or a point after lysis)",
                                             tags$br(),
                                             "►'Min abs' finds and subtracts the lowest absorbance for each curve. An offset value can be added"
                                        )),
                         tags$h4("Maximum options"),                 
                         tags$blockquote(h5 ("►You can specify a global maximum or use the maximum clotting absorbance for each curve",
                                             tags$br(),
                                             "►Global maximum may be useful if not all curves clot fully but you know what the final absorbance should be",
                                             tags$br(),
                                             "►Global maximum is useful for Halo assays where some wells are used to provide a maximum absorbance for complete lysis 
                                             using blood that is not clotted"
                                             
                         )),
                         
                         tags$h4("Raw data or spline fitted data"),                 
                         tags$blockquote(h5 ("►If data points are very sparse or you are using a high 'Threshold' value, it may be helpful to generate extra points using the spline fitting option ",
                                             tags$br(),
                                             "►Options for curve start, number of additional points and data trunction are provided for spline curve generation",
                                             tags$br(),
                                             "►If you have good data with reasonable density of points simply use the defaults of raw data and interpolation",
                                             tags$br(),
                                             "►Input data (raw or with extra fitted points, as selected), are shown in the plots in the Curve tab and data in the Raw Data tab"
                                        )),
                         
                         tags$h4("Graph options"),                 
                         tags$blockquote(h5 ("►Displayed graphs and results table are presented in an adjustable arrangement",
                                             tags$br(),
                                             "►Input the number of rows to display the plots and results as required (e.g. 8 rows for a full 96 well plate)",
                                             tags$br(),
                                             "►Individual curves can be scrutinised on the next tab (Curve)",
                                             tags$br(),
                                             "►These options do not affect data or calculated results"
                                        )),
                         
                         
                         tags$h4("Select a parameter"),                 
                         tags$blockquote(h5 ("►Select what should be shown in the results table",
                                             tags$br(),
                                             "►The results table is organised like the graphical output above",
                                             tags$br(),
                                             "►Values can be copied and pasted for further analysis",
                                             tags$br()
                                        )),
                        
                         tags$h4("Other tabs"),                 
                         tags$blockquote(h5 ("►'Curve' focuses on an individual well",
                                             tags$br(),
                                             "►'All Results' is a table showing all results for each curve",
                                             tags$br(),
                                             "►'Raw Data' is the time and absorbance data set.  Time interval changes if spline fitted curves are used",
                                             tags$br(),
                                             "►'Explore' presents graphical outputs of analysed results",
                                             tags$br(),
                                             "►'Settings' generates a table of critical setings that can be recorded to help reproducibility",
                                             tags$br(),
                                             "►Code files and detailed help notes are available in a github repository", 
                                             #tags$a (href="https://github.com/drclongstaff/Clotlysis_2019/blob/master/Clotlysis_CL_help_notes_2019.pdf", "Here")
                                             tags$a (href="https://github.com/drclongstaff/clot_or_Halo/", "Here")
                                        )),
                         tags$h4("Halo assays"),                 
                         tags$blockquote(h5 ("►Settings on the app allow it to be used for halo assay data analysis [Bonnard et al, Sci Rep 7, 2346, 2017]",
                                             tags$br(),
                                             "►Typically in this assay wells are included that give estimates for baseline and global maximum, which can be input into the app",
                                             tags$br(),
                                             "►The high absorbance in the global max wells may mean the threshold needs to be increased and/or spline fitted data used",
                                             tags$br(),
                                             "►The resulting lysis curves resemble clotting curves and results can be reported as time to 50% lysis (equivalent to 50% clotting in the app) ",
                                             tags$br()
                                            
                                        )),
                  
                  #a png of the first page plots (is it necessary?)
                  #tags$img(src="PlotTable.png", width=600, height=500),
                  #tags$img(src="screenCapS6.png", width=600, height=700)
                  
                        tags$h5("Changes in version 0.3"),                 
                        tags$blockquote(h5 ("►Several changes have been made to the underlying code and appearance of Clot or Halo",
                                      tags$br(),
                                      "►The code dealing with curve analysis should be more robust and efficient",
                                      tags$br(),
                                      "►A new function is provided for reading user data which identifies the correct format of csv, txt, xlsx",
                                      tags$br(),
                                      "►The default settings are to use interpolation of curves, but the option to avoid interpolation is set by increasing the 'Threshold' value",
                                      tags$br(),
                                      "►Times at maximum rate increase are provided and shown on the clotting traces as a magenta line",
                                      tags$br(),
                                      "►These results are calculated from smoothed first derivative curves",
                                      tags$br(),
                                      tags$br(),
                                      tags$i("Please contact me with issues relating to:"),
                                      helpText(h5("Clot or Halo", vers,
                                                  " last accessed", Sys.Date()),
                                               tags$a(href="mailto: drclongstaff@gmail.com", "drclongstaff@gmail.com")),
                                      tags$a(href="https://drclongstaff.github.io/shiny-clots/", "Links to other apps")
                                      
                                      
                                  )),
                  
                  ),  #end of panel 7
                
               id = "tabselected"
                
            )
         )
      )
  )


