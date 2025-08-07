
options(shiny.maxRequestSize=30*1024^2)

function(input, output){
 
  #Load data using load_file function 
  readData <- reactive({
    
    inputFile <- input$data
    if (is.null(inputFile)) 
      read.csv("./Data/selected ancrod data.csv")
    else(
      load_file(input$data$name, input$data$datapath, input$sheet) %>% 
        as.data.frame()
               )
         })
  
  
  
  #displays the Excel sheet in the settings table
    whichsheet <- reactive({
       input$sheet
        })
  
    output$whichsheet <- renderUI({
       whichsheet()
        })
 
  #processed data dealing with background subtraction
    procdat0<- reactive({
    
                switch(input$abini,
                       "global zero"=readData(),
                       #"global zero"=as_tibble(cbind("Time"=readData()[[1]], readData()[,-1]-input$back)),
                       "nth absorbance"=map_df(readData()[,-1], ~BaselineNP(.x, input$arow)) %>% #BaselineNP is a function defined above
                         add_column("Time"=readData()[[1]], .before = TRUE),
                       "min+offset"=map_df(readData()[,-1], ~BaselineOff(.x, input$off)) %>% 
                         add_column("Time"=readData()[[1]], .before = TRUE)#,
                       #"raw data"=readData()
                      )
          })
  
  #For spline fitting and adjusting number of points
    procdat <- reactive({
                switch(input$spln,
                       "raw"= procdat0(),
                       #"spline"=procdat0()
                       "spline"=allSpline(NULL, input$npoints, input$zero, input$trunc, procdat0())
                        )
      #use clipr for pasting to clipboard locally only
      #clipr::write_clip(procdat0)
          })
  
      output$text3<-renderText({
          #text3 is a header for the results table to show what's been displayed
          paste(names(TabRes()[as.numeric(input$tabRes)]),"for ",input$ini, "%", "of maximum")
          })
  
      #Whichfile is for identifying the file and size loaded to include in the settings table
      whichfile <- reactive({
        if (is.null(input$data)) 
          return("Data/selected ancrod data.csv")
        else(input$data[,1:2])
           })
      
      #whichfile 2 is only the filename when needed
      whichfile2 <- reactive({
        inputFile <- input$data
        if (is.null(inputFile)) 
          return(file.path("Data/selected ancrod data.csv"))
        else(input$data$name)
           })
      
      #The following operations return file names to go into the ui    
      output$which <- renderUI({
        whichfile2()
          })
      
      
      output$which1 <- renderUI({
                whichfile()
          })
  
  
      output$which2 <- renderUI({
                whichfile2()
          })
  
      output$whichraw <- renderUI({
                whichfile2()
          })
  
      output$setfile <- renderUI({
              whichfile()
          })
  
      #This is to collect the column names without the first column of Time
      var<- reactive({
              colnames(readData()[,-1])
          })
  
      #This is for tab 2 to select the individual curve
      output$what<-renderUI({
              selectInput("colmnames",
                label= h5("Select a column of absorbance data"),
                choices = var())
           }) 
  
      #generates the table of all results with some rounding and re-ordering
      #output$contents<-renderDataTable({
      output$contents<-DT::renderDT({
    
              TabRes() %>% 
                    mutate(across(where(is.numeric), \(x) round(x, digits=3)))
    
            })
  
  
          ###Table of results
      #There are some functions to do the calculations
      #Then purrr makes the final table of results
      
      #first min and max for the basics of the curve analysis
      MaxandMin<-function(m, Time, ini, thresh){
    
            minAbs <- switch(input$abini,
                     "global zero"=input$back,
                     "nth absorbance"=m[input$arow],
                     "min+offset"=min(m, na.rm=TRUE)#+input$off
                          )
            switch(input$abfin,
                   "global max"=maxA<-input$finabs,
                   "each curve max"=maxA<-max(m, na.rm = TRUE))
            if(maxA>max(m, na.rm = TRUE)){pointmax<-length(m)}  else{pointmax<-which(abs(m-maxA)==min(abs(m-maxA)))[1]}
            #pointmax<-which.max(m)
            
            firstA=m[1]

            #maxA <- max(m, na.rm = TRUE)
            
            maxT <- Time[which.max(m)]#,
            changeA=maxA-minAbs
            
        #Here list the variables you want to collect
        Aminmax<-c(firstA, minAbs, maxA, changeA, maxT, pointmax)
    
                    }
  
      #uppity analysis the clotting part
      uppity<-function(u, Time, ini, thresh){
    
          minAbs <- switch(input$abini,
                     "global zero"=input$back,
                     "nth absorbance"=u[input$arow],
                     "min+offset"=min(u, na.rm=TRUE)+input$off
                        )
          
          switch(input$abfin,
                 "global max"=maxAbs<-input$finabs,
                 "each curve max"=maxAbs<-max(u, na.rm = TRUE))
          if(maxAbs>max(u, na.rm = TRUE)){pointmax<-length(u)}  else{pointmax<-which(abs(u-maxAbs)==min(abs(u-maxAbs)))[1]}
            
            #maxAbs <- max(u, na.rm = TRUE)
            #pointmax<-which.max(u)
            
            upTime<-Time[c(1:pointmax)] #vector of time to max
            upAbs<-u[c(1:pointmax)]  #vector of absorbances to max
            pcChange<-ini*(maxAbs-minAbs)+minAbs#this minAbs is determined in shiny, may be set or calculated
            #need to account for flat response
            if(pcChange>max(u, na.rm=TRUE)){startPoint <- length(u)} else{startPoint<-which(abs(upAbs-pcChange)==min(abs(upAbs-pcChange)))[1]}
            #startPoint<-which(abs(upAbs-pcChange)==min(abs(upAbs-pcChange)))[1]
            #StartAbs is fitted if abs > threshold, otherwise is closest point
            #This prevents crashing if there are blank wells
            ifelse(max(u)-min(u)<thresh,
                   startAbs<-upAbs[startPoint],
                   startAbs<-round(approx(upTime, upAbs, xout = pcChange, ties = mean)$x,3)
                  )
            
            #StartTime is fitted if abs > threshold, otherwise is closest point
            #This prevents crashing if there are blank wells
            ifelse(max(u)-min(u)<thresh,
                   startTime<-upTime[startPoint],
                   startTime<-round(approx(upAbs, upTime, xout = startAbs, ties = mean)$y,3)
                  )
    
        #Collect the variables
        upcurve<-c(minAbs, startTime,startAbs, startAbs, startPoint, pointmax)
    
        } 
  
      
  
      #Using Purrr to generate the table of results
      TabRes<- reactive ({
    
          ini <- input$ini*.01
          thresh <- input$thresh
          
          Time <- procdat()[[1]]
          #Time <- readData()[[1]]
          TabRes <- procdat()[-1] %>% 
            map_df(~ data.frame(firstAbs=MaxandMin(.x, Time, ini, thresh)[1], 
                                min.abs=MaxandMin(.x, Time, ini, thresh)[2], 
                                max.abs=MaxandMin(.x, Time, ini, thresh)[3], 
                                delta.abs=MaxandMin(.x, Time, ini, thresh)[4], 
                                max.time=MaxandMin(.x, Time, ini, thresh)[5], 
                                pointmax=MaxandMin(.x, Time, ini, thresh)[6], 
                                clot.time=uppity(.x, Time, ini, thresh )[2], 
                                clot.abs=uppity(.x, Time, ini, thresh)[4], 
                                startPoint=uppity(.x, Time, ini, thresh)[5])) %>% 
                                #lys.abs=downy(.x, Time, ini, thresh)[1],
                                #lys.time=downy(.x, Time, ini, thresh)[2],
                                #decayPoint=downy(.x, Time, ini, thresh)[3],
                                #endPoint=downy(.x, Time, ini, thresh)[4],
                                #end.time=downy(.x, Time, ini, thresh)[5])) %>% 
            #AUC=downy(.x, Time, ini, thresh)[6])) %>% 
            #mutate(clotTolys.time=lys.time-clot.time) %>% 
            add_column(Wells=colnames(readData()[-1]), .before = TRUE) %>% 
            add_column(maxrate.time=derivRes()[[1]]) %>% 
            select(Wells, min.abs, clot.time, clot.abs, max.abs,  delta.abs, max.time, maxrate.time, #max.abs,cut until plots functions sorted
                   # lys.time, lys.abs,  clotTolys.time, 
                  startPoint, pointmax) %>%  
                  #decayPoint, endPoint, end.time) %>% 
            mutate(across(where(is.numeric), \(x) round(x, digits=4))) 
          #No clipboard for online app
          #clipr::write_clip(TabRes)
          TabRes
    
      })
 
      #for multiple and single plots
      plot <- reactive({multi_plotFun(procdat(), input$numrows, TabRes())})
                   output$plot <- renderPlot({plot()})
 
      myplot <- reactive({one_plotFun(procdat(), input$colmnames, TabRes())})
                    output$myplot <- renderPlot({myplot()})
                    
                    
      derivRes <- reactive({
        #Make the derivative curves
        Tdif<-readData()[,1][-1]
        Time <- readData()[,1]
        
        fun_diff<-function(d1){
          res<-diff(d1)/diff(Time) #need to deal with the time column here in denominator
          res
        }
        
        #In case of smoothing
        #all smoothing
        
        fun_allsmooth <- function(as){
          smoothfit <- supsmu(Tdif, as)
          smoothfit$y
        }
        
        derivDat<-readData()[-1] %>%  map_df(~fun_diff(.x)) %>% 
          map_df(~fun_allsmooth(.x)) %>%
          add_column(Tdif, .before = 1) %>%
          as.data.frame()
        #derivDats <- derivDat[-1] %>% map_df(~fun_allsmooth(.x)) %>% add_column(Tdif, .before = 1) %>% as.data.frame()
        
        derivDat
        
        ini <- input$ini*.01
        thresh <- input$thresh
        
        Time <- readData()[[1]]
        #Time <- readData()[[1]]
        derivsT <- derivDat[-1] %>% 
          map_df(~ data.frame(max.dtime=MaxandMin(.x, Tdif, ini, thresh)[5] ))
        
        derivsT
        
      })
 
      
      #table of results
      output$head <- renderTable({TabRes()})
  
      #raw or processed data
      output$raw<-renderTable({procdat()})
 
      #table of results to match the multiple plots
      plotTabl<-reactive({
      if(is.null(input$colmnames)){return(NULL)} # To stop this section running and producing an error before the data has uploaded
      
        RowNum<-input$numrows
        
         data<-switch(input$tabRes, 
                    "1"=matrix(TabRes()[,1], byrow=TRUE, nrow=RowNum),
                    "2"=matrix(TabRes()[,2], byrow=TRUE, nrow=RowNum),
                    "3"=matrix(TabRes()[,3], byrow=TRUE, nrow=RowNum),
                    "4"=matrix(TabRes()[,4], byrow=TRUE, nrow=RowNum),
                    "5"=matrix(TabRes()[,5], byrow=TRUE, nrow=RowNum), 
                    "6"=matrix(TabRes()[,6], byrow=TRUE, nrow=RowNum),
                    "7"=matrix(TabRes()[,7], byrow=TRUE, nrow=RowNum), 
                    "8"=matrix(TabRes()[,8], byrow=TRUE, nrow=RowNum),
                   # "9"=matrix(TabRes()[,9], byrow=TRUE, nrow=RowNum),
                   # "10"=matrix(TabRes()[,10], byrow=TRUE, nrow=RowNum),
                   # "11"=matrix(TabRes()[,11], byrow=TRUE, nrow=RowNum),
                   # "12"=matrix(TabRes()[,12], byrow=TRUE, nrow=RowNum),
                   # "13"=matrix(TabRes()[,13], byrow=TRUE, nrow=RowNum),
                   # "14"=matrix(TabRes()[,14], byrow=TRUE, nrow=RowNum),
                   # "15"=matrix(TabRes()[,15], byrow=TRUE, nrow=RowNum),
                    #"16"=matrix(TabRes()[,16], byrow=TRUE, nrow=RowNum),
                   # "17"=matrix(TabRes()[,17], byrow=TRUE, nrow=RowNum),
                   # "18"=matrix(TabRes()[,18], byrow=TRUE, nrow=RowNum)
                
                    )
   
   colnames(data) =as.character(1:(length(data)/RowNum))
   data
   
             })
 
    output$resultsTable<-renderTable({plotTabl()})
    #output$resultsTable<-renderTable({derivRes()})
    
    
 
    #results to go with single plot
    output$curveTable<-renderTable({
        if(is.null(input$colmnames)){return(NULL)} # To stop this section running and producing an error before the data has uploaded
        TabNames<-colnames(TabRes())
        ColNam<-c("Parameter", "Result")
   
         All.Res<-TabRes() %>% 
           mutate(across(where(is.numeric), \(x) round(x, digits=3))) %>% 
           filter(Wells == input$colmnames) %>% select(1:10)
         All.Res.Tab<-cbind(TabNames[1:10], t(All.Res))
         colnames(All.Res.Tab)<-ColNam
         
         All.Res.Tab
         
              })
      #Get the names from results table
      varnames<- reactive({
              mynames<-colnames(TabRes()[1:10])
              })
      
      #Find data for graphs in Explore tab
       output$whatx<-renderUI({
             if(is.null(input$colmnames)){return(NULL)} # To stop this section running and producing an error before the data has uploaded
              selectInput("mynamesx",
                     label= h5("On X axis"),
                     choices = varnames(), selected = colnames(TabRes()[1]))
              })
       #Find data for graphs in Explore tab
       output$whaty<-renderUI({
              if(is.null(input$colmnames)){return(NULL)} # To stop this section running and producing an error before the data has uploaded
              selectInput("mynamesy",
                     label= h5("On Y axis"),
                     choices = varnames(), selected = colnames(TabRes()[3]))
              })
      
       #select data for heatmap
      output$whatheat<-renderUI({
             if(is.null(input$colmnames)){return(NULL)} # To stop this section running and producing an error before the data has uploaded
             selectInput("mynamesheat",
                      label= h5("Choose data for heatmap"),
                      choices = varnames(), selected = colnames(TabRes()[3]))
             })
 
 
        #generate the heatmap or scatterplot in the Explore tab
       output$exploreplot<-renderPlotly({
             if(is.null(input$colmnames)){return(NULL)} # To stop this section running and producing an error before the data has uploaded
            TabRes_r <- TabRes() %>% mutate(across(where(is.numeric), \(x) round(x, digits=2)))
            N<-TabRes()[[1]]
            # X<-TabRes()[,input$mynamesx]
            # Y<-TabRes()[,input$mynamesy]
             X<-TabRes_r[,input$mynamesx]
             Y<-TabRes_r[,input$mynamesy]
             Z<-matrix(TabRes()[,input$mynamesheat], byrow=TRUE, nrow=input$numrows)
             switch(input$heat,
                    "heat" = plot_ly(z= Z, colors = colorRamp(c("grey", "pink", "red")),  type = "heatmap") %>% layout(yaxis = list(autorange = "reversed")), 
                    "scatter" = plot_ly(type = "scatter", mode = "markers") %>% 
                      add_trace(x=X, y=Y, text = N, marker = list(size = 10, color = "pink", line = list(color="red", width = 2)), hoverinfo = N) )
             
             })
 
      #Collection of setting for records in the settings tab
        setsTab<-reactive({
          
          #make a matrix 4 columns, 10 rows to show as a table
           setTable<-matrix(c( 
             "Date", format(Sys.Date(), "%d %b %Y"), "Time", format(Sys.time(), "%X"),
             
             "Read interval s", readData()[2,1]-readData()[1,1],  "Chosen % clotting ", input$ini,
             
             "Threshold value", input$thresh, "", "",
             
             "Raw or spline fitted",  input$spln,  "",  "",
             
             "Fit start", input$zero,  "Truncate for fit", input$trunc, 
             
             "Read interval with fit", round(procdat()[2,1]-procdat()[1,1], 4), "Number of points in fit", input$npoints,
             
             "Zero method", input$abini,  "", "",
             
             "", "","Global zero", input$back, 
             
             "", "","nth point", input$arow,
             
             "", "", "offset with min abs", input$off,
             
             "Curve max method", input$abfin, "", "",
             
             "", "", "global max value", input$finabs
             
                        ),
         byrow=TRUE, nrow=12)
   
        colnames(setTable)<-c("Parameter", "Value", "Parameter", "Value")
        setTable
             })
 
 
      #generate the table 
      output$settings<-renderTable({
            if(is.null(input$colmnames)){return(NULL)}
            setTable<-setsTab()
      setTable
             })
 
    }
