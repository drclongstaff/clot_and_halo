
#This function plots data from the selected well 'WELLNUM'
#The table of results, TabRes, is used to draw the lines and pick out curve features
one_plotFun <- function(PLATE, WELLNUM, TABRES) {
  #define some parameters for the plot      
        TabRes <- TABRES
        #k<-WELLNUM #the selected well
        k<-which(colnames(PLATE)==WELLNUM)[1]-1
        #k <- 4
        RowNum <- 1 #here there is only one plot, so only one row is needed
        Time<-PLATE[[1]] #time in the first column of the plate data
        plateData<-PLATE[,-1] #The platedata without the time column
        absWells <- 1 #This is for single plot
        mint<-min(Time, na.rm = TRUE)
        maxt<-max(Time, na.rm = TRUE)
        maxy<-max(plateData, na.rm = TRUE) #These min and max values are for scaling the plots
        samples <- colnames(plateData) #for writing the well name
        par(mfrow=c(RowNum,(absWells/RowNum)))
        par(mar=c(4,4,1,1)) # dimensions for figure
  
        yi<- plateData[[k]]
        #yi <- plateData[[which(colnames(plateData)==WELLNUM)[1]]]
        plot<-plot(Time, yi, type = "l", col= "steelblue", lwd = 2, xlim= c(0, maxt), 
                   ylim=c(0, maxy*1.2), ylab="Absorbance")
        
        points(Time, yi, pch=21, col = "slategrey")
        #different coloured lines for start to clotting to max to %lysis and to full lysis
        lines(Time[1:TabRes[[k,9]]], yi[1:TabRes[[k,9]]],col="red", lwd=3)
        lines(Time[TabRes[[k,9]]:TabRes[[k,10]]], yi[TabRes[[k,9]]:TabRes[[k,10]]],col="green", lwd=3)
       # lines(Time[TabRes[[k,12]]:TabRes[[k,13]]], yi[TabRes[[k,12]]:TabRes[[k,13]]],col="blue", lwd=3)
       # lines(Time[TabRes[[k,13]]:TabRes[[k,14]]], yi[TabRes[[k,13]]:TabRes[[k,14]]],col="darkorange2", lwd=3)
        
        #legend including the well number and name
        legend(TabRes[k,3], maxy*1.2, xjust=TabRes[k,3]/maxt, bty="n", paste0(samples[k],"=",k), cex=2)
    
        #dashed lines showing the baseline and times
        abline("v"= TabRes[k,3], lty=2)
        abline("v"= TabRes[k,7], lty=2 )
        #abline("v"= TabRes[k,8], lty=1, col="magenta")
        abline("v"= Time[1], lty=2)
        abline("h"= TabRes[k,2], lty = 2)
        abline("v"= TabRes[k,8], col="magenta")
        abline("h"= TabRes[k,4], lty = 2)
        abline("h"= TabRes[k,9], lty = 2)
        abline("h"= TabRes[k,5], lty = 2)

    }


  
