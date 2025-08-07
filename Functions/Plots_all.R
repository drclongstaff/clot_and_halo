
#This function plots data from all the wells
#The ROWNUM argument organises the grid and TABRES is results table for lines and colors
multi_plotFun <- function(PLATE, ROWNUM, TABRES) {
  #define some parameters
        TabRes <- TABRES #results table
        RowNum<-ROWNUM #number of rows to plot
        Time<-PLATE[[1]] #time in the first column
        plateData<-PLATE[,-1] #the absorbance data without the time column
        absWells <- length(plateData[1,]) #the no. of columns of the absorbance data
        mint<-min(Time, na.rm = TRUE)
        maxt<-max(Time, na.rm = TRUE)
        maxy<-max(plateData, na.rm = TRUE) #min and max values for scaling the plots
        samples <- colnames(plateData) #names of columns to go on the plots
        par(mfrow=c(RowNum,(absWells/RowNum)))
        par(mar=c(0.2,0.2,0.2,0.2)) # dimensions for figure
        
        #a loop to plot all the wells
        for(k in seq_along(plateData) ) {
        yi<- plateData[[k]]
        
        plot<-plot(Time, yi, type = "l", col= "steelblue", lwd = 2, xlim= c(0, maxt), 
                   ylim=c(0, maxy*1.2), xaxt="n", yaxt="n", xlab = "" , ylab="")
        
        #points(Time, yi, pch=21, col = "slategrey") # better to omit points for many plots
        #lines to show clotting, lysis etc
        lines(Time[1:TabRes[[k,9]]], yi[1:TabRes[[k,9]]],col="red", lwd=3)
        lines(Time[TabRes[[k,9]]:TabRes[[k,10]]], yi[TabRes[[k,9]]:TabRes[[k,10]]],col="green", lwd=3)
       # lines(Time[TabRes[[k,12]]:TabRes[[k,13]]], yi[TabRes[[k,12]]:TabRes[[k,13]]],col="blue", lwd=3)
       # lines(Time[TabRes[[k,13]]:TabRes[[k,14]]], yi[TabRes[[k,13]]:TabRes[[k,14]]],col="darkorange2", lwd=3)
        
        #placing the legend is fiddly with different data
        #legend(TabRes[k,3],maxy*1.2, bty="n", paste0(samples[k],"=",k), cex=1.5)
        #legend(TabRes[k,3], maxy*1.2, xjust=TabRes[k,3]/maxt, bty="n", paste0(samples[k],"=",k), cex=1.5)
        #legend(TabRes[k,3],TabRes[k,5]*1.2, bty="n", paste0(samples[k],"=",k), cex=1.5)
        legend("topright", bty="n", paste0(samples[k],"=",k), cex=1.5 )
        
        #results lines for baseline and times
        abline("v"= TabRes[k,3], lty=2)
        abline("v"= TabRes[k,7], lty=2)
        abline("v"= Time[1], lty=2)
        abline("h"= TabRes[k,2], lty = 2)
        abline("v"= TabRes[k,8], col="magenta")
        abline("h"= TabRes[k,4], lty = 2)
        abline("h"= TabRes[k,10], lty = 2)
        abline("h"= TabRes[k,5], lty = 2)
      }
  
  }


