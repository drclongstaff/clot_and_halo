#This file will hold the functions for all the plotting
# This is the recipe for making a plot
plotmake_fun <- function(myPlate, TabRes, k, axx, axy) {
  Time <- myPlate[[1]] # Time is first column
  plateData <- myPlate[, -1] # All absorbance data without time column
  absWells <- length(plateData[1, ])
  mint <- min(Time, na.rm = TRUE)
  maxt <- max(Time, na.rm = TRUE)
  miny <- min(plateData, na.rm = TRUE)
  maxy <- max(plateData, na.rm = TRUE) # To set the maximum absorbance for plotting
  samples <- colnames(plateData)
  
  y.k <- plateData[[k]]
  
  plot <- plot(Time, y.k,
               type = "l", col = "blue", lwd = 3, xlim = c(0, maxt),
               ylim = c(0, maxy), ylab = "Abs" , xaxt = axx, yaxt = axy
  )
  points(Time, y.k, pch = 21, col = "gray70", cex = .8)
  
  lines(Time[TabRes[k, 9]:TabRes[k, 10]], y.k[TabRes[k, 9]:TabRes[k, 10]], col = "green", lwd = 3)
  lines(Time[1:TabRes[k, 9]], y.k[1:TabRes[k, 9]], col = "red", lwd = 3)
  
  #legend(TabRes[k, 8], maxy * 0.95, xjust = TabRes[k, 8] / maxt, bty = "n", paste0(samples[k], "=", k), cex = 1.5)
  legend("topright", bty = "n", paste0(samples[k], "=", k), cex = 1.5)
  
  abline("v" = TabRes[k, 3], lty = 2)
  abline("v" = TabRes[k, 7], lty = 2)
  abline("h" = TabRes[k, 2], lty = 2)
  abline("v"= TabRes[k,8], col="magenta")
  arrows(mint, TabRes[k, 4], TabRes[k, 3], TabRes[k, 4], length = 0.075, angle = 10)
  arrows(mint, TabRes[k, 5], TabRes[k, 7], TabRes[k, 5], length = 0.075, angle = 10)
  abline("h"= TabRes[k,9], lty = 2)
}

# This function plots data from a chosen well 
one_plotFun <- function(PLATE, WELLNUM, TABRES) {
  # Set plotting parameters
  TabRes <- TABRES
  
  #k <- WELLNUM
  k<-which(colnames(PLATE)==WELLNUM)[1]-1
  
  par(mar = c(4, 4, 1, 1)) # dimensions for figure
  
  # Only one plot in this case so no looping needed
  plotmake_fun(PLATE, TabRes, k, axx="s", axy="s")
}

#plotting all wells
multi_plotFun <- function(PLATE, ROWNUM, TABRES) {
  # Set up some plotting parameters
  TabRes <- TABRES
  RowNum <- ROWNUM
  
  plateData <- PLATE[, -1]
  absWells <- length(plateData[1, ])
  par(mfrow = c(RowNum, (absWells / RowNum))) # Organisation of multiple plots
  par(mar = c(0.2, 0.2, 0.2, 0.2)) # Dimensions for figure
  
  # Generate the plots from plotmake_fun
  lapply(seq_along(plateData), function(k) {
    # for(k in seq_along(plateData)){ #This is an alternative loop
    
    plotmake_fun(PLATE, TabRes, k, axx="n", axy="n")
  }) # remove this ')' with alternative loop
}

