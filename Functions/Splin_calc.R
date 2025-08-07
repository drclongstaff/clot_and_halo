##Spline fitting

library(dplyr)
library(purrr)
library(tibble)

#myPlate <- read.csv("Data/ClotLysisDoses.csv")
#myPlate <- read.csv("Data/Clotlysistrim.csv")


allSpline <- function(S,  N, MIN, MAX, PLATE){
Ts <- PLATE[[1]]
aSpline <- function(S, TS, N, MIN, MAX){
  maxx<-max(Ts)-MAX
  DatSP <- spline(TS,S, N, MIN, maxx,   method = "natural") 
        }
DatSP <- map_df(PLATE[,-1], ~aSpline(.x, Ts, N, MIN, MAX))
DatSPdf <- data.frame(split(DatSP$y,   ceiling(seq_along(DatSP$y) / N)))
colnames(DatSPdf) <- colnames(PLATE[,-1])
DatSPdf <- cbind("Ts"=DatSP$x[1:N], DatSPdf)
 
 }


#myPlateS <- allSpline(NULL, 50, 0, 5000, myPlate)

#multi_plotsimple(myPlateS, 3)

#plot(myPlateS$Ts, myPlateS$F5, type = "l")
