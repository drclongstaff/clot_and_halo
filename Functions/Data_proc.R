#Initial data processing

#Some functions for data processing

#This is minus a fixed background value

BaselineBack <- function(n, Back){
  n-Back
}

#This function subtracts a background from all data defined by row in all data
BaselineNP<-function(n, NP){
  n <- n-n[NP]
}

#Subtract the min abs of each column and include an offset
BaselineOff <- function(n, OF) {
  n <- n-min(n)-OF
    }
 
 



