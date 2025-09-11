#' @title Add carryover dummy variables
#'
#' @description Create dummy variables associated with first-order
#' carryover effect in a Crossover Design
#'
#' @param data A data frame with the variables of
#' the crossover experimental design
#' @param period A character string specifying the name of vector with the
#' observation period of the responses of  the crossover experimental design
#' @param treatment A character string specifying the name of vector with the
#' treatment applied at each observation of the crossover experimental design
#' @param id A  character string specifying the name of vector which identifies
#' the experimental units.
#' @param carrySimple \code{TRUE} = simple carry-over, where the residual effect of a
#' treatment affects equally each of the treatments that are preceded by it,
#'  and \code{FALSE} = carry-over complex, where the residual effect of the
#'   treatment affects each of the other treatments differently.
#' @return A list with:
#' \item{data}{A data frame with all the variables of the crossover
#' experimental design and the carryover variables}
#' \item{carryover}{The names of new carryover variables.}
#'
#' @examples
#' data(Water)
#' carryover <- createCarry(data=Water,
#'                          treatment = "Treatment", id = "ID",
#'                          period = "Period", carrySimple = FALSE)
#' carryover$carryover
#' carryover$data
#' @export
#' @importFrom dplyr "%>%"


createCarry <- function(data, treatment, period, id, carrySimple=TRUE){
  data_temp <- data %>% dplyr::arrange_at(c(id, period)) %>% data.frame()
  data_temp[,period] <- factor(data_temp[,period])
  data_temp[,treatment] <- factor(data_temp[,treatment])
  treat <- levels(data_temp[,treatment])
  per <- levels(data_temp[,period])
  addCarry <- character()
  if(carrySimple){
    for ( trett in treat[-1]){
      vari <- paste("Carry", trett,sep="_")
      addCarry <- c(addCarry, vari)
      data_temp[, vari] <- 0
      for (ids in unique(data_temp[, id])){
        for ( peri in 2:length(per)){
          treatAnt <- data_temp[data_temp[,id]==ids &
                                  data_temp[,period]==per[peri-1],treatment][1]
          data_temp[data_temp[,id]==ids &
           data_temp[,period]==per[peri],vari] <-ifelse(trett==treatAnt,1,0)
        }
      }
    }
  }else{
    for ( trettA in treat){
      for ( trettB in treat[treat!=trettA]){
        vari = paste("Carry", trettA, "over", trettB,sep="_")
        addCarry =c(addCarry, vari)
        data_temp[, vari] <- 0
        for (ids in unique(data_temp[, id])){
          for ( peri in 2:length(per)){
            treatAnt <- data_temp[data_temp[,id]==ids &
                                  data_temp[,period]==per[peri-1],treatment][1]
            data_temp[data_temp[,id]==ids & data_temp[,period]==per[peri] &
                        data_temp[, treatment]==trettB,vari] <-
              ifelse(trettA==treatAnt,1,0)
      }
    }

  }
    }
  }
  print(paste("This variable was added to data:", addCarry))
return(list(data=data_temp, carryover = addCarry))
}


