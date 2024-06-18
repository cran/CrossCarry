#' @title Run a GEE model for data from a crossover experiment with
#'  repeated measures
#'
#' @description Provides a GEE model for the data of a crossover design with S
#' sequences of T periods. There must be at least two observations of each
#' experimental unit in each period.
#'
#'
#' @param response A character string specifying the name of the
#' response variable of the crossover experimental design
#' @param period A character string specifying the name of vector with the
#' observation period of the responses of  the crossover experimental design
#' @param treatment A character string specifying the name of vector with the
#' treatment applied at each observation of the crossover experimental design
#' @param id A  character string specifying the name of vector which identifies
#' the experimental units. The length of ‘id’
#' should be the same as the number of observations. Data are assumed to be sorted so
#' that observations on each cluster appear as contiguous rows in data. If data
#' is not sorted this way, the function will not identify the clusters correctly.
#'  If data is not sorted this way, a warning will be issued.
#' @param time A character string specifying the name of the vector with the
#'  measurement time within each period
#' @param carry A vector of  character string specifying the name set of dummy
#'  variables that indicates the treatment applied
#' in the previous period of each experimental unit. They must be 0 in period 1
#' @param covar A vector of  character string specifying the name of possible
#' covariates of the crossover experimental design
#' @param data A data frame with all the variables of the crossover experimental design
#' @param family See corresponding documentation to \code{glm}.
#' @param correlation character string specifying the correlation within periods
#'  structure.
#'  The following are permitted: "independence", "exchangeable", "ar1" and
#'  "unstructured".
#' @param formula A formula related the response variable with the explanatory
#'  variables. If it is \code{NULL} the formula,
#' \code{response~period+treatment+carry+time+covar} will be evaluated
#' @param tol 	the tolerance used in the fitting algorithm.
#' @param niter the maximum number of iterations.
#'  \code{response~period+treatment+carry+time+covar} will be evaluated
#' @return \code{QIC} The QIC of the model: The model are fitted by \code{geeglm}
#' @return \code{model} The model fitted by \code{geeglm}.
#' @return \code{Within} The estimated correlation matrix within the period
#' with the structure determined by \code{correlation}.
#' @return \code{Between} The estimated correlation matrix between periods
#' @references Cruz, N.A., Melo, O.O. & Martinez, C.A.
#'  A correlation structure for the analysis of Gaussian and non-Gaussian
#'   responses in crossover experimental designs with repeated measures.
#'    Statistical Papers 65, 263–290 (2024)
#' @source https://doi.org/10.1007/s00362-022-01391-z
#' @examples
#' data(Arterial)
#'
#' carrydata <- createCarry(data=Arterial, treatment = "Treatment",
#'  period = "Period",id="Subject")
#'
#' data <- carrydata$data
#' carry <- carrydata$carryover
#' model <- CrossGEEKron(response = "Pressure", treatment = "Treatment",
#' period = "Period", id="Subject", time="Time",
#'  carry=c("Carry_B","Carry_C"),data=data, correlation = "ar1")
#'
#' model$QIC
#' model$Within
#' model$Between
#' summary(model$model)
#' model2 <- CrossGEEKron(response = "Pressure", treatment = "Treatment",
#'  period = "Period", id="Subject", time="Time",
#'  carry=c("Carry_B","Carry_C"), data=data,
#'  correlation = "ar1",formula=Pressure ~ Treatment+
#'  Period+ Carry_B+Carry_C)
#'
#' model2$QIC
#' model2$Within
#' model2$Between
#' summary(model2$model)

#' @export
#' @importFrom dplyr "%>%"
#' @importFrom stats "gaussian"


CrossGEEKron <- function(response,period,treatment,id,
                         time,carry, covar=NULL,data,
                         family=gaussian(), correlation="independence",
                         formula =NULL,tol = 1e-4, niter=100){
  data["Per_id"]=as.numeric(as.factor(paste(data[,id], data[,period])))
  totalVar <- c(response, period, treatment,carry, covar, id)
  if(sum(totalVar %in% names(data))!=length(totalVar)){
    stop("Some variables are not present in data")
  }

  data <- data %>% dplyr::select_at(c(totalVar,id, period, time, "Per_id" ))  %>%
    dplyr::arrange_at(c(id, period, time)) %>%
    data.frame() %>% stats::na.exclude() %>% data.frame()
  data["id"] <- data[id]
  data <- data %>% dplyr::mutate(id = as.numeric(factor(id)))
  data[id] <- data["id"]
  data["id"] <- as.numeric(as.factor(data[,id]))
  if(is.null(formula)){
    form1 <- stats::as.formula(paste(response,
                              paste(c(period, treatment,carry, time, covar),
                                    collapse=" + "), sep=" ~ "))

  }else{
    form1 <- formula}
  if(correlation=="ar1"){
    correlation <- "AR-M"
  }

  init=gee::gee(form1,family=family,
                data=data, id=id, corstr = correlation)
  pp <- length(init$coefficients)
  resid <- init$residuals
  phi <- sum(resid^2)/(length(resid)-pp)
  PP <- length(table(data[period]))
  LL <- max(table(data["Per_id"]))
  n_suj <- length(table(data[id]))
  if(correlation=="AR-M"){
    K1 <- sum(table(data["Per_id"])-1)
    n <- length(table(data["Per_id"]))
    suma <- 0
    for(i in 1:n){
      residu_temp <- resid[data["Per_id"]==i]
      n_temp <- length(residu_temp)
      suma <- suma + sum(residu_temp[-n_temp]*residu_temp[-1])
    }
    alpha <- suma/((K1-pp)*phi)
    Ralpha <- matrix(0, nrow=LL, ncol=LL)
    for(i in 1:LL){
      for(j in 1:LL){
        Ralpha[i,j] <- alpha^(abs(i-j))
      }
    }
  }else if(correlation=="exchangeable"){
    K1 <- sum((table(data["Per_id"])-1)*(table(data["Per_id"])))
    n <- length(table(data["Per_id"]))
    suma <- 0
    for(i in 1:n){
      residu_temp <- resid[data["Per_id"]==i]
      a_temp <- outer(residu_temp, residu_temp)
      suma <- suma + sum(a_temp[lower.tri(a_temp)])
    }
    alpha <- suma/((K1-pp)*phi)
    Ralpha <- stats::toeplitz(c(1, rep(alpha, LL-1)))
  }else if(correlation=="independence"){
    alpha <- 0
    Ralpha <- stats::toeplitz(c(1, rep(alpha, LL-1)))
  }else if(correlation=="unstructured"){
    alpha <- 0.5
    Ralpha <- stats::toeplitz(c(1, rep(alpha, LL-1)))
  }


  data["residu"] <- resid/sqrt(phi)
  Psi <- matrix(0, PP, PP)
  r_media <- matrix(0, PP, LL)
  for(ii in 1:PP){
    media_temp <- rep(0, LL)
    data_temp <- data[data[period]==ii,]
    for(jj in 1:n_suj){
      media_temp <- media_temp + data_temp[data_temp[id]==jj, "residu"]
    }
    r_media[ii,] <- media_temp/(n_suj)
  }
  for(ii in 1:PP){
    for(jj in 1:PP){
      suma_temp <- 0
      for( kk in 1:n_suj){
        data_temp1 <- data[data[period]==ii,]
        data_temp2 <- data[data[period]==jj,]
        residu1 <- data_temp1[data_temp1[id]==kk, "residu"]
        residu2 <- data_temp2[data_temp2[id]==kk, "residu"]
        AAA <- solve(Ralpha) %*%((residu1- r_media[ii,])%*%
                                  t((residu2- r_media[jj,])))
        suma_temp  <- suma_temp  + sum(diag(AAA))
      }
      Psi[ii,jj] <- suma_temp /n_suj
    }
  }
  a=Psi
  Psi=(diag(1/sqrt(diag(a))))%*% a %*% (diag(1/sqrt(diag(a))))

  RR <- kronecker(Psi, Ralpha)
  beta_init <- init$coefficients
  diference <- 1
  counter <- 0
  while(diference>tol & counter < niter){
    model <- suppressWarnings(gee::gee(form1,family=family,
                 data=data, id=id, corstr = "fixed",
                 R = RR, b = beta_init, maxiter = 1, tol=100))
    diference <- sum((beta_init-model$coefficients)^2)
    beta_init <- model$coefficients
    resid <- model$residuals
    phi <- sum(resid^2)/(length(resid)-pp)
    PP <- length(table(data[period]))
    LL <- max(table(data["Per_id"]))
    n_suj <- length(table(data[id]))
    if(correlation=="AR-M"){
      K1 <- sum(table(data["Per_id"])-1)
      n <- length(table(data["Per_id"]))
      suma <- 0
      for(i in 1:n){
        residu_temp <- resid[data["Per_id"]==i]
        n_temp <- length(residu_temp)
        suma <- suma + sum(residu_temp[-n_temp]*residu_temp[-1])
      }
      alpha <- suma/((K1-pp)*phi)
      Ralpha <- matrix(0, nrow=LL, ncol=LL)
      for(i in 1:LL){
        for(j in 1:LL){
          Ralpha[i,j] <- alpha^(abs(i-j))
        }
      }
    }else if(correlation=="exchangeable"){
      K1 <- sum((table(data["Per_id"])-1)*(table(data["Per_id"])))
      n <- length(table(data["Per_id"]))
      suma <- 0
      for(i in 1:n){
        residu_temp <- resid[data["Per_id"]==i]
        a_temp <- outer(residu_temp, residu_temp)
        suma <- suma + sum(a_temp[lower.tri(a_temp)])
      }
      alpha <- suma/((K1-pp)*phi)
      Ralpha <- stats::toeplitz(c(1, rep(alpha, LL-1)))
    }else if(correlation=="independence"){
      alpha <- 0
      Ralpha <- stats::toeplitz(c(1, rep(alpha, LL-1)))
    }else if(correlation=="unstructured"){
      alpha <- 0.5
      Ralpha <- stats::toeplitz(c(1, rep(alpha, LL-1)))
    }
    data["residu"] <- resid/sqrt(phi)
    Psi <- matrix(0, PP, PP)
    r_media <- matrix(0, PP, LL)
    for(ii in 1:PP){
      media_temp <- rep(0, LL)
      data_temp <- data[data[period]==ii,]
      for(jj in 1:n_suj){
        media_temp <- media_temp + data_temp[data_temp[id]==jj, "residu"]
      }
      r_media[ii,] <- media_temp/(n_suj)
    }
    for(ii in 1:PP){
      for(jj in 1:PP){
        suma_temp <- 0
        for( kk in 1:n_suj){
          data_temp1 <- data[data[period]==ii,]
          data_temp2 <- data[data[period]==jj,]
          residu1 <- data_temp1[data_temp1[id]==kk, "residu"]
          residu2 <- data_temp2[data_temp2[id]==kk, "residu"]
          AAA <- solve(Ralpha) %*%((residu1- r_media[ii,])%*%
                                     t((residu2- r_media[jj,])))
          suma_temp  <- suma_temp  + sum(diag(AAA))
        }
        Psi[ii,jj] <- suma_temp /n_suj
      }
    }
    a <- Psi
    Psi <- (diag(1/sqrt(diag(a))))%*% a %*% (diag(1/sqrt(diag(a))))
    RR <- kronecker(Psi, Ralpha)
    counter <- counter+1
  }

  if(counter >= niter){
    stop(("Convergence not achieved, change the formula"))
  }
  model1 <- geepack::geeglm(form1, data=data, family=gaussian,
                            corstr = "userdefined",
                            zcor=rep(RR[upper.tri(RR)],n_suj),
                            id=id)
  QICmodels <- data.frame(geepack::QIC(model1))
  return(list(QIC =QICmodels, model=model1, Within=Ralpha, Between=Psi))

}
