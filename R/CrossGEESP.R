#' @title Run a semi-parametric GEE model for data from a crossover experiment with
#'  repeated measures
#'
#' @description Provides a GEE model for the data of a crossover design with S
#' sequences of T periods. There must be at least two observations of each
#' experimental unit in each period. The effect of time within period and
#' the possible carryover effects are modeled by means of splines.
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
#' @param family See corresponding documentation to \code{glm}
#' @param correlation 	a character string specifying the correlation structure.
#'   The following are permitted: "independence", "fixed", "stat_M_dep",
#'  "non_stat_M_dep", "exchangeable", "AR-M" and "unstructured"
#' @param Mv When correlation is "stat_M_dep", "non_stat_M_dep", or "AR-M"
#'   then Mv must be specified.
#' @param formula A formula related the response variable with the explanatory
#'  variables. If it is \code{NULL}, formula
#'  \code{response~period+treatment+carry+time+covar} will be evaluated.
#' @param tol 	the tolerance used in the fitting algorithm.
#' @param niter the maximum number of iterations.
#' @param nodes Number of nodes in the estimation of the splines.
#' By default, the base 2 logarithm of the number of observations per period is used.
#' @return \code{QIC} The QIC of the model: The model are fitted by \code{geeglm}
#' @return \code{model} The model fitted by \code{geeglm}.
#' @return \code{graphs} The graphs estimated by splines.
#' In position 1 the graph of the effect of time appears and from then on,
#' it appears one for each carryover effect declared in the \code{carry} option.
#' The graphs are built with \code{ggplot2}, therefore they allow manipulation
#' of axes and other graphic parameters of that library.
#' @references Cruz Gutierrez NA, Melo OO, Martinez CA. Semiparametric generalized
#'  estimating equations for repeated measurements in cross-over designs.
#'   Statistical Methods in Medical Research, 2023;32(5):1033-1050.
#' @source https://doi.org/10.1177/09622802231158736
#' @examples
#' data(Arterial)
#'
#' carrydata <- createCarry(data=Arterial, treatment = "Treatment",
#'                          period = "Period",id="Subject", carrySimple = FALSE)
#' data <- carrydata$data
#' carry <- carrydata$carryover
#' model1 <- CrossGEESP(response = "Pressure", treatment = "Treatment",
#'                     period = "Period", id="Subject", time="Time",
#'                     carry=carrydata$carryover,data=data,
#'                     correlation = "exchangeable")
#'
#'
#' model2 <- CrossGEESP(response = "Pressure", treatment = "Treatment",
#'                      period = "Period", id="Subject", time="Time",
#'                      carry=carrydata$carryover,data=data, correlation = "AR-M")
#'
#'
#' model1$QIC
#' model2$QIC
#' summary(model1$model)
#' summary(model2$model)
#' model1$graph[[1]]
#' model1$graph[[2]]
#' plot <- model1$graph[[1]] + ggplot2::xlab("Time in minutes")+
#' ggplot2::ylab("Change in systolic blood pressure")
#' plot
#' @export
#' @importFrom dplyr "%>%"


CrossGEESP <- function(response,period,treatment,id,
                         time,carry, covar=NULL,data,
                         family=gaussian, correlation="independence",
                         formula =NULL,tol = 1e-4, niter=100, nodes=NULL, Mv=1){
  data["Per_id"]=as.numeric(as.factor(paste(data[,id], data[,period])))
  totalVar <- c(response, period, treatment,carry, covar, id)
  if(sum(totalVar %in% names(data))!=length(totalVar)){
    stop("Some variables are not present in data")
  }
  data[,period] <- factor(data[,period])

  data <- data %>% dplyr::select_at(c(totalVar,id, period, time, "Per_id" ))  %>%
    dplyr::arrange_at(c(id, period, time)) %>%
    data.frame() %>% stats::na.exclude() %>% data.frame()

  data["id"] <- as.numeric(as.factor(data[,id]))
  if(is.null(formula)){
    form1 <- stats::as.formula(paste(response,
                              paste(c(period, treatment, covar),
                                    collapse=" + "), sep=" ~ "))

  }else{
    form1 <- formula}
  if(correlation=="ar1"){
    correlation1 <- "AR-M"
  }else{
    correlation1 <- correlation}

  Time_total <- data[,time]
  if(is.null(nodes)){
    nodes <- max(1, floor(log(length(unique(Time_total)),2)))
    if(length(unique(Time_total))<8){
      nodes=1
    }
    if(length(unique(Time_total))<5){
      stop("Too few repetitions per period to be able to estimate splines")
    }
  }
  splines1 <- splines::bs(Time_total, knots = stats::quantile(Time_total, 1:nodes/(nodes+1)))
  splines2 <- splines1
  for(carryover in carry){
    TimeTemp <- data %>% dplyr::filter_at(carryover,~.x !=0) %>%
      dplyr::select_at(time) %>%  data.frame()
    carryTemp <- splines::bs(TimeTemp[,1], knots = stats::quantile(TimeTemp[,1], 1:nodes/(nodes+1)))
    splinesTemp <- matrix(0, ncol=ncol(splines1), nrow=nrow(splines1))
    j=1
    for(i in 1:nrow(data)){
      if(data[,carryover][i]!=0){
        splinesTemp[i,]<- carryTemp[j,]
        j=j+1
      }
    }
    splines2 <- cbind(splines2, splinesTemp)
  }

  varTempSp <- paste("Var", 1:ncol(splines2), sep="")
  colnames(splines2) <- varTempSp

  diference <- 1
  counter <- 0
  modTemp <- suppressMessages(gee::gee(form1,data=data,
                                  corstr = correlation1, id=Per_id))

  while(diference>tol & counter < niter){
    beta0 <- stats::coef(modTemp)
    Xbeta <- stats::model.matrix(stats::lm(form1, data=data))
    R_alpha <- modTemp$working.correlation
    dataTemp1 <- data.frame(Xbeta=Xbeta %*%beta0, splines2,
                            data %>% dplyr::select_at(c(response,period,
                                                        treatment, covar,"Per_id")) %>%
                              data.frame())
    form2 <- stats::as.formula(paste(response,
                              paste(c(varTempSp, "Xbeta"),
                                    collapse=" + "), sep=" ~ "))
    modTemp <- suppressMessages(gee::gee(form2, data=dataTemp1,
                                    R=R_alpha, maxiter = 10,
                   corstr = "fixed", id=Per_id))
    alpha <- stats::coef(modTemp)[2:(ncol(splines1)+1)]
    ZZ <- dataTemp1[,2:(ncol(splines1)+1)]
    TT <- as.matrix(ZZ)%*%alpha

    for(ii in 1:(length(carry))){
      alphaTemp <- stats::coef(modTemp)[(2+ncol(splines1)*ii):(1+ncol(splines1)*(ii+1))]
      ZTemp <- dataTemp1[,(2+ncol(splines1)*ii):(1+ncol(splines1)*(ii+1))]
      TTemp=as.matrix(ZTemp)%*%alphaTemp
      TT <- cbind(TT, TTemp)
    }
    colnames(TT) <- paste("Eff", c(time, carry), sep="")
    dataTemp2 <- data.frame(data,TT)
    form3 <- stats::as.formula(paste(response,
                              paste(c(period, treatment, covar,
                                      paste("offset(",colnames(TT),")",sep="" )),
                                    collapse=" + "), sep=" ~ "))

    modTemp <- suppressMessages(gee::gee(form3,data=dataTemp2,
                                    corstr =correlation1, id=Per_id))
    summary(modTemp)
    diference <- sum((beta0-stats::coef(modTemp))^2)
    counter <- counter+1
  }

  if(counter >= niter){
    stop(("Convergence not achieved, change the formula"))
  }

  modTemp2=suppressMessages(gee::gee(form2, data=dataTemp1, R=R_alpha, maxiter = 10,
               corstr = "fixed", id=Per_id))
  graphs <- list()
  alpha <- stats::coef(modTemp2)[2:(ncol(splines1)+1)]
  ZZ <- as.matrix(dataTemp1[,2:(ncol(splines1)+1)])
  varTemp <- modTemp2$robust.variance[2:(ncol(splines1)+1), 2:(ncol(splines1)+1)]
  TT <- as.matrix(ZZ)%*%alpha
  varTT <- ZZ %*% varTemp %*% t(ZZ)
  dataTemp5 <- data.frame(Time = Time_total, Y_hat = TT, VarY_hat =diag(varTT)) %>%
    base::unique() %>% dplyr::mutate(liminf=Y_hat - 1.96*sqrt(VarY_hat),
                                     limsup=Y_hat + 1.96*sqrt(VarY_hat)) %>%
    dplyr::select(Time, Y_hat, liminf, limsup)
  colnames(dataTemp5)=c("Time", "Estimate_time_effect", "lower_band","upper_band")
  plot1 <- dataTemp5 %>% ggplot2::ggplot( ggplot2::aes(x = Time, y =Estimate_time_effect)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_band, ymax = upper_band), alpha=0.4) +
    ggplot2::geom_line(color="black",lwd=.8, linetype='solid') +
    ggplot2::xlab("Time")+ggplot2::ylab(paste("Change in", response))+
    ggplot2::geom_hline(yintercept = 0)+ggplot2::geom_vline(xintercept = 0)

  graphs[[1]] <- plot1

  for(ii in 1:(length(carry))){
    alpha <- stats::coef(modTemp2)[(2+ncol(splines1)*ii):(1+ncol(splines1)*(ii+1))]
    ZZ <- as.matrix(dataTemp1[data[, carry[ii]]==1,
                    (2+ncol(splines1)*ii):(1+ncol(splines1)*(ii+1))])
    varTemp <- modTemp2$robust.variance[(2+ncol(splines1)*ii):(1+ncol(splines1)*(ii+1)),
                                        (2+ncol(splines1)*ii):(1+ncol(splines1)*(ii+1))]
    TT <- ZZ%*%alpha
    varTT <- ZZ %*% varTemp %*% t(ZZ)
    dataTemp5 <- data.frame(Time = Time_total[data[, carry[ii]]==1],
                            Y_hat = TT, VarY_hat =diag(varTT)) %>%
      base::unique() %>% dplyr::mutate(liminf=Y_hat - 1.96*sqrt(VarY_hat),
                                       limsup=Y_hat + 1.96*sqrt(VarY_hat)) %>%
      dplyr::select(Time, Y_hat, liminf, limsup)
    colnames(dataTemp5)=c("Time", "Estimate_time_effect", "lower_band","upper_band")
    plot1 <- dataTemp5 %>% ggplot2::ggplot(ggplot2::aes(x = Time, y =Estimate_time_effect)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_band, ymax = upper_band), alpha=0.4) +
      ggplot2::geom_line(color="black",lwd=.8, linetype='solid') +
      ggplot2::xlab("Time")+ggplot2::ylab(carry[ii])+
      ggplot2::geom_hline(yintercept = 0)+ggplot2::geom_vline(xintercept = 0)

    graphs[[ii+1]] <- plot1
  }

  names(graphs) <- c(time, carry)

  model1 <- suppressMessages(gee::gee(form3,data=dataTemp2,
                                             corstr =correlation, id=Per_id))
  QICmodels <- data.frame(computeqic(model1))
  names(QICmodels) <- "QICs"
  return(list(QIC =QICmodels, model=model1, graphs=graphs))
}
