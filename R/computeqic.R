#' Quasi Information Criterion
#'
#' Function for calculating the quasi-likelihood under the independence model
#' information criterion (QIC), quasi-likelihood, correlation information
#' criterion (CIC), and corrected QIC for one or several fitted geeglm model
#' object from the geepack package.
#'
#' QIC is used to select a correlation structure. The QICu is used to compare
#' models that have the same working correlation matrix and the same
#' quasi-likelihood form but different mean specifications. CIC has been
#' suggested as a more robust alternative to QIC when the model for the mean
#' may not fit the data very well and when models with different correlation
#' structures are compared.
#'
#' Models with smaller values of QIC, CIC, QICu, or QICC are preferred.
#' @param object a fitted GEE model from the gee package.
#' @return A vector or matrix with the QIC, QICu, quasi likelihood,
#'     CIC, the number of mean effect parameters, and the corrected
#'     QIC for each GEE object
#' @author Alirio Cruz \email{nelson-alirio.cruz@uib.es},
#'         Claus Ekstrom \email{claus@rprimer.dk},
#'         Brian McLoone \email{bmcloone@pdx.edu} and
#'          Steven Orzack \email{orzack@freshpond.org}
#' @references Pan, W. (2001). Akaike's information criterion in
#'     generalized estimating equations. Biometrics, 57, 120-125.
#'     Hardin, J.W.  and Hilbe, J.M. (2012). Generalized
#'     Estimating Equations, 2nd Edition, Chapman and Hall/CRC: New
#'     York.
#'
#'     Hin, L.-Y. and Wang, Y-G.
#'     (2009). Working-correlation-structure identification in
#'     generalized estimating equations, Statistics in Medicine 28:
#'     generalized estimating equations, Statistics in Medicine 28:
#'     642-658.  Thall, P.F.  and Vail, S.C. (1990).
#'     Covariance Models for Longitudinal Count Data with
#'     Overdispersion.  Biometrics, 46, 657-671.
#' @examples
#'
#' library(gee)
#' data(Arterial)
#' fit <- gee(Pressure ~ Time + Treatment, id=Subject,
#'        data=Arterial, family=gaussian, corstr="AR-M")
#' computeqic(fit)
#'
#' @export
#' @importFrom MASS ginv
#' @importFrom stats coef
#' @importFrom stats family

computeqic <- function(object) {
  ## Fitted and observed values for quasi likelihood
  mu <- object$fitted.values
  y  <- object$y

  ## Quasi Likelihood for Poisson
  ## quasi.R <- sum((y*log(mu.R)) - mu.R) # poisson()$dev.resids - scale and weights = 1
  type <- family(object)$family
  quasi <- switch(type,
                  poisson  = sum((y * log(mu)) - mu),
                  gaussian = sum(((y - mu)^2)/-2),
                  binomial = sum(y * log(mu / (1 - mu)) + log(1 - mu)),
                  Gamma    = sum(-y/(mu - log(mu))),
                  stop("Error: distribution not recognized"))

  ## Fit model with independence correlation structure
  object$call$corstr <- "independence"

  model.indep <- eval(object$call, parent.frame())
  ## model.indep <- eval(object$call, parent.frame()) ## FIXME parent.frame() is wrong...
  ## model.indep <- update(object, corstr="independence",zcorr=NULL)

  # Trace term (penalty for model complexity)
  AIinverse <- MASS::ginv(model.indep$naive.variance)
  Vr <- object$robust.variance
  trace <- sum(diag(AIinverse %*% Vr))
  params <- length(coef(object)) # Mean parameters in the model

  kpm <- params+length(object$scale)

  # QIC
  QIC  <- -2*(quasi - trace)
  QICu <- -2*(quasi - params)
  QICC <- QIC + (2 * kpm * (kpm + 1)) / (length(unique(object$id)) - kpm - 1)
  output <- c(QIC, QICu, quasi, trace, params, QICC)
  names(output) <- c("QIC", "QICu", "Quasi Lik", "CIC", "params", "QICC")
  output
}

