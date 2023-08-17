# Rfun_szgaViz
# based on Rcode_112_design, Rcode_112_design_func
# R function: marginalPower1, marginalPower2

#' @name szgaViz
#' @title Sample size optimization using graphical approach in clinical trial design with two hypotheses
#' @description This function computes the optimal design using graphical approach along with the minimum sample size when two hypotheses are considered in a clinical trial.
#'
#' @param alpha a value of overall type I error rate
#' @param beta1 a value of one minus marginal powers for testing H1
#' @param beta2 a value of one minus marginal powers for testing H2
#' @param deltaVec a vector of effect sizes for testing H1 nd H2, respectively
#' @param cVec a vector of coefficients. When testing continuous endpoints, these coefficients are exactly one. When testing binary endpoints, the values are roughly one but not exactly one
#' @param rho a value of correlation coefficients between two hypotheses
#' @param wunit a value of initial weight on H1 for grid search and visualization
#' @param initIntvl a vector of lower and upper limits for searching optimal sample size
#' @param visualization a logical value, indicating whether a visualization is needed
#' @return a vector of three numbers: the optimal weight on H1 \code{w1}, and optimal sample size \code{n1} (based on H1) and \code{n2} (based on H2), where \code{n1} and \code{n2} should be roughly the same
#' @export
#' @import stats
#' @import mvtnorm
#' @import graphics
#' @author Jiangtao Gou
#' @author Fengqing (Zoe) Zhang
#' @references
#'  Zhang, F. and Gou, J. (2023). Sample size optimization for clinical trials using graphical approaches for multiplicity adjustment, Technical Report.
#'  Gou, J. (2022). Sample size optimization and initial allocation of the significance levels in group sequential trials with multiple endpoints. \emph{Biometrical Journal}, 64(2), 301-311.
#' @examples
#' szgaViz(alpha = 0.05, beta1 = 0.20, beta2 = 0.20,
#'          deltaVec = c(0.3,0.3), cVec = c(1,1), rho = 0.0,
#'          wunit= 0.01, initIntvl = c(1,1000),
#'          visualization = FALSE)
#'
szgaViz <- function (alpha, beta1, beta2, deltaVec, cVec, rho, wunit, initIntvl, visualization = TRUE) {
  #
  # Rfun_pwr2
  # based on Rcode_119_func
  # package: mvtnorm
  # R function: none
  # Marginal Power Functions for two hypotheses
  #
  marginalPower1 <- function (w, n, alpha, rho, deltaVec, cVec) {
    c1 <- cVec[1]
    c2 <- cVec[2]
    Delta1 <- deltaVec[1] * sqrt(n)
    Delta2 <- deltaVec[2] * sqrt(n)
    #
    zcv1 <- stats::qnorm(p = w*alpha, mean = 0, sd = 1, lower.tail = FALSE)
    zcv2 <- stats::qnorm(p = (1-w)*alpha, mean = 0, sd = 1, lower.tail = FALSE)
    zcv <- stats::qnorm(p = alpha, mean = 0, sd = 1, lower.tail = FALSE)
    #
    partA <- stats::pnorm(q = c1*zcv1 - Delta1, mean = 0, sd = 1, lower.tail = FALSE)
    #
    lowerB <- c(c1*zcv, c2*zcv2)
    upperB <- c(c1*zcv1, +1000)
    meanV <- c(Delta1, Delta2)
    corrM <- matrix(c(1,rho,rho,1), nrow = 2)
    resultIntgl <- mvtnorm::pmvnorm(lowerB, upperB, meanV, corrM, algorithm=Miwa(steps=128))
    partB <- resultIntgl[1]
    #
    return(partA+partB)
  }
  #
  marginalPower2 <- function (w, n, alpha, rho, deltaVec, cVec) {
    c1 <- cVec[1]
    c2 <- cVec[2]
    Delta1 <- deltaVec[1] * sqrt(n)
    Delta2 <- deltaVec[2] * sqrt(n)
    #
    zcv1 <- stats::qnorm(p = w*alpha, mean = 0, sd = 1, lower.tail = FALSE)
    zcv2 <- stats::qnorm(p = (1-w)*alpha, mean = 0, sd = 1, lower.tail = FALSE)
    zcv <- stats::qnorm(p = alpha, mean = 0, sd = 1, lower.tail = FALSE)
    #
    partA <- stats::pnorm(q = c2*zcv2 - Delta2, mean = 0, sd = 1, lower.tail = FALSE)
    #
    lowerB <- c(c1*zcv1, c2*zcv)
    upperB <- c(1000, c2*zcv2)
    meanV <- c(Delta1, Delta2)
    corrM <- matrix(c(1,rho,rho,1), nrow = 2)
    resultIntgl <- mvtnorm::pmvnorm(lowerB, upperB, meanV, corrM, algorithm=Miwa(steps=128))
    partB <- resultIntgl[1]
    #
    return(partA+partB)
  }
  #
  # target1 and target2
  target1 <- function (n, w, alpha, beta1, rho, deltaVec, cVec) {
    mp1 <- marginalPower1(w, n, alpha, rho, deltaVec, cVec)
    return (mp1 - 1 + beta1)
  }
  #
  target2 <- function (n, w, alpha, beta2, rho, deltaVec, cVec) {
    mp2 <- marginalPower2(w, n, alpha, rho, deltaVec, cVec)
    return (mp2 - 1 + beta2)
  }
  #
  wvec <- seq(from = wunit, to = 1 - wunit, by = wunit)
  wlength <- length(wvec)
  #
  n1vec <- rep(0, times = wlength)
  n2vec <- rep(0, times = wlength)
  #
  for (i in 1:wlength) {
    result1 <- stats::uniroot(target1, interval=initIntvl, tol=2.5e-16, w = wvec[i], alpha=alpha, beta1 = beta1, rho = rho, deltaVec=deltaVec, cVec=cVec)
    result2 <- stats::uniroot(target2, interval=initIntvl, tol=2.5e-16, w = wvec[i], alpha=alpha, beta2 = beta2, rho = rho, deltaVec=deltaVec, cVec=cVec)
    n1vec[i] <- result1$root
    n2vec[i] <- result2$root
  }
  if (visualization) {
    graphics::plot(wvec, n1vec, type = "l", lty = 1)
    graphics::lines(wvec, n2vec, lty = 4)
  }
  idx <- which(abs(n1vec-n2vec) == min(abs(n1vec-n2vec)))
  output <- c(wvec[idx], n1vec[idx], n2vec[idx])
  names(output) <- c("w1", "n1", "n2")
  return(output)
}
