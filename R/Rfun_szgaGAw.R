# Rfun_szgaGAw
# based on Rcode_223_optGAw
# package: GA
# R function: marginalPower31, marginalPower32, marginalPower33

#' @name szgaGAw
#' @title Sample size optimization using graphical approach in clinical trial design with three hypotheses when the transition matrix is pre-specified
#' @description This function computes the optimal design using graphical approach along with the minimum sample size when three hypotheses are considered in a clinical trial. The transition matrix is pre-specified and fixed.
#'
#' @param alpha a value of overall type I error rate
#' @param betaVec a vector of one minus marginal powers for testing H1, H2 and H3, respectively
#' @param deltaVec a vector of effect sizes for testing H1, H2 and H3, respectively
#' @param cVec a vector of coefficients. When testing continuous endpoints, these coefficients are exactly one. When testing binary endpoints, the values are roughly one but not exactly one
#' @param rhoMat a matrix of the correlation coefficients among three hypotheses
#' @param transMat a matrix of the fixed transition probabilities among three hypotheses
#' @param lower a vector of lower limit of sample size n, and initial weights w1 and w2, where w3 is computed by 1 - w1 - w2
#' @param upper a vector of upper limit of sample size n, and initial weights w1 and w2, where w3 is computed by 1 - w1 - w2
#' @param gaIter a vector of two numbers. The first one is the parameter maxiter of the ga function, and the second one is the parameter run of the ga function
#' @param penPara a number of penalization parameter for optimization to balance the sample size requirement and the power requirement
#' @param seed a number of the seed of the random number generator
#' @return a vector of three numbers: the optimal sample size \code{n}, and initial weights \code{w1} and \code{w2}
#' @export
#' @importFrom GA ga
#' @import stats
#' @import mvtnorm
#' @author Jiangtao Gou
#' @details R package \code{GA} is used for Genetic Algorithms.
#' @references
#'  Zhang, F. and Gou, J. (2023). Sample size optimization for clinical trials using graphical approaches for multiplicity adjustment, Technical Report.
#'  Gou, J. (2022). Sample size optimization and initial allocation of the significance levels in group sequential trials with multiple endpoints. \emph{Biometrical Journal}, 64(2), 301-311.
#' @examples
#' start <- Sys.time()
#' szgaGAw(alpha = 0.025, betaVec = c(0.15, 0.20, 0.10),
#'         deltaVec = c(0.1111952, 0.1037179, 0.1335865),
#'         cVec = c(1.003086, 1.002686, 1.004451),
#'         rhoMat = matrix(c(1,0.5,0.8, 0.5,1,0.6, 0.8,0.6,1), nrow = 3, byrow = TRUE),
#'         transMat = matrix(c(0,0.50,0.50, 0.5,0,0.5, 0.5,0.5,0), nrow = 3, byrow = TRUE),
#'         lower = c(700, rep(0.05, 2)),
#'         upper = c(900, rep(0.95, 2)),
#'         gaIter = c(10, 5),
#'         penPara = 0.0135,
#'         seed = 234)
#' end <- Sys.time()
#' data.frame(time = end - start)
#'
szgaGAw <- function (alpha, betaVec, deltaVec, cVec, rhoMat, transMat, lower = c(1, rep(1e-6, 2)), upper = c(1e4, rep(1-1e-6, 2)), gaIter = c(20, 20), penPara = 0.1, seed = 2022) {
  #
  # Rfun_pwr
  # based on Rcode_211_pwr
  # package: mvtnorm
  # R function: none
  # Marginal Power Functions for three hypotheses
  #
  marginalPower31 <- function (wVec, transMat, n, alpha, rhoMat, deltaVec, cVec){
    DeltaVec <- deltaVec * sqrt(n)
    #
    PartA <- stats::pnorm(cVec[1] * stats::qnorm(p = wVec[1]*alpha, lower.tail = FALSE),
                          mean = DeltaVec[1], sd = 1, lower.tail = FALSE)
    #
    lowerB <- c(
      cVec[1] * stats::qnorm(p = (wVec[1] + wVec[2]*transMat[2,1])*alpha, lower.tail = FALSE),
      cVec[2] * stats::qnorm(p = wVec[2]*alpha, lower.tail = FALSE)
    )
    upperB <- c(
      cVec[1] * stats::qnorm(p = wVec[1]*alpha, lower.tail = FALSE),
      +1000
    )
    meanV <- DeltaVec[c(1,2)]
    corrM <- rhoMat[c(1,2), c(1,2)]
    resultIntgl <- mvtnorm::pmvnorm(lowerB, upperB, meanV, corrM, algorithm = Miwa(steps=128))
    PartB <- resultIntgl[1]
    #
    lowerB <- c(
      cVec[1] * stats::qnorm(p = alpha, lower.tail = FALSE),
      cVec[2] * stats::qnorm(p = wVec[2]*alpha, lower.tail = FALSE),
      cVec[3] * stats::qnorm(p = (wVec[3] + wVec[2]*transMat[2,3])*alpha, lower.tail = FALSE)
    )
    upperB <- c(
      cVec[1] * stats::qnorm(p = (wVec[1] + wVec[2]*transMat[2,1])*alpha, lower.tail = FALSE),
      +1000,
      +1000
    )
    meanV <- DeltaVec
    corrM <- rhoMat
    resultIntgl <- mvtnorm::pmvnorm(lowerB, upperB, meanV, corrM, algorithm = Miwa(steps=128))
    PartC <- resultIntgl[1]
    #
    lowerB <- c(
      cVec[1] * stats::qnorm(p = (wVec[1] + wVec[3]*transMat[3,1])*alpha, lower.tail = FALSE),
      -1000,
      cVec[3] * stats::qnorm(p = wVec[3]*alpha, lower.tail = FALSE)
    )
    upperB <- c(
      cVec[1] * stats::qnorm(p = wVec[1]*alpha, lower.tail = FALSE),
      cVec[2] * stats::qnorm(p = wVec[2]*alpha, lower.tail = FALSE),
      +1000
    )
    resultIntgl <- mvtnorm::pmvnorm(lowerB, upperB, meanV, corrM, algorithm = Miwa(steps=128))
    PartD <- resultIntgl[1]
    #
    lowerB <- c(
      cVec[1] * stats::qnorm(p = alpha, lower.tail = FALSE),
      cVec[2] * stats::qnorm(p = (wVec[2] + wVec[3]*transMat[3,2])*alpha, lower.tail = FALSE),
      cVec[3] * stats::qnorm(p = wVec[3]*alpha, lower.tail = FALSE)
    )
    upperB <- c(
      cVec[1] * stats::qnorm(p = (wVec[1] + wVec[3]*transMat[3,1])*alpha, lower.tail = FALSE),
      cVec[2] * stats::qnorm(p = wVec[2]*alpha, lower.tail = FALSE),
      +1000
    )
    resultIntgl <- mvtnorm::pmvnorm(lowerB, upperB, meanV, corrM, algorithm = Miwa(steps=128))
    PartE <- resultIntgl[1]
    #
    return(sum(c(PartA, PartB, PartC, PartD, PartE)))
  }
  #
  marginalPower32 <- function (wVec, transMat, n, alpha, rhoMat, deltaVec, cVec){
    #twoFirst <- c(2,1,3)
    twoFirst <- c(2,3,1)
    wVec <- wVec[twoFirst]
    transMat <- transMat[twoFirst,twoFirst]
    rhoMat <- rhoMat[twoFirst,twoFirst]
    deltaVec <- deltaVec[twoFirst]
    cVec <- cVec[twoFirst]
    result <- marginalPower31(wVec, transMat, n, alpha, rhoMat, deltaVec, cVec)
    return(result)
  }
  #
  marginalPower33 <- function (wVec, transMat, n, alpha, rhoMat, deltaVec, cVec){
    twoFirst <- c(3,1,2)
    #twoFirst <- c(3,2,1)
    wVec <- wVec[twoFirst]
    transMat <- transMat[twoFirst,twoFirst]
    rhoMat <- rhoMat[twoFirst,twoFirst]
    deltaVec <- deltaVec[twoFirst]
    cVec <- cVec[twoFirst]
    result <- marginalPower31(wVec, transMat, n, alpha, rhoMat, deltaVec, cVec)
    return(result)
  }
  #
  target <- function (x, alpha, betaVec, deltaVec, cVec, rhoMat, transMat) {
    #
    # we need to maximize the function
    pen <- .Machine$double.xmax^penPara  # penalty term
    #
    n <- floor(x[1])
    wVec <- c(x[2], x[3], 1 - x[2] - x[3])
    penaltyWt <- max(-wVec[3], 0)*pen*100+1000 # 100 measure the difference of penalties on weights and powers
    if (wVec[3] < 0) {
      return (- penaltyWt)
    } # Three weights should all be positive
    #
    pwrCalc <- c(marginalPower31(wVec, transMat, n, alpha, rhoMat, deltaVec, cVec),
                 marginalPower32(wVec, transMat, n, alpha, rhoMat, deltaVec, cVec),
                 marginalPower33(wVec, transMat, n, alpha, rhoMat, deltaVec, cVec))
    pwrConstraint <- 1 - betaVec - pwrCalc
    penaltyPwr <- max(pwrConstraint,0)*pen # penalisation
    return(-n - penaltyPwr)
  } # End of function target
  #
  resultGA <- GA::ga("real-valued", fitness = target,
                     lower = lower, upper = upper,
                     maxiter = gaIter[1], run = gaIter[2], seed = seed,
                     alpha = alpha, betaVec = betaVec,
                     deltaVec = deltaVec, cVec = cVec, rhoMat = rhoMat, transMat = transMat)
  return(resultGA@solution)
}
