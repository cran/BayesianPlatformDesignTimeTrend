#' @title demo_Cutoffscreening
#' @description This function does a cutoff screening for trial simulation.
#' @param ntrials A numeric variable indicating how many trial replicates you want to run
#' @param trial.fun The function of trial simulation, related to MainFunction.R
#' @param grid.inf A list of grid information to create start grid and extend grid for cutoff screening.
#' @param input.info A list of input information including all information required for trial simulation.
#' @param cl A numeric variable indicating how many cores you want to use in parallel programming.
#'
#' @return A vector of recommended cutoff. The final value is the latest recommended value. A plot for all tested cutoff and error rate
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom iterators icount
#' @importFrom stats lm
#' @importFrom stats predict
#' @importFrom rstan rstan_options
#' @importFrom RcppParallel RcppParallelLibs CxxFlags
#' @importFrom graphics lines
#' @import BiocManager
#' @export
#'
#' @examples
#' \donttest{demo_Cutoffscreening(ntrials = 2, cl = 2,
#'     grid.inf = list(start = c(0.9, 0.95, 1), extendlength = 2))}
#' @author Ziyan Wang
demo_Cutoffscreening = function(ntrials = 1000,
                                trial.fun = simulatetrial,
                                grid.inf = list(start = c(0.9, 0.95, 1), extendlength =
                                                  15),
                                input.info = list(
                                  response.probs = c(0.4, 0.4),
                                  ns = c(30, 60, 90, 120, 150),
                                  max.ar = 0.75,
                                  rand.algo = "Urn",
                                  max.deviation = 3,
                                  model.inf = list(
                                    model = "tlr",
                                    ibb.inf = list(
                                      pi.star = 0.5,
                                      pess = 2,
                                      betabinomialmodel = ibetabinomial.post
                                    ),
                                    tlr.inf = list(
                                      beta0_prior_mu = 0,
                                      beta1_prior_mu = 0,
                                      beta0_prior_sigma = 2.5,
                                      beta1_prior_sigma = 2.5,
                                      beta0_df = 7,
                                      beta1_df = 7,
                                      reg.inf =  "main",
                                      variable.inf = "Fixeffect"
                                    )
                                  ),
                                  Stop.type = "Early-Pocock",
                                  Boundary.type = "Symmetric",
                                  Random.inf = list(
                                    Fixratio = FALSE,
                                    Fixratiocontrol = NA,
                                    BARmethod = "Thall",
                                    Thall.tuning.inf = list(tuningparameter = "Fixed",  fixvalue = 1)
                                  ),
                                  trend.inf = list(
                                    trend.type = "step",
                                    trend.effect = c(0, 0),
                                    trend_add_or_multip = "mult"
                                  )
                                ),
                                cl = 2) {
  old <- options()# code line i
  on.exit(options(old))
  #Set start grid of screening
  startgrid <-
    data.frame(tpIE = rep(NA, length(grid.inf$start)), cutoff = grid.inf$start)

  extendgrid <-
    data.frame(
      tpIE = rep(NA, grid.inf$extendlength),
      cutoff = rep(NA, grid.inf$extendlength)
    )

  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores(logical = FALSE))

  registerDoParallel(cores = cl)
  message("Start the start grid screening")
  for (j in 1:dim(startgrid)[1]) {
    #Construct the stop.inf list
    Stopbound.inf = Stopboundinf(
      Stop.type = input.info$Stop.type,
      Boundary.type = input.info$Boundary.type,
      cutoff = c(startgrid[j, 2], 1 - startgrid[j, 2])
    )

    result = foreach(icount(ntrials)) %dopar% trial.fun(
      response.probs = input.info$response.probs,
      ns = input.info$ns,
      max.ar = input.info$max.ar,
      rand.algo = input.info$rand.algo,
      max.deviation = input.info$max.deviation,
      model.inf = input.info$model.inf,
      Stopbound.inf = Stopbound.inf,
      Random.inf = input.info$Random.inf,
      trend.inf = input.info$trend.inf
    )

    # perHtypeIerror=mean(perHtypeIerrorfunc(result))
    FWER = FWER_disconjunctivepowerfunc(result)
    startgrid[j, 1] = FWER
  }

  startgrid$cutoff2 <- startgrid$cutoff ^ 2
  quadratic.model <-
    lm(tpIE ~ cutoff + cutoff2, data = data.frame(startgrid))
  cutoffgrid <- seq(0.9, 1, 0.0001)
  predictedtpIE <-
    predict(quadratic.model,
            list(cutoff = cutoffgrid, cutoff2 = cutoffgrid ^ 2))
  # plot(tpIE~cutoff, pch=16, xlab = "cutoff", ylab = "tpIE", cex.lab = 1.3, col = "blue",data = data.frame(startgrid))
  # lines(cutoffgrid, predictedtpIE, col = "darkgreen", lwd = 3)
  potentialcutoff = cutoffgrid[abs(predictedtpIE - 0.05) <= 0.0025]
  e = 1e-10
  randomprobability = (1 / (abs(predictedtpIE[abs(predictedtpIE - 0.05) <=
                                                0.0025] - 0.05) + e)) / sum(1 / (abs(predictedtpIE[abs(predictedtpIE - 0.05) <=
                                                                                                     0.0025] - 0.05) + e))
  nextcutoff = sample(potentialcutoff, 1, replace = T, prob = randomprobability)
  extendgrid[1, 2] = nextcutoff
  recommand = {

  }
  message(paste("Start the extend grid screening.","There are", grid.inf$extendlength ,"cutoff values under investigation in the extend grid"))
  for (cutoffindex in 1:(dim(extendgrid)[1])) {
    #Construct the stop.inf list
    Stopbound.inf = Stopboundinf(
      Stop.type = input.info$Stop.type,
      Boundary.type = input.info$Boundary.type,
      cutoff = c(extendgrid[cutoffindex, 2], 1 - extendgrid[cutoffindex, 2])
    )
    restlr090five = foreach(icount(ntrials)) %dopar% trial.fun(
      response.probs = input.info$response.probs,
      ns = input.info$ns,
      max.ar = input.info$max.ar,
      rand.algo = input.info$rand.algo,
      max.deviation = input.info$max.deviation,
      model.inf = input.info$model.inf,
      Stopbound.inf = Stopbound.inf,
      Random.inf = input.info$Random.inf,
      trend.inf = input.info$trend.inf
    )
    FWER = FWER_disconjunctivepowerfunc(restlr090five)
    extendgrid[cutoffindex, 1] = FWER
    extendgrid$cutoff2 <- extendgrid$cutoff ^ 2
    quadratic.model <-
      lm(tpIE ~ cutoff + cutoff2, data = data.frame(rbind(startgrid, extendgrid)))
    cutoffgrid <- seq(0.9, 1, 0.0001)
    predictedtpIE <-
      predict(quadratic.model,
              list(cutoff = cutoffgrid, cutoff2 = cutoffgrid ^ 2))
    # plot(tpIE~cutoff, pch=16, xlab = "cutoff", ylab = "tpIE", cex.lab = 1.3, col = "blue",data = data.frame(rbind(startgrid,extendgrid)))
    # lines(cutoffgrid, predictedtpIE, col = "darkgreen", lwd = 3)
    potentialcutoff = cutoffgrid[abs(predictedtpIE - 0.05) <= 0.0025]
    randomprobability = (1 / (abs(predictedtpIE[abs(predictedtpIE - 0.05) <=
                                                  0.0025] - 0.05) + e)) / sum(1 / (abs(predictedtpIE[abs(predictedtpIE - 0.05) <=
                                                                                                       0.0025] - 0.05) + e))
    if (length(potentialcutoff) == 0) {
      randomprobability = 1
      potentialcutoff = extendgrid[cutoffindex, 2]
    }
    extendgrid[cutoffindex + 1, 2] = sample(potentialcutoff, 1, replace = T, prob = randomprobability)
    recommand = c(recommand, cutoffgrid[as.numeric(names(which.max(randomprobability)))])
    message(paste("Finished extend grid screening round", cutoffindex))
  }
  message("Output data recording")
  dataloginformd = data.frame(rbind(startgrid, extendgrid))
  recommandloginformd = recommand
  quadratic.model <-
    lm(tpIE ~ cutoff + cutoff2, data = dataloginformd)
  cutoffgrid <- seq(0.9, 1, 0.0001)
  predictedtpIEinformd <-
    predict(quadratic.model,
            list(cutoff = cutoffgrid, cutoff2 = cutoffgrid ^ 2))
  doParallel::stopImplicitCluster()
  return(
    list(
      detailsforgrid = dataloginformd,
      recommandcutoff = recommandloginformd,
      predictedtpIEinformd = predictedtpIEinformd
    )
  )
}
