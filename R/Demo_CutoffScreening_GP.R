#' @title A demo for cutoff screening using Bayesian optimisation
#' @description This function does a cutoff screening for trial simulation using Bayesian optimisation.
#' @param ntrials A numeric variable indicating how many trial replicates you want to run
#' @param trial.fun The function of trial simulation, related to MainFunction.R
#' @param grid.inf A list of grid information to create start grid and extend grid for cutoff screening.
#'   'start.length' is the size of start grid. Default is 10. 'grid.min' is the lower bound for screening grid.
#'   'grid.max' is the upper bound for screening grid. 'errorrate' refers to the target of type I error rate or family-wise error rate.
#'   'confidence.level' is a numeric value indicating the confidence level of estimate. Default is 0.95.
#'   'grid.length' is the accuracy of grid. Default is 5000.
#'   'change.scale' is a logic value indicating whether we want to change scale when doing Gaussian process. Default is FALSE.
#'   'noise' is a logic value indicating whether the input x is noisy. Default is TRUE.
#'   'simulationerror' is a numeric value indicating the tolerable error for simulated type I error rate. Default is 0.01,
#'   'iter.max' is a numeric value indicating the maximum number of evaluations. Default is 15.
#'   'plotornot' is a logic value indicating whether the errorrate vs grid plot needed to be generated at each iteration. Default is FALSE.
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
#' \donttest{demo_Cutoffscreening.GP(ntrials = 2, cl = 2,
#'     grid.inf = list(
#'     start.length = 10,
#'     confidence.level = 0.95,
#'     grid.length = 5000,
#'     change.scale = FALSE,
#'     noise = TRUE,
#'     errorrate = 0.1,
#'     simulationerror = 0.01,
#'     iter.max = 15,
#'     plotornot = FALSE))}
#' @author Ziyan Wang
demo_Cutoffscreening.GP = function(ntrials = 1000,
                                   trial.fun = simulatetrial,
                                   grid.inf = list(
                                     start.length = 10,
                                     grid.min = NULL,
                                     grid.max = NULL,
                                     confidence.level = 0.95,
                                     grid.length =
                                       5000,
                                     change.scale = FALSE,
                                     noise = TRUE,
                                     errorrate = 0.1,
                                     simulationerror = 0.01,
                                     iter.max = 15,
                                     plotornot = FALSE
                                   ),
                                   input.info = list(
                                     response.probs = c(0.15, 0.15, 0.15, 0.15),
                                     ns = c(120, 240, 360, 480, 600),
                                     max.ar = 0.85,
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
                                     Stop.type = "Early-OBF",
                                     Boundary.type = "Symmetric",
                                     Random.inf = list(
                                       Fixratio = FALSE,
                                       Fixratiocontrol = NA,
                                       BARmethod = "Thall",
                                       Thall.tuning.inf = list(tuningparameter = "Unfixed",  fixvalue = 1)
                                     ),
                                     trend.inf = list(
                                       trend.type = "step",
                                       trend.effect = c(0, 0, 0, 0),
                                       trend_add_or_multip = "mult"
                                     )
                                   ),
                                   cl = 2) {
  old <- options()# code line i
  on.exit(options(old))
  Boundary.type=input.info$Boundary.type
  #Set start grid of screening
  if (is.null(grid.inf$grid.min) | is.null(grid.inf$grid.max)) {
    if (input.info$Stop.type == "Early-OBF") {
      grid.min = 1
      grid.max = 8
    }
    else{
      grid.min = 0.9
      grid.max = 1
    }
  }
    else if (grid.inf$grid.min >= grid.inf$grid.max) {
      stop("Error: grid.min should be greater tha grid.max")
    }
    else{
      grid.min = grid.inf$grid.min
      grid.max = grid.inf$grid.max
    }

  if (Boundary.type == "Symmetric") {
    start = lhs::maximinLHS(grid.inf$start.length, 1) * (grid.max - grid.min) +
      grid.min
    startgrid <-
      data.frame(tpIE = rep(NA, grid.inf$start.length),
                 cutoff = start)

    extendgrid <-
      data.frame(
        tpIE = rep(NA, grid.inf$iter.max),
        cutoff = rep(NA, grid.inf$iter.max)
      )
  }
  else{
    start = lhs::maximinLHS(grid.inf$start.length, 2) * (grid.max - grid.min) +
      grid.min
    start = start[start[, 1] < start[, 2],]
    startgrid <-
      data.frame(
        tpIE = rep(NA, nrow(start)),
        cutoffdown = start[, 1],
        cutoffup = start[, 2]
      )

    extendgrid <-
      data.frame(
        tpIE = rep(NA, grid.inf$iter.max),
        cutoffdown = rep(NA, grid.inf$iter.max),
        cutoffup = rep(NA, grid.inf$iter.max)
      )
  }
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores(logical = FALSE))

  registerDoParallel(cores = cl)
  message("Start the start grid screening")

  for (j in 1:dim(startgrid)[1]) {
    #Construct the stop.inf list
    if (Boundary.type == "Symmetric") {
      if (input.info$Stop.type == "Early-OBF") {
        cutoff = c(startgrid[j, 2], startgrid[j, 2])
      }
      else{
        cutoff = c(startgrid[j, 2], 1 - startgrid[j, 2])
      }
    }
    else{
      cutoff = c(startgrid[j, 2], startgrid[j, 3])
    }
    Stopbound.inf = Stopboundinf(
      Stop.type = input.info$Stop.type,
      Boundary.type = input.info$Boundary.type,
      cutoff = cutoff
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

  GP.res = GP.optim(startgrid$cutoff, startgrid$tpIE, grid.min=grid.min, grid.max=grid.max)
  nextcutoff = GP.res$next.cutoff
  if (Boundary.type == "Symmetric") {
    extendgrid[1, 2] = nextcutoff
  }
  else{
    extendgrid[1, 2] = nextcutoff[1]
    extendgrid[1, 3] = nextcutoff[2]
  }
  message(
    paste(
      "Start the extend grid screening.",
      "There are at most",
      grid.inf$iter.max ,
      "cutoff values under evaluation in the extend grid"
    )
  )
  for (cutoffindex in 1:(dim(extendgrid)[1])) {
    if (Boundary.type == "Symmetric") {
      if (input.info$Stop.type == "Early-OBF") {
        cutoff = c(extendgrid[cutoffindex, 2], extendgrid[cutoffindex, 2])
      }
      else{
        cutoff = c(extendgrid[cutoffindex, 2], 1 - extendgrid[cutoffindex, 2])
      }
    }
    else{
      cutoff = c(extendgrid[cutoffindex, 2], extendgrid[cutoffindex, 3])
    }
    #Construct the stop.inf list
    Stopbound.inf = Stopboundinf(
      Stop.type = input.info$Stop.type,
      Boundary.type = input.info$Boundary.type,
      cutoff = cutoff
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

    if (FWER <= grid.inf$errorrate * (1 + grid.inf$simulationerror) &
        FWER >= grid.inf$errorrate * (1 - grid.inf$simulationerror)) {
      earlyend = TRUE
    }
    else{
      earlyend = FALSE
    }
    extendgrid[cutoffindex, 1] = FWER

    GP.res = GP.optim(startgrid$cutoff, startgrid$tpIE, grid.min=grid.min, grid.max=grid.max)
    nextcutoff = GP.res$next.cutoff
    if (Boundary.type == "Symmetric") {
      extendgrid[cutoffindex + 1, 2] = nextcutoff
    }
    else{
      extendgrid[cutoffindex + 1, 2] = nextcutoff[1]
      extendgrid[cutoffindex + 1, 3] = nextcutoff[2]
    }
    prediction = data.frame(
      yhat = GP.res$prediction$yhat,
      sd = GP.res$prediction$sd,
      qup = GP.res$prediction$qup,
      qdown = GP.res$prediction$qdown,
      xgrid = GP.res$prediction$xgrid
    )
    # Plot or not
    if (isTRUE(grid.inf$plotornot)) {
      if (Boundary.type == "Symmetric") {
        GPplot = ggplot(data = prediction) +
          geom_ribbon(aes(
            x = prediction$xgrid,
            ymin = prediction$qdown,
            ymax = prediction$qup
          ), alpha = 0.5) +
          geom_line(ggplot2::aes(prediction$xgrid, prediction$yhat)) +
          geom_point(ggplot2::aes(cutoff[1:sum(!is.na(prediction$tpIE))], prediction$tpIE[1:sum(!is.na(prediction$tpIE))]),
                              data = data.frame(rbind(startgrid, extendgrid))) +
          geom_hline(yintercept = grid.inf$errorrate) +
          geom_text(aes(x=grid.min,y=grid.inf$errorrate+0.05,label=paste0("FWER target is ",grid.inf$errorrate)),hjust=0,vjust=1)+
          geom_vline(xintercept = nextcutoff, linetype = 2) +
          geom_text(aes(x=nextcutoff,y=grid.inf$errorrate*2,label=paste0("Next cutoff value is ",round(nextcutoff,3))))+
          theme_minimal() +ylab("FWER")+xlab("Cutoff grid")+
          labs(title = paste0("Iteration", cutoffindex))
      }
      # Debug here for asymmetric boundary on April 30, 2023
      else if (Boundary.type != "Symmetric") {
        stop("Error: Plot is not available for asymmetric boundary now")
      }
    }
    # Stop the iteration
    if (isTRUE(earlyend)) {
      message(paste("Finished extend grid screening round", cutoffindex))
      message("The iteration stopped early because the optimal cutoff value is found.")
      break
    }
    else{
      message(paste("Finished extend grid screening round", cutoffindex))
    }
  }
  message("Output data recording")
  dataloginformd = data.frame(rbind(startgrid, extendgrid))
  doParallel::stopImplicitCluster()
  return(list(
    detailsforgrid = dataloginformd,
    predictedtpIEinformd = GP.res$prediction
  ))
}
