#' @title A demo for cutoff screening using Bayesian optimisation
#' @description This function does a cutoff screening for trial simulation using Bayesian optimisation.
#'
#' @param ntrials A numeric variable indicating how many trial replicates you want to run
#' @param trial.fun The function of trial simulation, related to MainFunction.R
#' @param grid.inf A list of grid information to create start grid and extend grid for cutoff screening.
#'   'start.length' is the size of start grid. Default is 10.
#'   'grid.min' A numeric value or vector (for asymmetric boundary) indicating the lower bound of the grid for screening. For asymmetric boundary, the first value is efficacy minimum value and the second value is futility minimum value.
#'   'grid.max' A numeric value or vector (for asymmetric boundary) indicating the upper bound of the grid for screening. For asymmetric boundary, the first value is efficacy maximum value and the second value is futility maximum value.
#'   'errorrate' refers to the target of type I error rate or family-wise error rate. Default is 0.05. User can change it to 0.1 for FWER if they think 0.05 is too conservative. The per-hypothesis type I error equals errorrate / (K-1) where (K-1) is the number of treatment arms.
#'   'confidence.level' is a numeric value indicating the confidence level of estimate. Default is 0.95.
#'   'grid.length' A numeric value indicating the grid resolution. Default is 5000 for symmetric boundary. For asymmetric boundary, the length of grid is 101 for both efficacy grid and futility grid. A numeric value indicating the grid resolution. Default is 5000 for symmetric boundary. For asymmetric boundary, the length of grid is 101 for both efficacy grid and futility grid.
#'   'change.scale' is a logic value indicating whether we want to change scale when doing Gaussian process. Default is FALSE.
#'   'noise' is a logic value indicating whether the input x is noisy. Default is TRUE.
#'   'simulationerror' is a numeric value indicating the tolerable error for simulated type I error rate. Default is 0.01,
#'   'iter.max' is a numeric value indicating the maximum number of evaluations. Default is 15.
#'   'plotornot' is a logic value indicating whether the errorrate vs grid plot needed to be generated at each iteration. Default is FALSE.
#' @param input.info A list of input information including all information required for trial simulation.
#' @param cl A numeric variable indicating how many cores you want to use in parallel programming.
#' @param power.type A indicator of which type of power we need to optimise when tuning the cutoff value for asymmetric boundary. Default is NA (Symmetric boundary).
#'   The choice of power type is Conjunctive power ("Conjunctive") and Disconjunctive power ("Disconjunctive"). In a two arm trial design, these power type are the same.
#' @param response.probs.alt A vector of response probability of each arm under the alternative scenario. This is used for power optimisation when tuning the cutoff values for asymmetric boundary. Default is NA.
#'
#' @return A vector of recommended cutoff. The final value is the latest recommended value. A plot for all tested cutoff and error rate
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom iterators icount
#' @importFrom stats lm dist runif
#' @importFrom stats predict
#' @importFrom rstan rstan_options
#' @importFrom RcppParallel RcppParallelLibs CxxFlags
#' @importFrom graphics lines
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics image par points
#' @importFrom ggpubr ggarrange
#' @importFrom RColorBrewer brewer.pal
#' @import BiocManager
#' @export
#'
#' @examples
#' \donttest{
#' #Two arm asymmetric boundary screening. Default is OBF boundary.
#' demo_Cutoffscreening.GP(ntrials = 2, cl = 2,
#'  power.type = NA,
#'  response.probs.alt = NA,
#'  grid.inf = list(
#'  start.length = 10,
#'  confidence.level = 0.95,
#'  grid.length = 5000,
#'  change.scale = FALSE,
#'  noise = TRUE,
#'  errorrate = 0.1,
#'  simulationerror = 0.01,
#'  iter.max = 15,
#'  plotornot = FALSE))
#'
#'  #Four arm asymmetric OBF boundary screening where conjunctive power is optimised.
#'  demo_Cutoffscreening.GP(ntrials = 2, cl = 2,
#'  power.type = "Conjunctive",
#'  response.probs.alt = c(0.4,0.6,0.6,0.4),
#'  grid.inf = list(
#'  start.length = 10,
#'  confidence.level = 0.95,
#'  grid.length = 101,
#'  change.scale = FALSE,
#'  noise = TRUE,
#'  errorrate = 0.1,
#'  simulationerror = 0.01,
#'  iter.max = 15,
#'  plotornot = FALSE))
#'  input.info = list(
#'  response.probs.null = c(0.4,0.4,0.4,0.4),
#'  ns = c(120, 240, 360, 480, 600),
#'  max.ar = 0.85,
#'  rand.algo = "Urn",
#'  max.deviation = 3,
#'  test.type = "Twoside",
#'  model.inf = list(
#'  model = "tlr",
#'  ibb.inf = list(
#'  pi.star = 0.5,
#'  pess = 2,
#'  betabinomialmodel = ibetabinomial.post
#'  ),
#'  tlr.inf = list(
#'  beta0_prior_mu = 0,
#'  beta1_prior_mu = 0,
#'  beta0_prior_sigma = 2.5,
#'  beta1_prior_sigma = 2.5,
#'  beta0_df = 7,
#'  beta1_df = 7,
#'  reg.inf =  "main",
#'  variable.inf = "Fixeffect"
#'  )
#'  ),
#'  Stop.type = "Early-OBF",
#'  Boundary.type = "Asymmetric",
#'  Random.inf = list(
#'  Fixratio = FALSE,
#'  Fixratiocontrol = NA,
#'  BARmethod = "Thall",
#'  Thall.tuning.inf = list(tuningparameter = "Unfixed",  fixvalue = 1)
#'  ),
#'  trend.inf = list(
#'  trend.type = "step",
#'  trend.effect = c(0, 0, 0, 0),
#'  trend_add_or_multip = "mult"
#'  )
#'  )
#'  }
#'
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
                                   power.type = NA,
                                   response.probs.alt = NA,
                                   input.info = list(
                                     response.probs.null = c(0.4,0.4,0.4,0.4),
                                     ns = c(120, 240, 360, 480, 600),
                                     max.ar = 0.85,
                                     rand.algo = "Urn",
                                     max.deviation = 3,
                                     test.type = "Twoside",
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
  grid.length=grid.inf$grid.length
  #Set start grid of screening
  if (Boundary.type == "Symmetric") {
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
    if (is.na(response.probs.alt) | input.info$response.probs.null[1] != response.probs.alt[1]) {
      stop("Error: For asymmetric boundary, please input one null scenario and one of alternative scenarios related to the null for type I error control and power optimisation")
    }

    # Define the number of samples and the number of variables
    num_samples <- grid.inf$start.length  # Change this to the desired number of samples
    num_variables <- 2  # Change this to the number of variables you have

    # Define the number of optimization iterations
    num_iterations <- 10000
    dis_vec<-{}
    # Create a function to calculate the minimum pairwise distance in the LHS sample
    calculate_min_pairwise_distance <- function(lhs_sample) {
      dist_matrix <- as.matrix(dist(lhs_sample))
      min_distance <- min(dist_matrix[upper.tri(dist_matrix)])
      return(min_distance)
    }

    # Initialize an LHS sample
    lhs_sample <- matrix(NA, nrow = num_samples, ncol = num_variables)

    if (is.null(grid.inf$grid.min[1]) | is.null(grid.inf$grid.max[1])|is.null(grid.inf$grid.min[1]) | is.null(grid.inf$grid.max[2])) {
      if (input.info$Stop.type == "Early-OBF") {
        grideff.min = 1
        grideff.max = 8
        gridfut.min = 1
        gridfut.max = 8
      }
      else{
        grideff.min = 0.9
        grideff.max = 1
        gridfut.min = 0
        gridfut.max = 0.1
      }
    }
    else if ((grid.inf$grid.min[1] >= grid.inf$grid.max[1]|
             grid.inf$grid.min[2] >= grid.inf$grid.max[2]) & input.info$Stop.type != "Early-OBF") {
      stop("Error: grid.min should be greater tha grid.max for boundary excepting for OBF")
    }
    else{
      grideff.min = grid.inf$grid.min[1]
      grideff.max = grid.inf$grid.max[1]
      gridfut.min = grid.inf$grid.min[2]
      gridfut.max = grid.inf$grid.max[2]
    }

    # Generate a set of random permutations for each variable within their respective ranges
    for (i in 1:num_variables) {
      if (i == 1) {
        lhs_sample[, i] <- runif(num_samples, min = grid.inf$grid.min[1], max = grid.inf$grid.max[1])
      } else if (i == 2) {
        lhs_sample[, i] <- runif(num_samples, min = grid.inf$grid.min[2], max = grid.inf$grid.max[2])
      }
    }

    # Calculate the initial minimum pairwise distance
    best_min_distance <- calculate_min_pairwise_distance(lhs_sample)
    best_lhs_sample <- lhs_sample
    lhs_sample.temp=lhs_sample
    # Perform optimization
    for (iteration in 1:num_iterations) {
      # Shuffle the columns of the LHS sample
      lhs_sample.temp[,1] <- lhs_sample[sample(num_samples), 1]
      lhs_sample.temp[,2] <- lhs_sample[sample(num_samples), 2]

      # Calculate the minimum pairwise distance of the shuffled sample
      current_min_distance <- calculate_min_pairwise_distance(lhs_sample.temp)
      dis_vec<-c(dis_vec,current_min_distance)
      # If the current minimum distance is greater, update the best sample
      if (current_min_distance > best_min_distance) {
        best_min_distance <- current_min_distance
        best_lhs_sample <- lhs_sample.temp
      }
    }
    startgrid.asy <- data.frame(tpIE =rep(NA,15),
                                cutoff.eff = best_lhs_sample[,1],
                                cutoff.fut = best_lhs_sample[,2])

    extendgrid <-
      data.frame(
        tpIE = rep(NA, grid.inf$iter.max),
        cutoff.eff = rep(NA, grid.inf$iter.max),
        cutoff.fut = rep(NA, grid.inf$iter.max)
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
      Stopbound.inf = Stopboundinf(
        Stop.type = input.info$Stop.type,
        Boundary.type = input.info$Boundary.type,
        cutoff = cutoff
      )
      result = foreach(icount(ntrials)) %dopar% trial.fun(
        response.probs = input.info$response.probs.null,
        ns = input.info$ns,
        max.ar = input.info$max.ar,
        rand.algo = input.info$rand.algo,
        max.deviation = input.info$max.deviation,
        model.inf = input.info$model.inf,
        Stopbound.inf = Stopbound.inf,
        Random.inf = input.info$Random.inf,
        trend.inf = input.info$trend.inf
      )
      FWER = conjuncativepower_or_FWER(result,input.info$response.probs.null,test.type = input.info$test.type)
      startgrid[j, 1] = FWER

    }
    else{
      # Record effective sample size data for each cutoff pair
      samplesize.start.null={}
      samplesize.start.alt={}
      samplesize.ext.null={}
      samplesize.ext.alt={}
      if (power.type %in% c("Conjunctive", "Disconjunctive")){
        conj.pow = {}
        disconj.pow = {}
        marg.pow = {}
      }
      else {
        stop("Error: Please specify a correct power type for asymmetric cutoff screening.")
      }
      cutoff = c(startgrid.asy[j, 2], startgrid.asy[j, 3])

    Stopbound.inf = Stopboundinf(
      Stop.type = input.info$Stop.type,
      Boundary.type = input.info$Boundary.type,
      cutoff = cutoff
    )
      result.null = foreach(icount(ntrials)) %dopar% trial.fun(
        response.probs = input.info$response.probs.null,
        ns = input.info$ns,
        max.ar = input.info$max.ar,
        rand.algo = input.info$rand.algo,
        max.deviation = input.info$max.deviation,
        model.inf = input.info$model.inf,
        Stopbound.inf = Stopbound.inf,
        Random.inf = input.info$Random.inf,
        trend.inf = input.info$trend.inf
      )
      FWER = conjuncativepower_or_FWER(result.null,input.info$response.probs.null,test.type = input.info$test.type)
      startgrid.asy[j, 1] = FWER
      samplesize.start.null = c(samplesize.start.null, sum(Nfunc(result.null)))

      result.alt = foreach(icount(ntrials)) %dopar% trial.fun(
        response.probs = response.probs.alt,
        ns = input.info$ns,
        max.ar = input.info$max.ar,
        rand.algo = input.info$rand.algo,
        max.deviation = input.info$max.deviation,
        model.inf = input.info$model.inf,
        Stopbound.inf = Stopbound.inf,
        Random.inf = input.info$Random.inf,
        trend.inf = input.info$trend.inf
      )
      conj.pow = c(conj.pow,conjuncativepower_or_FWER(result.alt,response.probs.alt,test.type = input.info$test.type))
      disconj.pow = c(disconj.pow, disconjunctivepowerfunc(result.alt))
      samplesize.start.alt = c(samplesize.start.alt, sum(Nfunc(result.alt)))
      # marg.pow = c(marg.pow)
    }
  }

  if (Boundary.type == "Symmetric") {
    GP.res = GP.optim(x = matrix(startgrid$cutoff,ncol=1), y.t1E = startgrid$tpIE,
                      grid.min=grid.min, grid.max=grid.max,errorrate = grid.inf$errorrate,
                      grid.length=grid.length,
                      Boundary.type = Boundary.type)
    nextcutoff = GP.res$next.cutoff
    extendgrid[1, 2] = nextcutoff
  }
  else{
    if (power.type == "Conjunctive"){
      startgrid.asy$pow = conj.pow
    }
    else{
      startgrid.asy$pow = disconj.pow
    }
    GP.res = GP.optim(x = cbind(startgrid.asy$cutoff.eff, startgrid.asy$cutoff.fut),
                      y.t1E = startgrid.asy$tpIE, y.pow = startgrid.asy$pow,
                      errorrate = grid.inf$errorrate,
                      grid.min=c(grideff.min,gridfut.min), grid.max=c(grideff.max,gridfut.max),
                      grid.length=grid.length,
                      Boundary.type = Boundary.type)
    nextcutoff = GP.res$next.cutoff
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

    #Construct the stop.inf list
    Stopbound.inf = Stopboundinf(
      Stop.type = input.info$Stop.type,
      Boundary.type = input.info$Boundary.type,
      cutoff = cutoff
    )
    result.null = foreach(icount(ntrials)) %dopar% trial.fun(
      response.probs = input.info$response.probs.null,
      ns = input.info$ns,
      max.ar = input.info$max.ar,
      rand.algo = input.info$rand.algo,
      max.deviation = input.info$max.deviation,
      model.inf = input.info$model.inf,
      Stopbound.inf = Stopbound.inf,
      Random.inf = input.info$Random.inf,
      trend.inf = input.info$trend.inf
    )
    FWER = conjuncativepower_or_FWER(result.null,input.info$response.probs.null,test.type = input.info$test.type)

    if (FWER <= grid.inf$errorrate * (1 + grid.inf$simulationerror) &
        FWER >= grid.inf$errorrate * (1 - grid.inf$simulationerror)) {
      earlyend = TRUE
    }
    else{
      earlyend = FALSE
    }
    extendgrid[cutoffindex, 1] = FWER

    dnew=data.frame(rbind(startgrid,extendgrid))
    dnew=dnew[!is.na(dnew$cutoff),]
    GP.res = GP.optim(dnew$cutoff, dnew$tpIE, grid.min=grid.min, grid.max=grid.max,
                      grid.length=grid.length,
                      errorrate = grid.inf$errorrate)
    nextcutoff = GP.res$next.cutoff
    extendgrid[cutoffindex + 1, 2] = nextcutoff
    prediction = data.frame(
      yhat.t1E = GP.res$prediction$yhat.t1E,
      sd.t1E = GP.res$prediction$sd.t1E,
      qup.t1E = GP.res$prediction$qup.t1E,
      qdown.t1E = GP.res$prediction$qdown.t1E,
      xgrid = GP.res$prediction$xgrid
    )
    # Plot or not
    if (isTRUE(grid.inf$plotornot)) {
        GPplot = ggplot(data = prediction) +
          geom_ribbon(aes(
            x = prediction$xgrid,
            ymin = prediction$qdown.t1E,
            ymax = prediction$qup.t1E
          ),col="#f8766d", alpha = 0.5,linetype = 2) +
          geom_line(ggplot2::aes(prediction$xgrid, prediction$yhat.t1E),col = "#f8766d") +
          geom_point(ggplot2::aes(cutoff[1:sum(!is.na(prediction$tpIE))], prediction$tpIE[1:sum(!is.na(prediction$tpIE))]),
                     data = data.frame(rbind(startgrid, extendgrid)),col = "#00bfc4") +
          geom_hline(yintercept = grid.inf$errorrate) +
          geom_text(aes(x=grid.min,y=grid.inf$errorrate+0.05,label=paste0("FWER target is ",grid.inf$errorrate)),hjust=0,vjust=1)+
          geom_vline(xintercept = nextcutoff, linetype = 2) +
          geom_text(aes(x=nextcutoff,y=grid.inf$errorrate*2,label=paste0("Next cutoff value is ",round(nextcutoff,3))))+
          theme_minimal() +ylab("FWER")+xlab("Cutoff grid")+
          geom_point(aes(nextcutoff, 0.1),
                     data = data.frame(tpIE=prediction$tpIE,cutoff=cutoff),col = "#f8766d") +
          theme(plot.background = element_rect(fill = "#e6dfba"))
        labs(title = paste0("Iteration", cutoffindex))
      }
    }
    else{
      cutoff = c(extendgrid[cutoffindex, 2], extendgrid[cutoffindex, 3])
      #Construct the stop.inf list
      Stopbound.inf = Stopboundinf(
        Stop.type = input.info$Stop.type,
        Boundary.type = input.info$Boundary.type,
        cutoff = cutoff
      )
      result.null = foreach(icount(ntrials)) %dopar% trial.fun(
        response.probs = input.info$response.probs.null,
        ns = input.info$ns,
        max.ar = input.info$max.ar,
        rand.algo = input.info$rand.algo,
        max.deviation = input.info$max.deviation,
        model.inf = input.info$model.inf,
        Stopbound.inf = Stopbound.inf,
        Random.inf = input.info$Random.inf,
        trend.inf = input.info$trend.inf
      )
      FWER = conjuncativepower_or_FWER(result.null,input.info$response.probs.null,test.type = input.info$test.type)
      extendgrid[cutoffindex, 1] = FWER
      samplesize.ext.null = c(samplesize.ext.null, sum(Nfunc(result.null)))

      result.alt = foreach(icount(ntrials)) %dopar% trial.fun(
        response.probs = response.probs.alt,
        ns = input.info$ns,
        max.ar = input.info$max.ar,
        rand.algo = input.info$rand.algo,
        max.deviation = input.info$max.deviation,
        model.inf = input.info$model.inf,
        Stopbound.inf = Stopbound.inf,
        Random.inf = input.info$Random.inf,
        trend.inf = input.info$trend.inf
      )
      conj.pow = conjuncativepower_or_FWER(result.alt,response.probs.alt,test.type = input.info$test.type)
      disconj.pow = disconjunctivepowerfunc(result.alt)
      samplesize.ext.alt = c(samplesize.ext.alt, sum(Nfunc(result.alt)))
      # marg.pow = c(marg.pow)

      if (power.type == "Conjunctive"){
        extendgrid$pow[cutoffindex] = conj.pow
      }
      else{
        extendgrid$pow[cutoffindex] = disconj.pow
      }
      dnew=data.frame(rbind(startgrid.asy,extendgrid))
      dnew=dnew[!is.na(dnew$cutoff.eff),]
      GP.res = GP.optim(x = cbind(dnew$cutoff.eff, dnew$cutoff.fut),
                        y.t1E = dnew$tpIE, y.pow = dnew$pow,grid.length=grid.length,
                        errorrate = grid.inf$errorrate,
                        grid.min=c(grideff.min,gridfut.min), grid.max=c(grideff.max,gridfut.max),
                        Boundary.type = Boundary.type)
      nextcutoff = GP.res$next.cutoff

      extendgrid[cutoffindex + 1, 2] = nextcutoff[1]
      extendgrid[cutoffindex + 1, 3] = nextcutoff[2]

      prediction = data.frame(
        yhat.t1E = GP.res$prediction$yhat.t1E,
        yhat.pow = GP.res$prediction$yhat.pow,
        yhat.ESS.null = GP.res$prediction$yhat.ESS.null,
        yhat.ESS.alt = GP.res$prediction$yhat.ESS.alt,
        sd.t1E = GP.res$prediction$sd.t1E,
        sd.pow = GP.res$prediction$sd.pow,
        sd.ESS.null = GP.res$prediction$sd.ESS.null,
        sd.ESS.alt = GP.res$prediction$sd.ESS.alt,
        qup.ESS.null = GP.res$prediction$qup.ESS.null,
        qup.ESS.alt = GP.res$prediction$qup.ESS.alt,
        qdown.ESS.null = GP.res$prediction$qdown.ESS.null,
        qdown.ESS.alt = GP.res$prediction$qdown.ESS.alt,
        potentialcutoff = GP.res$prediction$potentialcutoff,
        qup.t1E = GP.res$prediction$qup.t1E,
        qdown.t1E = GP.res$prediction$qdown.t1E,
        qup.pow = GP.res$prediction$qup.pow,
        qdown.pow = GP.res$prediction$qdown.pow,
        xgrid = GP.res$prediction$xgrid
      )
      # Construct the ESS data frame for plot or analysis
      sample.mat=cbind(c(samplesize.start.null,samplesize.ext.null),c(samplesize.start.alt,samplesize.ext.alt))
      # Plot or not
      if (isTRUE(grid.inf$plotornot)) {
        prediction=GP.res$prediction
        colormap=colorRampPalette(rev(brewer.pal(11,'Spectral')))(32)
        target_line=grid.inf$errorrate
        xgrid.eff=prediction$xgrid[,1]
        xgrid.fut=prediction$xgrid[,2]
        nextcutoff.predict = nextcutoff
        colnames(nextcutoff.predict)=c("eff","fut","FWER")
        cleandata = dnew
        colnames(cleandata)=c("tpIE","eff","fut","pow")
        df=data.frame(FWER=prediction$yhat.t1E,eff=xgrid.eff,fut=xgrid.fut)
        Contour.tIE<-ggplot(df,aes(x=df$eff,y=df$fut,z=FWER))+
          scale_fill_gradientn(colors = colormap)+geom_tile(aes(fill=FWER))+
          geom_contour(breaks=c(target_line, seq(min(df$FWER),max(df$FWER),by=(max(df$FWER)-min(df$FWER))/10)),color="black")+
          geom_contour(breaks=target_line,color="white",linewidth=1.1)+
          labs(title="Mean type I error rate (FWER)", x="Cutoff value for efficacy",y="Cutoff value for futility")
        Contour.tIE=Contour.tIE+geom_point(data=cleandata,aes(x=cleandata$eff,y=cleandata$fut),color="black")+
          geom_point(data=nextcutoff.predict,aes(x=nextcutoff.predict$eff,y=nextcutoff.predict$fut),color="pink")

        # Extract the contour data
        contour_data_tIE <- ggplot_build(Contour.tIE)$data[[2]]
        # Record the contour that has FWER equal to the target
        contour_data_tIE_subset <- contour_data_tIE[contour_data_tIE$level == target_line, ]
        # Order and split the data to ensure the plot is drawn correctly
        contour_data_tIE_subset=contour_data_tIE_subset[order(contour_data_tIE_subset$piece,contour_data_tIE_subset$x),]
        contour_data_tIE_subset_1=contour_data_tIE_subset[contour_data_tIE_subset$piece==1,]
        contour_data_tIE_subset_2=contour_data_tIE_subset[contour_data_tIE_subset$piece==2,]
        # To make sure the data frame is not empty
        if (nrow(contour_data_tIE_subset_1) == 0){
          contour_data_tIE_subset_1[1,]=(rep(NA,dim(contour_data_tIE_subset_1)[[2]]))
        } else if (nrow(contour_data_tIE_subset_2) == 0){
          contour_data_tIE_subset_2[1,]=(rep(NA,dim(contour_data_tIE_subset_2)[[2]]))
        }


        df=data.frame(precision=prediction$sd.t1E,eff=xgrid.eff,fut=xgrid.fut)
        Contour.sd<-ggplot(df,aes(x=df$eff,y=df$fut,z=df$precision))+
          scale_fill_gradientn(colors = colormap)+geom_tile(aes(fill=df$precision))+
          geom_contour(breaks=seq(min(df$precision),max(df$precision),by=(max(df$precision)-min(df$precision))/10),color="black")+labs(title="sd of each contour plot", x="Cutoff value for efficacy",y="Cutoff value for futility")
        Contour.sd=Contour.sd+
          geom_path(data = contour_data_tIE_subset_1,
                    aes(contour_data_tIE_subset_2$x,contour_data_tIE_subset_2$y,z=NA),color="white",linewidth=1.1)+
          geom_path(data = contour_data_tIE_subset_2, aes(contour_data_tIE_subset_2$x,contour_data_tIE_subset_2$y,z=NA),color="white",linewidth=1.1)+
          geom_point(data=cleandata,aes(cleandata$eff,cleandata$fut,z=NA),color="black")+
          geom_point(data=nextcutoff.predict,aes(x=nextcutoff.predict$eff,y=nextcutoff.predict$fut,z=NA),color="pink")

        df=data.frame(Power=prediction$yhat.pow,eff=xgrid.eff,fut=xgrid.fut)
        Contour.pow<-ggplot(df,aes(x=df$eff,y=df$fut,z=df$Power))+
          scale_fill_gradientn(colors = colormap)+geom_tile(aes(fill=df$Power))+
          geom_contour(breaks=seq(min(df$Power),max(df$Power),by=(max(df$Power)-min(df$Power))/10),color="black")+labs(title="Mean power", x="Cutoff value for efficacy",y="Cutoff value for futility")
        Contour.pow=Contour.pow+
          geom_path(data = contour_data_tIE_subset_1,
                    aes(contour_data_tIE_subset_2$x,contour_data_tIE_subset_2$y,z=NA),color="white",linewidth=1.1)+
          geom_path(data = contour_data_tIE_subset_2, aes(contour_data_tIE_subset_2$x,contour_data_tIE_subset_2$y,z=NA),color="white",linewidth=1.1)+
          geom_point(data=cleandata,aes(cleandata$eff,cleandata$fut,z=NA),color="black")+
          geom_point(data=nextcutoff.predict,aes(x=nextcutoff.predict$eff,y=nextcutoff.predict$fut,z=NA),color="pink")

        df=data.frame(NullESS=prediction$yhat.ESS.null,eff=xgrid.eff,fut=xgrid.fut)
        Contour.nullESS<-ggplot(df,aes(x=df$eff,y=df$fut,z=df$NullESS))+
          scale_fill_gradientn(colors = colormap)+geom_tile(aes(fill=df$NullESS))+
          geom_contour(breaks=seq(min(df$NullESS),max(df$NullESS),by=(max(df$NullESS)-min(df$NullESS))/10),color="black")+labs(title="Mean ESS under null", x="Cutoff value for efficacy",y="Cutoff value for futility")
        Contour.nullESS=Contour.nullESS+
          geom_path(data = contour_data_tIE_subset_1,
                    aes(contour_data_tIE_subset_2$x,contour_data_tIE_subset_2$y,z=NA),color="white",linewidth=1.1)+
          geom_path(data = contour_data_tIE_subset_2, aes(contour_data_tIE_subset_2$x,contour_data_tIE_subset_2$y,z=NA),color="white",linewidth=1.1)+
          geom_point(data=cleandata,aes(cleandata$eff,cleandata$fut,z=NA),color="black")+
          geom_point(data=nextcutoff.predict,aes(x=nextcutoff.predict$eff,y=nextcutoff.predict$fut,z=NA),color="pink")

        df=data.frame(AltESS=prediction$yhat.ESS.alt,eff=xgrid.eff,fut=xgrid.fut)
        Contour.altESS<-ggplot(df,aes(x=df$eff,y=df$fut,z=df$AltESS))+
          scale_fill_gradientn(colors = colormap)+geom_tile(aes(fill=df$AltESS))+
          geom_contour(breaks=seq(min(df$AltESS),max(df$AltESS),by=(max(df$AltESS)-min(df$AltESS))/10),color="black")+labs(title="Mean ESS under alternative", x="Cutoff value for efficacy",y="Cutoff value for futility")
        Contour.altESS=Contour.altESS+
          geom_path(data = contour_data_tIE_subset_1,
                    aes(contour_data_tIE_subset_2$x,contour_data_tIE_subset_2$y,z=NA),color="white",linewidth=1.1)+
          geom_path(data = contour_data_tIE_subset_2, aes(contour_data_tIE_subset_2$x,contour_data_tIE_subset_2$y,z=NA),color="white",linewidth=1.1)+
          geom_point(data=cleandata,aes(cleandata$eff,cleandata$fut,z=NA),color="black")+
          geom_point(data=nextcutoff.predict,aes(x=nextcutoff.predict$eff,y=nextcutoff.predict$fut,z=NA),color="pink")
        # Plot these figures
        ggarrange(Contour.tIE,Contour.pow,Contour.nullESS,Contour.altESS,Contour.sd,ncol = 2,nrow=3)
      }
      # Need to add a early stop criteria discuss with Dave on 05/09/2023
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
