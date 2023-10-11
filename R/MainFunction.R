#' @title simulatetrial
#' @description This function simulates a MAMS trial applying adaptive methods where the time trend effect can be studied.
#' @param ii Meaning less parameter but required for foreach function in doParallel package
#' @param response.probs A vector of true response probability for each arm. Default response.probs = c(0.4, 0.4).
#' @param test.type A indicator of whether to use one side test or two side test for each treatment-control comparison.
#' @param ns  A vector of accumulated number of patient at each stage. Default is ns = c(30, 60, 90, 120, 150).
#' @param max.ar The upper boundary for randomisation ratio for each arm. Default is 0.75 for a two arm trial. The minimum value depends on K where 1 - max.ar <= 1/K
#' @param rand.algo The method of applying patient allocation with a given randomisation probability vector. Default is "Urn".
#' @param max.deviation The tuning parameter for Urn randomisation method. Default is 3.
#' @param Stopbound.inf The list of stop boundary information for more see \code{\link{Stopboundinf}}
#' @param Random.inf The list of Adaptive randomisation information for more see \code{\link{Randomisation.inf}}
#' @param trend.inf The list of time trend information
#' @param model.inf The list of interim data analysis model information for more see \code{\link{modelinf.fun}}
#'
#' @return A matrix including all evaluation metrics
#' @export
#'
#' @examples
#' set.seed(1)
#' simulatetrial(response.probs = c(0.4, 0.4),
#'                ns = c(30, 60, 90, 120, 150),
#'                max.ar = 0.75,
#'                test.type = "Twoside",
#'                rand.algo = "Urn",
#'                max.deviation = 3,
#'                model.inf = list(
#'                model = "tlr",
#'                ibb.inf = list(
#'                pi.star = 0.5,
#'                pess = 2,
#'                betabinomialmodel = ibetabinomial.post
#'                ),
#'                tlr.inf = list(
#'                beta0_prior_mu = 0,
#'                beta1_prior_mu = 0,
#'                beta0_prior_sigma = 2.5,
#'                beta1_prior_sigma = 2.5,
#'                beta0_df = 7,
#'                beta1_df = 7,
#'                reg.inf =  "main",
#'                variable.inf = "Fixeffect"
#'                )
#'                ),
#'                Stopbound.inf = Stopboundinf(
#'                Stop.type = "Early-Pocock",
#'                Boundary.type = "Symmetric",
#'                cutoff = c(0.99,0.01)
#'                ),
#'                Random.inf = list(
#'                Fixratio = FALSE,
#'                Fixratiocontrol = NA,
#'                BARmethod = "Thall",
#'                Thall.tuning.inf = list(tuningparameter = "Fixed",  fixvalue = 1),
#'                Trippa.tuning.inf = list(a = 10, b = 0.75)
#'                ),
#'                trend.inf = list(
#'                trend.type = "step",
#'                trend.effect = c(0, 0),
#'                trend_add_or_multip = "mult"
#'                ))
#' @author Ziyan Wang
simulatetrial <- function(ii,
                           response.probs = c(0.4, 0.4),
                           ns = c(30, 60, 90, 120, 150),
                           test.type = "Twoside",
                           max.ar = 0.75,
                           rand.algo = "Urn",
                           max.deviation = 3,
                           model.inf  = list(
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
                           Stopbound.inf = Stopbound.inf,
                           Random.inf = Random.inf,
                           trend.inf = trend.inf) {
  #-Boundary construction-
  Boundary = Boundaryconstruction(Stopbound.inf, ns = ns)
  cutoffeff = Boundary$Efficacy.boundary
  cutoffful = Boundary$Fultility.boundary

  #-Randomisation inf check-
  Random.inf = Randomisation.inf(Random.inf)
  Fixratio = Random.inf$Fixratio
  Fixratiocontrol = Random.inf$Fixratiocontrol
  if (!is.logical(Fixratio)) stop("Error: Fixratio should be a logical value (TRUE/FALSE)")
  if (isTRUE(Fixratio) & !is.numeric(Fixratiocontrol)) stop("Error: Fixratiocontrol argument should be numeric for fix ratio approach")
  BARmethod = Random.inf$BARmethod
  #List of information required for Thall's approach
  Thall.tuning.inf = Random.inf$Thall.tuning.inf
  #List of information required for Trippa's approach
  Trippa.tuning.inf = Random.inf$Trippa.tuning.inf
  #Identify whether the tuning parameter for Thall's approach is fixed or not
  tuningparameter = Thall.tuning.inf$tuningparameter
  #Fixed tuning parameter for Thall's approach if existing
  c = Thall.tuning.inf$c
  #Fixed tuning parameter for Trippa's approach if existing
  a = Trippa.tuning.inf$a
  b = Trippa.tuning.inf$b
  #-Simulation setting-
  #Initialize Data
  initialised.par = Initializetrialparameter(response.probs, ns)
  #Extract the list of initialized parameter
  K = initialised.par$K
  armleft = initialised.par$armleft
  treatmentindex = initialised.par$treatmentindex
  n = initialised.par$n
  y1 = initialised.par$y1
  groupwise.response.probs = initialised.par$groupwise.response.probs
  randomprob = initialised.par$randomprob
  z = initialised.par$z
  y = initialised.par$y
  group_indicator = initialised.par$group_indicator
  post.prob.best.mat = initialised.par$post.prob.best.mat

  if (test.type == "Twoside" & isFALSE(model.inf$Random.inf$Fixratio)){
    warning("Adaptive randomisation aims to allocate to superior arms,
            while two side test conclude both superiority and inferiority.
            Therefore, these two approaches have conflicting purposes.")
  }

  #-max.ar check
  if (1 - max.ar > 1/K){
    stop("Error: The lower allocation ratio should be at least 1/K. Please check the number of arm at the beginning and the max.ar")
  }

  #Initialize the output data frame: stats
  stats = OutputStats.initialising(model.inf$tlr.inf$variable.inf,
                                   model.inf$tlr.inf$reg.inf,
                                   ns,
                                   K)
  # ----For random effect initialisation
  ntemp=matrix(rep(NA,nrow(stats)*K),nrow=K)
  ytemp=matrix(rep(NA,nrow(stats)*K),nrow=K)
  # -----
  #Generating time trend function ("trend.function") based on input trend information
  #and check whether the input is reasonable using the "Timeindicator" variable
  Timetrendfunctionlist = Timetrend.fun(trend.inf)
  trend.function = Timetrendfunctionlist$trend.function
  Timeindicator = Timetrendfunctionlist$timetrendornot
  trend_add_or_multip = Timetrendfunctionlist$trend_add_or_multip
  trend.effect = Timetrendfunctionlist$trend.effect
  #-Start data generation and analysing stage by stage-
  for (group in 1:nrow(stats)) {
    # Number of new enrolls during current group
    n.new = c(0, ns)[group + 1] - c(0, ns)[group]

    # Patient randomisation based on given randomization ratio
    # Parameter checking
    if (rand.algo %in% c("Coin", "Urn")) {
      random.output = AdaptiveRandomisation(
        Fixratio,
        rand.algo,
        K,
        n.new,
        randomprob,
        treatmentindex,
        groupwise.response.probs,
        group,
        armleft,
        max.deviation,
        trend_add_or_multip,
        trend.function,
        trend.effect,
        ns,
        Fixratiocontrol
      )
    }
    else {
      stop("Error: randomisation type wrong")
    }

    # Extract the output list
    nstage = random.output$nstage
    ystage = random.output$ystage
    znew = random.output$znew
    ynew = random.output$ynew

    # Update parameters n: accumulated patients to each arm; y1: accumulated outcomes
    # z: treatment index vector; y: outcome index vector; N: current total number of patients
    # group_indicator: stage index vector indicating the stage at which each patient is treated
    n = n + nstage
    y1 = y1 + ystage
    ntemp[,group]=nstage
    ytemp[,group]=ystage
    stats2 = as.vector(matrix(c(n, y1), nrow = 2, byrow = TRUE))
    z = c(z, znew)
    y = c(y, ynew)
    N = length(z)
    group_indicator = c(group_indicator, rep(group, n.new))

    #-Modelling-
    if (model.inf$model == "ibb") {
      # Beta-binomial model can not estimate mean treatment effect (stats4);
      #    estimate treatment effect variance (stats5);
      #    can not deal with time trend effect (stats6);
      #    and time trend treatment effect interaction (stats7);
      #    Therefore stats4,5,6,7 are set to be blank to generate output matrix.
      stats4 = {

      }
      stats5 = {

      }
      stats6 = {

      }
      stats7 = {

      }
      if (BARmethod == "Thall") {
        #Debugged for K arm by Ziyan Wang on 18:21 11/08/2022. Used to be stoperror for K arm)
        zdropped = z[z %in% as.numeric(names(randomprob[, randomprob != 0]))]
        ydropped = y[z %in% as.numeric(names(randomprob[, randomprob != 0]))]
        Ndropped = length(zdropped)
        group_indicator_dropped = group_indicator[z %in% as.numeric(names(randomprob[, randomprob != 0]))]
        zlevel = as.numeric(levels(factor(zdropped)))

        for (zindex in 1:armleft) {
          zdropped[zdropped == levels(factor(zdropped))[zindex]] = zindex
        }

        data = list(
          K = armleft,
          N = Ndropped,
          y = array(ydropped, dim = Ndropped),
          z = array(zdropped, dim = Ndropped),
          group = group,
          pistar = model.inf$ibb.inf$pi.star,
          pess = model.inf$ibb.inf$pess
        )

        fit <- rstan::sampling(
          stanmodels$betabinom,
          data = data,
          chains = 1,
          refresh = 0,
          warmup = 2500,
          iter = 5000
        )
        sampeff = rstan::extract(fit, 'theta')[[1]]
        colnames(sampeff) <- c(0, treatmentindex)
        control = matrix(sampeff[, 1])
        colnames(control) <- 0
        treatment = matrix(sampeff[, -1], ncol = dim(sampeff)[2] - 1)
        colnames(treatment) <- treatmentindex
        #posterior of each treatment better than control
        post.prob.btcontrol = colMeans(sapply(treatmentindex, function(treatmentindex) {
          postprob = treatment[, as.numeric(colnames(treatment)) == treatmentindex] >
            control
          return(postprob)
        }))

        stats1 = rep(NA, K - 1)
        names(stats1) = seq(1, K - 1)
        stats1[treatmentindex] = post.prob.btcontrol

        post.prob.best = colMeans(rstan::extract(fit, 'times_to_be_best')[[1]])

        for (q in 1:armleft) {
          post.prob.best.mat[group, zlevel[q]] = post.prob.best[q]
        }

        post.prob.best = post.prob.best.mat[group,]
        #Normalizing in case any value equals zero
        post.prob.best = post.prob.best + 1e-7
        post.prob.best = post.prob.best / sum(post.prob.best)
        # Drop both superior and inferior arm and make hypothesis testing
        test_drop.inf = testing_and_armdropping(
          K = K,
          armleft = armleft,
          post.prob.btcontrol = post.prob.btcontrol,
          group = group,
          cutoffeff = cutoffeff,
          cutoffful = cutoffful,
          treatmentindex = treatmentindex,
          test.type = test.type
        )
        stats3 = test_drop.inf$stats3
        armleft = test_drop.inf$armleft
        treatmentindex = test_drop.inf$treatmentindex
      }
      else {
        nibb = n[c(1, treatmentindex + 1)]
        y1ibb = y1[c(1, treatmentindex + 1)]
        # Analytic result of posterior probability
        # More than one posterior probability for more than one treatment arm situation
        # (Much faster than stan,
        # however it can not generate posterior probability of each arm to be the best).
        resultibb = model.inf$ibb.inf$ibetabinomial.post(
          n = nibb,
          y = y1ibb,
          pi.star = model.inf$ibb.inf$pi.star,
          pess = model.inf$ibb.inf$pess
        )
        post.prob.btcontrol = resultibb
        stats1 = rep(NA, K - 1)
        names(stats1) = seq(1, K - 1)
        stats1[treatmentindex] = post.prob.btcontrol
        # Drop both superior and inferior arm and make hypothesis testing
        test_drop.inf = testing_and_armdropping(
          K = K,
          armleft = armleft,
          post.prob.btcontrol = post.prob.btcontrol,
          group = group,
          cutoffeff = cutoffeff,
          cutoffful = cutoffful,
          treatmentindex = treatmentindex,
          test.type = test.type
        )
        stats3 = test_drop.inf$stats3
        armleft = test_drop.inf$armleft
        treatmentindex = test_drop.inf$treatmentindex
      }
      if (isFALSE(Fixratio)) {
        #-Adjust the posterior randomisation ratio-
        randomprob = ARmethod(
          BARmethod,
          group,
          stats,
          post.prob.btcontrol,
          K,
          n,
          tuningparameter,
          c,
          a,
          b,
          post.prob.best,
          max.ar,
          armleft,
          treatmentindex
        )
      }
    }
    else if (model.inf$model == "tlr") {
      stan.data.temp = stan.logisticmodeltrans(
        z,
        y,
        randomprob,
        group_indicator,
        armleft,
        group,
        model.inf$tlr.inf$variable.inf,
        model.inf$tlr.inf$reg.inf
      )
      zdropped = stan.data.temp$zdropped
      ydropped = stan.data.temp$ydropped
      Ndropped = stan.data.temp$Ndropped
      group_indicator_dropped = stan.data.temp$group_indicator_dropped
      zlevel = stan.data.temp$zlevel
      xdummy = stan.data.temp$xdummy
      if (model.inf$tlr.inf$variable.inf == "Fixeffect" |
          group == 1) {
        data = list(
          K = dim(xdummy)[2],
          N = Ndropped,
          y = array(ydropped, dim = Ndropped),
          z = array(zdropped, dim = Ndropped),
          x = xdummy,
          group = group_indicator_dropped,
          beta0_prior_mu = model.inf$tlr.inf$beta0_prior_mu,
          beta1_prior_mu = model.inf$tlr.inf$beta1_prior_mu,
          beta0_prior_sigma = model.inf$tlr.inf$beta0_prior_sigma,
          beta1_prior_sigma = model.inf$tlr.inf$beta1_prior_sigma,
          beta0_nu = model.inf$tlr.inf$beta0_df,
          beta1_nu = model.inf$tlr.inf$beta1_df
        )
        fit <- rstan::sampling(
          stanmodels$logisticdummy,
          data = data,
          chains = 1,
          refresh = 0,
          warmup = 2500,
          iter = 5000
        )
        beta0 = matrix(rstan::extract(fit, 'b_Intercept')[[1]], ncol = 1)
        statsbeta0 = colMeans(beta0)
        processedfitresult = resultstantoRfunc(
          group = group,
          reg.inf = model.inf$tlr.inf$reg.inf,
          variable.inf = model.inf$tlr.inf$variable.inf,
          fit = fit,
          armleft = armleft,
          treatmentindex = treatmentindex,
          K = K,
          ns = ns
        )
        stats4 = processedfitresult$stats4
        stats5 = processedfitresult$stats5
        stats6 = processedfitresult$stats6
        stats7 = processedfitresult$stats7
        stats1 = processedfitresult$stats1
        sampefftotal = processedfitresult$sampefftotal
        post.prob.btcontrol = processedfitresult$post.prob.btcontrol

        #-Calculating posterior probability of each arm (including control) to be the best arm-
        # post.prob.best: The posterior probability of each arm (including control) to be the best arm
        # This is required for Thall's randomisation approach
        for (q in 1:armleft) {
          post.prob.best.mat[group, zlevel[q]] = (sum(max.col(sampefftotal) == q)) /
            2500
        }
        post.prob.best = post.prob.best.mat[group,]
        #Normalizing in case any value equals zero
        post.prob.best = post.prob.best + 1e-7
        post.prob.best = post.prob.best / sum(post.prob.best)

        # Drop both superior and inferior arm and make hypothesis testing
        test_drop.inf = testing_and_armdropping(
          K = K,
          armleft = armleft,
          post.prob.btcontrol = post.prob.btcontrol,
          group = group,
          cutoffeff = cutoffeff,
          cutoffful = cutoffful,
          treatmentindex = treatmentindex,
          test.type = test.type
        )
        stats3 = test_drop.inf$stats3
        armleft = test_drop.inf$armleft
        treatmentindex = test_drop.inf$treatmentindex
      }
      else if (model.inf$tlr.inf$variable.inf == "Mixeffect.stan") {
        dataran = list(
          K = armleft,
          N = Ndropped,
          Y = array(ydropped, dim = Ndropped),
          z = array(zdropped, dim = Ndropped),
          X = xdummy,
          groupmax = max(group_indicator_dropped),
          group = group_indicator_dropped,
          beta0_prior_mu = model.inf$tlr.inf$beta0_prior_mu,
          beta1_prior_mu = model.inf$tlr.inf$beta1_prior_mu,
          beta0_prior_sigma = model.inf$tlr.inf$beta0_prior_sigma,
          beta1_prior_sigma = model.inf$tlr.inf$beta1_prior_sigma,
          beta0_nu = model.inf$tlr.inf$beta0_df,
          beta1_nu = model.inf$tlr.inf$beta1_df
        )

        fit <- rstan::sampling(
          stanmodels$randomeffect,
          data = dataran,
          chains = 1,
          refresh = 0,
          warmup = 2500,
          iter = 5000
        )

        beta0 = matrix(rstan::extract(fit, 'b_Intercept')[[1]], ncol = 1)
        beta1 = rstan::extract(fit, "beta")[[1]]
        statsbeta0 = mean(beta0+beta1[,1])

        processedfitresult.rand = resultstantoRfunc.rand(
          group = group,
          fit = fit,
          armleft = armleft,
          treatmentindex = treatmentindex,
          K = K,
          ns = ns
        )
        stats4 = processedfitresult.rand$stats4
        stats5 = processedfitresult.rand$stats5
        stats6 = processedfitresult.rand$stats6
        names(stats6) = max(group_indicator_dropped) - as.numeric(names(stats6)) +
          2
        stats6 = stats6[order(as.numeric(names(stats6)))]
        stats7 = processedfitresult.rand$stats7
        stats1 = processedfitresult.rand$stats1
        sampefftotal = processedfitresult.rand$sampefftotal
        post.prob.btcontrol = processedfitresult.rand$post.prob.btcontrol

        #-Calculating posterior probability of each arm (including control) to be the best arm-
        # post.prob.best: The posterior probability of each arm (including control) to be the best arm
        # This is required for Thall's randomisation approach
        for (q in 1:armleft) {
          post.prob.best.mat[group, zlevel[q]] = (sum(max.col(sampefftotal) == q)) /
            2500
        }
        post.prob.best = post.prob.best.mat[group,]
        #Normalizing in case any value equals zero
        post.prob.best = post.prob.best + 1e-7
        post.prob.best = post.prob.best / sum(post.prob.best)
        # Drop both superior and inferior arm and make hypothesis testing
        test_drop.inf = testing_and_armdropping(
          K = K,
          armleft = armleft,
          post.prob.btcontrol = post.prob.btcontrol,
          group = group,
          cutoffeff = cutoffeff,
          cutoffful = cutoffful,
          treatmentindex = treatmentindex,
          test.type = test.type
        )
        stats3 = test_drop.inf$stats3
        armleft = test_drop.inf$armleft
        treatmentindex = test_drop.inf$treatmentindex
      }
      if (isFALSE(Fixratio)) {
        #-Adjust the posterior randomisation ratio-
        randomprob = ARmethod(
          BARmethod,
          group,
          stats,
          post.prob.btcontrol,
          K,
          n,
          tuningparameter,
          c,
          a,
          b,
          post.prob.best,
          max.ar,
          armleft,
          treatmentindex
        )
      }
    }
    stats[group,] = c(stats1, stats2, stats3, round(statsbeta0, 3), stats4, stats5, stats6, stats7)
    if (armleft == 1) {
      break
    }
  }
  return(stats)
}
