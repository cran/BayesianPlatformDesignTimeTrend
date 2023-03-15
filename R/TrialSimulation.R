#' @title Trial simulation
#' @description This function simulates and does final analysis of a trial with one scenario.
#'     The time cost of this function depend on the cpu cores of the user's cpmputer.
#' @param ntrials A numeric variable indicating how many trial replicates you want to run. Default is 5000.
#' @param trial.fun The function of trial simulation for more see \code{\link{simulatetrial}}
#' @param input.info A list of input information including all information required for trial simulation.
#' @param cl A numeric variable indicating how many cores you want to use in parallel programming.
#'
#' @return A list of output including the final output of each trial replicates called 'result'
#'     The analysis result table of the specific trial called 'OPC' and the file name for saving these output on the computer
#' @importFrom foreach foreach
#' @importFrom matrixStats colVars
#' @import ggplot2
#' @import rstantools
#' @export
#'
#' @examples
#' set.seed(1)
#' \donttest{Trial.simulation(ntrials = 2, cl = 2)}
#' @author Ziyan Wang
Trial.simulation = function(ntrials = 5000,
                            trial.fun = simulatetrial,
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
                              Stopbound.inf = Stopboundinf(
                                Stop.type = "Early-Pocock",
                                Boundary.type = "Symmetric",
                                cutoff = c(0.99, 0.01)
                              ),
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
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores(logical = FALSE))

  registerDoParallel(cores = cl)

  message("Start trial information initialisation")
  #-Initialising evaluation metrics-
  bias = {

  }

  TIE_POWER = {

  }

  rMSE = {

  }

  totalsample = {

  }

  totaleachN = {

  }

  totaleachS = {

  }

  result = {

  }
  #-----------------------------------------
  result = foreach(icount(ntrials)) %dopar% trial.fun(
    response.probs = input.info$response.probs,
    ns = input.info$ns,
    max.ar = input.info$max.ar,
    rand.algo = input.info$rand.algo,
    max.deviation = input.info$max.deviation,
    model.inf = input.info$model.inf,
    Stopbound.inf = input.info$Stopbound.inf,
    Random.inf = input.info$Random.inf,
    trend.inf = input.info$trend.inf
  )
  #Save output data
  FWER = FWER_disconjunctivepowerfunc(result)
  TIE_POWER = rbind(TIE_POWER, FWER)
  meanres = Meanfunc(result) - (logit(input.info$response.probs[-1]) - logit(input.info$response.probs[1]))
  # varres=Varfunc(result)
  # bias=rbind(bias,paste0(meanres," (",varres,")"))
  bias = rbind(bias, meanres)
  var = varfunc(result)
  rmse = sqrt(meanres ^ 2 + var)
  rMSE = rbind(rMSE, rmse)

  Nperarm = Nfunc(result)
  totaleachN = rbind(totaleachN, Nperarm)
  Sperarm = Sperarmfunc(result)
  totaleachS = rbind(totaleachS, Sperarm)
  totalN = sum(Nperarm)
  totalsample = rbind(totalsample, totalN)
  OPC = data.frame(
    "Type I Error or Power" = TIE_POWER,
    Bias = bias,
    rMSE = rMSE,
    "N per arm" = totaleachN,
    "Survive per arm" = totaleachS,
    N = totalsample
  )

  Nameofsaveddata = Save.resulttoRDatafile(input.info)
  result = list(result)
  names(result) = paste0(
    paste0(paste(
      as.character(
        stringr::str_replace_all(input.info$response.probs, "[^[:alnum:]]", "")
      )
    ), collapse = ""),
    "TimeTrend",
    paste0(paste(
      as.character(
        stringr::str_replace_all(input.info$trend.inf$trend.effect, "[^[:alnum:]]", "")
      ), collapse = ""
    )),
    "stage",
    length(input.info$ns),
    input.info$model.inf$tlr.inf$reg.inf
  )
  rownames(OPC) =
    paste0(
      paste0(paste(
        as.character(
          stringr::str_replace_all(input.info$response.probs, "[^[:alnum:]]", "")
        )
      ), collapse = ""),
      "TimeTrend",
      paste0(paste(
        as.character(
          stringr::str_replace_all(input.info$trend.inf$trend.effect, "[^[:alnum:]]", "")
        ), collapse = ""
      )),
      "stage",
      length(input.info$ns),
      input.info$model.inf$tlr.inf$reg.inf
    )
  doParallel::stopImplicitCluster()
  return(list(
    result = result,
    OPC = OPC,
    Nameofsaveddata = Nameofsaveddata
  ))
}
