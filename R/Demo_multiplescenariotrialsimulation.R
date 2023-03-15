#' @title demo_multscenario
#' @description This is a demo function simulating multi-arm multi-stage design with two different null scenarios where the response probability of control is 0.15 and 0.4, respectively.
#'    The clinically meaningful increment on probability scale is 0.2. The stopping boundary is the OBF. The cutoff vector in the demo is tuned to keep Type I error rate to be 0.05. The output data can be saved as .RData file
#'
#' @param ntrials A numeric value. The number of total trail replicates for each scenario.
#' @param cl A numeric variable indicating how many cores you want to use in parallel programming.
#' @param save_data An indicator of whether the output data need to be saved. Default is FALSE.
#'
#' @return A list of data for plotting. One is results of trial replicates for all scenarios. The other one is a data frame containing all summarised evaluation metrics for all scenarios
#' @export
#'
#' @examples
#' \donttest{demo_multscenario(ntrials = 2, cl = 2, save_data = FALSE)}
#' @author Ziyan Wang
demo_multscenario = function(ntrials = 1000,
                             cl = 2,
                             save_data = FALSE) {
  message("Start trial information initialisation")
  # ns = list(seq(60, 300, 60), seq(30, 300, 30))
  ns = list(seq(60, 300, 60))
  null.response.probs1 = 0.15
  alt.response.probs1 = 0.35
  null.response.probs2 = 0.4
  alt.response.probs2 = 0.6
  scenario = matrix(
    c(
      null.response.probs1,
      null.response.probs1,
      null.response.probs1,
      alt.response.probs1,
      null.response.probs2,
      null.response.probs2,
      null.response.probs2,
      alt.response.probs2
    ),
    ncol = 2,
    nrow = 4,
    byrow = T
  )

  cutoffearlyOBF = c(4.391, 4.661, 4.281, 4.512)

  result = {

  }
  OPC = {

  }
  cutoffindex = 1
  message(
    "Start trial simulation. This is a two arm trial simulation. There are two null scenarios and two alternative scenarios and for each scenario there are two vectors of number of patients at each stage in this demo. There are 8 rounds."
  )
  for (i in 1:dim(scenario)[1]) {
    for (z in 1:length(ns)) {
      message(paste(
        "Scenario",
        i,
        "with patient number sequence",
        ns[z],
        "under simulation"
      ))
      restlr = Trial.simulation(
        ntrials = ntrials,
        trial.fun = simulatetrial,
        input.info = list(
          response.probs = scenario[i,],
          ns = ns[[z]],
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
            Stop.type = "Early-OBF",
            Boundary.type = "Symmetric",
            cutoff = c(cutoffearlyOBF[cutoffindex], cutoffearlyOBF[cutoffindex])
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
        cl = cl
      )
      cutoffindex = cutoffindex + 1
      result = c(result, restlr$result)
      OPC = rbind(OPC, restlr$OPC)
      message(paste("Finished round", cutoffindex))
    }
  }
  if (isTRUE(save_data)) {
    save(result, file = restlr$Nameofsaveddata$nameData)
    save(OPC, file = restlr$Nameofsaveddata$nameTable)
  }
  return(list(result = result, OPC = OPC))
}
