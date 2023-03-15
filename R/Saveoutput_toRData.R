#' @title Save.resulttoRDatafile
#' @description This function generates the name of file for output table and dataset
#'
#' @param input.info A list of input information require for trial simulation
#'
#' @return A list of name for table and dataset
#' @export
#'
#' @examples
#' Save.resulttoRDatafile(
#' input.info = list(
#'   response.probs = c(0.4, 0.4),
#'   ns = c(30, 60, 90, 120, 150),
#'   max.ar = 0.75,
#'   rand.type = "Urn",
#'   max.deviation = 3,
#'   model.inf = list(
#'     model = "tlr",
#'     ibb.inf = list(
#'       pi.star = 0.5,
#'       pess = 2,
#'       betabinomialmodel = ibetabinomial.post
#'     ),
#'     tlr.inf = list(
#'       beta0_prior_mu = 0,
#'       beta1_prior_mu = 0,
#'       beta0_prior_sigma = 2.5,
#'       beta1_prior_sigma = 2.5,
#'       beta0_df = 7,
#'       beta1_df = 7,
#'       reg.inf =  "main",
#'       variable.inf = "Fixeffect"
#'     )
#'   ),
#'   Stop.type = "Early-Pocock",
#'   Boundary.type = "Symmetric",
#'   Random.inf = list(
#'     Fixratio = FALSE,
#'     Fixratiocontrol = NA,
#'     BARmethod = "Thall",
#'     Thall.tuning.inf = list(tuningparameter = "Fixed",  fixvalue = 1)
#'   ),
#'   trend.inf = list(
#'     trend.type = "step",
#'     trend.effect = c(0, 0),
#'     trend_add_or_multip = "mult"
#'   )
#' ))
#' @author Ziyan Wang
Save.resulttoRDatafile = function(input.info = list(
  response.probs = c(0.4, 0.4),
  ns = c(30, 60, 90, 120, 150),
  max.ar = 0.75,
  rand.type = "Urn",
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
)) {
  if (sum(input.info$trend.inf$trend.effect != 0) > 0) {
    trendornot = "TREND"
  }
  else{
    trendornot = "NOTREND"
  }
  nameTable = paste0(
    "TABLE",
    trendornot,
    toupper(input.info$Stopbound.inf$Stop.type),
    toupper(input.info$Stopbound.inf$Boundary.type),
    toupper(input.info$Random.inf$BARmethod),
    ".RData"
  )
  nameData = paste0(
    "DATA",
    trendornot,
    toupper(input.info$Stopbound.inf$Stop.type),
    toupper(input.info$Stopbound.inf$Boundary.type),
    toupper(input.info$Random.inf$BARmethod),
    ".RData"
  )
  return(list(nameTable = nameTable, nameData = nameData))
}
