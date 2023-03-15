#' @title modelinf.fun
#' @description This function summarize the input parameters describing the model for analysis and transfer them into a list
#' @param model The statistical model. ibb: betabinomial model / tlr: logistic model
#' @param ibb.inf The list of information for betabinomial model including:
#'     betabinomialmodel: The betabinomial model, pi.star: prior response rate,
#'     pess: prior effective sample size
#' @param tlr.inf The list of information for logistic model including:
#'     The mean (mu), variance (sigma), degree of freedom (df) of the intercept
#'     and the main effect of the linear terms in logistic model.
#'     reg.inf: The type of linear function in logistic model.
#'     variable.inf: Fixeffect/Mixeffect. Indicating whether a mix effect model
#'     is used (for time trend effect modelling)
#'
#' @return A list of model information including model, ibb.inf and tlr.inf
#' @export
#'
#' @examples
#' modelinf.fun(model = "tlr",
#'   ibb.inf = list(pi.star = 0.5,
#'                  pess = 2,
#'                  betabinomialmodel = ibetabinomial.post),
#'   tlr.inf = list(beta0_prior_mu = 0,
#'                  beta1_prior_mu = 0,
#'                  beta0_prior_sigma = 2.5,
#'                  beta1_prior_sigma = 2.5,
#'                  beta0_df = 7,
#'                  beta1_df = 7,
#'                  reg.inf =  "main",
#'                  variable.inf = "Fixeffect"
#'                  ))
#' @author Ziyan Wang
modelinf.fun <- function(model = "tlr",
                      ibb.inf = list(pi.star = 0.5,
                                     pess = 2,
                                     betabinomialmodel = ibetabinomial.post),
                      tlr.inf = list(
                        beta0_prior_mu = 0,
                        beta1_prior_mu = 0,
                        beta0_prior_sigma = 2.5,
                        beta1_prior_sigma = 2.5,
                        beta0_df = 7,
                        beta1_df = 7,
                        reg.inf =  "main",
                        variable.inf = "Fixeffect"
                      )) {
  return(list(model,
              ibb.inf = ibb.inf,
              tlr.inf = tlr.inf))
}
