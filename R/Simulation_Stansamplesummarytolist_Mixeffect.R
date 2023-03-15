#' @title resultstantoRfunc.rand
#' @description This function summarise the mix effect stan output data to and transform them to be readable.

#' @param group The current stage
#' @param fit The stan output
#' @param armleft The number of treatment left in the platform (>2)
#' @param treatmentindex A vector of treatment index at the beginning of a trial
#' @param K Total number of arms at the beginning
#' @param ns A vector of accumulated number of patient at each stage
#'
#' @return A list of stan result inference
#'     stats1: A vector of posterior probability for all treatment arms
#'         including dropped and active treatment arm
#'     stats4: The mean treatment effect estimate of each treatment compared to control
#'     stats5: The variance of treatment effect estimate of each treatment compared to control
#'     post.prob.btcontrol: a vector including Posterior probability
#'         of each active treatment arm better than control
#'     stats6: A vector of mean stage (time trend) effect estimate
#'     stats7: A vector of mean treatment - stage (time trend) interaction effect estimate
#'     sampefftotal: The posterior samples of response probability
#'         of each active arm on logit scale.
#'         This can be transformed to probit scale by using inv.logit() function.
#' @export
#'
#' @examples
#' \dontrun{resultstantoRfunc.rand(group, fit, armleft, treatmentindex, K, ns)}
#' @author Ziyan Wang
resultstantoRfunc.rand = function(group, fit, armleft, treatmentindex, K, ns) {
  stats6 = rep(NA, length(ns) - 1)
  stats7 = {
  }
  if (group == 1) {
    stop("Error: Random effect model is not used for the first stage")
  }

  else if (group > 1) {
    beta0 = matrix(rstan::extract(fit, 'b_Intercept')[[1]], ncol = 1)
    beta1 = rstan::extract(fit, "beta")[[1]]
    statsbeta0 = mean(beta0+beta1[,1])
    alpha = rstan::extract(fit, "alpha")[[1]][, -1]
    sampeff = cbind(beta1[,-1]-beta1[,1], alpha)
    trteff = matrix(sampeff[, 1:(armleft - 1)], ncol = armleft - 1)
    resulttrt = resultrtostats.rand(
      trteff = trteff,
      treatmentindex = treatmentindex,
      armleft = armleft,
      K = K,
      fit = fit,
      group = group,
      ns = ns
    )
    stats4 = resulttrt$stats4
    stats5 = resulttrt$stats5
    stats1 = resulttrt$stats1
    post.prob.btcontrol = resulttrt$post.prob.btcontrol

    stageeff = matrix(sampeff[,-(1:(armleft - 1))], ncol = group - 1)
    stats6 = rep(NA, length(ns) - 1)
    names(stats6) = seq(2, length(ns))
    stats6[1:group - 1] = round(colMeans(stageeff), 3)
    #Sample distribution of reference in logistic regression
    sampefftotal = beta0 + beta1[,1]
    #Sample distribution of treatment in logistic regression
    for (temp in 1:dim(trteff)[2]) {
      sampefftotal = cbind(sampefftotal, beta0 + trteff[, temp] + 0)
    }
    #Transfer from logit scale to probability scale
    sampoutcome = inv.logit(sampefftotal)
  }
  return(
    list(
      stats1 = stats1,
      stats4 = stats4,
      stats5 = stats5,
      stats6 = stats6,
      stats7 = stats7,
      sampefftotal = sampefftotal,
      post.prob.btcontrol = resulttrt$post.prob.btcontrol
    )
  )
}
