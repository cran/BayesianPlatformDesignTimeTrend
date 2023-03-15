#' @title resultstantoRfunc
#' @description This function summarise the fix effect stan output data to and transform them to be readable.

#' @param group The current stage
#' @param reg.inf The information of how much accumulated information will be used
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
#' \dontrun{resultstantoRfunc(group, reg.inf, fit, armleft, treatmentindex, K, ns)}
#' @author Ziyan Wang
resultstantoRfunc = function(group,
                             reg.inf,
                             fit,
                             armleft,
                             treatmentindex,
                             K,
                             ns) {
  sampeff = rstan::extract(fit, 'beta1')[[1]]
  sampeff0 = matrix(rstan::extract(fit, 'b_Intercept')[[1]], ncol = 1)
  if (reg.inf == "main") {
    stats6 = {
    }
    stats7 = {
    }
  }
  else if (reg.inf == "main + stage_continuous") {
    stats6 = rep(NA, 1)
    stats7 = {
    }
  }
  else if (reg.inf == "main * stage_continuous") {
    stats6 = rep(NA, 1)
    stats7 = rep(NA, K - 1)
  }
  else if (reg.inf == "main + stage_discrete") {
    stats6 = rep(NA, length(ns) - 1)
    stats7 = {
    }
  }
  else if (reg.inf == "main * stage_discrete") {
    stats6 = rep(NA, length(ns) - 1)
    stats7 = rep(NA, (K - 1) * (length(ns) - 1))
  }
  #--------------------Only Main effect----------------
  if (group == 1 | reg.inf == "main") {
    trteff = matrix(sampeff[, 1:(armleft - 1)], ncol = armleft - 1)
    resulttrt = resultrtostats(
      trteff = trteff,
      treatmentindex = treatmentindex,
      armleft = armleft,
      K = K,
      fit = fit,
      reg.inf = reg.inf,
      group = group,
      ns = ns
    )
    #stats4: Treatmenteffect mean
    #stats5: Treatmenteffect variance
    #stats1: probability of treatment better than control
    stats4 = resulttrt$stats4
    stats5 = resulttrt$stats5
    stats1 = resulttrt$stats1
    post.prob.btcontrol = resulttrt$post.prob.btcontrol

    #Sample distribution of reference in logistic regression
    sampefftotal = sampeff0
    #Sample distribution of treatment in logistic regression
    for (temp in 1:dim(trteff)[2]) {
      sampefftotal = cbind(sampefftotal, sampeff0 + trteff[, temp])
    }
    #Transfer from logit scale to probability scale
    sampoutcome = inv.logit(sampefftotal)
  }
  #--------------------Main effect and stage effect independent----------------
  else if (group > 1 & stringr::str_detect(reg.inf, "\\+")) {
    trteff = matrix(sampeff[, 1:(armleft - 1)], ncol = armleft - 1)
    resulttrt = resultrtostats(
      trteff = trteff,
      treatmentindex = treatmentindex,
      armleft = armleft,
      K = K,
      fit = fit,
      reg.inf = reg.inf,
      group = group,
      ns = ns
    )
    stats4 = resulttrt$stats4
    stats5 = resulttrt$stats5
    stats1 = resulttrt$stats1
    post.prob.btcontrol = resulttrt$post.prob.btcontrol
    if (stringr::str_detect(reg.inf, "continuous")) {
      stageeff = matrix(sampeff[,-(1:(armleft - 1))], ncol = 1)
      stats6 = rep(NA, 1)
      names(stats6) = "stageeffect"
      stats6 = round(colMeans(stageeff), 3)
      #Sample distribution of reference in logistic regression
      sampefftotal = sampeff0 + stageeff * group
      #Sample distribution of treatment in logistic regression
      for (temp in 1:dim(trteff)[2]) {
        sampefftotal = cbind(sampefftotal, sampeff0 + trteff[, temp] + stageeff * group)
      }
      #Transfer from logit scale to probability scale
      sampoutcome = inv.logit(sampefftotal)
    }
    else {
      stageeff = matrix(sampeff[,-(1:(armleft - 1))], ncol = group - 1)
      stats6 = rep(NA, length(ns) - 1)
      names(stats6) = seq(2, length(ns))
      stats6[1:group - 1] = round(colMeans(stageeff), 3)
      #Sample distribution of reference in logistic regression
      sampefftotal = sampeff0 + stageeff[, group - 1]
      #Sample distribution of treatment in logistic regression
      for (temp in 1:dim(trteff)[2]) {
        sampefftotal = cbind(sampefftotal, sampeff0 + trteff[, temp] + stageeff[, group - 1])
      }
      #Transfer from logit scale to probability scale
      sampoutcome = inv.logit(sampefftotal)
    }
  }
  #--------------------Main effect and stage effect dependent----------------
  else if (group > 1 & stringr::str_detect(reg.inf, "\\*")) {
    trteff = matrix(sampeff[, 1:(armleft - 1)], ncol = armleft - 1)
    resulttrt = resultrtostats(
      trteff = trteff,
      treatmentindex = treatmentindex,
      armleft = armleft,
      K = K,
      fit = fit,
      reg.inf = reg.inf,
      group = group,
      ns = ns
    )
    stats4 = resulttrt$stats4
    stats5 = resulttrt$stats5
    # post.prob.btcontrol = resulttrt$post.prob.btcontrol
    if (stringr::str_detect(reg.inf, "continuous")) {
      stageeff = matrix(sampeff[, armleft], ncol = 1)
      stats6 = rep(NA, 1)
      names(stats6) = "stageeffect"
      stats6 = round(colMeans(stageeff), 3)

      interactioneff = matrix(sampeff[, (armleft + 1):(2 * armleft - 1)], ncol = armleft - 1)
      stats7 = rep(NA, K - 1)
      names(stats7) = seq(1, K - 1)
      stats7[treatmentindex] = round(colMeans(interactioneff), 3)

      trteffect_with_int = {
      }
      #Sample distribution of reference in logistic regression
      sampefftotal = sampeff0 + stageeff * group
      #Sample distribution of treatment in logistic regression
      for (temp in 1:dim(trteff)[2]) {
        sampefftotal = cbind(sampefftotal,
                             sampeff0 + trteff[, temp] + stageeff * group + interactioneff[, temp] * group)
        trteffect_with_int = cbind(trteffect_with_int, trteff[, temp] + interactioneff[, temp] * group)
      }

      #posterior probability of each arm including control to be the best.
      post.prob.btcontrol = colMeans(trteffect_with_int > 0)
      stats1 = rep(NA, K - 1)
      names(stats1) = seq(1, K - 1)
      stats1[treatmentindex] = post.prob.btcontrol
      #Transfer from logit scale to probability scale
      sampoutcome = inv.logit(sampefftotal)
    }
    else {
      stageeff = matrix(sampeff[, (armleft - 1 + 1):(armleft - 1 + group - 1)], ncol = group - 1)
      stats6 = rep(NA, length(ns) - 1)
      names(stats6) = seq(2, length(ns))
      stats6[1:group - 1] = round(colMeans(stageeff), 3)

      interactioneff = matrix(sampeff[,-(1:(armleft - 1 + group - 1))], ncol = (group - 1) * (armleft - 1))
      stats7 = rep(NA, (K - 1) * (length(ns) - 1))
      names(stats7) = seq(1, (K - 1) * (length(ns) - 1))
      for (stat7temp in 1:length(treatmentindex)) {
        stats7[(1 + (length(ns) - 1) * (treatmentindex[stat7temp] - 1)):((group -
                                                                            1) + (length(ns) - 1) * (treatmentindex[stat7temp] - 1))] =
          round(colMeans(interactioneff)[(1 + (group - 1) * (stat7temp -
                                                               1)):((group - 1) * stat7temp)], 3)
      }

      trteffect_with_int = {
      }
      #Sample distribution of reference in logistic regression
      sampefftotal = sampeff0 + stageeff[, group - 1]
      #Sample distribution of treatment in logistic regression
      for (temp in 1:dim(trteff)[2]) {
        sampefftotal = cbind(sampefftotal,
                             sampeff0 + trteff[, temp] + stageeff[, group - 1] + interactioneff[, (group - 1) * temp])
        trteffect_with_int = cbind(trteffect_with_int, trteff[, temp] + interactioneff[, (group - 1) * temp])
      }

      #posterior probability of each arm including control to be the best.
      post.prob.btcontrol = colMeans(trteffect_with_int > 0)
      stats1 = rep(NA, K - 1)
      names(stats1) = seq(1, K - 1)
      stats1[treatmentindex] = post.prob.btcontrol
      #Transfer from logit scale to probability scale
      sampoutcome = inv.logit(sampefftotal)
    }
  }
  else if (group > 1 & reg.inf != "main") {
    stop("Regression information inputted wrong")
  }
  return(
    list(
      stats1 = stats1,
      stats4 = stats4,
      stats5 = stats5,
      stats6 = stats6,
      stats7 = stats7,
      sampefftotal = sampefftotal,
      post.prob.btcontrol = post.prob.btcontrol
    )
  )
}
