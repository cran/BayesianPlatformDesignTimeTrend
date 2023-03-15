#' @title resultrtostats
#' @description This is an inner function of the function \code{\link{resultstantoRfunc}}

#' @param trteff Stan posterior samples of treatment effect sample distribution
#' @param treatmentindex A vector of treatment index at the beginning of a trial
#' @param armleft The number of treatment left in the platform (>2)
#' @param K Total number of arms at the beginning
#' @param group The current stage
#' @param reg.inf The information of how much accumulated information will be used
#' @param fit The stan output
#' @param ns A vector of accumulated number of patient at each stage
#'
#' @return A list of stan result inference
#'     stats1: A vector of posterior probability for all treatment arms
#'         including dropped and active treatment arm
#'     stats4: The mean treatment effect estimate of each treatment compared to control
#'     stats5: The variance of treatment effect estimate of each treatment compared to control
#'     post.prob.btcontrol: a vector including Posterior probability
#'         of each active treatment arm better than control
#' @export
#'
#' @examples
#' \dontrun{resultrtostats(trteff = NA, treatmentindex = NA, armleft, K, group, reg.inf, fit, ns)}
#' @author Ziyan Wang
resultrtostats = function(trteff = NA,
                          treatmentindex = NA,
                          armleft,
                          K,
                          group,
                          reg.inf,
                          fit,
                          ns) {
  sampeff = rstan::extract(fit, 'beta1')[[1]]
  if (group > 1 & stringr::str_detect(reg.inf, "\\*")) {
    if (stringr::str_detect(reg.inf, "continuous")) {
      interactioneff = matrix(sampeff[, (armleft + 1):(2 * armleft - 1)], ncol = armleft - 1)
      for (temp in 1:dim(trteff)[2]) {
        trteff[, temp] = trteff[, temp] + interactioneff[, temp] * group
      }
      #Mean estimates and variance of treatment effect
      beta1mean = matrix(colMeans(trteff), ncol = armleft - 1)
      colnames(beta1mean) = treatmentindex
      beta1var = matrix(sapply(data.frame(trteff), var), ncol = armleft - 1)
      colnames(beta1var) = treatmentindex

      #Record mean estimate and variance
      stats4 = rep(NA, K - 1)
      names(stats4) = seq(1, K - 1)
      stats4[treatmentindex] = round(beta1mean, 3)
      # stats4[treatmentindex] = paste0(round(beta1mean,3),"(", round(beta1var,3), ")")
      stats5 = rep(NA, K - 1)
      names(stats5) = seq(1, K - 1)
      stats5[treatmentindex] = round(beta1var, 3)
      # #posterior probability of each trt to be better than control.
      # post.prob.btcontrol = colMeans(trteff>0)
      # stats1 = rep(NA,K-1)
      # names(stats1) = seq(1,K-1)
      # stats1[treatmentindex] = post.prob.btcontrol

      #Debugged at 00:13 on 18/10/2022
      #Because there are interactions, stats1 and post.prob.btcontrol will be calculated in resultstantoRfunc
      return(list(stats4 = stats4, stats5 = stats5))
    }
    else {
      interactioneff = matrix(sampeff[,-(1:(armleft - 1 + group - 1))], ncol = (group - 1) * (armleft - 1))
      for (temp in 1:dim(trteff)[2]) {
        trteff[, temp] = trteff[, temp] + interactioneff[, (group - 1) * temp]
      }
      #Mean estimates and variance of treatment effect
      beta1mean = matrix(colMeans(trteff), ncol = armleft - 1)
      colnames(beta1mean) = treatmentindex
      beta1var = matrix(sapply(data.frame(trteff), var), ncol = armleft - 1)
      colnames(beta1var) = treatmentindex

      #Record mean estimate and variance
      stats4 = rep(NA, K - 1)
      names(stats4) = seq(1, K - 1)
      stats4[treatmentindex] = round(beta1mean, 3)
      # stats4[treatmentindex] = paste0(round(beta1mean,3),"(", round(beta1var,3), ")")
      stats5 = rep(NA, K - 1)
      names(stats5) = seq(1, K - 1)
      stats5[treatmentindex] = round(beta1var, 3)
      # #posterior probability of each arm including control to be the best.
      # post.prob.btcontrol = colMeans(trteff>0)
      # stats1 = rep(NA,K-1)
      # names(stats1) = seq(1,K-1)
      # stats1[treatmentindex] = post.prob.btcontrol

      #Debugged at 00:13 on 18/10/2022
      #Because there are interactions, stats1 and post.prob.btcontrol will be calculated in resultstantoRfunc
      return(list(stats4 = stats4, stats5 = stats5))
    }
  }
  else {
    #Mean estimates and variance of treatment effect
    beta1mean = matrix(colMeans(trteff), ncol = armleft - 1)
    colnames(beta1mean) = treatmentindex
    beta1var = matrix(sapply(data.frame(trteff), var), ncol = armleft - 1)
    colnames(beta1var) = treatmentindex

    #Record mean estimate and variance
    stats4 = rep(NA, K - 1)
    names(stats4) = seq(1, K - 1)
    stats4[treatmentindex] = round(beta1mean, 3)
    # stats4[treatmentindex] = paste0(round(beta1mean,3),"(", round(beta1var,3), ")")
    stats5 = rep(NA, K - 1)
    names(stats5) = seq(1, K - 1)
    stats5[treatmentindex] = round(beta1var, 3)
    #posterior probability of each arm including control to be the best.
    post.prob.btcontrol = colMeans(trteff > 0)
    stats1 = rep(NA, K - 1)
    names(stats1) = seq(1, K - 1)
    stats1[treatmentindex] = post.prob.btcontrol
    return(
      list(
        stats1 = stats1,
        stats4 = stats4,
        stats5 = stats5,
        post.prob.btcontrol = post.prob.btcontrol
      )
    )
  }
}
