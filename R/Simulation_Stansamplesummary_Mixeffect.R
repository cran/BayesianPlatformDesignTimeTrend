#' @title resultrtostats.rand
#' @description The inner function of function \code{\link{resultstantoRfunc.rand}}

#' @param trteff Stan posterior samples of treatment effect sample distribution
#' @param treatmentindex A vector of treatment index at the beginning of a trial
#' @param armleft The number of treatment left in the platform (>2)
#' @param K Total number of arms at the beginning
#' @param group The current stage
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
#' \dontrun{resultrtostats.rand(trteff = NA, treatmentindex = NA, armleft, K, group, fit, ns)}
#' @author Ziyan Wang
resultrtostats.rand = function(trteff = NA,
                               treatmentindex = NA,
                               armleft,
                               K,
                               group,
                               fit,
                               ns) {
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
