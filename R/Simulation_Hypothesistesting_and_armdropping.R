#' @title testing_and_armdropping
#' @description This function makes a decision on whether any active arm should be dropped based on posterior probability and
#'     return the vector of decision on each arm, the vector of active arms index and the number of arms left for further study.
#'
#' @param post.prob.btcontrol A numeric vector of posterior probability of each treatment arm better than control
#' @param group A numeric value. The current stage index.
#' @param cutoffeff A numeric vector of the cutoff value at each stage for efficacy boundary.
#' @param cutoffful A numeric vector of the cutoff value at each stage for futility boundary.
#' @param treatmentindex A numeric vector of the current active treatment arm index
#' @param test.type A character indicating which hypothesis testing we are use.
#'     "Oneside": H_0: \\pi_k \\leq \\pi_0; H_0: \\pi_k \> \\pi_0
#'     "Twoside": H_0: \\pi_k \\eq \\pi_0; H_0: \\pi_k \\neq \\pi_0
#' @param K A numeric value indicating the total number of arm at the beginning of trial including both control and treatment.
#' @param armleft A numeric vector indicating the number of active arms before this interim analysis;
#'
#' @return A list of information including armleft: the number of active arms after this interim analysis;
#'     treatmentindex: the index vector of active arm after this interim analysis;
#'     stats3: the vector of conclusion on whether null hypothesis is rejected
#' @export
#'
#' @examples
#' testing_and_armdropping(
#' K = 4,
#' armleft = 4,
#' post.prob.btcontrol = c(0.5,0.99,0.02),
#' group = 3,
#' cutoffeff = c(1, 0.99, 0.975, 0.96, 0.95),
#' cutoffful = c(0, 0.01, 0.025, 0.04, 0.05),
#' treatmentindex = c(1,2,3),
#' test.type = "Oneside")
#'
#' testing_and_armdropping(
#' K = 4,
#' armleft = 4,
#' post.prob.btcontrol = c(0.5,0.99,0.02),
#' group = 3,
#' cutoffeff = c(1, 0.99, 0.975, 0.96, 0.95),
#' cutoffful = c(0, 0.01, 0.025, 0.04, 0.05),
#' treatmentindex = c(1,2,3),
#' test.type = "Twoside")
testing_and_armdropping = function(K,
                                   armleft,
                                   post.prob.btcontrol,
                                   group,
                                   cutoffeff,
                                   cutoffful,
                                   treatmentindex,
                                   test.type) {
  #----Justify if type I error was made for each arm----
  # post.prob.btcontrol>cutoffeff[group]: Efficacy boundary is hit at this stage
  # post.prob.btcontrol<cutoffful[group]: Fultility boundary is hit at this stage
  Justify = post.prob.btcontrol > cutoffeff[group] |
    post.prob.btcontrol < cutoffful[group]
  # If one side test, we would only make conclusion on arm for superiority.
  # However, if two side test, we could make conclusion on arm for either inferiority and superiority.
  Conclude.efficacy = post.prob.btcontrol > cutoffeff[group]
  Conclude.twoside = Justify
  #Identify which active arm should be dropped at current stage
  treatmentdrop = treatmentindex[post.prob.btcontrol > cutoffeff[group] |
                                   post.prob.btcontrol < cutoffful[group]]
  # Delete the arm dropped in this round information
  post.prob.btcontrol = post.prob.btcontrol[!Justify]

  stats3 = rep(NA, K - 1)
  names(stats3) = seq(1, K - 1)
  if (test.type == "Oneside") {
    stats3[treatmentindex] = Conclude.efficacy
  }
  else if (test.type == "Twoside") {
    stats3[treatmentindex] = Conclude.twoside
  }
  else{
    stop("The hypothesis testing type should be specified")
  }
  if (sum(Justify) > 0) {
    armleft = armleft - sum(Justify)
    #Debugged for K arm by Ziyan Wang on 12:00 26/07/2022 for three arm. Used to be treatmentindex = treatmentindex[-treatmentdrop]
    #Debugged for K arm by Ziyan Wang on 18:58 26/07/2022 for more than 3 arm. Used to be treatmentindex = treatmentindex[!(treatmentindex==treatmentdrop)]
    treatmentindex = treatmentindex[is.na(match(treatmentindex, treatmentdrop))]
  }

  return(list(
    stats3 = stats3,
    armleft = armleft,
    treatmentindex = treatmentindex
  ))
}
