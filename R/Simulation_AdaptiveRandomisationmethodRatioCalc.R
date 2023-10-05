#' @title ARmethod
#' @description This function adjusts the posterior randomisation probability for each arm using many approaches.
#'     Currently Thall's approach and Trippa's approach are used.
#'     Double biased coin and other method will be added in the next version.
#' @param Fixratio A indicator TRUE/FALSE
#' @param BARmethod The indicator of which adaptive randomisation method is used
#' @param group The current stage
#' @param stats The output matrix
#' @param post.prob.btcontrol The vector of posterior probability of each active treatment arm better than control
#' @param K Total number of arms at the beginning
#' @param n The vector of sample size for each arm
#' @param tuningparameter The tuning parameter indicator for Thall's approach
#' @param c The tuning parameter for Thall's approach
#' @param a The hyperparamter parameter for Trippa's approach
#' @param b The hyperparamter parameter for Trippa's approach
#' @param post.prob.best Posterior probability of each arm to be the best
#' @param max.ar The upper boundary for randomisation ratio for each arm, which is used in Thall's approach since Trippa's approach has protection on control arm.
#' @param armleft The number of treatment left in the platform (>2)
#' @param treatmentindex The vector of treatment arm index excluding the control arm whose index is 0
#'
#'
#' @return randomprob: The vector of adjusted randomisation probability to each arm
#' @export
#'
#' @examples
#' ARmethod(Fixratio = FALSE,
#' BARmethod = "Thall",
#' group = 1,
#' stats = matrix(rep(NA, 40), ncol = 8, nrow = 5),
#' post.prob.btcontrol = 0.5,
#' K = 2,
#' n = c(30, 30),
#' tuningparameter = "fixed",
#' c = 1,
#' post.prob.best = c(0.5, 0.5),
#' max.ar = 0.75,
#' armleft = 2,
#' treatmentindex = 1)
#'
#' ARmethod(Fixratio = FALSE,
#' BARmethod = "Trippa",
#' group = 1,
#' stats = matrix(rep(NA, 40), ncol = 8, nrow = 5),
#' post.prob.btcontrol = c(0.5, 0.6),
#' K = 3,
#' n = c(30, 30, 40),
#' tuningparameter = NA,
#' c = NA,
#' a = 3,
#' b = 0.75,
#' post.prob.best = c(0.3, 0.3, 0.4),
#' max.ar = NA,
#' armleft = 3,
#' treatmentindex = c(1, 2))
#'
#' @references Bayesian adaptive randomized trial design for patients with recurrent glioblastoma. Trippa, Lorenzo, Eudocia Q. Lee, Patrick Y. Wen, Tracy T. Batchelor, Timothy Cloughesy, Giovanni Parmigiani, and Brian M. Alexander. Journal of Clinical Oncology 30, no. 26 (2012): 3258.
#'     A simulation study of outcome adaptive randomization in multi-arm clinical trials. Wathen, J. Kyle, and Peter F. Thall. Clinical Trials 14, no. 5 (2017): 432-440.
#'
#' @author Ziyan Wang
ARmethod = function(Fixratio,
                    BARmethod,
                    group,
                    stats,
                    post.prob.btcontrol,
                    K,
                    n,
                    tuningparameter = NA,
                    c = NA,
                    a = NA,
                    b = NA,
                    post.prob.best,
                    max.ar = NA,
                    armleft,
                    treatmentindex) {
  # Validate inputs
  if (!is.logical(Fixratio)) stop("Error: Fixratio should be a logical value (TRUE/FALSE)")
  if (!is.character(BARmethod)) stop("Error: BARmethod should be a character value")
  if (!is.numeric(group) || group <= 0) stop("Error: group should be a positive numeric value")
  if (BARmethod == "Thall" & (!is.numeric(max.ar) || max.ar <= 0 || max.ar >= 1)) stop("max.ar should be a numeric value between 0 and 1 for Thall's approach")
  if (BARmethod == "Trippa" & (!is.numeric(a) || !is.numeric(b))) stop("hyperparameters a and b should be a numeric value for Trippa's approach")
  if (K < 2) stop("Error: K should be an integer value greater than or equal to 2")

  #---------------------Trippa's approach---------------------
  if (Fixratio == F & BARmethod == "Trippa") {
    ##Tuning the paprameter using method mentioned in Trippa's paper (2014)
    gamma_stage = a * ((group / dim(stats)[1])) ^ b
    eta_stage = 0.25 * (group / dim(stats)[1])
    ##Reweigh the allocation probability
    ###K >= 1, treatment group

    allocate_trt = post.prob.btcontrol ^ gamma_stage / sum(post.prob.btcontrol ^
                                                             gamma_stage)
    ###k = 0, control group
    allocate_control = 1 / (armleft - 1) * (exp(max(n[-1]) - n[1])) ^ eta_stage

    sum_pi = allocate_control + sum(allocate_trt)
    alloc.prob.btcontrol = c(allocate_control / sum_pi, allocate_trt / sum_pi)

    rpk = matrix(rep(0,armleft),ncol = armleft)
    randomprob = matrix(rep(0,K),ncol = K)
    colnames(rpk) = c(1,treatmentindex+1)
    colnames(randomprob) = seq(1,K)
    rpk[1] = alloc.prob.btcontrol[1]
    rpk[-1] = alloc.prob.btcontrol[-1]
    rpk = rpk/sum(rpk)
    randomprob[as.numeric(colnames(rpk))] = randomprob[as.numeric(colnames(rpk))]+rpk
  }

  #---------------------Thall's approach---------------------
  else if (Fixratio == F & BARmethod == "Thall") {
    ##Tuning parameter c for Thall's approach

    if (tuningparameter == "Unfixed") {
      c = group / (2 * dim(stats)[1])
    }
    else {
      c = c
    }

    ##Reweigh the allocation probability
    alloc.prob.best = post.prob.best ^ c / sum(post.prob.best ^ c)
    rpblocktwoarm = min(max.ar, max(1 - max.ar, post.prob.btcontrol))

    randomprob = alloc.prob.best

    #-------------------Allocation bounds restriction (two arm)---------------
    if (K == 2) {
      lower = ifelse(alloc.prob.best < (1 - max.ar), 1 - max.ar, alloc.prob.best)
      upper = ifelse(lower > max.ar, max.ar, lower)
      randomprob = upper
      randomprob = matrix(randomprob, ncol = length(upper))
      colnames(randomprob) = seq(1, K)
    }
    #-------------------Allocation bounds restriction (K arm: restriction on control)---------------
    else if (K > 2) {
      rpk = matrix(rep(0,armleft),ncol = armleft)
      randomprob = matrix(rep(0,K),ncol = K)
      colnames(rpk) = c(1,treatmentindex+1)
      colnames(randomprob) = seq(1,K)
      rpk[1] = alloc.prob.best[1]
      rpk[-1] = alloc.prob.best[-1][treatmentindex]
      rpk = rpk/sum(rpk)
      lower = ifelse(rpk<(1-max.ar),1-max.ar,rpk)
      rpk = lower
      upper = ifelse(rpk>max.ar,max.ar,rpk)
      rpk = upper
      rpk[!(rpk==(1-max.ar))]=(1-sum(rpk[rpk==1-max.ar]))*(rpk[!(rpk==1-max.ar)]/sum(rpk[!(rpk==1-max.ar)]))
      randomprob[as.numeric(colnames(rpk))] = randomprob[as.numeric(colnames(rpk))]+rpk
    }

  }
  else {
    stop("Error: Please check the input of Fixratio and BARmethod")
  }
  return(randomprob)
}
