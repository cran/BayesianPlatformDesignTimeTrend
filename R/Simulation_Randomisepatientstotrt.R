#' @title AdaptiveRandomisation
#' @description This is a function doing the randomisation process. This Function generates the Sequence for patient allocation to each arm, patient outcomes.

#' @param Fixratio A indicator TRUE/FALSE
#' @param rand.algo Randomisation algorithm: "Coin": Biased coin; "Urn": Urn method
#' @param K Total number of arms at the beginning
#' @param n.new The cohort size
#' @param randomprob A named vector of randomisation probability to each arm
#' @param groupwise.response.probs A matrix of response probability of each arm
#' @param group The current stage
#' @param armleft The number of treatment left in the platform (>2)
#' @param max.deviation Tuning parameter of using urn randomisation method.
#' @param trend_add_or_multip How time trend affects the true response probability: "add" or "mult"
#' @param trend.function The function returns time trend effect regarding to different time trend pattern
#' @param trend.effect The strength of time trend effect as a parameter in trend.function()
#' @param treatmentindex The vector of treatment arm index excluding the control arm whose index is 0
#' @param ns A vector of accumulated number of patient at each stage
#' @param Fixratiocontrol A numeric value indicating the weight of control in randomisation.
#'     Eg. 1 means equal randomisation, 2 means thw number of patients allocated to control is twice as large as other treatment arm.
#'
#' @return A list of patient allocation and patient outcome
#'     nstage: A vector of the number of patients allocated to each arm
#'     ystage: A vector of the patients outcome after treating with each arm
#'     znew: A vector of treatment index assigned to each patient in the current cohort
#'     ynew: A vector of outcome index record for each patient after treatment in the current cohort
#' @importFrom stats rbinom
#' @importFrom boot logit
#' @importFrom boot inv.logit
#' @export
#'
#' @examples
#' AdaptiveRandomisation(
#' Fixratio = FALSE,
#' rand.algo = "Urn",
#' K = 2,
#' n.new = 30,
#' randomprob = matrix(c(0.5, 0.5), ncol = 2, dimnames = list(c(),c("1","2"))),
#' treatmentindex = 1,
#' groupwise.response.probs = matrix(rep(c(0.4, 0.4), 5), byrow = TRUE, ncol = 2, nrow = 5),
#' group = 1,
#' armleft = 2,
#' max.deviation = 3,
#' trend_add_or_multip = "mult",
#' trend.function = function(ns, group, i, trend.effect) {delta = 0; return(delta)},
#' trend.effect = c(0, 0),
#' ns = c(30, 60, 90, 120, 150),
#' Fixratiocontrol = NA)
#'
#' AdaptiveRandomisation(
#' Fixratio = TRUE,
#' rand.algo = "Urn",
#' K = 4,
#' n.new = 30,
#' randomprob = NA,
#' treatmentindex = c(1,3),
#' groupwise.response.probs = matrix(rep(c(0.4, 0.4,0.4, 0.4), 5), byrow = TRUE, ncol = 4, nrow = 5),
#' group = 1,
#' armleft = 3,
#' max.deviation = 3,
#' trend_add_or_multip = "mult",
#' trend.function = function(ns, group, i, trend.effect) {delta = 0; return(delta)},
#' trend.effect = c(0, 0),
#' ns = c(30, 60, 90, 120, 150),
#' Fixratiocontrol = 1)
#'
#' @references Mass weighted urn design—a new randomization algorithm for unequal allocations. Zhao, Wenle. Contemporary clinical trials 43 (2015): 209-216.
#' @author Ziyan Wang
AdaptiveRandomisation = function(Fixratio,
                                 rand.algo,
                                 K,
                                 n.new,
                                 randomprob,
                                 treatmentindex,
                                 groupwise.response.probs,
                                 group,
                                 armleft,
                                 max.deviation,
                                 trend_add_or_multip,
                                 trend.function,
                                 trend.effect,
                                 ns,
                                 Fixratiocontrol) {
  #unfix ratio
  if (Fixratio == FALSE) {
    #Generate assignments and outcomes for current interval
    if (rand.algo == "Coin") {
      #Assuming K arms including Control (K-1 treatment vs 1 Control)

      #Randomisation to each k arm including control arm, where rand.prob=c(AR1,AR2,AR3,...,ARK-1) for K arm allocation
      #rand.prob was suggested by Trippa et.al (2012; 2014)
      #The first element is the control group, the elements 2 to K are treatment groups

      randomsample = sample(K, n.new, replace = T, prob = randomprob)
      randcount = table(factor(randomsample, levels = seq(1, K)))
      nstage = randcount

      #Simulate the response of randomised patients
      ystage = rbinom(K, nstage, prob = groupwise.response.probs[group,])
      #Debugged due to adding time trend at 17:36 on 30/08/2022 by Ziyan Wang
      znew = randomsample
      ynew = rep(0, sum(nstage))

      m = seq(1, K)
      y.outcome.for.each.k = sapply(m, function(m) {
        ynew[which(znew == m)[sample(ystage[m])]] <- 1
        return(ynew)
      })
      ynew = rowSums(y.outcome.for.each.k)
    }
    if (rand.algo == "Urn") {
      #Data generation
      outcome = NULL
      allocation = NULL
      rand.prob.temp = randomprob
      count = matrix(rep(0, K), ncol = K)
      ycount = matrix(rep(0, K), ncol = K)
      for (i in 1:n.new) {
        urnprob = randomprob * max.deviation - as.vector(count) + (i - 1) * randomprob

        maxcon = ifelse(urnprob > 0, urnprob, 0)
        rand.prob.temp = maxcon / sum(maxcon)
        # allocation[i] = rbinom(1,1,min(max(rand.prob.temp,0),1))
        allocation[i] = sample(as.numeric(colnames(randomprob)), 1, T, prob = rand.prob.temp)
        count[allocation[i]] = count[allocation[i]] + 1


        #Time trend add to response probability. Coded due to adding time trend at 18:01 on 30/08/2022 by Ziyan Wang
        if (trend_add_or_multip == "add") {
          timetrend.response.probs = groupwise.response.probs[group,] + trend.function(ns, group , i, trend.effect)
        }
        else if (trend_add_or_multip == "mult") {
          timetrend.response.probs = inv.logit(
            logit(groupwise.response.probs[group,]) + trend.function(ns, group , i, trend.effect)
          )
        }
        else {
          stop("Error: Trend should be additive of multiplicative")
        }
        outcome[i] = rbinom(1, 1, prob = timetrend.response.probs[allocation[i]])
        ycount[allocation[i]] = ycount[allocation[i]] + outcome[i]
        # print(timetrend.response.probs)
      }

      #Data recording. Coded due to adding time trend at 18:01 on 30/08/2022 by Ziyan Wang
      randomsample = allocation
      randomoutcome = outcome
      nstage = as.vector(count)
      ystage = as.vector(ycount)
      #Debugged due to adding time trend at 17:36 on 30/08/2022 by Ziyan Wang
      znew = randomsample
      ynew = randomoutcome
    }
    #Update dataset
  }
  #Debugged at 0002 on 01112022 by Ziyan wang for adding additional fix ratio scenario
  #PS: Need additional code to be feasible to all fix ratio scenarios.
  else{
    # -----------------------Fix ratio----------------------
    rpk = matrix(rep(0, armleft), ncol = armleft)
    randomprob = matrix(rep(0, K), ncol = K)
    colnames(rpk) = c(1, treatmentindex + 1)
    colnames(randomprob) = seq(1, K)

    rpk[1] = Fixratiocontrol / (Fixratiocontrol + 1 * (armleft - 1))
    # Debugged on 02052023 by ziyan wang. The original code is not suitable for multiarm fix ratio.
    rpk[-1] = rep(1 / (armleft - 1 + Fixratiocontrol), armleft -
                    1)
    rpk = rpk / sum(rpk)
    randomprob[as.numeric(colnames(rpk))] = randomprob[as.numeric(colnames(rpk))] + rpk

    #Data generation
    outcome = NULL
    allocation = NULL
    rand.prob.temp = randomprob
    count = matrix(rep(0, K), ncol = K)
    ycount = matrix(rep(0, K), ncol = K)
    for (i in 1:n.new) {
      urnprob = randomprob * max.deviation - as.vector(count) + (i - 1) * randomprob

      maxcon = ifelse(urnprob > 0, urnprob, 0)
      rand.prob.temp = maxcon / sum(maxcon)
      # allocation[i] = rbinom(1,1,min(max(rand.prob.temp,0),1))
      allocation[i] = sample(as.numeric(colnames(randomprob)), 1, T, prob = rand.prob.temp)
      count[allocation[i]] = count[allocation[i]] + 1


      #Time trend add to response probability. Coded due to adding time trend at 18:01 on 30/08/2022 by Ziyan Wang
      if (trend_add_or_multip == "add") {
        timetrend.response.probs = groupwise.response.probs[group,] + trend.function(ns, group , i, trend.effect)
      }
      else if (trend_add_or_multip == "mult") {
        timetrend.response.probs = inv.logit(
          logit(groupwise.response.probs[group,]) + trend.function(ns, group , i, trend.effect)
        )
      }
      else {
        stop("Error: Trend should be additive of multiplicative")
      }
      outcome[i] = rbinom(1, 1, prob = timetrend.response.probs[allocation[i]])
      ycount[allocation[i]] = ycount[allocation[i]] + outcome[i]
      # print(timetrend.response.probs)
    }

    #Data recording. Coded due to adding time trend at 18:01 on 30/08/2022 by Ziyan Wang
    randomsample = allocation
    randomoutcome = outcome
    nstage = as.vector(count)
    ystage = as.vector(ycount)
    #Debugged due to adding time trend at 17:36 on 30/08/2022 by Ziyan Wang
    znew = randomsample
    ynew = randomoutcome

  }
  return(list(
    nstage = nstage,
    ystage = ystage,
    znew = znew,
    ynew = ynew
  ))
}
