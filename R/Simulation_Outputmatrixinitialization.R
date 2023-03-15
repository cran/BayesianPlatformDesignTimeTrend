#' @title OutputStats.initialising
#' @description This function initializes the output matrix including all evaluation metrics based on input information

#' @param variable.inf  The parameter information in the model
#' @param reg.inf  The model information. For the fixed effect model, the input of reg.inf can be main, main + stage_continuous, main * stage_continuous, main + stage_discrete,
#'      main * stage_discrete.
#'      For the mixed effect model, the reg.inf is invalid.
#' @param ns  A vector of accumulated number of patient at each stage
#' @param K  Total number of arm including control
#'
#' @return The empty output matrix including different evaluation metrics.
#' @export
#'
#' @examples
#' OutputStats.initialising(
#'     variable.inf = "Fixeffect",
#'     reg.inf = "main",
#'     ns = c(15, 30, 45, 60, 75),
#'     K = 2)
#' @author Ziyan Wang
OutputStats.initialising = function(variable.inf, reg.inf, ns, K) {
  #Storage object for posterior probabilities
  if (variable.inf == "Fixeffect") {
    if (reg.inf == "main") {
      stats = matrix(NA,
                     nrow = length(ns),
                     ncol = K - 1 + K * 2 + K - 1 + 1 + K - 1 + K - 1)
      #K-1 posterior probability better than control and K columns number of success + K columns number of patient
      #K-1 Indicator of K-1 hypothesis
      #1: control mean estimates Debugged at 02:25 10/10/2022 by ZIYAN WANG
      #K-1 Mean estimates of treatment effect at logit scale. Debugged at 14:52 05/09/2022 by ZIYAN WANG
      #K-1 Variance of treatment effect at logit scale. Debugged at 14:52 05/09/2022 by ZIYAN WANG
    }
    else if (reg.inf == "main + stage_continuous") {
      stats = matrix(NA,
                     nrow = length(ns),
                     ncol = K - 1 + K * 2 + K - 1 + 1 + K - 1 + K - 1 + 1)
      #K-1 posterior probability better than control and K columns number of success + K columns number of patient
      #K-1 Indicator of K-1 hypothesis
      #1: control mean estimates Debugged at 02:25 10/10/2022 by ZIYAN WANG
      #K-1 Mean estimates of treatment effect at logit scale. Debugged at 14:52 05/09/2022 by ZIYAN WANG
      #K-1 Variance of treatment effect at logit scale. Debugged at 14:52 05/09/2022 by ZIYAN WANG
      #1 Stage effect when stage is treated as continuous. Debugged at 22:00 16/09/2022 by ZIYAN WANG
    }
    else if (reg.inf == "main * stage_continuous") {
      stats = matrix(NA,
                     nrow = length(ns),
                     ncol = K - 1 + K * 2 + K - 1 + 1 + K - 1 + K - 1 + 1 + K - 1)
      #K-1 posterior probability better than control and K columns number of success + K columns number of patient
      #K-1 Indicator of K-1 hypothesis
      #1: control mean estimates Debugged at 02:25 10/10/2022 by ZIYAN WANG
      #K-1 Mean estimates of treatment effect at logit scale. Debugged at 14:52 05/09/2022 by ZIYAN WANG
      #K-1 Variance of treatment effect at logit scale. Debugged at 14:52 05/09/2022 by ZIYAN WANG
      #1 Stage effect when stage is treated as continuous. Debugged at 22:00 16/09/2022 by ZIYAN WANG
      #K-1 Interaction effect when stage is treated as continuous. Debugged at 22:00 16/09/2022 by ZIYAN WANG
    }
    else if (reg.inf == "main + stage_discrete") {
      stats = matrix(NA,
                     nrow = length(ns),
                     ncol = K - 1 + K * 2 + K - 1 + 1 + K - 1 + K - 1 + length(ns) - 1)
      #K-1 posterior probability better than control and K columns number of success + K columns number of patient
      #K-1 Indicator of K-1 hypothesis
      #1: control mean estimates Debugged at 02:25 10/10/2022 by ZIYAN WANG
      #K-1 Mean estimates of treatment effect at logit scale. Debugged at 14:52 05/09/2022 by ZIYAN WANG
      #K-1 Variance of treatment effect at logit scale. Debugged at 14:52 05/09/2022 by ZIYAN WANG
      #length(ns)-1 Total number of Stage effect when stage is treated as discrete Debugged at 22:00 16/09/2022 by ZIYAN WANG
    }
    else if (reg.inf == "main * stage_discrete") {
      stats = matrix(
        NA,
        nrow = length(ns),
        ncol = K - 1 + K * 2 + K - 1 + 1 + K - 1 + K - 1 + length(ns) - 1 + (K -
                                                                               1) * (length(ns) - 1)
      )
      #K-1 posterior probability better than control and K columns number of success + K columns number of patient
      #K-1 Indicator of K-1 hypothesis
      #1: control mean estimates Debugged at 02:25 10/10/2022 by ZIYAN WANG
      #K-1 Mean estimates of treatment effect at logit scale. Debugged at 14:52 05/09/2022 by ZIYAN WANG
      #K-1 Variance of treatment effect at logit scale. Debugged at 14:52 05/09/2022 by ZIYAN WANG
      #length(ns)-1 Total number of Stage effect when stage is treated as discrete Debugged at 22:00 16/09/2022 by ZIYAN WANG
      #(K-1)*(length(ns)-1) Total number of Interaction effect when stage is treated as discrete Debugged at 22:00 16/09/2022 by ZIYAN WANG
    }
    else{
      stop(
        "Error: reg.inf must be in c(main, main + stage_continuous, main * stage_continuous, main + stage_discrete, main * stage_discrete)"
      )
    }
    #Debugged at 02:01 17/09/2022 by ZIYAN WANG. For compressing different version of codes.
    #Use name function to assign the names of each column of stats matrix
    namefunc = function(K, ns) {
      ngroup = length(ns)
      #probability of each treatment better than control
      name1 = {
      }
      #number of patients in each arm and their outcomes
      name2 = {
      }
      #Hypothesis testing result
      name3 = {
      }
      #mean treatment effect
      name4 = {
      }
      #varance of treatment effect
      name5 = {
      }
      #mean Stage effect
      name6 = {
      }
      #mean Interaction effect
      name7 = {
      }
      #mean fix intercept
      name8 = paste0("Intercept")

      if (reg.inf == "main") {
        for (n in 1:(K - 1)) {
          name1 = c(name1, paste0("PP", n, "C"))
          name2 = c(name2, paste0("nE", n), paste0("yE", n))
          name3 = c(name3, paste0("H1^", n, "tpIE"))
          name4 = c(name4, paste0("Trt", n, "_Mean"))
          name5 = c(name5, paste0("Trt", n, "_Var"))
        }
      }
      else if (reg.inf == "main + stage_continuous") {
        for (n in 1:(K - 1)) {
          name1 = c(name1, paste0("PP", n, "C"))
          name2 = c(name2, paste0("nE", n), paste0("yE", n))
          name3 = c(name3, paste0("H1^", n, "tpIE"))
          name4 = c(name4, paste0("Trt", n, "_Mean"))
          name5 = c(name5, paste0("Trt", n, "_Var"))
        }
        name6 = c(name6, paste0("stageeffect"))
      }
      else if (reg.inf == "main * stage_continuous") {
        for (n in 1:(K - 1)) {
          name1 = c(name1, paste0("PP", n, "C"))
          name2 = c(name2, paste0("nE", n), paste0("yE", n))
          name3 = c(name3, paste0("H1^", n, "tpIE"))
          name4 = c(name4, paste0("Trt", n, "_Mean"))
          name5 = c(name5, paste0("Trt", n, "_Var"))
          name7 = c(name7, paste0("stage_treatment", n, "interaction"))
        }
        name6 = c(name6, paste0("stageeffect"))
      }
      else if (reg.inf == "main + stage_discrete") {
        for (n in 1:(K - 1)) {
          name1 = c(name1, paste0("PP", n, "C"))
          name2 = c(name2, paste0("nE", n), paste0("yE", n))
          name3 = c(name3, paste0("H1^", n, "tpIE"))
          name4 = c(name4, paste0("Trt", n, "_Mean"))
          name5 = c(name5, paste0("Trt", n, "_Var"))
        }
        for (grouptemp in 2:ngroup) {
          name6 = c(name6, paste0("stageeffect_", grouptemp))
        }
      }
      else if (reg.inf == "main * stage_discrete") {
        for (n in 1:(K - 1)) {
          name1 = c(name1, paste0("PP", n, "C"))
          name2 = c(name2, paste0("nE", n), paste0("yE", n))
          name3 = c(name3, paste0("H1^", n, "tpIE"))
          name4 = c(name4, paste0("Trt", n, "_Mean"))
          name5 = c(name5, paste0("Trt", n, "_Var"))
          for (grouptemp in 2:ngroup) {
            name7 = c(name7,
                      paste0("trt", n, "_stage", grouptemp, "interaction"))
          }
        }
        for (grouptemp in 2:ngroup) {
          name6 = c(name6, paste0("stageeffect_", grouptemp))
        }
      }
      return(c(
        name1,
        "nC",
        "yC",
        name2,
        name3,
        name8,
        name4,
        name5,
        name6,
        name7
      ))
    }
  }
  else if (variable.inf == "Mixeffect") {
    stats = matrix(NA,
                   nrow = length(ns),
                   ncol = K - 1 + K * 2 + K - 1 + 1 + K - 1 + K - 1 + length(ns) - 1)
    #K-1 posterior probability better than control and K columns number of success + K columns number of patient
    #K-1 Indicator of K-1 hypothesis
    #1: control mean estimates Debugged at 02:25 10/10/2022 by ZIYAN WANG
    #K-1 Mean estimates of treatment effect at logit scale. Debugged at 14:52 05/09/2022 by ZIYAN WANG
    #K-1 Variance of treatment effect at logit scale. Debugged at 14:52 05/09/2022 by ZIYAN WANG
    #length(ns)-1 Total number of Stage effect when stage is treated as discrete Debugged at 22:00 16/09/2022 by ZIYAN WANG
    namefunc = function(K, ns) {
      ngroup = length(ns)
      #probability of each treatment better than control
      name1 = {
      }
      #number of patients in each arm and their outcomes
      name2 = {
      }
      #Hypothesis testing result
      name3 = {
      }
      #mean treatment effect
      name4 = {
      }
      #varance of treatment effect
      name5 = {
      }
      #mean Stage effect
      name6 = {
      }
      #mean Interaction effect
      name7 = {
      }
      #mean fix intercept
      name8 = paste0("Intercept")

      for (n in 1:(K - 1)) {
        name1 = c(name1, paste0("PP", n, "C"))
        name2 = c(name2, paste0("nE", n), paste0("yE", n))
        name3 = c(name3, paste0("H1^", n, "tpIE"))
        name4 = c(name4, paste0("Trt", n, "_Mean"))
        name5 = c(name5, paste0("Trt", n, "_Var"))
      }
      for (grouptemp in 2:ngroup) {
        name6 = c(name6, paste0("Randomintercept", grouptemp))
      }
      return(c(
        name1,
        "nC",
        "yC",
        name2,
        name3,
        name8,
        name4,
        name5,
        name6,
        name7
      ))
    }
  }
  statsname = namefunc(K, ns)
  colnames(stats) = statsname
  rownames(stats) = 1:length(ns)
  return(stats)
}
