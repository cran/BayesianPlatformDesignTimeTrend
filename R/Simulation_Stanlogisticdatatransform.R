#' @title stan.logisticmodeltrans
#' @description This function transform the data in trial simulation to the data required for stan modelling

#' @param z A vector of all treatment index data from the beginning of a trial
#' @param y A vector of all outcome data from the beginning of a trial
#' @param randomprob A named vector of randomisation probability to each arm
#' @param group_indicator A vector for the stage at which each patient was treated
#' @param armleft The number of treatment left in the platform (>2)
#' @param group The current stage
#' @param variable.inf Fixeffect/Mixeffect for logistic model parameter
#' @param reg.inf The information of how much accumulated information will be used
#'
#' @return A list of information require for the stan model including:
#'     zdropped: The vector of treatment index for each patient
#'         whose treatment arm is active at current stage.
#'     ydropped: The vector of outcome index for each patient
#'         whose treatment arm is active at current stage.
#'     Ndropped: The total number of patients
#'         that are treated with active treatment arms at current stage.
#'     group_indicator_dropped: The vector of stage index for each patient
#'         whose treatment arm is active at current stage.
#'     zlevel: The active treatment arm index at current stage
#'     xdummy: A design matrix transformed from zdropped and group_indicator_dropped for modelling
#' @importFrom stats model.matrix
#' @export
#'
#' @examples
#' stan.logisticmodeltrans(
#' z = c(1,2,1,2,2,1,2,1),
#' y = c(0,0,0,0,1,1,1,1),
#' randomprob = matrix(c( 0.5, 0.5), ncol = 2, dimnames = list(c("Stage1"), c("1", "2"))),
#' group_indicator = c(1,1,1,1,1,1,1,1),
#' armleft = 2,
#' group = 1,
#' variable.inf = "Fixeffect",
#' reg.inf = "main")
#' @author Ziyan Wang
stan.logisticmodeltrans = function(z,
                                   y,
                                   randomprob,
                                   group_indicator,
                                   armleft,
                                   group,
                                   variable.inf,
                                   reg.inf) {
  zdropped = z[z %in% as.numeric(names(randomprob[, randomprob != 0]))]
  ydropped = y[z %in% as.numeric(names(randomprob[, randomprob != 0]))]
  Ndropped = length(zdropped)
  group_indicator_dropped = group_indicator[z %in% as.numeric(names(randomprob[, randomprob != 0]))]
  zlevel = as.numeric(levels(factor(zdropped)))

  for (zindex in 1:armleft) {
    zdropped[zdropped == levels(factor(zdropped))[zindex]] = zindex
  }
  #Construct data set x for different regression models because of rewriting stan model.
  #Debugged at 18:29 16/09/2022 by ZIYAN WANG.
  if (group == 1) {
    #Stage 1, there is no stage effect and interaction effect
    xdummy = model.matrix( ~ factor(zdropped))
  }
  else if (variable.inf == "Mixeffect") {
    xdummy = model.matrix( ~ factor(zdropped))
  }
  else if (group > 1 &
           reg.inf == "main" & variable.inf == "Fixeffect") {
    xdummy = model.matrix( ~ factor(zdropped))
  }
  else if (group > 1 &
           reg.inf == "main + stage_continuous" &
           variable.inf == "Fixeffect") {
    xdummy = model.matrix( ~ factor(zdropped) + group_indicator_dropped)
  }
  else if (group > 1 &
           reg.inf == "main * stage_continuous" &
           variable.inf == "Fixeffect") {
    xdummy = model.matrix( ~ factor(zdropped) * group_indicator_dropped)
  }
  else if (group > 1 &
           reg.inf == "main + stage_discrete" &
           variable.inf == "Fixeffect") {
    xdummy = model.matrix( ~ factor(zdropped) + as.factor(group_indicator_dropped))
  }
  else if (group > 1 &
           reg.inf == "main * stage_discrete" &
           variable.inf == "Fixeffect") {
    xdummy = model.matrix( ~ factor(zdropped) * as.factor(group_indicator_dropped))
  }
  return(
    list(
      zdropped = zdropped,
      ydropped = ydropped,
      Ndropped = Ndropped,
      group_indicator_dropped = group_indicator_dropped,
      zlevel = zlevel,
      xdummy = xdummy
    )
  )
}

