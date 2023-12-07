#' @title Timetrend.fun
#' @description This function generate the time trend function based on trend information. This function also check the validity of the input time trend information.

#' @param trend.inf The list of information for time trend effect including 'trend.type', 'trend.effect', 'trend_add_or_multip'.
#'     'trend.type' is the shape of time trend. Default is "step". Other types are "linear", "inverse.U.linear", "plateau".
#'     "trend.effect" the vector of the strength of time trend for each arm. The first element is for the control arm.
#'         The value of time trend is the gap between the start of the trial and the end of the trial. The change between each interim or each patient is calculated in the function.
#'         For example, for linear time trend with trend.effect = c(0.2, 0.2). The trend effect increment in control group for patient $i$ is 0.2(i-1)/(N_max-1), for stage $j$ is 0.2(j-1)/(length(ns)-1).
#'     "trend_add_or_multip" the pattern of time trend affecting the true response probability. Default is "mult".
#'
#' @return A list containing the time trend function according to input trend.type variable,
#'     and a indicator of whether there is a time trend in data generation
#'     based on input trend information
#' @export
#'
#' @examples
#' Timetrend.fun(trend.inf = list(
#'                trend.type = "step",
#'                trend.effect = c(0, 0),
#'                trend_add_or_multip = "mult"
#'                ))
#' @author Ziyan Wang
Timetrend.fun = function(trend.inf) {
  # Time trend pattern function
  trend.type = trend.inf$trend.type
  trend.effect = trend.inf$trend.effect
  trend_add_or_multip = trend.inf$trend_add_or_multip

  switch(
    trend.type,
    "step" = {
      if (sum(trend.effect != 0) > 0) {
        trend.function = function(ns, group , i, trend.effect) {
          delta = (group - 1) * trend.effect/(length(ns)-1)
          return(delta)
        }
        timetrendornot = c("There is time trend during data generation")
      }
    },
    "linear" = {
      if (sum(trend.effect != 0) > 0) {
        trend.function = function(ns, group , i, trend.effect) {
          delta = (i - 1 + ns[group] - ns[1]) * trend.effect /  (ns[length(ns)] - 1)
          return(delta)
        }
        timetrendornot = c("There is time trend during data generation")
      }
    },
    "inverse.U.linear" = {
      if (sum(trend.effect != 0) > 0) {
        trend.function = function(ns, group , i, trend.effect) {
          delta = ifelse(
            group <= round(length(ns) / 2),
            (i - 1 + ns[group] - ns[1]) * trend.effect /  (ns[length(ns)] - 1),
            (ns[1] - 1 + ns[round(length(ns) / 2)] - ns[1]) * trend.effect /  (ns[length(ns)] - 1) - (i - 1 +
                                                                                                        ns[group - round(length(ns) / 2)] - ns[1]) * trend.effect /  (ns[length(ns)] - 1)
          )
          return(delta)
        }
        timetrendornot = c("There is time trend during data generation")
      }
    },
    "plateau" = {
      if (sum(trend.effect != 0) > 0) {
        trend.function = function(ns, group , i, trend.effect) {
          delta = trend.effect * (i - 1 + ns[group] - ns[1]) / (max(ns) / 10 + (i -
                                                                                  1 + ns[group] - ns[1]))
          return(delta)
        }
        timetrendornot = c("There is time trend during data generation")
      }
    },
    stop(
      "Error: Wrong trend type or strength of time effect for data generation"
    )
  )

  if (sum(trend.effect != 0) == 0) {
    trend.function = function(ns, group, i, trend.effect) {
      delta = 0
      return(delta)
    }
    timetrendornot = c("There is no time trend during data generation")
  }

  return(
    list(
      trend.function = trend.function,
      timetrendornot = timetrendornot,
      trend_add_or_multip = trend_add_or_multip,
      trend.effect = trend.effect
    )
  )
}
