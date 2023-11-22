#' @title Timetrend.fun
#' @description This function generate the time trend function based on trend information. This function also check the validity of the input time trend information.

#' @param trend.inf The list of information for time trend effect including 'trend.type', 'trend.effect', 'trend_add_or_multip'.
#'     'trend.type' is the shape of time trend. Default is "step". Other types are "linear", "inverse.U.linear", "plateau".
#'     "trend.effect" the vector of the strength of time trend for each arm. The first element is for the control arm.
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
          delta = (group - 1) * trend.effect
          return(delta)
        }
        timetrendornot = c("There is time trend during data generation")
      }
    },
    "linear" = {
      if (sum(trend.effect != 0) > 0) {
        trend.function = function(ns, group , i, trend.effect) {
          delta = (i - 1) * trend.effect /  (ns[length(ns)] - 1)
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
            (i - 1) * trend.effect /  (ns[length(ns)] - 1),
            (- 1 + ns[round(length(ns) / 2)]) * trend.effect /  (ns[length(ns)] - 1) - (i - 1 - ns[round(length(ns) / 2)]) * trend.effect /  (ns[length(ns)] - 1)
          )
          return(delta)
        }
        timetrendornot = c("There is time trend during data generation")
      }
    },
    "plateau" = {
      if (sum(trend.effect != 0) > 0) {
        trend.function = function(ns, group , i, trend.effect) {
          delta = trend.effect * (i - 1) / (max(ns) / 10 + (i - 1))
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
