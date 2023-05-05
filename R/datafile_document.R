#' @title Operation characteristic table for alternative scenario
#' @description Operation characteristic table for alternative scenario using main + continuousstage model. The time trend pattern is step.
#'     The strength of time trend is 0.1 equally for all arm. The effect of time trend on true response probability is multiplicative.
#' @format A data frame with 3 rows and 16 variables:
#' \describe{
#'   \item{\code{Type.I.Error.or.Power}}{Power}
#'   \item{\code{Bias.1}}{Treatment effect bias for treatment 1}
#'   \item{\code{Bias.2}}{Treatment effect bias for treatment 2}
#'   \item{\code{Bias.3}}{Treatment effect bias for treatment 3}
#'   \item{\code{rMSE.1}}{Rooted mean squared error for treatment 1}
#'   \item{\code{rMSE.2}}{Rooted mean squared error for treatment 2}
#'   \item{\code{rMSE.3}}{Rooted mean squared error for treatment 3}
#'   \item{\code{N.per.arm.1}}{Mean total number of patient allocated to control}
#'   \item{\code{N.per.arm.2}}{Mean total number of patient allocated to treatment 1}
#'   \item{\code{N.per.arm.3}}{Mean total number of patient allocated to treatment 2}
#'   \item{\code{N.per.arm.4}}{Mean total number of patient allocated to treatment 3}
#'   \item{\code{Survive.per.arm.1}}{Mean total number of patient allocated to control}
#'   \item{\code{Survive.per.arm.2}}{Mean total number of patient survived when using treatment 1}
#'   \item{\code{Survive.per.arm.3}}{Mean total number of patient survived when usin treatment 2}
#'   \item{\code{Survive.per.arm.4}}{Mean total number of patient survived when usin treatment 3}
#'   \item{\code{N}}{Mean total number of patient in a trial}
#'}

"OPC_alt"

#' @title Operation characteristic table for null scenario
#' @description Operation characteristic table for null scenario using main and main + continuousstage model. The main effect model was run for a null scenario with and without time trend.
#'  The time trend pattern is step. The strength of time trend is 0.1 equally for all arm. The effect of time trend on true response probability is multiplicative.
#' @format A data frame with 3 rows and 16 variables:
#' \describe{
#'   \item{\code{Type.I.Error.or.Power}}{Family wise error rate}
#'   \item{\code{Bias.1}}{Treatment effect bias for treatment 1}
#'   \item{\code{Bias.2}}{Treatment effect bias for treatment 2}
#'   \item{\code{Bias.3}}{Treatment effect bias for treatment 3}
#'   \item{\code{rMSE.1}}{Rooted mean squared error for treatment 1}
#'   \item{\code{rMSE.2}}{Rooted mean squared error for treatment 2}
#'   \item{\code{rMSE.3}}{Rooted mean squared error for treatment 3}
#'   \item{\code{N.per.arm.1}}{Mean total number of patient allocated to control}
#'   \item{\code{N.per.arm.2}}{Mean total number of patient allocated to treatment 1}
#'   \item{\code{N.per.arm.3}}{Mean total number of patient allocated to treatment 2}
#'   \item{\code{N.per.arm.4}}{Mean total number of patient allocated to treatment 3}
#'   \item{\code{Survive.per.arm.1}}{Mean total number of patient allocated to control}
#'   \item{\code{Survive.per.arm.2}}{Mean total number of patient survived when using treatment 1}
#'   \item{\code{Survive.per.arm.3}}{Mean total number of patient survived when usin treatment 2}
#'   \item{\code{Survive.per.arm.4}}{Mean total number of patient survived when usin treatment 3}
#'   \item{\code{N}}{Mean total number of patient in a trial}
#'}

"OPC_null"

#' @title Cutoff screening example: the recommended grid value at each time point
#' @description he recommended grid value at each time point. There are 20 cutoff value explored
#' @format A data frame with 20 rows and 1 variables:
#' \describe{
#'   \item{\code{recommandloginformd}}{The cutoff value at each time point}
#'}

"recommandloginformd"

#' @title Cutoff screening example: the predicted value from quadratic model
#' @description The predicted value from quadratic model for famliy wise error rate vs cutoff value plotting
#' @format A data frame with 1001 rows and 1 variables:
#' \describe{
#'   \item{\code{predictedtpIEinformd}}{The predicted FWER value of a large grid}
#'}

"predictedtpIEinformd"

#' @title Cutoff screening example: the details of grid
#' @description Details of grid including famliy wise error rate of a cutoff value, the cutoff value and the square of cutoff value for modelling and prediction
#' @format A data frame with 24 rows and 3 variables:
#' \describe{
#'   \item{\code{tpIE}}{FWER}
#'   \item{\code{cutoff}}{Cutoff value}
#'   \item{\code{cutoff2}}{Square of cutoff value}
#'}

"dataloginformd"

#' @title Operation characteristic table for Trial.simulation() null scenario
#' @description Operation characteristic table for null scenario using main model.
#' @format A data frame with 1 rows and 8 variables:
#' \describe{
#'   \item{\code{Type.I.Error.or.Power}}{Power}
#'   \item{\code{Bias}}{Treatment effect bias for treatment 1}
#'   \item{\code{rMSE}}{Rooted mean squared error for treatment effect 1}
#'   \item{\code{N.per.arm.1}}{Mean total number of patient allocated to control}
#'   \item{\code{N.per.arm.2}}{Mean total number of patient allocated to treatment 1}
#'   \item{\code{Survive.per.arm.1}}{Mean total number of patient allocated to control}
#'   \item{\code{Survive.per.arm.2}}{Mean total number of patient survived when using treatment 1}
#'   \item{\code{N}}{Mean total number of patient in a trial}
#'}

"OPC_Trial.simulation"

#' @title A list of data from Gaussian process for cutoff screening.
#' @description A list of data from Gaussian process for cutoff screening.
#' @format A list with two element:
#' \describe{
#'   \item{\code{next.cutoff}}{The cutoff value for the next evaluation}
#'   \item{\code{prediction}}{A list of values from Gaussian process model}
#'   \item{\code{tpIE}}{A vector of type I error rate data}
#'   \item{\code{cutoff}}{A vector of cutoff data}
#'}

"optimdata"
