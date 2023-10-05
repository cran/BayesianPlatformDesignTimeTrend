#' GP.optim: optimiser to give the next cutoff for evaluation
#' @description A function to predict the next cutoff value for evaluation.
#'
#' @param x A numeric vector of cutoff data
#' @param confidence.level A numeric value indicating the confidence level of estimate. Default is 0.95
#' @param grid.length A numeric value indicating the grid resolution. Default is 1000.
#' @param change.scale A logic value indicating whether we want to change scale when doing Gaussian process. Default is FALSE.
#' @param noise A logic value indicating whether the input x is noisy. Default is TRUE.
#' @param grid.min A numeric value or vector (for asymmetric boundary) indicating the lower bound of the grid for screening. For asymmetric boundary, the first value is efficacy minimum value and the second value is futility minimum value.
#' @param grid.max A numeric value or vector (for asymmetric boundary) indicating the upper bound of the grid for screening. For asymmetric boundary, the first value is efficacy maximum value and the second value is futility maximum value.
#' @param y.t1E A numeric vector of type I error rate data
#' @param y.pow A numeric vector of power data. You can input conjucntive, disconjunctive and marginal power data. Default is NA. Only used when Boundary.type == "Asymmetric"
#' @param Boundary.type A text indicating what type of boundary used. Default is "Symmetric"
#' @param ESS A matrix of effective sample size. This is only called for asymmetric boundary cutoff screening. Default is NA for symmetric boundary.
#'     The first column is the ESS for different cutoff pair under the null scenario, the second column is the ESS for different cutoff pair under the alternative scenario.
#' @param errorrate 'errorrate' refers to the target of type I error rate or family-wise error rate. Default is 0.05. User can change it to 0.1 for FWER if they think 0.05 is too conservative. The per-hypothesis type I error equals errorrate / (K-1) where (K-1) is the number of treatment arms.
#'
#' @return A list including the next cutoff value for evaluation `next.cutoff` and a list of predictions for screening grid.
#' @importFrom  laGP distance
#' @importFrom  stats qnorm
#' @importFrom  stats optimize
#' @export
#'
#' @examples
#' x = c(7.123968, 6.449631, 1.984406,
#' 3.507463, 4.972510, 2.925768,
#' 5.816682, 4.367796,
#' 7.349160, 1.113648)
#' y.t1E = c(0.0396, 0.0450,
#' 0.5116, 0.2172,
#' 0.1040, 0.3058,
#' 0.0592, 0.1384,
#' 0.0296, 0.7936)
#' grid.min=1
#' grid.max=8
#' GP.res=GP.optim(x=x, y.t1E=y.t1E, errorrate = 0.1, grid.min = grid.min, grid.max = grid.max)
#' GP.res$next.cutoff
#'
#' x =  data.frame(matrix(c(
#' 0.9563408, 0.006295626,
#' 0.9669739, 0.014395030,
#' 0.9959410, 0.034858339,
#' 0.9635357, 0.048435579,
#' 0.9794314, 0.021659226,
#' 0.9552201, 0.018442535,
#' 0.9655374, 0.035281833,
#' 0.9837123, 0.010656442,
#' 0.9974910, 0.047741842,
#' 0.9989172, 0.012982826), byrow=TRUE, ncol = 2))
#' y.t1E = c(0.3044, 0.2938, 0.2573, 0.4780, 0.2923, 0.3733, 0.4263, 0.1962, 0.2941, 0.1131)
#' y.pow = c(0.8300, 0.8239, 0.7102, 0.7291, 0.8205, 0.7984, 0.7709, 0.8418, 0.6359, 0.5609)
#' ESS = data.frame(matrix(c(
#' 594.672, 560.580,
#' 596.148, 566.328,
#' 597.840, 590.124,
#' 590.052, 544.800,
#' 597.024, 574.716,
#' 593.952, 554.580,
#' 593.676, 554.400,
#' 598.500, 583.896,
#' 595.740, 590.520,
#' 599.580, 598.644),byrow=TRUE,ncol=2))
#' grid.min=c(0.95,0)
#' grid.max=c(1,0.05)
#' GP.res_asy=GP.optim(x=x, y.t1E=y.t1E, y.pow=y.pow, ESS=ESS,errorrate = 0.1,
#' grid.min = grid.min, grid.max = grid.max, Boundary.type="Asymmetric")
#' GP.res_asy$next.cutoff
#' @references Surrogates: Gaussian process modeling, design, and optimization for the applied sciences. CRC press. Gramacy, R.B., 2020.
#'    Bayesian optimization for adaptive experimental design: A review. IEEE access, 8, 13937-13948. Greenhill, S., Rana, S., Gupta, S., Vellanki, P., & Venkatesh, S. (2020).
#' @author Ziyan Wang
GP.optim = function(x,
                    y.t1E,
                    y.pow=NA,
                    ESS=NA,
                    errorrate = 0.05,
                    confidence.level = 0.95,
                    grid.length = 1000,
                    change.scale = FALSE,
                    noise = T,
                    grid.min,
                    grid.max,
                    Boundary.type = "Symmetric") {
  # Debug here for GP model with Bayesian optimisation on April 25, 2023
  eps = .Machine$double.eps
  # loglikelihood function
  nlg = function(g, D, y) {
    n = length(y)
    K = exp(-D) + diag(g, n)
    Ki = solve(K)
    logdetK = determinant(K, logarithm = T)$modulus
    loglike = -(n / 2) * log(t(y) %*% Ki %*% y) - 0.5 * logdetK
    return(-loglike)
  }

  D = distance(x)
  if (noise) {
    g = optimize(nlg,
                 interval = c(eps ^ 0.5, var(y.t1E)),
                 D = D,
                 y = y.t1E)$minimum
  }
  else{
    g = eps
  }
if (Boundary.type == "Symmetric"){
  gridx = matrix(seq(grid.min, grid.max, length.out = grid.length), ncol=1)
  if (change.scale) {
    grid.new = (gridx - min(x)) / (max(x) - min(x))
    x.new = (x - min(x)) / (max(x) - min(x))
  }
  else{
    x.new = x
    grid.new = gridx
  }
}
  else{
#     if (grid.length > 101) {
#       warning(
#         "The grid length(accuracy) is too large for asymmetric boundary screening.
# This may lead to the error message:'Error: vector memory exhausted (limit reached?)'.
# This is due to the calculation of distance between each grid point (generating a (grid.length^2) rows with (grid.length^2) columns matrix).
# I recommand use grid.length = 101 which is enough for the accuracy."
#       )
#     }
    gridx.eff = seq(grid.min[1], grid.max[1], length.out = grid.length)
    gridx.fut = seq(grid.min[2], grid.max[2], length.out = grid.length)
    if (change.scale) {
      grid.new.eff = (gridx.eff - min(x[,1])) / (max(x[,1]) - min(x[,1]))
      grid.new.fut = (gridx.fut - min(x[,2])) / (max(x[,2]) - min(x[,2]))
      grid.new = expand.grid(grid.new.eff, grid.new.fut)
      x.new.eff = (x[,1] - min(x[,1])) / (max(x[,1]) - min(x[,1]))
      x.new.fut = (x[,2] - min(x[,2])) / (max(x[,2]) - min(x[,2]))
      x.new = cbind(x.new.eff, x.new.fut)
    }
    else{
      x.new = x
      grid.new = expand.grid(gridx.eff, gridx.fut)
    }
  }
  D = distance(x.new)
  K = exp(-D) + diag(g, ncol(D))
  # DXX = distance(grid.new)
  # Only take the diagonal of distance matrix which is zero for further use to save memory.
  # If user want to use the whole distance matrix which can be used to generate a full variance covariance matrix for further use,
  # Please use DXX and KXX instead as marked by #. However, the covariance is very huge if you want a fine grid.
  DXX_diag = rep(0, dim(grid.new)[[1]])
  KXX_diag = exp(-DXX_diag) + rep(g, length(DXX_diag))
  # KXX = exp(-DXX) + diag(g, ncol(DXX))
  DX = distance(grid.new, x.new)
  KX = exp(-DX)
  Ki = solve(K)

  if (Boundary.type == "Symmetric"){
    yhat.t1E = KX %*% Ki %*% y.t1E
    yhat.pow = NA
    tau_squared.t1E = drop(t(y.t1E) %*% Ki %*% y.t1E / length(y.t1E))
    sigma.t1E_diag = tau_squared.t1E * (KXX_diag - rowSums((KX %*% Ki) * KX))
    sd.t1E = sqrt(abs(sigma.t1E_diag))
    sd.pow = NA
    # construct confidence interval
    qup.t1E = yhat.t1E + qnorm(confidence.level, 0, sd.t1E)
    qdown.t1E = yhat.t1E - qnorm(confidence.level, 0, sd.t1E)
    qdown.pow = NA
    qup.pow = NA
    yhat.ESS.null = NA
    yhat.ESS.alt = NA
    sd.ESS.null = NA
    sd.ESS.alt = NA
    qup.ESS.null = NA
    qup.ESS.alt = NA
    qdown.ESS.null = NA
    qdown.ESS.alt = NA
    # GP finished
    # Bayesian optimisation start
    # construct the target. The ideal value should locate in (errorrate - 1%, errorrate + 1%)
    target = abs(yhat.t1E - errorrate) <= errorrate / 100
    # construct the potential cutoff set
    potentialcutoff = grid.new[which(target)]
    e = 1e-10
    weighs = 1 / sqrt((abs(yhat.t1E[target] - errorrate) + e) * abs(sigma.t1E_diag)[which(target)])
    randomprobability = weighs / sum(weighs)
    # GP finished
    # randomise the next value from the potential set
    # Debugged on 11/06/2023 by Ziyan wang. Cran check and find one error due to the use of sample()
    if (length(potentialcutoff) >= 1){
      next.cutoff = potentialcutoff[sample(length(potentialcutoff), 1, replace = T, prob = randomprobability)]
    }
    else {
      next.cutoff = grid.new[which.min(abs(yhat.t1E - errorrate))]
    }
  }
  else{
  # Model Type I error rate
  yhat.t1E = KX %*% Ki %*% y.t1E
  tau_squared.t1E = drop(t(y.t1E) %*% Ki %*% y.t1E / length(y.t1E))
  sigma.t1E_diag = tau_squared.t1E * (KXX_diag - rowSums((KX %*% Ki) * KX))
  sd.t1E = sqrt(abs(sigma.t1E_diag))
  # construct confidence interval
  qup.t1E = yhat.t1E + qnorm(confidence.level, 0, sd.t1E)
  qdown.t1E = yhat.t1E - qnorm(confidence.level, 0, sd.t1E)
  # Model power
  yhat.pow = KX %*% Ki %*% y.pow
  tau_squared.pow = drop(t(y.pow) %*% Ki %*% y.pow / length(y.pow))
  sigma.pow_diag = tau_squared.pow * (KXX_diag - rowSums((KX %*% Ki) * KX))
  sd.pow = sqrt(abs(sigma.pow_diag))
  # construct confidence interval
  qup.pow = yhat.pow + qnorm(confidence.level, 0, sd.pow)
  qdown.pow = yhat.pow - qnorm(confidence.level, 0, sd.pow)

  # Model ESS under null
  y.ESS.null = ESS[,1]
  yhat.ESS.null = KX %*% Ki %*% y.ESS.null
  tau_squared.ESS.null = drop(t(y.ESS.null) %*% Ki %*% y.ESS.null / length(y.ESS.null))
  sigma.ESS.null_diag = tau_squared.ESS.null * (KXX_diag - rowSums((KX %*% Ki) * KX))
  sd.ESS.null= sqrt(abs(sigma.ESS.null_diag))
  # construct confidence interval
  qup.ESS.null = yhat.ESS.null + qnorm(confidence.level, 0, sd.ESS.null)
  qdown.ESS.null = yhat.ESS.null - qnorm(confidence.level, 0, sd.ESS.null)

  # Model ESS under alternative
  y.ESS.alt = ESS[,2]
  yhat.ESS.alt = KX %*% Ki %*% y.ESS.alt
  tau_squared.ESS.alt = drop(t(y.ESS.alt) %*% Ki %*% y.ESS.alt / length(y.ESS.alt))
  sigma.ESS.alt_diag = tau_squared.ESS.alt * (KXX_diag - rowSums((KX %*% Ki) * KX))
  sd.ESS.alt = sqrt(abs(sigma.ESS.alt_diag))
  # construct confidence interval
  qup.ESS.alt = yhat.ESS.alt + qnorm(confidence.level, 0, sd.ESS.alt)
  qdown.ESS.alt = yhat.ESS.alt - qnorm(confidence.level, 0, sd.ESS.alt)
  # GP finished
  # Bayesian optimisation start
  # construct the target. The ideal value should locate in (errorrate - 1%, errorrate + 1%)
  target = abs(yhat.t1E - errorrate) <= errorrate / 100
  # construct the potential cutoff set
  potentialcutoff = grid.new[which(target),]

  # We optimise power for asymmetric boundary so we dont do weighed randomisation.
  # Instead, we pick the cutoff pair that has highest predicted power when type I error is under control.

  if (dim(potentialcutoff)[1] >= 1){
    next.cutoff = potentialcutoff[which.max(yhat.pow[as.numeric(row.names(potentialcutoff))]),]
  }
  else {
    next.cutoff = grid.new[which.min(abs(yhat.t1E - errorrate)),]
  }
}
  return(list(
    next.cutoff = next.cutoff,
    prediction = list(
      yhat.t1E = yhat.t1E,
      yhat.pow = yhat.pow,
      yhat.ESS.null = yhat.ESS.null,
      yhat.ESS.alt = yhat.ESS.alt,
      sd.t1E = sd.t1E,
      sd.pow = sd.pow,
      sd.ESS.null = sd.ESS.null,
      sd.ESS.alt = sd.ESS.alt,
      qup.t1E = qup.t1E,
      qdown.t1E = qdown.t1E,
      qup.pow = qup.pow,
      qdown.pow= qdown.pow,
      qup.ESS.null = qup.ESS.null,
      qup.ESS.alt = qup.ESS.alt,
      qdown.ESS.null = qdown.ESS.null,
      qdown.ESS.alt = qdown.ESS.alt,
      potentialcutoff = potentialcutoff,
      xgrid = grid.new
    )
  ))
}
