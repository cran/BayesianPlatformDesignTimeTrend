#' GP.optim: optimiser to give the next cutoff for evaluation
#' @description A function to predict the next cutoff value for evaluation.
#' @param x A numeric vector of cutoff data
#' @param y A numeric vector of error rate data
#' @param errorrate A numeric value. The error rate we want to achieve. Error rate here means type I error rate or family-wise error rate. Default is 0.05.
#' @param confidence.level A numeric value indicating the confidence level of estimate. Default is 0.95
#' @param grid.length A numeric value indicating the grid resolution. Default is 5000.
#' @param change.scale A logic value indicating whether we want to change scale when doing Gaussian process. Default is FALSE.
#' @param noise A logic value indicating whether the input x is noisy. Default is TRUE.
#' @param grid.min A numeric value indicating the lower bound of the grid for screening.
#' @param grid.max A numeric value indicating the upper bound of the grid for screening.
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
#' y = c(0.0396, 0.0450,
#' 0.5116, 0.2172,
#' 0.1040, 0.3058,
#' 0.0592, 0.1384,
#' 0.0296, 0.7936)
#' grid.min=1
#' grid.max=8
#' GP.res=GP.optim(x=x, y=y, errorrate = 0.1, grid.min = grid.min, grid.max = grid.max)
#' GP.res$next.cutoff
#' @references Surrogates: Gaussian process modeling, design, and optimization for the applied sciences. CRC press. Gramacy, R.B., 2020.
#'    Bayesian optimization for adaptive experimental design: A review. IEEE access, 8, 13937-13948. Greenhill, S., Rana, S., Gupta, S., Vellanki, P., & Venkatesh, S. (2020).
#' @author Ziyan Wang
GP.optim = function(x,
                    y,
                    errorrate = 0.05,
                    confidence.level = 0.95,
                    grid.length = 5000,
                    change.scale = FALSE,
                    noise = T,
                    grid.min,
                    grid.max) {
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
                 interval = c(eps ^ 0.5, var(y)),
                 D = D,
                 y = y)$minimum
  }

  gridx = seq(grid.min, grid.max, length.out = grid.length)
  if (change.scale) {
    grid.new = (gridx - min(x)) / (max(x) - min(x))
    x.new = (x - min(x)) / (max(x) - min(x))
  }
  else{
    x.new = x
    grid.new = gridx
  }
  D = distance(x.new)
  K = exp(-D) + diag(g, ncol(D))
  DXX = distance(grid.new)
  KXX = exp(-DXX) + diag(g, ncol(DXX))
  DX = distance(grid.new, x.new)
  KX = exp(-DX)
  Ki = solve(K)

  yhat = KX %*% Ki %*% y
  tau_squared = drop(t(y) %*% Ki %*% y / length(y))
  sigma = tau_squared * (KXX - KX %*% Ki %*% t(KX))
  sd = sqrt(abs(diag(sigma)))
  # construct confidence interval
  qup = yhat + qnorm(confidence.level, 0, sd)
  qdown = yhat - qnorm(confidence.level, 0, sd)

  # GP finished
  # Bayesian optimisation start
  # construct the target. The ideal value should locate in (errorrate - 1%, errorrate + 1%)
  target = abs(yhat - errorrate) <= errorrate / 100
  # construct the potential cutoff set
  potentialcutoff = grid.new[which(target)]
  e = 1e-10
  weighs = sqrt(1 / (abs(grid.new[target] - errorrate) + e) * abs(diag(sigma))[which(target)])
  randomprobability = weighs / sum(weighs)
  # randomise the next value from the potential set
  next.cutoff = sample(potentialcutoff, 1, replace = T, prob = randomprobability)
  return(list(
    next.cutoff = next.cutoff,
    prediction = list(
      yhat = yhat,
      sd = sd,
      qup = qup,
      qdown = qdown,
      xgrid = grid.new
    )
  ))
}
