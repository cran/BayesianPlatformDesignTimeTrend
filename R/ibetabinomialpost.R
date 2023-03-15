#' @title ibetabinomial.post
#' @description This function calculates the posterior probability of each active treatment arm better than control using betabinomial model
#'
#' @param n A vector of treated patients for each arm (The first element is for control)
#' @param y A vector of treated patient outcomes for each arm (The first element is for control)
#' @param pi.star The prior response probability. The default is 0.5
#' @param pess The effective sample size of beta prior. The default is 2
#'
#' @return A vector posterior probability of each active treatment arm better than control
#' @importFrom stats dbeta
#' @importFrom stats integrate
#' @importFrom stats lm
#' @importFrom stats pbeta
#' @export
#'
#' @examples
#' n <- c(20,20,20,20)
#' y <- c(12,12,12,6)
#' ibetabinomial.post(n, y, pi.star = 0.5, pess = 2)
#' #[1] 0.5000000 0.5000000 0.0308018
#' @author Ziyan Wang
ibetabinomial.post = function(n, y, pi.star = 0.5, pess = 2) {
  #First element of n and y are from control
  K = length(n)
  # Posterior Probability of each arm better than the same control arm
  # See https://www.evanmiller.org/bayesian-ab-testing.html#cite1
  # Treatment_k~beta(a_k,b_k), Control~beta(a_1,b_1)

  # p.prior*ess.prior: prior success
  # (1-p.prior)*ess.prior: prior failure


  post.prob = {
  }
  for (k in 2:K) {
    post.prob[k - 1] <- unlist(integrate(
      function(x)
        pbeta(x, y[k] + pi.star * pess, (n[k] - y[k]) + (1 - pi.star) * pess, lower.tail =
                FALSE) * dbeta(x, y[1] + pi.star * pess, (n[1] - y[1]) + (1 - pi.star) *
                                 pess),
      lower = 0,
      upper = 1
    ))$value
    names(post.prob[k - 1]) = paste("Treatment", k - 1, "vs", "Control", sep = " ")
  }

  #Slower version
  # # p.prior*ess.prior: prior success
  # # (1-p.prior)*ess.prior: prior failure
  # narm=length(n)
  # rn=matrix(rbeta(random.number*narm,
  #                 y+p.prior*ess.prior,
  #                 n-y+(1-p.prior)*ess.prior),
  #           random.number, byrow = TRUE)
  #
  # controlrn=rn[,1]
  # treatmentrn=as.matrix(rn[,-1])
  #
  # # rnT=rbeta(random.number,y[1]+p.prior*ess.prior,n[1]-y[1]+(1-p.prior)*ess.prior)
  # # rnC=rbeta(random.number,y[2]+p.prior*ess.prior,n[2]-y[2]+(1-p.prior)*ess.prior)
  # # postprob<-mean(rnT>rnC)
  #
  # post.prob=colMeans(treatmentrn>controlrn)

  return(post.prob)
}
