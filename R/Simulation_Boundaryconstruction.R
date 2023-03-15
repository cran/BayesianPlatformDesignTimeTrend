#' @title Boundaryconstruction
#' @description This function constructs the stopping boundary based on input information

#' @param Stopbound.inf The list of stop boundary information for more see \code{\link{Stopboundinf}}
#' @param ns  A vector of accumulated number of patient at each stage
#'
#' @return A list of the futility boundary and the efficacy boundary
#' @importFrom stats pnorm
#' @importFrom stats var
#' @export
#'
#' @examples
#'  Stopbound.inf=list(Stop.type="Early-Pocock",Boundary.type="Symmetric",cutoff=c(0.9928,0.0072))
#'  ns=c(60,120,180,240,300)
#'  Boundaryconstruction(Stopbound.inf, ns)
#' @author Ziyan Wang
Boundaryconstruction = function(Stopbound.inf = Stopbound.inf, ns = ns) {
  cutoff.temp = Stopbound.inf$cutoff #(cutoff1, cutoff2): cutoff1 is efficacy cutoff; cutoff2 is futility cutoff
  if (length(cutoff.temp) != 2 | sum(is.na(cutoff.temp)) != 0) {
    stop(
      "Error: Please input the cutoff value as a vector with two values for both symmetric and asymmetric boundary"
    )
  }
  boundary = Stopbound.inf$Boundary.type
  #------------Need to Check--------------
  stage = length(ns)
  #---------------------------------------
  Stoptype = Stopbound.inf$Stop.type
  if (boundary == "Symmetric") {
    if (Stoptype == "Noearly") {
      cutoffeff = c(rep(1.01, stage - 1), cutoff.temp[1]) #Ensure no stop during the trial (The efficacy cutoff before the final stage > 1)
      cutoffful = c(rep(-0.01, stage - 1), cutoff.temp[2])  #Ensure no stop during the trial (The futility cutoff before the final stage < 0)
      if (sum(cutoff.temp) != 1) {
        stop(
          "Error: The cutoff inputted is Asymmetric. Please revise the Boundary.type in Stopbound.inf list. The sum of them for symmetric Noearly should be 1"
        )
      }
      if ((cutoff.temp[1] > 1) > 0 |
          (cutoff.temp[1] < cutoff.temp[2]) > 0 | cutoff.temp[2] < 0) {
        stop(
          "Error: The fultility cutoff should be smaller than efficacy at the final stage. Also, both cutoffs should be <= 1 except for OBF boundary"
        )
      }
    }
    else if (Stoptype == "Early-Pocock") {
      cutoffeff = rep(cutoff.temp[1], stage)
      cutoffful = rep(cutoff.temp[2], stage)
      if (sum(cutoff.temp) != 1) {
        stop(
          "Error: The cutoff inputted is Asymmetric. Please revise the Boundary.type in Stopbound.inf list. The sum of them for symmetric Noearly should be 1"
        )
      }
      if (sum(cutoffeff > 1) > 0 |
          sum(cutoffeff < cutoffful) > 0 | cutoff.temp[2] < 0) {
        stop(
          "Error: The cutoff of Pocock boundary should be smaller or equal to 1 (<= 1), The fultility cutoff should be smaller than efficacy at each stage"
        )
      }
    }
    else if (Stoptype == "Early-OBF") {
      cutoffeff = pnorm(sqrt(stage / seq(1:stage) * cutoff.temp[1]))
      cutoffful = pnorm(-sqrt(stage / seq(1:stage) * cutoff.temp[2]))
      if (var(cutoff.temp) != 0) {
        stop(
          "Error: For symmetric OBF boundary, the cutoff for efficacy and fultility should be the same. The function will atuomatically construct an symmetric boundary set. cutoff can be greater than 1"
        )
      }
    }
    else {
      stop(
        "Error: Please input the stopping boundary type correctly should be in c(Noearly,Early-Pocock,Early-OBF)"
      )
    }
  }
  else if (boundary == "Asymmetric") {
    if (Stoptype == "Noearly") {
      cutoffeff = c(rep(1.01, stage - 1), cutoff.temp[1]) #Ensure no stop during the trial (The efficacy cutoff before the final stage > 1)
      cutoffful = c(rep(-0.01, stage - 1), cutoff.temp[2])  #Ensure no stop during the trial (The futility cutoff before the final stage < 0)
      if (sum(cutoff.temp) == 1) {
        stop(
          "Error: The cutoff inputted is symmetric. Please revise the Boundary.type in Stopbound.inf list."
        )
      }
      if (cutoff.temp[1] > 1 |
          (cutoff.temp[1] < cutoff.temp[2]) > 0 | cutoff.temp[2] < 0) {
        stop(
          "Error: The fultility cutoff should be smaller than efficacy at the final stage. Also, both cutoffs should be <= 1 except for OBF boundary"
        )
      }
    }
    else if (Stoptype == "Early-Pocock") {
      cutoffeff = rep(cutoff.temp[1], stage)
      cutoffful = rep(cutoff.temp[2], stage)
      if (sum(cutoff.temp) == 1) {
        stop(
          "Error: The cutoff inputted is symmetric. Please revise the Boundary.type in Stopbound.inf list."
        )
      }
      if (sum(cutoffeff > 1) > 0 |
          sum(cutoffeff < cutoffful) > 0 | cutoff.temp[2] < 0) {
        stop(
          "Error: The cutoff of Pocock boundary should be smaller or equal to 1 (<= 1), The fultility cutoff should be smaller than efficacy at each stage"
        )
      }
    }
    else if (Stoptype == "Early-OBF") {
      cutoffeff = pnorm(sqrt(stage / seq(1:stage) * cutoff.temp[1]))
      cutoffful = pnorm(-sqrt(stage / seq(1:stage) * cutoff.temp[2]))
      if (var(cutoff.temp) == 0) {
        stop(
          "Error: For Asymmetric OBF boundary, the cutoff for efficacy and fultility should be the different. The function will atuomatically construct a asymmetric boundary set. cutoff can be greater than 1"
        )
      }
    }
    else {
      stop(
        "Error: Please input the stopping boundary type correctly should be in c(Noearly,Early-Pocock,Early-OBF)"
      )
    }
  }
  else {
    stop("Error: The boundary type is invalid")
  }
  return(list(
    Efficacy.boundary = cutoffeff,
    Fultility.boundary = cutoffful
  ))
}
#---------------------------------------------------------------------
