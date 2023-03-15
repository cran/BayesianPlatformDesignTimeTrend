#' @title perHtypeIerror_powerfunc
#' @description This function reads in the output matrix of a number of trial replicates to calculate the error rate
#' @param res A list of output matrix of a number of trial replicates
#'
#' @return Type I error rate or power of each treatment - control comparison
#' @export
#'
#' @examples
#' \dontrun{perHtypeIerror_powerfunc(res)}
#' @author Ziyan Wang
perHtypeIerror_powerfunc = function(res) {
  colMeans(t(sapply(res, function(x) {
    resname = colnames(x)
    K = sum(stringr::str_detect(colnames(x), "H")) + 1
    #Indentify which hypothesis is rejected
    reject = which(matrix(x[, (K - 1 + 2 * K + 1):(K - 1 + 2 * K + K - 1)] %in% 1, ncol =
                            K - 1), arr.ind = TRUE)[, 2]
    if (length(reject) >= 1) {
      rejectres = rep(0, K - 1)
      rejectres[reject] = 1
      return(rejectres)
    }
    else{
      return(rep(0, K - 1))
    }
  })))
}
#' @title FWER_disconjunctivepowerfunc
#' @description This function reads in the output matrix of a number of trial replicates to calculate the Family wise error rate or disconjunctive power
#' @param res A list of output matrix of a number of trial replicates
#'
#' @return Family wise error rate or disconjunctive power
#' @export
#'
#' @examples
#' \dontrun{FWER_disconjunctivepowerfunc(res)}
FWER_disconjunctivepowerfunc = function(res) {
  mean(sapply(res, function(x) {
    resname = colnames(x)
    K = sum(stringr::str_detect(colnames(x), "H")) + 1
    if (sum(x[, (K - 1 + 2 * K + 1):(K - 1 + 2 * K + K - 1)] %in% 1) >= 1) {
      return(1)
    }
    else{
      return(0)
    }
  }))
}

#' @title Meanfunc
#' @description This function reads in the output matrix of a number of trial replicates to calculate mean treatment effect estimate.
#' @param res A list of output matrix of a number of trial replicates
#'
#' @return Mean treatment effect estimates of each treatment arm
#' @export
#'
#' @examples
#' \dontrun{Meanfunc(res)}
Meanfunc = function(res) {
  K = mean(sapply(res, function(x) {
    K = sum(stringr::str_detect(colnames(x), "H")) + 1
    return(K)
  }))
  meaneffect = colMeans(matrix(t(sapply(res, function(x) {
    stage = dim(x)[1]
    resname = colnames(x)
    K = sum(stringr::str_detect(colnames(x), "H")) + 1
    reject = which(matrix(x[, (K - 1 + 2 * K + 1):(K - 1 + 2 * K + K - 1)] %in% 1, ncol =
                            K - 1), arr.ind = TRUE)[,2]
    if (length(reject) >= 1) {
      drop.at = which(matrix(x[, (K - 1 + 2 * K + 1):(K - 1 + 2 * K + K - 1)] %in% 1, ncol =
                               K - 1), arr.ind = TRUE)[, 1]
      drop.at.all = rep(stage, K - 1)
      drop.at.all[reject] = drop.at
      treatmentindex = seq(1, K - 1)
      trtmean.loc = cbind(drop.at.all, treatmentindex)
      meanres = matrix(x[, (K - 1 + 2 * K + K - 1 + 1 + 1):(K - 1 + 2 *
                                                              K + K - 1 + 1 + K - 1)], ncol = K - 1)
      result = rep(NA, K - 1)
      for (i in 1:(K - 1)) {
        result[i] = meanres[trtmean.loc[i, 1], trtmean.loc[i, 2]]
      }
      return(result)
    }
    else{
      drop.at.all = rep(stage, K - 1)
      treatmentindex = seq(1, K - 1)
      trtmean.loc = cbind(drop.at.all, treatmentindex)
      meanres = matrix(x[, (K - 1 + 2 * K + K - 1 + 1 + 1):(K - 1 + 2 *
                                                              K + K - 1 + 1 + K - 1)], ncol = K - 1)
      result = rep(NA, K - 1)
      for (i in 1:(K - 1)) {
        result[i] = meanres[trtmean.loc[i, 1], trtmean.loc[i, 2]]
      }
      return(result)
    }
  })), ncol = K - 1))
  return(meaneffect)
}
#' @title varfunc
#' @description This function reads in the output matrix of a number of trial replicates to calculate the variance of treatment effect estimate.
#' @param res A list of output matrix of a number of trial replicates
#'
#' @return The variance of Treatment effect estimates of each treatment arm
#' @export
#'
#' @examples
#' \dontrun{varfunc(res)}
varfunc = function(res) {
  K = mean(sapply(res, function(x) {
    K = sum(stringr::str_detect(colnames(x), "H")) + 1
    return(K)
  }))
  meaneffect = matrixStats::colVars(matrix(t(sapply(res, function(x) {
    stage = dim(x)[1]
    resname = colnames(x)
    K = sum(stringr::str_detect(colnames(x), "H")) + 1
    reject = which(matrix(x[, (K - 1 + 2 * K + 1):(K - 1 + 2 * K + K - 1)] %in% 1, ncol =
                            K - 1), arr.ind = TRUE)[,2]
    if (length(reject) >= 1) {
      drop.at = which(matrix(x[, (K - 1 + 2 * K + 1):(K - 1 + 2 * K + K - 1)] %in% 1, ncol =
                               K - 1), arr.ind = TRUE)[, 1]
      drop.at.all = rep(stage, K - 1)
      drop.at.all[reject] = drop.at
      treatmentindex = seq(1, K - 1)
      trtmean.loc = cbind(drop.at.all, treatmentindex)
      meanres = matrix(x[, (K - 1 + 2 * K + K - 1 + 1 + 1):(K - 1 + 2 *
                                                              K + K - 1 + 1 + K - 1)], ncol = K - 1)
      result = rep(NA, K - 1)
      for (i in 1:(K - 1)) {
        result[i] = meanres[trtmean.loc[i, 1], trtmean.loc[i, 2]]
      }
      return(result)
    }
    else{
      drop.at.all = rep(stage, K - 1)
      treatmentindex = seq(1, K - 1)
      trtmean.loc = cbind(drop.at.all, treatmentindex)
      meanres = matrix(x[, (K - 1 + 2 * K + K - 1 + 1 + 1):(K - 1 + 2 *
                                                              K + K - 1 + 1 + K - 1)], ncol = K - 1)
      result = rep(NA, K - 1)
      for (i in 1:(K - 1)) {
        result[i] = meanres[trtmean.loc[i, 1], trtmean.loc[i, 2]]
      }
      return(result)
    }
  })), ncol = K - 1))
  return(meaneffect)
}

#' @title Nfunc
#' @description This function reads in the output matrix of a number of trial replicates to calculate mean estimate of total number of patients allocated to each arm
#' @param res A list of output matrix of a number of trial replicates
#'
#' @return The mean estimate of total number of patients allocated to each arm
#' @export
#'
#' @examples
#' \dontrun{Nfunc(res)}
Nfunc = function(res) {
  K = mean(sapply(res, function(x) {
    K = sum(stringr::str_detect(colnames(x), "H")) + 1
    return(K)
  }))
  Nmean = colMeans(matrix(t(sapply(res, function(x) {
    stage = dim(x)[1]
    resname = colnames(x)
    K = sum(stringr::str_detect(colnames(x), "H")) + 1
    reject = which(matrix(x[, (K - 1 + 2 * K + 1):(K - 1 + 2 * K + K - 1)] %in% 1, ncol =
                            K - 1), arr.ind = TRUE)[,2]
    if (length(reject) >= 1) {
      drop.at = which(matrix(x[, (K - 1 + 2 * K + 1):(K - 1 + 2 * K + K - 1)] %in% 1, ncol =
                               K - 1), arr.ind = TRUE)[, 1]
      drop.at.all = rep(stage, K - 1)
      drop.at.all[reject] = drop.at
      treatmentindex = seq(1, K - 1)
      trtmean.loc = cbind(c(max(drop.at.all), drop.at.all), c(1, treatmentindex +
                                                                1))
      Nres = matrix(x[, seq(K, K - 1 + 2 * K - 1, 2)], ncol = K)
      result = rep(NA, K)
      for (i in 1:K) {
        result[i] = Nres[trtmean.loc[i, 1], trtmean.loc[i, 2]]
      }
      return(result)
    }
    else{
      drop.at.all = rep(stage, K)
      treatmentindex = seq(1, K - 1)
      trtmean.loc = cbind(drop.at.all, c(1, treatmentindex + 1))
      Nres = matrix(x[, seq(K, K - 1 + 2 * K - 1, 2)], ncol = K)
      result = rep(NA, K)
      for (i in 1:K) {
        result[i] = Nres[trtmean.loc[i, 1], trtmean.loc[i, 2]]
      }
      return(result)
    }
  })), ncol = K))
  return(Nmean)
}

#' @title Sperarmfunc
#' @description This function reads in the output matrix of a number of trial replicates to calculate mean total number of survived patients of each arm
#' @param res A list of output matrix of a number of trial replicates
#'
#' @return The mean total number of survived patients of each arm
#' @export
#'
#' @examples
#' \dontrun{Sperarmfunc(res)}
Sperarmfunc = function(res) {
  K = mean(sapply(res, function(x) {
    K = sum(stringr::str_detect(colnames(x), "H")) + 1
    return(K)
  }))
  Smean = colMeans(matrix(t(sapply(res, function(x) {
    stage = dim(x)[1]
    resname = colnames(x)
    K = sum(stringr::str_detect(colnames(x), "H")) + 1
    reject = which(matrix(x[, (K - 1 + 2 * K + 1):(K - 1 + 2 * K + K - 1)] %in% 1, ncol =
                            K - 1), arr.ind = TRUE)[,2]
    if (length(reject) >= 1) {
      drop.at = which(matrix(x[, (K - 1 + 2 * K + 1):(K - 1 + 2 * K + K - 1)] %in% 1, ncol =
                               K - 1), arr.ind = TRUE)[, 1]
      drop.at.all = rep(stage, K - 1)
      drop.at.all[reject] = drop.at
      treatmentindex = seq(1, K - 1)
      trtmean.loc = cbind(c(max(drop.at.all), drop.at.all), c(1, treatmentindex +
                                                                1))
      Nres = matrix(x[, seq(K + 1, K - 1 + 2 * K, 2)], ncol = K)
      result = rep(NA, K)
      for (i in 1:K) {
        result[i] = Nres[trtmean.loc[i, 1], trtmean.loc[i, 2]]
      }
      return(result)
    }
    else{
      drop.at.all = rep(stage, K)
      treatmentindex = seq(1, K - 1)
      trtmean.loc = cbind(drop.at.all, c(1, treatmentindex + 1))
      Nres = matrix(x[, seq(K + 1, K - 1 + 2 * K, 2)], ncol = K)
      result = rep(NA, K)
      for (i in 1:K) {
        result[i] = Nres[trtmean.loc[i, 1], trtmean.loc[i, 2]]
      }
      return(result)
    }
  })), ncol = K))
  return(Smean)
}

# list.of.analysisfunction <-
#   list(
#     perHtypeIerrorfunc = perHtypeIerrorfunc,
#     FWERfunc = FWERfunc,
#     Meanfunc = Meanfunc,
#     varfunc = varfunc,
#     Nfunc = Nfunc,
#     Sperarmfunc = Sperarmfunc
#   )
