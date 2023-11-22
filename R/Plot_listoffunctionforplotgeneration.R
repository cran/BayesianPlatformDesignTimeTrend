#' @title alphaspending
#' @description This function estimates the mean error rate spent at each interim analysis for a trial
#' Example usage:
#' 1. sapply(res = result, fun = alphaspending) will generate list of the proportion of trial replicates are stopped at each stage for all scenarios in result where result is a list containing output data for different scenario
#' 2. sapply(sapply(result,FUN = alphaspending),sum) will generate the type I error rate or power for all scenario on the result list
#' 3. alpha(result) generate the proportion of trial replicates are stopped at each stage where result is the output data for one specific scenario
#' 4. sum(alpha(result)) will generate the type I error rate or power for a specific scenario
#' @param res A list of output matrix of a number of trial replicates
#'
#' @return The error rate at each interim analysis
#' @export
#'
#' @examples
#' \dontrun{alphaspending(res)}
#' @author Ziyan Wang
alphaspending = function(res) {
  K = mean(sapply(res, function(x) {
    K = sum(stringr::str_detect(colnames(x), "H")) + 1
    return(K)
  }))
  ntrials = length(res)
  round(colSums(t(sapply(res, function(x) {
    rejectres = matrix(x[, (K - 1 + 2 * K + 1):(K - 1 + 2 * K + K - 1)], ncol =
                         K - 1)
    return(rejectres)
  })), na.rm = T) / ntrials, 4)
}
# Example:
# 1. sapply(res = result, fun = alphaspending) will generate list of the proportion of trial replicates are stopped at each stage for all scenarios in result
# where result is a list containing output data for different scenario
# 2. sapply(sapply(result,FUN = alphaspending),sum) will generate the type I error rate or power for all scenario on the result list
# 3. alpha(result) generate the proportion of trial replicates are stopped at each stage
# where result is the output data for one specific scenario
# 4. sum(alpha(result)) will generate the type I error rate or power for a specific scenario

#' @title trtbias
#' @description This function estimates the mean bias of treatment effect
#' @param res A list of output matrix of a number of trial replicates
#'
#' @param trueeffect A vector of true treatment effect in each scenario
#'
#' @return A matrix of mean treatment effect bias
#' @export
#'
#' @examples
#' \dontrun{trtbias(res, trueeffect)}
trtbias = function(res, trueeffect) {
  namedata = names(res)
  datatempmeanout = {

  }
  count = 1
  for (temp in 1:length(res)) {
    K = mean(sapply(res[[temp]], function(x) {
      K = sum(stringr::str_detect(colnames(x), "H")) + 1
      return(K)
    }))
    ntrials = length(res[[temp]])
    datatempmean = matrix(rep(NA, ntrials * (K - 1)), ncol = (K - 1))
    datatempmean = matrix(t(sapply(res[[temp]], function(x) {
      stage = dim(x)[1]
      resname = colnames(x)
      K = sum(stringr::str_detect(colnames(x), "H")) + 1
      reject = which(matrix(x[, (K - 1 + 2 * K + 1):(K - 1 + 2 * K + K -
                                                       1)] %in% 1, ncol = K - 1), arr.ind = TRUE)[, 2]
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
        for (i in 1:K - 1) {
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
        for (i in 1:K - 1) {
          result[i] = meanres[trtmean.loc[i, 1], trtmean.loc[i, 2]]
        }
        return(result)
      }
    })), ncol = K - 1)
    datatempmeanout = cbind(datatempmeanout, datatempmean - trueeffect[count])
    count = count + 1
  }
  tempname = rep(namedata, each = K - 1)
  colnames(datatempmeanout) = tempname
  # colnames(datatempmeanout) = paste0(tempname,rep(c(1,2),length(res)))
  datatemp = reshape::melt(datatempmeanout)
  names(datatemp)[names(datatemp) == "value"] = "Treatmenteffect"
  names(datatemp)[names(datatemp) == "X2"] = "Model"
  return(datatemp)
}

# intdataout = function(res) {
#   namedata = names(res)
#   datatempmean = {
#
#   }
#   for (temp in 1:length(res)) {
#     K = mean(sapply(res[[temp]], function(x) {
#       K = sum(stringr::str_detect(colnames(x), "H")) + 1
#       return(K)
#     }))
#     datatempmean = cbind(datatempmean, matrix(t(sapply(res[[temp]], function(x) {
#       stage = dim(x)[1]
#       resname = colnames(x)
#       K = sum(stringr::str_detect(colnames(x), "H")) + 1
#       reject = which(matrix(x[, (K - 1 + 2 * K + 1):(K - 1 + 2 * K + K -
#                                                        1)] %in% 1, ncol = K - 1), arr.ind = TRUE)[, 2]
#       if (length(reject) >= 1) {
#         drop.at = which(matrix(x[, (K - 1 + 2 * K + 1):(K - 1 + 2 * K + K - 1)] %in% 1, ncol =
#                                  K - 1), arr.ind = TRUE)[, 1]
#         drop.at.all = rep(stage, K - 1)
#         drop.at.all[reject] = drop.at
#         treatmentindex = seq(1, K - 1)
#         trtmean.loc = cbind(drop.at.all, treatmentindex)
#         meanres = matrix(x[, (K - 1 + 2 * K + K - 1 + 1 + K - 1 + K - 1 +
#                                 (stage - 1) + 1):(K - 1 + 2 * K + K - 1 + 1 + K - 1 + K - 1 + (stage - 1) +
#                                                     (stage - 1) * (K - 1))], ncol = (stage - 1) * (K - 1))
#         result = rep(NA, K - 1)
#         for (i in 1:(K - 1)) {
#           if (trtmean.loc[i, 1] == 1) {
#             result[i] = NA
#           }
#           else{
#             result[i] = meanres[trtmean.loc[i, 1], trtmean.loc[i, 2] * (trtmean.loc[i, 1] -
#                                                                           1)]
#           }
#         }
#         return(c(result, drop.at.all))
#       }
#       else{
#         drop.at.all = rep(stage, K - 1)
#         treatmentindex = seq(1, K - 1)
#         trtmean.loc = cbind(drop.at.all, treatmentindex)
#         meanres = matrix(x[, (K - 1 + 2 * K + K - 1 + 1 + K - 1 + K - 1 +
#                                 (stage - 1) + 1):(K - 1 + 2 * K + K - 1 + 1 + K - 1 + K - 1 + (stage - 1) +
#                                                     (stage - 1) * (K - 1))], ncol = (stage - 1) * (K - 1))
#         result = rep(NA, K - 1)
#         for (i in 1:K - 1) {
#           result[i] = meanres[trtmean.loc[i, 1], trtmean.loc[i, 2] * (stage - 1)]
#         }
#         return(c(result, drop.at.all))
#       }
#     })), ncol = 2 * (K - 1)))
#   }
#   dataname = {
#
#   }
#   for (nameind in 1:length(namedata)) {
#     dataname = c(dataname, namedata[nameind], paste0(namedata[nameind], "_stage"))
#   }
#   colnames(datatempmean) = dataname
#   # datatemp=reshape::melt(datatempmean)
#   # names(datatemp)[names(datatemp)=="value"]="interactioneffect"
#   # names(datatemp)[names(datatemp)=="X2"]="Model"
#   return(datatempmean)
# }

#' @title intbias
#' @description This function estimates the mean bias of treatment - stage interaction effect
#' @param res A list of output matrix of a number of trial replicates
#'
#' @return A matrix of mean treatment - stage interaction effect bias
#' @export
#'
#' @examples
#' \dontrun{intbias(res)}
intbias = function(res) {
  namedata = names(res)
  datatempmean = {

  }
  for (temp in 1:length(res)) {
    K = mean(sapply(res[[temp]], function(x) {
      K = sum(stringr::str_detect(colnames(x), "H")) + 1
      return(K)
    }))
    datatempmean = cbind(datatempmean, matrix(t(sapply(res[[temp]], function(x) {
      stage = dim(x)[1]
      resname = colnames(x)
      K = sum(stringr::str_detect(colnames(x), "H")) + 1
      reject = which(matrix(x[, (K - 1 + 2 * K + 1):(K - 1 + 2 * K + K -
                                                       1)] %in% 1, ncol = K - 1), arr.ind = TRUE)[, 2]
      if (length(reject) >= 1) {
        drop.at = which(matrix(x[, (K - 1 + 2 * K + 1):(K - 1 + 2 * K + K - 1)] %in% 1, ncol =
                                 K - 1), arr.ind = TRUE)[, 1]
        drop.at.all = rep(stage, K - 1)
        drop.at.all[reject] = drop.at
        treatmentindex = seq(1, K - 1)
        trtmean.loc = cbind(drop.at.all, treatmentindex)
        meanres = matrix(x[, (K - 1 + 2 * K + K - 1 + 1 + K - 1 + K - 1 +
                                (stage - 1) + 1):(K - 1 + 2 * K + K - 1 + 1 + K - 1 + K - 1 + (stage - 1) +
                                                    (stage - 1) * (K - 1))], ncol = (stage - 1) * (K - 1))
        result = rep(NA, K - 1)
        for (i in 1:(K - 1)) {
          if (trtmean.loc[i, 1] == 1) {
            result[i] = NA
          }
          else{
            result[i] = meanres[trtmean.loc[i, 1], trtmean.loc[i, 2] * (trtmean.loc[i, 1] -
                                                                          1)]
          }
        }
        return(result)
      }
      else{
        drop.at.all = rep(stage, K - 1)
        treatmentindex = seq(1, K - 1)
        trtmean.loc = cbind(drop.at.all, treatmentindex)
        meanres = matrix(x[, (K - 1 + 2 * K + K - 1 + 1 + K - 1 + K - 1 +
                                (stage - 1) + 1):(K - 1 + 2 * K + K - 1 + 1 + K - 1 + K - 1 + (stage - 1) +
                                                    (stage - 1) * (K - 1))], ncol = (stage - 1) * (K - 1))
        result = rep(NA, K - 1)
        for (i in 1:K - 1) {
          result[i] = meanres[trtmean.loc[i, 1], trtmean.loc[i, 2] * (stage - 1)]
        }
        return(result)
      }
    })), ncol = K - 1))
  }
  colnames(datatempmean) = namedata
  datatemp = reshape::melt(datatempmean)
  names(datatemp)[names(datatemp) == "value"] = "interactioneffect"
  names(datatemp)[names(datatemp) == "X2"] = "Model"
  return(datatemp)
}

#' trteffect
#' @description This function estimates the mean treatment effect bias and its rooted mean squared error
#' @param res A list of output matrix of a number of trial replicates
#' @param trueeff A vector of true treatment effect in each scenario
#'
#' @return A vector of mean treatment effect bias and its rooted mean squared error
#' @export
#'
#' @examples
#' \dontrun{trteffect(res, trueeff)}
trteffect = function(res, trueeff) {
  K = sapply(res, function(x) {
    K = sum(stringr::str_detect(colnames(x), "H")) + 1
    return(K)
  })
  ntrials = length(res)
  samp = t(sapply(res, function(x) {
    stage = dim(x)[1]
    resname = colnames(x)
    K = sum(stringr::str_detect(colnames(x), "H")) + 1
    reject = which(matrix(x[, (K - 1 + 2 * K + 1):(K - 1 + 2 * K + K - 1)] %in% 1, ncol =
                            K - 1), arr.ind = TRUE)[, 2]

    meanres = matrix(x[, (K - 1 + 2 * K + K - 1 + 1 + 1):(K - 1 + 2 * K +
                                                            K - 1 + 1 + K - 1)], ncol = K - 1)
  })) - trueeff


  meaneffect = colMeans(samp, na.rm = T)
  sdeffect = apply(samp, 2, stats::sd, na.rm = T)

  replicate = {

  }
  for (i in 1:dim(samp)[2]) {
    replicate[i] = ntrials - sum(is.na(samp[, i]))
  }
  seeffect = sdeffect / sqrt(replicate)
  return(rbind(meaneffect, seeffect))
}

# list.of.Plotfunction <-
#   list(
#     alphaspending = alphaspending,
#     trtbias = trtbias,
#     intdataout = intdataout,
#     intbias = intbias,
#     trteffect = trteffect,
#     Nfunc = Nfunc,
#     Sperarmfunc = Sperarmfunc
#   )
