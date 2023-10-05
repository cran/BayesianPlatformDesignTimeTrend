#' @title Randomisation.inf
#' @description This function checks the validity of the randomisation information input

#' @param Random.inf A list of adaptive randomisation information.
#'     'Fixratio' a indicator of whether the randomisation process uses fix ratio. Default is FALSE.
#'     'Fixratiocontrol' the numerical value indicating the randomisation weight of the control arm compared to the treatment arms. Default is NA for Fixratio = FALSE.
#'     'BARmethod' the Bayesian adaptive randomisation type. Default is "Thall" indicating the use of Thall's approach in the randomisation process. The other value is 'Trippa'.
#'     'Thall.tuning.inf' the list of tuning parameter for Thall's approach including 'tuningparameter' (Default is "Fixed" indicating that the tuning paramter is fixed for all stages) and 'fixvalue' (Default is 1).
#' @return A list of input randomisation information
#' @export
#'
#' @examples
#' Randomisation.inf(Random.inf = list(
#' Fixratio = FALSE,
#' Fixratiocontrol = NA,
#' BARmethod = "Thall",
#' Thall.tuning.inf = list(tuningparameter = "Fixed",  fixvalue = 1),
#' Trippa.tuning.inf = list(a = NA, b = NA)
#' ))
#' @author Ziyan Wang
Randomisation.inf = function(Random.inf = list(
  Fixratio = FALSE,
  Fixratiocontrol = NA,
  BARmethod = "Thall",
  Thall.tuning.inf = list(tuningparameter = "Fixed",  fixvalue = 1),
  Trippa.tuning.inf = list(a = 10, b = 0.75)
)) {
  Fixratio = Random.inf$Fixratio
  if (Fixratio == T) {
    if (is.na(Random.inf$Fixratiocontrol) |
        Random.inf$Fixratiocontrol <= 0) {
      stop(
        "Error: The value R > 0 for fix randomisation (R:1:1:1:......) should be input which is Fixratiocontrol"
      )
    }
    else{
      Fixratiocontrol = Random.inf$Fixratiocontrol
    }
    return(
      list(
        Fixratio = Fixratio,
        Fixratiocontrol = Fixratiocontrol,
        BARmethod = NA,
        Thall.tuning.inf = list(tuningparameter = NA,  c = NA),
        Trippa.tuning.inf = list(a = NA, b = NA)
      )
    )
  }
  else {
    BARmethod = Random.inf$BARmethod
    if (BARmethod == "Thall") {
      tuning.inf = Random.inf$Thall.tuning.inf$tuningparameter
      if (tuning.inf == "Fixed") {
        tuningparameter = "Fixed"
        c = Random.inf$Thall.tuning.inf$fixvalue
        if (is.na(Random.inf$Thall.tuning.inf$fixvalue)) {
          stop("Error: The value of tuning parameter in Thall's approach should be specified.")
        }
        return(
          list(
            BARmethod = BARmethod,
            Thall.tuning.inf = list(tuningparameter = tuningparameter, c = c),
            Fixratio = Fixratio,
            Fixratiocontrol = NA
          )
        )
      }
      else {
        tuningparameter = "Unfixed"
      }
      return(
        list(
          BARmethod = BARmethod,
          Thall.tuning.inf = list(tuningparameter = tuningparameter,  c = NA),
          Fixratio = Fixratio,
          Fixratiocontrol = NA
        )
      )
    }
    else{
      BARmethod = "Trippa"
      a = Random.inf$Trippa.tuning.inf$a
      b = Random.inf$Trippa.tuning.inf$b
      return(
        list(
          BARmethod = BARmethod,
          Trippa.tuning.inf = list(a = a, b = b),
          Thall.tuning.inf = list(tuningparameter = NA,  c = NA),
          Fixratio = Fixratio,
          Fixratiocontrol = NA
        )
      )
    }
  }
}
