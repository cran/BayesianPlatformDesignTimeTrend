## ----include = FALSE----------------------------------------------------------
 knitr::opts_chunk$set(
   collapse = TRUE,
   comment = "#>"
 )

## ----setup--------------------------------------------------------------------
library(BayesianPlatformDesignTimeTrend)

## -----------------------------------------------------------------------------

ntrials = 1000 # Number of trial replicates
ns = seq(120,600,120) # Sequence of total number of accrued patients at each interim analysis
null.reponse.prob = 0.4
alt.response.prob = 0.6

# We investigate the type I error rate for different time trend strength
null.scenario = matrix(
  c(
    null.reponse.prob,
    null.reponse.prob,
    null.reponse.prob,
    null.reponse.prob
  ),
  nrow = 1,
  ncol = 4,
  byrow = T
)
# alt.scenario = matrix(c(null.reponse.prob,null.reponse.prob,null.reponse.prob,null.reponse.prob,
#                     null.reponse.prob,alt.response.prob,null.reponse.prob,null.reponse.prob,
#                     null.reponse.prob,alt.response.prob,alt.response.prob,null.reponse.prob,
#                     null.reponse.prob,alt.response.prob,alt.response.prob,alt.response.prob), nrow=3, ncol = 4,byrow=T)
model = "tlr" #logistic model
max.ar = 0.75  #limit the allocation ratio for the control group (1-max.ar < r_control < max.ar)
#------------Select the data generation randomisation methods-------
rand.type = "Urn" # Urn design
max.deviation = 3 # The recommended value for the tuning parameter in the Urn design

# Require multiple cores for parallel running
cl = 2

# Set the model we want to use and the time trend effect for each model used.
# Here the main model will be used twice for two different strength of time trend c(0,0,0,0) and c(1,1,1,1) to investigate how time trend affect the evaluation metrics in BAR setting.
# Then the main + stage_continuous model which is the treatment effect + stage effect model will be applied for strength equal c(1,1,1,1) to investigate how the main + stage effect model improve the evaluation metrics.
reg.inf = c("main", "main", "main + stage_continuous")
trend.effect = matrix(
  c(0, 0, 0, 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
  ncol = 4,
  nrow = 3,
  byrow = T
)

#
cutoffearly = matrix(rep(0.994, dim(null.scenario)[1]), ncol = 1)

K = dim(null.scenario)[2]
print(
  paste0(
    "Start trial simulation. This is a ",
    K,
    "-arm trial simulation. There are one null scenario and ",
    K - 1 ,
    " alternative scenarios. There are ",
    K ,
    " rounds."
  )
)
cutoffindex = 1

## ----eval=FALSE---------------------------------------------------------------
#  result = {
#  
#  }
#  OPC_null = {
#  
#  }
#  for (i in 1:dim(null.scenario)[1]) {
#    trendindex = 1
#    for (j in 1:length(reg.inf)){
#      restlr = Trial.simulation(
#        ntrials = ntrials,
#        # Number of trial replicates
#        trial.fun = simulatetrial,
#        # Call the main function
#        input.info = list(
#          response.probs = null.scenario[cutoffindex, ],
#          #The scenario vector in this round
#          ns = ns,
#          # Sequence of total number of accrued patients at each interim analysis
#          max.ar =  max.ar,
#          #limit the allocation ratio for the control group (1-max.ar < r_control < max.ar)
#          rand.type = rand.type,
#          # Which randomisation methods in data generation.
#          max.deviation = max.deviation,
#          # The recommended value for the tuning parameter in the Urn design
#          model.inf = list(
#            model = model,
#            #Use which model?
#            ibb.inf = list(
#              #independent beta-binomial model which can be used only for no time trend simulation
#              pi.star = 0.5,
#              # beta prior mean
#              pess = 2,
#              # beta prior effective sample size
#              betabinomialmodel = ibetabinomial.post # beta-binomial model for posterior estimation
#            ),
#            tlr.inf = list(
#              beta0_prior_mu = 0,
#              # Stan logistic model t prior location
#              beta1_prior_mu = 0,
#              # Stan logistic model t prior location
#              beta0_prior_sigma = 2.5,
#              # Stan logistic model t prior sigma
#              beta1_prior_sigma = 2.5,
#              # Stan logistic model t prior sigma
#              beta0_df = 7,
#              # Stan logistic model t prior degree of freedom
#              beta1_df = 7,
#              # Stan logistic model t prior degree of freedom
#              reg.inf =  reg.inf[trendindex],
#              # The model we want to use
#              variable.inf = "Fixeffect" # Use fix effect logistic model
#            )
#          ),
#          Stopbound.inf = Stopboundinf(
#            Stop.type = "Early-Pocock",
#            # Use Pocock like early stopping boundary
#            Boundary.type = "Symmetric",
#            # Use Symmetric boundary where cutoff value for efficacy boundary and futility boundary sum up to 1
#            cutoff = c(cutoffearly[cutoffindex, 1], 1 - cutoffearly[cutoffindex, 1]) # The cutoff value for stopping boundary
#          ),
#          Random.inf = list(
#            Fixratio = FALSE,
#            # Do not use fix ratio allocation
#            Fixratiocontrol = NA,
#            # Do not use fix ratio allocation
#            BARmethod = "Thall",
#            # Use Thall's Bayesian adaptive randomisation approach
#            Thall.tuning.inf = list(tuningparameter = "Fixed",  fixvalue = 1) # Specified the tunning parameter value for fixed tuning parameter
#          ),
#          trend.inf = list(
#            trend.type = "step",
#            # Linear time trend pattern
#            trend.effect = trend.effect[trendindex, ],
#            # Stength of time trend effect
#            trend_add_or_multip = "mult" # Multiplicative time trend effect on response probability
#          )
#        ),
#        cl = 2 # 2 cores required
#      )
#  
#    trendindex = trendindex + 1
#    # The result list can be used for plotting and the OPC table is the summary evaluaton metrics for each scenario
#    result = c(result, restlr$result)
#    OPC_null = rbind(OPC_null, restlr$OPC)
#    }
#    cutoffindex = cutoffindex + 1
#  }

## -----------------------------------------------------------------------------
print("Finished null scenario study")
save_data = FALSE
if (isTRUE(save_data)) {
  save(result, file = restlr$Nameofsaveddata$nameData)
  save(OPC_null, file = restlr$Nameofsaveddata$nameTable)
}


## -----------------------------------------------------------------------------
# Characteristic table
print(OPC_null)

## -----------------------------------------------------------------------------

ntrials = 1000 # Number of trial replicates
ns = seq(120,600,120) # Sequence of total number of accrued patients at each interim analysis
null.reponse.prob = 0.4
alt.response.prob = 0.6

# We investigate the type I error rate for different time trend strength
alt.scenario = matrix(
  c(
    null.reponse.prob,
    alt.response.prob,
    null.reponse.prob,
    null.reponse.prob,
    null.reponse.prob,
    alt.response.prob,
    alt.response.prob,
    null.reponse.prob,
    null.reponse.prob,
    alt.response.prob,
    alt.response.prob,
    alt.response.prob
  ),
  nrow = 3,
  ncol = 4,
  byrow = T
)
model = "tlr" #logistic model
max.ar = 0.75  #limit the allocation ratio for the control group (1-max.ar < r_control < max.ar)
#------------Select the data generation randomisation methods-------
rand.type = "Urn" # Urn design
max.deviation = 3 # The recommended value for the tuning parameter in the Urn design

# Require multiple cores for parallel running
cl = 2

# Set the model we want to use and the time trend effect for each model used.
# Here the main model will be used twice for two different strength of time trend c(0,0,0,0) and c(1,1,1,1) to investigate how time trend affect the evaluation metrics in BAR setting.
# Then the main + stage_continuous model which is the treatment effect + stage effect model will be applied for strength equal c(1,1,1,1) to investigate how the main + stage effect model improve the evaluation metrics.
reg.inf = c("main + stage_continuous")
trend.effect = matrix(c(0.1, 0.1, 0.1, 0.1),
                      ncol = 4,
                      nrow = 1,
                      byrow = T)

#
cutoffearly = matrix(rep(0.994, dim(alt.scenario)[1]), ncol = 1)

K = dim(alt.scenario)[2]
print(
  paste0(
    "Start trial simulation. This is a ",
    K,
    "-arm trial simulation. There are one null scenario and ",
    K - 1 ,
    " alternative scenarios. There are ",
    K ,
    " rounds."
  )
)
cutoffindex = 1

## ----eval=FALSE---------------------------------------------------------------
#  
#  result = {
#  
#  }
#  OPCalt = {
#  
#  }
#  for (i in 1:dim(alt.scenario)[1]) {
#    trendindex = 1
#    for (j in 1:length(reg.inf)){
#      restlr = Trial.simulation(
#        ntrials = ntrials,
#        # Number of trial replicates
#        trial.fun = simulatetrial,
#        # Call the main function
#        input.info = list(
#          response.probs = alt.scenario[cutoffindex, ],
#          #The scenario vector in this round
#          ns = ns,
#          # Sequence of total number of accrued patients at each interim analysis
#          max.ar =  max.ar,
#          #limit the allocation ratio for the control group (1-max.ar < r_control < max.ar)
#          rand.type = rand.type,
#          # Which randomisation methods in data generation.
#          max.deviation = max.deviation,
#          # The recommended value for the tuning parameter in the Urn design
#          model.inf = list(
#            model = model,
#            #Use which model?
#            ibb.inf = list(
#              #independent beta-binomial model which can be used only for no time trend simulation
#              pi.star = 0.5,
#              # beta prior mean
#              pess = 2,
#              # beta prior effective sample size
#              betabinomialmodel = ibetabinomial.post # beta-binomial model for posterior estimation
#            ),
#            tlr.inf = list(
#              beta0_prior_mu = 0,
#              # Stan logistic model t prior location
#              beta1_prior_mu = 0,
#              # Stan logistic model t prior location
#              beta0_prior_sigma = 2.5,
#              # Stan logistic model t prior sigma
#              beta1_prior_sigma = 2.5,
#              # Stan logistic model t prior sigma
#              beta0_df = 7,
#              # Stan logistic model t prior degree of freedom
#              beta1_df = 7,
#              # Stan logistic model t prior degree of freedom
#              reg.inf =  reg.inf[trendindex],
#              # The model we want to use
#              variable.inf = "Fixeffect" # Use fix effect logistic model
#            )
#          ),
#          Stopbound.inf = Stopboundinf(
#            Stop.type = "Early-Pocock",
#            # Use Pocock like early stopping boundary
#            Boundary.type = "Symmetric",
#            # Use Symmetric boundary where cutoff value for efficacy boundary and futility boundary sum up to 1
#            cutoff = c(cutoffearly[cutoffindex, 1], 1 - cutoffearly[cutoffindex, 1]) # The cutoff value for stopping boundary
#          ),
#          Random.inf = list(
#            Fixratio = FALSE,
#            # Do not use fix ratio allocation
#            Fixratiocontrol = NA,
#            # Do not use fix ratio allocation
#            BARmethod = "Thall",
#            # Use Thall's Bayesian adaptive randomisation approach
#            Thall.tuning.inf = list(tuningparameter = "Fixed",  fixvalue = 1) # Specified the tunning parameter value for fixed tuning parameter
#          ),
#          trend.inf = list(
#            trend.type = "step",
#            # Linear time trend pattern
#            trend.effect = trend.effect[trendindex, ],
#            # Stength of time trend effect
#            trend_add_or_multip = "mult" # Multiplicative time trend effect on response probability
#          )
#        ),
#        cl = 2 # 2 cores required
#      )
#    trendindex = trendindex + 1
#    # The result list can be used for plotting and the OPC table is the summary evaluaton metrics for each scenario
#    result = c(result, restlr$result)
#    OPC_alt = rbind(OPC_alt, restlr$OPC)
#    }
#    cutoffindex = cutoffindex + 1
#  }

## -----------------------------------------------------------------------------
print("Finished alternative scenario study")
save_data = FALSE
if (isTRUE(save_data)) {
  save(result, file = restlr$Nameofsaveddata$nameData)
  save(OPC_alt, file = restlr$Nameofsaveddata$nameTable)
}


## -----------------------------------------------------------------------------
# Characteristic table
print(OPC_alt)

