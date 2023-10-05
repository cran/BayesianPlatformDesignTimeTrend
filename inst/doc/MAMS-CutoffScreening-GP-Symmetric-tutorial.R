## ---- include = FALSE---------------------------------------------------------
 knitr::opts_chunk$set(
   collapse = TRUE,
   comment = "#>"
 )

## ----setup--------------------------------------------------------------------
library(BayesianPlatformDesignTimeTrend)

## ----eval=FALSE---------------------------------------------------------------
#  
#  ntrials = 1000 # Number of trial replicates
#  ns = seq(120, 600, 60) # Sequence of total number of accrued patients at each interim analysis
#  null.reponse.prob = 0.15
#  alt.response.prob = 0.35
#  
#  # We investigate the type I error rate for different time trend strength
#  null.scenario = matrix(
#    c(
#      null.reponse.prob,
#      null.reponse.prob,
#      null.reponse.prob,
#      null.reponse.prob
#    ),
#    nrow = 1,
#    ncol = 4,
#    byrow = T
#  )
#  # alt.scenario = matrix(c(null.reponse.prob,null.reponse.prob,null.reponse.prob,null.reponse.prob,
#  #                     null.reponse.prob,alt.response.prob,null.reponse.prob,null.reponse.prob,
#  #                     null.reponse.prob,alt.response.prob,alt.response.prob,null.reponse.prob,
#  #                     null.reponse.prob,alt.response.prob,alt.response.prob,alt.response.prob), nrow=3, ncol = 4,byrow=T)
#  model = "tlr" #logistic model
#  max.ar = 0.85  #limit the allocation ratio for the control group (1-max.ar < r_control < max.ar)
#  #------------Select the data generation randomisation methods-------
#  rand.type = "Urn" # Urn design
#  max.deviation = 3 # The recommended value for the tuning parameter in the Urn design
#  
#  # Require multiple cores for parallel running
#  cl = 2
#  
#  # Set the model we want to use and the time trend effect for each model used.
#  # Here the main model will be used twice for two different strength of time trend c(0,0,0,0) and c(1,1,1,1) to investigate how time trend affect the evaluation metrics in BAR setting.
#  # Then the main + stage_continuous model which is the treatment effect + stage effect model will be applied for strength equal c(1,1,1,1) to investigate how the main + stage effect model improve the evaluation metrics.
#  reg.inf = "main"
#  trend.effect = c(0,0,0,0)
#  
#  result = {
#  
#  }
#  OPC = {
#  
#  }
#  K = dim(null.scenario)[2]
#  cutoffindex = 1
#  trendindex = 1
#  
#  cutoff.information=demo_Cutoffscreening.GP (
#    ntrials = ntrials,
#    # Number of trial replicates
#    trial.fun = simulatetrial,
#    # Call the main function
#    grid.inf = list(
#      start.length = 10,
#      grid.min = NULL,
#      grid.max = NULL,
#      confidence.level = 0.95,
#      grid.length = 5000,
#      change.scale = FALSE,
#      noise = T,
#      errorrate = 0.1,
#      simulationerror = 0.01,
#      iter.max = 15,
#      plotornot = FALSE),
#    # Set up the cutoff grid for screening. The start grid has three elements. The extended grid has fifteen cutoff value under investigation
#    input.info = list(
#      response.probs = null.scenario[1,],
#      #The scenario vector in this round
#      ns = ns,
#      # Sequence of total number of accrued patients at each interim analysis
#      max.ar =  max.ar,
#      #limit the allocation ratio for the control group (1-max.ar < r_control < max.ar)
#      rand.type = rand.type,
#      # Which randomisation methods in data generation.
#      max.deviation = max.deviation,
#      # The recommended value for the tuning parameter in the Urn design
#      model.inf = list(
#        model = model,
#        #Use which model?
#        ibb.inf = list(
#          #independent beta-binomial model which can be used only for no time trend simulation
#          pi.star = 0.5,
#          # beta prior mean
#          pess = 2,
#          # beta prior effective sample size
#          betabinomialmodel = ibetabinomial.post # beta-binomial model for posterior estimation
#        ),
#        tlr.inf = list(
#          beta0_prior_mu = 0,
#          # Stan logistic model t prior location
#          beta1_prior_mu = 0,
#          # Stan logistic model t prior location
#          beta0_prior_sigma = 2.5,
#          # Stan logistic model t prior sigma
#          beta1_prior_sigma = 2.5,
#          # Stan logistic model t prior sigma
#          beta0_df = 7,
#          # Stan logistic model t prior degree of freedom
#          beta1_df = 7,
#          # Stan logistic model t prior degree of freedom
#          reg.inf =  reg.inf,
#          # The model we want to use
#          variable.inf = "Fixeffect" # Use fix effect logistic model
#        )
#      ),
#      Stop.type = "Early-OBF",
#      # Use Pocock like early stopping boundary
#      Boundary.type = "Symmetric",
#      # Use Symmetric boundary where cutoff value for efficacy boundary and futility boundary sum up to 1
#      Random.inf = list(
#        Fixratio = FALSE,
#        # Do not use fix ratio allocation
#        Fixratiocontrol = NA,
#        # Do not use fix ratio allocation
#        BARmethod = "Thall",
#        # Use Thall's Bayesian adaptive randomisation approach
#        Thall.tuning.inf = list(tuningparameter = "Unfixed",  fixvalue = 1) # Specified the tunning parameter value for fixed tuning parameter
#      ),
#      trend.inf = list(
#        trend.type = "linear",
#        # Linear time trend pattern
#        trend.effect = trend.effect,
#        # Stength of time trend effect
#        trend_add_or_multip = "mult" # Multiplicative time trend effect on response probability
#      )
#    ),
#    cl = 2
#  )
#  

## -----------------------------------------------------------------------------
library(ggplot2)
# Details of grid
optimdata=optimdata_sym
# Recommend cutoff at each screening round
nextcutoff = optimdata$next.cutoff
prediction = optimdata$prediction
cutoff=optimdata$cutoff
tpIE=optimdata$tpIE
cutoff=cutoff[1:sum(!is.na(tpIE))]
tpIE=tpIE[1:sum(!is.na(tpIE))]
GP.res = optimdata
prediction = data.frame(yhat = GP.res$prediction$yhat.t1E,
                        sd = matrix(GP.res$prediction$sd.t1E,ncol=1),
                        qup = GP.res$prediction$qup.t1E,
                        qdown = GP.res$prediction$qdown.t1E,
                        xgrid = GP.res$prediction$xgrid)
GPplot=ggplot(data = prediction) +
        geom_ribbon(aes(x = xgrid, ymin = qdown, ymax = qup),col="#f8766d", alpha = 0.5,linetype = 2) +
        geom_line(aes(xgrid, yhat),col = "#f8766d") +
        geom_point(aes(cutoff[1:sum(!is.na(tpIE))], tpIE[1:sum(!is.na(tpIE))]),
                   data = data.frame(tpIE=tpIE,cutoff=cutoff),col = "#00bfc4") +
        geom_point(aes(nextcutoff, 0.1),
                   data = data.frame(tpIE=tpIE,cutoff=cutoff),col = "#f8766d") +
        geom_hline(yintercept = 0.1,linetype = 2) +
        geom_text(aes(x=1,y=0.15,label=paste0("FWER target is 0.1")),hjust=0,vjust=1)+
        geom_vline(xintercept = nextcutoff, linetype = 2) +
        geom_text(aes(x=6,y=0.8,label=paste0("Next cutoff value is ",round(nextcutoff,3))))+
        theme_minimal()+ylab("FWER")+xlab("Cutoff value of the OBF boundary (c*)")+
  geom_point(aes(nextcutoff, 0.1),
                   data = data.frame(tpIE=tpIE,cutoff=cutoff),col = "#f8766d") +
  theme(plot.background = element_rect(fill = "#e6dfba"))
print(GPplot)

