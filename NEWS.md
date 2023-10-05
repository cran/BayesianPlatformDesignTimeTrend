# Sept 30, 2023, BayesianPlatformDesignTimeTrend version 1.2.0
## Summary of changes
- This release corrects the boundary construction formula of OBF boundary; 
- This release adds internal functions of calculating conjunctive power, disconjunctive power and marginal power.
- This release extend the type of hypothesis testing which are Oneside test and the two side test now; 
- This release implements method of cutoff screening for asymmetric boundary using Bayesian optimization. The cutoff screening for asymmetric boundary is released in this version (including No early stopping boundary, Pocock boundary and OBF boundary). The plot of boundary shape refers to the plot on github: https://github.com/ZXW834/BayesianPlatformDesignTimeTrend.
- This release implements method of hyperparameter tuning for Trippa's adaptive randomisation approach. There are some flow chats in the tutorial and on the github.
- This release adds my supervisor David Woods as the contributor to the package.

## Bugs
- In OBF boundary construction, the boundary formula was $\theta_j=\phi \left(\sqrt{\frac{J}{j}c} \right)$ which is not the same as that in the reference which is $\theta_j=\phi \left(\sqrt{\frac{J}{j}}c \right)$. Although the original boundary formula works in trial simulation, the correction is still made.

## New feature
- In previous version, the power function was only for disconjunctive power and FWER. In this version, the conjunctive power `conjuncativepower_or_FWER` and marginal power (per-hypothesis type I error)  `perHtypeIerror_marginalpowerfunc` can also be calculated for the used of design evaluation and power optimisation in asymmetric boundary cutoff screening. The function to calculate disconjunctive power is renamed as `disconjunctivepowerfunc`.

- In previous version, the hypothesis testing during the interim analysis is $H_0$: $\pi_k = \pi_0$; $H_1$: $\pi_k \neq \pi_0$. Now we add one side test to only conclude for superiority: $H_0$: $\pi_k \leq \pi_0$; $H_1$: $\pi_k > \pi_0$.

- In previous version, the cutoff value can only be tuned for symmetric boundary. In this version, asymmetric boundary cutoff screening is added so that user can decide whether to use a more aggressive efficacy boundary with a more aggressive futility boundary (or a more conservative efficacy boundary with more conservative futility boundary). The asymmetric boundary screening process first select the type I error rate contour equal to the target, then pick a cutoff pair (eff, fut) that maximizes the power.

- In previous version, the Trippa's adaptive randomisation approach used fixed hyperparameters suggested in their paper. In this version, the hyperparameters a and b are tuned to optimise power when type I error (FWER) is controlled at target. This idea can be found on Trippa's paper (2012). We used a gaussian process model to smooth the contour plot. The pair (a, b) with highest predicted power will be exploit for hyperparameter screening until there are enough point overlapped on the contour plot. In addition, the hyperparameter point (a, b) can be pick in a crude way to save computation time via setting the maximum number of points to be looked at.

# Aug 22, 2023, BayesianPlatformDesignTimeTrend version 1.1.3
This release fix one command in fix effect model analysis of main effect that is confusing.

## Bugs
- I fixed the bug in function `resultrtostats`. The calculation of "stats4" was the sum of main effect and interaction effect when using the model with interaction term (main * stage_continuous / main * stage_discrete). The sum was used to calculate the posterior probability of superior or inferior. However, "stats4" only represents the mean values of main effect in the final output matrix. The wrong value was given to "stats4". This is the same for "stats5" which is the standard errors of the main effect in the model

# Jun 25, 2023, BayesianPlatformDesignTimeTrend version 1.1.2
This release fix two bugs.

## Bugs
- I fixed the bug in function `GP.optim` where the formula of information weighed randomisation is wrong.
- I fixed the bug in function `demo_Cutoffscreening` where the nextcutoff vector for sample may have only one element. This will lead to the error in function `sample` when you only want to sample one value greater than 1. The argument 'ntrials' in each example should be large (> 100) instead 2 to make the example more like an actual simulation example. I use ntrials = 2 in the example to speed up the check process.

# Jun 11, 2023, BayesianPlatformDesignTimeTrend version 1.1.1
This release fix one bug reported by the cran team.

## Bugs
- I fixed the bug in function `GP.optim` where the nextcutoff vector for sample may have only one element. This will lead to the error in function `sample` when you only want to sample one value greater than 1.
- The argument 'ntrials' in each example should be large (> 100) instead 2 to make the example more like an actual simulation example. I use ntrials = 2 in the example to speed up the check process.

# May 2, 2023, BayesianPlatformDesignTimeTrend version 1.1.0
This release implements additional method of cutoff screening using Bayesian optimization.

## New feature
- Add a new Demo function called `demo_Cutoffscreening_GP`. The function indicates how to use Gaussian process-based Bayesian optimisation algorithm to calibrate the critical value for stopping boundary under the null scenario to control the Type I error rate or FWER before trial simulation. 

- Add a function `GP_optim` which returns a cutoff value for the next evaluation during the cutoff screening process. The function is used in demo `demo_Cutoffscreening_GP`.

## Bugs
- I found that the output data matrix can not be generated  when there are four arm, the 'variable.inf' = "Mixeffect" and 'reg.inf' is "main"/"main + stage_continuous"/"main * stage_discrete".
It only work for 'variable.inf'="Mixeffect" and 'reg.inf = main + stage_discrete'. The reason why this happen is that the length of stats6 and stats7 variable are not generated correctly at stage 1 since the fixed effect model is used for 'variable.inf'="Mixeffect" at stage 1 to speed up the simulation.

## Debug
- I add one ifelse command to enforce the generation of stats6 and stats7 when 'variable.inf'="Mixeffect". It does not depend on the value of 'reg.inf' now.

# May 1, 2023, BayesianPlatformDesignTimeTrend version 1.0.1
## Major changes
- I fixed a bug where the fix ratio allocation method code is not suitable for the multi-arm trial

- Old code: rpk[-1] = rep((armleft - 1) / (armleft - 1 + Fixratiocontrol), armleft - 1)
- New code: rpk[-1] = rep(1 / (armleft - 1 + Fixratiocontrol), armleft - 1)

# Apr 25, 2023, BayesianPlatformDesignTimeTrend version 1.0.1
## Major changes
- I added a check command for max.ar which is the upper boundary for randomisation ratio for each arm. The command is added in the MainFunction.R. Details are shown below.
-   #-max.ar check
  if (1 - max.ar > 1/K){
    stop("Error: The lower allocation ratio should be at least 1/K. Please check the number of arm at the beginning and the max.ar")
  }
  
- Another change is the randomisation ratio adjustment in the Simulation_AdaptiveRandomisationmethodRatioCalc.R.

- I modified the randomisation ratio adjustment command for Thall's approach.

- ##----Original code------ (This code only protects the control arm's allocation ratio)
- rpk = matrix(rep(0, armleft), ncol = armleft)
- randomprob = matrix(rep(0, K), ncol = K)
- colnames(rpk) = c(1, treatmentindex + 1)
- colnames(randomprob) = seq(1, K)
- rpk[1] = alloc.prob.best[1]
- rpk[-1] = alloc.prob.best[-1][treatmentindex]
- rpk = rpk / sum(rpk)
- lower = ifelse(rpk[1] < (1 - max.ar), 1 - max.ar, rpk[1])
- rpk[1] = lower
- upper = ifelse(rpk[1] > max.ar, max.ar, rpk[1])
- rpk[1] = upper
- rpk[-1] = (1 - rpk[1]) * rpk[-1] / sum(rpk[-1])
- randomprob[as.numeric(colnames(rpk))] = randomprob[as.numeric(colnames(rpk))] + rpk

- ##----New code------ (This code protects all arms' allocation ratio. Also, this code work together with max.ar check code)
- rpk = matrix(rep(0,armleft),ncol = armleft)
- randomprob = matrix(rep(0,K),ncol = K)
- colnames(rpk) = c(1,treatmentindex+1)
- colnames(randomprob) = seq(1,K)
- rpk[1] = alloc.prob.best[1]
- rpk[-1] = alloc.prob.best[-1][treatmentindex]
- rpk = rpk/sum(rpk)
- lower = ifelse(rpk<(1-max.ar),1-max.ar,rpk)
- rpk = lower
- upper = ifelse(rpk>max.ar,max.ar,rpk)
- rpk = upper
- rpk[!(rpk==(1-max.ar))]=(1-sum(rpk[rpk==1-max.ar]))*(rpk[!(rpk==1-max.ar)]/sum(rpk[!(rpk==1-max.ar)]))
- randomprob[as.numeric(colnames(rpk))] = randomprob[as.numeric(colnames(rpk))]+rpk

