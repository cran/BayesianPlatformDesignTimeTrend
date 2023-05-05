# May 2, 2023, BayesianPlatformDesignTimeTrend version 1.1.0
This release implements additional method of cutoff screening using Bayesian optimization.
## New feature
- Add a new Demo function called `demo_Cutoffscreening_GP`. The function indicates how to use Gaussian process-based Bayesian optimisation algorithm to calibrate the critical value for stopping boundary under the null scenario to control the Type I error rate or FWER before trial simulation. 

- Add a function `GP_optim` which returns a cutoff value for the next evaluation during the cutoff screening process. The function is used in demo `demo_Cutoffscreening_GP`.

## Bugs
- I found that the output data matrix can not be generated  when there are four arm, the 'variable.inf'="Mixeffect" and 'reg.inf' is "main"/"main + stage_continuous"/"main * stage_discrete".
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

