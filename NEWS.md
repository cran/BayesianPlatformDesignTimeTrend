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

# May 1, 2023, BayesianPlatformDesignTimeTrend version 1.0.1
## Major changes
- I fixed a bug where the fix ratio allocation method code is not suitable for the multi-arm trial

- Old code: rpk[-1] = rep((armleft - 1) / (armleft - 1 + Fixratiocontrol), armleft - 1)
- New code: rpk[-1] = rep(1 / (armleft - 1 + Fixratiocontrol), armleft - 1)
