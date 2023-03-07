##### Game Plan ####
# Decide deltas/hypothesis situations to test
# This will allow us to decide correct n
# determine ideal number of simulations
# run simulations
# win.

# ---- Determining S ------
m = 0.02
S.power = (1.96/m)^2 *.8*(1-.8); S.power
S.alpha = (1.96/m)^2 *.95*(1-.95); S.alpha
# implies we need at least 1537 montecarlo runs.


#### Alternative Settings ####
######## Setting 1 #########
#p1 = 0, p2 = .15, sig.level = 0.05, power = 0.80
power.prop.test(n= NULL, p1 = 0, p2 = 0.15,
                sig.level = 0.05, power= 0.80,
                alternative = "two.sided")
# initial n is 48 in each group, 96 total.
# this 24 in each group at the interim, 48 total





######## Setting 2 #########
#p1 = .05, p2 = .20, sig.level = 0.05, power = 0.80
power.prop.test(n= NULL, p1 = 0.05, p2 = 0.2,
                sig.level = 0.05, power= 0.80,
                alternative = "two.sided")
# initial n is 76 in each group, 152 total.
# this 38 in each group at the interim, 76 total


######## Setting 3 #########

#p1 = .1, p2 = .25, sig.level = 0.05, power = 0.80
power.prop.test(n= NULL, p1 = 0.1, p2 = 0.25,
                sig.level = 0.05, power= 0.80,
                alternative = "two.sided")
# initial n is 100 in each group, 200 total.
# this 50 in each group at the interim, 100 total


######## Setting 4 #########
#p1 = .15, p2 = .3, sig.level = 0.05, power = 0.80
power.prop.test(n= NULL, p1 = 0.15, p2 = 0.3,
                sig.level = 0.05, power= 0.80,
                alternative = "two.sided")
# initial n is 120 in each group, 240 total.
# this 60 in each group at the interim, 120 total


#### Null Settings ####

######## Setting 1 #########



p.t = 0.2
p.c = 0.3
n = 200
S = 50

power.prop.test(n = 100, p1 = .3, p2 = .2, sig.level = 0.05, 
                alternative = "two.sided")

simulationSettingName = "test"

trial.sim = function(p.t, p.c, n, S, simulationSettingName, boundary1 = 0.005, boundary2 = 0.048){
  # p.t is proportion of successes we see on treatment
  # p.c is proprotion fo successes we see on control.
  # n is the total number of people that we need to enroll across 
  # both arms of the trial
  # S is the number of monte carlo replicates.
  # simulationSetting name is a string that allows us to save the result of the simulation
  results = matrix(NA, nrow= S, ncol = 9)
  colnames(results) = c("replicate", "p.11", "p.12", "T.1",
                        "p.21", "p.22", "T.2", "decision", "overall.n")
  for(r in 1:S){
    # interim analysis
    #get number of success for the controls
    p.11 = sum(rbinom(n/4, 1, prob = p.c))
    #get number of success for the treatment
    p.12 = sum(rbinom(n/4, 1, prob = p.t))
    T.1 = ifelse(prop.test(x = c(p.11, p.12), n = c(n/4, n/4))$p.value<boundary1,
                 1, 0)
    # using the O'Brien-Fleming boundary
    # 1 for rejecting the null, 0 for failing to reject the null.
    if(T.1 == 0){
      # if we fail to reject the null during the interim analysis:
      p.21 = p.11 + sum(rbinom(n/4, 1, prob = p.c))
      p.22 = p.12 + sum(rbinom(n/4, 1, prob = p.t))
      overall.n = n
      T.2 = ifelse(prop.test(x = c(p.21, p.22), n = c(n/2, n/2))$p.value<boundary2,
                   1, 0)
    }else{
      #rejected the null, found a difference in the interim analysis:
      p.21 = 0
      p.22 = 0
      overall.n = n/2
      T.2 = 0
    }

  results[r, ] = c(r, p.11, p.12, T.1, p.21, p.22, T.2, T.1+T.2, overall.n)  
  }
}

