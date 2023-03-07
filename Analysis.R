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





trial.sim = function(p.t, p.c, n, S, simulationSettingName){
  # p.t is proportion of successes we see on treatment
  # p.c is proprotion fo successes we see on control.
  # n is the total number of people that we need to enroll across 
  # both arms of the trial
  # S is the number of monte carlo replicates.
  # simulationSetting name is a string that allows us to save the result of the simulation
  results = matrix(NA, nrow= S, ncol = 8)
  colnames(results) = c("replicate", "p.11", "p.12", "T.1",
                        "p.21", "p.22", "T.2", "decision")
  for(r in 1:S){
    # interim analysis
    #get number of success for the controls
    p.11 = sum(rbinom(n/4, 1, prob = p.c))
    #get number of success for the treatment
    p.12 = sum(rbinom(n/4, 1, prob = p.t))
    T.1 = ifelse(prop.test(x = c(p.11, p.12), n = c(n/4, n/4))$p.value<0.005,
                 1, 0)
    
    
    
  }
}

