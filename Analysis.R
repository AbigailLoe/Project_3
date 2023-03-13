##### Game Plan ####
# Decide deltas/hypothesis situations to test
# This will allow us to decide correct n
# determine ideal number of simulations
# run simulations
# win.

# ---- Determining S ------
m = 0.01
S.power = (1.96/m)^2 *.8*(1-.8); S.power
S.alpha = (1.96/m)^2 *.95*(1-.95); S.alpha
# implies we need at least 1537 montecarlo runs.



trial.sim = function(p.t, p.c, n, S, simulationSettingName= NA, 
                     boundary1 = 0.005, boundary2 = 0.048){
  # p.t is proportion of successes we see on treatment
  # p.c is proportion of successes we see on control.
  # n is the number in each arm (assuming equal arms)
  # S is the number of monte carlo replicates.
  # simulationSetting name is a string that allows us to save the result of the simulation
  results = matrix(NA, nrow= S, ncol = 11)
  colnames(results) = c("replicate", "p.11", "p.12", "Z.1","T.1",
                        "p.21", "p.22","Z.2", "T.2", "reject", "overall.n")
  for(r in 1:S){
    # interim analysis
    #get number of success for the controls
    p.11 = sum(rbinom(n/2, 1, prob = p.c))
    #get number of success for the treatment
    p.12 = sum(rbinom(n/2, 1, prob = p.t))
    neum = asin(sqrt(p.11/(n/2)))-asin(sqrt(p.12/(n/2)))
    denom = sqrt(1/(4* n/2)+ 1/(4*n/2))
   
     Z.1 = abs(neum/denom)
    #need to perform a two-sided test.
    pval1 = 2*(1-pnorm(Z.1))
    
    T.1 = ifelse(pval1< boundary1, 1, 0 )
    # problem with the T.1 being NA right now.
    # using the O'Brien-Fleming boundary
    # 1 for rejecting the null, 0 for failing to reject the null.
    if(T.1 == 0){
      # if we fail to reject the null during the interim analysis:
      p.21 = p.11 + sum(rbinom(n/2, 1, prob = p.c))
      p.22 = p.12 + sum(rbinom(n/2, 1, prob = p.t))
      
      overall.n = 2*n
      
      neum = asin(sqrt(p.21/n))-asin(sqrt(p.22/n))
      denom =sqrt(1/(4* n)+ 1/(4*n))
      
      Z.2 = abs(neum/denom)
      #need to perform a two-sided test.
      pval2 = 2*(1-pnorm(Z.2))
      
      T.1 = ifelse(pval2< boundary2, 1, 0 )
      
      # T.2 = ifelse(prop.test(x = c(p.21, p.22), n = c(n/2, n/2))$p.value<boundary2,
      #              1, 0)
    }else{
      #rejected the null, found a difference in the interim analysis:
      p.21 = 0
      p.22 = 0
      overall.n = n
      Z.2 = 0
      T.2 = 0
    }
    
    results[r, ] = c(r, p.11, p.12, Z.1, T.1, p.21, p.22, Z.2, T.2, T.1+T.2, overall.n)  
  }
  if(!is.na(simulationSettingName)){
    write.csv(paste0(simulationSettingName, ".csv"))
  }
  return(results)
}




######## Setting 1 #########
#p1 = 0, p2 = .15, sig.level = 0.05, power = 0.80
power.prop.test(n= NULL, p1 = 0, p2 = 0.15,
                sig.level = 0.05, power= 0.80,
                alternative = "two.sided")
# initial n is 48 in each group, 96 total.
# this 24 in each group at the interim, 48 total
set.seed(16)
test1 = ( trial.sim(p.t = 0.15, p.c = 0.00, n = 30, S = round(S.power), 
                  boundary1 = 0.005, boundary2 = 0.04))
colMeans(test1[, c("reject", "overall.n")])

null1 = trial.sim(p.t = 0, p.c = 0.00, n = 1000, S = round(S.power), 
                  boundary1 = 0.005, boundary2 = 0.04)
colMeans(null1[, c("reject", "overall.n")])
# type 1 error is miniscule here! this feels off.

######## Setting 2 #########
#p1 = .05, p2 = .20, sig.level = 0.05, power = 0.80
power.prop.test(n= NULL, p1 = 0.05, p2 = 0.2,
                sig.level = 0.05, power= 0.80,
                alternative = "two.sided")
# initial n is 76 in each group, 152 total.
# this 38 in each group at the interim, 76 total

set.seed(16)

test2 = ( trial.sim(p.t = 0.2, p.c = 0.05, n = 70, S = round(S.power), 
                    boundary1 = 0.005, boundary2 = 0.048 ))
colMeans(test2[, c("reject", "overall.n")])
#`168 or 172 are sufficient

set.seed(16); null2 = ( trial.sim(p.t = 0.05, p.c = 0.05, n = 70, S = round(S.power), 
                    boundary1 = 0.005, boundary2 = 0.048))
colMeans(null2[, c("reject", "overall.n")])
# really large type I error here...

set.seed(16); null2 = ( trial.sim(p.t = 0.05, p.c = 0.05, n = 700, S = round(S.power), 
                                  boundary1 = 0.005, boundary2 = 0.048))
colMeans(null2[, c("reject", "overall.n")])
# better if we make the sample size much much larger.

`
######## Setting 3 #########

#p1 = .1, p2 = .25, sig.level = 0.05, power = 0.80
power.prop.test(n= NULL, p1 = 0.1, p2 = 0.25,
                sig.level = 0.05, power= 0.80,
                alternative = "two.sided")
# initial n is 100 in each group, 200 total.
# this 50 in each group at the interim, 100 total

set.seed(16); test3 = ( trial.sim(p.t = 0.25, p.c = 0.1, n = 96, S = round(S.power), 
                    boundary1 = 0.005, boundary2 = 0.048 ))
colMeans(test3[, c("reject", "overall.n")])

set.seed(16); null3 = ( trial.sim(p.t = 0.1, p.c = 0.1, n = 100, S = round(S.power), 
                                  boundary1 = 0.005, boundary2 = 0.04 ))
colMeans(null3[, c("reject", "overall.n")])

set.seed(16); null3 = ( trial.sim(p.t = 0.1, p.c = 0.1, n = 110, S = round(S.power), 
                                  boundary1 = 0.005, boundary2 = 0.04 ))
colMeans(null3[, c("reject", "overall.n")])

######## Setting 4 #########
#p1 = .15, p2 = .3, sig.level = 0.05, power = 0.80
power.prop.test(n= NULL, p1 = 0.15, p2 = 0.3,
                sig.level = 0.05, power= 0.80,
                alternative = "two.sided")
# initial n is 120 in each group, 240 total.
# this 60 in each group at the interim, 120 total
set.seed(16); test4 = ( trial.sim(p.t = 0.3, p.c = 0.15, n = 120, S = round(S.power), 
                    boundary1 = 0.005, boundary2 = 0.045 ))
colMeans(test4[, c("reject", "overall.n")])


set.seed(16); null4 = ( trial.sim(p.t = 0.3, p.c = 0.3, n = 120, S = round(S.power), 
                                  boundary1 = 0.005, boundary2 = 0.045 ))
colMeans(null4[, c("reject", "overall.n")])

