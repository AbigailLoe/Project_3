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
# implies we need at least 6147 montecarlo runs for power, and 1825 for alpha

# remains to be done: solve for the true margin of error for alpha!!!



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
    
    response.1 = matrix(c(p.11, n/2-p.11, p.12, n/2-p.12 ), nrow = 2)
    
    test.1 = fisher.test(response.1, alternative = "two.sided")
    Z.1 = test.1$estimate
    pval1= test.1$p.value
    T.1 = ifelse(pval1< boundary1, 1, 0 )
    # problem with the T.1 being NA right now.
    # using the O'Brien-Fleming boundary
    # 1 for rejecting the null, 0 for failing to reject the null.
    if(T.1 == 0){
      # if we fail to reject the null during the interim analysis:
      p.21 = p.11 + sum(rbinom(n/2, 1, prob = p.c))
      p.22 = p.12 + sum(rbinom(n/2, 1, prob = p.t))
      
      overall.n = 2*n
      response.2 = matrix(c(p.21, n-p.21, p.22, n-p.22 ), nrow = 2)
      test.2 = fisher.test(response.2, alternative = "two.sided")
      
      Z.2 = test.2$estimate
      pval2= test.2$p.value
      T.2 = ifelse(pval2< boundary2, 1, 0 )
      
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
set.seed(16); test1 = ( trial.sim(p.t = 0.15, p.c = 0.0, n = 80, S = ceiling(S.power), 
                    boundary1 = 0.01, boundary2 = 0.04 ))
colMeans(test1[, c("reject", "overall.n")])

set.seed(16); null1 = ( trial.sim(p.t = 0.0, p.c = 0.0, n = 80, S = ceiling(S.power), 
                                  boundary1 = 0.01, boundary2 = 0.04 ))
colMeans(null1[, c("reject", "overall.n")])


######## Setting 2 #########
#p1 = .05, p2 = .20, sig.level = 0.05, power = 0.80
power.prop.test(n= NULL, p1 = 0.05, p2 = 0.2,
                sig.level = 0.05, power= 0.80,
                alternative = "two.sided")
# initial n is 76 in each group, 152 total.
# this 38 in each group at the interim, 76 total

set.seed(16); test2 = ( trial.sim(p.t = 0.2, p.c = 0.05, n = 90, S = ceiling(S.power), 
                    boundary1 = 0.01, boundary2 = 0.04 ))
colMeans(test2[, c("reject", "overall.n")])
#140

set.seed(16); null2 = ( trial.sim(p.t = 0.05, p.c = 0.05, n = 76, S = ceiling(S.alpha), 
                    boundary1 = 0.01, boundary2 = 0.04))
colMeans(null2[, c("reject", "overall.n")])



######## Setting 3 #########

#p1 = .1, p2 = .25, sig.level = 0.05, power = 0.80
power.prop.test(n= NULL, p1 = 0.1, p2 = 0.25,
                sig.level = 0.05, power= 0.80,
                alternative = "two.sided")
# initial n is 100 in each group, 200 total.
# this 50 in each group at the interim, 100 total

set.seed(16); test3 = ( trial.sim(p.t = 0.25, p.c = 0.1, n = 116, S = ceiling(S.power), 
                    boundary1 = 0.01, boundary2 = 0.04 ))
colMeans(test3[, c("reject", "overall.n")])

set.seed(16); null3 = ( trial.sim(p.t = 0.1, p.c = 0.1, n = 116, S = ceiling(S.alpha), 
                                  boundary1 = 0.01, boundary2 = 0.04 ))
colMeans(null3[, c("reject", "overall.n")])



######## Setting 4 #########
#p1 = .15, p2 = .3, sig.level = 0.05, power = 0.80
power.prop.test(n= NULL, p1 = 0.15, p2 = 0.3,
                sig.level = 0.05, power= 0.80,
                alternative = "two.sided")
# initial n is 120 in each group, 240 total.
# this 60 in each group at the interim, 120 total
set.seed(16); test4 = ( trial.sim(p.t = 0.3, p.c = 0.15, n = 138, S = ceiling(S.power), 
                    boundary1 = 0.01, boundary2 = 0.04 ))
colMeans(test4[, c("reject", "overall.n")])


set.seed(16); null4 = ( trial.sim(p.t = 0.3, p.c = 0.3, n = 138, S = ceiling(S.alpha), 
                                  boundary1 = 0.01, boundary2 = 0.04 ))
colMeans(null4[, c("reject", "overall.n")])





######### In conclusion ########
# want to creat a graph that has the various powers and alpha levels
# also varying n 

pcDat = seq(from = 0.0, to = 0.15, by = 0.03)

nChoices = seq(from = 122, to = 140, by = 2); nChoices


finalDataFrame = as.data.frame(matrix(ncol = 6, nrow = length(pcDat)*length(nChoices)))
colnames(finalDataFrame) = c("propNull","expectedSampleAlt","expectedSampleNull", 
                             "typeIError", "power", "trialSize")

iter = 1
set.seed(16)
for (j in 1:length(nChoices)){
  print(j)
  nSize = nChoices[j]
  for(i in 1:length(pcDat)){
    phat = pcDat[i]
    
    alt = trial.sim(p.t = phat+0.15, p.c = phat, n = nSize, S = ceiling(S.power), 
                    boundary1 = 0.01, boundary2 = 0.04 )
    finalDataFrame[iter, c("power", "expectedSampleAlt", "trialSize")]=  
      c(colMeans(alt[, c("reject", "overall.n")]), 2*nSize)
    
    pointNull = trial.sim(p.t = phat, p.c = phat, n = nSize, S = ceiling(S.alpha), 
                          boundary1 = 0.01, boundary2 = 0.04)
    finalDataFrame[iter, c("propNull", "typeIError", "expectedSampleNull")] = 
      c(phat, colMeans(pointNull[,c("reject", "overall.n")]))
    iter = iter+1
  }
}


# write.csv(finalDataFrame, "simulationResultsFinal.csv")
# 
# 
library(ggplot2)
colors <- c("Type I" = "blue", "Type II" = "orange")

ggplot(data = finalDataFrame, aes(x = propNull))+
  geom_ribbon(aes(ymin = 0.2-m,  ymax = 0.2+m), fill = "lightgrey")+
  geom_abline(intercept = 0.2, slope = 0, color = "orange", linetype = 2) +
  geom_ribbon(aes(ymin = 0.05-m,  ymax = 0.05+m), fill = "lightgrey")+
  geom_line(aes(y = typeIError, color = "Type I"))+
  geom_line(aes(y = 1-power, color = "Type II"))+
  geom_abline(intercept = 0.05, slope = 0, color = "blue", linetype = 3) +
  labs(x = "Response P(MRB) under the Null",
       y = "Error Probability",
       color = "Error Types",
       title = "Type I and Type 2 Errors",
       subtitle = "(Stratified by Total Sample Size)")+
  scale_color_manual(values = colors)+
  facet_wrap(~trialSize)

# ggsave("Errors.png")


colors2 = c("Alternative" = "blue", "Null" = "orange")
ggplot(data = finalDataFrame, aes(x = propNull))+
  geom_line(aes(y = expectedSampleAlt, color = "Alternative"))+
  geom_line(aes(y = expectedSampleNull, color = "Null"))+
  labs(x = "Response P(MRB) under the Null",
       y = "Average Number of People Enrolled",
       color = "Settings",
       title = "Expected Sample Sizes",
       subtitle = "(Stratified by Total Possible Sample Size)")+
  scale_color_manual(values = colors2)+
  facet_wrap(~trialSize)

# ggsave("Expected.sample.sizes.png")
