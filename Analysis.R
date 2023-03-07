##### Game Plan ####
# Decide deltas/hypothesis situations to test
# This will allow us to decide correct n
# determine ideal number of simulations
# run simulations
# win.


#### Null Settings ####
######## Setting 1 #########
#p1 = 0, p2 = .15, sig.level = 0.05, power = 0.80
power.prop.test(n= NULL, p1 = 0, p2 = 0.15,
                sig.level = 0.05, power= 0.80,
                alternative = "two.sided")
# initial n is 48 in each group, 96 total.
# this 24 in each group at the interim, 48 total

# ---- Determining S ------
m = 0.02
S.power = (1.96/m)^2 *.8*(1-.8); S.power
S.alpha = (1.96/m)^2 *.95*(1-.95); S.alpha
# implies we need at least 1537 montecarlo runs.





######## Setting 2 #########


######## Setting 3 #########


######## Setting 4 #########


#### Alternative Settings ####