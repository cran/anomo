## ----fig.width = 7, fig.height = 3.5------------------------------------------
 library(anomo)
myci <- mcci(d = c(.1, .15), se = c(.01, .01))
# Note. Effect difference (the black square representing d1 - d2), 90% MCCI 
# (the thick horizontal line) for the test of equivalence, and 95% MCCI 
# (the thin horizontal line) for the test of moderation 
# (or difference in effects).


## -----------------------------------------------------------------------------
# Adjust the plot
myci <- mcci(d = c(.1, .15), se = c(.01, .01),
             eq.bd = c(-0.2, 0.2), xlim = c(-.2, .7))

## ----fig.width = 7, fig.height = 3.5------------------------------------------
MyCI.Mediation <- mcci(d = c(.60, .40, .60, .80), 
             se = c(.019, .025, .016, .023), mediation = TRUE)
#Note. The order of d is a1, b1, a2, and b2 (e.g., treatment-mediator
#   and mediator-outcome path in group/study 1 and 2, respectively). 
#   se is in the same order for the standard errors.

## ----conventional.power.analysis----------------------------------------------
 # 1. Conventional Power Analyses from Difference Perspectives
 # Calculate the required sample size to achieve certain level of power
 mysample <- power.1.eq(d = .1, eq.dis = 0.1,  p =.5,
                             r12 = .5, q = 1, power = .8)
 mysample$out
 # Calculate power provided by a sample size allocation
 mypower <- power.1.eq(d = 0.1, eq.dis = 0.1, n = 1238, p =.5,
                            r12 = .5, q = 1)
 mypower$out
 # Calculate minimum detectable distance a given sample size allocation can achieve
 myeq.dis <- power.1.eq(d = .1, n = 1238, p =.5,
                            r12 = .5, q = 1, power = .8)
 myeq.dis$out

## ----power.analysis.with.costs------------------------------------------------
 # 2. Power Analyses Using Optimal Sample Allocation
 # Optimal sample allocation identification
 od <- od.1.eq(r12 = 0.5, c1 = 1, c1t = 10)
 # Required budget and sample size at the optimal allocation
 budget <- power.1.eq(expr = od, d = 0.1, eq.dis = 0.1, 
                         power = .8)  
 # Required budget and sample size by an balanced design with p = .50
 budget.balanced <- power.1.eq(expr = od, d = 0.1, eq.dis = 0.1,
                                    power = .8,
                                    constraint = list(p = .50))
 # 27% more budget required from the balanced design with p = 0.50.
 (budget.balanced$out$m-budget$out$m)/budget$out$m *100

## ----power.curve--------------------------------------------------------------
plot.power.eq(expr = od, d = 0.1, eq.dis = 0.1) 

