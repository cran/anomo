## ----fig.width = 7, fig.height = 3.5------------------------------------------
 library(anomo)
myci <- mcci(d1 = .1, se1 = .1, d2 = .3, se2 = .1)
# Note. Effect difference (the black square representing d2 - d1), 90% MCCI 
# (the thick horizontal line) for the test of equivalence, and 95% MCCI 
# (the thin horizontal line) for the test of moderation 
# (or difference in effects).


## -----------------------------------------------------------------------------
# Adjust the plot
myci <- mcci(d1 = .1, se1 = .1, d2 = .4, se2 = .1,
             eq.bd = c(-0.2, 0.2), xlim = c(-.2, .7))

## ----fig.width = 7, fig.height = 3.5------------------------------------------
MyCI.Mediation <- mcci(d1 = c(.60, .40), 
             se1 = c(.019, .025),
             d2 = c(.60, .80),
             se2 = c(.016, .023))

## ----conventional.power.analysis----------------------------------------------
 # 1. Conventional Power Analyses from Difference Perspectives
 # Calculate the required sample size to achieve certain level of power
 mysample <- power.eq.2group(d = .1, eq.dis = 0.1,  p =.5,
                             r12 = .5, q = 1, power = .8)
 mysample$out
 # Calculate power provided by a sample size allocation
 mypower <- power.eq.2group(d = 0.1, eq.dis = 0.1, n = 1238, p =.5,
                            r12 = .5, q = 1)
 mypower$out
 # Calculate minimum detectable distance a given sample size allocation can achieve
 myeq.dis <- power.eq.2group(d = .1, n = 1238, p =.5,
                            r12 = .5, q = 1, power = .8)
 myeq.dis$out

## ----power.analysis.with.costs------------------------------------------------
 # 2. Power Analyses Using Optimal Sample Allocation
 # Optimal sample allocation identification
 od <- od.eq.2group(r12 = 0.5, c1 = 1, c1t = 10)
 # Required budget and sample size at the optimal allocation
 budget <- power.eq.2group(expr = od, d = 0.1, eq.dis = 0.1,
                            q = 1, power = .8)  
 # Required budget and sample size by an balanced design with p = .50
 budget.balanced <- power.eq.2group(expr = od, d = 0.1, eq.dis = 0.1,
                                    q = 1, power = .8,
                                    constraint = list(p = .50))
 # 27% more budget required from the balanced design with p = 0.50.
 (budget.balanced$out$m-budget$out$m)/budget$out$m *100

## ----power.curve--------------------------------------------------------------
pwr <- NULL
p.range <- seq(0.01, 0.99, 0.01)
for(p in p.range){
  pwr <- c(pwr, power.eq.2group(expr = od, constraint = list(p = p),
                                m = budget$out$m, d = 0.1, eq.dis = 0.1,
                                q = 1, verbose = FALSE)$out$power)
}
plot(p.range*100,
     pwr*100,
     type = "l", lty = 1,
     xlim = c(0, 100), ylim = c(0, 100),
     xlab = "Proportion of Units in Treated (%)", ylab = "Power (%)",
     main = "", col = "black")
abline(v=od$out$p*100, lty = 2, col = "black")

