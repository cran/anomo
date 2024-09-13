## ----fig.width = 7, fig.height = 3.5------------------------------------------
 library(anomo)
myci <- mcci(d1 = .1, se1 = .1, d2 = .3, se2 = .1)
# Note. Effect difference (the black square representing d2 - d1), 90% MCCI (the thick horizontal line) for the test of equivalence, and 95% MCCI (the thin horizontal line) for the test of moderation (or difference in effects).


## ----setup--------------------------------------------------------------------
# Adjust the plot
myci <- mcci(d1 = .1, se1 = .1, d2 = .4, se2 = .1,
             bound.eq = c(-0.1, 0.2), xlim = c(-.2, .7))

