# dataGenCD
# 
# This script generates budget set and hours data, assuming at most one kink in the budget set  
# and Cobb-Douglas preferences. The kink is due to an income tax effective above a certain 
# threshold income (no tax below the kink).
#
# Authors: Jackson Bunting, Attila Gyetvai
# jackson.bunting@duke.edu, attila.gyetvai@duke.edu
#
# Feb 19, 2017


# Initializing variables --------------------------------------------

n <- 100    # Number of obs

# w: Wage rate; slope
w     <- matrix(NA, nrow=n, ncol=2)
w.m   <- 40    # Mean of w
w.v   <- 20    # Var of w
w.min <- 10    # Minimum wage rate
w[,1] <- pmax(w.min, rnorm(n, w.m, w.v))

# y: Nonlabor income; intercept
y     <- matrix(NA, nrow=n, ncol=2)
y.m   <- 4     # Mean of y
y.v   <- .5    # Var of y
y.min <- 0     # Minimum nonlabor income
y[,1] <- pmax(y.min, rnorm(n, y.m, y.v))

# Kink
b.k <- 20                 # Income level at kink
t.k <- .3                 # Tax rate above kink
b.f <- y[,1] + w[,1]*1    # Budget if use entire time endowment for work
i.k <- b.f>=b.k           # Indicator of kink

# Hours at the kink
l      <- rep(1, n)
l[i.k] <- (b.k-y[i.k,1])/w[i.k,1]

# Wage rate above the kink
w[,2]     <- w[,1]
w[i.k,2]  <- (1-t.k) * w[i.k,1]

# Hypothetical nonlabor income for segment above the kink
y[,2]    <- y[,1]
y[i.k,2] <- b.k - w[i.k,2]*l[i.k]


# Choosing labor supply ---------------------------------------------

# v: Cobb-Douglas parameter, drawn from uniform
v.m <- 0
v.v <- 1
v   <- runif(n, v.m, v.v)

# Functions from Blomquist and Newey (2002)
pi.fn     <- function(y,w,v){ 1 - v - v*y/w }         # Labor supply
pi.bar.fn <- function(y,w)  { 1/2 - y/(2*w) }         # Integrating v out
pi.inv.fn <- function(y,w,l){ w*(1-l)/(w+y) }         # Inverse of pi wrt v
mu.fn     <- function(y,w,l){ w*(1-l)^2/(2*(w+y)) }

# Actual hours
h.1 <- pmin(pmax(0, pi.fn(y[,1], w[,1], v)), 1)
h.2 <- pmin(pmax(0, pi.fn(y[,2], w[,2], v)), 1)
c.1 <- y[,1] + w[,1]*h.1
c.2 <- y[,2] + w[,2]*h.2
U.1 <- c.1^(1-v)*(1-h.1)^v
U.2 <- c.2^(1-v)*(1-h.2)^v
i.U <- U.2>=U.1

h       <- matrix(NA, nrow=n, ncol=1)
h[i.U]  <- h.2[i.U]
h[!i.U] <- h.1[!i.U]

# Conditional expectation of hours
pi.bar <- pi.bar.fn(y[,2], w[,2])
mu.1   <- mu.fn(y[,1], w[,1], l)
mu.2   <- mu.fn(y[,2], w[,2], l)
h.bar  <- pi.bar + mu.1 - mu.2


