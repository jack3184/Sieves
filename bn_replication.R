# bn_replication
# 
# This script replicates the study by Blomquist and Newey (2002 ECMA). We simulate budget set
# and hours data assuming one kink in the budget set. The kink is due to an income tax effective 
# above a certain threshold income (no tax below the kink). First we assume Cobb-Douglas preferences,
# second we assume isoelastic utility. Then we estimate the labor supply function via sieve LS. We 
# then estimate the supply response to a counterfactual tax rate change policy.
#
# Authors: Jackson Bunting, Attila Gyetvai
# jackson.bunting@duke.edu, attila.gyetvai@duke.edu
#
# Feb 28, 2017


#### Initialization ####
rm(list = ls())

#### Simulating budget set data ####

n <- 1000    # Number of obs

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


#### Labor supply -- Cobb-Douglas preferences ####

# v: Cobb-Douglas parameter, drawn from uniform
v.m <- 0
v.v <- 1
v   <- runif(n, v.m, v.v)

# Functions from Blomquist and Newey (2002)
pi.fn     <- function(y,w,v){ 1 - v - v*y/w }                  # Labor supply
pi.bar.fn <- function(y,w)  { 1/2 - y/(2*w) - w/(2*(w+y)) }    # Integrating v out
pi.inv.fn <- function(y,w,l){ w*(1-l)/(w+y) }                  # Inverse of pi wrt v
mu.fn     <- function(y,w,l){ w/(w+y)*(1-l-(1-l)^2/2 - 1/2) }

# Actual hours
h.1  <- pmin(pmax(0, pi.fn(y[,1], w[,1], v)), 1)
h.2  <- pmin(pmax(0, pi.fn(y[,2], w[,2], v)), 1)
# Bunching at the kink
i.l1      <- h.1 >= l
i.l2      <- h.2 <= l
h.1[i.l1] <- l[i.l1]
h.2[i.l2] <- l[i.l2]

c.1  <- y[,1] + w[,1]*h.1
c.2  <- y[,2] + w[,2]*h.2
U.1  <- c.1^(1-v)*(1-h.1)^v
U.2  <- c.2^(1-v)*(1-h.2)^v
i.U  <- U.2>=U.1

h       <- matrix(NA, nrow=n, ncol=1)
h[i.U]  <- h.2[i.U]
h[!i.U] <- h.1[!i.U]
eps     <- rnorm(n, 0, 0.5)
h       <- pmax(pmin(h+eps, 1), 0)

# Conditional expectation of hours
pi.bar <- pi.bar.fn(y[,2], w[,2])
mu.1   <- mu.fn(y[,1], w[,1], l)
mu.2   <- mu.fn(y[,2], w[,2], l)
h.bar  <- pi.bar + mu.1 - mu.2


##### Estimation routine #####

# Generating basis
d.max <- 20    # Maximum degree of polynomial considered

P.fn.base <- function(y, w, l){
  P <- rep(1, n)
  # Adds columns in order of degree
  for (i in 1:d.max){
    P.temp <- cbind(poly(cbind(y[,1], w[,1]), degree=i, raw=T)[1:n,], 
              - poly(cbind(y[,2], w[,2], l), degree=i, raw=T)[1:n,]
              + poly(cbind(y[,1], w[,1], l), degree=i, raw=T)[1:n,])
    head(P.temp)
    P <- cbind(P, P.temp[,!(colnames(P.temp) %in% colnames(P))])
  }
  P
}

# Constructing lin ind P matrix (li is index lin ind columns)
P.fn <- function(y, w, l){
  P    <- P.fn.base(y, w, l)
  li   <- 1
  P.li <- matrix(P[,1])
  
  for (i in 2:ncol(P)){
    P.temp <- cbind(P.li[,which(li<i)], P[,i], P.li[,which(li>i)])
    colnames(P.temp)[max(max(which(li<i) + 1),1)] <- colnames(P)[i]
    # Near-multicollinearity is a problem. The first 'if' ensures matrix invertibility; 
    # the second 'if' prevents the values from being too big for the computer. More here:
    # http://stackoverflow.com/questions/38553336/get-wrong-inverse-matrix-from-both-solve-and-ginv-functions
    if (ncol(P.temp) == qr(P.temp)$rank ) {
    # if ( ncol(t(P.temp) %*% P.temp) == qr(t(P.temp) %*%P.temp)$rank ) {
      P.li <- P.temp
      li <- sort(c(li,i))
    }
  }
  P.li
}
P <- P.fn(y,w,l)
ncol(P)

# Sieve regression
library(MASS)
beta <- function(K) {ginv(t(P[,1:K]) %*% P[,1:K]) %*% t(P[,1:K]) %*% h}
# beta <- function(K) {matrix(lm(h ~ P[,1:K]-1)$coefficients,ncol=1)}    # Alternative OLS coefs
hhat <- function(K) {P[,1:K] %*% beta(K)}

# Leave-one-out CV
betaLOO <- function(K, i) {ginv(t(P[-i,1:K]) %*% P[-i,1:K]) %*% t(P[-i,1:K]) %*% h[-i]}
ghatLOO <- function(K, i) {P[i,1:K] %*% betaLOO(K, i)}

SSE <- function(K){
  temp <- matrix(NA, nrow=n, ncol=1)
  P_l  <- P[,1:K]
  ghat <- matrix(sapply(1:n, ghatLOO, K=K), ncol=1)
  temp <- (h-ghat)^2 / (1-diag(P_l %*% ginv(t(P_l) %*% P_l) %*% t(P_l)))^2
  sum(temp)
}

CV <- function(K){1 - SSE(K) / sum((h-mean(h))^2)}


#### Running estimation ####

SieveReg <- function (K.max, plotfilename) {
  # This function implements the sieve regression routine outlined above. K.max is the
  # max number of polynomials (not degree!) considered, plotfilename is the name of the
  # saved .pdf file of the CDF of results.
  t.start <- Sys.time ()    # Timer on
  CV.mat  <- matrix(sapply(1:K.max,CV), ncol=1)
  CV.opt  <- max(CV.mat)
  K.opt   <- which.max(CV.mat)
  t.end   <- Sys.time () - t.start    # Timer off
  
  pdf(file = paste(plotfilename, '.pdf', sep=''))
  plot(ecdf(h), col="GREEN")                            # CDF of h
  lines(ecdf(P[,1:K.opt] %*% beta(K.opt)), col="BLUE")  # CDF of h.hat
  lines(ecdf(h.bar), col="RED")                         # CDF of h.bar
  dev.off()
  
  return(list(CV.mat, CV.opt, K.opt, t.end))
}
results1 <- SieveReg(20, 'fig1')


#### Labor supply -- isoelastic utility ####
# v now represents measurement error

a  <- 4/3    # Utility fn param; alpha in the write-up is 1/a = 3/4

# Functions from Blomquist and Newey (2002)
pi        <- function(y,w,v) { pmax(0, ((w)^a - (y*v)) / (w^a+w)) }    # Labor supply
pi.bar.fn <- function(y,w)   { 1/(w^a+w)*((w^a-y/2)^(w^a/y >= 1))*
                              ((w^(a*2)/(y*2))^(w^a/y < 1)) }          # Integrating v out
#pi.inv.fn <- function(y,w,l) { (w^2 - l*(w^2+w)/y) * (l>0)}            # Inverse of pi wrt v
pi.inv.fn <- function(y,w,l) { ((w^a - l*(w^a+w))/y) * (l>0)}            # Inverse of pi wrt v
# Notice: inverse not well defined on some part of domain, mu.fn accounts for this
mu.fn     <- function(y,w,l) { 
  0 * (0 >= (w^a-l*(w^a+w))/y ) +
  (w^(2*a)-2*l*w^a*(w^a+w) + l^2*(w^a+w)^2)/(2*y*(w^a+w)) * (1 >= (w^a-l*(w^a+w))/y &(w^a-l*(w^a+w))/y >= 0) +
  ((w^(a)-1/2*y)/(w^a+w)-l) * ((w^a-l*(w^a+w))/y >= 1) 
}

# Simulated hours
h.bar.fn <- function(y,w,l){ pi.bar.fn(y[,2],w[,2]) + mu.fn(y,w,rep(l))[,1] - mu.fn(y,w,rep(l))[,2]}
h.bar    <- h.bar.fn(y,w,l)
h        <- pi(y[,1],w[,1],v)

# Basis
P <- P.fn(y,w,l)
ncol(P)

# Estimation
results2 <- SieveReg(20, 'fig2')

#### Policy simulation ####
# .a is for `after', .b is for `before'
P.b <- P
y.a <- y
w.a <- w
l.a <- l

# Minimum wage: new budget sets
w.a[w.a<=(mean(w.a)*.4)] <- mean(w.a)*.4
l.a                      <- pmin(1, (b.k)/w.a[,1])

# Generating new basis
K.opt <- as.integer(results2[3])
K.pol <- K.opt+5

P.a <- P.fn.base(y.a,w.a,l.a)
P.a <- cbind(rep(1,n), P.a[,colnames(P.a) %in% colnames(P)])

P.b.mean <- apply(P.b, 2, mean)
P.a.mean <- apply(P.a, 2, mean)

# These functions follow B&N section 3.4
B            <- function(K){ matrix(P.b.mean[1:K], nrow=K) }
D            <- function(K){ matrix(P.a.mean[1:K], nrow=K) - B(K)}
theta.fn     <- function(K){ (t(D(K)) %*% beta(K)) / (t(B(K)) %*% beta(K)) }
C            <- function(K){ (D(K) - rep(theta.fn(K))* B(K)) / rep(t(B(K)) %*% beta(K)) }
Q.hat.fn     <- function(K){ t(P.b[,1:K]) %*% P.b[,1:K]/n }
Sigma.hat.fn <- function(K){
  x <- array(NA, dim=c(K,K,n))
  for(i in 1:n){ t(t(P[i,1:K])) %*% t(P[i,1:K]) -> x[,,i] }
  x <- x * rep( (h-P.b[,1:K] %*% beta(K))^2, each=K*K)/rep(n)
  apply(x, c(1,2), sum)
}
SE.m.fn       <- function(K){
  sqrt( t(C(K)) %*% ginv(Q.hat.fn(K)) %*% Sigma.hat.fn(K) %*% ginv(Q.hat.fn(K)) %*% C(K) /rep(n))
}

M.hat <- sapply(3:K.pol, theta.fn)
M.SE  <-sapply(3:K.pol, SE.m.fn)
plot(M.hat)
plot(M.SE)

Additional.term <- colnames(P)[3:K.pol]
CV <- round(results2[[1]][3:K.pol],4)
K <- 3:K.pol

library(xtable)
table1 <- xtable(cbind(K,Additional.term,CV,round(cbind(M.hat,M.SE),4)),label=NULL)
print.xtable(table1, type="latex", file="table1.tex", include.rownames=FALSE )


#### END OF CODE ####


