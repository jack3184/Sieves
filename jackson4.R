#Very simple
# 12 Feb, with correct(ing) (y,w)
n <- 10

# EAch indv has one kink in their budget constraint (including linear as a special case)
j <- round(rep(2, n) )
y <- matrix(NA, nrow=n, ncol=j[1])
w <- matrix(NA, nrow=n, ncol=j[1])
y[,1] <- pmax(0, rnorm(n, 4, .5))
w[,1] <- pmax(10, rnorm(n, 40, 20))
# kink is when labor tax kicks in
kink <- 20
tax.r <- .3
l <- pmin(1, kink/w[,1])
w[(l<1),2] <- w[(l<1),1]*(1-tax.r)
w[(l==1),2] <- 1
y[,2] <- l*(w[,1]-w[,2]) + y[,1]

cbind(y,w,l)


# V is the Cobb-douglas parameter (1-\alpha), drawn here from uniform
v.mean <- 0
v.sigm <- 1
v <- round(runif(n, 0, 1),1)

# # as mentioned above, pi is cobb-douglas. Allows pi.bar (mean over V), pi.inv and mu to be easily computed
pi <- function(y,w,v){ (v*(w+y)-y)/w }

pi.bar.fn <- function(y,w){  (1/2*(w+y)-y)/w   }
pi.inv.fn <- function(y,w,l) { (w*l+y)/(w+y)}
mu.fn <- function(y,w,l){ -1/2*(w*l+y)^2/(w*(w+y))}

# pi.bar <- pi.bar.fn(y[,2],w[,2])
# mu <- mu.fn(y,w,rep(l))

h.bar.fn <- function(y,w,l){ pi.bar.fn(y[,2],w[,2]) + mu.fn(y,w,rep(l))[,1] - mu.fn(y,w,rep(l))[,2]}
h.bar <- h.bar.fn(y,w,l)
plot(sort(h.bar))

# # here are some plots for interest
plot(l,h.bar)
plot(y[,1],h.bar)
plot(y[,2],h.bar)
plot(w[,1],h.bar)
plot(w[,2],h.bar)
# 
# 
# # Estimation
# h data is generated now
h <- pmax(.5, rnorm(n, .58, .1))
# # plot(sort(h-h.bar))
# # plot(sort(h))
# 
# # basis has just increasing powers, up to k
# k <- 2
# P <- matrix(NA, nrow=n, ncol= k*2)
# for (i in 1:k){
# P[,i] <- (y*w)[,2]^(i)
# P[,k+i] <- (l*(y[,1]*w[,1]-y[,2]*w[,2]))^(i)
# }
# 
# beta = ginv(t(P) %*% P) %*% t(P) %*% h
# 
# plot(sort(P %*% beta- h.bar))
# plot(sort(P %*% beta- h))
# 
# lines(sort((P %*% beta- h.bar)[h.bar>.51]))
# lines(sort((P %*% beta- h)[h.bar>.51]))
# 
# lines(sort((P %*% beta- h.bar)[h>.51]),col="blue")
# lines(sort((P %*% beta- h)[h>.51]), col="blue")
# # some observations: 
# # the estimation does better! if the basis is limited (see what happens when k = 8), why is this?
# # having a look at Table I and section 5.3, it sems like 'additional terms' are added v arbitrarily. Must read more on this
# # How to pick power functions? SHould we try non-uniformly increasing functions (like y^3*w^2) ?
# 
# 
