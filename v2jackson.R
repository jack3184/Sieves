#Very simple
# 14 Feb, with correct(ing) (y,w)
n <- 100

# EAch indv has one kink in their budget constraint (including linear as a special case)
j <- round(rep(2, n) )
y <- matrix(NA, nrow=n, ncol=j[1])
w <- matrix(NA, nrow=n, ncol=j[1])
y[,1] <- pmax(0, rnorm(n, 10, .5))
# y[,1] <- rep(0,n)
w[,1] <- pmax(1, rnorm(n, 15, 10))
# kink is when labor tax kicks in
kink <- 20
tax.r <- .3
l <- pmin(1, kink/w[,1])
w[(l<1),2] <- w[(l<1),1]*(1-tax.r)
w[(l==1),2] <- 1
y[,2] <- l*(w[,1]-w[,2]) + y[,1]

head(cbind(y,w,l))


# V is noise on measurement of y
v.mean <- 0
v.sigm <- 1
v <- (runif(n, 0, 1))

a <- 8/6

# pi is from u=(c^(1-a)+(1-h)^(1-a))/(1-a)
pi <- function(y,w,v){ pmax(0, ((w)^a-(y*v))/(w^a+w)) }
# pi <- function(y,w,v){ pmax((v*(w+y)-y)/w,0) }
mean(pi(y,w,v))
sum(pi(y,w,v))
plot(sort(pi(y[,1],w[,1],rep(v))))


pi.bar.fn <- function(y,w){ 1/(w^a+w)*((w^a-y/2)^(w^a/y >= 1))*
    ((w^(a*2)/(y*2))^(w^a/y < 1))}
lines(sort(pi.bar.fn(y[,1],w[,1])))

pi.inv.fn <- function(y,w,l) { 
  (w^2-l*(w^2+w)/y)*(l>0)}

mu.fn <- function(y,w,l){ 
  0*(0 >= (w^a-l*(w^a+w))/y ) +
    (w^(2*a)-2*l*w^a*(w^a+w) +l^2*(w^a+w)^2)/(2*y*(w^a+w))*(1 >= (w^a-l*(w^a+w))/y &(w^a-l*(w^a+w))/y >= 0) +
    ((w^(a)-1/2*y)/(w^a+w)-l)*((w^a-l*(w^a+w))/y >= 1) 
}
plot(((w^a-l*(w^a+w))/y),mu.fn(y,w,l))
head(mu.fn(y,w,l))
head((w^a-l*(w^a+w))/y)

# pi.bar <- pi.bar.fn(y[,2],w[,2])
# mu <- mu.fn(y,w,rep(l))

h.bar.fn <- function(y,w,l){ pi.bar.fn(y[,2],w[,2]) + mu.fn(y,w,rep(l))[,1] - mu.fn(y,w,rep(l))[,2]}
h.bar <- h.bar.fn(y,w,l)
plot(sort(h.bar))
# 
# # here are some plots for interest
# plot(l,h.bar)
# plot(w[,1],pi(y[,1],w[,1],rep(v)))
# plot(w[,1],h.bar)
# plot(w[,2],h.bar)
# plot(y[,1],h.bar)
# plot(y[,2],h.bar)
#
#
# # Estimation
# h data is generated now
h <- pi(y[,1],w[,1],v)
plot(sort(h))

#
# # basis has just increasing powers, up to k
K.sq <- 7
K.mu <- 3

P.sq <- cbind(rep(1,n),poly(cbind(y[,1],w[,1]), degree=K.sq, raw=T)[1:n,])
P.mu <- -poly(cbind(l,y[,2],w[,2]), degree=K.mu, raw=T)[1:n,]+
  poly(cbind(l,y[,1],w[,1]), degree=K.mu, raw=T)[1:n,]
P <- cbind(P.sq, P.mu)

# check for linear dependence: ginv makes this redundant, but interesting to see
rankifremoved <- sapply(1:ncol(P), function (x) qr(P[,-x])$rank)
print(which(rankifremoved != max(rankifremoved)))
head(P)

library(MASS)
beta = ginv(t(P) %*% P) %*% t(P) %*% h
# beta = solve(t(P) %*% P)

plot(sort(P %*% beta- h.bar))
lines(sort(P %*% beta- h))

# sum(P %*% beta - h)

plot(ecdf(h),col="GREEN")
lines(ecdf(P%*% beta), col="BLUE")
lines(ecdf(h.bar),col="RED")

## OLD:
# some observations:
# the estimation does better! if the basis is limited (see what happens when k = 8), why is this?
# having a look at Table I and section 5.3, it sems like 'additional terms' are added v arbitrarily. Must read more on this
# How to pick power functions? SHould we try non-uniformly increasing functions (like y^3*w^2) ?
# # Next steps:
# Why these basis functions? They seem to suck for CD
# How pick the basis? Do his CV?
# gaussian, explain why screw up with big k, run another function, do 2/3-variable taylor polynomial
