# estim
# 
# TODO describe estimation routine
#
# Authors: Jackson Bunting, Attila Gyetvai
# jackson.bunting@duke.edu, attila.gyetvai@duke.edu
#
# Feb 19, 2017

P.sq <- function(K){ cbind(rep(1, n), poly(cbind(y[,1], w[,1]), degree=K, raw=T)[1:n,]) }
P.mu <- function(K){ -poly(cbind(l, y[,2], w[,2]), degree=K, raw=T)[1:n,]
  + poly(cbind(l, y[,1], w[,1]), degree=K, raw=T)[1:n,]}
P    <- function(K){ cbind(P.sq(K), P.mu(K)) }

library(MASS)
beta <- function(K) {ginv(t(P(K)) %*% P(K)) %*% t(P(K)) %*% h}
hhat <- function(K) {P(K) %*% beta(K)}

# Minimize CV
SSE <- function(K){
  temp = matrix(NA, nrow=n, ncol=1)
  hhat_l = hhat(K)
  P_l   = P(K)
  for (i in 1:n) {
    temp[i] = (h[i]-hhat_l[i])^2 / (1-t(P_l[i,]) %*% ginv(t(P_l) %*% P_l) %*% P_l[i,] )^2
  }
  sum(temp)
}

CV <- function(K){1 - SSE(K) / sum((h-h.bar)^2)}

Kmax = 4
CV.mat = matrix(NA, nrow=Kmax, ncol=1)
for (i in 1:Kmax) {
  CV.mat[i] = CV(i)
}

CV.opt = min(CV.mat)
K.opt = which.min(CV.mat)
