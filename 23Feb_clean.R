
##### Generate data #####
n <- 1000
### Budget sets: EAch indv has one kink in their budget constraint (including linear as a special case)
j <- round(rep(2, n) )
y <- matrix(NA, nrow=n, ncol=j[1])
w <- matrix(NA, nrow=n, ncol=j[1])
y[,1] <- pmax(0, rnorm(n, 10, .5))
# y[,1] <- rep(0,n)
w[,1] <- pmax(1, rnorm(n, 15, 10))
# kink is when labor tax kicks in
kink <- 10
l <- pmin(1, (kink)/w[,1])
tax.r <- .3
w[(l<1),2] <- w[(l<1),1]*(1-tax.r)
w[(l==1),2] <- 1
y[,2] <- 1
y[(l<1),2] <-  l[(l<1)]*(w[(l<1),1]-w[(l<1),2]) + y[(l<1),1]

head(cbind(y,w,l))


### V is noise on measurement of y
v.mean <- 0
v.sigm <- 1
v <- runif(n)

### pi is from u=(c^(1-a)+(1-h)^(1-a))/(1-a) (?CES)
a <- 4/3
pi <- function(y,w,v){ pmax(0, ((w)^a-(y*v))/(w^a+w)) }
# mean(pi(y,w,v))
# sum(pi(y,w,v))
# plot(sort(pi(y[,1],w[,1],rep(v))))

pi.bar.fn <- function(y,w){ 1/(w^a+w)*((w^a-y/2)^(w^a/y >= 1))*
    ((w^(a*2)/(y*2))^(w^a/y < 1))}
# lines(sort(pi.bar.fn(y[,1],w[,1])))

pi.inv.fn <- function(y,w,l) { (w^2-l*(w^2+w)/y)*(l>0)}
### notice inverse not well defined on some part of domain, mu.fn accounts for this

mu.fn <- function(y,w,l){ 
  0 * (0 >= (w^a-l*(w^a+w))/y ) +
  (w^(2*a)-2*l*w^a*(w^a+w) +l^2*(w^a+w)^2)/(2*y*(w^a+w))  * (1 >= (w^a-l*(w^a+w))/y &(w^a-l*(w^a+w))/y >= 0) +
  ((w^(a)-1/2*y)/(w^a+w)-l) * ((w^a-l*(w^a+w))/y >= 1) 
}
# plot(((w^a-l*(w^a+w))/y),mu.fn(y,w,l))
# head(mu.fn(y,w,l))
# head((w^a-l*(w^a+w))/y)

### h.bar.fn would have to be re-written for >2 segments
h.bar.fn <- function(y,w,l){ pi.bar.fn(y[,2],w[,2]) + mu.fn(y,w,rep(l))[,1] - mu.fn(y,w,rep(l))[,2]}
h.bar <- h.bar.fn(y,w,l)
plot(sort(h.bar))

### here are some plots for interest
# plot(l,h.bar)
# plot(w[,1],pi(y[,1],w[,1],rep(v)))
# plot(w[,1],h.bar)
# plot(w[,2],h.bar)
# plot(y[,1],h.bar)
# plot(y[,2],h.bar)

####### Estimation #####
### h data is generated now
h <- pi(y[,1],w[,1],v)
# plot(sort(h))

##### Generate Basis ####
## I need the .base function for the policy simulation
P.fn.base <- function(y,w,l){
P <- rep(1,n)
### adds columns in order of degree (try with something other than 10)
for (i in 1:20){
P.temp <- cbind(poly(cbind(y[,1],w[,1]), degree=i, raw=T)[1:n,], -poly(cbind(y[,2],w[,2],l), degree=i, raw=T)[1:n,]+
  poly(cbind(y[,1],w[,1],l), degree=i, raw=T)[1:n,])
head(P.temp)
P <- cbind(P,P.temp[,!(colnames(P.temp) %in% colnames(P))])
}
P
}
# head(P)
# print(ncol(P))

### now construct lin ind P matrix (li is index lin ind columns)
P.fn <- function(y,w,l){
P <- P.fn.base(y,w,l)
li <- 1
P.li <- matrix(P[,1])

for (i in 2:ncol(P)){
  # head(P.li)
  P.temp <- cbind(P.li[,which(li<i)], P[,i], P.li[,which(li > i)])
  colnames(P.temp)[max(max(which(li<i)+1),1)] <- colnames(P)[i]
  # head(P.temp)
  ## near-multicollinearity is a problem, the first 'if' ensure matrix invertibility; the second 'if' prevents the values from being too big for the computer
  ## http://stackoverflow.com/questions/38553336/get-wrong-inverse-matrix-from-both-solve-and-ginv-functions
  if ( ncol(P.temp) == qr(P.temp)$rank ) {
  # if ( ncol(t(P.temp) %*% P.temp) == qr(t(P.temp) %*%P.temp)$rank ) {
    P.li <- P.temp
    li <- sort(c(li,i))
  }
}
P.li
}
P <- P.fn(y,w,l)
ncol(P)


##### Do (sieve) regression ####
library(MASS)
beta <- function(K) {ginv(t(P[,1:K]) %*% P[,1:K]) %*% t(P[,1:K]) %*% h}
# beta <- function(K) {matrix(lm(h ~ P[,1:K]-1)$coefficients,ncol=1)}
hhat <- function(K) {P[,1:K] %*% beta(K)}

##### CV (vectorization is quicker than loops, change to `leave one out') ####
SSE <- function(K){
  temp = matrix(NA, nrow=n, ncol=1)
  hhat_l = hhat(K)
  P_l   = P[,1:K]
  temp <- (h-hhat_l)^2 / (1-diag(P_l %*% ginv(t(P_l) %*% P_l) %*% t(P_l)))^2
  sum(temp)
}

CV <- function(K){1 - SSE(K) / sum((h-h.bar)^2)}

Kmax = ncol(P)
CV.mat <- matrix(sapply(1:Kmax,CV), ncol=1)

CV.opt = max(CV.mat)
K.opt = which.max(CV.mat)

K.opt=22

pdf(file = 'plot1.pdf')
plot(ecdf(h),col="GREEN")
lines(ecdf(P[,1:K.opt]%*% beta(K.opt)), col="BLUE")
lines(ecdf(h.bar),col="RED")
dev.off()

#### policy simulation (\hat{M})
### .a is for `after', .b is for before
### seems pretty sensitive to changing parameters: it makes sense, since taking ^p of data, small initial changes mean big changes in basis
P.b <- P
y.a <- y
w.a <- w
l.a <- l

## min.wage
w.a[w.a<=(mean(w.a)*.4)] <- mean(w.a)*.4
l.a <- pmin(1, (kink)/w.a[,1])

## increase welfare
# y.a[,1] <- 15

## change tax
# tax.r <- .2
# w.a[(l.a<1),2] <- w.a[(l.a<1),1]*(1-tax.r)
# w.a[(l.a==1),2] <- 1
# y.a[,2] <- l.a*(w.a[,1]-w.a[,2]) + y.a[,1]

K.pol <- K.opt+5

P.a <- P.fn.base(y.a,w.a,l.a)
P.a <- cbind( rep(1,n), P.a[,colnames(P.a) %in% colnames(P)])

P.b.mean <- apply(P.b, 2, mean)
P.a.mean <- apply(P.a, 2, mean)
## functions
{
B <- function(K){ matrix(P.b.mean[1:K],nrow=K) }
D <- function(K) {matrix(P.a.mean[1:K],nrow=K)-B(K)}
theta.fn <- function(K){ (t(D(K)) %*% beta(K)) / (t(B(K)) %*% beta(K))}
C <- function(K) {(D(K) - rep(theta.fn(K))* B(K)) / rep(t(B(K)) %*% beta(K))}
Q.hat.fn <- function(K){ t(P.b[,1:K])%*% P.b[,1:K]/n}
Sigma.hat.fn <- function(K){x <- array(NA, dim=c(K,K,n))
  for(i in 1:n){t(t(P[i,1:K])) %*% t(P[i,1:K]) -> x[,,i]}
  x <- x * rep( (h-P.b[,1:K] %*% beta(K))^2, each=K*K)/rep(n)
  apply(x, c(1,2), sum)
}
SE.m.fn <- function(K){
  sqrt( t(C(K)) %*% ginv(Q.hat.fn(K)) %*% Sigma.hat.fn(K) %*% ginv(Q.hat.fn(K)) %*% C(K) /rep(n))
}
}

M.hat <- sapply(3:K.pol, theta.fn)
M.SE <-sapply(3:K.pol, SE.m.fn)
plot(M.hat)
plot(M.SE)

Additional.term <- colnames(P)[3:K.pol]

write.table(cbind(name,round(cbind(M.hat,M.SE),4)), 'table1.txt')

### (22 FEb)
# fix basis generation
# do M: done 
# standard errors?
# produce output for .tex
# 
