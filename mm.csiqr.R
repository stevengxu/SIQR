library(Rlab)
set.seed(1)

n <- 500
beta <- c(1, 1, 1)/sqrt(3)
X <- matrix(stats::runif(length(beta)*n,-1,1), ncol=length(beta))
beta <- sign(beta[1])*beta/sqrt(sum(beta^2));
si <- X%*%beta
y <- sin(pi*si)+rnorm(n,0,0.5)
XX <- rowMeans(X)*0.5
expp=exp(XX)/(1+exp(XX))
t <- rbern(n,expp)

## y
y[t==0] <- 0

# PS
pp <- as.vector(glm(t~X, family = "binomial")$fitted.values)

mm.csiqr <- function(X, y, t, pp, tau = 0.50, isControl = FALSE, max.iter = 500, tol = 1e-5){
  
  p <- NCOL(X)
  n <- length(y)
  eps <- tol/n
  if (isControl) {
    t  <- 1-t
    pp <- 1-pp
  }
  yo  <- y[t==1]   # observed response
  no  <- length(yo)
  Xo  <- X[t==1,]  # observed covariates
  ppo <- pp[t==1]  # observed ps
  
  # initialize
  fit0 <- rq(y~X,weights=t,tau=tau)
  beta <- coef(fit0)[-1]
  beta <- sign(beta[1])*beta/sqrt(sum(beta^2))
  
  alfa <- matrix(0,nrow=n,ncol=2)
  alfa[,1] <- fit0$fitted.values # g
  alfa[,2] <- 1 # g.prime
  
  si  <- X%*%beta # single index
  hm  <- KernSmooth::dpill(si[t==1], yo);
  hp  <- h <- hm*(tau*(1-tau)/(stats::dnorm(stats::qnorm(tau)))^2)^.2
  
  z <- mapply(function(j,i){
    si[j]-si[i] # z_ij = si_j - si_i
  },seq_along(si))
  z <- z[,t==1]
  wp <- w <- dnorm(z/h)
  wp <- wp/rowSums(wp)
  
  for (iter in 1:max.iter) {
    plot(si[t==1],yo)
    lines(si[order(si)], alfa[order(si),1], col=2)
    AA <- bb <- 0
    alfa.new <- alfa
    for (i in 1:n) {
      ## updating a_i and b_i
      v <- cbind(1, z[i,])
      s <- 1/(eps + as.vector(abs(yo-v%*%alfa[i,])))
      #A <- t(v)%*%diag(s*w[i,]/ppo)%*%v
      A <- crossprod(v*s*w[i,]/ppo,v)
      r <- s*yo+2*tau-1
      b <- crossprod(v, r*w[i,]/ppo)
      #for (ii in 1:n) {
      #  A <- A + t(v)%*%diag(s*wp[ii,]*w[i,]*(1-t[ii]/pp[ii]))%*%v
      #  b <- b + crossprod(v, (s*yo+2*tau-1)*wp[ii,]*w[i,]*(1-t[ii]/pp[ii]))
      #}
      for (j in 1:no) {
        ww <- sum((1-t/pp)*wp[,j])*w[i,j]
        A  <- A + ww*s[j]*tcrossprod(v[j,],v[j,])
        b  <- b + ww*r[j]*v[j,]
      }
      alfa.new[i,] <- solve(A,b)
      ## update residuals
      s  <- 1/(eps + as.vector(abs(yo-v%*%alfa[i,])))
      u  <- alfa[i,2]*sweep(Xo,2,X[i,],"-")
      AA <- AA + t(u)%*%diag(s*w[i,]/ppo)%*%u
      AA <- AA + crossprod(u*s*w[i,]/ppo,u)
      r  <- s*(yo-alfa[i,1])+2*tau-1
      bb <- bb + crossprod(u, r*w[i,]/ppo)
      for (j in 1:no) {
        ww <- sum((1-t/pp)*wp[,j])*w[i,j]
        AA <- AA + ww*s[j]*tcrossprod(u[j,],u[j,])
        bb <- bb + ww*r[j]*u[j,]
      }
    }
    if (sum((alfa.new[,1]-alfa[,1])^2)/sum(alfa[,1]^2) < tol) break
    alfa <- alfa.new
    ## update beta
    beta.new <- solve(AA,bb)
    beta.new <- sign(beta.new[1])*beta.new/sqrt(sum(beta.new^2))
    beta <- beta.new
    # update weights
    si <- X%*%beta
    #if (iter%%10==0) {
    #  hm  <- KernSmooth::dpill(si[t==1], yo);
    #  hp  <- h <- hm*(tau*(1-tau)/(stats::dnorm(stats::qnorm(tau)))^2)^.2
    #}
    z <- mapply(function(j,i){
      si[j]-si[i] # z_ij = si_j - si_i
    },seq_along(si))
    z <- z[,t==1]
    wp <- w <- dnorm(z/h)
    wp <- wp/rowSums(wp)
  }
  return(list(alfa=alfa, beta=beta, iter=iter))
}


fit <- mm.csiqr(X,y,t,pp,tau=0.75)
lines(si[order(si)],sin(pi*si[order(si)])+qnorm(0.75,0,0.5))






qt.ipw <- mm.csiqr(X,y,t,pp,0.5,"IPW",FALSE)
qt.eep <- mm.csiqr(X,y,t,pp,0.5,"EEP",FALSE)
qt.aipw <- mm.csiqr(X,y,t,pp,0.5,"AIPW",FALSE)

qc.ipw <- mm.csiqr(X,y,t,pp,0.5,"IPW",TRUE)
qc.eep <- mm.csiqr(X,y,t,pp,0.5,"EEP",TRUE)
qc.aipw <- mm.csiqr(X,y,t,pp,0.5,"AIPW",TRUE)

qt.qr <- mm.siqr(X[t==1,],y[t==1],max.iter=100)
qc.qr <- mm.siqr(X[t==0,],y[t==0],max.iter=100)


qt.ans <- rowMeans(X[,1:4]) + abs(X[,1])*qnorm(0.5)
qc.ans <- sin(2*pi*X[,1])+qnorm(0.5,0,0.3)


par(mfrow=c(1,4))
si <- X%*%qt.ipw$beta
plot(si,y1)
ipw.fit <- qt.ipw$alfa[,1]
lines(si[order(si)],ipw.fit[order(si)],col=2,lwd=2)
lines(si[order(si)],qt.ans[order(si)],col=3,lwd=2)

si <- X%*%qt.eep$beta
plot(si,y1)
eep.fit <- qt.eep$alfa[,1]
lines(si[order(si)],eep.fit[order(si)],col=2,lwd=2)
lines(si[order(si)],qt.ans[order(si)],col=3,lwd=2)

si <- X%*%qt.aipw$beta
plot(si,y1)
aipw.fit <- qt.aipw$alfa[,1]
lines(si[order(si)],aipw.fit[order(si)],col=2,lwd=2)
lines(si[order(si)],qt.ans[order(si)],col=3,lwd=2)

si <- X%*%qt.qr$beta
plot(si,y1)
qr.fit <- qt.qr$alfa[,1]
lines(si[order(si)],qr.fit[order(si)],col=2,lwd=2)
lines(si[order(si)],qt.ans[order(si)],col=3,lwd=2)



a <- rnorm(3)
b <- rnorm(3)
z <- cbind(1,rnorm(3))

A <- 0
for (i in 1:3) {
  A <- A + a[i]*b[i]*z[i,]%*%t(z[i,])
}

B <- t(z)%*%diag(a*b)%*%z
