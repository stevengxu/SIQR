#' Data generation function for simulation and demonstration
#' There are three settings.
#'
#' @param n sample size
#' @param true.theta true single-index coefficients,
#' default is c(1,1,1)/sqrt(3) for setting 1 and c(1,2)/sqrt(5) for other settings
#' @param sigma the standard deviation of error term
#' @param setting chose from three settings
#' @param ncopy generates multiple copies of data for Monte Carlo simulations
#'
#' @return X  predictors
#' @return Y  response variables
#' @return single.index.values single index term
#' @export
generate.data <- function(n,true.theta=NULL,sigma=0.1,setting="setting1",ncopy=1){
  
  if(setting == "setting1"){
    
    #parameter setting
    true.theta = if(is.null(true.theta)) c(1, 1, 1)/sqrt(3) else true.theta
    c1 = sqrt(3)/2-1.645/sqrt(12) #0.3912
    c2 = sqrt(3)/2+1.645/sqrt(12)#1.3409
    
    X = matrix(stats::runif(length(true.theta)*n), ncol=length(true.theta))
    true.theta = sign(true.theta[1])*true.theta/sqrt(sum(true.theta^2)) 
    U = X%*%true.theta
    si = sin( (U-c1)*pi/(c2 -c1) )
    y = si + stats::rnorm(length(si),0,sigma)
    if(ncopy>1){
      ylist <- lapply(vector(mode = "list", length = ncopy),function(x){si + stats::rnorm(n,0,sigma)})
    }
  }else if(setting == "setting2"){
    
  }else if(setting == "setting3"){
    true.theta = if(is.null(true.theta)) c(1, 2)/sqrt(5) else true.theta
    X = matrix(stats::rnorm(length(true.theta)*n), ncol=length(true.theta))
    U = X%*%true.theta
    si = 5*cos(U)+exp(-U^2)
    e = stats::rexp(n,rate=.5)
    y = si+e 
    if(ncopy>1){
      ylist <- lapply(vector(mode = "list", length = ncopy),function(x){si + stats::rexp(n,rate=.5)})
    }
  }
  
  if(ncopy>1){
    return(list("X" = X, "Y" = ylist,"single.index.values"=si))
  }else{
    return(list("X" = X, "Y" = y,"single.index.values"=si))
  }
}


#' A supporting function that return the local polynomial regression quantile.
#' This estimates the quantile and its derivative at the point x.0
#'
#' @param x covariate sequence 
#' @param y response sequence 
#' @param t treatment vector 
#' @param pp propensity score vector 
#' @param h bandwidth(scalar) 
#' @param tau - left-tail probability
#' @param x0 point at which the quantile is estimated
#' @param t0 treatment of x0
#' @param pp0 propsensity score of x0
#' @param isControl whether conditional quantile of control is estimated
#'
#' @return x0  a scalar
#' @return fv  quantile est  dv - quantile derivative est
clprq0 <- function (x, y, t, pp, h, tau = 0.5, x0)  #used in step 1 of the algorithm
{
  z  <- x - x0
  wx <- stats::dnorm(z/h)*t
  px <- t/pp
  wx <- px*wx
  r  <- quantreg::rq(y ~ z, weights = wx, tau = tau, ci = FALSE)
  fv <- r$coef[1]
  dv <- r$coef[2]
  return(list(x0 = x0, fv = fv, dv = dv))
}


#' Main estimation function of single index quantile regression model.
#' a two step method.
#'
#' @param y response vector 
#' @param X covariate matrix 
#' @param t treatment/exposure vector 
#' @param pp propensity score vector 
#' @param tau left-tail probability (quantile index), scalar
#' @param beta.initial starting value of beta, the single index coefficients
#' @param h user-defined bandwidth
#' @param maxiter max iteration number
#' @param tol toleration for convergence
#'
#' @return a siqr object, which includes:
#'         beta - the fitted single index coefficients with unit norm and first component being non negative
#'         flag.conv  - whether the iterations converge
#' @examples
#' #generate data
#' set.seed(2021)
#' data <- generate.data(50)
#' X <- data$X
#' y0<- data$Y
#'
#' #initials
#' beta0 <- NULL
#' #quantile
#' tau = 0.75
#' siqr.result <- siqr(y0,X,beta.initial = beta0, tau=tau)
#' summary(siqr.result)
#'
#' @export
csiqr <- function (y, X, t, pp, tau=0.5, isControl=FALSE, 
                   beta.initial=NULL, h=NULL, maxiter=100, tol=1e-8)
{
  n <- length(y)
  d <- ncol(X) 
  
  if (isControl) {
    t <- 1-t
    pp <- 1-pp
  }
  
  ii <- rep(1:n, each=n)
  jj <- rep(1:n, n)
  pwts <- t[ii]/pp[ii]
  diff <- X[ii,,drop=FALSE]-X[jj,,drop=FALSE]
  
  if (is.null(beta.initial)) {
    wts.initial <- ifelse(isControl, t, 1-t)
    beta.initial <- stats::coef(quantreg::rq(y~X, weights=t, tau=tau))[-1]
    beta.initial[1] <- abs(beta.initial[1])
  }
  flag.conv <- FALSE  #flag whether maximum iteration is achieved
  
  beta.new<-sign(beta.initial[1])*beta.initial/sqrt(sum(beta.initial^2)) 
  
  a <- rep(0,n)
  b <- rep(0,n)
  
  iter <- 0 
  beta.old <- 2*beta.new 
  
  
  while ((iter < maxiter) & (sum((beta.new-beta.old)^2)>tol))

  {
    #print(iter)
    #print(beta.new)
    
    beta.old <- beta.new 
    iter <- iter+1 
    ####################################
    #  step 1: compute a.j,b.j  j=1:n  #
    ####################################
    si <- X%*%beta.old #n-sequence, dim=null
    
    if ((iter-1) %% 20 == 0) {
      hm <- KernSmooth::dpill(si, y) 
      h <- hm*(tau*(1-tau)/(stats::dnorm(stats::qnorm(tau)))^2)^.2 
    }

    x0 <- 0 
    for (j in 1:n) {
      x0 <- si[j] 
      fit <- clprq0(si, y, t, pp, h, tau, x0)
      a[j] <- fit$fv 
      b[j] <- fit$dv 
    }
    
    #############################
    # step 2: compute beta.new #
    #############################
    # here, let v.j=1/n 
    ynew <- y[ii] - a[jj]
    xnew <- b[jj]*diff
    kwts <- stats::dnorm(diff%*%beta.old/h)*t[jj]
    wts <- kwts*pwts
    
    #fit<-quantreg::rq(ynew ~ 0+ xnew, weights = wts, tau = tau, method="pfn")   #pfn for very large problems
    fit <- quantreg::rq(ynew ~ 0 + xnew, weights = wts, tau = tau, ci = FALSE)   #default: br method, for several thousand obs
    # 0, to exclude intercept
    beta.new <- fit$coef 
    beta.new <- sign(beta.new[1])*beta.new/sqrt(sum(beta.new^2))    #normalize
    
  } #end iterates over iter 
  
  
  flag.conv <- (iter < maxiter)
  
  beta <- beta.new 
  names(beta) <- colnames(X)
  
  si <- X%*%beta
  hm <- KernSmooth::dpill(si,y) 
  if (is.null(h)) {
    h <- hm*(tau*(1-tau)/(stats::dnorm(stats::qnorm(tau)))^2)^.2 
  } else {
    h <- h 
  }
  
  yhat <- rep(0,n) 
  for (j in 1:n) {
    x0 <- si[j]
    local.fit <- clprq0(si, y, t, pp, h, tau, x0) 
    yhat[j] <- local.fit$fv 
  }
  
  err <- y-yhat 
  R <- sum(abs(err)+(2*tau-1)*err)/n 
  
  csiqr_ojb <- list(beta=beta, flag.conv=flag.conv, X=X, y=y, yhat=yhat, tau=tau,
                    rqfit=fit, MSAE = R)
  
  class(csiqr_ojb) <- "csiqr"
  return(csiqr_ojb)
}


#' plot function of siqr
#'
#' @param x The SIQR model object
#' @param ... optional arguments
#' @param bootstrap.interval whether to calculate and plot bootstrap interval
#'
#' @return None
#' @export plot.siqr
#' @exportS3Method
plot.siqr <- function(x, ..., bootstrap.interval = FALSE) {
  si <- x$X%*%x$beta
  y <- x$y
  plot(si,y,xlab = "Single Index", ylab = "Predicted Y", col="gray", main="Fitted Quantile Plot") 
  graphics::lines(sort(si),x$yhat[order(si)],lty=1,lwd=1.5,col="red") 
  
  if(bootstrap.interval){
    tau <- x$tau
    hm <- KernSmooth::dpill(si,y)
    h <- hm*(tau*(1-tau)/(stats::dnorm(stats::qnorm(tau)))^2)^.2
    
    #get residual
    res <- y-x$yhat
    n <- length(res)
    
    #get bootstrap y.hat
    #v1
    # B=100
    # y.hat.B <- matrix(NA,length(y.B),B)
    # for(b in 1:B){
    #   #get residual bootstrap data
    #   bs.index<-sample(n,replace=T)
    #   res.B<-res[bs.index]
    #   y.B<-x$yhat+res.B
    #   fit.B <- siqr(y.B, X, beta.initial = beta0, tau=tau,maxiter = 20,tol = 1e-6, method = "Wu")
    #   y.hat.B[,b] <- fit.B$yhat
    # }
    
    
    #v2
    B=100
    y.hat.B <- matrix(NA,length(y),B)
    for(b in 1:B){
      for(i in 1:length(y)){
        #get residual bootstrap data
        bs.index<-sample(n,replace=T)
        res.B<-res[bs.index]
        y.B<-x$yhat+res.B
        fit.B <- clprq0(si, y.B, h, tau=tau, si[i])
        y.hat.B[i,b] <- fit.B$fv
      }
    }
    
    #get stats::sd of bootstrap Y.hat
    se.yhat <- apply(y.hat.B,1,stats::sd)
    #2*stats::sd +/- original y.hat to form the interval
    yhat.B.025 <- x$yhat - 2 * se.yhat
    yhat.B.975 <- x$yhat + 2 * se.yhat
    #plot
    #plot.si(x = x)
    graphics::lines(sort(si),yhat.B.025[order(si)],lty=6,lwd=1.5,col="blue")
    graphics::lines(sort(si),yhat.B.975[order(si)],lty=6,lwd=1.5,col="blue")
  }
}


#' Function to print summary
#'
#' @param object the single index quantile regression model object
#' @param digits controls digits in output
#' @param signif.stars whether show the significance stars
#' @param ... extra arguments
#'
#' @return the summarized information object
#' @exportS3Method
summary.siqr <- function(object, digits = max(5, getOption("digits") - 3),
                         signif.stars = getOption("show.signif.stars"), ...)
  
{
  cat("Tau: ",object$tau,"\n",sep="")
  
  if (length(object$beta)>0)
  { cat("\nsingle index coefficients:\n")
    stats::printCoefmat(data.frame(Coefficients=object$beta), digits = digits, 
                        signif.stars = signif.stars, na.print = "NA", ...)
  }
  cat("\n")
  cat("Model MSAE: ",object$MSAE,"\n",sep="")
  cat("Model convergence: ",object$flag.conv,"\n",sep="")
  invisible(object)
}
