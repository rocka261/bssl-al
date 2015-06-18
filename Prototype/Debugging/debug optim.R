
#  ------------------------------------------------------------------------

#  Functions for MCMC sampling with the Metrolopis Hastings algorithm
#  with part of the resonse with missing labels
#  The missing labels for the reonse will be simultaniously estimates
#  while simulating from the posterior of the parameters

#  Included are some evaluation and plotting functions

#  ------------------------------------------------------------------------

library(mvtnorm)
library(nnet)   # Multinomial Regression
library(caret)
library(GGally)
library(optimx)
library(ggplot2)
library(plyr)
require(gridExtra)
library(reshape2)
library(glmnet)
library(Matrix)
require(foreign)
require(nnet)
require(reshape2)



# DATA --------------------------------------------------------------------

ml <- read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")



y <- ml$prog
x <- ml[,c('ses' , 'write')]
X <- model.matrix(~.,x)
y.mm <- model.matrix(~y-1,y)




# Functions ---------------------------------------------------------------



calc_prob <- function(beta, X){
  # Calcs the probability of the class y_j for j in 1:J
  # beta: dim = c(q, (J-1)) 
  # X: dim = c(n, q)
  exp.eta <-  t(exp(t(beta) %*% t(X)))
  exp.eta[exp.eta == Inf] <- 10^100
  p1 <- 1/(1+rowSums(exp.eta))
  p <- cbind(p1, p1*exp.eta)
  return(p)
}
logPost <- function(betaVect, y, X, mu, Sigma,...){
  # Calculates the log posterior of a multinomial logistic regression
  # given a y vector with dim = c(n, J)
  # X design matrix with dim = c(n,q)
  # And a beta vector with dim = c((q, (J-1))
  q <- ncol(X)
  J <- ncol(y)
  p <- calc_prob(matrix(betaVect, q), X)
  logLik <- sum(y * log(p))
  if (abs(logLik) == Inf | is.na(logLik)) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(betaVect, mu, Sigma, log=TRUE)
  return(logLik + logPrior)
}


logPostx <- function(betaVect, opt = list(y=y.mm,X=X,mu=mu,Sigma=Sigma)){
  # Calculates the log posterior of a multinomial logistic regression
  # given a y vector with dim = c(n, J)
  # X design matrix with dim = c(n,q)
  # And a beta vector with dim = c((q, (J-1))
  y <- opt$y
  X <- opt$X
  mu <- opt$mu
  Sigma <- opt$Sigma
  q <- ncol(X)
  J <- ncol(y)
  p <- calc_prob(matrix(betaVect, q), X)
  logLik <- sum(y * log(p))
  if (abs(logLik) == Inf | is.na(logLik)) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(betaVect, mu, Sigma, log=TRUE)
  return(logLik + logPrior)
}


lp.gr <- function(betaVect, y, X, mu, Sigma, ...){
  # Gives the gradient of the log posterior of the the multinomial logistic function
  # betaVect: dim = c((J-1)*q, 1)
  # y: dim = c(n, J)
  # X: dim = c(n, q)
  J <- ncol(y)
  q <- ncol(X)
  p <- calc_prob(matrix(betaVect, q), X)
  likGradient <- as.vector(t(X) %*% (y[ ,2:J] - p[ ,2:J]))
  priorGradient <- -betaVect %*% solve(Sigma)   # Correct?
  return(likGradient + priorGradient)
}

lp.grx <- function(betaVect, opt = list(y=y.mm,X=X,mu=mu,Sigma=Sigma)){
  # Gives the gradient of the log posterior of the the multinomial logistic function
  # betaVect: dim = c((J-1)*q, 1)
  # y: dim = c(n, J)
  # X: dim = c(n, q)
  y <- opt$y
  X <- opt$X
  mu <- opt$mu
  Sigma <- opt$Sigma
  J <- ncol(y)
  q <- ncol(X)
  p <- calc_prob(matrix(betaVect, q), X)
  likGradient <- as.vector(t(X) %*% (y[ ,2:J] - p[ ,2:J]))
  priorGradient <- -betaVect %*% solve(Sigma)   # Correct?
  return(likGradient + priorGradient)
}


# Gives the hessian for the logposterior of the multinomial logistic regression
lp.hesx <- function(betaVect, opt = list(y=y.mm,X=X,mu=mu,Sigma=Sigma)){
  # betaVect: dim = c((J-1)*q, 1)
  # y: dim = c(n, J)
  # X: dim = c(n, q)
  y <- opt$y
  X <- opt$X
  mu <- opt$mu
  Sigma <- opt$Sigma
  p <- calc_prob(beta = matrix(betaVect, q), X)[,-1]   # Leave the first column out
  J <- ncol(y)
  xbar <- bdiag(rep(list(X), (J-1)))
  wbar <- Matrix(0,ncol=(J-1)*n, nrow=(J-1)*n)
  if(J>2){
    for(k in 1:(J-1)){
      for(t in 1:(J-1)){
        if(k==t){
          wbar[(1+(k-1)*n):(k*n),(1+(t-1)*n):(t*n)] <- diag(p[,k]*(1-p[,t]))
        }
        else{
          wbar[(1+(k-1)*n):(k*n),(1+(t-1)*n):(t*n)]  <- diag(p[,k]*(0-p[,t]))
        }
      }
    }
  }
  else{
    wbar <- diag(p*(1-p)) 
  }
  hessian <- -t(xbar) %*% wbar %*% xbar - solve(Sigma)
  return(as.matrix(hessian))
}

lp.hesx(betaVect=as.vector(initVal), opt = list(y=y.mm,X=X,mu=mu,Sigma=Sigma))
betaVect=as.vector(initVal)
opt = list(y=y.mm,X=X,mu=mu,Sigma=Sigma)
y <- opt$y
X <- opt$X
mu <- opt$mu
Sigma <- opt$Sigma
p <- calc_prob(beta = matrix(betaVect, q), X)[,-1]
mult <- array(0, dim=c(4,4,nrow(p)))
for(i in 1:nrow(p)){
  mult[,,i] <- -p[i,1]*(1-p[i,1])*X[i,]%*%t(X[i,])
}

apply(mult, c(1,2), sum) - solve(diag(tau^2,4))


mult <- array(0, dim=c(4,4,nrow(p)))
for(i in 1:nrow(p)){
  mult[,,i] <- p[i,1]*(p[i,2])*X[i,]%*%t(X[i,])
}
apply(mult, c(1,2), sum) - solve(diag(tau^2,4))



# Test --------------------------------------------------------------------

q <- ncol(X)
J <- ncol(y.mm)
n <- nrow(X)
tau <- 10
# Priors 
mu <- rep(0, (J-1)*q)
Sigma <- tau^2*diag(q*(J-1))

initVal <- solve(crossprod(X,X))%*%t(X)%*%y.mm[,2:J]


# Test optim and optimx
OptimResults<-optim(as.vector(initVal),logPost,gr=lp.gr ,
                    y.mm,X,mu,Sigma,
                    method=c("BFGS"),control=list(maxit = 30000,fnscale=-1),hessian=TRUE)

optimxResults <- optimx(as.vector(initVal),fn=logPostx, gr=lp.grx, hess = lp.hesx, 
                        method = c("BFGS"), hessian = TRUE, lower=-Inf, upper=Inf,
                        control = list(fnscale=-1),opt = list(y=y.mm,X=X,mu=mu,Sigma=Sigma))
OptimResults$hessian
attr(optimxResults,'details')[[3]]
