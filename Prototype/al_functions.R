
#  ------------------------------------------------------------------------

#  Different active learning methods

#  ------------------------------------------------------------------------

# Load user defined config file
source('./config.R')


# Add local path to user installed packages
.libPaths( c( .libPaths(), localLibs))


# Load functions from helpers.R script
if(!exists("entr", mode="function")) source("./helpers.R")


# Passive  ----------------------------------------------------------------

# Random query
random <- function(y.p,...){
  sample(1:nrow(y.p), 1)  
}



# Uncertainty strategies --------------------------------------------------

### Least confident
LC <- function(y.p,...){
  max.prob <- apply(y.p,1,max)
  return(which.min(max.prob))
}
# LC(y.p)


### Margin
marg <- function(y.p,...){
  max.prob <- apply(y.p,1,max)
  second.max.prob <- apply(y.p,1,function(x) max(x[-which.max(x)]))
  marg <- max.prob - second.max.prob
  return(which.min(marg))
}
# marg(y.p)


### Entropy
maxEntropy <- function(y.p,...){
 rowEntropy <- apply(y.p, 1, entr) 
 return(which.max(rowEntropy))
}
# maxEntropy(y.p)



# Expected model change ---------------------------------------------------

lp.gr <- function(betaVect, y, X, mu, Sigma, ...){
  # Gives the gradient of the log posterior of the the multinomial logistic function
  # betaVect: dim = c((J-1)*q, 1)
  # y: dim = c(n, J)
  # X: dim = c(n, q)
  if(!exists('calc_prob', mode="function")) source("./bssl.R")
  J <- ncol(y)
  q <- ncol(X)
  p <- calc_prob(matrix(betaVect, q), X)
  likGradient <- as.vector(t(X) %*% (y[ ,2:J] - p[ ,2:J]))
  priorGradient <- -betaVect %*% solve(Sigma)   # Correct?
  return(likGradient + priorGradient)
}

ml.gr <- function(betaVect, y, X,...){
  # Gives the gradient of the log posterior of the the multinomial logistic function
  # betaVect: dim = c((J-1)*q, 1)
  # y: dim = c(n, J)
  # X: dim = c(n, q)
  if(!exists('calc_prob', mode="function")) source("./bssl.R")
  J <- ncol(y)
  q <- ncol(X)
  p <- calc_prob(matrix(betaVect, q), X)
  likGradient <- as.vector(t(X) %*% (y[ ,2:J] - p[ ,2:J]))
  return(likGradient)
}

norm_vec <- function(x) sqrt(sum(x^2))

exp.model.change.post <- function(yl.mm, Xl, Xu,  tau = 10, postMean.betaVect,...){
  # yl.mm: dim=c(nl, J)
  # Xl: dim=c(nl, q)
  # Xu: dim=c(nu,q)
  J <- ncol(yl.mm)
  q <- ncol(Xl)
  nu <- nrow(Xu)
  mu <- rep(0, (J-1)*q)
  Sigma <- tau^2*diag(q*(J-1))
  if(!exists('calc_prob', mode="function")) source("./bssl.R")

  gr0 <- lp.gr(postMean.betaVect, yl.mm, Xl, mu, Sigma)  
  gr.change <- matrix(0, nrow=nu, ncol=J)
  is.y <- diag(J)
  for(i in 1:nu){
    for(j in 1:J){
      gr.change[i,j] <- y.p[i,j] * norm_vec(lp.gr(postMean.betaVect, rbind(yl.mm,is.y[j,]), rbind(Xl, Xu[i,]), mu, Sigma))
    }
  }
  which.max(rowSums(gr.change))
}


exp.model.change.ml <- function(y.p, y, x, param, ...){
  # yl.mm: dim=c(nl, J)
  # Xl: dim=c(nl, q)
  # Xu: dim=c(nu,q)
  yl <- y[!is.na(y)]
  yl.mm <- model.matrix(~yl-1,yl)
  X <- model.matrix(~., data.frame(x))
  Xl <- X[!is.na(y),]
  Xu <- X[is.na(y),]
  J <- ncol(yl.mm)
  q <- ncol(Xl)
  nu <- nrow(Xu)
  if(!exists('ml.gr', mode="function")) source("./bssl.R")
  
  gr0 <- ml.gr(param, yl.mm, Xl)  
  gr.change <- matrix(0, nrow=nu, ncol=J)
  is.y <- diag(J)
  for(i in 1:nu){
    for(j in 1:J){
      gr.change[i,j] <- y.p[i,j] * norm_vec(ml.gr(param, rbind(yl.mm,is.y[j,]), rbind(Xl, Xu[i,])))
    }
  }
  which.max(rowSums(gr.change))
}



# Expected error reduction ------------------------------------------------

# För alla klasser J
# För alla observationer i valideringssetted V

# Intractable to compute?
# Use ml-version of multilogit
exp.error.reduction <- function(y.p, y, x, param, loss='0-1',...){
  # yl.mm: dim=c(nl, J)
  # Xl: dim=c(nl, q)
  # Xu: dim=c(nu,q)
  if(!exists('calc_prob', mode="function")) source("./bssl.R")
  yl <- y[!is.na(y)]
  yl.mm <- model.matrix(~yl-1,yl)
  X <- model.matrix(~., data.frame(x))
  Xl <- X[!is.na(y),]
  Xu <- X[is.na(y),]
  J <- ncol(yl.mm)
  q <- ncol(Xl)
  nu <- nrow(Xu)
  nl <- nrow(Xl)
  labeled <- data.frame(yl, Xl[,-1])
  errors <- matrix(0, nrow=nu, ncol=J)
  levels <- levels(yl)
  for(i in 1:nu){
    for(j in 1:J){
      model <- multinom(yl~., rbind(labeled,data.frame(yl=levels[j],t(Xu[i,-1]))))
      if(loss == '0-1'){
        # 0-1 loss
        if(J>2){
          maxpred <- apply(calc_prob(t(coef(model)), Xu[-i,]),1, max)
          errors[i,j] <- y.p[i,j] * sum(1-maxpred)
        }
        else{
          maxpred <- apply(calc_prob(coef(model), Xu[-i,]),1, max)
          errors[i,j] <- y.p[i,j] * sum(1-maxpred)
        }
        }
        
      if(loss == 'log'){
        # log-loss
        errors[i,j] <- y.p[i,j] * sum(apply(calc_prob(matrix(coef(model),q), Xu[-i,]),1, entr))  
      }
    }
  }
  which.min(rowSums(errors))
}

# # start <- proc.time()
# exp.error.reduction(y.p = ssl.pred$prob, y = y, x =  x, 
#                     param=apply(ssl_model.UL.update$beta, 2, mean), loss='0-1')
# proc.time() - start
# # 0-1: 378
# 
# model <- multinom(yl~.,labeled)
# predict(model, Xu)
# colnames(testdata)


