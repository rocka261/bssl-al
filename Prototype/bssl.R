
#  ------------------------------------------------------------------------

#  Functions for MCMC sampling with the Metrolopis Hastings algorithm
#  with part of the resonse with missing labels
#  The missing labels for the response will be simultaniously estimated
#  while simulating from the posterior of the parameters

#  Includes some evaluation and plotting functions

#  ------------------------------------------------------------------------

# Load user defined config file
source('./config.R')


# Add local path to user installed packages
.libPaths( c( .libPaths(), localLibs))

library(mvtnorm)
library(nnet)   
# library(optimx)   # For optimization with the hessian
library(ggplot2)
library(plyr)
require(gridExtra)
library(reshape2)
library(glmnet)
library(Matrix)



calc_prob <- function(beta, X){
  #  Calcs the probability of the classes y = j for j in 1:J
  #  beta: dim = c(q, (J-1)) 
  #  X: dim = c(n, q)
  exp.eta <-  t(exp(t(beta) %*% t(X)))
  exp.eta[exp.eta == Inf] <- 10^100   # To handle Inf values
  p1 <- 1/(1+rowSums(exp.eta))
  p <- cbind(p1, p1*exp.eta)
  return(p)
}


logPost <- function(betaVect, y, X, mu, Sigma,...){
  #   Calculates the log posterior predictive distribution of a 
  # multinomial logistic regression model given a y vector 
  # with dim = c(n, J)
  #   X design matrix with dim = c(n,q)
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
  #   This is just a version of logPost() where the input parameters are given
  # in a list. Suits optimx() better than optim() for optimization using the 
  # a given hessian as input. 
  #   X design matrix with dim = c(n,q)
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




init <- function(X, y, initVal=NULL, mu, Sigma){
  # Init the parameter values with the OLS estimate, if the X-matrix sufficient
  # otherwise with the multinom() function from the nnet package
  tryCatch(
{
  y.mm <- model.matrix(~y-1, y)
  J <- ncol(y.mm)
  q <- ncol(X)
  if(is.null(initVal)){
    initVal <- solve(crossprod(X,X))%*%t(X)%*%y.mm[,2:J]
  }
  OptimResults<-optim(as.vector(initVal),logPost,gr=lp.gr ,
                      y.mm,X,mu,Sigma,
                      method=c("BFGS"),control=list(maxit = 30000,fnscale=-1),hessian=TRUE)
  inh <- -solve(OptimResults$hessian)
  betaVect <- OptimResults$par 
  conv <- OptimResults$convergence
  return(list(betaVect=betaVect, inh=inh, conv=conv))   # Initial values by OLS
  
#   optimResults <- optimx(as.vector(initVal),fn=logPostx, gr=lp.grx, hess = lp.hesx, 
#                          method = c("BFGS"), hessian = TRUE, lower=-Inf, upper=Inf,
#                          control = list(fnscale=-1),opt = list(y=y.mm,X=X,mu=mu,Sigma=Sigma))
#   inh <- forceSymmetric(-solve(attr(optimResults, 'details')[[3]]))
#   betaVect <- optimResults[1:((J-1)*q)]
  return(list(betaVect=betaVect, inh=inh))   # Initial values by OLS
}
,error=function(cond){
  message('An exception has been made due to an error:')
  message(cond)
  initVal <- coef(multinom(y~X[,-1]))    # Init with multinom estimate
  if(J>3) initVal <- t(initVal)
  OptimResults<-optim(as.vector(initVal),logPost,gr=lp.gr ,
                      y.mm,X,mu,Sigma,
                      method=c("BFGS"),control=list(maxit = 30000,fnscale=-1),hessian=TRUE)
  inh <- -solve(OptimResults$hessian)
  betaVect <- OptimResults$par 
  conv <- OptimResults$convergence
  return(list(betaVect=betaVect, inh=inh, conv=conv))
# 
#   optimResults <- optimx(as.vector(initVal),fn=logPostx, gr=lp.grx, hess = lp.hesx, 
#                          method = c("BFGS"), hessian = TRUE, lower=-Inf, upper=Inf,
#                          control = list(fnscale=-1),opt = list(y=y.mm,X=X,mu=mu,Sigma=Sigma))
#   inh <- -solve(attr(optimResults, 'details')[[3]])
#   betaVect <- optimResults[1:((J-1)*q)]
  return(list(betaVect=betaVect, inh=inh))
}
,finally=function(cond){
  message('This is the origional message')
  message(cond)
} 
  )
}



ssl_MH <- function(y, x,                       # Data
                   burn = 500, samp = 1000,    # Sampling
                   kappa = 0.3, tau = 10,     # Constant for acceptance rate and prior scale
                   proposal = NULL,            # An optional list containing a betaVect and inh
                   scale = FALSE, initVal = 'multinom', adapt = FALSE,  # Other control parameters
                   init.yu = 'none', sim.yu = TRUE, sim.yu.partial=NULL)              
{
  #   Metropolice hastings sampler for multinomial logistic regression with 
  #   possability to use unlabeled instances
  #  scale=TRUE scales the covariates to have mean=0 and variance=1 
  #  adapt=TRUE reestimates the proposal distribution after half of the samples 
  # init.yu: takes values 'none' meaning no inital values are estimated for yu for the proposal distribution
  # 'rand' gives a label to yu by random
  # 'multinom' mean that the labels are estimated from the multinomial function in the nnet package
  # 'ridge' means that the labels for yu are estimated with ridge regression
  # 'input' mean that proposal is provided as input, and not needs to be estimated
  #  sim.yu=TRUE lets the unlabeled instances to be estimated
  #  sim.yu.partial=<some number> lets yu be updated at fewer times each mcmc iteration 
  
  n <- length(y)
  y <- as.factor(y)
  labels <- levels(y)
  J <- nlevels(y)
  if(scale) {
    X <- model.matrix(~., data.frame(scale(x)))
  }else{
    X <- model.matrix(~., data.frame(x))
  }
  q <- ncol(X)
  
  
  isna <- is.na(y)   # Instances which labels has to be estimated
  yl <- y[!isna]
  yu <- y[isna]
  Xl <- X[!isna,]
  Xu <- X[isna,]
  
  # Priors 
  mu <- rep(0, (J-1)*q)
  Sigma <- tau^2*diag(q*(J-1))
  
  #### Proposal
  y.temp <- y  
  
  # Inital values for yu - uniform prior
  if(init.yu == 'rand'){
    if(sum(isna)>0){
      y.temp[isna] <- replicate(length(yu), labels[rmultinom(1,1,rep(1, J) * 1/J)==1])
    }
    proposal <- init(X=X, y=y.temp, initVal = initVal, mu=mu, Sigma=Sigma)  
  }
  # Inital values for yu - prediction with multinom
  if(init.yu == 'multinom'){
    if(sum(isna)>0){
      mmodel <- multinom(y.temp~ X[,-1]) 
      y.temp[isna] <- apply(calc_prob(matrix(coef(mmodel), q),X[isna,]), 1, 
                            function(x) labels[rmultinom(1,1,x)==1])
    }
    proposal <- init(X=X, y=y.temp, initVal = initVal, mu=mu, Sigma=Sigma)  
    
  }
  # Initial values yu - prediction with ridge regression
  if(init.yu == 'ridge'){
    if(sum(isna)>0)   {
      ridge <- cv.glmnet(Xl, yl, family = "multinomial", alpha = 0)
      y.temp[isna] <- predict(ridge, newx=Xu, type = "class")
    }
    proposal <- init(X=X, y=y.temp, initVal = initVal, mu=mu, Sigma=Sigma)
    
  }
  # Proposal just based on the known labels
  if(init.yu == 'none'){
    proposal <- init(X=Xl, y=yl, initVal = initVal, mu=mu, Sigma=Sigma)
    # Update y with the chosen beta
    if(sum(isna)>0) {
      y.temp[isna] <- apply(calc_prob(matrix(proposal$betaVect,q),X[isna,]), 1, 
                              function(x) labels[rmultinom(1,1,x)==1])
    }
  }
  
  betaVect <- as.vector(t(proposal$betaVect))
  inh <- proposal$inh 
  
#   if(sum(isna)>1){
#     # Update y with the chosen beta
#     y.temp[isna] <- apply(calc_prob(matrix(betaVect, q),X[isna,]), 1, 
#                           function(x) labels[rmultinom(1,1,x)==1])
#   }


  # Additional Vectors and storing 
  beta.new <- mu  
  nsim <- burn + samp 
  beta.store <- array(0, dim=c(nsim,q, J-1))
  dimnames(beta.store) <- list(NULL, colnames(X), labels[-1])
  beta.store[1, , ] <- matrix(betaVect, q)
  if(!is.null(sim.yu.partial)) update.points <- round(seq(2,nsim, length.out = sim.yu.partial))
  y.sim <- matrix(0, n, nsim)
  A<-0   # Counter for accepted draws
  A.ad <- 0   # Counter for accepted draws after adaption
  adapt.made = FALSE
  update <- 0
  
  
  # MCMC
  start <- Sys.time()
  for (i in 2:nsim) {
    
    # Update y model matrix
    y.model.matrix <- model.matrix(~y.temp-1, y.temp)
    
    # Draw Candidates by random walk
    betaVect.cand <- rmvnorm(1, mean=betaVect, sigma=kappa^2*inh)
    #betaVect.cand <- betaVect + rmv(1,sigma= kappa*cov, q)    # By multivariate t-distribution. How many df are reasonable?
    
    # Acceptance prob on log scale
    r <- logPost(betaVect.cand, y.model.matrix, X, mu, Sigma) - logPost(betaVect, y.model.matrix, X, mu, Sigma)
    
    # Accept with a probability
    if(log(runif(1))<r){
      betaVect <- betaVect.cand
      if (i>burn) A<-A+1
      if(adapt){
        if (i>(burn + round(samp/2))) A.ad<-A.ad+1
      }
    }
    # Update yu
    if(sum(isna)>0 & sim.yu){
      # For partial update
      if(!is.null(sim.yu.partial)){
        if(i %in% update.points){
          # Update y with the chosen beta
          update <- update + 1
          # Update with beta draws since last yu update
          if(update>1){
            pred.y <- predict.with.posterior(beta.store[update.points[update-1]:i,,], 
                                             X[isna,], levels=labels)
            y.temp[isna] <- apply(pred.y$prob, 1, 
                                  function(x) labels[rmultinom(1,1,x)==1])
          }else{
            y.temp[isna] <- apply(calc_prob(matrix(betaVect, q),X[isna,]), 1, 
                                  function(x) labels[rmultinom(1,1,x)==1])
          }
        }
      # For update every iteration  
      }else{
        # Update y with the chosen beta
        y.temp[isna] <- apply(calc_prob(matrix(betaVect, q),X[isna,]), 1, 
                              function(x) labels[rmultinom(1,1,x)==1])
      }
    }
    
    
    # Adaptions step
    if(adapt){
      adapt.made = FALSE
      if(i == burn + round(samp/2)){
        if(A/(samp/2) > 0.5 | A/(samp/2) < 0.1){
          kappa = ifelse(A/(samp/2) > 0.5, kappa*2, kappa/2)
          adapt.made = TRUE
        }
      }
      if (i>burn + round(samp/2)) A.ad<-A.ad+1
    }
    
    # Store Results
    beta.store[i, , ] <- matrix(betaVect, q)
    y.sim[,i] <- y.temp
    if (i%%100==0) print(i)
    
  } # End MCMC
  print(Sys.time() - start)
  print(paste('Acceptance rate:', A/samp))
  if(adapt.made) {
    print('Adaption mode during sample')
    print(paste('Acceptance rate after adaption:', A.ad/(burn + round(samp/2))))
  }
  return(list(beta = beta.store[(burn+1):nsim, , ], 
              ysim = y.sim[, (burn + 1):nsim], 
              isna = isna, a.rate = A/samp, levels = labels))
}





# Evaluation --------------------------------------------------------------



# For prediction with the posterior sample on new X
predict.with.posterior <- function(post.beta, X, levels, y=NULL){
  # Post beta should have dimensions: dim=c(nSim, q, J-1)
  # X: dim=c(n, q)
  # y: dim=c(1): should be a factor with J levels
  n <- nrow(X)
  J <- length(levels)
  nsim <- dim(post.beta)[1]
  y.pred <- matrix(0, n, nsim)
  p <- apply(post.beta, 1, function(x) calc_prob(x, X))
  if(is.null(dim(post.beta)[3])| is.na(dim(post.beta)[3])) {
    for(i in 1:nsim){
      .p <- matrix(p[,i], ncol = J)
      y.pred[,i] <- apply(.p, 1, function(x) levels[rmultinom(1,1,x)==1])
    }
  }else{
    for(i in 1:nsim){
      .p <- matrix(p[,i], ncol = J)
      y.pred[,i] <- apply(.p, 1, function(x) levels[rmultinom(1,1,x)==1])
    }
  }
  
  y.freq <- t(apply(y.pred, 1, function(x) table(x)))
  if(typeof(y.freq)=='list') y.freq<-ldply(y.freq, rbind)
  y.freq[is.na(y.freq )] <- 0
  y.prob <- t(apply(y.freq, 1, function(x) x/sum(x)))
  entropy <- -y.prob*log(y.prob)
  entropy[is.na(entropy)] <- 0
  row.entropy <- rowSums(entropy)
  y.pred.prob <- apply(y.prob, 1, function(x) levels[rmultinom(1,1,x)==1])
  acc.prob <- sum(y.pred.prob == y)/n
  y.pred.max <- apply(y.prob, 1, function(x) levels[which.max(x)])
  acc.max <- sum(y.pred.max == y)/n
  return(list(freq=y.freq, prob=y.prob, pred=y.pred, entr=row.entropy, 
              a.p = acc.prob, a.m = acc.max))
}




# Plotting ----------------------------------------------------------------

mcmc.boxplots <- function(beta.draws, title=NULL, pdf = FALSE, filename = 'mcmcmboxplot.pdf'){
  # beta.draws: dim=c(samp, q, J-1)
  melted <- melt(beta.draws, row.names = 0)
  if(is.na(dim(beta.draws)[3])){
    melted[,2] <- as.factor(melted[,2])
    colnames(melted) <- c('Iteration', 'Covariate',  'Coefficient')  
    plots <- ggplot(melted, aes(y=Coefficient, x=Covariate)) + geom_boxplot() + ggtitle(title) + theme_bw() +
      theme(plot.title = element_text(size=20)) + coord_flip() + 
      scale_x_discrete(limits = rev(dimnames(beta.draws)[[2]]))
  }else{
    melted[,2:3] <- apply(melted[,2:3], 2, function(x) as.factor(x))
    colnames(melted) <- c('Iteration', 'Covariate', 'Class', 'Coefficient')  
    plots <- ggplot(melted, aes(y=Coefficient, x=Covariate, fill=Class)) + geom_boxplot() + ggtitle(title) + theme_bw() +
      theme(plot.title = element_text(size=20))+ coord_flip() + 
      scale_x_discrete(limits = rev(dimnames(beta.draws)[[2]]))
  }
  if(pdf)
  {
    ggsave(filename, plots)
  }else{
    plots
  }
}


mcmc.plots <- function(beta.draws, title = NULL, pdf=FALSE, filename = 'mcmcplots.pdf'){
  # beta.draws: dim=c(samp, q, J-1)
  q <- dim(beta.draws)[2]
  melted <- melt(beta.draws, row.names = 0)
  coefs <- dimnames(beta.draws)[[2]]
  # Automaically prints to a pdf if the number of covariates is higher than 6 
  if (q > 6){
    pages <- ceiling(q/6)
    coef.split <- split(coefs, ceiling(seq_along(coefs)/6))
    pdf(filename)
    for(page in 1:pages){
      if(is.na(dim(beta.draws)[3])){
        colnames(melted) <- c('Iteration', 'Covariate', 'Coefficient')
        melted.page <- melted[melted$Covariate %in% coef.split[[page]],]
        melted.page[,2] <- as.factor(melted.page[,2])
        t <- ggplot(melted.page, aes(y=Coefficient, x=Iteration)) +
          geom_line() +  facet_grid(Covariate~. , scales = 'free') + theme_bw() + theme(legend.position="none")
        d <- ggplot(melted.page, aes(x=Coefficient)) + 
          geom_density(alpha=.5) + facet_grid(Covariate~. , scales = 'free') + theme_bw() + theme(legend.position="none")
      }
      else{
        colnames(melted) <- c('Iteration', 'Covariate', 'Class', 'Coefficient') 
        melted.page <- melted[melted$Covariate %in% coef.split[[page]],]
        melted.page[,2:3] <- apply(melted.page[,2:3], 2, function(x) as.factor(x))
        t <- ggplot(melted.page, aes(y=Coefficient, x=Iteration, color=Class)) +
          geom_line() +  facet_grid(Covariate~. , scales = 'free') + theme_bw() + theme(legend.position="none")
        d <- ggplot(melted.page, aes(x=Coefficient, fill=Class)) + 
          geom_density(alpha=.5) + facet_grid(Covariate~. , scales = 'free') + theme_bw() + theme(legend.position="none")
      }
      do.call("grid.arrange", c(list(t, d), list(ncol = 2, main=textGrob(title,gp=gpar(fontsize=16)))))
    }
    dev.off()
  }
  else{
    if(is.na(dim(beta.draws)[3])){
      melted[,2] <- as.factor(melted[,2])
      colnames(melted) <- c('Iteration', 'Covariate', 'Coefficient')  
      t <- ggplot(melted, aes(y=Coefficient, x=Iteration)) +
        geom_line() +  facet_grid(Covariate~. , scales = 'free') + theme_bw() + theme(legend.position="none")
      d <- ggplot(melted, aes(x=Coefficient)) + 
        geom_density(alpha=.5) + facet_grid(Covariate~. , scales = 'free') + theme_bw() + theme(legend.position="none")
    }
    else{
      melted[,2:3] <- apply(melted[,2:3], 2, function(x) as.factor(x))
      colnames(melted) <- c('Iteration', 'Covariate', 'Class', 'Coefficient') 
      t <- ggplot(melted, aes(y=Coefficient, x=Iteration, color=Class)) +
        geom_line() +  facet_grid(Covariate~. , scales = 'free') + theme_bw() + theme(legend.position="none")
      d <- ggplot(melted, aes(x=Coefficient, fill=Class)) + 
        geom_density(alpha=.5) + facet_grid(Covariate~. , scales = 'free') + theme_bw() + theme(legend.position="none")
    }
    
    do.call("grid.arrange", c(list(t, d), list(ncol = 2, main=textGrob(title,gp=gpar(fontsize=16)))))
  }  
}





#### Taken from Winston Chang
# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper functions

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}






