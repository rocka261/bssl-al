
#  ------------------------------------------------------------------------

#  SSL test on two gaussians

#  ------------------------------------------------------------------------

require(ggplot2)
require(sampling)

if(!exists("calc_prob", mode="function")) source("C:/Users/erogeka/Documents/R/util.R")

x1 <- c(rnorm(100,-2, 1),rnorm(100,2, 1)) 
x2 <- rnorm(200, 0, 1)
x <- data.frame(x1,x2)
X <- model.matrix(~., x)

beta <- c(1,2,0)
levels <- c('class.1', 'class.2')
y <- apply(calc_prob(beta, X), 1, function(x) levels[rmultinom(1,1,x)==1])
y <- factor(y)
data <- data.frame(y,x)[order(y),]

ggplot(data, aes(x1,x2, color=y)) + geom_point()

# Removing some observations
# How many labeles in total?
nlab <- 10
pick <- strata(data, 'y', round(table(y)*(sum(table(y))-nlab)/sum(table(y))), method = 'srswor')$ID_unit
y.na <- y
y.na[pick] <- NA



# Modeling ----------------------------------------------------------------

# Labeled and unlabeled
ssl_model.UL.update <-  ssl_MH(y =  y.na, x = x, 
                               burn = 1000, samp = 2000,
                               kappa = 0.7, tau=10, scale = FALSE, initVal = NULL, 
                               adapt = FALSE, init.yu = 'multinom', sim.yu = TRUE)
mcmc.plots(ssl_model.UL.update$beta, pdf = FALSE,filename = 'trace-UL-update.pdf', title = 'Trace and density plots from MCMC simulation on both \nlabeled and unlabeled data when Yu is updated each iteration')
# Partial update
ssl_model.UL.update <-  ssl_MH(y =  y.na, x = x, 
                               burn = 1000, samp = 2000,
                               kappa = 0.1, tau=10, scale = FALSE, initVal = NULL, 
                               adapt = FALSE, init.yu = 'multinom', sim.yu = TRUE, 
                               sim.yu.partial = 150)
mcmc.plots(ssl_model.UL.update$beta, pdf = FALSE,filename = 'trace-UL-update.pdf', title = 'Trace and density plots from MCMC simulation on both \nlabeled and unlabeled data when Yu is updated each iteration')




apply(ssl_model.UL.update$beta, 2, mean)
pred.UL.upd <- predict.with.posterior(post.beta=ssl_model.UL.update$beta, X=X, y=y, levels=levels(y))
pred.UL.upd$a.m
ggplot(data, aes(x1, x2, color = factor(ssl_model.UL.update$ysim[,1]))) + geom_point()
ggplot(data, aes(x1, x2, color = factor(ssl_model.UL.update$ysim[,1000]))) + geom_point()
ggplot(data, aes(x1, x2, color = factor(ssl_model.UL.update$ysim[,2000]))) + geom_point()
ggplot(data, aes(x1, x2, color = factor(ssl_model.UL.update$ysim[,3000]))) + geom_point()
ggplot(data, aes(x1, x2, color = factor(ssl_model.UL.update$ysim[,4000]))) + geom_point()
ggplot(data, aes(x1, x2, color = factor(ssl_model.UL.update$ysim[,5000]))) + geom_point()
ggplot(data, aes(x1, x2, color = factor(ssl_model.UL.update$ysim[,6000]))) + geom_point()
ggplot(data, aes(x1, x2, color = factor(ssl_model.UL.update$ysim[,7000]))) + geom_point()
ggplot(data, aes(x1, x2, color = factor(ssl_model.UL.update$ysim[,8000]))) + geom_point()
ggplot(data, aes(x1, x2, color = factor(ssl_model.UL.update$ysim[,9000]))) + geom_point()
ggplot(data, aes(x1, x2, color = factor(ssl_model.UL.update$ysim[,10000]))) + geom_point()



ssl_model.UL.no.update <-  ssl_MH(y =  y.na, x = x, 
                                  burn = 5000, samp = 10000,
                                  kappa = 1, tau=10, scale = FALSE, initVal = NULL, 
                                  adapt = FALSE, init.yu = 'multinom', sim.yu = FALSE)
mcmc.plots(ssl_model.UL.no.update$beta, pdf = FALSE,filename = 'trace-UL-no-update.pdf', title = 'Trace and density plots from MCMC simulation on both labeled \nand unlabeled data when Yu not is updated each iteration')
apply(ssl_model.UL.update$beta, 2, mean)
pred.UL <- predict.with.posterior(post.beta=ssl_model.UL.no.update$beta, X=X, y=y, levels=levels(y))
pred.UL$a.m

# Labeled and unlabeled
ssl_model.L <-  ssl_MH(y =  y[!is.na(y)], x = x[!is.na(y),], 
                       burn = 5000, samp = 10000,
                       kappa = 1.1, tau=10, scale = FALSE, initVal = NULL, 
                       adapt = FALSE, init.yu = 'multinom', sim.yu = FALSE)

mcmc.plots(ssl_model.L$beta, pdf = FALSE,filename = 'trace-L.pdf',  title = 'Trace and density plots from MCMC simulation on only labeled data')
apply(ssl_model.L$beta, 2, mean)
pred.L <- predict.with.posterior(ssl_model.L$beta, X=X, y=y, levels=levels(y))
pred.L$a.m


# Reference
ssl_model.ref <-  ssl_MH(y =  y, x = x, 
                       burn = 5000, samp = 10000,
                       kappa = 1.2, tau=10, scale = FALSE, initVal = NULL, 
                       adapt = FALSE, init.yu = 'multinom', sim.yu = FALSE)

mcmc.plots(ssl_model.ref$beta, pdf = FALSE,filename = 'trace-L.pdf',  title = 'Trace and density plots from MCMC simulation on only labeled data')
apply(ssl_model.ref$beta, 2, mean)
pred.ref <- predict.with.posterior(ssl_model.ref$beta, X=X, y=y, levels=levels(y))
pred.ref$a.m









