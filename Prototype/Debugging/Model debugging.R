
# Test pred ---------------------------------------------------------------


# Load functions from util.R script
if(!exists("SSL_MH", mode="function")) source("C:/Users/erogeka/Documents/R/util.R")


require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)
require(sampling)

ml <- read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")

ml$prog2 <- relevel(ml$prog, ref = "academic")
test <- multinom(prog2 ~ ses + write, data = ml)

head(pp <- fitted(test))

p <- calc_prob(t(coef(test)), model.matrix(prog2 ~ ses + write, data = ml))
p ==pp
any(abs(p-pp ) > 10e-10)   # Diff but margin small

# Test with bayes ---------------------------------------------------------

y <- ml$prog2
x <- ml[,c('ses' , 'write')]

ssl.model <-  ssl_MH(y = y, x= x, burn = 1000, samp = 2000, kappa = 0.6, 
                        tau=10, scale = FALSE, initVal = NULL, 
                        adapt = FALSE, init.yu = 'none')
mcmc.plots(ssl.model$beta, pdf = FALSE)



t(apply(ssl.model$beta, c(2,3), mean))
coef(test)

mcmc.boxplots(beta.draws = ssl.model$beta, pdf = FALSE)
mcmc.plots(beta.draws = ssl.model$beta, pdf = FALSE)


pred <- predict.with.posterior(ssl.model$beta, X = model.matrix(~., x), y = y, levels = levels(y))
pred$a.m

sum(predict(test)  == y) / length(y)


# Verkar fungera som det ska.



# Iris data ---------------------------------------------------------------

data(iris)
x <- iris[,1:4]
y <- iris[,5]
# Removing some observations
# How many labelein total?
nlab <- 10
pick <- strata(data.frame(y, x), 'y', 
               round(table(y)*(sum(table(y))-nlab)/sum(table(y))), method = 'srswor')$ID_unit
y.na <- y
y.na[pick] <- NA


## multinom()
mmlogit <- multinom(Species ~., data = iris)


ssl.model.ref <-  ssl_MH(y = y, x= x, burn = 2000, samp = 4000, kappa = 0.6, 
                     tau=20, scale = FALSE, initVal = NULL, 
                     adapt = FALSE, init.yu = 'none')


ssl.model.U <-  ssl_MH(y = y.na, x= x, burn = 20000, samp = 10000, kappa = 0.5, 
                     tau=20, scale = FALSE, initVal = NULL, 
                     adapt = FALSE, init.yu = 'multinom', sim.yu  = TRUE)

mcmc.plots(beta.draws = ssl.model.U$beta, pdf = FALSE)

ssl.model.L <-  ssl_MH(y = y[!is.na(y.na)], x= x[!is.na(y.na),], burn = 10000, samp = 10000, kappa = 0.5, 
                       tau=20, scale = FALSE, initVal = NULL, 
                       adapt = FALSE, init.yu = 'multinom', sim.yu  = FALSE)

mcmc.plots(beta.draws = ssl.model.L$beta, pdf = FALSE)

t(apply(ssl.model$beta, c(2,3), mean))
t(apply(ssl.model.U$beta, c(2,3), mean))
t(apply(ssl.model.L$beta, c(2,3), mean))

coef(mmlogit)

mcmc.boxplots(beta.draws = ssl.model$beta, pdf = FALSE)
mcmc.plots(beta.draws = ssl.model$beta, pdf = FALSE)



pred <- predict.with.posterior(ssl.model$beta, X = model.matrix(~., x), y = y, levels=levels(y))
pred$a.m
pred <- predict.with.posterior(ssl.model.U$beta, X = model.matrix(~., x), y = y, levels=levels(y))
pred$a.m
pred <- predict.with.posterior(ssl.model.L$beta, X = model.matrix(~., x), y = y, levels=levels(y))
pred$a.m




sum(predict(mmlogit)  == y) / length(y)



