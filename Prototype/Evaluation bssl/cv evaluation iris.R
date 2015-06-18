
#  ------------------------------------------------------------------------

#  Cross validation of SSL-MH classification on different amount of
#  labeled data
#  On newsgroups

#  ------------------------------------------------------------------------

library(sampling)
library(topicmodels)
library(nnet)

# Load functions from util.R script
if(!exists("SSL_MH", mode="function")) source("C:/Users/erogeka/Documents/R/util.R")


# Data --------------------------------------------------------------------

data(iris)
n <- nrow(iris)
x <- iris[,1:4]
X <- model.matrix(~., x) 
X <- as.matrix(X)
q <- ncol(X)
# Response
y <- iris[,5]
J <- nlevels(y)

testmodel <- multinom(y~X[,-1])
summary(testmodel)
sum(predict(testmodel) == y) / length(y)


# Iterate through different proportions of missing data -------------------


#################  User settings: #########################################
# Missing
# missing <- 1-c(.1, .25, .5, 0.75)   # Ratio of missing labels tested
nlabels <- c(10, 25, 50, 75, 100)
# CV
folds <- 10   # Number of cross validation sets
set.seed(987654321)
folding <- sample(rep(seq_len(folds), ceiling(n))[seq_len(n)])

# Number of models:
nmod <- 3

# Constant trimming the acceptance rate of the candidates
ref.kappa <- 0.75
UL.kappa <- 0.3
L.kappa <- 0.75
burn <- 1000
samp <- 2000

###########################################################################



# Storing
acc.prob <- array(0, dim=c(length(nlabels), nmod, folds))
acc.max <- array(0, dim=c(length(nlabels), nmod, folds))
mean.entropy <- array(0, dim=c(length(nlabels), nmod, folds))


# # To test
m <- 1
k <- 1

for(m in seq_along(nlabels)){
  for(k in seq_len(folds)){
    
    # Pick out training set
    y.train <- y[folding != k]
    x.train <- x[folding != k, ]
    
    # Draws a stratified random sample to "unlabel" a part of the instances
    pick <- strata(data.frame(y.train, x.train), 'y.train', 
                   round(table(y.train)*(sum(table(y.train))-nlabels[m])/sum(table(y.train))), method = 'srswor')$ID_unit
    y.na <- y.train
    y.na[pick] <- NA
    
    #### Modeling
    
    ### Reference model on all of the data, nothing unlabeled
    ref <- ssl_MH(y = y.train, x = x.train, 
                  burn = burn, samp = samp, 
                  kappa = ref.kappa, tau=10)
    # Checks that acceptance rate is acceptable
    while(ref$a.rate < 0.1 | ref$a.rate > 0.5){
      print(paste('Acceptance rate on fold', k, 
                  'and label ratio', m, 'for the ref model is bad:', 
                  ref$a.rate))
      ref.kappa <- ifelse(ref$a.rate < 0.1 , ref.kappa/2, ref.kappa*2)
      ref <- ssl_MH(y = y.train, x = x.train, 
                    burn = burn, samp = samp, 
                    kappa = ref.kappa, tau = 10)
    }
    
    ### Model on both labeled and unlabeled instances
    UL <- ssl_MH(y = y.na, x = x.train, 
                 burn = burn, samp = samp, 
                 kappa = UL.kappa, tau = 10, init.yu = 'none')
    # Checks that acceptance rate is acceptable
    while(UL$a.rate < 0.1 | UL$a.rate > 0.5){
      print(paste('Acceptance rate on fold', k, 
                  'and label ratio', m, 'for the UL model is bad:', 
                  UL$a.rate))
      UL.kappa <- ifelse(UL$a.rate < 0.1 , UL.kappa/2, UL.kappa*2)
      UL <- ssl_MH(y = y.na, x = x.train, 
                   burn = burn, samp = samp, 
                   kappa = UL.kappa, tau = 10)
    }
    
    ### Model on just labeled instances
    L <- ssl_MH(y = y.train[!(is.na(y.na))], 
                x = x.train[!(is.na(y.na)),], 
                burn = burn, samp = samp, 
                kappa = L.kappa, tau = 10)
    # Checks that acceptance rate is acceptable
    while(L$a.rate < 0.1 | L$a.rate > 0.5){
      print(paste('Acceptance rate on fold', k, 
                  'and label ratio', m, 'for the L model is bad:', 
                  L$a.rate))
      L.kappa <- ifelse(L$a.rate < 0.1 , L.kappa/2, L.kappa*2)
      L <- ssl_MH(y = y.train[!(is.na(y.na))], 
                  x = x.train[!(is.na(y.na)),], 
                  burn = burn, samp = samp, 
                  kappa = L.kappa, tau = 10)
    }
    
    #### Predictions
    
    pred.ref <- predict.with.posterior(post.beta = ref$beta, X = X[folding == k,], 
                                       levels = levels(y), y =y[folding == k])
    acc.prob[m,1,k] <- pred.ref$a.p
    acc.max[m,1,k] <- pred.ref$a.m
    mean.entropy[m,1,k] <- mean(pred.ref$entr)
    
    pred.UL <- predict.with.posterior(UL$beta, X[folding == k,], 
                                      levels = levels(y), y =y[folding == k])
    acc.prob[m,2,k] <- pred.UL$a.p
    acc.max[m,2,k] <- pred.UL$a.m
    mean.entropy[m,2,k] <- mean(pred.UL$entr)
    
    pred.L <- predict.with.posterior(L$beta, X[folding == k,], 
                                     levels = levels(y), y =y[folding == k])
    acc.prob[m,3,k] <- pred.L$a.p
    acc.max[m,3,k] <- pred.L$a.m
    mean.entropy[m,3,k] <- mean(pred.L$entr)
  }
}

# # Just a check
# mcmc.plots(beta.draws = ref$beta, pdf = TRUE, filename = 'ref-trace.pdf')
# mcmc.plots(UL$beta, pdf = TRUE, filename = 'UL-trace.pdf')
# mcmc.plots(L$beta, pdf = TRUE, filename = 'L-trace.pdf')
# 

save(acc.prob, acc.max, mean.entropy, file='C:/Users/erogeka/Documents/Analysis/cv-ssl/result-iris.rda')
# load(file='C:/Users/erogeka/Documents/Analysis/cv-ssl/result-iris.rda')


# Evaluation --------------------------------------------------------------

library(ggplot2)

## Accuracy prob
dimnames(acc.prob) <- list(nlabels, c('Reference', 'Labeled-Unlabeled', 'Labeled'), 1:10)
melted.acc.prob <- melt(acc.prob)
colnames(melted.acc.prob) <- c('Labeled', 'Model', 'Fold', 'Accuracy') 
melted.acc.prob$Fold <- as.factor(melted.acc.prob$Fold)

# plot.a.prob <- ggplot(melted.acc.prob, aes(x=Missing, y=Accuracy, color=Model)) +
#   geom_line(stat="summary", fun.y = "mean") + scale_x_reverse() + theme_bw()
# ggsave(plot.a.prob, file='C:/Users/erogeka/Documents/Analysis/cv-ssl/aprob.pdf',
#        width=3, height=2, scale =2)
# Lines with error bars
pd <- position_dodge(1) # move them .05 to the left and right
acc.p.sum <- summarySE(melted.acc.prob, measurevar="Accuracy", groupvars=c("Labeled","Model"))
plot.a.prob <- ggplot(acc.p.sum, aes(x=Labeled, y=Accuracy, colour=Model, group=Model)) + 
  geom_errorbar(aes(ymin=Accuracy-ci, ymax=Accuracy+ci), width=1,  position=pd) +
  geom_line(position=pd) + geom_point(position=pd) + theme_bw() 
ggsave(plot.a.prob, file='C:/Users/erogeka/Documents/Analysis/cv-ssl/aprob-iris.pdf',
       width=3, height=2, scale =2)


## Accuracy max
dimnames(acc.max) <- list(nlabels, c('Reference', 'Labeled-Unlabeled', 'Labeled'), 1:10)
melted.acc.max <- melt(acc.max)
colnames(melted.acc.max) <- c('Labeled', 'Model', 'Fold', 'Accuracy') 
melted.acc.max$Fold <- as.factor(melted.acc.max$Fold)

# plot.a.max <- ggplot(melted.acc.max, aes(x=Missing, y=Accuracy, color = Model)) +
#   geom_line(stat="summary", fun.y = "mean") +  scale_x_reverse() + theme_bw()
# ggsave(plot.a.max, file='C:/Users/erogeka/Documents/Analysis/cv-ssl/amax.pdf',
#        width=3, height=2, scale =2)
# Lines with error bars
acc.m.sum <- summarySE(melted.acc.max, measurevar="Accuracy", groupvars=c("Labeled","Model"))
plot.a.max <- ggplot(acc.m.sum, aes(x=Labeled, y=Accuracy, colour=Model, group=Model)) + 
  geom_errorbar(aes(ymin=Accuracy-ci, ymax=Accuracy+ci), width=1,  position=pd) +
  geom_line(position=pd) + geom_point(position=pd) + theme_bw()
ggsave(plot.a.max, file='C:/Users/erogeka/Documents/Analysis/cv-ssl/amax-iris.pdf',
       width=3, height=2, scale =2)


## Entropy
dimnames(mean.entropy) <- list(nlabels, c('Reference', 'Labeled-Unlabeled', 'Labeled'), 1:10)
melted.entropy <- melt(mean.entropy)
colnames(melted.entropy) <- c('Labeled', 'Model', 'Fold', 'Entropy') 
melted.entropy$Fold <- as.factor(melted.entropy$Fold)


# plot.entropy <- ggplot(melted.entropy, aes(x=Missing, y=Entropy, color = Model)) +
#   geom_line(stat="summary", fun.y = "mean") +  scale_x_reverse() + theme_bw() + geom_point()
# ggsave(plot.entropy, file='C:/Users/erogeka/Documents/Analysis/cv-ssl/entropy.pdf',
#        width=3, height=2, scale =2)
# Lines with error bars
entr.sum <- summarySE(melted.entropy, measurevar="Entropy", groupvars=c("Labeled","Model"))
plot.entropy <- ggplot(entr.sum, aes(x=Labeled, y=Entropy, colour=Model, group=Model)) + 
  geom_errorbar(aes(ymin=Entropy-ci, ymax=Entropy+ci), width=1,  position=pd) +
  geom_line(position=pd) + geom_point(position=pd) + theme_bw()
ggsave(plot.entropy, file='C:/Users/erogeka/Documents/Analysis/cv-ssl/entropy-iris.pdf',
       width=3, height=2, scale =2)


# Combined acc and entr for the report
plot.acc <- plot.a.max + theme(legend.position="none")
pdf('C:/Users/erogeka/Documents/Analysis/cv-ssl/combo-iris.pdf', width = 8, height = 3)
do.call("grid.arrange", c(list(plot.acc, plot.entropy), list(ncol = 2, widths = c(2,3.2))))
dev.off()



# To do: 
# - Forsätt utvärdera ssl-modeleringen.
# - Finns det andra sätt för utvärdering?
# - Är nuvarande sätt rimligt och rätt?

