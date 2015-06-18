
#  ------------------------------------------------------------------------

#  Predict new data 

#  ------------------------------------------------------------------------


# # For test
# setwd('C:/Users/erogeka/bla/Protoype/Analysis')
# model <- 'bsl-model'

# User input 
args <- commandArgs(TRUE)
modelSelection <- args[1]

# Load user defined config file
source('./config.R')

# Load other functions
source('./helpers.R')
source('./bssl.R')


# Add local path to user installed packages
.libPaths( c( .libPaths(), localLibs))

library(topicmodels)

# Data --------------------------------------------------------------------

# Load the new data
newDTM <- getDTM('data/newData/')

# import list of words after tfidf reduction and lda model
load('./results/lda.rda')
load('./results/tfidfFeatures.rda')

# Select the features left after tfidf reduction
newDTM.tfidf <- rep(0, length(features))
names(newDTM.tfidf) <- features
newDTM.tfidf[features[features %in% newDTM@Dimnames[[2]]]] <- 
  newDTM[,features[features %in% newDTM@Dimnames[[2]]]]
newDTM.tfidf <- Matrix(newDTM.tfidf, nrow = 1, ncol = length(features))
newDTM.tfidf@Dimnames <- list(newDTM@Dimnames[[1]], features)

# Precidt the topic proportions 
lda.pred <- LDA(newDTM.tfidf, model = lda.model, 
                control = list(estimate.beta = FALSE, 
                               burnin = 1000, thin = 100, iter = 1000, 
                               best = TRUE))

x <- data.frame(lda.pred@gamma)
names(x) <- paste('t', 1:ncol(x), sep='')
X <- model.matrix(~., data=x)


# Prediction --------------------------------------------------------------

# Bayesian supervised learning
if(modelSelection == 'bsl'){
  load('./results/bsl-model.rda')
  pred <- predict.with.posterior(model$beta, X=X, levels = model$levels)
  pred$prob  
}

# Bayesian semi-supervised learning
if(modelSelection == 'bssl'){
  load('./results/bssl-model.rda')
  pred <- predict.with.posterior(model$beta, X=X, levels = model$levels)
  pred$prob  
}

# Supervised learning (non-Bayesian)
if(modelSelection == 'sl'){
  load('./results/sl-model.rda')
  pred <- calc_prob(matrix(coef(model)), X)
  colnames(pred) <- model$lev
  pred
}

