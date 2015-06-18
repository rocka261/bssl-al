
#  ------------------------------------------------------------------------

#  Non-Bayesian supervised learning

#  ------------------------------------------------------------------------

# # For test
# setwd('C:/Users/erogeka/bla/Protoype/Analysis')


# User input 
args <- commandArgs(TRUE)

# Load user defined config file
source('./config.R')


# Add local path to user installed packages
.libPaths( c( .libPaths(), localLibs))

library(nnet)   



# Data --------------------------------------------------------------------


# Topic proportions
load('./results/topics.rda')

labeledLogs <- read.csv2('./data/labeled_logs.csv', row.names = 1)

# Merge topic distributions and info of labels
data <- merge(labeledLogs[, c('general.failure.cause', 
                              'specific.failure.cause')],
              topic.proportions, by=0)
row.names(data) <- data$Row.names; data$Row.names <- NULL

# Set "design/environment" as unknown
data$general.failure.cause[data$general.failure.cause == 'design/environment'] <- NA

# Divide the data on labeled and unlabeled sets
yl <- as.character(data[!is.na(data[,response]),response])
yu <- as.character(data[is.na(data[,response]), response])
Xl <- model.matrix(~. -1, data[!is.na(data[,response]), covariates])
Xu <- model.matrix(~. -1, data[is.na(data[,response]), covariates])
x <- data.frame(rbind(Xl, Xu))
X <- model.matrix(~.,x)
y <- factor(c(yl, yu))


# Modeling ----------------------------------------------------------------

model <- multinom(yl ~ Xl[,])


save(model, file='./results/sl-model.rda')




