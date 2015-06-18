
#  ------------------------------------------------------------------------

#  Active learning for multinomial logistic regression


#  ------------------------------------------------------------------------

library(nnet)

# User input 
args <- commandArgs(TRUE)
alMethod <- args[1]

# # For test 
# alMethod <- 'exp.model.change.ml'


# Load user defined config file
source('./config.R')


# Add local path to user installed packages
.libPaths( c( .libPaths(), localLibs))


source("./helpers.R")
source("./al_functions.R")
source("./bssl.R")





# Data --------------------------------------------------------------------


# Topic proportions
load('./results/topics.rda')

labeledLogs <- read.csv2( labeledLogsDir, row.names = 1)

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




# Model -------------------------------------------------------------------


# fit model
model <- multinom(yl ~ Xl)


# Predict unlabeled data
yu.probs <- calc_prob(matrix(coef(model)), cbind(1,Xu))

query <- do.call('exp.model.change.ml', args = 
                   list(yu.probs, y[is.na(y)], Xu, 
                        as.vector(t(coef(model))), loss='0-1'))

names(query) <- row.names(Xu)[query]

print(paste("The selected query with", alMethod, "is", names(query)))









