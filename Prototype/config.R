
#  ------------------------------------------------------------------------

#  User defined parameters for build log analysis

#  ------------------------------------------------------------------------



# Set this to the directory of your user installed R-packages
localLibs <- "C:/Users/erogeka/Documents/R/win-library/3.1"
# This is the file contining labeled logs
labeledLogsDir <- 'C:/Users/erogeka/bla/Protoype/Analysis/data/labeled_logs.csv'

# Set resonse for modeling and choose covariates
response <- 'general.failure.cause'
covariates <- paste('t', 1:20, sep='')

# Parameter settings for modeling 
tau <- 10
kappa <- 0.5   # Alter this parameter if the acceptance rate is to high or low
burn <- 2000
samp <- 4000
