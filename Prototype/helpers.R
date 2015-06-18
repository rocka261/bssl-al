
# Helpers -----------------------------------------------------------------



entr <- function(p){
  # p is a probability distribution
  return(-sum(p*log(p)))
}


# Jensen-Shannon divergence - Symmetric and smoothed KL-divergence. 
JSD  <- function(P, w = 1/2){
  # P is a matrix of the probability distributions, rows sum to 1
  jsd.out <- matrix(0, nrow=nrow(P), ncol=nrow(P))
  for(i in 1:nrow(P)){
    for(j in 1:nrow(P)){
      jsd.out[i,j] <- entr(1/2*(P[i,] + P[j,])) - sum(1/2*(entr(P[i,]) + entr(P[j,])))
    }
  }
  return(jsd.out)
}


getDTM <- function(directory){
  # To load dtm with feature names and log ids
  require(Matrix)
  # Load matrix in MatrixMarket format. 
  dtm <- readMM(paste(directory, 'dtm.mtx', sep=''))
  # Get id and features
  rnames <- scan(paste(directory, 'fileID.txt.txt', sep=''), 
                 what = 'character', nmax = dim(dtm)[1])
  row.names(dtm) <- sapply(strsplit(rnames, '_'), 
                           function(x) paste(x[1:2], collapse = '_'))
  colnames(dtm) <- scan(paste(directory, 'featureNames.txt', sep=''), 
                        what = 'character', nmax = dim(dtm)[2])
  
  return(dtm)
}