#' Use lasso to estimate graphical Granger causality
#' @param X input array
#' @param d number of time lags to consider
#' @param group vector of group indices of length p*d; if null, no group structure
#' @param penCoefMethod penalty method for lasso; options are "errBased" and "direct"
#' @param typeIerr acceptable type I error rate
#' @param lambda penalty coefficient for lasso
#' @param weights weights for adaptive lasso
#' @param Alasso.power power for adaptive lasso
#' @param eps threshold to set small values to 0
#' @return a list including a matrix of estimated coefficients and final lambda value
grangerLasso <-
  function(
    X, 				  #input array dim=(n,p,T) (longitudinal), or (p,T) (time series); last time=Y
    d = NULL, 	  #number of time lags to consider
    group = NULL, #group indices
    penCoefMethod = 'errBased', #choose between errBased and direct
    typeIerr = 0.10, 		  #sig level for lambda (...Method=errBased)
    lambda = NULL,		  #value of lambda (...Method=direct)
    weights = NULL,			  #matrix of weights for Alasso
    Alasso.power = 1,		  #power for Adaptive lasso
    eps = 1e-8			  #used for setting small values to zero
  ){
    ####### START OF FN #######
    n <- dim(X)[1]
    p <- dim(X)[2]
    tp <- dim(X)[3]

    useAlasso <- !is.null(weights)
    Onep = matrix(1,p,(p*d))
    estMat <- array( 0, c(p, p, d) )

    ##scale the X matrix
    for (i in 1:(tp-1))
    {
      X[,,i] <- scale( X[,,i] )*sqrt(n/(n-1))
    }

    ##first put all X matrices into one big matrix
    XX <- array2mat(X[,,(tp-d):(tp-1)])
    YY <- X[,,tp]

    if (!useAlasso)
    {
      temp = pldag.set(XX, YY, group = group, sigLevel=typeIerr, wantScale=TRUE)
    }
    else
    {
      W = abs( array2mat(weights[,,1:d]) )
      W = ( W + eps*Onep ) ^ (-Alasso.power)
      if (sum(W < 1) > 0)
      {
        W[(W < 1)] = 1
      }

      temp = pldag.set(XX, YY, group = group, sigLevel=typeIerr,
                      useWghts=TRUE, wghts=W, wantScale=TRUE)
    }
    AA <- as.matrix(temp$AA)
    lambda <- temp$lambda
    sigma <- temp$sigma
    intercepts <- temp$intercepts
    ##Put the matrix output of pldag.set into an array to make it
    ##compatible with other parts of the code
    for(i in 1:p)
    {
      estMat[i,,] = AA[i,]
    }
    rm(temp)

    return(list(estMat = estMat, lambda = lambda, sigma = sigma, intercepts = intercepts))
  }

