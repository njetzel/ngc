#' Use adaptively thresholded lasso to estimate graphical Granger causality
#' @param X input array
#' @param d number of time lags to consider
#' @param group group indices
#' @param typeIerr acceptable type I error rate
#' @param typeIIerr acceptable type II error rate
#' @param weights weights for adaptive lasso
#' @param thresholdConstant constant used to compute threshold value
#' @param refit whether to refit a linear regression after initial thresholding
#' @return a list including a matrix of estimated coefficients, final lambda value, and time series order estimate
grangerThrLasso <-
  function(
    X, #input array dim=(n,p,T), last time=Y
    d = NULL, #number of time lags to consider
    group = NULL, #group indices
    typeIerr = NULL, #acceptable type I error rate
    typeIIerr = 0.10, #acceptable type II error rate
    weights = NULL, #weights for granger lasso
    thresholdConstant = NULL, #constant used to compute threshold value. If null, use CV
    refit = FALSE #whether to refit a linear regression after initial thresholding
  )
  {
    ####### START OF FN #######
    n <- dim(X)[1]
    p <- dim(X)[2]
    tp <- dim(X)[3]

    #initial granger lasso fit and refit
    grangerLassoFit <- grangerLasso(X, d = d, group = group, typeIerr = typeIerr, weights = weights)
    lambda <- grangerLassoFit$lambda
    sigma <- grangerLassoFit$sigma
    intercepts <- grangerLassoFit$intercepts
    lambda_not <- sqrt(2*log((tp-1)*p)/n)

    #set thresholding constants
    edgeThreshold <- p^2*typeIIerr/(tp-1)
    #set threshold value
    #if thresholdConstant isn't provided, use CV
    if (!is.null(thresholdConstant))
    {
      thresholdValue <- thresholdConstant*lambda_not*sigma
    }
    else #min(10,n)-fold cross-validation
    {
      nfolds <- min(10,n)
      thresholdValues <- seq(0.2, 0.9, 0.1)*lambda_not*sigma
      cvMseMin <- Inf
      thresholdValue <- NA
      #guarantee each fold has roughly equal size (some may have 1 more observation than others)
      foldid <- rep(1:nfolds, floor(n/nfolds))
      #reshuffle
      foldid <- sample(c(foldid, sample(1:nfolds, n-length(foldid))))
      for (c2 in thresholdValues)
      {
        #threshold coefficients
        thresholdEst <- thresholdGgcFit(grangerLassoFit, edgeThreshold, c2)
        estMat <- thresholdEst$estMat
        mseSum <- 0
        #iterate over folds and perform CV
        #fit linear regression on nonzero coefficients and validate on held-out data
        #MSE for the fit on the held-out fold
        for (i in 1:nfolds)
        {
          mseSum <- mseSum + fitLm(X, estMat, d, p, tp, foldid != i)$mse
        }
        #compare mean MSE over each held-out fold with the current minimum
        #if we have found a new minimum, set c2 as the threshold value
        if (mseSum/nfolds < cvMseMin)
        {
          cvMseMin <- mseSum/nfolds
          thresholdValue <- c2
        }
      }
    }

    thresholdEst <- thresholdGgcFit(grangerLassoFit, edgeThreshold, thresholdValue)
    estMat <- thresholdEst$estMat
    tsOrder <- thresholdEst$tsOrder
    #Refit using standard linear regression on nonzero entries
    if (refit==TRUE)
    {
      lmFit <- fitLm(X, estMat, d, p, tp)
      estMat <- lmFit$estMat
      intercepts <- lmFit$intercepts
    }

    return(list(estMat = estMat, lambda = lambda, tsOrder = tsOrder, intercepts = intercepts))
  }

  thresholdGgcFit <- function(ggcFit, edgeThreshold, thresholdValue)
  {
    grangerLassoEst = adj_left2right(ggcFit$estMat)
    d <- dim(grangerLassoEst)[1]
    p <- dim(grangerLassoEst)[2]
    thresholdEst = array(0, c(d, p, p))
    tsOrder <- 0
    #Threshold the adjacency matrices
    for (i in 1:d)
    {
      adjacencyMat <- grangerLassoEst[i,,]
      nonzeroEdgeCount <- sum(adjacencyMat!=0)
      if (i > 1 & nonzeroEdgeCount < edgeThreshold)
      {
        next
      }
      else
      {
        thresholdEst[i,,] <- adjacencyMat*1*(abs(adjacencyMat)>thresholdValue)
        if (sum(thresholdEst[i,,]!=0)>=edgeThreshold)
        {
          tsOrder <- i
        }
      }
    }
    return(list(estMat = adj_right2left(thresholdEst), tsOrder = tsOrder))
  }

  #fits a linear regression on nonzero coefficients
  #if cvIx is not null, use to construct train and test set
  fitLm <- function(X, estMat, d, p, tp, cvIx = NULL)
  {
    # create n x (p*d) matrix of predictors
    predictorMatrix <- array2mat(X[,,(tp-d):(tp-1)])
    mse <- rep(NA, p)
    intercepts <- rep(0, p)
    for (i in 1:p)
    {
      #get response for ith covariate
      response <- X[,i,tp]
      #set nonzero predictors for ith response variable
      nonZeroIx <- which(estMat[i,,]!=0)
      if (length(nonZeroIx)>0)
      {
        predictors <- predictorMatrix[,nonZeroIx, drop = FALSE]
        #make sure predictors is in matrix format for compatibility
        if (is.null(cvIx))
        {
          lmfit <- lm(response~predictors)
          coefs <- coef(lmfit)
          intercepts <- coefs[1]
          coefs <- coefs[-1]
          coefs[is.na(coefs)] <- 0
        }
        else #split into train and test set
        {
          trainPredictors <- predictors[cvIx,, drop=FALSE]
          testPredictors <- predictors[!cvIx,, drop=FALSE]
          trainResponse <- response[cvIx]
          testResponse <- response[!cvIx]
          lmfit <- lm(trainResponse~trainPredictors)
          coefs <- coef(lmfit)[-1]
          coefs[is.na(coefs)] <- 0
          # compute mean squared difference between response values and fitted values in the test data set
          mse[i] <- sum((testResponse-(coef(lmfit)[1] + as.vector(testPredictors%*%coefs)))^2)/length(testResponse)
        }
        betas <- rep(0,d*p)
        betas[nonZeroIx] <- coefs
        estMat[i,,] = betas
      }
      else if (!is.null(cvIx))
      #If no nonzero predictors, we can't fit a regression
      #Instead compute sum of squared deviations for test data
      {
        testResponse <- response[!cvIx]
        mse[i] <- sum((testResponse-mean(testResponse))^2)/length(testResponse)
      }
    }
    return(list(estMat = estMat, mse = mean(mse), intercepts = intercepts))
  }
