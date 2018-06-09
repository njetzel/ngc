#' Use sparse additive modeling to estimate graphical Granger causality
#' @param X input array
#' @param d number of time lags to consider
#' @param deg degree of spline basis; cannot currently be set in ngc function
#' @return a list including a matrix of estimated coefficients
grangerSam <- 
  function(
    X,
    d,
    deg = 3
  )
  {
    n <- dim(X)[1]
    p <- dim(X)[2]
    tp <- dim(X)[3]
    estMat <- array( 0, c(p, p, d) )
    
    XX <- array2mat(X[,,(tp-d):(tp-1)])
    YY <- X[,,tp]
    
    lambdas <- rep(NA, p)
    
    for (i in 1:p)
    {
      y <- YY[,i]
      fit <- samQL(XX, y, p = deg)
      lambda <- fit$lambda
      #choose best lambda with cross-validation
      nfolds <- max(min(10, ceiling(n/5)), 3)
      samp <- sample(n)%%nfolds + 1
      #get lambda*nfolds matrix of validation errors from nfolds-fold cross-validation
      valErr <- sapply(1:nfolds, function(k)
      {
        valIds <- which(samp == k)
        xtr <- XX[-valIds,]
        ytr <- y[-valIds]
        xtest <- XX[valIds,]
        ytest <- y[valIds]
        trFit <- samQL(xtr, ytr, p = deg, lambda = lambda)
        #Fix ridiculous error in samQL output
        trFit$knots = trFit$nkots
        yhatTr <- tryCatch({
          predict.samQL(trFit, xtest)$values
        },
        error = function(cond){
          save(trFit, file="trFit.rdata")
          save(xtest, file="xtest.rdata")
          stop('oops')
        })
        
        return(apply((ytest - yhatTr)^2, 2, mean))
      })
      
      cvm <- apply(valErr, 1, mean)
      cvsd <- apply(valErr, 1, sd)
      ix_best_lambda <- order(cvm)[1]
      ix_lambda_1se <- min(which(cvm < cvm[ix_best_lambda] + cvsd[ix_best_lambda]))
      lambdas[i] <- lambda[ix_lambda_1se]
      #estimate edge matrix
      w <- fit$w[,ix_lambda_1se]
      estMat[i,,] <- sqrt(sapply(seq(1, length(w), deg), function(j)
      {
        sum(w[j:(j+deg-1)]^2)
      }))
    }

    return(list(estMat = estMat, lambdas = lambdas))
  }

#' Use sparse additive modeling to estimate graphical Granger causality
#' @param X input array
#' @param d number of time lags to consider
#' @param nbasis max number of basis functions; cannot currently be set in ngc function
#' @return a list including a matrix of estimated coefficients
grangerHierBasis <- 
  function(
    X,
    d,
    nbasis = 10
  )
  {
    n <- dim(X)[1]
    p <- dim(X)[2]
    tp <- dim(X)[3]
    estMat <- array( 0, c(p, p, d) )
    
    XX <- array2mat(X[,,(tp-d):(tp-1)])
    YY <- X[,,tp]
    
    lambdas <- rep(NA, p)
    
    for (i in 1:p)
    {
      y <- YY[,i]
      fit <- AdditiveHierBasis(XX, y, nbasis = nbasis)
      max.lambda <- max(fit$lam)
      #choose best lambda with cross-validation
      nfolds <- max(min(10, ceiling(n/5)), 3)
      samp <- sample(n)%%nfolds + 1
      #get lambda*nfolds matrix of validation errors from nfolds-fold cross-validation
      valErr <- sapply(1:nfolds, function(k)
      {
        valIds <- which(samp == k)
        xtr <- XX[-valIds,]
        ytr <- y[-valIds]
        xtest <- XX[valIds,]
        ytest <- y[valIds]
        #set max lambda in order to get the same sequence of lambda each time
        trFit <- AdditiveHierBasis(xtr, ytr, nbasis = nbasis, max.lambda = max.lambda)
        yhatTr <- predict(trFit, new.x = xtest)
        
        return(apply((ytest - yhatTr)^2, 2, mean))
      })
      
      cvm <- apply(valErr, 1, mean)
      cvsd <- apply(valErr, 1, sd)
      ix_best_lambda <- order(cvm)[1]
      ix_lambda_1se <- min(which(cvm < cvm[ix_best_lambda] + cvsd[ix_best_lambda]))
      lambdas[i] <- fit$lam[ix_lambda_1se]
      #estimate edge matrix
      beta <- fit$beta[,ix_lambda_1se]
      estMat[i,,] <- sqrt(sapply(seq(1, length(beta), nbasis), function(j)
      {
        sum(beta[j:(j+nbasis-1)]^2)
      }))
    }
    return(list(estMat = estMat, lambdas = lambdas))
  }

