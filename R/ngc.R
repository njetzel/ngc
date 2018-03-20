#' Estimate graphical Granger causality
#' @param X input array
#' @param d number of time lags to consider
#' @param method estimation method to use; options are "regular", "truncate", or "threshold"
#' @param grpIndex vector of group indices of length p or p*d; if null, no group structure
#' @param typeIerr acceptable type I error rate
#' @param typeIIerr acceptable type II error rate
#' @param weights weights for adaptive lasso
#' @param thresholdConstant constant used to compute threshold
#' @param refit whether to refit a linear regression after initial thresholding
#' @param edgeThreshold absolute value threshold for including edge in graph
#' @return a list including a matrix of estimated coefficients, final lambda value, and time series order estimate
#' @examples
#' p <- 15
#' tp <- 20
#' d_actual <- 3
#' d <- 10
#' nsample <- 25
#' sigma <- 0.3
#' edge <- defn_net(d = d_actual, p = p, n = nsample)
#' X <- simulate_data(nsample, edge[1:d_actual,,], tp, error_sd = sigma)
#' fit1 <- ngc(X, method = 'truncate')
#' fit2 <- ngc(X, d = d, method = 'threshold', refit = TRUE)
#' @export ngc
ngc <-
  function(
    X, #input array dim=(n,p,T) (longitudinal), or (p,T) (time series); last time=Y
    d = NULL, #number of time lags to consider
    method = 'regular', #method to use. options are "regular", "truncate", and "threshold"
    grpIndex = NULL, #vector of group indices of length p or p*d; if null, no group structure
    typeIerr = NULL, #acceptable type I error rate; if provided, error-based lasso is fitted
    typeIIerr = 0.1, #acceptable type II error rate
    weights = NULL, #wmatrix of weights for Alasso. If no weights are provided, use regular lasso.
    thresholdConstant = NULL, #constant used for calculating threshold value
    refit = FALSE, #whether to refit a linear regression after initial thresholding
    edgeThreshold = 1e-6 #absolute value threshold for including edge in graph
  ){
    ####### START OF FN #######
    if (method != 'regular' & method != 'truncate' & method != 'threshold')
    {
      stop('Invalid estimation method specified')
    }

    if (!is.array(X) | length(dim(X)) <2 | length(dim(X)) > 3)
    {
      stop('Invalid X')
    }

    if (length(dim(X))==2)
    {
      n <- 1
      p <- dim(X)[1]
      tp <- dim(X)[2]
    }
    else
    {
      n <- dim(X)[1]
      p <- dim(X)[2]
      tp <- dim(X)[3]
    }

    if (is.null(d))
    {
      d <- tp-1
    }

    if (d >= tp)
    {
      stop('Number of time lags to consider cannot exceed number of time points')
    }

    #Set up replicates for the time series case
    #Put X into array format
    #The transformed matrix has (tp-d) replicates over d+1 time points
    if (n == 1)
    {
      if (d >= tp-1)
      {
        stop('Number of time lags to consider must be restricted in order to fit time series')
      }
      cat('Warning: stationarity assumption is required for time series data')
      xMat <- X
      n <- tp-d
      tp <- d+1
      X <- array(0, c(n,p,tp))
      for (i in 1:n)
      {
        X[i,,] <- xMat[,i:(d+i)]
      }
    }
    
    if (!is.null(grpIndex))
    {
      if (length(grpIndex)!=p & length(grpIndex)!=p*d)
      {
        stop('Invalid group specification')
      }
      if (!is.numeric(grpIndex))
      {
        stop('Groups must be specified with consecutive integers')
      }
      ngrp = length(unique(grpIndex))
      if (!all.equal(sort(unique(grp)), 1:ngrp))
      {
        stop("Groups must be specified with consecutive integers")
      }
      # apply groups across time points if group vector length equals p
      if (length(grpIndex) == p)
      {
        grpIndex <- grpIndex + rep(seq(0, (d-1)*ngrp, by = ngrp), each = p)
      }
    }

    if (method == 'regular')
    {
      fit <- grangerLasso(X, d = d, grpIndex = grpIndex, typeIerr = typeIerr,
                          weights = weights)
    }

    else if (method == 'truncate')
    {
      fit <- grangerTlasso(X, d = d, typeIerr = typeIerr,
                           typeIIerr = typeIIerr, weights = weights)
    }

    else #threshold
    {
      fit <- grangerThrLasso(X, d = d, typeIerr = typeIerr,
                             typeIIerr = typeIIerr, weights = weights,
                             thresholdConstant = thresholdConstant,
                             refit = refit)
    }

    fit$estMat <- fit$estMat*(abs(fit$estMat)>edgeThreshold)
    dagMat <- Matrix(0, nrow=p*(d+1), ncol=p*(d+1), sparse = TRUE)
    ringMat <- Matrix(0, nrow=p, ncol=p)
    edgeIx <- which(fit$estMat != 0, arr.ind = T)
    edgeCount <- dim(edgeIx)[1]
    if (is.null(fit$tsOrder))
    {
      tsOrder <- ifelse(edgeCount > 0, max(edgeIx[,3]), 0)
      fit$tsOrder <- ifelse(!is.null(tsOrder), tsOrder, 0)
    }
    if (edgeCount > 0)
    {
      for (i in 1:edgeCount)
      {
        edge <- edgeIx[i,]
        pStart <- edge[2]
        pEnd <- edge[1]
        lag <- edge[3]
        dagMat[((lag-1)*p + pStart),(d*p + pEnd)] <- fit$estMat[pEnd, pStart, lag]
        ringMat[pStart, pEnd] <- ringMat[pStart, pEnd] + abs(fit$estMat[pEnd, pStart, lag])
      } 
    }
    fit$dag <- graph_from_adjacency_matrix(dagMat, mode = 'directed', weighted = TRUE)
    fit$ring <- graph_from_adjacency_matrix(ringMat, mode = 'directed', weighted = TRUE)
    fit$method <- method
    fit$n <- n
    fit$p <- p
    fit$d <- d
    fit$grp <- grpIndex
    class(fit) = "ngc"
    return(fit)
  }

#' Plot DAG of network or ring graph showing Granger causality
#' @param fit object of class ngc
#' @param ngc.ring whether to plot ring graph
plot.ngc <- 
  function(
    fit, #object of class ngc
    ngc.ring = FALSE #whether to plot ring graph
  ){
    if (class(fit) != "ngc")
    {
      stop("Class of argument must be ngc")
    }
    p <- fit$p
    d <- fit$d
    grp <- fit$grp
    if (ngc.ring)
    {
      g <- fit$ring
      edgeThickness = E(g)$weight/mean(E(g)$weight)
      plot(g, layout = layout_in_circle(g), 
           vertex.shape = "none", edge.width = edgeThickness)
    }
    else
    {
      xcoords = rep(1:(d+1), each=p)
      ycoords = rep(p:1, d+1)
      layout_matrix = matrix(c(xcoords, ycoords), ncol=2)
      g <- fit$dag
      groupList <- NULL
      if (!is.null(grp))
      {
        groupList <- lapply(unique(grp),function(x){which(grp==x)})
      }
      par(mar=c(2.5, 2.5, 2.5, 2.5))
      edgeColor = ifelse(E(g)$weight > 0, "green", "red")
      edgeThickness = abs(E(g)$weight/mean(abs(E(g)$weight)))
      plot(g, asp = 0.3, layout=layout_matrix,
           mark.groups = groupList, mark.border = NA,
           vertex.label = rep(1:p, d+1), mark.expand = 20, vertex.shape="none",
           edge.color = edgeColor, edge.width = edgeThickness,
           rescale = FALSE, xlim = c(1, d+1), ylim = c(0, p))
      text(0, -0.5, "Lag")
      for (i in 1:d)
      {
        text(i, -0.5, d-i+1)
      } 
    }
  }

#' Predict covariate values at a given time point
#' @param fit object of class ngc from which to predict
#' @param X input array of size n x p x T
#' @param tp time point at which to predict covariate values 
#' e.g. if tp = 2, output is fitted covariates at time T+2
#' @return n x p matrix of fitted covariates for each replicate
predict.ngc <- 
  function(
    fit, #object of class ngc
    X, #input array dim = (n, p, T) from which to predict
    tp #time point at which to predict covariate value
  ){
    if (class(fit) != "ngc")
    {
      stop("Class of argument must be ngc")
    }
    if (tp < 1)
    {
      stop("Specify positive time point")
    }
    n <- dim(X)[1]
    p <- dim(X)[2]
    d1 <- dim(X)[3]
    if (fit$p != p)
    {
      stop("Incompatible dimensions between fit and input array")
    }
    estMat <- fit$estMat
    d2 <- dim(estMat)[3]
    tsOrder <- fit$tsOrder
    if (tsOrder > d1)
    {
      cat("Warning: X is shorter than estimated time series order")
    }
    i = 1
    while (i <= tp)
    {
      Y <- matrix(rep(fit$intercepts, each = n), nrow = n)
      d1 <- dim(X)[3]
      nlags <- min(d1, tsOrder)
      for (j in 1:nlags)
      {
        Y <- Y + X[,,(d1-j+1)]%*%t(estMat[,,(d2-j+1)])
      }
      X2 <- array(0, c(n, p, d1+1))
      X2[,,1:d1] <- X
      X2[,,(d1+1)] <- Y
      X <- X2
      rm(X2)
      i <- i+1
    }
    return(Y)
  }
