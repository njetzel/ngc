#' Estimate graphical Granger causality
#' @param X input array
#' @param d number of time lags to consider
#' @param fitMethod "lasso" for Granger lasso fit; "sam" or "hierbasis" for sparse additive modeling
#' @param estMethod estimation method to use with the lasso; options are "regular", "truncate", or "threshold"
#' @param group vector of group indices of length p; if null, no group structure
#' @param groupByTime whether to group covariates across time points
#' @param typeIerr acceptable type I error rate
#' @param typeIIerr acceptable type II error rate
#' @param weights weights for adaptive lasso
#' @param thresholdConstant constant used to compute threshold
#' @param refit whether to refit a linear regression after initial thresholding
#' @param edgeThreshold absolute value threshold for including edge in graph
#' @param covNames covariate names
#' @return a list including a matrix of estimated coefficients, final lambda value, and time series order estimate
#' @examples
#' p <- 9
#' len <- 20
#' d_actual <- 3
#' d <- 6
#' n <- 30
#' sigma <- 0.3
#' edge <- defn_net(d = d_actual, p = p, n = n)
#' X <- simulate_data(n, edge, len, error_sd = sigma)
#' fit1 <- ngc(X, d = d, typeIerr = 0.05)
#' fit2 <- ngc(X, d = d, estMethod = "threshold", refit = TRUE)
#' @export ngc
ngc <-
  function(
    X, #input array dim=(n,p,T) (longitudinal), or (p,T) (time series); last time=Y
    d = NULL, #number of time lags to consider
    fitMethod = "lasso", #"lasso" for Granger lasso fit; "sam" or "hierbasis" for sparse additive modeling
    estMethod = "regular", #method to use. options are "regular", "truncate", and "threshold"
    group = NULL, #vector of group indices of length p; if null, no group structure
    groupByTime = FALSE, #whether to group covariates across time points
    typeIerr = NULL, #acceptable type I error rate; if provided, error-based lasso is fitted
    typeIIerr = 0.1, #acceptable type II error rate
    weights = NULL, #wmatrix of weights for Alasso. If no weights are provided, use regular lasso.
    thresholdConstant = NULL, #constant used for calculating threshold value
    refit = FALSE, #whether to refit a linear regression after initial thresholding
    covNames = NULL #covariate names
  ){
    ####### START OF FN #######
    if (fitMethod != "lasso" & fitMethod != "sam" & fitMethod != "hierbasis")
    {
      stop("Invalid fit method specified")
    }
    
    if (estMethod != "regular" & estMethod != "truncate" & estMethod != "threshold")
    {
      stop("Invalid estimation method specified")
    }
    
    if (fitMethod != "lasso" & estMethod != "regular")
    {
      stop("This estimation method is not supported for additive modeling")
    }

    if (!is.array(X) | length(dim(X)) <2 | length(dim(X)) > 3)
    {
      stop("Invalid X")
    }

    if (length(dim(X))==2)
    {
      n <- 1
      p <- dim(X)[1]
      len <- dim(X)[2]
    }
    else
    {
      n <- dim(X)[1]
      p <- dim(X)[2]
      len <- dim(X)[3]
    }

    if (is.null(d))
    {
      d <- len-1
    }

    if (d >= len)
    {
      stop("Number of time lags to consider cannot exceed number of time points")
    }

    #Set up replicates for the time series case
    #Put X into array format
    #The transformed matrix has (len-d) replicates over d+1 time points
    if (n == 1)
    {
      if (d >= len-1)
      {
        stop("Number of time lags to consider must be restricted in order to fit time series")
      }
      cat("Warning: stationarity assumption is required for time series data")
      xMat <- X
      n <- len-d
      len <- d+1
      X <- array(0, c(n,p,len))
      for (i in 1:n)
      {
        X[i,,] <- xMat[,i:(d+i)]
      }
    }
    
    if (!is.null(covNames))
    {
      if (length(covNames) != p)
      {
        stop("Number of covariate names must match number of covariates")
      }
    }
    
    if (!is.null(group))
    {
      if (length(group)!=p)
      {
        stop("Invalid group specification")
      }
      if (!is.numeric(group))
      {
        stop("Groups must be specified with consecutive integers")
      }
      if (!all.equal(order(group), 1:p))
      {
        stop("Groups must be specified with consecutive integers")
      }
      # apply groups across time points
      if (groupByTime)
      {
        group <- rep(group, d)
      }
      else
      {
        ngrp = length(unique(group))
        group <- group + rep(seq(0, (d-1)*ngrp, by = ngrp), each = p)
      }
    }

    if (fitMethod == "lasso")
    {
      if (estMethod == "regular")
      {
        fit <- grangerLasso(X, d = d, group = group, typeIerr = typeIerr,
                            weights = weights) 
      }
      else if (estMethod == "truncate")
      {
        fit <- grangerTlasso(X, d = d, group = group, typeIerr = typeIerr,
                             typeIIerr = typeIIerr, weights = weights)
      }
      
      else #threshold
      {
        fit <- grangerThrLasso(X, d = d, group = group, typeIerr = typeIerr,
                               typeIIerr = typeIIerr, weights = weights,
                               thresholdConstant = thresholdConstant,
                               refit = refit)
      }
    }
    else if (fitMethod == "sam")
    {
      fit <- grangerSam(X, d = d)
    }
    else if (fitMethod == "hierbasis")
    {
      fit <- grangerHierBasis(X, d = d)
    }

    dagMat <- Matrix(0, nrow=p*(d+1), ncol=p*(d+1), sparse = TRUE)
    ringMat <- Matrix(0, nrow=p, ncol=p)
    edgeIx <- which(fit$estMat != 0, arr.ind = T)
    edgeCount <- dim(edgeIx)[1]
    if (is.null(fit$tsOrder))
    {
      tsOrder <- ifelse(edgeCount > 0, max(d-edgeIx[,3]+1), 0)
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
        ringMat[pStart, pEnd] <- ringMat[pStart, pEnd] + fit$estMat[pEnd, pStart, lag]^2
      } 
    }
    fit$dag <- graph_from_adjacency_matrix(dagMat, mode = "directed", weighted = TRUE)
    fit$ring <- graph_from_adjacency_matrix(sqrt(ringMat), mode = "directed", weighted = TRUE)
    fit$fitMethod <- fitMethod
    fit$estMethod <- estMethod
    fit$n <- n
    fit$p <- p
    fit$len <- len
    fit$d <- d
    fit$group <- group
    fit$X <- X
    fit$covNames <- covNames
    class(fit) <- "ngc"
    return(fit)
  }

#' Plot DAG of network or ring graph showing Granger causality
#' @param fit object of class ngc
#' @param ngc.type type of plot to create
plot.ngc <- 
  function(
    fit, #object of class ngc
    ngc.type = "dag", #"dag" or "granger"
    sparsify = FALSE #whether to remove covariates with no edges in the high-dimensional case
  ){
    if (class(fit) != "ngc")
    {
      stop("Class of argument must be ngc")
    }
    p <- fit$p
    d <- fit$d
    if (ngc.type == "dag" & p > 20)
    {
      cat("Warning: plot.ngc is not designed for plotting networks of more than 20 covariates.")
      if (sparsify)
      {
        cat("plot.ngc will remove covariates with no edges in order to plot the network.")
      }
    }
    covNames <- fit$covNames
    group <- fit$group
    if (ngc.type == "granger")
    {
      plot_ring(fit$ring, p, d)
    }
    else
    {
      plot_network(fit$dag, p, d, group, fit$fitMethod == "lasso", sparsify)
    }
    if (!is.null(covNames))
    {
      legend(d+1.5, p, paste(1:p, covNames, sep = " - "), cex = labelCex, ncol = p%/%10 + 1, title = "Legend")
    }
  }

plot_ring <- function(g, p, d)
{
  if (is.null(E(g)$weight))
  {
    edgeThickness = 0
  }
  else
  {
    edgeThickness = E(g)$weight^2/mean(E(g)$weight^2)
  }
  #control maximum and minimum thickness
  edgeThickness <- ifelse(edgeThickness > 0.2, edgeThickness, 0.2)
  edgeThickness <- ifelse(edgeThickness < 5, edgeThickness, 5)
  labelCex <- max(min(10/p, 1), 0.3)
  arrowSize <- 0.5*labelCex
  plot(g, layout = layout_in_circle(g), vertex.label.cex = labelCex,
       edge.arrow.size = arrowSize, vertex.shape = "none", edge.width = edgeThickness)
}

plot_network <- function(g, p, d, group =  NULL, signed = TRUE, sparsify = FALSE, title = NULL)
{
  edgeTails <- tail_of(g, E(g))
  edgeHeads <- head_of(g, E(g))
  if (sparsify)
  {
    nodes <- as.numeric(unique(c(edgeHeads%%p, edgeTails%%p)))
    nodes[nodes==0] <- p
    nodes <- sort(nodes)
    toRemove <- (1:p)[-nodes] + rep(seq(0, p*d, p), each = p - length(nodes))
    g <- delete_vertices(g, toRemove)
    p <- length(nodes)
    vertexLabels <- rep(nodes, d+1)
  }
  else
  {
    vertexLabels <- rep(1:p, d+1)
  }
  xcoords = rep(1:(d+1), each=p)
  ycoords = rep(p:1, d+1) 
  layout_matrix = matrix(c(xcoords, ycoords), ncol=2)
  groupList <- NULL
  if (!is.null(group))
  {
    groupList <- lapply(unique(group),function(x){which(group==x)})
  }
  par(mar=c(2.5, 2.5, 2.5, 2.5))
  if (signed)
  {
    edgeColor = ifelse(E(g)$weight > 0, "blue", "red")
  }
  else
  {
    edgeColor = NULL
  }
  if (is.null(E(g)$weight))
  {
    edgeThickness = 0
  }
  else
  {
    edgeThickness <- E(g)$weight^2/mean(E(g)$weight^2)
  }
  #control maximum and minimum thickness
  edgeThickness <- ifelse(edgeThickness > 0.2, edgeThickness, 0.2)
  edgeThickness <- ifelse(edgeThickness < 5, edgeThickness, 5)
  labelCex <- max(min(10/p, 1), 0.3)
  arrowSize <- 0.5*labelCex
  #curve edges that are more than 1 lag
  edgeCurvature <- (edgeTails <= p*(d-1))*0.25
  edgeCurvature <- edgeCurvature*(-1)^((head_of(g, E(g)) %% p) < (edgeTails %% p))
  aRatio <- ((d+3)/p)/2
  plot(g, asp = aRatio, layout = layout_matrix, main = title,
       mark.groups = groupList, mark.border = NA,
       vertex.label.cex = labelCex,
       vertex.label = vertexLabels, vertex.shape = "none",
       edge.color = edgeColor, edge.width = edgeThickness,
       edge.arrow.size = arrowSize, edge.curved = edgeCurvature,
       rescale = FALSE, xlim = c(0, d+2), ylim = c(0, p))
  text(0, -0.25*labelCex, "Lag", cex = labelCex)
  lagStep <- ifelse(d < 10, 1, 5)
  for (i in seq(lagStep, d, lagStep))
  {
    text(i, -0.25*labelCex, d-i+1, cex = labelCex)
  }
}

#' Predict covariate values at a given time point
#' @param fit object of class ngc from which to predict
#' @param X input array of size n x p x T
#' @param len time point at which to predict covariate values 
#' e.g. if len = 2, oulenut is fitted covariates at time T+2
#' @return n x p matrix of fitted covariates for each replicate
predict.ngc <- 
  function(
    fit, #object of class ngc
    tp = 0#time point at which to predict covariate value
  ){
    if (class(fit) != "ngc")
    {
      stop("Class of argument must be ngc")
    }
    if (fit$fitMethod != "lasso")
    {
      stop("Predict method is not currently implemented for this fit method")
    }
    if (tp < 0)
    {
      stop("Specify current or future time point")
    }
    n <- fit$n
    p <- fit$p
    d <- fit$d
    len <- fit$len
    X <- fit$X
    estMat <- fit$estMat
    tsOrder <- fit$tsOrder
    #adjust scaled parameters back to original scale
    intercepts <- fit$intercepts

    i <- 0
    while (i <= tp)
    {
      Y <- matrix(rep(intercepts, each = n), nrow = n, ncol = p)
      d1 <- dim(X)[3]
      if (tsOrder >= 1)
      {
        for (j in 1:tsOrder)
        {
          Y <- Y + X[,,(d1-j)]%*%t(estMat[,,(d-j+1)])
        }
      }
      X2 <- array(0, c(n, p, d1+1))
      X2[,,1:d1] <- X
      scaledY <- scale(Y)*sqrt(n/(n-1))
      X2[,,(d1+1)] <- ifelse(is.nan(scaledY), 0, scaledY)
      X <- X2
      rm(X2)
      i <- i+1
    }
    
    return(Y)
  }
