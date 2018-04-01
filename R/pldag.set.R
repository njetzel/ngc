#' Fit the lasso with glmnet
#' @param X1 predictor matrix
#' @param X2 response matrix
#' @param group vector of group indices matching number of columns in X1; if null, no group structure
#' @param sigLevel significance level for penalties
#' @param useWghts whether to use weights for the penalty
#' @param wghts weights for penalty
#' @param excld which variables to exclude
#' @param wantScale whether to scale
#' @param useNaive whether to use naive or covariance method for glmnet
#' @param useTLASSO whether to use the truncating lasso
#' @param d number of time lags to consider
#' @return a list including fitted coefficients and the final lambda value
pldag.set <-
function(
	X1,				##predictor set, nxp1 matrix
	X2, 				##outcome set, nxp2 matrix (X2~X1)
	group = NULL, #group indices
	sigLevel = NULL, 		##siglevel for MB & BdE penalties
	useWghts = FALSE,		##use wghts for the penalty?
					##if TRUE, wghts should be provided
	wghts=NULL, 		##wghts for penalty, e.g. Alasso
	excld=NULL,			##which variabels to exclude (this
					##is equiv to setting wght=Inf)
	wantScale = FALSE,	##whether to scale
	useNaive=TRUE,		##see below!
	useTLASSO=FALSE,		##is truncating lasso used?
	d = NULL			##number of time lags to cosider
){
####### START OF FN #######
	method <- "naive"		#method used in glmnet (see glmnet help)
	if (!useNaive)
  {
		method <- "covariance"
	}

	if (dim(X1)[1] != dim(X2)[1])
  {
		stop("Number of observations does not  match")
	}
	n <- dim(X2)[1]

	p1 <- dim(X1)[2]
	p2 <- dim(X2)[2]
	prdctrIndx = 1:p1

	if (is.null(wghts))
  {
		if (useWghts)
	  {
			cat("WARNING: No weights provided, using weights=1", "\n")
		}
	  if (is.null(group))
	  {
	    wghts <- matrix(1,p2,p1)
	  }
	  else
	  {
	    wghts <- matrix(rep(as.vector(sqrt(table(group))), p2), nrow=p2)
	  }
	}

	##calculate penalty coefficient (lambda)
	nvar <- p2		#no of variables, to use in calculation of lambda
	ncov <- p1		#no of covariates, to use in calculation of lambda

	#in case of TLASSO, one set of values of X is given at each time, but
	#ncov needs to be adjusted based on the truncating pealty
	if (useTLASSO)
  {
		if (is.null(d))
	  {
			stop('Number of effective time lags needed for TLASSO')
		}
	  else
    {
			ncov <- d*p1
		}
	}

	AA <- matrix(0, p2, p1)

	if ( is.null(excld) )
  {
		excld <- matrix(FALSE,p2,p1)
	}

	if ( (dim(excld)[1]!= dim(AA)[1]) || (dim(excld)[2]!= dim(AA)[2]) )
  {
		stop("Wrong dimension for variables to exclude")
	}

	lambda <- NULL
	if(!is.null(sigLevel))
	{
    lambda <- (1/sqrt(n))*qnorm(1-sigLevel/(2*ncov*nvar))
	}

	lambdas <- rep(NA, p2)
	sigmas <- rep(NA, p2)
	sigmas2 <- rep(NA, p2)
	intercepts <- rep(0, p2)
	##main estimation loop
	for (i in 1:p2)
  {
		y <- X2[ , i]
		ww <- wghts[i, ]

		temp <- excld[i, ]
		excldIndx <-prdctrIndx[temp]

		if (length(excldIndx) < p1-1)
	  {
		  #estimate sigma with deviance
		  if (!is.null(lambda))
		  {
		    if (is.null(group))
		    {
		      fit1 <- glmnet(X1, y, lambda = lambda, penalty.factor = ww,
		                     standardize = wantScale, exclude = excldIndx)
		    }
		    else
		    {
		      fit1 <- gglasso(X1, y, group = sort(group), loss="ls", 
		                      lambda = lambda)
		    }
		    betas <- coef(fit1)
		    intercepts[i] <- betas[1]
		    betas <- betas[-1]
		    dev <- deviance(fit1)
		    if (!is.null(dev))
		    {
		      residDf <- max(n-1-sum(betas!=0), 1)
		      sigmas[i] <- sqrt(dev/residDf)
		    }
		    else
		    {
		      sigmas[i] = 1
		    }
		  }
		  else #estimate sigma with cvm
		  {
		    if (is.null(group))
		    {
		      fit1 <- cv.glmnet(X1, y, penalty.factor = ww,
		                        standardize = wantScale, exclude = excldIndx)
		    }
		    else
		    {
		      fit1 <- cv.gglasso(X1, y, group = sort(group), pred.loss = "L1", pf  = ww)
		    }
		    lambdas[i] <- fit1$lambda.1se
		    betas <- coef(fit1, s="lambda.1se")
		    intercepts[i] <- betas[1]
		    betas <- betas[-1]
		    sigmas[i] <- mean(sqrt(fit1$cvm))
		  }
			if (length(betas) > 0)
		  {
				AA[i, ] <- betas
			}
			rm(fit1)
			rm(betas)
		}
	}
	if (is.null(lambda))
	{
	  lambda <- lambdas
	}
  return(list(AA = AA, lambda = lambda, sigma = mean(sigmas, na.rm = TRUE), intercepts = intercepts))
}

