## run GDM models on inverts

options(warn=1)

# require(SDMTools)
# require(raster)
##source("Y:\\MOD\\G\\GDM\\_tools\\obs.gdm.plot_logit.R")
##source("Y:\\MOD\\G\\GDM\\_tools\\create_5splineOutput2.R")
##source("Y:\\MOD\\G\\GDM\\_tools\\gdm.spline.plot.R")

##############################################################################################################################
##
## FUNCTIONS TO RUN AN OBSERVATION-PAIRS GDM
## 
## These functions are developed to support the fitting of an observation pairs gdm as described in:
## 	Hoskins et al. (2019) BILBI: supporting global biodiversity assessment through high-resolution macroecological modelling.
##
## The model fitting procedure 


##############################################################################################################################
##
## The inverse logit function for converting 0-1 logit predictions into logit space
##
## Inputs:
## 		x - a value between 0 and 1 representing 
##
## Outputs:
## 		x transformed into its' representation in logit space
##
inv.logit <- function(x){exp(x)/(1+exp(x))}
##
##############################################################################################################################

##############################################################################################################################
##
## Transform predictions from an observation pair GDM into estimates of 
##
## Inputs:
## 		p0 - intercept value 
## 		w  - precalculated weighting value. See eq. 4 and supporting text for details
## 		p  - 0-1 representing a vector of predicted values from an obs-pair GDM 
##
## Outputs:
## 		prw - p reweighted to represent estimated p if downweighting of 0's hadn't occurred
## 		out - 0-1 representing the transformation of p into an estimated pairwise dissimilarity value
##
ObsTrans <- function(p0,w,p){
	prw <- (p*w) / ((1-p) + (p*w))
	p0w <- (p0*w) / ((1-p0) + (p0*w))
	out <- 1 - ((1-prw) / (1-p0w))
	return(list(prw=prw,out=out))
}
##
##############################################################################################################################

##############################################################################################################################
##
## I-spline transformation of intput vector 
##		- following Ramsay (1988) Monoton Regression Splines in Action. Statistical science.
##
## Inputs:
## 		predVal - vector of values to be splined 
## 		q1 - first quantile for transformation
## 		q2 - second quantile for transformation
## 		q3 - third quantile for transformation
##
## Outputs:
## 		vector representing the I-spline transformed values
##
# Calculate the I-Spline value for predVal given quantiles q1, q2, q3
I_spline <- function(predVal, q1, q2, q3 ){
	outVal <- rep(NA,length(predVal))
	outVal[predVal <= q1] <- 0
	outVal[predVal >= q3] <- 1
	outVal[predVal > q1 & predVal <= q2] <- ( ( (predVal[predVal > q1 & predVal <= q2] - q1) * (predVal[predVal > q1 & predVal <= q2] - q1) ) / ( (q2-q1) * (q3-q1) ) )
	outVal[predVal > q2 & predVal < q3] <- ( 1.0 - ( ( (q3 - predVal[predVal > q2 & predVal < q3]) * (q3 - predVal[predVal > q2 & predVal < q3]) ) / ( (q3-q2) * (q3-q1) ) ) )	
	return(outVal)     
}
##
##############################################################################################################################


plotSpline <- function(env, nSplines, pCoeffs, pQuants,xlab=NA,ylab=NA){
	dCalc = 0.0
	for(s in 1:nSplines){
		if(pCoeffs[s] == 0.0 ){
		# Skip calcs for zero coefficient
		} else {
		# This is the first knot
		if( s == 1 ){
			d0 = pCoeffs[s] * sapply(env, I_spline,
			pQuants[1], pQuants[1], pQuants[2] )
			dCalc = dCalc + d0
		} else {
			# This is the last knot
			if( s == nSplines ){
				d2 = pCoeffs[s] * sapply(env, I_spline,
				pQuants[nSplines-1], pQuants[nSplines], pQuants[nSplines] )
				dCalc = dCalc + d2
			} else {
			# This is a middle knot
					d1 = pCoeffs[s] * sapply(env, I_spline,
					pQuants[s-1], pQuants[s], pQuants[s+1] )
					dCalc = dCalc + d1
				}
			}
		}
	}
	plot(dCalc~env,type="l",xlab=xlab,ylab=ylab)
}

prepSpline <- function(env, nSplines, pCoeffs, pQuants, nmn, zero){
  dCalc = 0.0
  for(s in 1:nSplines){
    if(pCoeffs[s] == 0.0 ){
      # Skip calcs for zero coefficient
    } else {
      # This is the first knot
      if( s == 1 ){
        d0 = pCoeffs[s] * sapply(env, I_spline,
                                 pQuants[1], pQuants[1], pQuants[2] )
        dCalc = dCalc + d0
      } else {
        # This is the last knot
        if( s == nSplines ){
          d2 = pCoeffs[s] * sapply(env, I_spline,
                                   pQuants[nSplines-1], pQuants[nSplines], pQuants[nSplines] )
          dCalc = dCalc + d2
        } else {
          # This is a middle knot
          d1 = pCoeffs[s] * sapply(env, I_spline,
                                   pQuants[s-1], pQuants[s], pQuants[s+1] )
          dCalc = dCalc + d1
        }
      }
    }
  }
  return(data.frame(dCalc = dCalc, env = env, nmn = nmn, zero = zero))
}

##############################################################################################################################
##
## Take gdm formatted table and create spline distance table
##
## Inputs:
## 		X          - data.frame or matrix containing environmental values for pairs fo sites 
## 		splines    - number of splines to be used for each environmental covariate
## 		quantiles  - quantile values for each set fo splines 
##
## Outputs:
## 		data.frame of the distances in spline transformed space of all input site pairs
##
##  
splineData <- function(X,splines=NULL,quantiles=NULL){

	## fold X and create site vector
	nc <- ncol(X)
	nc2 <- nc/2
	if(nc %% 2 != 0){stop("X must be a matrix with even columns")}
	X1 <- X[,1:nc2]
	X2 <- X[,(nc2+1):nc]
	nms <- colnames(X1)
	colnames(X2) <- nms
	## site vector
	sv <- c(rep(1,nrow(X1)),rep(2,nrow(X2)))
	XX <- rbind(X1,X2)

	
	## error checking
	if(is.null(splines)){
		message("No splines specified. Using 3 splines (0%, 50%, 100% quantiles)")
		splines <- 	rep(3,ncol(XX))
	}
	if(is.null(quantiles)){
		if(all(splines != 3)){stop("Must specify quantile positions if all(splines) != 3")}
		quantiles <- unlist(lapply(1:ncol(XX),function(x){quantile(XX[,x],c(0,0.5,1))}))
	}
	if(length(quantiles) != sum(splines)){stop("Number of quantiles must equal number of splines")}
	
	## spline data
	csp <- c(0,cumsum(splines))
	out.tab <- c()
	
	for(c in 1:ncol(XX)){
		ns <- splines[c]
		predVal <- XX[,c]
		quan <- quantiles[(csp[c]+1):(csp[c]+ns)]
		for(sp  in 1:ns){
			if(sp == 1){spl <- I_spline(predVal,quan[1],quan[1],quan[2])}
			if(sp == ns){spl <- I_spline(predVal,quan[ns-1],quan[ns],quan[ns])}
			if(sp != 1 & sp != ns){spl <- I_spline(predVal,quan[sp-1],quan[sp],quan[sp+1])}
			out.tab <- cbind(out.tab,spl)
			##if(anyNA(spl)){print(c);print(sp);break} ## error catch, unhash to use
		}
	}
	NMS <- rep(nms,splines)
	SPNMS <- paste("spl",unlist(lapply(splines,function(x){1:x})),sep="")
	NMS <- paste(NMS,SPNMS,sep="_")
	colnames(out.tab) <- NMS
	XX1 <- out.tab[sv == 1,]
	XX2 <- out.tab[sv == 2,]
	return(abs(XX1 - XX2))
}

splineDataNew <- function(X,splines=NULL,quantiles=NULL){
  
  ## fold X and create site vector
  nc <- ncol(X)
  nc2 <- nc/2
  if(nc %% 2 != 0){stop("X must be a matrix with even columns")}
  X1 <- X[,1:nc2]
  X2 <- X[,(nc2+1):nc]
  nms <- colnames(X1)
  colnames(X2) <- nms
  ## site vector
  sv <- c(rep(1,nrow(X1)),rep(2,nrow(X2)))
  XX <- rbind(X1,X2)
  
  
  ## error checking
  if(is.null(splines)){
    message("No splines specified. Using 3 splines (0%, 50%, 100% quantiles)")
    splines <- 	rep(3,ncol(XX))
  }
  if(is.null(quantiles)){
    if(all(splines != 3)){stop("Must specify quantile positions if all(splines) != 3")}
    quantiles <- unlist(lapply(1:ncol(XX),function(x){quantile(XX[,x],c(0,0.5,1))}))
  }
  if(length(quantiles) != sum(splines)){stop("Number of quantiles must equal number of splines")}
  
  ## spline data
  csp <- c(0,cumsum(splines))
  out.tab <- c()
  
  for(c in 1:ncol(XX)){
    ns <- splines[c]
    predVal <- XX[,c]
    quan <- quantiles[(csp[c]+1):(csp[c]+ns)]
    for(sp  in 1:ns){
      if(sp == 1){spl <- I_spline(predVal,quan[1],quan[1],quan[2])}
      if(sp == ns){spl <- I_spline(predVal,quan[ns-1],quan[ns],quan[ns])}
      if(sp != 1 & sp != ns){spl <- I_spline(predVal,quan[sp-1],quan[sp],quan[sp+1])}
      out.tab <- cbind(out.tab,spl)
      ##if(anyNA(spl)){print(c);print(sp);break} ## error catch, unhash to use
    }
  }
  NMS <- rep(nms,splines)
  SPNMS <- paste("spl",unlist(lapply(splines,function(x){1:x})),sep="")
  NMS <- paste(NMS,SPNMS,sep="_")
  colnames(out.tab) <- NMS
  XX1 <- out.tab[sv == 1,]
  XX2 <- out.tab[sv == 2,]
  return(list(abs(XX1 - XX2), quantiles))
}

##
##############################################################################################################################

##############################################################################################################################
##
## Wrapper function for fitting GDM using base R glm functions. 
##
## Inputs:
## 		formula - intercept value 
## 		data    - precalculated weighting value. See eq. 4 and supporting text for details
##
## Outputs:
## 		standard glm model object fitted to splined data
##
fitGDM <- function(formula=NULL,data=gdm.tab,family=binomial(link=negexp()),method='nnls.fit',...){

	if(is.null(formula)){
		f1 <- paste(colnames(data)[-1],collapse="+")
		formula <- as.formula(paste(colnames(data)[1],"~",f1,sep=""))
		}

	glm(formula,family=family,data=data,control=list(maxit=500),method=method)
	
}
##
##############################################################################################################################

##############################################################################################################################
##
## Negative exponential link funciton for fitting GDM models
##
negexp <- function()
{
	## link
    linkfun <- function(mu) -log(1-mu)
	## inverse link
    linkinv <- function(eta) 1-exp(-eta)
	## derivative of inverse link wrt eta
    mu.eta <- function(eta) exp(-eta)
    valideta <- function(eta) all(is.finite(eta))
    link <- paste0("negexp")
    structure(list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta, name = link),
              class = "link-glm")
}
##
##############################################################################################################################

##############################################################################################################################
##
## glm.fit function modified to allow forced non-negative acoefficient estimates
## 		see ?glm.fit for detail on glm.fit and ?nnnpls for details on nnnpls
##
nnls.fit <- function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
	mustart = NULL, offset = rep(0, nobs), family = gaussian(),
	control = list(), intercept = TRUE, singular.ok = TRUE){
    control <- do.call("glm.control", control)
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y)) 
        rownames(y)
    else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights)) 
        weights <- rep.int(1, nobs)
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("'family' argument seems not to be a valid family object", 
            call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null) if (is.null(x)) 
        if.null
    else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta)) 
            stop("invalid linear predictor values in empty model", 
                call. = FALSE)
        mu <- linkinv(eta)
        if (!validmu(mu)) 
            stop("invalid fitted means in empty model", call. = FALSE)
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep_len(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric()
        iter <- 0L
    }
    else {
        coefold <- NULL
        eta <- if (!is.null(etastart)) 
            etastart
        else if (!is.null(start)) 
            if (length(start) != nvars) 
                stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", 
                  nvars, paste(deparse(xnames), collapse = ", ")), 
                  domain = NA)
            else {
                coefold <- start
                offset + as.vector(if (NCOL(x) == 1L) 
                  x * start
                else x %*% start)
            }
        else family$linkfun(mustart)
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta))) 
            stop("cannot find valid starting values: please specify some 1", 
                call. = FALSE)
        devold <- sum(dev.resids(y, mu, weights))
        boundary <- conv <- FALSE
        for (iter in 1L:control$maxit) {
            good <- weights > 0
            varmu <- variance(mu)[good]
            if (anyNA(varmu)) 
                stop("NAs in V(mu)")
            if (any(varmu == 0)) 
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
            good <- (weights > 0) & (mu.eta.val != 0)
            if (all(!good)) {
                conv <- FALSE
                warning(gettextf("no observations informative at iteration %d", 
                  iter), domain = NA)
                break
            }
            z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
            ngoodobs <- as.integer(nobs - sum(!good))
			
			## USE NNLS instead of least squares
			require(nnls)
            fit <- nnls(x[good, , drop = FALSE] *  w, z * w)
			fit$coefficients <- fit$x
			## calculate qr decomposition..
			QR <- qr(x[good, , drop = FALSE] *  w)
			fit$qr <- QR$qr
			fit$rank <- QR$rank
			fit$pivot <- QR$pivot
			fit$qraux <- QR$qraux	
			fit$effects <- fit$fitted
			### END NNLS
			
            if (any(!is.finite(fit$coefficients))) {
                conv <- FALSE
                warning(gettextf("non-finite coefficients at iteration %d", 
                  iter), domain = NA)
                break
            }
            if (nobs < fit$rank) 
                stop(sprintf(ngettext(nobs, "X matrix has rank %d, but only %d observation", 
                  "X matrix has rank %d, but only %d observations"), 
                  fit$rank, nobs), domain = NA)
            start[fit$pivot] <- fit$coefficients
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace) 
                cat("Deviance = ", dev, " Iterations - ", iter, 
                  "\n", sep = "")
            boundary <- FALSE
            if (!is.finite(dev)) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found: please supply starting values 2", 
                    call. = FALSE)
                warning("step size truncated due to divergence", 
                  call. = FALSE)
                ii <- 1
                while (!is.finite(dev)) {
                  if (ii > control$maxit) 
                    stop("inner loop 1; cannot correct step size", 
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace) 
                  cat("Step halved: new deviance = ", dev, "\n", 
                    sep = "")
            }
            if (!(valideta(eta) && validmu(mu))) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found: please supply starting values 3", 
                    call. = FALSE)
                warning("step size truncated: out of bounds", 
                  call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxit) 
                    stop("inner loop 2; cannot correct step size", 
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace) 
                  cat("Step halved: new deviance = ", dev, "\n", 
                    sep = "")
            }
            if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
                conv <- TRUE
                coef <- start
                break
            }
            else {
                devold <- dev
                coef <- coefold <- start
            }
        }
        if (!conv) 
            warning("glm.fit: algorithm did not converge", call. = FALSE)
        if (boundary) 
            warning("glm.fit: algorithm stopped at boundary value", 
                call. = FALSE)
        eps <- 10 * .Machine$double.eps
        if (family$family == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps)) 
                warning("glm.fit: fitted probabilities numerically 0 or 1 occurred", 
                  call. = FALSE)
        }
        if (family$family == "poisson") {
            if (any(mu < eps)) 
                warning("glm.fit: fitted rates numerically 0 occurred", 
                  call. = FALSE)
        }
        if (fit$rank < nvars) 
            coef[fit$pivot][seq.int(fit$rank + 1, nvars)] <- NA
        xxnames <- xnames[fit$pivot]
        residuals <- (y - mu)/mu.eta(eta)
        fit$qr <- as.matrix(fit$qr)
        nr <- min(sum(good), nvars)
        if (nr < nvars) {
            Rmat <- diag(nvars)
            Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars]
        }
        else Rmat <- fit$qr[1L:nvars, 1L:nvars]
        Rmat <- as.matrix(Rmat)
        Rmat[row(Rmat) > col(Rmat)] <- 0
        names(coef) <- xnames
        colnames(fit$qr) <- xxnames
        dimnames(Rmat) <- list(xxnames, xxnames)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    if (!EMPTY) 
        names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("", 
            sum(good) - fit$rank))
    wtdmu <- if (intercept) 
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY) 
        0
    else fit$rank
    resdf <- n.ok - rank
    aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
    list(coefficients = coef, residuals = residuals, fitted.values = mu, 
        effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat, 
        rank = rank, qr = if (!EMPTY) structure(fit[c("qr", "rank", 
            "qraux", "pivot", "tol")], class = "qr"), family = family, 
        linear.predictors = eta, deviance = dev, aic = aic.model, 
        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights, 
        df.residual = resdf, df.null = nulldf, y = y, converged = conv, 
        boundary = boundary)
}

##
##############################################################################################################################

##############################################################################################################################
##
## glm.fit function modified to allow forced non-negative and non-positive coefficient estimates
## 		see ?glm.fit for detail on glm.fit and ?nnnpls for details on nnnpls
##
nnnpls.fit <- function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL, 
    mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
    control = list(), intercept = TRUE, singular.ok = TRUE) 
{
    control <- do.call("glm.control", control)
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y)) 
        rownames(y)
    else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights)) 
        weights <- rep.int(1, nobs)
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("'family' argument seems not to be a valid family object", 
            call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null) if (is.null(x)) 
        if.null
    else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta)) 
            stop("invalid linear predictor values in empty model", 
                call. = FALSE)
        mu <- linkinv(eta)
        if (!validmu(mu)) 
            stop("invalid fitted means in empty model", call. = FALSE)
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric()
        iter <- 0L
    }
    else {
        coefold <- NULL
        eta <- if (!is.null(etastart)) 
            etastart
        else if (!is.null(start)) 
            if (length(start) != nvars) 
                stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", 
                  nvars, paste(deparse(xnames), collapse = ", ")), 
                  domain = NA)
            else {
                coefold <- start
                offset + as.vector(if (NCOL(x) == 1L) 
                  x * start
                else x %*% start)
            }
        else family$linkfun(mustart)
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta))) 
            stop("cannot find valid starting values: please specify some", 
                call. = FALSE)
        devold <- sum(dev.resids(y, mu, weights))
        boundary <- conv <- FALSE
        for (iter in 1L:control$maxit) {
            good <- weights > 0
            varmu <- variance(mu)[good]
            if (any(is.na(varmu))) 
                stop("NAs in V(mu)")
            if (any(varmu == 0)) 
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
            good <- (weights > 0) & (mu.eta.val != 0)
            if (all(!good)) {
                conv <- FALSE
                warning("no observations informative at iteration ", 
                  iter)
                break
            }
            z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
            ngoodobs <- as.integer(nobs - sum(!good))
			
			## USE NNNPLS instead of least squares
			require(nnls)
            fit <- nnnpls(x[good, , drop = FALSE] *  w, z * w,con=c(-1,rep(1,ncol(x)-1)))
			fit$coefficients <- fit$x
			## calculate qr decomposition..
			QR <- qr(x[good, , drop = FALSE] *  w,tol=min(1e-07, control$epsilon/1000))
			fit$qr <- QR$qr
			fit$rank <- QR$rank
			fit$pivot <- QR$pivot
			fit$qraux <- QR$qraux	
			fit$effects <- fit$fitted
			### END NNNPLS
			
            if (any(!is.finite(fit$coefficients))) {
                conv <- FALSE
                warning(gettextf("non-finite coefficients at iteration %d", 
                  iter), domain = NA)
                break
            }
            if (nobs < fit$rank) 
                stop(gettextf("X matrix has rank %d, but only %d observations", 
                  fit$rank, nobs), domain = NA)
            start[fit$pivot] <- fit$coefficients ## change here - remove $qr$
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace) 
                cat("Deviance =", dev, "Iterations -", iter, 
                  "\n")
            boundary <- FALSE
            if (!is.finite(dev)) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found: please supply starting values", 
                    call. = FALSE)
                warning("step size truncated due to divergence", 
                  call. = FALSE)
                ii <- 1
                while (!is.finite(dev)) {
                  if (ii > control$maxit) 
                    stop("inner loop 1; cannot correct step size", 
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
            if (!(valideta(eta) && validmu(mu))) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found: please supply starting values", 
                    call. = FALSE)
                warning("step size truncated: out of bounds", 
                  call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxit) 
                    stop("inner loop 2; cannot correct step size", 
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
            if (((dev - devold)/(0.1 + abs(dev)) >= control$epsilon) & 
                (iter > 1)) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found: please supply starting values", 
                    call. = FALSE)
                warning("step size truncated due to increasing deviance", 
                  call. = FALSE)
                ii <- 1
                while ((dev - devold)/(0.1 + abs(dev)) > -control$epsilon) {
                  if (ii > control$maxit) 
                    break
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                }
                if (ii > control$maxit) 
                  warning("inner loop 3; cannot correct step size")
                else if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
            if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
                conv <- TRUE
                coef <- start
                break
            }
            else {
                devold <- dev
                coef <- coefold <- start
            }
        }
        if (!conv) 
            warning("glm.fit2: algorithm did not converge. Try increasing the maximum iterations", 
                call. = FALSE)
        if (boundary) 
            warning("glm.fit2: algorithm stopped at boundary value", 
                call. = FALSE)
        eps <- 10 * .Machine$double.eps
        if (family$family == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps)) 
                warning("glm.fit2: fitted probabilities numerically 0 or 1 occurred", 
                  call. = FALSE)
        }
        if (family$family == "poisson") {
            if (any(mu < eps)) 
                warning("glm.fit2: fitted rates numerically 0 occurred", 
                  call. = FALSE)
        }
        if (fit$rank < nvars) 
            coef[fit$pivot][seq.int(fit$rank + 1, nvars)] <- NA 
        xxnames <- xnames[fit$pivot] 
        residuals <- (y - mu)/mu.eta(eta)
        fit$qr <- as.matrix(fit$qr) 
        nr <- min(sum(good), nvars)
        if (nr < nvars) {
            Rmat <- diag(nvars)
            Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars] 
        }
        else Rmat <- fit$qr[1L:nvars, 1L:nvars] 
        Rmat <- as.matrix(Rmat)
        Rmat[row(Rmat) > col(Rmat)] <- 0
        names(coef) <- xnames
        colnames(fit$qr) <- xxnames 
        dimnames(Rmat) <- list(xxnames, xxnames)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    if (!EMPTY) 
        names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("", 
            sum(good) - fit$rank))
    wtdmu <- if (intercept) 
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY) 
        0
    else fit$rank
    resdf <- n.ok - rank
    aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
    list(coefficients = coef, residuals = residuals, fitted.values = mu, 
        effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat, 
        rank = rank, qr = if (!EMPTY) structure(fit[c("qr", 
            "rank", "qraux", "pivot", "tol")], class = "qr"), 
        family = family, linear.predictors = eta, deviance = dev, 
        aic = aic.model, null.deviance = nulldev, iter = iter, 
        weights = wt, prior.weights = weights, df.residual = resdf, 
        df.null = nulldf, y = y, converged = conv, boundary = boundary)  
}

##
##############################################################################################################################

##############################################################################################################################
##
## Plotting function to produce visualisations of model fit 
##
## Inputs:
## 		formula - intercept value 
## 		data    - precalculated weighting value. See eq. 4 and supporting text for details
##
## Outputs:
## 		standard glm model object fitted to splined data
##
obs.gdm.plot <- function(model,title,w,Is,raw){

	## setup plotting parameters
	par(mfcol=c(2,2))
	## binned observed vs predicted
	plot(fitted(model), model$y, xlab = "Predicted proportion of species mismatch",ylab = "Observed proportion of species mismatch", type = "n",ylim=c(0,1),xlim=c(0,1))	
	sqs <- seq(0.05,0.95,by=0.1)
	rat <-unlist(lapply(sqs,function(x){
		data <- model$y[fitted(model) >= (x-0.05) & fitted(model) <= (x+0.05)]
		length(data[data==1]) / length(data)
	}))	
	lgth <-unlist(lapply(sqs,function(x){
		data <- model$y[fitted(model) >= (x-0.05) & fitted(model) <= (x+0.05)]
		length(data)
	}))
	rat_var <- rat*(1-rat)
	points(sqs, rat, pch = 20, cex = 2,col = rgb(0,0,1))
	segments(x0=sqs,y0=rat+rat_var,x1=sqs,y1=rat-rat_var,col=rgb(0,0,1))
	overlayX <- overlayY <- seq(from = min(fitted(model)),to = max(fitted(model)), length = 200)
	lines(overlayX, overlayY, lwd = 1,lty="dashed")
	## observed vs predicted
	plot(model$linear.predictors, model$y, xlab = "Predicted Ecological Distance",ylab = "Proportion of species mismatch", type = "n")
   	overlayX <- seq(from = min(model$linear.predictors), to = max(model$linear.predictors),length = 200)
    overlayY <- inv.logit(overlayX) 	
	ecoR <- range(model$linear.predictors)

	tt <- try(seq(ecoR[1],ecoR[2],by=0.1),silent=T)
	if(!inherits(tt, "try-error")){
		sqs <- seq(ecoR[1],ecoR[2],by=0.1)
	}
	if(inherits(tt, "try-error")){
		sqs <- seq(ecoR[1],10,by=0.1)
	}	
	
	rat <-unlist(lapply(sqs,function(x){
		data <- model$y[model$linear.predictors >= (x-0.05) & model$linear.predictors <= (x+0.05)]
		length(data[data==1]) / length(data)
	}))
	rat_var <- rat*(1-rat)	
	points(sqs, rat, pch = 20, cex = 1,col = rgb(0,0,1))
	segments(x0=sqs,y0=rat+rat_var,x1=sqs,y1=rat-rat_var,col=rgb(0,0,1))
	lines(overlayX, overlayY, lwd = 2,lty="dashed",col="green")
	legend("bottomright",legend=c("Observed","Predicted"),col=c("blue","green"),lty=c("solid","dashed"),pch=c(16,NA))
	## density of observed
	match <- density(x=model$linear.predictors[model$y == 0])
	miss <- density(x=model$linear.predictors[model$y == 1])
	mx <- max(c(match$y,miss$y))
	xrng <- range(model$linear.predictors)
	plot(1,1,type="n",xlim=xrng,ylim=c(0,mx),xlab="Predicted Ecological Distance",ylab="Density")
	polygon(match,col=rgb(1,0,0,0.3))
	polygon(miss,col=rgb(0,0,1,0.3))
	lines(match,col="red")
	lines(miss,col="blue")
	legend("topright",legend=c("Species match","Species mismatch"),col=c("red","blue"),lty=1)
	## calculate SORENSON vs binned obs match ratio
	brk <- seq(0.0,1,by=0.01)
	
	## crude richness threshold on sites for calculating SORENSON.
	raw2 <- raw
	r_thr <- 15
	raw2 <- raw[raw$Richness.S1 > r_thr & raw$Richness.S2 > r_thr,]
	binnedSor <- lapply(brk,function(x){length(raw2$Match[raw2$Match == 1 & raw2$SORENSON >= (x-0.05) & raw2$SORENSON <= (x+0.05)]) / length(raw2$Match[raw2$SORENSON >= (x-0.05) & raw2$SORENSON <= (x+0.05)])})
	raw2 <- raw[raw$TARGET,]
	raw2 <- raw2[raw2$Richness.S1 > r_thr & raw2$Richness.S2 > r_thr,]
	binnedSor2 <- lapply(brk,function(x){length(raw2$Match[raw2$Match == 1 & raw2$SORENSON >= (x-0.05) & raw2$SORENSON <= (x+0.05)]) / length(raw2$Match[raw2$SORENSON >= (x-0.05) & raw2$SORENSON <= (x+0.05)])})
	dat2 <- data.frame(y=brk,x=unlist(binnedSor),x2=unlist(binnedSor2))
	dat2 <- dat2[dat2$x != 0,]

	if(!(all(is.na(dat2$x)) | all(is.na(dat2$y)))){
		plot(dat2$y~dat2$x,ylab="Dissimilarity",xlab="Proportion of mismatches",type="n",xlim=c(0,1),ylim=c(0,1))
		points(dat2$x, dat2$y, pch = 20, cex = 1,col = rgb(0,0,1))
		rat_var <- dat2$x*(1-dat2$x)
		segments(y0=dat2$y,x0=dat2$x+rat_var,y1=dat2$y,x1=dat2$x-rat_var,col=rgb(0,0,1))
		cols <- c(rgb(0,0,1),"green")

		p0 <- inv.logit(Is)
		modOut <- seq(p0,1,length.out=100)
		tran <- ObsTrans(p0,w,modOut)
		lines(modOut,tran$out,col=cols[2],lty="dashed",lwd=2)
		legend("bottomright",legend=c("Observed","Predicted"),col=cols,lty=c("solid","dashed"),pch=c(16,NA))

	}

}


##############################################################################################################################
##
## Plotting function to produce visualizations of model fit 
##
## Inputs:
## 		formula - intercept value 
## 		data    - precalculated weighting value. See eq. 4 and supporting text for details
##
## Outputs:
## 		standard glm model object fitted to splined data
##
obs.gdm.plot.new <- function(model, title, w, Is, r_thr, raw){
  
  ## setup plotting parameters
  par(mfcol=c(2,2))
  ## binned observed vs predicted
  plot(fitted(model), model$y, xlab = "Predicted proportion of species mismatch",ylab = "Observed proportion of species mismatch", type = "n",ylim=c(0,1),xlim=c(0,1))	
  sqs <- seq(0.05,0.95,by=0.1)
  rat <-unlist(lapply(sqs,function(x){
    data <- model$y[fitted(model) >= (x-0.05) & fitted(model) <= (x+0.05)]
    length(data[data==1]) / length(data)
  }))	
  lgth <-unlist(lapply(sqs,function(x){
    data <- model$y[fitted(model) >= (x-0.05) & fitted(model) <= (x+0.05)]
    length(data)
  }))
  rat_var <- rat*(1-rat)
  points(sqs, rat, pch = 20, cex = 2,col = rgb(0,0,1))
  segments(x0=sqs,y0=rat+rat_var,x1=sqs,y1=rat-rat_var,col=rgb(0,0,1))
  overlayX <- overlayY <- seq(from = min(fitted(model)),to = max(fitted(model)), length = 200)
  lines(overlayX, overlayY, lwd = 1,lty="dashed")
  ## observed vs predicted
  plot(model$linear.predictors, model$y, xlab = "Predicted Ecological Distance",ylab = "Proportion of species mismatch", type = "n")
  overlayX <- seq(from = min(model$linear.predictors), to = max(model$linear.predictors),length = 200)
  overlayY <- inv.logit(overlayX) 	
  ecoR <- range(model$linear.predictors)
  
  tt <- try(seq(ecoR[1],ecoR[2],by=0.1),silent=T)
  if(!inherits(tt, "try-error")){
    sqs <- seq(ecoR[1],ecoR[2],by=0.1)
  }
  if(inherits(tt, "try-error")){
    sqs <- seq(ecoR[1],10,by=0.1)
  }	
  
  rat <-unlist(lapply(sqs,function(x){
    data <- model$y[model$linear.predictors >= (x-0.05) & model$linear.predictors <= (x+0.05)]
    length(data[data==1]) / length(data)
  }))
  rat_var <- rat*(1-rat)	
  points(sqs, rat, pch = 20, cex = 1,col = rgb(0,0,1))
  segments(x0=sqs,y0=rat+rat_var,x1=sqs,y1=rat-rat_var,col=rgb(0,0,1))
  lines(overlayX, overlayY, lwd = 2,lty="dashed",col="green")
  legend("bottomright",legend=c("Observed","Predicted"),col=c("blue","green"),lty=c("solid","dashed"),pch=c(16,NA))
  ## density of observed
  match <- density(x=model$linear.predictors[model$y == 0])
  miss <- density(x=model$linear.predictors[model$y == 1])
  mx <- max(c(match$y,miss$y))
  xrng <- range(model$linear.predictors)
  plot(1,1,type="n",xlim=xrng,ylim=c(0,mx),xlab="Predicted Ecological Distance",ylab="Density")
  polygon(match,col=rgb(1,0,0,0.3))
  polygon(miss,col=rgb(0,0,1,0.3))
  lines(match,col="red")
  lines(miss,col="blue")
  legend("topright",legend=c("Species match","Species mismatch"),col=c("red","blue"),lty=1)
  ## calculate SORENSON vs binned obs match ratio
  brk <- seq(0.0,1,by=0.01)
  
  ## crude richness threshold on sites for calculating SORENSON.
  raw2 <- raw
  r_thr <- r_thr
  raw2 <- raw[raw$Richness.S1 > r_thr & raw$Richness.S2 > r_thr,]
  binnedSor <- lapply(brk,function(x){length(raw2$Match[raw2$Match == 1 & raw2$SORENSON >= (x-0.05) & raw2$SORENSON <= (x+0.05)]) / length(raw2$Match[raw2$SORENSON >= (x-0.05) & raw2$SORENSON <= (x+0.05)])})
  raw2 <- raw[raw$TARGET,]
  raw2 <- raw2[raw2$Richness.S1 > r_thr & raw2$Richness.S2 > r_thr,]
  binnedSor2 <- lapply(brk,function(x){length(raw2$Match[raw2$Match == 1 & raw2$SORENSON >= (x-0.05) & raw2$SORENSON <= (x+0.05)]) / length(raw2$Match[raw2$SORENSON >= (x-0.05) & raw2$SORENSON <= (x+0.05)])})
  dat2 <- data.frame(y=brk,x=unlist(binnedSor),x2=unlist(binnedSor2))
  dat2 <- dat2[dat2$x != 0,]
  
  if(!(all(is.na(dat2$x)) | all(is.na(dat2$y)))){
    plot(dat2$y~dat2$x,ylab="Dissimilarity",xlab="Proportion of mismatches",type="n",xlim=c(0,1),ylim=c(0,1))
    points(dat2$x, dat2$y, pch = 20, cex = 1,col = rgb(0,0,1))
    rat_var <- dat2$x*(1-dat2$x)
    segments(y0=dat2$y,x0=dat2$x+rat_var,y1=dat2$y,x1=dat2$x-rat_var,col=rgb(0,0,1))
    cols <- c(rgb(0,0,1),"green")
    
    p0 <- inv.logit(Is)
    modOut <- seq(p0,1,length.out=100)
    tran <- ObsTrans(p0,w,modOut)
    lines(modOut,tran$out,col=cols[2],lty="dashed",lwd=2)
    legend("bottomright",legend=c("Observed","Predicted"),col=cols,lty=c("solid","dashed"),pch=c(16,NA))
    
  }
  
}

##
##############################################################################################################################
