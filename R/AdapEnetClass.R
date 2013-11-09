
plotObsEst <- function(yObs, yEst, delta, xlab = NULL,
ylab = NULL, title = NULL, legendplot = TRUE,
legendpos = "topleft", maxvalue = NULL, minvalue = NULL)
{
if(is.null(maxvalue))
maxvalue <- max(yObs, yEst)
if(is.null(minvalue))
minvalue <- min(yObs, yEst)
lim <- c(minvalue, maxvalue)
plot(yObs[delta == 1], yEst[delta == 1], pch = 1, xlim =
lim, ylim = lim, col = 2, main = title, xlab=xlab, ylab=ylab)
abline(0, 1)
points(yObs[delta == 0], yEst[delta == 0], pch = 4, col
= 4)
if(legendplot)
{
leg.txt <- c("Uncensor", "Censor")
legend(legendpos, leg.txt, pch = c(1,4), col = c(2, 4), cex = 3/4)
}
}


Enet.wls<-function(X, Y, delta)
{
n <- nrow(X) # number of samples
p <- ncol(X) # number of predictors
if(n != length(delta) || n != length(Y))
stop("dimensions of X, Y and delta don't match!")
kw <- aft.kmweight(Y,delta)$kmwts
XW <- apply(as.matrix(X[delta == 1,] * kw[delta == 1]), 2,
sum) / sum(kw[delta == 1])
YW <- sum(Y[delta == 1] * kw[delta == 1]) /
sum(kw[delta == 1])
for(i in 1:n)
X[i,] <- X[i,] - XW
X <- as.matrix(sqrt(kw) * X)
Y <- sqrt(kw) * (Y - YW)
fit<-glmnet(X, Y)
beta<-fit$beta[,dim(fit$beta)[2]]
return(list(beta=beta, fit=fit))
}


mrbj<-function(formula, data, subset, trace=FALSE, gehanonly=FALSE, cov=FALSE, 
na.action=na.exclude, residue=FALSE, mcsize=100)
    {
	lss.betag<-function(x,y,delta,z)
	{
	row=ncol(x)
	col=ncol(z)
	betagm<-matrix(0,ncol=col,nrow=row)
	dimnum<-dim(x)
	n1<-dimnum[1]
	n2<-dimnum[2]	
	yy0<-rep(y,rep(n1,n1))
	delta1<-rep(delta,rep(n1,n1))
	yy1<-rep(y,n1)
	yy2<-delta1*(yy0-yy1)	
	xx0<-matrix(rep(as.vector(x),rep(n1,n1*n2)),nrow=n1*n1)
	xx1<-t(matrix(rep(as.vector(t(x)),n1),nrow=n2))
	xx2<-xx0-xx1
	
	for(i in 1:col)
	{
		zz=rep(z[,i],rep(n1,n1))*rep(z[,i],n1)
		xxdif<-xx2*zz*delta1
		xnew<-apply(xxdif,2,sum)
		xnew<-rbind(xxdif)
		yynew<-c(yy2*zz)
            fit <- Enet.wls(xnew, yynew, delta1)$beta
    		betagm[,i] <- fit
	}
	betagm
	}

	eps <- .Machine$double.eps^(2/3)
	call <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
	mf <- mf[c(1, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1]] <- as.name("model.frame")
	mf <- eval(mf, sys.parent())
	Terms <- attr(mf, "terms")
	xvars <- as.character(attr(Terms, "variables"))
	yvar <- attr(Terms, "response")
	if((yvar <- attr(Terms, "response")) > 0)
		xvars <- xvars[ - yvar]
	else xlevels <- NULL
	y <- model.extract(mf, "response")
	x <- model.matrix(Terms, mf)

	if(all(x[, 1] == 1))
		x <- x[, -1]	
		
	if(ncol(as.matrix(y)) != 2)
		stop("Response must be a right-censored survival object!")

	nobs <- nrow(y)
	#nvar <- ncol(x)
      nvar1 <- ncol(x)

	fit <- list(converged = FALSE, gehanonly=gehanonly, cov=cov, mcsize=mcsize)

	fit$call <- call
	fit$nobs <- nobs
	fit$censored <- nobs - sum(y[,2])
	fit$niter <- 0
	fit$printkm <- residue
	if(gehanonly)
		{
	z <- matrix(rexp(nobs*mcsize), ncol=mcsize)
	zdummy <- matrix(rep(1,nobs), ncol=1)

	beta <- lss.betag(x, y[,1], y[,2], zdummy)
	betastar <- lss.betag(x, y[,1], y[,2], z)

	fit$betag <- beta
      beta <- lss.betag(x, y[,1], y[,2], zdummy)
	betastar <- lss.betag(x, y[,1], y[,2], z)
	fit$cnames <- dimnames(x)[[2]]
      nvar <- ncol(x)	
      bbar <- apply(betastar, 1, mean)
	tmp <- betastar - bbar
	fit$gehancov <- tmp %*% t(tmp)/(mcsize - 1)
	fit$gehansd <- sqrt(diag(fit$gehancov))
	fit$gehanzvalue <- beta/fit$gehansd
	fit$gehanpvalue <- (1 - pnorm(abs(fit$gehanzvalue))) * 2
	dimnames(fit$gehancov) <- list(fit$cnames,fit$cnames)
	if(trace)
		cat("\nbetag: ", format(beta), "\n\n")
       	}
	gnet<-Enet.wls(x, y[,1], y[,2])
	fit$enet <- gnet$beta
	fit$fit <- gnet$fit
	fit
     }



AEnet.aft<-function(X, Y, delta, weight, lambda2, maxit=10)
{
	n <- nrow(X) # number of samples
	p <- ncol(X) # number of predictors
	if(n != length(delta) || n != length(Y))
	stop("dimensions of X, Y and censorship don't match!")
	weight[weight==0]<-0.001
	w <-1/weight
	kw <- aft.kmweight(Y,delta)$kmwts
	XW <- apply(as.matrix(X[delta == 1,] * kw[delta == 1]), 2,
	sum) / sum(kw[delta == 1])
	YW <- sum(Y[delta == 1] * kw[delta == 1]) /
	sum(kw[delta == 1])
	for(i in 1:n)
	X[i,] <- X[i,] - XW
	X <- as.matrix(sqrt(kw) * X)
	Y <- sqrt(kw) * (Y - YW)
	Xextra<-diag(sqrt(lambda2), p)
	Yextra<-c(rep(0,p))
	delta.extra<-c(rep(1,p))

	X<-rbind(X, Xextra)
	for (j in 1:p)
	X[,j] <-X[,j]/w[j]
	Y<-c(Y, Yextra)
	delta<-c(delta, delta.extra)
	meanx <- apply(X,2,mean)
	normx <- sqrt(apply(X^2,2,sum))
	meany <- mean(Y)
	Y <- Y-meany
	X <- t(t(X)/normx)
    XX <- t(X)%*%X
    h <- rep(0,p)
    beta <- rbind(rep(0,p))
    cor <- t(X)%*%Y
    it <- 1
    potential <- active <- old.active <- rep(F,p)
    NZ <- rep(F,2*p)

get.h <-
		function(X,Y,beta,cor,active,it,XX,NZ){
    p <- ncol(X)
    K <- sum(active)
    XX <- matrix((XX*sign(cor))[active,],nrow=K)
    A.mat <-  cbind(-XX,XX,diag(K))
    NZ <- c(NZ,rep(F,K))
    exclude <- c(NZ[((p+1):(2*p))],NZ[1:p],rep(F,K))
    cor.sign <- sign(cor[which.min(-abs(cor))])
    if (K>1){
      B1.inv <- solve(A.mat[-K,NZ])
      B2 <- A.mat[K,NZ]
      B2.B1inv <- B2%*%B1.inv
      alpha <- as.vector(B2.B1inv%*%rep(1,K-1)-1)
      B1.BK <- B1.inv%*%A.mat[-K,!NZ & !exclude]
      q <- A.mat[K,!NZ& ! exclude]-B2%*%B1.BK
      cost <- c(rep(1,2*p),rep(0,K))
      eval <- alpha*(t(cost[NZ])%*%B1.BK-cost[!NZ & !exclude])/q
      eval[alpha*q<0] <- -10^100
      best <- which.min(-eval)
      new.col <- (1:(2*p+K))[!NZ & !exclude][best]
      old.cols <- (1:(2*p+K))[NZ]
      NZ[new.col] <- T
      q <- q[best]
      B.inv <- B1.inv*q+B1.BK[,best]%*%B2.B1inv
      B.inv <- rbind(cbind(B.inv,-B1.BK[,best]),cbind(-B2.B1inv,1))/q
      B.inv <- B.inv[order(c(old.cols,new.col)),]
    }
    else{
      NZ[which.min(-abs(cor))-p*(cor.sign-1)/2] <- T
      B.inv <- 1/A.mat[,NZ]}
    active.beta <- NZ[1:p] | NZ[(p+1):(2*p)]
    h <- rep(0,p)
    beta.pos <- (1:p)[NZ[1:p]]
    beta.neg <- (1:p)[NZ[(p+1):(2*p)]]
    h[c(beta.pos,beta.neg)] <- -(B.inv%*%rep(1,K))[1:sum(active.beta)]
    h[beta.neg] <- -h[beta.neg]
    slack <- 0
    if (sum(NZ[-(1:(2*p))])>0)
      slack <- (1:p)[active][NZ[-(1:(2*p))]]
    list(h=h,slack=slack,NZ=NZ[1:(2*p)])}

    active[which.max(abs(t(X)%*%Y))] <- T
    while (max(abs(cor))>10^-5 & it<=maxit){
      res <- Y-X%*%beta[it,]
      cor <- as.vector(t(X)%*%res)
      if (max(abs(cor))>10^-8){
        S <- sign(cor)
        max.cor <- max(abs(cor))
        obj <- get.h(X,Y,beta[it,],as.vector(cor),active,it,XX,NZ)
        NZ <- obj$NZ
        h <- obj$h
        active <- active & (1:p)!=obj$slack
        act <- min((1:p)[active])
        direc.1 <- t(h)%*%(XX[act,]-XX)
        direc.2 <- t(h)%*%(XX[act,]+XX)
        dist.hitzero <- -beta[it,]/h
        index2 <- (1:p)[h!=0][dist.hitzero[h!=0]>10^-11]
        dist1 <- (cor[act]-cor)/direc.1
        dist2 <- (cor[act]+cor)/direc.2
        index <- c((1:p)[!active][dist1[!active]>10^-11],
        (1:p)[!active][dist2[!active]>10^-11])
        index3 <- c((1:p)[!active][dist1[!active]>10^-11], ((p+1):(2*p))[!active][dist2[! active]>10^-11])
        dist <- c(dist1,dist2)[index3]
        if (length(index)==0)
          gamma1 <- 1/0
        else
          gamma1 <- min(dist)
        if (length(index2)==0)
          gamma2 <- 1/0
        else
          gamma2 <- min(dist.hitzero[index2])
        gamma3 <- min(c(dist1[act],dist2[act]),na.rm=T)
        gamma <- min(gamma1,gamma2,gamma3)
        if (gamma1<gamma2)
          active[index[which.min(dist)]] <- T
        else{NNZ <- (1:p)[index2][which.min(dist.hitzero[index2])]
              NZ[c(NNZ,NNZ+p)] <- F}
        beta <- rbind(beta,beta[it,]+h*gamma)
        it <- it+1
      }
    }
    beta=t(t(beta)/normx)
    for (j in 1:p)
    beta[,j] <-beta[,j]/w[j]
    beta <-(1+lambda2)*(beta)
    obj <- list(beta=beta,mu=meany,meanx=meanx,normx=normx,type="LASSO")
    class(obj) <- "lars"
    obj
  }


AEnetCC.aft<-function(X, Y, delta, weight, C=1,s = 1, lambda2)
{
n <- nrow(X) # number of samples
p <- ncol(X) # number of predictors
if(n != length(delta) || n != length(Y))
stop("dimensions of X, Y and censorship don't match!")
weight[weight==0]<-0.001
w <-1/weight
kw <- aft.kmweight(Y,delta)$kmwts
XW <- apply(as.matrix(X[delta == 1,] * kw[delta == 1]), 2,
sum) / sum(kw[delta == 1])
YW <- sum(Y[delta == 1] * kw[delta == 1]) /
sum(kw[delta == 1])
for(i in 1:n)
X[i,] <- X[i,] - XW
X <- as.matrix(sqrt(kw) * X)
Y <- sqrt(kw) * (Y - YW)

Xextra<-diag(sqrt(lambda2), p)
Yextra<-c(rep(0,p))
delta.extra<-c(rep(1,p))

X<-rbind(X, Xextra)
for (j in 1:p)
X[,j] <-X[,j]/w[j]
Y<-c(Y, Yextra)
delta<-c(delta, delta.extra)
meanx <- apply(X,2,mean)
normx <- sqrt(apply(X^2,2,sum))
meany <- mean(Y)
Y <- Y-meany
X <- t(t(X)/normx)

ncensor <- n+p - sum(delta)
nB <- 2 * p + ncensor
Dmat <- diag(C, nB)
DTmp <- cbind(X[delta==1, ], -X[delta==1, ])
Dmat[1:(2*p), 1:(2*p)] <- t(DTmp) %*% DTmp + diag(lambda2, 2 * p)
dvec <- c( drop(Y[delta==1]) %*% DTmp, numeric(ncensor) )
Aind <- matrix(0, 2*p + 2, nB+1)
Aind[1, 1:ncensor] <- rep(2*p+1, ncensor)
Aind[1, (ncensor+1):nB] <- rep(1, 2*p)
Aind[1, (nB+1)] <- 2*p
for(i in 1:(2*p))
Aind[(i+1), 1:ncensor] <- rep(i, ncensor)
Aind[(2*p+2), 1:ncensor] <- seq(2*p + 1, nB)
Aind[2, (ncensor+1):(ncensor+2*p)] <- seq(1, 2*p)
Aind[2:(2*p+1), nB+1] <- seq(1, 2*p)
Amat <- cbind( rbind(t(X[delta==0,]), t(-X[delta==0,]),
rep(1, ncensor)),
rbind(rep(1, 2*p), matrix(0, 2*p, 2*p)),
c(rep(-1, 2*p), 0) )
bvec <- c(Y[delta==0], numeric(2*p), -s)
qpobj <- solve.QP.compact(Dmat,dvec,Amat,Aind,bvec)
beta <- (qpobj$solution)[1:p] - (qpobj$solution)[(p+1):(2*p)]
beta=t(t(beta)/normx)
beta <-(1+lambda2)*(beta/w)
return( c(beta) )
}


WEnet.aft<-function(X, Y, delta, weight, lambda2, maxit=10)
{
n <- nrow(X) # number of samples
p <- ncol(X) # number of predictors
if(n != length(delta) || n != length(Y))
stop("dimensions of X, Y and censorship don't match!")
weight[weight==0]<-0.0001
w <-weight
kw <- aft.kmweight(Y,delta)$kmwts
XW <- apply(as.matrix(X[delta == 1,] * kw[delta == 1]), 2,
sum) / sum(kw[delta == 1])
YW <- sum(Y[delta == 1] * kw[delta == 1]) /
sum(kw[delta == 1])
for(i in 1:n)
X[i,] <- X[i,] - XW
X <- as.matrix(sqrt(kw) * X)
Y <- sqrt(kw) * (Y - YW)
Xextra<-diag(w*sqrt(lambda2), p)
Yextra<-c(rep(0,p))
delta.extra<-c(rep(1,p))
X<-(1+lambda2)^(-0.5)*rbind(X, Xextra)
for (j in 1:p)
X[,j] <-X[,j]/w[j]
Y<-c(Y, Yextra)
delta<-c(delta, delta.extra)
meanx <- apply(X,2,mean)
normx <- sqrt(apply(X^2,2,sum))
meany <- mean(Y)
Y <- Y-meany
X <- t(t(X)/normx)
    XX <- t(X)%*%X
    h <- rep(0,p)
    beta <- rbind(rep(0,p))
    cor <- t(X)%*%Y
    it <- 1
    potential <- active <- old.active <- rep(F,p)
    NZ <- rep(F,2*p)
		get.h <-
		function(X,Y,beta,cor,active,it,XX,NZ){
    p <- ncol(X)
    K <- sum(active)
    XX <- matrix((XX*sign(cor))[active,],nrow=K)
    A.mat <-  cbind(-XX,XX,diag(K))
    NZ <- c(NZ,rep(F,K))
    exclude <- c(NZ[((p+1):(2*p))],NZ[1:p],rep(F,K))
    cor.sign <- sign(cor[which.min(-abs(cor))])
    if (K>1){
      B1.inv <- solve(A.mat[-K,NZ])
      B2 <- A.mat[K,NZ]
      B2.B1inv <- B2%*%B1.inv
      alpha <- as.vector(B2.B1inv%*%rep(1,K-1)-1)
      B1.BK <- B1.inv%*%A.mat[-K,!NZ & !exclude]
      q <- A.mat[K,!NZ& ! exclude]-B2%*%B1.BK
      cost <- c(rep(1,2*p),rep(0,K))
      eval <- alpha*(t(cost[NZ])%*%B1.BK-cost[!NZ & !exclude])/q
      eval[alpha*q<0] <- -10^100
      best <- which.min(-eval)
      new.col <- (1:(2*p+K))[!NZ & !exclude][best]
      old.cols <- (1:(2*p+K))[NZ]
      NZ[new.col] <- T
      q <- q[best]
      B.inv <- B1.inv*q+B1.BK[,best]%*%B2.B1inv
      B.inv <- rbind(cbind(B.inv,-B1.BK[,best]),cbind(-B2.B1inv,1))/q
      B.inv <- B.inv[order(c(old.cols,new.col)),]
    }
    else{
      NZ[which.min(-abs(cor))-p*(cor.sign-1)/2] <- T
      B.inv <- 1/A.mat[,NZ]}
    active.beta <- NZ[1:p] | NZ[(p+1):(2*p)]
    h <- rep(0,p)
    beta.pos <- (1:p)[NZ[1:p]]
    beta.neg <- (1:p)[NZ[(p+1):(2*p)]]
    h[c(beta.pos,beta.neg)] <- -(B.inv%*%rep(1,K))[1:sum(active.beta)]
    h[beta.neg] <- -h[beta.neg]
    slack <- 0
    if (sum(NZ[-(1:(2*p))])>0)
      slack <- (1:p)[active][NZ[-(1:(2*p))]]
    list(h=h,slack=slack,NZ=NZ[1:(2*p)])}

    active[which.max(abs(t(X)%*%Y))] <- T
    while (max(abs(cor))>10^-5 & it<=maxit){
      res <- Y-X%*%beta[it,]
      cor <- as.vector(t(X)%*%res)
      if (max(abs(cor))>10^-8){
        S <- sign(cor)
        max.cor <- max(abs(cor))
        obj <- get.h(X,Y,beta[it,],as.vector(cor),active,it,XX,NZ)
        NZ <- obj$NZ
        h <- obj$h
        active <- active & (1:p)!=obj$slack
        act <- min((1:p)[active])
        direc.1 <- t(h)%*%(XX[act,]-XX)
        direc.2 <- t(h)%*%(XX[act,]+XX)
        dist.hitzero <- -beta[it,]/h
        index2 <- (1:p)[h!=0][dist.hitzero[h!=0]>10^-11]
        dist1 <- (cor[act]-cor)/direc.1
        dist2 <- (cor[act]+cor)/direc.2
        index <- c((1:p)[!active][dist1[!active]>10^-11],
        (1:p)[!active][dist2[!active]>10^-11])
        index3 <- c((1:p)[!active][dist1[!active]>10^-11], ((p+1):(2*p))[!active][dist2[! active]>10^-11])
        dist <- c(dist1,dist2)[index3]
        if (length(index)==0)
          gamma1 <- 1/0
        else
          gamma1 <- min(dist)
        if (length(index2)==0)
          gamma2 <- 1/0
        else
          gamma2 <- min(dist.hitzero[index2])
        gamma3 <- min(c(dist1[act],dist2[act]),na.rm=T)
        gamma <- min(gamma1,gamma2,gamma3)
        if (gamma1<gamma2)
          active[index[which.min(dist)]] <- T
        else{NNZ <- (1:p)[index2][which.min(dist.hitzero[index2])]
              NZ[c(NNZ,NNZ+p)] <- F}
        beta <- rbind(beta,beta[it,]+h*gamma)
        it <- it+1
      }
    }
    beta=t(t(beta)/normx)
    for (j in 1:p)
    beta[,j] <-beta[,j]/w[j]
    beta <-(1+lambda2)*beta
    obj <- list(beta=beta,mu=meany,meanx=meanx,normx=normx,type="LASSO")
    class(obj) <- "lars"
    obj
  }


WEnetCC.aft<-function(X, Y, delta, weight, C, s , lambda2)
{
n <- nrow(X) # number of samples
p <- ncol(X) # number of predictors
if(n != length(delta) || n != length(Y))
stop("dimensions of X, Y and censorship don't match!")

weight[weight==0]<-0.0001
w <-weight
gamma<-s/(sqrt(1+lambda2))

eta<- 0.01*sqrt(2*log(p))
kw <- aft.kmweight(Y,delta)$kmwts
XW <- apply(as.matrix(X[delta == 1,] * kw[delta == 1]), 2,
sum) / sum(kw[delta == 1])
YW <- sum(Y[delta == 1] * kw[delta == 1]) /
sum(kw[delta == 1])
for(i in 1:n)
X[i,] <- X[i,] - XW
X <- as.matrix(sqrt(kw) * X)
Y <- sqrt(kw) * (Y - YW)

Xextra<-diag(w*sqrt(lambda2), p)
Yextra<-c(rep(0,p))
delta.extra<-c(rep(1,p))

X<-(1+lambda2)^(-0.5)*rbind(X, Xextra)
for (j in 1:p)
X[,j] <-X[,j]/w[j]
Y<-c(Y, Yextra)
delta<-c(delta, delta.extra)
meanx <- apply(X,2,mean)
normx <- sqrt(apply(X^2,2,sum))
meany <- mean(Y)
Y <- Y-meany
X <- t(t(X)/normx)

ncensor <- n+p - sum(delta)
nB <- 2 * p + ncensor
Dmat <- diag(C, nB)
DTmp <- cbind(X[delta==1, ], -X[delta==1, ])
Dmat[1:(2*p), 1:(2*p)] <- t(DTmp) %*% DTmp + diag(lambda2, 2 * p)
#Dmat[1:(2*p), 1:(2*p)] <- t(DTmp) %*% DTmp
dvec <- c( drop(Y[delta==1]) %*% DTmp, numeric(ncensor) )
Aind <- matrix(0, 2*p + 2, nB+1)
Aind[1, 1:ncensor] <- rep(2*p+1, ncensor)
Aind[1, (ncensor+1):nB] <- rep(1, 2*p)
Aind[1, (nB+1)] <- 2*p
for(i in 1:(2*p))
Aind[(i+1), 1:ncensor] <- rep(i, ncensor)
Aind[(2*p+2), 1:ncensor] <- seq(2*p + 1, nB)
Aind[2, (ncensor+1):(ncensor+2*p)] <- seq(1, 2*p)
Aind[2:(2*p+1), nB+1] <- seq(1, 2*p)
Amat <- cbind( rbind(t(X[delta==0,]), t(-X[delta==0,]),
rep(1, ncensor)),
rbind(rep(1, 2*p), matrix(0, 2*p, 2*p)),
c(rep(-1, 2*p), 0) )
bvec <- c(Y[delta==0], numeric(2*p), -gamma)
qpobj <- solve.QP.compact(Dmat,dvec,Amat,Aind,bvec)
beta <- (qpobj$solution)[1:p] - (qpobj$solution)[(p+1):(2*p)]
beta=t(t(beta)/normx[1:p])
beta <-(1+lambda2)*(beta/w)
return( c(beta) )
}


cv.AWEnet <-
function (X, Y, delta, weight, lambda2, maxit, K = 10, fraction = seq(from = 0, to = 1, length = 100), 
            plot.it = F, se = TRUE, AEnet=T, all.folds=NULL) 
{
  #require(lars)
   bds<-sort(lambda2)
	cv.Enet <-
	function (X, Y, delta, weight, lambda2, maxit, K, fraction, AEnet) 
	{
  	if (is.null(all.folds))
    	all.folds <- cv.folds(length(Y), K)
  	residmat <- matrix(0, length(fraction), K)
  	for (i in seq(K)) {
    	omit <- all.folds[[i]]
    	if (AEnet)
      fit <- AEnet.aft(X, Y, delta, weight, lambda2, maxit)
    	else
      fit <- WEnet.aft(X, Y, delta, weight, lambda2, maxit)
    	fit <- predict(fit, X[omit,, drop = FALSE], mode = "fraction",s = fraction)$fit
    	if (length(omit) == 1) 
      fit <- matrix(fit, nrow = 1)
    	residmat[, i] <- apply((Y[omit] - fit)^2, 2, mean)
  	}
  	cv <- apply(residmat, 1, mean)
  	cv.error <- sqrt(apply(residmat, 1, var)/K)
  	object <- list(index = fraction, cv = cv, cv.error =  cv.error,all.folds=all.folds, mode="fraction")
  	if (plot.it) 
    	plotCVLars(object, se = se)
  	invisible(object)
		}

   index<-NULL
	for (lambda in bds)
	{
      if (AEnet){
	cvEnet<-cv.Enet(X, Y, delta, weight, lambda, maxit, K , fraction, AEnet=T) 
	s<-cvEnet$index[which.min(cvEnet$cv)]
	cv.mse<-which.min(cvEnet$cv)
	cv.error<-which.min(cvEnet$cv.error)
	index<-rbind(index, c(lambda, s, cv.mse ,cv.error))}
	else{
	cvEnet<-cv.Enet(X, Y, delta, weight, lambda, maxit, K , fraction, AEnet=F) 
	gama<-cvEnet$index[which.min(cvEnet$cv)]
	s<-gama*sqrt(1+lambda)
	cv.mse<-which.min(cvEnet$cv)
	cv.error<-which.min(cvEnet$cv.error)
	index<-rbind(index, c(lambda, s, cv.mse ,cv.error))}
	}
list(index=index)
}



cv.AWEnetCC<-function(X, Y, delta, weight, kFold = 10, C, s, lambda2, AEnetCC=T)
{
n <- nrow(X) # number of samples
p <- ncol(X) # number of predictors
if(n != length(delta) || n != length(Y))
stop("dimensions of X, Y and censorship don't match!")
aft.cv.censorkfold <- function(n, delta, kfold, maxTry = 10)
{
for (i in 1:maxTry) # try at most maxTry times
{
cvFold <- split(sample(1:n), rep(1:kfold, length = n))
kfoldsuccess <- T
for (k in 1:kfold)
{
if (sum(delta[cvFold[[k]]] == 1) <= 1)
{
kfoldsuccess <- F
break
}
}
if (kfoldsuccess)
break
}
if(!kfoldsuccess)

stop(message="Too small number of samples in the given number of
folds! Try decreasing number of folds!")
return(cvFold)
}

cvFold <- aft.cv.censorkfold(n, delta, kFold)
Betas <- numeric(p)
cvscore <- numeric(n)
for (i in 1:kFold)
{
foldIndex <- cvFold[[i]]
if(AEnetCC)
beta<-AEnetCC.aft(as.matrix(X[-foldIndex, ]), Y[-foldIndex], delta[-foldIndex], weight, C, s, lambda2)
else 
beta<-WEnetCC.aft(as.matrix(X[-foldIndex, ]), Y[-foldIndex], delta[-foldIndex], weight, C, s, lambda2)

Betas <- cbind(Betas, beta)
cvscore[foldIndex] <- Y[foldIndex] - as.matrix(X[foldIndex, ]) %*% beta
}
kw <- aft.kmweight(Y, delta)$kmwts * delta
cvscore <- sum((cvscore^2) * kw) / 2
list(beta = Re(apply(Betas[,], 1, mean)), betavar =
Re(apply(Betas[,], 1, var)), cvscore = Re(cvscore) )
}


