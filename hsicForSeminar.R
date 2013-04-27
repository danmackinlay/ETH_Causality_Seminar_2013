# Copyright (c) 2010-2012  Jonas Peters  [jopeters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms.

library(kernlab)

hsicTest <- function(x,y,alpha=0.05, pars = list())
{
    # outputs the test statistic (N*HSIC) and the critical value (according to alpha). If the test statistic is 
    # larger than the critical value, 
    # H_0 (X and Y are independent) is rejected.

    if(is.matrix(x)==FALSE){
    x<-as.matrix(x)}
    if(is.matrix(y)==FALSE){
    y<-as.matrix(y)}
    len <- dim(x)[1]

    # compute distance matrices
    if(len>1000)
    {
        xhilf<-sample(x,1000)
        yhilf<-sample(y,1000)
    }
    else
    {
        xhilf<-x
        yhilf<-y
    }
    xnorm<-as.matrix(dist(xhilf,method="euclidean",diag=TRUE,upper=TRUE))
    xnorm<-xnorm^2
    ynorm<-as.matrix(dist(yhilf,method="euclidean",diag=TRUE,upper=TRUE))
    ynorm<-ynorm^2

    # choose median heuristic for bandwidth
    sigmax<-sqrt(0.5*median(xnorm[lower.tri(xnorm,diag=FALSE)]))
    sigmay<-sqrt(0.5*median(ynorm[lower.tri(xnorm,diag=FALSE)]))

    if(sigmax == 0)
    {
      sigmax<-sqrt(0.5*mean(xnorm[lower.tri(xnorm,diag=FALSE)]))
    }
    if(sigmay == 0)
    {
      sigmay<-sqrt(0.5*mean(ynorm[lower.tri(xnorm,diag=FALSE)]))
    }
    
    ###
    # Compute HSIC
    ###

    ## exact, but slow HSIC
    # KX<-exp(-xnorm/(2*sigmax^2))
    # KY<-exp(-ynorm/(2*sigmay^2))
    # H<-diag(1,len)-1/len*matrix(1,len,len)
    # HSIC_slow<-1/(len^2)*sum(diag(KX%*%H%*%KY%*%H))

    ## incomplete cholesky decomposition
    LX<-inchol(x, kernel="rbfdot", kpar=list(sigma=1/(2*sigmax^2)), tol = 0.0001, maxiter = 300)
    LX<-matrix(LX,nrow=dim(LX)[1], ncol=dim(LX)[2])
    LY<-inchol(y, kernel="rbfdot", kpar=list(sigma=1/(2*sigmay^2)), tol = 0.0001, maxiter = 300)
    LY<-matrix(LY,nrow=dim(LY)[1], ncol=dim(LY)[2])

    LXc<-LX-1/len*(as.matrix(rep(1,len))%*%colSums(LX))
    LYc<-LY-1/len*(as.matrix(rep(1,len))%*%colSums(LY))

    # tr(H*LX*LX'*H*LY*LY')
    #=tr(LXc*LX'*LYc*LY')
    #=tr((LX'*LYc)*(LY'*LXc*))
    
    HSIC<-1/(len^2)*sum(diag((t(LX)%*%LYc)%*%(t(LY)%*%LXc)))


    ###
    # Compute Gamma Approximation
    ###

    mux<-(crossprod(colSums(LX))-len)/(len*(len-1))
    muy<-(crossprod(colSums(LY))-len)/(len*(len-1))

    mean_h0<-1/len*(1+mux*muy-mux-muy)
    var_h0<-(2*(len-4)*(len-5))/(len*(len-1)*(len-2)*(len-3))*1/((len-1)^2)*sum(diag((t(LX)%*%LXc)%*%(t(LX)%*%LXc)))*1/((len-1)^2)*sum(diag((t(LY)%*%LYc)%*%(t(LY)%*%LYc)))

    a<-(mean_h0^2)/var_h0
    b<-len*var_h0/mean_h0

    ress <- list()
    ress$statistic <- len*HSIC
    ress$criticalValue <- qgamma(1-alpha,shape=a,scale=b)
    ress$pValue <- 1-pgamma(len*HSIC,shape=a,scale=b)
    return(ress)
}
