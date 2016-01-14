#help=Efficient Global Optimization (EGO)
#type=Optimization
#output=Optimum
#options=initBatchSize=4,batchSize=4,iterations=10,bounds=true,trend=y~1,covtype=matern3_2,liar=max,search_min=false
#require=lhs,DiceKriging,DiceView,pso

#' constructor and initializer of R session
EGO <- function(options) {
  library(lhs)
  library(DiceKriging)
  library(DiceView)
  library(pso)

  # all parameters are initialy strings, so you have to put as global non-string values
  options$initBatchSize <- as.integer(options$initBatchSize)
  options$batchSize <- as.integer(options$batchSize)
  options$iterations <- as.integer(options$iterations)
  options$bounds <- as.logical(options$bounds)
  options$trend <- as.formula(options$trend)
  options$search_min <- as.logical(options$search_min)

  ego = new.env()
  ego$i = 0
  lapply(names(options),function(x) assign(x,options[[x]],ego))
  return(ego)
}

#' first design building. All variables are set in [0,1].
#' @param d number of variables
#' @return next design of experiments
getInitialDesign <- function(ego,d) {
  set.seed(1)
  if(ego$initBatchSize < 100){
    lhs <- optimumLHS(n=ego$initBatchSize,k=d)
  }else{
    lhs <- maximinLHS(n=ego$initBatchSize,k=d)
  }
  if (ego$bounds) {
    e=c(0,1)
    id=1
    while(id<d){
      e=rbind(cbind(e,0),cbind(e,1))
      id=id+1
    }
    Xinit=rbind(as.matrix(e),as.matrix(lhs))
  } else {
    Xinit=as.matrix(lhs)
  }
  return(Xinit)
}

#' iterated design building.
#' @param X data frame of current doe variables (in [0,1])
#' @param Y data frame of current results
#' @return  next doe step
getNextDesign <- function(ego,X,Y) {
  if (ego$i > ego$iterations) return();

  d = dim(X)[2]
  if (dim(Y)[2] == 2) {
    noise.var <- as.array(Y[,2])^2
  } else {
    noise.var <- NULL
  }

  if (ego$search_min) {y=Y[,1]} else {y=-Y[,1]}

  ego$kmi <- km(control=list(trace=FALSE),ego$trend,optim.method='BFGS',covtype=ego$covtype, noise.var = noise.var,design=X,response=y)

  EGOi <- max_qEI(model=ego$kmi,npoints=ego$batchSize,L=ego$liar,lower=rep(0,d),upper=rep(1,d),control=list(trace=FALSE))
  if (is.null(EGOi)) return()

  Xnext <- EGOi$par

  ego$i <- ego$i + 1
  return(as.matrix(Xnext))
}

#' final analysis.
#' @param X data frame of doe variables (in [0,1])
#' @param Y data frame of results
#' @return HTML string of analysis
displayResults <- function(ego,X,Y) {
  analyse.files <- paste("sectionview_",ego$i-1,".png",sep="")
  resolution <- 600

  if (dim(Y)[2] == 2) {
    noise.var <<- as.array(Y[,2])^2
    yname=paste0("N(",names(Y)[1],",",names(Y)[2])
  } else {
    noise.var <<- NULL
    yname=names(Y)
  }

  if (ego$search_min) {
    m = min(Y[,1])
    x = as.matrix(X)[which(Y[,1]==m),]
    html=paste(sep="<br/>",paste("<HTML>minimum is ",m),paste(sep="","found at ",paste(collapse="<br/>",paste(sep="= ",names(x),x)),"<br/><img src='",analyse.files,"' width='",resolution,"' height='",resolution,"'/></HTML>"))
  } else {
    m = max(Y[,1])
    x = as.matrix(X)[which(Y[,1]==m),]
    html=paste(sep="<br/>",paste("<HTML>maximum is ",m),paste(sep="","found at ",paste(collapse="<br/>",paste(sep="= ",names(x),x)),"<br/><img src='",analyse.files,"' width='",resolution,"' height='",resolution,"'/></HTML>"))
  }

  png(file=analyse.files,bg="transparent",height=resolution,width = resolution)
  try(sectionview.km(ego$kmi,center=x,Xname=names(X),yname=yname))
  dev.off()

  return(html)
}

################################################################################

#' @test X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); EI(runif(100),kmi)
#' @test X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); DiceView::sectionview.fun(function(x)EI(x,kmi),dim=1)
#' @test X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); DiceView::contourview.fun(function(x)EI(x,kmi),dim=2)
EI <- function (x, model, plugin=NULL) {
  if (is.null(plugin)){ if (model@noise.flag) plugin <- min(model@y-2*sqrt(model@noise.var)) else plugin <- min(model@y) }
  m <- plugin

  ########################################################################################
  # Convert x in proper format(s)
  if (!is.matrix(x)) x <- matrix(x,ncol= model@d)
  d <- ncol(x)
  if (d != model@d){ stop("x does not have the right number of columns (",d," instead of ",model@d,")") }
  newdata <- x
  colnames(newdata) = colnames(model@X)

  ########################################################################################
  #cat("predict...")
  predx <- predict.km(object=model, newdata=newdata, type="UK", checkNames = FALSE)
  #cat(" done.\n")
  kriging.mean <- predx$mean
  kriging.sd   <- predx$sd

  xcr <- (m - kriging.mean)/kriging.sd

  xcr.prob <- pnorm(xcr)
  xcr.dens <- dnorm(xcr)
  res <- (m - kriging.mean) * xcr.prob + kriging.sd * xcr.dens

  too.close = which(kriging.sd/sqrt(model@covariance@sd2) < 1e-06)
  res[too.close] <- max(0,m - kriging.mean)

  return(res)
}

#' @test set.seed(1); X=matrix(runif(20),ncol=2); y=branin(X); kmi <- km(design=X,response=y); kmi=km(design=X,response=y); DiceView::contourview.fun(function(x)EI(x,kmi),dim=2); points(max_EI(kmi,lower=c(0,0),upper=c(1,1))$par)
max_EI <-function(model,  lower, upper, control=NULL) {

  d <- ncol(model@X)

  if (is.null(control$print.level)) control$print.level <- 1
  if (is.null(control$max.parinit.iter)) control$max.parinit.iter <- 10^d
  if(d<=6) N <- 10*2^d else N <- 100*d
  if (is.null(control$pop.size))  control$pop.size <- N
  if (is.null(control$solution.tolerance))  control$solution.tolerance <- 1e-15

  pars=NULL
  for (i in 1:d) pars=cbind(pars,matrix(runif(N,lower[i],upper[i]),ncol=1))

  #t=Sys.time()
  ei <- EI(pars,model)
  #print(capture.output(Sys.time()-t))
  print(cbind(pars,ei))

  good_start = which(ei==max(ei,na.rm=T))
  par0=matrix(pars[good_start[sample(1:length(good_start),1)],],nrow=1)

  o <- psoptim(par=par0,fn=function(x){
    EI(x,model)
  },lower=lower,upper=upper,
  control=list( fnscale=-1, trace=control$print.level,maxit=10*d))

  o$par <- t(as.matrix(o$par))
  colnames(o$par) <- colnames(model@X)
  o$value <- as.matrix(o$value)
  colnames(o$value) <- "EI"

  return(list(par=o$par, value=o$value, counts=o$counts,par.all=o$par.all))
}

#' @test set.seed(1); X=matrix(runif(20),ncol=2); y=apply(FUN=branin,X,1); kmi <- km(design=X,response=y);  kmi=km(design=X,response=y); DiceView::contourview.fun(function(x)EI(x,kmi),dim=2); points(max_qEI(kmi,npoints=5,L="upper95",lower=c(0,0),upper=c(1,1))$par)
max_qEI <- function(model, npoints, L,  lower, upper,  control=NULL, ...) {
  n1 <- nrow(model@X)
  for (s in 1:npoints) {
    oEGO <- max_EI(model=model, lower=lower, upper=upper, control, ...)

    if (distXmin(oEGO$par,model@X)<=prod(upper-lower)*1E-10) {
        warning("Proposed a point already in design !");
        npoints=s-1;
        break;
    }

    model@X <- rbind(model@X, oEGO$par)

    if (L=="min")
      l = min(model@y)
    else if (L=="max")
      l = max(model@y)
    else if (L=="upper95")
      l = predict.km(object = model,newdata = oEGO$par,type="UK",light.return = TRUE)$upper95
    else if (L=="lower95")
      l = predict.km(object = model,newdata = oEGO$par,type="UK",light.return = TRUE)$lower95
    else l = L

    model@y <- rbind(model@y, l, deparse.level=0)

    model@F <- trendMatrix.update(model, Xnew=data.frame(oEGO$par))
    if (model@noise.flag) {
      model@noise.var = c(model@noise.var, 0) # here is the fix!
    }
    newmodel = NULL
    try(newmodel <- computeAuxVariables(model))
    if (is.null(newmodel)) {warning("Unable to update model !");npoints=s-1;break;}
    model = newmodel
  }

  if (npoints==0) return()
  return(list(par = model@X[(n1+1):(n1+npoints),, drop=FALSE], value = model@y[(n1+1):(n1+npoints),, drop=FALSE]))
}

distXmin <- function(x,Xmin) {
  return(min(sqrt(rowSums((Xmin-matrix(x,nrow=nrow(Xmin),ncol=ncol(Xmin),byrow=TRUE))^2))))
}

