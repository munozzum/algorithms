#help=Random Embedding Bayesian Optimization (REMBO)
#type=Optimization
#output=Optimum
#options=initBatchSize=4,batchSize=4,iterations=10,bounds=true,trend=y~1,covtype=matern3_2,liar=max,search_min=false,
#d=<Low dimension (integer)>,D=<Problem dimension (integer)>,low_size=sqrt(d),kernel=Warped,A_norm=true  
#require=lhs,DiceKriging,DiceView,pso,MASS

#' constructor and initializer of R session
REMBO <- function(options) {
  library(lhs)
  library(DiceKriging)
  library(DiceView)
  library(pso)
  library(MASS)
  
  # Rembo options: d (low dimension for Rembo),
  #                D (high dimension/ Initial problem dimension)
  #                low_size=[-sqrt(d), sqrt(d)]^d (default box size for the low dimensional search)
  #                kernel = "Warped" (kernel for Rembo: "lowDim", "highDim", "Warped" (default))
  #                A_norm = true (should the row of the embedding matrix be normalized)
  
  
  # all parameters are initialy strings, so you have to put as global non-string values
  options$initBatchSize <- as.integer(options$initBatchSize)
  options$batchSize <- as.integer(options$batchSize)
  options$iterations <- as.integer(options$iterations)
  options$bounds <- as.logical(options$bounds)
  options$trend <- as.formula(options$trend)
  options$search_min <- as.logical(options$search_min)
  options$d <- as.integer(options$d)
  options$D <- as.integer(options$D)
  options$low_size <- as.double(options$low_size)
  options$kernel <- as.character(options$kernel)
  options$A_norm <- as.logical(options$A_norm)
  
  A <- matrix(rnorm(options$D*options$d), options$D, options$d)
  if(options$A_norm)
    A <- A/sqrt(rowSums(A^2))
  
  options$A <- A
  
  if(options$kernel == "Warped"){
    options$pA <- A %*% ginv(t(A) %*% A) %*% t(A)
    options$WObs <- NULL # for storing warped dimensional points
  }
  
  options$lowObs <- numeric(0)  # for storing low dimensional points
  
  rembo = new.env()
  rembo$i = 0
  lapply(names(options),function(x) assign(x,options[[x]],rembo))
  return(rembo)
}

#' first design building. All variables are set in [0,1]. 
#' @param rembo a rembo env
#' @return next design of experiments
getInitialDesign <- function(rembo) {
  set.seed(1)
  if(rembo$initBatchSize < 100){
    lhs <- optimumLHS(n = rembo$initBatchSize, k = rembo$d)
  }else{
    lhs <- maximinLHS(n = rembo$initBatchSize, k = rembo$d)
  }
  if (rembo$bounds) {
    e=c(0,1)
    id=1
    while(id<rembo$d){
      e=rbind(cbind(e,0),cbind(e,1))
      id=id+1
    }
    Xinit=rbind(as.matrix(e),as.matrix(lhs))
  } else {
    Xinit=as.matrix(lhs)
  }
  
  ## Specific Rembo part
  #1) First resize to [-boxsize,boxsize]^d
  Yinit <- 2*Xinit - 1# design in the low dimensional space
  Xinit <- t(apply(Yinit, 1, mapping_to_X, rembo$A)) # design is the high dimensional one
  
  #2) Check that no replicates are present (coming from the convex projection) and eventually replace them
  if(any(duplicated(Xinit))){
    size <- nrow(Xinit) # required size
    Ytemp <- Yinit[!duplicated(Xinit),]
    while(nrow(Ytemp) < size){
      tmp <- size - nrow(Ytemp)
      Yinit <- augmentLHS(Yinit, tmp)
      Xinit <- rbind(Xinit, t(apply(Yinit[-c(1:nrow(Xinit)),,drop = FALSE], 1, mapping_to_X, rembo$A)))
      Ytemp <- Yinit[!duplicated(Xinit),]
    }
    Yinit <- Ytemp
    Xinit <- Xinit[!duplicated(Xinit),]
  }
  
  # assign("lowObs", Yinit, envir = rembo)
  rembo$lowObs <- Yinit
  if(rembo$kernel=="Warped"){
    rembo$WObs <- t(apply(Yinit, 1, Psi, A = rembo$A, pA = rembo$pA))
    # assign("WObs", WObs, envir = rembo)
  }
  return(Xinit)
}

#' iterated design building.
#' @param rembo a REMBO env
#' @param X data frame of current doe variables (in [0,1])
#' @param Y data frame of current results
#' @return  next doe step
getNextDesign <- function(rembo,X,Y) {
  if (rembo$i > rembo$iterations) return();
  
  d = rembo$d
  if (dim(Y)[2] == 2) {
    noise.var <- as.array(Y[,2])^2
  } else {
    noise.var <- NULL
  }
  
  if (rembo$search_min) {y=Y[,1]} else {y=-Y[,1]}
  
  # Different models depending on the kernel
  if(rembo$kernel == "lowDim"){
    rembo$kmi <- km(control = list(trace=FALSE), formula = rembo$trend, optim.method='BFGS',
                    covtype = rembo$covtype, noise.var = noise.var, design = rembo$lowObs, response=y,
                    iso = T)
  }else{
    if(rembo$kernel == "highDim"){
      rembo$kmi <- km(control = list(trace=FALSE), formula = rembo$trend, optim.method='BFGS',
                      covtype = rembo$covtype, noise.var = noise.var, design = X, response=y,
                      iso = T)
    }
    if(rembo$kernel == "Warped"){
      rembo$kmi <- km(control = list(trace=FALSE), formula = rembo$trend, optim.method='BFGS',
                      covtype = rembo$covtype, noise.var = noise.var, design = rembo$WObs, response=y,
                      iso = T)
    }
  }
  
  
  EGOi <- max_qEI_REMBO(model = rembo$kmi, npoints = rembo$batchSize, rembo = rembo, L = rembo$liar, 
                        lower = rep(-rembo$low_size, d),
                        upper = rep(rembo$low_size, d), control=list(trace=FALSE))
  if (is.null(EGOi)) return()
  
  rembo$lowObs <- rbind(rembo$lowObs, EGOi$par)
  if(rembo$kernel == "Warped"){
    rembo$WObs <- rbind(rembo$WObs, t(apply(EGOi$par, 1, Psi, A = rembo$A, pA = rembo$pA)))
  }
  
  Xnext <- t(apply(EGOi$par, 1, mapping_to_X, A = rembo$A))
  
  rembo$i <- rembo$i + 1
  return(as.matrix(Xnext))
}

#' final analysis.
#' @param rembo a rembo env
#' @param X data frame of doe variables (in [0,1])
#' @param Y data frame of results
#' @return HTML string of analysis
displayResults <- function(rembo,X,Y) {
  analyse.files <- paste("sectionview_",rembo$i-1,".png",sep="")
  resolution <- 600
  
  if (dim(Y)[2] == 2) {
    noise.var <<- as.array(Y[,2])^2
    yname=paste0("N(",names(Y)[1],",",names(Y)[2])
  } else {
    noise.var <<- NULL
    yname=names(Y)
  }
  
  if (rembo$search_min) {
    m = min(Y[,1])
    x = as.matrix(X)[which(Y[,1]==m),]
    if(rembo$kernel == "lowDim"){
      xr = as.matrix(rembo$lowObs)[which(Y[,1]==m),]
    }
    if(rembo$kernel == "highDim"){
      xr = as.matrix(X)[which(Y[,1]==m),]
    }
    if(rembo$kernel == "Warped"){
      xr = as.matrix(rembo$WObs)[which(Y[,1]==m),]
    }
    
    html=paste(sep="<br/>",paste("<HTML>minimum is ",m),paste(sep="","found at ",paste(collapse="<br/>",paste(sep="= ",names(x),x)),"<br/><img src='",analyse.files,"' width='",resolution,"' height='",resolution,"'/></HTML>"))
  } else {
    m = max(Y[,1])
    x = as.matrix(X)[which(Y[,1]==m),]
    if(rembo$kernel == "lowDim"){
      xr = as.matrix(rembo$lowObs)[which(Y[,1]==m),]
    }
    if(rembo$kernel == "highDim"){
      xr = as.matrix(X)[which(Y[,1]==m),]
    }
    if(rembo$kernel == "Warped"){
      xr = as.matrix(rembo$WObs)[which(Y[,1]==m),]
    }
    html=paste(sep="<br/>",paste("<HTML>maximum is ",m),paste(sep="","found at ",paste(collapse="<br/>",paste(sep="= ",names(x),x)),"<br/><img src='",analyse.files,"' width='",resolution,"' height='",resolution,"'/></HTML>"))
  }
  
  png(file=analyse.files,bg="transparent",height=resolution,width = resolution)

  try(sectionview.km(rembo$kmi,center=xr,Xname=names(rembo$kmi@X),yname=yname))
  dev.off()
  
  return(html)
}

################################################################################

#' @test X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); EI(runif(100),kmi)
#' @test X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); DiceView::sectionview.fun(function(x)EI(x,kmi),dim=1)
# #' @test X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); DiceView::contourview.fun(function(x)EI(x,kmi),dim=2)
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

# (not done yet)#' @test set.seed(1); X=matrix(runif(20),ncol=2); y=branin(X); kmi <- km(design=X,response=y); kmi=km(design=X,response=y); DiceView::contourview.fun(function(x)EI(x,kmi),dim=2); points(max_EI(kmi,lower=c(0,0),upper=c(1,1))$par)
max_EI_REMBO <-function(model, rembo, lower, upper, control=NULL) {
  
  d <- ncol(rembo$lowObs)
  
  if (is.null(control$print.level)) control$print.level <- 1
  if (is.null(control$max.parinit.iter)) control$max.parinit.iter <- 10^d
  if(d<=6) N <- 10*2^d else N <- 100*d
  if (is.null(control$pop.size))  control$pop.size <- N
  if (is.null(control$solution.tolerance))  control$solution.tolerance <- 1e-15
  
  pars=NULL
  for (i in 1:d) pars=cbind(pars,matrix(runif(N,lower[i],upper[i]),ncol=1))
  
  #t=Sys.time()
  if(rembo$kernel == "highDim"){
    par_W <- t(apply(pars, 1, mapping_to_X, rembo$A))
  }
  if(rembo$kernel == "Warped"){
    par_W <- t(apply(pars, 1, Psi, A = rembo$A, pA = rembo$pA))
  }
  if(rembo$kernel == "lowDim")
    par_W <- pars
  ei <- EI(par_W, model)
  
  #print(capture.output(Sys.time()-t))
  # print(cbind(pars,ei))
  
  good_start = which(ei==max(ei,na.rm=T))
  par0=matrix(pars[good_start[sample(1:length(good_start),1)],],nrow=1)
  
  o <- psoptim(par=par0, fn=function(x){
    if(rembo$kernel == "highDim"){
      x <- t(mapping_to_X(x, rembo$A))
    }
    if(rembo$kernel == "Warped"){
      x <- Psi(x, A = rembo$A, pA = rembo$pA)
    }
    EI(x,model)
  },lower = lower, upper = upper,
  control = list(fnscale = -1, trace = control$print.level, maxit = 10*d))
  
  o$par <- t(as.matrix(o$par))
  colnames(o$par) <- colnames(rembo$lowObs)
  o$value <- as.matrix(o$value)
  colnames(o$value) <- "EI"
  
  return(list(par=o$par, value=o$value, counts=o$counts,par.all=o$par.all))
}

# #' @test set.seed(1); X=matrix(runif(20),ncol=2); y=apply(FUN=branin,X,1); kmi <- km(design=X,response=y);  kmi=km(design=X,response=y); DiceView::contourview.fun(function(x)EI(x,kmi),dim=2); points(max_qEI(kmi,npoints=5,L="upper95",lower=c(0,0),upper=c(1,1))$par)
max_qEI_REMBO <- function(model, npoints, L, rembo,  lower, upper,  control=NULL, ...) {
  n1 <- nrow(model@X)
  
  newlowDimPoints <- NULL
  
  for (s in 1:npoints) {
    oEGO <- max_EI_REMBO(model=model, rembo = rembo, lower=lower, upper=upper, control, ...)
    
    newlowDimPoints <- rbind(newlowDimPoints, oEGO$par)
    newX <- t(mapping_to_X(as.vector(oEGO$par), rembo$A)) # replace by a apply?
    
    if(rembo$kernel == "lowDim"){
      if(distXmin(oEGO$par,model@X)<=prod(upper-lower)*1E-10) {
        warning("Proposed a point already in design !");
        npoints=s-1;
        break;
      }
      newPointKrig <- oEGO$par
    }
    
    if(rembo$kernel == "highDim"){
      if (distXmin(newX,model@X)<=prod(upper-lower)*1E-10) {
        warning("Proposed a point already in design !");
        npoints=s-1;
        break;
      }
      newPointKrig <- newX
    }
    
    if(rembo$kernel == "Warped"){
      newW <- Psi(as.vector(oEGO$par), A = rembo$A, rembo$pA)
      if (distXmin(newW,model@X)<=prod(upper-lower)*1E-10) {
        warning("Proposed a point already in design !");
        npoints=s-1;
        break;
      }
      newPointKrig <- newW
    }
    
    model@X <- rbind(model@X, newPointKrig)
    
    if (L=="min")
      l = min(model@y)
    else if (L=="max")
      l = max(model@y)
    else if (L=="upper95")
      l = predict.km(object = model,newdata = newPointKrig,type="UK",light.return = TRUE)$upper95
    else if (L=="lower95")
      l = predict.km(object = model,newdata = newPointKrig,type="UK",light.return = TRUE)$lower95
    else l = L
    
    model@y <- rbind(model@y, l, deparse.level=0)
    
    model@F <- trendMatrix.update(model, Xnew=data.frame(newPointKrig))
    if (model@noise.flag) {
      model@noise.var = c(model@noise.var, 0) # here is the fix!
    }
    newmodel = NULL
    try(newmodel <- computeAuxVariables(model))
    if (is.null(newmodel)) {warning("Unable to update model !");npoints=s-1;break;}
    model = newmodel
  }
  
  if (npoints==0) return()
  return(list(par = newlowDimPoints, value = model@y[(n1+1):(n1+npoints),, drop=FALSE]))
}

distXmin <- function(x,Xmin) {
  return(min(sqrt(rowSums((Xmin-matrix(x,nrow=nrow(Xmin),ncol=ncol(Xmin),byrow=TRUE))^2))))
}

# Bounds are supposed to be [0,1]
mapping_to_X <- function(y, A){
  Xmap <- A %*% y
  Xmap = pmin(Xmap, 1)
  Xmap = pmax(Xmap, -1)
  Xmap <- (Xmap + 1)/2
  return(Xmap)
}

# Warped kernel
Psi <- function(y, A, pA){
  
  Xmap <- t(A %*% y)
  Xmap <- pmin(Xmap, 1)
  Xmap <- pmax(Xmap, -1) # convex projection

  Xo <- t(pA %*% t(Xmap)) # orthogonal projection
  
  if(max(abs(Xo)) > 1){
    pivot <- Xo/max(abs(Xo))
    tmp <- dist(rbind(Xo, Xmap))[1] # distance between Xo and Xmap
    tmp2 <- sqrt(sum(pivot^2)) # norm of pivot
    Xo <- Xo * (tmp2 + tmp)/tmp2 # distortion
  }
  return(Xo) 
}

# type : type of warping, 'none' for low_dim kernel, 'convex' for high_dim kernel, 'psi' for warped kernel
warping <- function(y, A, type){
  if(type == 'none')
    return(y)
  if(type == 'convex')
    return(mapping_to_X(y, A))
  if(type== 'warped')
    return(Psi(y, A))
}



