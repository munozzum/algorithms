#help=Morris's Elementary Effects Screening Method
#type=Sensitivity analysis
#output=Sensitivities
#options=r=10,levels=4
#require=sensitivity

#' constructor and initializer of R session
Morris <- function(options) {
  library(sensitivity)

  # all parameters are initialy strings, so you have to put as global non-string values
  options$r <- as.integer(options$r)
  options$levels <- as.integer(options$levels)

  algorithm = new.env()
  lapply(names(options),function(x) assign(x,options[[x]],algorithm))
  return(algorithm)
}

#' first design building. All variables are set in [0,1].
#' @param d number of variables
#' @return next design of experiments
getInitialDesign <- function(algorithm,d) {
  set.seed(1)
  algorithm$m <- morris(model=NULL,factors=d,r=algorithm$r,design=list(type="oat",levels=algorithm$levels))
  return(algorithm$m$X)
}

#' iterated design building.
#' @param X data frame of current doe variables (in [0,1])
#' @param Y data frame of current results
#' @return  next doe step
getNextDesign <- function(algorithm,X,Y) {
  return()
}

#' final analysis.
#' @param X data frame of doe variables (in [0,1])
#' @param Y data frame of results
#' @return HTML string of analysis
displayResults <- function(algorithm,X,Y) {
  eval(expr=parse(text="tell(m,Yi)"),envir=algorithm) #tell(algorithm$m,Y)
  png(file="plot.png",bg="transparent",height=600,width=600)
  plot(algorithm$m)
  dev.off()
  
  mu=colMeans(algorithm$m$ee)
  mu.star=colMeans(abs(algorithm$m$ee))
  sig=sqrt((colSums(algorithm$m$ee^2)-mu^2*4)/(algorithm$r-1))
  effect=rep("Weak effect",length(mu))
  for (i in 1:length(mu)) if (mu.star[i]>0.2*max(mu.star)) if (sig[i]>0.2*max(sig)) effect[i]="Non linear or interactive effect" else effect[i]="Linear effect"
  
  html=paste0("<HTML name='sensitivity'>",paste0(collapse="\n",capture.output(print(xtable::xtable(  cbind(mu,mu.star,sig,effect)),type="html"))),"<br/><img src='plot.png'/></HTML>")
  return(html)
}

