cv4postpr_nnet = function (index, sumstat, postpr.out = NULL, nval, tols, method, 
          subset = NULL, kernel = "epanechnikov", numnet = 10, sizenet = 5, 
          lambda = c(1e-04, 0.001, 0.01), trace = FALSE, maxit = 500, 
          ...) 
{
  linout <- TRUE
  if (missing(nval)) 
    stop("'nval' must be supplied.", call. = F)
  if (is.null(postpr.out) && missing(method)) 
    stop("Method must be supplied when 'postpr.out' is NULL.", 
         call. = F)
  if (length(index) != na.omit(length(index))) 
    stop("'index' contains missing values. Models must be specified for each simulation.", 
         call. = F)
  if (!prod(table(index) > nval)) 
    stop("'nval' has to be smaller or equal to number of simulations for any of the models. Choose a smaller 'nval'.", 
         call. = F)
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (!is.null(postpr.out)) {
    subset <- postpr.out$na.action
    method <- postpr.out$method
    kernel <- "epanechnikov"
  }
  if (is.vector(sumstat)) {
    numstat <- 1
    sumstat <- matrix(sumstat, ncol = 1)
  }
  else numstat <- dim(sumstat)[2]
  numsim <- length(index)
  if (!is.null(postpr.out)) {
    if (numstat != postpr.out$numstat || numsim != length(postpr.out$na.action)) {
      stop("The number of summary statistics, or simulations provided in 'sumstat' are not the same as in 'postpr.out'.", 
           call. = F)
    }
    else if (!prod(unique(index) %in% postpr.out$names$models)) {
      stop("Models in 'index' are not the same as in 'postpr.out', or different names are used.", 
           call. = F)
    }
    else if (!prod(colnames(sumstat) %in% postpr.out$names$statistics.names)) {
      stop("Summary statistics in 'sumstat' are not the same as in 'postpr.out', or different names are used.", 
           call. = F)
    }
    else {
      mymodels <- postpr.out$names$models
      statnames <- postpr.out$names$statistics.names
    }
  }
  else {
    mymodels <- levels(factor(index))
    if (length(colnames(sumstat))) {
      statnames <- colnames(sumstat)
    }
    else {
      warning("No statistics names are given, using S1, S2, ...", 
              call. = F, immediate = T)
      statnames <- paste("S", 1:numstat, sep = "")
    }
  }
  gwt <- rep(TRUE, length(sumstat[, 1]))
  gwt[attributes(na.omit(sumstat))$na.action] <- FALSE
  if (is.null(subset)) 
    subset <- rep(TRUE, length(sumstat[, 1]))
  gwt <- as.logical(gwt * subset)
  cvsamp <- unlist(tapply(c(1:length(index))[gwt], index[gwt], 
                          sample, nval))
  tols <- sort(tols)
  allprobs <- list()
  mycall <- list()
  for (mytol in tols) {
    res <- matrix(ncol = length(unique(index)), nrow = length(cvsamp))
    for (i in 1:length(cvsamp)) {
      mysamp <- cvsamp[i]
      mytrue <- index[mysamp]
      mytarget <- sumstat[mysamp, ]
      myindex <- index[-mysamp]
      mysumstat <- sumstat[-mysamp, ]
      mysubset <- subset[-mysamp]
      # 'abc' R package version 2.1
      # Error in 'cv4postpr' due to 'too much weights' i.e. too much variables/hidden layers (nnet is limited to 1000 weights by default)
      # -> cv4postpr() not passing supplementary arguments to nnet at contrary to postpr()
      # at line 114 in 'cv4postpr.R' (referred as version 1.7 of CRAN/ABC) https://github.com/cran/abc/blob/master/R/cv4postpr.R
      # Fixed by adding ', ...' at the end of arguments passed to postpr()
      # in order to successfully pass supplementary arguments (MaxNWts = 100000 in my case) up to nnet()
      subres <- postpr(target = mytarget, index = myindex, 
                       sumstat = mysumstat, tol = mytol, subset = mysubset, 
                       method = method, kernel = kernel, ...)
      if (subres$method == "rejection") 
        res[i, ] <- summary.postpr(subres, print = F, 
                                   ...)$Prob
      if (subres$method == "mnlogistic") 
        res[i, ] <- summary.postpr(subres, print = F, 
                                   ...)$mnlogistic$Prob
      if (subres$method == "neuralnet") 
        res[i, ] <- summary.postpr(subres, print = F, 
                                   ...)$neuralnet$Prob
    }
    colnames(res) <- mymodels
    rownames(res) <- index[cvsamp]
    allprobs[[paste("tol", mytol, sep = "")]] <- res
    allnames <- lapply(allprobs, apply, 1, function(xx) mymodels[which(xx == 
                                                                         max(xx))])
    mycall[[paste("tol", mytol, sep = "")]] <- call("postpr", 
                                                    target = quote(target), index = quote(index), sumstat = quote(sumstat), 
                                                    tol = mytol, subset = quote(subset), method = subres$method, 
                                                    kernel = subres$kernel)
  }
  cv4postpr.out <- list(calls = mycall, cvsamples = cvsamp, 
                        tols = tols, true = index[cvsamp], estim = allnames, 
                        model.probs = allprobs, method = method, names = list(models = mymodels, 
                                                                              statistics.names = statnames), seed = seed)
  class(cv4postpr.out) <- "cv4postpr"
  invisible(cv4postpr.out)
}

is.cv4postpr <- function(x){
  if (inherits(x, "cv4postpr")) TRUE
  else FALSE
}

summary.cv4postpr <- function(object, probs=TRUE, print = TRUE, digits = max(3, getOption("digits")-3), ...){
  
  if (!inherits(object, "cv4postpr")) 
    stop("Use only with objects of class \"cv4postpr\".", call.=F)
  
  cv4postpr.out <- object
  tols <- cv4postpr.out$tols
  numtols <- length(tols)
  true <- cv4postpr.out$true
  estim <- cv4postpr.out$estim
  method <- cv4postpr.out$method
  model.probs <- cv4postpr.out$model.probs
  nmodels <- length(cv4postpr.out$names$models)
  nval <- length(true)/nmodels
  
  if(print) cat("Confusion matrix based on ", nval, " samples for each model.\n\n", sep="")
  cm <- lapply(estim, function(x) table(true, x))
  cm <- lapply(cm, function(x) {attributes(dimnames(x))$names <- NULL; x})
  if(print) print(cm); cat("\n")
  
  if(probs){
    if(print) cat(paste("Mean model posterior probabilities (", method ,")\n\n", sep="")) 
    myprs <- lapply(model.probs, apply, 2, tapply, true, mean)
    if(print) print(lapply(myprs, round, digits=digits))
    out <- list(conf.matrix=cm, probs=myprs)
  }
  else out <- cm
  
  invisible(out)
}

plot.cv4postpr <- function(x, probs=FALSE, file = NULL, postscript = FALSE, onefile = TRUE, ask = !is.null(deviceIsInteractive()), caption = NULL, ...){
  
  if (!inherits(x, "cv4postpr")) 
    stop("Use only with objects of class \"cv4postpr\".", call.=F)
  
  cv4postpr.out <- x
  tols <- cv4postpr.out$tols
  numtols <- length(tols)
  true <- cv4postpr.out$true
  estim <- cv4postpr.out$estim
  method <- cv4postpr.out$method
  model.probs <- cv4postpr.out$model.probs
  nmodels <- length(cv4postpr.out$names$models)
  nval <- length(true)/nmodels
  
  if(is.null(caption)) caption <- "Confusion matrix"
  
  ## Devices
  save.devAskNewPage <- devAskNewPage()
  if(!is.null(file)){
    file <- substitute(file)
    if(!postscript) pdf(file = paste(file, "pdf", sep="."), onefile=onefile)
    if(postscript) postscript(file = paste(file, "ps", sep="."), onefile=onefile)
  }
  else{
    if (ask && 1 < numtols) {
      devAskNewPage(TRUE)
    }
  }
  
  par(cex = 1, cex.main = 1.2, cex.lab = 1.1)
  for(i in 1:numtols){
    if(probs){
      mym <- lapply(model.probs, apply, 2, tapply, true, mean)
      barplot(t(mym[[paste("tol", tols[i], sep="")]]), ylab="Mean model probability", ...)
    }
    else{
      mym <- lapply(estim, table, true)
      barplot(mym[[paste("tol", tols[i], sep="")]], ylab="Frequency", ...)
    }
    title(caption, sub=paste("Tolerance rate = ", tols[i], sep=""))
  }
  
  if(!is.null(file)){
    dev.off()
  }
  else devAskNewPage(save.devAskNewPage)
  invisible()
  
}


######################################################################
#
# postpr.R
#
# copyright (c) 2011-05-30, Katalin Csillery, Olivier Francois and
# Michael GB Blum with some initial code from Mark Beaumont
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
# 
# Part of the R/abc package
# Contains: postpr, summary.postpr
#
######################################################################

postpr <- function(target, index, sumstat, tol, subset=NULL, method, corr=TRUE, kernel="epanechnikov",
                   numnet = 10, sizenet = 5, lambda = c(0.0001,0.001,0.01), trace = TRUE, maxit = 500, ...){
  
  linout <- FALSE
  call <- match.call()
  
  ## general checks that the function is used correctly
  ## ###################################################
  
  if(missing(target)) stop("'target' is missing with no default", call.=F)
  if(missing(index)) stop("'index' is missing with no default", call.=F)
  if(missing(sumstat)) stop("'sumstat' is missing with no default", call.=F)
  if(!is.vector(index)) stop("'index' has to be a vector.", call.=F)
  if(!is.matrix(sumstat) && !is.data.frame(sumstat) && !is.vector(sumstat)) stop("'sumstat' has to be a matrix, data.frame or vector.", call.=F)
  if(missing(tol)) stop("'tol' is missing with no default", call.=F)
  if(missing(method)) stop("'method' is missing with no default", call.=F)
  if(!any(method == c("rejection", "mnlogistic", "neuralnet")))
    stop("Method must be 'rejection', 'mnlogistic' or 'neuralnet'.", call.=F)
  if(length(unique(index)) == 1)
    stop("At least two different models must be given.", call.=F)
  if(method == "rejection") rejmethod <- TRUE
  else rejmethod <- FALSE
  
  if(is.data.frame(sumstat)) sumstat <- as.matrix(sumstat)
  if(is.vector(sumstat)) sumstat <- matrix(sumstat, ncol=1)
  if(is.list(target)) target <- unlist(target)
  if(length(target)!=dim(sumstat)[2]) stop("Number of summary statistics in 'target' has to be the same as in 'sumstat'.", call.=F)
  if(length(index) != length(sumstat[,1]))
    stop("'index' must be the same length as the number of rows in 'sumstat'.", call.=F)
  
  if(!any(kernel == c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine"))){
    kernel <- "epanechnikov"
    warning("Kernel is incorrectly defined. Setting to default kernel (Epanechnikov)", call.=F, immediate=T)
  }
  
  if(is.vector(sumstat)) sumstat <- matrix(sumstat, ncol=1)
  if(length(target)!=dim(sumstat)[2]) stop("Number of summary statistics in 'target' has to be the same as in 'sumstat'.", call.=F)
  
  index <- factor(index)
  mymodels <- levels(index)
  if(!is.numeric(target)) target <- as.numeric(target)
  
  ## parameter and/or sumstat values that are to be excluded
  ## #######################################################
  gwt <- rep(TRUE,length(sumstat[,1]))
  gwt[attributes(na.omit(sumstat))$na.action] <- FALSE
  if(is.null(subset)) subset <- rep(TRUE,length(sumstat[,1]))
  gwt <- as.logical(gwt*subset)
  
  sumstat <- as.data.frame(sumstat) ## ?????????????????
  ## extract names of statistics if given
  ## ####################################
  nss <- length(sumstat[1,])
  if(!length(colnames(sumstat))){
    warning("No summary statistics names are given, using S1, S2, ...", call.=F, immediate=T)
    statnames <- paste("S", 1:nss, sep="")
  }
  else statnames <- colnames(sumstat)
  
  ## stop if zero var in sumstat
  ## ###########################
  cond1 <- as.logical(apply(sumstat, 2, function(x) length(unique(x))-1))
  if(!all(cond1)) stop("Summary statistic(s) have zero variance.", call.=F)
  if(!any(cond1)){
    warning("Statistic(s) ", statnames[!cond1], " have zero variance. Excluding from estimation....", sep="\t", call.=F, immediate=T)
    sumstat <- sumstat[,cond1]
    nss <- length(sumstat[1,])
    statnames <- colnames(sumstat)
    target <- target[cond1]
  }
  
  ## scale everything
  ## ################
  scaled.sumstat <- sumstat
  for(j in 1:nss){
    scaled.sumstat[,j] <- normalise(sumstat[,j],sumstat[,j][gwt])
  }
  
  for(j in 1:nss){
    target[j] <- normalise(target[j],sumstat[,j][gwt])
  }
  
  ## calculate euclidean distance
  ## ############################
  sum1 <- 0
  for(j in 1:nss){
    sum1 <- sum1 + (scaled.sumstat[,j]-target[j])^2
  }
  dist <- sqrt(sum1)
  
  # includes the effect of gwt in the tolerance
  dist[!gwt] <- floor(max(dist[gwt])+10)
  
  # wt1 defines the region we're interested in
  abstol <- quantile(dist,tol)
  if(kernel == "gaussian") wt1 <- rep(TRUE, length(dist)) ## ???????????????????????????
  else{
    ceiling(length(dist)*tol)->nacc
    sort(dist)[nacc]->ds
    wt1 <- (dist <= ds)
    aux<-cumsum(wt1)
    wt1 <- wt1 & (aux<=nacc)
  }
  
  ## select summary statistics in region
  ## ##################################
  ss <- scaled.sumstat[wt1,]
  values <- index[wt1]
  pred <- table(values)/length(values)
  
  statvar <- as.logical(apply(scaled.sumstat[wt1, , drop=FALSE], 2, function(x) length(unique(x))-1))
  cond2 <- !any(statvar)
  
  if(cond2 && !rejmethod)
    stop("Zero variance in the summary statistics in the selected region.\nTry: checking summary statistics, choosing larger tolerance, or rejection method.", call.=F)
  
  
  ## if simple rejection or in the selected region there is no var in sumstat
  ## #########################################################################
  if(rejmethod){
    if(cond2) warning("Zero variance in the summary statistics in the selected region. Check summary statistics, consider larger tolerance.", call.=F, immediate=T)
    weights <- NULL
    pred.logit <- NULL
  }
  
  ## regression correction
  ## ######################
  else{
    if(cond2) cat("Warning messages:\nStatistic(s)",
                  statnames[!statvar],
                  "has/have zero variance in the selected region.\nConsider using larger tolerance or the rejection method or discard this/these statistics.\n", sep="\t")
    
    ## weights
    if(kernel == "epanechnikov") weights <- 1 - (dist[wt1]/abstol)^2
    if(kernel == "rectangular") weights <- dist[wt1]/abstol
    if(kernel == "gaussian") weights <- 1/sqrt(2*pi)*exp(-0.5*(dist/abstol)^2)
    if(kernel == "triangular") weights <- 1 - abs(dist[wt1]/abstol)
    if(kernel == "biweight") weights <- (1 - (dist[wt1]/abstol)^2)^2
    if(kernel == "cosine") weights <- cos(pi/2*dist[wt1]/abstol)
    
    ok <- index[wt1] # models accepted
    fml <- as.formula(paste("ok ~ ", paste(statnames, collapse= "+")))
    
    if(length(unique(ok)) < length(mymodels)) {
      warning(paste("There are",length(mymodels),"models but only",length(unique(ok)), "for which simulations have been accepted.\nNo regression is performed, method is set to rejection.\nConsider increasing the tolerance rate."),sep="", call.=F, immediate=T)
      weights <- NULL
      pred.logit <- NULL
      method <- "rejection"
    }
    
    ## calculating the number of weights for multinom
    mymnw <- (nss+2) * length(mymodels)
    
    if(method == "mnlogistic"){
      ss<-data.frame(ss)
      colnames(ss)<-statnames
      fit1 <- multinom(fml, data = ss, weigths = weights, trace=F, MaxNWts = mymnw + 1, ...)
      target <- as.data.frame(matrix(target, nrow=1))
      names(target) <- statnames
      pred <- predict(fit1, target, type="probs")
      if(length(pred) == 1){
        pred <- c(1-pred,pred)
        names(pred) <- levels(ok)
      }
    }
    
    else if(method == "neuralnet"){
      ss<-data.frame(ss)
      colnames(ss)<-statnames
      lambda <- sample(lambda, numnet, replace=T)
      target <- as.data.frame(matrix(target, nrow=1))
      names(target) <- statnames
      pred <- 0
      for(i in 1:numnet){
        fit1 <- nnet(fml, data=ss, weights = weights, decay = lambda[i],
                     size = sizenet, trace = trace, linout = linout, maxit = maxit, ...)
        if(length(mymodels)==2)
        {
          auxm<-predict(fit1, target, type="raw")
          pred <- pred + c(1-auxm,auxm)
        }
        else
          pred <- pred + predict(fit1, target, type="raw")
      }
      pred <- pred/numnet
      if(length(mymodels)!=2)
      {
        temp <- rep(0, length(mymodels))
        names(temp) <- mymodels
        temp[match(colnames(pred), mymodels)] <- pred
        pred <- temp
      }
      else names(pred) <- levels(ok)
    }
    ## correction for potentially different numbers of simulations per models
    ratio <- (pred*length(index)*tol) / table(index)
    pred <- ratio/sum(ratio)
    attributes(dimnames(pred)) <- NULL
  }
  
  if(rejmethod){
    postpr.out <- list(values=values, ss=ss, call=call, na.action=gwt, method=method, corr=corr, nmodels=table(index),
                       numstat=nss, names=list(models=mymodels, statistics.names=statnames))
  }
  else{
    postpr.out <- list(values=values, pred=pred, ss=ss, weights=weights, call=call, na.action=gwt, method=method, corr=corr, nmodels=c(table(index)),
                       numstat=nss, names=list(models=mymodels, statistics.names=statnames))
  }
  class(postpr.out) <- "postpr"
  invisible(postpr.out)
}

summary.postpr <- function(object, rejection = TRUE, print = TRUE, digits = max(3, getOption("digits")-3), ...){
  
  if (!inherits(object, "postpr")) 
    stop("Use only with objects of class \"postpr\".", call.=F)
  
  postpr.out <- object
  cl <- postpr.out$call
  npost <- length(postpr.out$values)
  pred <- postpr.out$pred
  allvals <- postpr.out$values
  postmod <- levels(postpr.out$values)
  nmod <- length(postmod)
  method <- postpr.out$method
  corr <- postpr.out$corr
  nmodels <- postpr.out$nmodels
  
  if(print){
    cat("Call: \n")
    dput(cl, control=NULL)  
    cat(paste("Data:\n postpr.out$values (",npost," posterior samples)\n", sep=""))
    cat(paste("Models a priori:\n "))
    cat(postpr.out$names$models, sep=", ")
    cat(paste("\nModels a posteriori:\n "))
    cat(postmod, sep=", ")
    if(corr & length(unique(nmodels))>1){
      cat("\n")
      warning("Posterior model probabilities are corrected for unequal number of simulations per models.", immediate.=T, call.=F)
    }
    cat("\n\n")
  }
  
  if(rejection || method == "rejection"){
    
    if(print) cat("Proportion of accepted simulations (rejection):\n")
    allpr <- table(allvals)/length(allvals)
    if(corr){
      ratio <- (allpr*npost) / nmodels
      allpr <- ratio/sum(ratio)
    }
    prnames <- dimnames(allpr)$allvals
    allpr <- c(allpr); names(allpr) <- prnames
    if(print) print(round(allpr, digits=digits))
    
    if(nmod>1){
      pr.rej <- table(allvals)/length(allvals)
      bf.rej <- t(matrix(pr.rej, nmod, nmod, byrow=T)/matrix(pr.rej, nmod, nmod, byrow=F))
      colnames(bf.rej) <- postmod
      rownames(bf.rej) <- postmod
      bf.rej <- as.table(bf.rej)
      if(print){
        cat("\nBayes factors:\n")
        print(round(bf.rej, digits=digits))
        cat("\n\n")
      }
    }
    else bf.rej <- NA
    
  }
  
  if(method == "mnlogistic" | method == "neuralnet"){
    
    if(print){
      cat(paste("Posterior model probabilities (", method, "):\n", sep=""))
      print(round(pred, digits=digits))
    }
    if(nmod>1){
      bf.reg <- t(matrix(pred[pred!=0], nmod, nmod, byrow=T)/matrix(pred[pred!=0], nmod, nmod, byrow=F))
      colnames(bf.reg) <- postmod
      rownames(bf.reg) <- postmod
      bf.reg <- as.table(bf.reg)
      if(print){
        cat("\nBayes factors:\n")
        print(round(bf.reg, digits=digits))
        cat("\n")
      }
    }
    else bf.reg <- NA
    
    if(rejection){
      if(method == "mnlogistic")
        out <- list(rejection=list(Prob=allpr, BayesF=bf.rej), mnlogistic=list(Prob=pred, BayesF=bf.reg))
      if(method == "neuralnet")
        out <- list(rejection=list(Prob=allpr, BayesF=bf.rej), neuralnet=list(Prob=pred, BayesF=bf.reg))
    }
    else{
      out <- list(Prob=pred, BayesF=bf.reg)
    }
  }
  else out <- list(Prob=allpr, BayesF=bf.rej)
  invisible(out)
}

######################################################################
#
# misc.R
#
# copyright (c) 2011-05-30, Katalin Csillery, Olivier Francois and
# Michael GB Blum with some initial code from Mark Beaumont
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
# 
# Part of the R/abc package
# Contains: normalise,namesWarningFilter
#
######################################################################


normalise <- function(x,y){
  if(mad(y) == 0)
    return (x)
  else
    return (x/mad(y))
}


namesWarningFilter <- function(x){
  if( any( grepl( "No parameter names are given", x) ) ) invokeRestart( "muffleWarning" )
  if( any( grepl( "No summary statistics names are given", x) ) ) invokeRestart( "muffleWarning" )
}