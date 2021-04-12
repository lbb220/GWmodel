##Model selection,specification by a step-wise like AIC procedure in a "forward" direction
gwr.model.selection<-function(DeVar=NULL,InDeVars=NULL, data=list(),bw=NULL,approach="CV",
                     adaptive=F,kernel="bisquare",dMat=NULL,p=2, theta=0, longlat=F,
                     parallel.method=F,parallel.arg=NULL)
{
  if (is.null(DeVar) || !is.character(DeVar) || is.null(InDeVars) || !is.character(InDeVars))
    stop("Input are not correct, please recheck!")
   ##Data points
  spdf <- data
  if (!is.null(data)) {
    if (is(data, "Spatial")) {
      p4s <- proj4string(data)
      dp.locat <- coordinates(data)
      data <- as(data, "data.frame")
    } else {
      if (!is(data, "data.frame"))
        stop("Given regression data must be data.frame or Spatial*DataFrame")
    }
  }
  else stop("No regression data frame is avaiable!")
  #################################################################
  vars.df <- names(data)
  dp.n <- nrow(data)
  var.n <- length(InDeVars)
  InDeVars.Sub <- InDeVars
  model.list <- list()
  GWR.df <- c()
   
  if (missing(dMat)) {
    DM.given <- F
    dMat <- matrix(0, 0, 0)
  }
  else if (is.null(dMat) || !is.matrix(dMat)) {
    DM.given<-F
    dMat <- matrix(0, 0, 0)
  }
  else {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
      stop ("Dimensions of dMat are not correct")
  }
  #vars.idxs<-c()###Record indices used for each model
  varsindx.list <- list()
  level.vars <- c()
  tag <- 1
  adapt <- NULL
  if (parallel.method == "cluster") {
     if (missing(parallel.arg)) {
        parallel.arg.n <- max(detectCores() - 4, 2)
        parallel.arg <- makeCluster(parallel.arg.n)
      } else parallel.arg.n <- length(parallel.arg)
      clusterCall(parallel.arg, function() { library(GWmodel) })
  }
  for (i in 1:var.n) {
    AICcs <- c()
    for (j in 1:(var.n - i + 1)) {
      vars.j <- c(level.vars, InDeVars.Sub[j])
      fml <- Generate.formula(DeVar, vars.j)
	    cat("Now calbrating the model: \n", fml, "\n")
	    matL <- extract.mat(fml, data)
	    y <- matL[[1]]
	    x <- matL[[2]]
      if (is.null(bw)) {
        part1 <- paste("bw<-bw.gwr(", fml, sep = "")
        part2 <- "data=spdf,kernel=kernel,approach=approach,dMat=dMat, parallel.method=parallel.method,parallel.arg=parallel.arg)"
        expression <- paste(part1, part2, sep = ",")
        print(expression)
        eval(parse(text = expression))
      } else {
        if (adaptive) {
          stopifnot(is.numeric(bw))
          stopifnot((bw >= 0))
        } else {
          stopifnot(is.numeric(bw))
          stopifnot((bw > min(dMat)))
        }
      } 

      ##############Calibrate the GWR model
      betas <- matrix(0, nrow = dp.n, ncol = ncol(x))
      s_hat <- numeric(2)
      if (parallel.method == FALSE) {
        res <- gw_reg_all(x, y, dp.locat, FALSE, dp.locat, DM.given, dMat, TRUE, p, theta, longlat, bw, kernel, adaptive)
        betas <- res$betas
        s_hat <- res$s_hat
      } else if (parallel.method == "omp") {
        if (missing(parallel.arg)) { threads <- 0 } else {
          threads <- ifelse(is(parallel.arg, "numeric"), parallel.arg, 0)
        }
        res <- gw_reg_all_omp(x, y, dp.locat, FALSE, dp.locat, DM.given, dMat, TRUE, p, theta, longlat, bw, kernel, adaptive, threads)
        betas <- res$betas
        s_hat <- res$s_hat
      } else if (parallel.method == "cuda") {
        if (missing(parallel.arg)) { groupl <- 16 } else {
          groupl <- ifelse(is(parallel.arg, "numeric"), parallel.arg, 16)
        }
        res <- gw_reg_all_cuda(x, y, dp.locat, FALSE, dp.locat, DM.given, dMat, TRUE, p, theta, longlat, bw, kernel, adaptive, groupl)
        betas <- res$betas
        s_hat <- res$s_hat
      } else if (parallel.method == "cluster") {
        parallel.arg.results <- clusterApplyLB(parallel.arg, 1:parallel.arg.n, function(group.i, parallel.arg.n, x, y, dp.locat, DM.given, dMat, p, theta, longlat, bw, kernel, adaptive) {
          res <- gw_reg_all(x, y, dp.locat, FALSE, dp.locat, DM.given, dMat, TRUE, p, theta, longlat, bw, kernel, adaptive, parallel.arg.n, group.i)
          return(res)
        }, parallel.arg.n, x, y, dp.locat, DM.given, dMat, p, theta, longlat, bw, kernel, adaptive)
        for (i in 1:parallel.arg.n) {
          res <- parallel.arg.results[[i]]
          betas = betas + res$betas
          s_hat = s_hat + res$s_hat
        }
      } else {
        for (i in 1:dp.n)
        {
          dist.vi<-dMat[,i]
          W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
          gw.resi<-gw_reg(x,y,W.i,hatmatrix=T,i)
          betas[i,]<-as.numeric(gw.resi[[1]])
          S[i,]<-gw.resi[[2]]
          #Ci<-gw.resi[[3]]
        }
      }
      aic.rss <- as.numeric(AICc_rss1(y, x, betas, s_hat))
      model.list[[tag]] <- list(fml, vars.j)
      GWR.df <- rbind(GWR.df, c(bw, aic.rss[2], aic.rss[3], aic.rss[1]))
      AICcs <- c(AICcs, aic.rss[3])
      tag <- tag + 1
    }
    idx <- which.min(AICcs)[1]
    level.vars <- c(level.vars, InDeVars.Sub[idx])
    InDeVars.Sub <- InDeVars.Sub[-idx]
  }
  if (parallel.method == "cluster") {
    if (missing(parallel.arg)) stopCluster(parallel.arg)     
  }
  res <- list(model.list, GWR.df)
  res
}

# This version of this function is kept to make the code work with the early versions of GWmodel (before 2.0-1)
model.selection.gwr <-function(DeVar=NULL,InDeVars=NULL, data=list(),bw=NULL,approach="CV",
                     adaptive=F,kernel="bisquare",dMat=NULL,p=2, theta=0, longlat=F,
                     parallel.method=F,parallel.arg=NULL)
{
  if (is.null(DeVar) || !is.character(DeVar) || is.null(InDeVars) || !is.character(InDeVars))
    stop("Input are not correct, please recheck!")
  ##Data points
  spdf <- data
  if (!is.null(data)) {
    if (is(data, "Spatial")) {
      p4s <- proj4string(data)
      dp.locat <- coordinates(data)
      data <- as(data, "data.frame")
    } else {
      if (!is(data, "data.frame"))
        stop("Given regression data must be data.frame or Spatial*DataFrame")
    }
  }
  else stop("No regression data frame is avaiable!")
  #################################################################
  vars.df <- names(data)
  dp.n <- nrow(data)
  var.n <- length(InDeVars)
  InDeVars.Sub <- InDeVars
  model.list <- list()
  GWR.df <- c()
  
  if (missing(dMat)) {
    dMat <- matrix(0, 0, 0)
    DM.given <- F
  }
  else if (is.null(dMat) || !is.matrix(dMat)) {
    DM.given<-F
    dMat <- matrix(0, 0, 0)
  }
  else {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
      stop ("Dimensions of dMat are not correct")
  }
  #vars.idxs<-c()###Record indices used for each model
  varsindx.list <- list()
  level.vars <- c()
  tag <- 1
  adapt <- NULL
  for (i in 1:var.n) {
    AICcs <- c()
    for (j in 1:(var.n - i + 1)) {
      vars.j <- c(level.vars, InDeVars.Sub[j])
      fml <- Generate.formula(DeVar, vars.j)
      cat("Now calbrating the model: \n", fml, "\n")
      matL <- extract.mat(fml, data)
      y <- matL[[1]]
      x <- matL[[2]]
      if (is.null(bw)) {
        part1 <- paste("bw<-bw.gwr(", fml, sep = "")
        part2 <- "data=spdf,kernel=kernel,approach=approach,dMat=dMat, parallel.method=parallel.method,parallel.arg=parallel.arg)"
        expression <- paste(part1, part2, sep = ",")
        print(expression)
        eval(parse(text = expression))
      } else {
        if (adaptive) {
          stopifnot(is.numeric(bw))
          stopifnot((bw >= 0))
        } else {
          stopifnot(is.numeric(bw))
          # stopifnot((bw > min(dMat)))
        }
      }

      ##############Calibrate the GWR model
      betas <- matrix(nrow = dp.n, ncol = var.n)
      s_hat <- numeric(2)
      if (parallel.method == FALSE) {
        res <- gw_reg_all(x, y, dp.locat, FALSE, dp.locat, DM.given, dMat, TRUE, p, theta, longlat, bw, kernel, adaptive)
        betas <- res$betas
        s_hat <- res$s_hat
      } else if (parallel.method == "omp") {
        if (missing(parallel.arg)) { threads <- 0 } else {
          threads <- ifelse(is(parallel.arg, "numeric"), parallel.arg, 0)
        }
        res <- gw_reg_all_omp(x, y, dp.locat, FALSE, dp.locat, DM.given, dMat, TRUE, p, theta, longlat, bw, kernel, adaptive, threads)
        betas <- res$betas
        s_hat <- res$s_hat
      } else if (parallel.method == "cuda") {
        if (missing(parallel.arg)) { groupl <- 0 } else {
          groupl <- ifelse(is(parallel.arg, "numeric"), parallel.arg, 0)
        }
        res <- gw_reg_all_cuda(x, y, dp.locat, FALSE, dp.locat, DM.given, dMat, TRUE, p, theta, longlat, bw, kernel, adaptive, groupl)
        betas <- res$betas
        s_hat <- res$s_hat
      } else if (parallel.method == "cluster") {
        if (missing(parallel.arg)) {
          parallel.arg.n <- max(detectCores() - 4, 2)
          parallel.arg <- makeCluster(parallel.arg.n)
        } else parallel.arg.n <- length(parallel.arg)
        clusterCall(parallel.arg, function() { library(GWmodel) })
        parallel.arg.results <- clusterApplyLB(parallel.arg, 1:parallel.arg.n, function(group.i, parallel.arg.n, x, y, dp.locat, DM.given, dMat, p, theta, longlat, bw, kernel, adaptive) {
          res <- gw_reg_all(x, y, dp.locat, FALSE, dp.locat, DM.given, dMat, TRUE, p, theta, longlat, bw, kernel, adaptive, parallel.arg.n, group.i)
          return(res)
        }, parallel.arg.n, x, y, dp.locat, DM.given, dMat, p, theta, longlat, bw, kernel, adaptive)
        for (i in 1:parallel.arg.n) {
          res <- parallel.arg.results[[i]]
          betas = betas + res$betas
          s_hat = s_hat + res$s_hat
        }
        if (missing(parallel.arg)) stopCluster(parallel.arg)
      } else {
        for (i in 1:dp.n) {
          dist.vi<-dMat[,i]
          W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
          gw.resi<-gw_reg(x,y,W.i,hatmatrix=T,i)
          betas[i,]<-as.numeric(gw.resi[[1]])
          S[i,]<-gw.resi[[2]]
          #Ci<-gw.resi[[3]]
        }
      }
      aic.rss <- as.numeric(AICc_rss1(y, x, betas, s_hat))
      model.list[[tag]] <- list(fml, vars.j)
      GWR.df <- rbind(GWR.df, c(bw, aic.rss[2], aic.rss[3], aic.rss[1]))
      AICcs <- c(AICcs, aic.rss[3])
      tag <- tag + 1
    }
    idx <- which.min(AICcs)[1]
    level.vars <- c(level.vars, InDeVars.Sub[idx])
    InDeVars.Sub <- InDeVars.Sub[-idx]
  }
  res <- list(model.list, GWR.df)
  res
}
#######Extract model matrix
extract.mat <- function(formula, data = list())
{
	this.call <- match.call()
	if (!is.null(data))
    {
       if (is(data, "Spatial"))
	   {
	     data <- as(data, "data.frame")
	   }
       else
       {
         if (!is(data, "data.frame"))
           stop("Given regression data must be data.frame or Spatial*DataFrame")
       }
    }
    else stop("No regression data frame is avaiable!")
	mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
  mat.list <- list(y, x)
	mat.list
}

gwr.model.view <- function(DeVar, InDeVars, model.list)
{
  n <- length(InDeVars)
  if (n > 10)
  {
    cex <- 10 / n
  }
  else
  {
    cex <- 1
  }

  #InDeVars<-sort(InDeVars)
  numModels<-length(model.list)
  alpha<-2*pi/numModels
  cols<-rainbow(n)
  pchs<-rep(c(8,9,10,15,16,17,18,23,24),length.out=n)
  plot(x=0,y=0,xlim=c(-3*n/4, n+6),ylim=c(-n-1, n+1), cex=2, axes=F, pch=22,xlab="",ylab="",main="View of GWR model selection with different variables")
  for (i in 1:numModels)
  {
    vars<-model.list[[i]][[2]]
    nvar<-length(vars)
    p1<-c(0,0)
    for (j in 1:nvar)
    {
      radius<-sqrt(n)*sqrt(j)
      var.idx<-which(InDeVars==vars[j])
      coord<-c(radius*cos((i-1)*alpha),radius*sin((i-1)*alpha))
      lines(x=c(p1[1], coord[1]),y=c(p1[2], coord[2]), col="grey",lwd=cex)
      points(x=coord[1], y=coord[2], col=cols[var.idx],pch=pchs[var.idx],cex=(cex*i/numModels+0.3))
      p1<-coord
    }
    text(x=(radius+0.5)*cos((i-1)*alpha),y=(radius+0.5)*sin((i-1)*alpha), as.character(i), cex=cex*0.6)
  }
  legend(x=n+2, y=n/2, col=c("black",cols),pch=c(22,pchs), c(DeVar, InDeVars),box.col="white")
}

# This version of this function is kept to make the code work with the early versions of GWmodel (before 2.0-1)
model.view.gwr <-function(DeVar, InDeVars, model.list)
{
  n<-length(InDeVars)
  if (n>10)
  {
    cex<-10/n
  }
  else
  {
    cex<-1
  }

  #InDeVars<-sort(InDeVars)
  numModels<-length(model.list)
  alpha<-2*pi/numModels
  cols<-rainbow(n)
  pchs<-rep(c(8,9,10,15,16,17,18,23,24),length.out=n)
  plot(x=0,y=0,xlim=c(-3*n/4, n+6),ylim=c(-n-1, n+1), cex=2, axes=F, pch=22,xlab="",ylab="",main="View of GWR model selection with different variables")
  for (i in 1:numModels)
  {
    vars<-model.list[[i]][[2]]
    nvar<-length(vars)
    p1<-c(0,0)
    for (j in 1:nvar)
    {
      radius<-sqrt(n)*sqrt(j)
      var.idx<-which(InDeVars==vars[j])
      coord<-c(radius*cos((i-1)*alpha),radius*sin((i-1)*alpha))
      lines(x=c(p1[1], coord[1]),y=c(p1[2], coord[2]), col="grey",lwd=cex)
      points(x=coord[1], y=coord[2], col=cols[var.idx],pch=pchs[var.idx],cex=(cex*i/numModels+0.3))
      p1<-coord
    }
    text(x=(radius+0.5)*cos((i-1)*alpha),y=(radius+0.5)*sin((i-1)*alpha), as.character(i), cex=cex*0.6)
  }
  legend(x=n+2, y=n/2, col=c("black",cols),pch=c(22,pchs), c(DeVar, InDeVars),box.col="white")
}
gwr.model.sort<-function(Sorting.list , numVars, ruler.vector)
{
  n<-length(Sorting.list)
  numMoldes<-length(ruler.vector)
  indxs<-c()
  tag<-0
  for (i in numVars:1)
  {
    tmpV<-ruler.vector[(tag+1):(tag+i)]
    indx<-sort(tmpV, decreasing=T, index=T)$ix
    indxs<-c(indxs, indx+tag)
    tag<-tag+i
  }
  res<-list()
  for (i in 1:n)
    {
      list.i<-Sorting.list[[i]]
      if (is.list(list.i))
      {
         tmp.L<-list()
         for (j in 1:numMoldes) tmp.L[[j]]<-list.i[[indxs[j]]]
         res[[i]]<-tmp.L
      }
      else
      {
        tmp.V<-c()
        for (j in 1:numMoldes) tmp.V<-rbind(tmp.V,list.i[indxs[j],])
        res[[i]]<-tmp.V
      }
    }
    res
}

# This version of this function is kept to make the code work with the early versions of GWmodel (before 2.0-1)
model.sort.gwr<-function(Sorting.list , numVars, ruler.vector)
{
  n<-length(Sorting.list)
  numMoldes<-length(ruler.vector)
  indxs<-c()
  tag<-0
  for (i in numVars:1)
  {
    tmpV<-ruler.vector[(tag+1):(tag+i)]
    indx<-sort(tmpV, decreasing=T, index=T)$ix
    indxs<-c(indxs, indx+tag)
    tag<-tag+i
  }
  res<-list()
  for (i in 1:n)
    {
      list.i<-Sorting.list[[i]]
      if (is.list(list.i))
      {
         tmp.L<-list()
         for (j in 1:numMoldes) tmp.L[[j]]<-list.i[[indxs[j]]]
         res[[i]]<-tmp.L
      }
      else
      {
        tmp.V<-c()
        for (j in 1:numMoldes) tmp.V<-rbind(tmp.V,list.i[indxs[j],])
        res[[i]]<-tmp.V
      }
    }
    res
}
####Generate formula based on the given depedent and indepednent variables
Generate.formula<-function(DeVar,InDeVars)
{
  fml<-paste(paste(DeVar, "~", sep=""), InDeVars[1],sep="")
  var.n<-length(InDeVars)
  if (var.n>1)
  {
    for (i in 2:var.n)
      fml<-paste(fml, InDeVars[i], sep="+")
  }
  fml
} 
