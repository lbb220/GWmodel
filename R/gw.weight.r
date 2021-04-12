################## NAME: gw.weight ###################
## ARGUMENTS IN:
#### vdist: numeric vector or matrix of distances (from gw.dist.r)
#### bw: scalar, bandwidth or number of nearest neighbours
#### kernel: text vector of function chosen
########## boxcar: wgt=1 if dist < bw, wgt=0 otherwise  
########## gaussian: wgt = exp(-.5*(vdist/bw)^2)
########## bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise   
########## tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise    
########## adaptive: if TRUE calulate the adaptive kernel, and bw correspond to the number of nearest neighbour
####
## ARGUMENTS OUT:
#### if vdist is a vector the output is a vector of weights of the length of vdist
#### if vdist is a matrix the output is a matrix of weights of the dim of vdist, containing in column i the weights for observation i
## REFERENCES: Book pg 56-57
########################

####### Fixed Kernel

# exponential kernel

gw.weight.exponential <- function(vdist,bw){
  if (is.matrix(vdist))
  {
    m <- ncol(vdist)
    bw <- rep(bw, m)
    wt <- exp_wt_mat(vdist, bw)
  }
  else
    wt <- as.numeric(exp_wt_vec(vdist, bw))
  wt
}


# Boxcar kernel

gw.weight.box  <- function(vdist,bw){
  {vdist<= bw}*1
}

# Gaussian kernel

gw.weight.gau  <- function(vdist,bw){
  if (is.matrix(vdist))
  {
    m <- ncol(vdist)
    bw <- rep(bw, m)
    wt <- gauss_wt_mat(vdist, bw)
  }
  else
    wt <- as.numeric(gauss_wt_vec(vdist, bw))
  wt
}

# Fixed Bisquare kernel

gw.weight.bis <- function(vdist,bw){
  if (is.matrix(vdist))
  {
    m <- ncol(vdist)
    bw <- rep(bw, m)
    wt <- bisq_wt_mat(vdist, bw)
  }
  else
    wt <- as.numeric(bisq_wt_vec(vdist, bw))
  wt
}

# Fixed Tricube kernel

gw.weight.tri <- function(vdist,bw){
  if (is.matrix(vdist))
  {
    m <- ncol(vdist)
    bw <- rep(bw, m)
    wt <- tri_wt_mat(vdist, bw)
  }
  else
    wt <- as.numeric(tri_wt_vec(vdist, bw))
  wt
}

####### Adaptive Kernel

# Adaptive Bisquare kernel 
## bw correspond to the number of nearest neighbours

#####

gw.weight.gau.ad <- function(vdist,bw){
  
  if (is.matrix(vdist)){
    nr <- nrow(vdist)
    dn <- bw/nr
    if(dn<=1)
    {
      rnk<-apply(vdist,2,rank,ties.method='first')
      bw<- vdist[rnk==bw]                  # bandwidth is at bw-th distance	
    }
    else
    {
      bw <- dn*apply(vdist,2,max)
    }
    if(length(bw)>0)
      wgt<- exp(vdist^2 / (-2 * bw^2))
    else
      wgt <- diag(1, dim(vdist)[1], dim(vdist)[2])
  }else{
    nr <- length(vdist)
    dn <- bw/nr
    if(dn<=1)
    {
      rnk<-rank(vdist,ties.method='first') # ranked distances
      cond<- which(rnk <= bw) 
      bw<- vdist[rnk==bw]                  # bandwidth is at bw-th distance
    }
    else
    {
      bw <- dn*max(vdist)
    }
    if(bw>0)
      wgt<- exp(vdist^2 / (-2 * bw^2))
    else
    {
      wgt <- numeric(length(vdist))
      wgt[cond] <- 1
    }  
  }
  wgt
}

#####

gw.weight.exp.ad <- function(vdist,bw){
  if (is.matrix(vdist)){
    nr <- nrow(vdist)
    dn <- bw/nr
    if(dn<=1)
    {
      rnk<-apply(vdist,2,rank,ties.method='first')
      bw<- vdist[rnk==bw]                  # bandwidth is at bw-th distance	
    }
    else
    {
      bw <- dn*apply(vdist,2,max)
    }
    if(length(bw)>0)
      wgt<-exp_wt_mat(vdist, bw)
    else
      wgt <- diag(1, dim(vdist)[1], dim(vdist)[2])
  }
  else
  {
    nr <- length(vdist)
    dn <- bw/nr
    if(dn<=1)
    {
      rnk<-rank(vdist,ties.method='first') # ranked distances
      cond<- which(rnk <= bw)
      bw<- vdist[rnk==bw]                  # bandwidth is at bw-th distance
    }
    else
    {
      bw <- dn*max(vdist)
    } 
    
    if(length(bw)>0)
      wgt<- as.numeric(exp_wt_vec(vdist, bw))
    else
    {
      wgt <- numeric(length(vdist))
      wgt[cond] <- 1
    }
  }
  wgt
}

#####


gw.weight.bis.ad <- function(vdist,bw){
  
  if (is.matrix(vdist))
  {
    nr <- nrow(vdist)
    dn <- bw/nr
    if(dn<=1)
    { 
      rnk<-apply(vdist,2,rank,ties.method='first')
      cond<- rnk <=bw  
      bw<- vdist[rnk == bw]
    }
    else
    {
      bw <- dn*apply(vdist,2,max)
    }	
    wgt<- matrix(0,nrow(vdist),ncol(vdist))
    if(length(bw)>0)
    {
      wgt<- bisq_wt_mat(vdist, bw)
    }
    else
      diag(wgt) <- 1
    
  }
  else
  {
    nr <- length(vdist)
    dn <- bw/nr
    if(dn<=1)
    {
      rnk<-rank(vdist,ties.method='first') # ranked distances
      cond<- which(rnk <= bw)               # condition for locations less than bw-th
      bw<- vdist[rnk==bw]                  # bandwidth is at bw-th distance
    }
    else
    {
      bw <- dn*max(vdist)
    } 
    wgt<-numeric(length(vdist))
    if(!is.na(bw) >0)
    {
      wgt<- as.numeric(bisq_wt_vec(vdist, bw))
    }
    else
    {
      wgt[cond]<- 1
    }
  }
  wgt
  
}


#####

gw.weight.tri.ad <- function(vdist,bw){
  
  if (is.matrix(vdist)){
    nr <- nrow(vdist)
    dn <- bw/nr
    if(dn<=1)
    { 
      rnk<-apply(vdist,2,rank,ties.method='first')
      cond<- rnk<= bw 
      bw<- vdist[rnk == bw]
    }
    else
    {
      bw <- dn*apply(vdist,2,max)
    }
    wgt<- matrix(0,nrow(vdist),ncol(vdist))
    if(length(bw)>0)
    {
      wgt <- tri_wt_mat(vdist, bw)
    }
    else
      diag(wgt) <- 1
    
  }else{
    nr <- length(vdist)
    dn <- bw/nr
    if(dn<=1)
    {
      rnk<-rank(vdist,ties.method='first') # ranked distances
      cond<- which(rnk <= bw)               # condition for locations less than bw-th
      bw<- vdist[rnk==bw]                  # bandwidth is at bw-th distance	
    }
    else
    {
      bw <- dn*max(vdist)
    }
    wgt<-numeric(length(vdist))
    if(!is.na(bw))
    {
      wgt <- as.numeric(tri_wt_vec(vdist, bw))
    }
    else
      wgt[cond]<-1
  }
  wgt
}


#####

gw.weight.box.ad <- function(vdist,bw){
  
  if (is.matrix(vdist)) 
  {
    nr <- nrow(vdist)
    rnk<-apply(vdist,2,rank,ties.method='first')
    if(bw>nr)
      bw <- nr
  }
  else 
  {
    nr <- length(vdist)
    rnk<-rank(vdist,ties.method='first') # ranked distances
    if(bw>nr)
      bw <- nr
  }
  {rnk <= bw}*1
}


# MAIN FUNCTION

gw.weight <-function(vdist,bw,kernel,adaptive=FALSE){
  
  if(adaptive==FALSE) switch(kernel,
                             gaussian = gw.weight.gau (vdist,bw),
                             bisquare = gw.weight.bis (vdist,bw),
                             tricube  = gw.weight.tri (vdist,bw),
                             boxcar   = gw.weight.box (vdist,bw),
                             exponential = gw.weight.exponential (vdist,bw))
  else switch(kernel,
              gaussian = gw.weight.gau.ad (vdist,bw),
              bisquare = gw.weight.bis.ad (vdist,bw),
              tricube  = gw.weight.tri.ad (vdist,bw),
              boxcar   = gw.weight.box.ad (vdist,bw),
              exponential = gw.weight.exp.ad (vdist,bw))
}


######### End of the Code
