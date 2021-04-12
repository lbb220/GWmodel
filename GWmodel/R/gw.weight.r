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
####C++ code embeded at 31/01/2021

# MAIN FUNCTION

gw.weight <-function(vdist,bw,kernel,adaptive=FALSE){
  if (is.matrix(vdist))
  {
    wgt <- gw_weight_mat(vdist,bw,kernel,adaptive)
  }
  else
  {
    wgt <- gw_weight_vec(vdist,bw,kernel,adaptive)
  }
    
  wgt
}
