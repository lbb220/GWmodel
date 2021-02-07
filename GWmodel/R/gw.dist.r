################################################################################
# gw.dist: caculate a distance vector or distance matrix between (a) regression point and data points
# dp.locat: numeric vector of two colunms, coordinates of data points
# rp.locat: numeric vector of two colunms, coordinates of regression points
# focus: an integer, indexing to the current regression point, if focus=0, all the distances between all the regression points and data points will be calculated and a distance matrix will be returned;
# if 0<focus<length(rp.locat), then the distances between the focusth regression points and data points will be calculated and a distance vector will be returned;
# p: the power of the Minkowski distance, default is 2, i.e. the Euclidean distance
# theta: an angle in radian to roate the coordinate system, default is 0
# longlat: if TRUE, great circle will be caculated

gw.dist<- function(dp.locat, rp.locat, focus = 0, p = 2, theta = 0, longlat = F)
{
  if (missing(dp.locat)||!is.numeric(dp.locat)||dim(dp.locat)[2]!=2)
  stop("Please input correct coordinates of data points")
  
  if (!missing(rp.locat)) {
    if (!is.numeric(rp.locat))
      stop("Please input correct coordinates of regression points")
    else
      rp.locat <- matrix(rp.locat, ncol = 2)
    if (focus < 0 || focus > length(rp.locat[,1]))
      stop("No regression point is fixed")
    dists <- gw_dist(dp.locat, rp.locat, focus - 1, ifelse(is.infinite(p), -1, p), theta, longlat, T)
  }
  else dists <- gw_dist(dp.locat, dp.locat, focus - 1, ifelse(is.infinite(p), -1, p), theta, longlat, F)
  dists
}
