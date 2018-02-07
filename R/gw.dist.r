################################################################################
# gw.dist: caculate a distance vector or distance matrix between (a) regression point and data points
# dp.locat: numeric vector of two colunms, coordinates of data points
# rp.locat: numeric vector of two colunms, coordinates of regression points
# focus: an integer, indexing to the current regression point, if focus=0, all the distances between all the regression points and data points will be calculated and a distance matrix will be returned;
# if 0<focus<length(rp.locat), then the distances between the focusth regression points and data points will be calculated and a distance vector will be returned;
# p: the power of the Minkowski distance, default is 2, i.e. the Euclidean distance
# theta: an angle in radian to roate the coordinate system, default is 0
# longlat: if TRUE, great circle will be caculated

 gw.dist<- function(dp.locat, rp.locat, focus=0, p=2, theta=0, longlat=F)
 {
   if (missing(dp.locat)||!is.numeric(dp.locat)||dim(dp.locat)[2]!=2)
      stop("Please input correct coordinates of data points")
  
   if (!missing(rp.locat)) 
     {
       rp.given<-T
     }
    else
     {
       rp.given<-F 
       rp.locat<- dp.locat
     } 
   
   if (!is.numeric(rp.locat))
      stop("Please input correct coordinates of regression points")
   else
      rp.locat <- matrix(rp.locat, ncol=2)
   if (focus<0||focus>length(rp.locat[,1]))
      stop("No regression point is fixed")

   n.rp<-length(rp.locat[,1])
   n.dp<-length(dp.locat[,1])
   if (focus>0)
       dists<-numeric(n.dp)
   else
       dists<-matrix(numeric(n.rp*n.dp),nrow=n.dp)
   ###Rotate the coordiante system
   if (p!=2&&theta!=0&&!longlat)
   {
     dp.locat<-coordinate_rotate(dp.locat, theta)
     rp.locat<-coordinate_rotate(rp.locat, theta)
   }

   ####Calculate a distance vector at regression point "focus"
   if (focus>0)
   {
      if (longlat)   ###Great circle distance, and FORTRAN or C++ code to be added
         dists<-spDistsN1(dp.locat, rp.locat[focus,],longlat=longlat)
      else
      { 
         if(p==2)
         {
             dists<-eu_dist_vec(dp.locat,rp.locat[focus,])
         }
         else if(is.infinite(p))
             dists<-cd_dist_vec(dp.locat,rp.locat[focus,])
         else if(p==1)
             dists<-md_dist_vec(dp.locat,rp.locat[focus,])
         else
             dists<-mk_dist_vec(dp.locat,rp.locat[focus,],p) 
      }              
   }
   else
   {
      if (longlat)
      {
         for (i in 1:n.rp)
           dists[,i]<-spDistsN1(dp.locat, matrix(rp.locat[i,],nrow=1),longlat=longlat)
      }
      else
      {
         if(p==2)
         {
             if(rp.given)
                dists<-eu_dist_mat(dp.locat,rp.locat)
             else
                dists<-eu_dist_smat(dp.locat)
         }
         else if(is.infinite(p))
         {
             if(rp.given)
                dists<-cd_dist_mat(dp.locat,rp.locat)
             else
                dists<-cd_dist_smat(dp.locat)
             
         }    
         else if(p==1)
         {
             if(rp.given)
                dists<-md_dist_mat(dp.locat,rp.locat)
             else
                dists<-md_dist_smat(dp.locat)
         }
         else
         {
             if(rp.given)
                dists<-mk_dist_mat(dp.locat,rp.locat,p)
             else
                dists<-mk_dist_smat(dp.locat,p)  
         }
      }     
   }
   dists     

 }
