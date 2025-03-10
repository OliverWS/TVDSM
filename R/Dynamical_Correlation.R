########################################################
### Estimating and Testing for Dynamical Correlation ###
########################################################

# Use the ind_DC function below to estimate dynamical correlations for each individual (Eq. 6)
# The sample-level dynamical correlation can then be computed by averaging the individual dynamical correlations (Eq. 7)

#install.packages("caTools")                             	### load package {caTools} to call 'trapz' ###                                      
require(caTools)   

ind_DC=function(x,y,t,na=FALSE){                     	### input x, y (data matrices with rows representing measurements and columns representing subjects) 
  #               t    (number of time points)
  #               na   (a logical variable indicating whether there are missing data)                        
  
  
  if (na) {                                   		  	 ### impute missing values by linear interpolation ###
    for (i in 1:dim(x)[2]){ 
      x[,i]=approx(t,x[,i],xout=t,rule=2)$y      
      y[,i]=approx(t,y[,i],xout=t,rule=2)$y
    }
  }
  
  temp1_x=temp1_y=matrix(0,nrow=dim(x)[1],ncol=dim(x)[2])
  temp2_x=temp2_y=matrix(0,nrow=dim(x)[1],ncol=dim(x)[2])
  M_x=M_y=z=numeric()
  
  for (i in 1:dim(x)[2]){
    aver_x=trapz(t,x[,i])/(tail(t,1)-head(t,1))     	### integrals approximated by trapezoidal rule ###
    aver_y=trapz(t,y[,i])/(tail(t,1)-head(t,1))
    temp1_x[,i]=x[,i]-aver_x                               	### center within subjects (vertical shift; Eq. 3) ###
    temp1_y[,i]=y[,i]-aver_y
  }
  
  M_x=rowMeans(temp1_x)                               	### the population-mean of vertically shifted curves (Eq. 4) ###
  M_y=rowMeans(temp1_y) 
  for (i in 1:dim(x)[2]){
    temp2_x[,i]=(temp1_x[,i]-M_x)/sqrt(trapz(t,(temp1_x[,i]-M_x)^2)/(tail(t,1)-head(t,1)))       	### population-centering and standardization (Eq. 5) ###
    temp2_y[,i]=(temp1_y[,i]-M_y)/sqrt(trapz(t,(temp1_y[,i]-M_y)^2)/(tail(t,1)-head(t,1)))
    z[i]=trapz(t,temp2_x[,i]*temp2_y[,i])/(tail(t,1)-head(t,1))   					### obtain individual dynamical correlation ###
  }
  return(z)
}
# Use the boot_test_DC function to conduct significance testing 

boot_test_DC=function(x1,y1,t1,x2,y2,t2,ms=c(FALSE,FALSE),B=1000,randPairs=F){     
  ### input x1,y1 (data matrices with rows representing measurements and columns representing subjects)
  #         x2,y2 (paired data matrices (if there are any) with rows representing measurements and columns representing subjects)
  #         t1,t2 (number of time points)
  #         ms    (a vector of logical variables indicating whether there are missings for paired data matrices) 
  #         B     (the number of bootstrap samples)                                                                            ###
  
  n=dim(x1)[2]
  if (missing(x2)) {                                           			### one-sample test ###
    dyncor_1=ind_DC(x1,y1,t1,na=ms[1])                                 	### observed DC ###
    obs_1_stud=mean(dyncor_1)*sqrt(n)/sd(dyncor_1)              	### observed standardized version ###
    boot_1_stud=numeric()
    boot_x1=boot_y1=matrix(0,nrow=length(t1),ncol=n)
    
    for(b in 1:B){ 
      idx=sample(c(1:n),replace = (randPairs==F))                         
      if(randPairs){
        boot_x1=x1[,1:n]  
        boot_y1=y1[,idx]
        boot_dyncor_1=ind_DC(boot_x1,boot_y1,t1,na=ms[1])                 	 ### DC based on bootstrap samples ###
        ### bootstrap replicates ###
        boot_1_stud[b]=mean(boot_dyncor_1)*sqrt(n)/sd(boot_dyncor_1)    	### bootstrap standardized version ###       
        
      }
      else {
        boot_x1=x1[,idx]    
        boot_y1=y1[,idx]
        boot_dyncor_1=ind_DC(boot_x1,boot_y1,t1,na=ms[1])                 	 ### DC based on bootstrap samples ###
        ### bootstrap replicates ###
        boot_1_stud[b]=mean(boot_dyncor_1-dyncor_1)*sqrt(n)/sd(boot_dyncor_1)    	### bootstrap standardized version ###       
        
      }
    }
    
    emp.stat=obs_1_stud                                          		### bootstrap test statistic ###
    emp.pval=length(boot_1_stud[boot_1_stud>abs(obs_1_stud) | boot_1_stud< -abs(obs_1_stud)])/B    ### p-value based on bootstrap null distribution ###
    outp=list(r=mean(dyncor_1),stats=emp.stat, pval=emp.pval,boot_1_stud=boot_1_stud,dyncor_1=dyncor_1,x1=x1,y1=y1,t=t1) 
    return(outp)
  }  else {                                                    			### two-sample paired test (similar as above) ###
    dyncor_1=ind_DC(x1,y1,t1,na=ms[1])
    dyncor_2=ind_DC(x2,y2,t2,na=ms[2])
    obs_1_stud=mean(dyncor_1)*sqrt(n)/sd(dyncor_1)
    obs_2_stud=mean(dyncor_2)*sqrt(n)/sd(dyncor_2)
    obsdiff=dyncor_2-dyncor_1
    obsdiff_stud=mean(obsdiff)*sqrt(n)/sd(obsdiff)
    boot_1_stud=boot_2_stud=boot_diff_stud=numeric()
    
    boot_x1=boot_y1=matrix(0,nrow=length(t1),ncol=n)
    boot_x2=boot_y2=matrix(0,nrow=length(t2),ncol=n)
    
    for(b in 1:B){
      idx=sample(c(1:n),replace = TRUE)
      boot_x1=x1[,idx]
      boot_y1=y1[,idx]
      boot_x2=x2[,idx]
      boot_y2=y2[,idx]
      boot_dyncor_1=ind_DC(boot_x1,boot_y1,t1)
      boot_dyncor_2=ind_DC(boot_x2,boot_y2,t2)
      boot_1_stud[b]=mean(boot_dyncor_1-dyncor_1)*sqrt(n)/sd(boot_dyncor_1)
      boot_2_stud[b]=mean(boot_dyncor_2-dyncor_2)*sqrt(n)/sd(boot_dyncor_2)
      boot_diff_stud[b]=mean(boot_dyncor_2-boot_dyncor_1-obsdiff)*sqrt(n)/sd(boot_dyncor_2-boot_dyncor_1)
    }
    
    emp.stat=obsdiff_stud
    emp.pval=length(boot_diff_stud[boot_diff_stud>abs(obsdiff_stud) | boot_diff_stud< -abs(obsdiff_stud)])/B
    outp=list(stats=emp.stat, pval=emp.pval,idc=dyncor_1) 
    return(outp)
  } 
}
# 
# ##########################
# ### Simulation Example###
# ##########################
# require(MASS)
# ###### DC=1 ######
# n=100             # sample size
# t=seq(0,1,length.out=100)       # length of data
# 
# # fixed effects models for x and y
# mu_quad_x=8*t^2-4*t+5
# mu_quad_y=8*t^2-12*t+6
# 
# ###### simulation ######
# fun=rbind(rep(1,length(t)),-t,t^2)
# z1=mvrnorm(n,rep(0,3),diag(c(2,16/3,4)))   #covariance matrix of random effects
# x1_quad=y1_quad=x1_quad_error=y1_quad_error=matrix(0,nrow=length(t),ncol=n)
# 
# 
# for (i in 1:n){
#   x1_quad[,i]=mu_quad_x+z1[i,]%*%fun           # x without error
#   y1_quad[,i]=mu_quad_y+2*z1[i,]%*%fun         # y without error
#   x1_quad_error[,i]=x1_quad[,i]+rnorm(length(t),0,0.01)      # x with error
#   y1_quad_error[,i]=y1_quad[,i]+rnorm(length(t),0,0.01)      # y  with error
# }
# 
# ###### estimating DC ######
# dyn1_quad=ind_DC(x1_quad,y1_quad,t)                     
# dyn1_quad         
# mean(dyn1_quad)    #DC=1 because there is no error
# 
# dyn1_quad_error=ind_DC(x1_quad_error,y1_quad_error,t)         
# dyn1_quad_error
# mean(dyn1_quad_error)    #DC close to 1 because of error
# 
# ###### test H0: DC=0 ######
# bt_DC=boot_test_DC(x1_quad_error,y1_quad_error,t,ms=FALSE,B=1000)        ### use bootstrap sample 1000 (default) for example ###
# bt_DC                                                        ### test results ###
# bt_DC$stats                                                  ### test statistic ###
# bt_DC$pval                                                   ### p-value ###
# 
# ###### test H0: DC1=DC2 ###### If there is another set of paired curves x2,y2,t2 (assuming no missing data)
# bt2_DC=boot_test_DC(x1_quad,y1_quad,t,x2,y2,t2,ms=c(FALSE,FALSE),B=1000)
# 
# a_mat = matrix(nrow=100,ncol=10)
# for(i in 1:100){a_mat[i,] <- 1:10}
# 
# b_mat = matrix(nrow=100,ncol=10)
# for(i in 1:100){b_mat[i,] <- 1:10}
# 
# n = 10
# for(b in 1:10){                                               
#   idx=sample(c(1:n),replace = TRUE)     
#   print(idx)
#   boot_x1=a_mat[,idx]                                        			### bootstrap replicates ###
#   boot_y1=b_mat[,idx]
#   #print(paste("Diff=",sum(boot_x1 - boot_y1)))
# }
