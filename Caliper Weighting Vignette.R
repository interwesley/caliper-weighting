library(sbw)


#######functions########
"transform"=function (formula, data){
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "numeric")
  T <-model.matrix(mt, mf)[,2]
  X <- model.matrix(mt, mf)[,-c(1,2)]
  data.frame(Y,T,X)
}



"caliper_w"=function(formula,data,w,h){
  dt=transform(formula,data)
  n=dim(dt)[1]
  
  k=length(w)
  y_hat=seq(k)
  t=matrix(0,n,k)
  ws=matrix(0,n,k)
  Ys=matrix(0,n,k)
  ord=matrix(0,n,k)
  
  for(i in 1:k){
    t[abs(dt[,2]-w[i])<h,i]=1
  }
  
  for(i in 1:k){
    Zt=cbind(t[,i],dt)
    Zt = Zt[order(t[,i]), ]
    ord[,i]=order(t[,i])
    X=Zt[,c(-1,-2,-3)]
    t[,i]=Zt[,1]
    Ys[,i]=Zt[,2]
    
    t[,i]=1-t[,i]
    
    
    data_frame = as.data.frame(cbind(t[,i], X,Ys[,i]))
    names(data_frame)[1] = "t_ind"
    names(data_frame)[dim(dt)[2]] = "Y"
    t_ind = "t_ind"
    # moment covariates
    bal = list()
    bal$bal_cov = colnames(X)
    bal$bal_tol = 0.02
    bal$bal_std = "group"
    
    tar = colMeans(X)
    names(tar) = bal$bal_cov
    sbwpop_object = sbw(dat = data_frame, ind = t_ind, out = "Y", bal = bal,
                        sol = list(sol_nam = "quadprog"), 
                        par = list(par_est = "pop"))
    
    ws[,i]=sbwpop_object[["dat_weights"]][["sbw_weights"]]
    t[,i]=1-t[,i]
    y_hat[i]=sum(ws[,i]*t[,i]*Ys[,i])
    
  }
  
  
  
  return(list(Ys=Ys,t=t,ws=ws,y_hat=y_hat,ord=ord))
}


cband=function(estimate,w){
  Ys=estimate[[1]]
  t=estimate[[2]]
  ws=estimate[[3]]
  y=estimate[[4]]
  
  
  
  sg=-Inf*seq(2000)
  for(i in 1:2000){
    kosi=rnorm(nrow(t),0,1)
    for(j in 1:length(w)){
      sgs=sum(kosi*Ys[,j]*t[,j]*ws[,j])-mean(kosi)*y[j];
      if(sgs>sg[i]){
        sg[i]=sgs
      }
    }}
  
  ssg=sort(sg)
  return(ssg[1900])
  
}

simulate_data=function(n){
  X1 = rbeta(n,2,2)+rnorm(n,0,1)
  X2=rgamma(n,3,1)/7
  X3= rexp(n)-runif(n,-1,1)
  W = log(abs(X1*X3))+X2+rnorm(n,0,2)+5
  Y=sin(W)+0.05*W^2+0.2*rnorm(n,0,1)-0.04*W*(X1+X2-0.5-3/7) +0.05*(X3^2-2-1/3)
  w=seq(0,9,0.2)
  k=length(w)
  Sample=data.frame(cbind(Y,W,X1,X2,X3))
  return(Sample)
}



########Instructions for the caliper_w function#######

#Input:
  
#formula: An expression for the relationship between outcome, treatment and covariates. The format is outcome~treatment+covariate_1+covariate_2+...+covariate_n.

#data: A data frame with an outcome variable, a treatment variable, and a series of covariates. It's typically in the form of a matrix.

#w: The treatment values at which we want to estimate the corresponding average treatment effect. It's typically a vector of values.

#h: The bandwidth.

#Output:
  
#y_hat is the estimated average treatment effect at each treatment level.



############Example############

n=5000
dataset=simulate_data(n)
w=seq(0,12,0.2)
h=12*n^(-5/22)

caliper_estimate=caliper_w(Y~W+X1+X2+X3,dataset,w,h)

y_hat=caliper_estimate$y_hat


##########Instructions for the cband function###########

#Input:
  
#estimate: The output of the caliper_w function

#w: The treatment values at which we want to estimate the corresponding average treatment effect. It should be the same as the input treatment values of the caliper_w function.

#Output: A single value. It's the width of the confidence band.

#################Example####################

alpha=cband(caliper_estimate,w)





#Plot true ADRF and estimated ADRF

plot(w,y_hat,type="l",ylim = c(-2,7.5),xlab="treatment levels",
     ylab="ADRF")
lines(w,y_hat+alpha,lty=2)
lines(w,y_hat-alpha,lty=2)


z=sin(w)+0.05*w^2
lines(w,z,col="green")

legend("bottomright", c("True ADRF","Estimated ADRF","95% confidence band"),
       lty=c(1,2,1),col=c("black","green","black"),cex=0.8)
