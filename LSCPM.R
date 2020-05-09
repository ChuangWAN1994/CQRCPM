
###########################
rm(list=ls())

LS.rho=function(t,y,x,z){
  u <- pmin(x-t,0)
  v <- pmax(x-t,0)
  XREG <- cbind(1,u,v,z)
  #fit <- rq.no.cross(y,XREG , taus=tao.v)
  fit <- lm.fit(XREG,y)
  #coef <- fit$coefficients
  return(sum(fit$residuals^2))
}

LS.fitting=function(y,x,z){
  n=length(y)
  resid=NULL
  x.sort=sort(x)
  x_l=x.sort[2]
  x_u=x.sort[n-1]
  minimization=optimize(f=LS.rho, interval=c(x_l,x_u),y=y,x=x,z=z)
  t.star=minimization$minimum
  u <- pmin(x-t.star,0)
  v <- pmax(x-t.star,0)
  #fit <- rq.no.cross(y,cbind(u,v,z),taus=tao.v)
  fit <- lm.fit(cbind(1,u,v,z),y)
  coef <- fit$coefficients
  vartheta <- c(coef,t.star)
  resid <- fit$residuals
  return(list(vartheta=vartheta,resid=resid))
}



#########Inference
#subject <- subject_id

LS.sym <- function(y,x,z,subject)
{
  obj <- LS.fitting(y,x,z)
  vartheta <- obj$vartheta
  res <- obj$resid
  q <- length(vartheta)
  t.est <- vartheta[q]
  bet.est <- vartheta[-q]
  b2 <- bet.est[2]
  b3 <- bet.est[3]
  n <- length(y)
  
  
  #p <- dim(cbind(1,x,z))+1
  #z <- cbind(z)
  #p <- dim(z)[2]
  
  design=model.matrix(~-1+factor(subject))
  N=dim(design)[2]
  x.lst=list()
  z.lst=list()
  resid.lst=list()
  nn=rep(0,N)
  for(i in 1:N){
    x.lst[[i]]=subset(x,design[,i]==1)
    z.lst[[i]]=subset(z,design[,i]==1)
    resid.lst[[i]]=subset(res,design[,i]==1)
    nn[i]=length(resid.lst[[i]])[1]
  }
  
  Vt <- cbind(1,pmin(x-t.est,0),pmax(x-t.est,0),z,
              -b2*(x<=t.est)-b3*(x>t.est))
  zmat0 <- 0*z
  Qn1 <- 1/n * t(Vt) %*% Vt
  
  Qn2 <- matrix(0,q,q)
  Qn2[,q] <- 1/n * apply(cbind(0,(x<=t.est),(x>t.est),zmat0,0)*(res), 2, sum)
  Qn2[q,] <- t(Qn2[,q])
  Qn <- Qn1+Qn2
  Qn.inv <- solve(Qn)
  
  Rn1 <- t(Vt) %*% diag(res^2) %*% Vt/n
  
  Rn2 <- matrix(0,q,q)
  
  for(i in 1:N){
    for(j in 1:nn[i]){
      for(j.pr in (1:nn[i])[-j]){
        mwj <- c(1,pmin(x.lst[[i]][j]-t.est,0),
                 pmax(x.lst[[i]][j]-t.est,0),z.lst[[i]][j],
                 -b2*(x.lst[[i]][j]<=t.est)-b3*(x.lst[[i]][j]>t.est))
        mwj.pr <- c(1,pmin(x.lst[[i]][j.pr]-t.est,0),
                    pmax(x.lst[[i]][j.pr]-t.est,0),z.lst[[i]][j.pr],
                    -b2*(x.lst[[i]][j.pr]<=t.est)-b3*(x.lst[[i]][j.pr]>t.est))
        Rn2 <- Rn2+mwj %*% t(mwj.pr)*(resid.lst[[i]][j]*resid.lst[[i]][j.pr])
      }
    }
  }
  Rn2 <- Rn2/n
  
  Rn <- Rn1 + Rn2
  
  #Cn <- Cn/n
  #diag(Cn)
  Vn <- Qn.inv%*%Rn%*%t(Qn.inv)
  #vtheta=array(vartheta.h,dim=c(length(vartheta.h),1))
  ese <- sqrt(diag(Vn/(n)))
  
  return(ese)
  
}



LS.Boot <- function(y,x,z,subject,n.boot)
{
    design=model.matrix(~-1+factor(subject))
    N=dim(design)[2]
    x.lst=list()
    z.lst=list()
    y.lst=list()
    for(i in 1:N){
      x.lst[[i]]=subset(x,design[,i]==1)
      z.lst[[i]]=subset(z,design[,i]==1)
      y.lst[[i]]=subset(y,design[,i]==1)
    }
    p <- dim(cbind(z))[2]
    par.boot <- matrix(NA,n.boot,5)
    for(i in 1:n.boot){
      N.boot <- sample(1:N,replace = T)
      x.boot <- unlist(x.lst[N.boot])
      z.boot <- unlist(z.lst[N.boot])
      y.boot <- unlist(y.lst[N.boot])
      obj <- LS.fitting(y.boot,x.boot,z.boot) #fitting.cqr(tao.v,y.boot,x.boot,z.boot)
      par.boot[i,] <- obj$vartheta
    }
    
    return(par.boot)
    
}



LAD.boot <- function(y,x,z,subject,tau,n.boot)
{
  
  design=model.matrix(~-1+factor(subject))
  N=dim(design)[2]
  x.lst=list()
  z.lst=list()
  y.lst=list()
  for(i in 1:N){
    x.lst[[i]]=subset(x,design[,i]==1)
    z.lst[[i]]=subset(z,design[,i]==1)
    y.lst[[i]]=subset(y,design[,i]==1)
  }
  p <- dim(cbind(z))[2]
  par.boot <- matrix(NA,n.boot,5)
  for(i in 1:n.boot){
    N.boot <- sample(1:N,replace = T)
    x.boot <- unlist(x.lst[N.boot])
    z.boot <- unlist(z.lst[N.boot])
    y.boot <- unlist(y.lst[N.boot])
    obj <- fitting(tau,y.boot,x.boot,z.boot) #fitting.cqr(tao.v,y.boot,x.boot,z.boot)
    par.boot[i,] <- as.vector(obj[[1]])
  }
  
  return(par.boot)
}


CQR.boot <- function(y,x,z,subject,tao.v,n.boot)
{
  
  design=model.matrix(~-1+factor(subject))
  N=dim(design)[2]
  x.lst=list()
  z.lst=list()
  y.lst=list()
  for(i in 1:N){
    x.lst[[i]]=subset(x,design[,i]==1)
    z.lst[[i]]=subset(z,design[,i]==1)
    y.lst[[i]]=subset(y,design[,i]==1)
  }
  p <- dim(cbind(z))[2]
  par.boot <- matrix(NA,n.boot,(4*length(tao.v)+1))
  for(i in 1:n.boot){
    N.boot <- sample(1:N,replace = T)
    x.boot <- unlist(x.lst[N.boot])
    z.boot <- unlist(z.lst[N.boot])
    y.boot <- unlist(y.lst[N.boot])
    obj <- fitting.cqr(tao.v,y.boot,x.boot,z.boot) #fitting.cqr(tao.v,y.boot,x.boot,z.boot)
    par.boot[i,] <- as.vector(obj$vartheta)
  }
  return(par.boot)
}




