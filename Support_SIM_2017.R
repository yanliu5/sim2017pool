
################################################################################################
# MC method; 2017 SIM Liu
# MCEM method; 2014 Biometrics Mitchell
################################################################################################

data_sim_para <- function(weight=TRUE,beta,sig,N,cj,J){
  
  # beta <- c(1,-0.5,1,0.5,-0.5,0.5,0.5);sig  <- 0.5
  # Generate 2 continuous covariates that are nonlinearly related to E(Y)
  
  x_1   <- rbinom(N,1,0.5)      # Gender
  race  <- t(rmultinom(N, 1, c(0.45,0.25,0.30)))[,-1]
  x_3   <- runif(N,-1,1)        # Age
  x_4   <- rnorm(N)             # BMI
  
  Wmat  <- cbind(1,x_1,race,x_3,x_4) # Design matrix
  
  Omat  <- matrix(NA,J,max(cj))                  # Weight matrix
  Vol   <- rep(N,1)                              # Generate volume of each specimen
  Wvec  <- rep(NA,J)
  
  # ====================================== Homogeneous pool
  Omat_HG  <- matrix(NA,J,max(cj))                  # Weight matrix of Homegeneous pooling
  oid      <- order(Wmat[,2],Wmat[,3],Wmat[,4],Wmat[,5])
  Wmat_HG  <- Wmat[oid,]
  Vol_HG   <- Vol[oid]
  Wvec_HG  <- rep(NA,J)
  #head(cbind(Wmat,x_1),20)
  #head(cbind(Wmat_HG,x_1_HG),20)
  
  csm <- 0
  for(jj in 1:J){
    if(weight){
      x    <- Vol[(csm+1):(csm+cj[jj])]
      x_HG <- Vol_HG[(csm+1):(csm+cj[jj])]
    }else{
      x    <- rep(1,cj[jj])
      x_HG <- rep(1,cj[jj])
    }
    Omat[jj,]    <- x#/sum(x)
    Omat_HG[jj,] <- x_HG#/sum(x_HG)
    Wvec[jj]     <- sum(x)
    Wvec_HG[jj]  <- sum(x_HG)
    csm <- csm + cj[jj]
  }
  
  Y_ind <- rlnorm(N,Wmat%*%beta,sig)
  Y_ind_HG <- Y_ind[oid]
  
  # create true pooled biomarker level
  Pool_id <- rep(1:J,cj)
  Imat    <- cbind(0,Pool_id,c(t(Omat)))
  Imat_HG <- cbind(0,Pool_id,c(t(Omat_HG)))
  
  Y_pool    <- aggregate(Y_ind*Imat[,3]~Imat[,2],FUN=sum)[,2]/Wvec
  Y_pool_HG <- aggregate(Y_ind_HG*Imat_HG[,3]~Imat_HG[,2],FUN=sum)[,2]/Wvec_HG
  # sum(abs(Y_pool1-Y_pool)) # should be 0
  Zmat    <- matrix(NA,J,2+max(cj))
  Zmat_HG <- matrix(NA,J,2+max(cj))
  j1 <- 0
  for(jj in 1:J){
    Zmat[jj,]    <- c(Y_pool[jj],cj[jj],(j1+1):(j1+cj[jj]))
    Zmat_HG[jj,] <- c(Y_pool_HG[jj],cj[jj],(j1+1):(j1+cj[jj]))
    j1 <- j1+cj[jj]
  }
  
  return(list(Y_ind=Y_ind,Y_pool=Y_pool,Wmat=Wmat,Omat=Omat,Imat=Imat,
              Zmat=Zmat,Y_ind_HG=Y_ind_HG,Y_pool_HG=Y_pool_HG,
              Wmat_HG=Wmat_HG,Omat_HG=Omat_HG,Imat_HG=Imat_HG,
              Zmat_HG=Zmat_HG,Wvec=Wvec,Wvec_HG=Wvec_HG))
}

###### Define Log-likelihood function for pooled data
llh<-function(theta,theta0,y,X,m,samp0){
  beta<-theta[-length(theta)]
  environment(mc)  <- environment()
  llmu<-X%*%beta  ## X design matrix
  llsig<-theta[length(theta)]
  beta0<-theta0[-length(theta0)]
  llmu0<-X%*%beta0  ## X design matrix
  llsig0<-theta0[length(theta0)]
  res<-sum(sapply(seq(1,gp),mc,y,llmu,llsig,llmu0,llsig0,m,samp0))
  return(res)
}
###### compute mc function
mc<-function(ind,y,llmu,llsig,llmu0,llsig0,m,samp0){
  index<-which(kmc$cluster==ind)
  kj<-length(index)
  mu.ind<-llmu[index]
  mu.ind0<-llmu0[index]
  y.p<-y[ind]
  y.mc0<-samp0[[ind]]
  if(kj==1){val<- -log(dlnorm(y.p,mu.ind,llsig))} else{
    z.mat<-y.mc0/rowSums(y.mc0)*kj*y.p
    w<-exp(rowSums(log(t(dlnorm(t(z.mat),mu.ind,llsig)/dlnorm(t(z.mat),mu.ind0,llsig0)))))
    z.sum<-rowSums(t(t(log(z.mat))-mu.ind0))
    cj<-1/y.p*sqrt(kj/(2*pi*llsig0^2))*exp(-z.sum^2/(2*kj*llsig0^2))
    cw<-cj*w
    val<- -log(mean(cw))
  }
  if (is.infinite(val)) {val<-2000000}
  return(val)
}
###### gradient theoretic formula
M.J<-function(id,llmu,llsig,llmu0,llsig0,X,y,samp0){
  index<-which(kmc$cluster==id)
  kj<-length(index)
  y.p<-y[id]
  y.mc0<-samp0[[id]]
  mu.ind<-llmu[index]
  mu.ind0<-llmu0[index]
  if(kj==1){
    w<-y.mc0
    cj<-dlnorm(y.p,mu.ind,llsig)
    a1<-t(1/llsig^2*(log(y.p)-mu.ind)*(X[index,])%*%t(w))
    a11<- -1/llsig+1/llsig^3*(log(y.p)-mu.ind)^2*w
  } else{
    z.mat<-y.mc0/rowSums(y.mc0)*kj*y.p
    w<-exp(rowSums(log(t(dlnorm(t(z.mat),mu.ind,llsig)/dlnorm(t(z.mat),mu.ind0,llsig0)))))
    z.sum<-rowSums(t(t(log(z.mat))-mu.ind0))
    cj<-1/y.p*sqrt(kj/(2*pi*llsig0^2))*exp(-z.sum^2/(2*kj*llsig0^2))
    u.sum<-(t(t(log(z.mat))-mu.ind))
    a1<-t(1/llsig^2*(t(X[index,])%*%t(u.sum)))
    a11<- -kj/llsig+1/(llsig^3)*rowSums(u.sum^2)
  }
  cw<-cj*w
  dcw<-(as.vector(cw)*cbind(a1,a11))
  res<-cbind(dcw,cw)
  return(res)
}
######
######
grr<-function(theta,theta0,y.pool,X.mat,m,samp0){
  p <- length(theta0)
  environment(M.J)  <- environment()
  beta<-theta[-length(theta)]
  llmu<-X.mat%*%beta  ## X design matrix
  llsig<-theta[length(theta)]
  beta0<-theta0[-length(theta0)]
  llmu0<-X.mat%*%beta0  ## X design matrix
  llsig0<-theta0[length(theta0)]
  mmat<-function(id){
    out<-M.J(id,llmu,llsig,llmu0,llsig0,X.mat,y.pool,samp0)
    return(colMeans(out))
  }
  res<-sapply(seq(1,gp),mmat)   
  return(-rowSums(t(t(res[1:p,])/res[(p+1),])))
}
############ Generate MC sample
samp<-function(theta0,X.mat,kmc,m){
  sam<-list()
  for (ind in seq(1,gp)){
    beta0<-theta0[-length(theta0)]
    llmu0<-X.mat%*%beta0  ## X design matrix
    llsig0<-theta0[length(theta0)]
    index<-which(kmc$cluster==ind)
    kj<-length(index)
    mu.ind0<-llmu0[index]
    if(kj==1){
      sam[[ind]]=rep(1,m)
    }else{
      sam[[ind]]<-matrix(rlnorm(m*kj,mu.ind0,llsig0),ncol=kj,byrow=T)
    }}
  return(sam)
}
######

############Variance Estimate 2##########
var.hat<-function(theta,theta0,X.mat,y.pool,mm,D2){
  environment(M.J)  <- environment()
  p<-length(theta)
  beta<-theta[-length(theta)]
  llmu<-X.mat%*%beta  ## X design matrix
  llsig<-theta[length(theta)]
  beta0<-theta0[-length(theta0)]
  llmu0<-X.mat%*%beta0  ## X design matrix
  llsig0<-theta0[length(theta0)]
  C.mat<-matrix(rep(0,p*gp*(p+1)),nrow=p)
  samp00<-samp(theta0,X.mat,kmc,mm)
  out.mat<-list()
  cov.mat<-list()
  for (id in seq(1,gp)){
    out<-M.J(id,llmu,llsig,llmu0,llsig0,X.mat,y.pool,samp00)
    out.mat[[id]]<-out
    cov.mat[[id]]<-cov(out)
    q.bar<-mean(out[,(p+1)])
    for (kk in seq(1,p)) {
      vkk.bar<-mean(out[,kk])
      C.mat[kk,(id-1)*(p+1)+kk]<- 1/q.bar
      C.mat[kk,id*(p+1)]<- -vkk.bar/q.bar^2
    }}
  H.mat<-bdiag(cov.mat)
  D1<-C.mat%*%H.mat%*%t(C.mat)
  Var.mat<-diag(solve(D2)%*%D1%*%solve(D2))
  return(Var.mat)
}
######################################
###### Define Log-likelihood function for pooled data
llh.em<-function(theta,y,X,m){
  beta<-theta[-length(theta)]
  llmu<-X%*%beta  ## X design matrix
  llsig<-theta[length(theta)]
  h.y<-sapply(seq(1,gp),mc.em,beta,y,llmu,llsig,m)
}
###### compute MC
mc.em<-function(ind,beta,y,llmu,llsig,m){
  index<-which(kmc$cluster==ind)
  kj<-length(index)
  mu.ind<-llmu[index]
  y.jp<-y[ind]
  yidv.mc<-matrix(rlnorm(m*kj,mu.ind,llsig),ncol=kj,byrow=T)
  z.mc<-yidv.mc/rowSums(yidv.mc)*kj*y.jp
  z.lsum<-rowSums(t(t(log(z.mc))-mu.ind))
  w.j<-exp(-z.lsum^2/(2*llsig^2*kj))
  botm<-sum(w.j)
  h1.y<-colSums(log(z.mc)*w.j)/botm
  h2.y<-colSums((log(z.mc))^2*w.j)/botm
  return(c(h1.y,h2.y))
}
######################## function for last step
llhl.em<-function(theta,y,X,m,samp0){
  beta<-theta[-length(theta)]
  llmu<-X%*%beta  ## X design matrix
  llsig<-theta[length(theta)]
  h.y<-sapply(seq(1,gp),mcl.em,beta,y,llmu,llsig,m,samp0)
}
###### compute MC
mcl.em<-function(ind,beta,y,llmu,llsig,m,samp0){
  index<-which(kmc$cluster==ind)
  kj<-length(index)
  mu.ind<-llmu[index]
  y.jp<-y[ind]
  yidv.mc<-samp0[[ind]]
  z.mc<-yidv.mc/rowSums(yidv.mc)*kj*y.jp
  z.lsum<-rowSums(t(t(log(z.mc))-mu.ind))
  w.j<-exp(-z.lsum^2/(2*llsig^2*kj))
  botm<-sum(w.j)
  h1.y<-colSums(log(z.mc)*w.j)/botm
  h2.y<-colSums((log(z.mc))^2*w.j)/botm
  return(c(h1.y,h2.y))
}
#############################
###### compute h1.y and h2.y
mat.ana<-function(mat){
  h1<-c()
  h2<-c()
  for (i in seq(1,gp)){
    if(class(mat)=="list"){
      a<-mat[[i]]
    }else{ a<-mat[,i]}
    ct<-length(a)
    h1<-c(h1,a[1:(ct/2)])
    h2<-c(h2,a[(ct/2+1):ct])
  }
  return(c(h1,h2))
}
###### initial fit for parameters
inifit<-function(X.pool,y.pool,kmc){
  gp <- length(y.pool)
  XX.pool<-c()
  for ( id in 1:gp){
    index<-which(kmc$cluster==id)
    if(length(index)==1){
      XX.pool<-rbind(XX.pool,X.pool[index,])
    }else{
      XX.pool<-rbind(XX.pool,colMeans(X.pool[index,]))}
  }
  
  K.mat<-diag(kmc$size)
  XX.mat<-cbind(rep(1,gp)+solve(K.mat)%*%rep(1,gp),XX.pool)
  b.hat<-solve(t(XX.mat)%*%K.mat%*%XX.mat)%*%t(XX.mat)%*%K.mat%*%log(y.pool)
  A.mat<-K.mat%*%XX.mat%*%solve(t(XX.mat)%*%K.mat%*%XX.mat)%*%t(XX.mat)%*%K.mat
  c.hat<-t(log(y.pool))%*%(K.mat-A.mat)%*%log(y.pool)/(gp-length(b.hat))
  b.hat[1]<-b.hat[1]+c.hat/2
  lsig.hat<-sqrt(log(1+c.hat))
  theta.ini<-c(b.hat,lsig.hat)
  
  return(theta.ini)
  
}

###################################
###################################
emfit<-function(m,cal,theta0,X.mat,y.pool,kmc){
  p         <- length(theta0)
  gp        <- length(y.pool)
  theta.ini <- theta0
  theta.rc  <- theta.ini
  environment(llhl.em)   <- environment()
  environment(llh.em)    <- environment()
  environment(mat.ana)    <- environment()
  
  t1<-proc.time()  ## record time1
  ic<-1
  while(ic <= cal){
    if(ic%%20==0) m<- ceiling(m + 0.25*m)
    if(ic==cal){
      beta0<-theta.ini[-length(theta.ini)]
      llmu0<-X.mat%*%beta0  ## X design matrix
      llsig0<-theta.ini[length(theta.ini)]
      samp0<-list()
      for (ii in seq(1,gp)){
        index<-which(kmc$cluster==ii)
        kj<-length(index)
        mu.ind0<-llmu0[index]
        samp0[[ii]]<-matrix(rlnorm(m*kj,mu.ind0,abs(llsig0)),ncol=kj,byrow=T)
      }
      
      res<-llhl.em(theta.ini,y.pool,X.mat,m,samp0)
      h.y<-mat.ana(res)
      h1.y<-h.y[1:N]
      h2.y<-h.y[(N+1):(2*N)]
    }else{
      
      res<-llh.em(theta.ini,y.pool,X.mat,m)
      h.y<-mat.ana(res)
      h1.y<-h.y[1:N]
      h2.y<-h.y[(N+1):(2*N)]}
    ###### reorder the design matrix
    X.n<-X.mat
    ###### Calculate parameter
    beta.ini<-solve(t(X.n)%*%X.n)%*%t(X.n)%*%h1.y
    mu.ini<-X.n%*%beta.ini
    lsig.ini<-sqrt(sum(h2.y-2*mu.ini*h1.y+mu.ini^2)/N)
    theta.ini<-c(beta.ini,lsig.ini)
    theta.rc<-c(theta.rc,theta.ini)
    ic<-ic+1
  }
  
  mat<-matrix(theta.rc,ncol=p,byrow=T)
  t2<-proc.time() ## record time3
  ################## Observed Fisher Information ################
  
  ###### compute Q_j
  qj<-function(theta,X,y,ind,yidv.mc){
    beta.t<-theta[-length(theta)]
    sig.t<-theta[length(theta)]
    mu.t<-X%*%beta.t
    theta.qj<-mat[cal,]
    beta.qj<-theta.qj[-length(theta.qj)]
    sig.qj<-theta.qj[length(theta.qj)]
    mu.qj<- X%*%beta.qj  ## X design matrix
    index<-which(kmc$cluster==ind)
    kj<-length(index)
    mu.ind<-mu.qj[index]
    y.jp<-y[ind]
    z.mc<-yidv.mc/rowSums(yidv.mc)*kj*y.jp
    z.lsum<-rowSums(t(t(log(z.mc))-mu.ind))
    w.j<-exp(-z.lsum^2/(2*sig.qj^2*kj))
    botm<-sum(w.j)
    h1.y<-colSums(log(z.mc)*w.j)/botm
    h2.y<-colSums((log(z.mc))^2*w.j)/botm
    val<- -kj*log(sig.t)-1/(2*sig.t^2)*sum(h2.y-2*mu.t[index]*h1.y+mu.t[index]^2)
    return(val)}
  ###### compute the last term of observed information
  ltj.j<-function(idx,theta,X,y.data,kj,y.jp,index){
    beta.t<-theta[-length(theta)]
    sig.t<-theta[length(theta)]
    mu.t<-X%*%beta.t
    yidv.mc<-y.data[idx,]
    z.mc<-yidv.mc/sum(yidv.mc)*kj*y.jp
    z.tsum<-(log(z.mc)-mu.t[index])
    if(kj==1){
      a1<-1/sig.t^2*(t(X[index,]*z.tsum))}else{
        a1<-1/sig.t^2*(t(X[index,]*z.tsum)%*%rep(1,kj))
      }
    a2<- -kj/sig.t+1/sig.t^3*sum(z.tsum^2)
    aa<-c(a1,a2)
    return(c(aa%*%t(aa)))
  }
  
  ltj<-function(theta,X,y,ind,samp0){
    theta.qj<-mat[cal,]
    beta.qj<-theta.qj[-length(theta.qj)]
    sig.qj<-theta.qj[length(theta.qj)]
    mu.qj<- X%*%beta.qj  ## X design matrix
    index<-which(kmc$cluster==ind)
    kj<-length(index)
    mu.ind<-mu.qj[index]
    y.jp<-y[ind]
    #y.data<-matrix(rlnorm(m*kj,mu.ind,sig.qj),ncol=kj,byrow=T)
    y.data<-samp0[[ind]]
    z.mc<-y.data/rowSums(y.data)*kj*y.jp
    test<-t(log(z.mc))-mu.ind
    z.lsum<-rowSums(t(test))
    w.j<-exp(-z.lsum^2/(2*sig.qj^2*kj))
    botm<-sum(w.j)
    res<-sapply(seq(1,m),ltj.j,mat[cal+1,],X.mat,y.data,kj,y.jp,index)
    res1<-colSums(t(res)*w.j)/botm
    val<-matrix(res1,ncol=length(theta.ini))
    return(val)
  }
  
  ress<-c()
  for ( ii in seq(1,gp)){
    QQ1<-function(theta){
      qj(theta,X.mat,y.pool,ind=ii,samp0[[ii]])
    }
    gad<-grad(QQ1,mat[cal+1,])
    hes<-hessian(QQ1,mat[cal+1,])
    lst<-ltj(mat[cal+1,],X.mat,y.pool,ii,samp0)
    ress<-cbind(ress,c(hes-gad%*%t(gad)+lst))
  }
  res<-rowSums(ress)
  mat.fish<-matrix(res,ncol=length(theta0))
  hh<-mat.fish
  t3<-proc.time() ## record time3
  
  t<-as.vector(t2-t1)[3]
  tt<-as.vector(t3-t2)[3]
  #res<-c(mat[cal+1,],hh,c(t,tt))
  
  return(list(theta=mat[cal+1,],hess=hh,clock=t+tt))
}

####################################
####################################
mcfit<-function(m,mm,theta0,X.mat,y.pool,kmc){
  status <- 1
  p  <- length(theta0)
  gp <- length(y.pool)
  
  environment(llh)      <- environment()
  environment(grr)      <- environment()
  environment(var.hat)  <- environment()
  environment(samp)     <- environment()
  
  t1<-proc.time()
  samp0<-samp(theta0,X.mat,kmc,m)
  tryCatch({r1<-optim(theta0,llh,grr,theta0=theta0,y=y.pool,X=X.mat,m=m,samp0=samp0,method="BFGS",hessian=T)}, error=function(e){status <<- 0})
  
  if(status==0){
    res<- c(status,rep(0,p^2+p*2+4))
    return(res)
  }
  
  
  
  theta00<-r1$par
  t2<-proc.time()
  t<-as.vector(t2-t1)[3]
  Var.mat<-var.hat(theta00,theta00,X.mat,y.pool,mm,r1$hessian)
  t3<-proc.time()
  
  tt<-as.vector(t3-t2)[3]
  if(sqrt(max(Var.mat/m)) < dd){
    m1<-max(Var.mat/dd^2)
    m2<-2000
    ttt<-0
    theta<-theta00
    hes<-r1$hessian
  }else{
    m1<-(max(Var.mat/dd^2)%/%10+1)*10
    Var.mat<-var.hat(theta00,theta00,X.mat,y.pool,mm,r1$hessian)
    m2<-(max(Var.mat/dd^2)%/%10+1)*10
    print(m2)
    samp0<-samp(theta00,X.mat,kmc,m2)
    tryCatch({r2<-optim(theta00,llh,grr,theta0=theta00,y=y.pool,X=X.mat,m=m2,samp0=samp0,method="BFGS",hessian=T)}, error=function(e){status <<- 0})
    
    if(status==0){
      res<- c(status,rep(0,p^2+p*2+4))
      return(res)
    }
    
    t4<-proc.time()
    ttt<-as.vector(t4-t3)[3]
    theta<-r2$par
    hes<-r2$hessian
  }
  
  #rres<-c(status,theta00,theta,hes,c(m1,m2,t,tt,ttt))
  return(list(status=status,theta1=theta00,theta2=theta,hess=hes,clock=t+tt+ttt))
}
