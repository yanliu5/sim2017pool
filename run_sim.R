###### estimate mle based on mc method
library(lattice)
library(Matrix)
library(numDeriv)

############################
###### Main Function

J    <- 200         # number of pools
cj   <- rep(2,J)    # pool size
N    <- sum(cj)     # total number of subjects
beta <- c(1,-0.5,1,0.5,-0.5,0.5) # True beta value
sig  <- 0.5                          # sd of lognormal distribution
dd   <- 0.01        # cutoff of Monte Carlo error    

###### Parameter for generating covariate information
mu <- 0
sig <- 1
p1<-0.5
###### True coefficients
sim_dat <- data_sim_para(weight=FALSE,beta,sig,N,cj,J)

Y_pool  <- sim_dat$Y_pool
Wmat    <- sim_dat$Wmat
Imat    <- sim_dat$Imat
Zmat    <- sim_dat$Zmat
Wvec    <- sim_dat$Wvec

###################################################################################
# Wmat: Design matrix, including intercept as the column 1
# Imat: a matrix of cj*J rows; 
#   col1=0, col2=pool id what ith individual was assigned to
#   col3=weight of ith individual in the pool 
# Zmat: a matrix of J rows
#   col1=Pooled biomarker level
#   col2=cj (pool size)
#   col3-col(2+cj)=indices of individuals assigned to jth pool
###################################################################################

######
m.em<-50
cal<-500
m<-2000
mm<-50000
######
tot<-500  # how many dataset to generate for simulation
for (rep in seq(1,tot)){
  # Monte Carlo method
  
  X.pool <- Wmat[,-1]
  X.mat  <- Wmat
  Y.pool <- Y_pool
  kmc<-list(cluster=rep(1:J,cj),size=cj)
  names(kmc) <- c("cluster","size")
  
  theta0 <- try(inifit(X.pool,Y.pool,kmc),silent=TRUE)

######
theta0 <- inifit(X.pool,Y.pool,kmc)
mc_res <- mcfit(m,mm,theta0,X.mat,Y.pool,kmc)
em_res <- emfit(m.em,cal,theta0,X.mat,Y.pool,kmc)

fit<-lm(log(y.obs)~X.pool)

write(c(theta0,res1),"emfit.txt",append=TRUE)
write(c(theta0,res2),"mcfit.txt",append=TRUE)
write(c(summary(fit)$coefficients[,1:2],summary(fit)$sigma),"indfit.txt",append=TRUE)
}
