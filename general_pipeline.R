################################################################################
## R script to draw plots for simulation on parametric vs emulation           ##
##  algorithms.                                                               ##
################################################################################
##
## James Liley
## 1 December 2021
##


################################################################################
## Setup                                                                      ##
################################################################################

# Set random seed
set.seed(21423)

# Save plot to file, or not
save_plot=TRUE
pdf_path="./"

# Load package
library(OptHoldoutSize)

# Suppose we have population size and cost-per-sample without a risk score as follows
N=100000
k1=0.4

# Suppose that true values of a,b,c are given by
theta_true=c(10000,1.2,0.2)
theta_lower=c(1,0.5,0.1) # lower bounds for estimating theta
theta_upper=c(20000,2,0.5) # upper bounds for estimating theta
theta_init=(theta_lower+theta_upper)/2 # We will start from this value when finding theta

# Kernel width and variance for Gaussian process
kw0=5000
vu0=1e7

# We will presume that these are the values of n for which cost can potentially be evaluated.
n=seq(1000,N,length=300)

# Sample variances var_w of d(n) uniformly between these values
vwmin=0.001; vwmax=0.02

## True form of k2, option 1: parametric assumptions satisfied (true)
true_k2_pTRUE=function(n) powerlaw(n,theta_true)

## True form of k2, option 2: parametric assumptions NOT satisfied (false)
true_k2_pFALSE=function(n) powerlaw(n,theta_true) + (1e4)*dnorm(n,mean=4e4,sd=8e3)



################################################################################
## Plot param./nonparam forms of k2                                           ##
################################################################################

## Plot
if (save_pdf) pdf(paste0(pdf_path,"k2_params.pdf"),width=4,height=4)
par(mar=c(6,4,1,1)) # for size consistency with later plots
plot(0,type="n",xlim=range(n),ylim=c(0,1.2),
  xlab="n",ylab=expression(paste("k"[2],"(n)")))

lines(n,true_k2_pTRUE(n),type="l",col="black")
lines(n,true_k2_pFALSE(n),type="l",col="red")
legend("topright",c("Param. satis.", "Param. not satis."),col=c("black","red"),
       lty=1,bty="n")
if (save_pdf) dev.off()

if (save_pdf) pdf(paste0(pdf_path,"cost_params.pdf"),width=4,height=4)
par(mar=c(6,4,1,1)) # for size consistency with later plots
plot(0,type="n",xlim=range(n),ylim=c(20000,80000),
     xlab="n",ylab="Total cost")

cost_TRUE=k1*n + true_k2_pTRUE(n)*(N-n)
cost_FALSE=k1*n + true_k2_pFALSE(n)*(N-n)
lines(n,cost_TRUE,type="l",col="black")
lines(n,cost_FALSE,type="l",col="red")

wTRUE=which.min(cost_TRUE)
wFALSE=which.min(cost_FALSE)
points(n[wTRUE],cost_TRUE[wTRUE],col="black",pch=16)
points(n[wFALSE],cost_FALSE[wFALSE],col="red",pch=16)

legend("topright",c("Param. satis.", "Param. not satis."),
       col=c("black","red"),lty=c(1,1),pch=c(16,16),bty="n")
if (save_pdf) dev.off()


################################################################################
## True optimal holdout sizes                                                 ##
################################################################################

nsamp=200 # Presume we have this many estimates of k_2(n), between 1000 and N

nset=round(runif(nsamp,1000,N))
var_w=runif(nsamp,vwmin,vwmax)
d_pTRUE=rnorm(nsamp,mean=true_k2_pTRUE(nset),sd=sqrt(var_w))
d_pFALSE=rnorm(nsamp,mean=true_k2_pFALSE(nset),sd=sqrt(var_w))


nc=1000:N

true_ohs_pTRUE=nc[which.min(k1*nc + true_k2_pTRUE(nc)*(N-nc))]
true_ohs_pFALSE=nc[which.min(k1*nc + true_k2_pFALSE(nc)*(N-nc))]

print(true_ohs_pTRUE)

print(true_ohs_pFALSE)

print("Line 62")


# Estimate a,b, and c from values nset and d
est_abc_pTRUE=powersolve(nset,d_pTRUE,
  lower=theta_lower,upper=theta_upper,init=theta_init)$par
est_abc_pFALSE=powersolve(nset,d_pFALSE,
  lower=theta_lower,upper=theta_upper,init=theta_init)$par

# Estimate optimal holdout sizes using parametric method
param_ohs_pTRUE=optimal_holdout_size(N,k1,theta=est_abc_pTRUE)$size
param_ohs_pFALSE=optimal_holdout_size(N,k1,theta=est_abc_pFALSE)$size

# Estimate optimal holdout sizes using semi-parametric (emulation) method
emul_ohs_pTRUE=optimal_holdout_size_emulation(nset,d_pTRUE,var_w,theta=est_abc_pTRUE,N=N,k1=k1)$size
emul_ohs_pFALSE=optimal_holdout_size_emulation(nset,d_pFALSE,var_w,theta=est_abc_pFALSE,N=N,k1=k1)$size

# In the parametrised case, the parametric model is better
print(true_ohs_pTRUE)
print(param_ohs_pTRUE)
print(emul_ohs_pTRUE)

# In the misparametrised case, the semi-parametric model is better
print(true_ohs_pFALSE)
print(param_ohs_pFALSE)
print(emul_ohs_pFALSE)





################################################################################
## Resample values of d and recalculate OHS                                   ##
################################################################################

n_var=1000 # Resample d this many times to estimate OHS variance

ohs_resample=matrix(NA,n_var,4) # This will be populated with OHS estimates

for (i in 1:n_var) {
  set.seed(36253 + i)

  # Resample values d
  d_pTRUE_r=rnorm(nsamp,mean=true_k2_pTRUE(nset),sd=sqrt(var_w))
  d_pFALSE_r=rnorm(nsamp,mean=true_k2_pFALSE(nset),sd=sqrt(var_w))


  # Estimate a,b, and c from values nset and d
  est_abc_pTRUE_r=powersolve(nset,d_pTRUE_r,y_var=var_w,
    lower=theta_lower,upper=theta_upper,init=theta_init)$par
  est_abc_pFALSE_r=powersolve(nset,d_pFALSE_r,y_var=var_w,
    lower=theta_lower,upper=theta_upper,init=theta_init)$par

  # Estimate optimal holdout sizes using parametric method
  param_ohs_pTRUE_r=optimal_holdout_size(N,k1,theta=est_abc_pTRUE_r)$size
  param_ohs_pFALSE_r=optimal_holdout_size(N,k1,theta=est_abc_pFALSE_r)$size

  # Estimate optimal holdout sizes using semi-parametric (emulation) method
  emul_ohs_pTRUE_r=optimal_holdout_size_emulation(nset,d_pTRUE_r,theta=est_abc_pTRUE_r,var_w,N,k1)$size
  emul_ohs_pFALSE_r=optimal_holdout_size_emulation(nset,d_pFALSE_r,theta=est_abc_pTRUE_r,var_w,N,k1)$size

  ohs_resample[i,]=c(param_ohs_pTRUE_r, param_ohs_pFALSE_r, emul_ohs_pTRUE_r, emul_ohs_pFALSE_r)

  print(i)
}

colnames(ohs_resample)=c("param_pTRUE","param_pFALSE","emul_pTRUE", "emul_pFALSE")
save(ohs_resample,file="data/ohs_resample.RData")




################################################################################
## Draw figure showing estimation of OHS                                      ##
################################################################################

d_pt=density(ohs_resample[,"param_pTRUE"])
d_et=density(ohs_resample[,"emul_pTRUE"])
d_pf=density(ohs_resample[,"param_pFALSE"])
d_ef=density(ohs_resample[,"emul_pFALSE"])

if (save_pdf) pdf(paste0(pdf_path,"ohs_distributions.pdf"),width=4,height=4)
par(mar=c(6,4,1,1))
plot(0,type="n",xlim=c(0,5),ylim=c(8000,45000),xaxt="n",ylab="OHS",xlab="")
axis(1,at=c(1.5,3.5),label=c("Par. satis.", "Par. unsatis."),las=2)

sc=1000; hsc=0.8

lines(c(0.5,2.5),rep(true_ohs_pTRUE,2),col="gray")
lines(c(2.5,4.5),rep(true_ohs_pFALSE,2),col="gray")

polygon(1+sc*c(d_pt$y,-rev(d_pt$y)),c(d_pt$x,rev(d_pt$x)),col=rgb(0,0,0,alpha=0.5),border=NA)
polygon(2+sc*c(d_et$y,-rev(d_et$y)),c(d_et$x,rev(d_et$x)),col=rgb(1,0,0,alpha=0.5),border=NA)
polygon(3+sc*c(d_pf$y,-rev(d_pf$y)),c(d_pf$x,rev(d_pf$x)),col=rgb(0,0,0,alpha=0.5),border=NA)
polygon(4+sc*c(d_ef$y,-rev(d_ef$y)),c(d_ef$x,rev(d_ef$x)),col=rgb(1,0,0,alpha=0.5),border=NA)

lines(1+hsc*c(-0.5,0.5),rep(median(ohs_resample[,"param_pTRUE"]),2),col="black",lty=2,lwd=2)
lines(2+hsc*c(-0.5,0.5),rep(median(ohs_resample[,"emul_pTRUE"]),2),col="red",lty=2,lwd=2)
lines(3+hsc*c(-0.5,0.5),rep(median(ohs_resample[,"param_pFALSE"]),2),col="black",lty=2,lwd=2)
lines(4+hsc*c(-0.5,0.5),rep(median(ohs_resample[,"emul_pFALSE"]),2),col="red",lty=2,lwd=2)


#legend("bottomleft",bty="n",
#       c("Dens. param.","Dens. emul","Med. param","Med. emul","True OHS"),
#       lty=c(NA,NA,2,2,1),lwd=c(NA,NA,2,2,1),
#       pch=c(16,16,NA,NA,NA),col=c(rgb(0,0,0,alpha=0.5),rgb(1,0,0,alpha=0.5),"black","red","gray"))

legend("bottomleft",bty="n",
       c("Param.","Emul","True OHS"),
       lty=c(2,2,1),lwd=c(2,2,1),
       pch=c(16,16,NA),col=c("black","red","gray"))


if (save_pdf) dev.off()




################################################################################
## Next point using parametrisation                                           ##
################################################################################

## Choose an initial five training sizes at which to evaluate k2
set.seed(32424) # start from same seed as before
nstart=5
nset0=round(runif(nstart,1000,N/2))
var_w0=runif(nstart,vwmin,vwmax)
d0_pTRUE=rnorm(nstart,mean=true_k2_pTRUE(nset0),sd=sqrt(var_w0))
d0_pFALSE=rnorm(nstart,mean=true_k2_pFALSE(nset0),sd=sqrt(var_w0))


# These are our sets of training sizes and k2 estimates, which will be built up.
nset_pTRUE=nset0
d_pTRUE=d0_pTRUE
var_w_pTRUE=var_w0

nset_pFALSE=nset0
d_pFALSE=d0_pFALSE
var_w_pFALSE=var_w0


# Go up to this many points
max_points=200

while(length(nset_pTRUE)<= max_points) {
  set.seed(46352 + length(nset_pTRUE))

  # Estimate parameters for parametric part of semi-parametric method
  theta_pTRUE=powersolve(nset_pTRUE,d_pTRUE,y_var=var_w_pTRUE,
    lower=theta_lower,upper=theta_upper,init=theta_init)$par
  theta_pFALSE=powersolve(nset_pTRUE,d_pFALSE,y_var=var_w_pTRUE,
    lower=theta_lower,upper=theta_upper,init=theta_init)$par

  # Mean and variance of emulator for cost function, parametric assumptions satisfied
  p_mu_pTRUE=mu_fn(n,nset=nset_pTRUE,d=d_pTRUE,var_w = var_w_pTRUE,theta=theta_pTRUE,N=N,k1=k1)
  p_var_pTRUE=psi_fn(n,nset=nset_pTRUE,N=N,var_w = var_w_pTRUE)

  # Mean and variance of emulator for cost function, parametric assumptions not satisfied
  p_mu_pFALSE=mu_fn(n,nset=nset_pFALSE,d=d_pFALSE,var_w = var_w_pFALSE,theta=theta_pFALSE,N=N,k1=k1)
  p_var_pFALSE=psi_fn(n,nset=nset_pFALSE,N=N,var_w = var_w_pFALSE)

  # Add vertical line at next suggested point
  exp_imp_em_pTRUE = exp_imp_fn(n,nset=nset_pTRUE,d=d_pTRUE,var_w = var_w_pTRUE, N=N,k1=k1)
  nextn_pTRUE = n[which.max(exp_imp_em_pTRUE)]

  # Find next suggested point, parametric assumptions not satisfied
  exp_imp_em_pFALSE = exp_imp_fn(n,nset=nset_pFALSE,d=d_pFALSE,var_w = var_w_pFALSE, N=N,k1=k1)
  nextn_pFALSE = n[which.max(exp_imp_em_pFALSE)]


  # New estimates of k2
  var_w_new_pTRUE=runif(1,vwmin,vwmax)
  d_new_pTRUE=rnorm(1,mean=true_k2_pTRUE(nextn_pTRUE),sd=sqrt(var_w_new_pTRUE))

  var_w_new_pFALSE=runif(1,vwmin,vwmax)
  d_new_pFALSE=rnorm(1,mean=true_k2_pFALSE(nextn_pFALSE),sd=sqrt(var_w_new_pFALSE))


  # Update data
  nset_pTRUE=c(nset_pTRUE,nextn_pTRUE)
  d_pTRUE=c(d_pTRUE,d_new_pTRUE)
  var_w_pTRUE=c(var_w_pTRUE,var_w_new_pTRUE)

  nset_pFALSE=c(nset_pFALSE,nextn_pFALSE)
  d_pFALSE=c(d_pFALSE,d_new_pFALSE)
  var_w_pFALSE=c(var_w_pFALSE,var_w_new_pFALSE)

  print(length(nset_pFALSE))


}

data_nextpoint_em=list(
  nset_pTRUE=nset_pTRUE,nset_pFALSE=nset_pFALSE,
  d_pTRUE=d_pTRUE,d_pFALSE=d_pFALSE,
  var_w_pTRUE=var_w_pTRUE,var_w_pFALSE=var_w_pFALSE)

save(data_nextpoint_em,file="data/data_nextpoint_em.RData")



## To draw plot with np points (np can be set using the button)

np=50 # or set using interactive session

par(mfrow=c(1,2))
yrange=c(0,100000)

# Mean and variance of emulator for cost function, parametric assumptions satisfied
p_mu_pTRUE=mu_fn(n,nset=nset_pTRUE[1:np],d=d_pTRUE[1:np],var_w = var_w_pTRUE[1:np],N=N,k1=k1)
p_var_pTRUE=psi_fn(n,nset=nset_pTRUE[1:np],N=N,var_w = var_w_pTRUE[1:np])

# Mean and variance of emulator for cost function, parametric assumptions not satisfied
p_mu_pFALSE=mu_fn(n,nset=nset_pFALSE[1:np],d=d_pFALSE[1:np],var_w = var_w_pFALSE[1:np],N=N,k1=k1)
p_var_pFALSE=psi_fn(n,nset=nset_pFALSE[1:np],N=N,var_w = var_w_pFALSE[1:np])


# Estimate parameters for parametric part of semi-parametric method
theta_pTRUE=powersolve(nset_pTRUE[1:np],d_pTRUE[1:np],y_var=var_w_pTRUE[1:np],lower=theta_lower,upper=theta_upper,init=theta_init)$par
theta_pFALSE=powersolve(nset_pFALSE[1:np],d_pFALSE[1:np],y_var=var_w_pFALSE[1:np],lower=theta_lower,upper=theta_upper,init=theta_init)$par

## First panel
plot(0,xlim=range(n),ylim=yrange,type="n",
  xlab="Training/holdout set size",
  ylab="Total cost (= num. cases)")
lines(n,p_mu_pTRUE,col="blue")
lines(n,p_mu_pTRUE - 3*sqrt(pmax(0,p_var_pTRUE)),col="red")
lines(n,p_mu_pTRUE + 3*sqrt(pmax(0,p_var_pTRUE)),col="red")
points(nset_pTRUE[1:np],k1*nset_pTRUE[1:np] + d_pTRUE[1:np]*(N-nset_pTRUE[1:np]),pch=16,cex=1,col="purple")
lines(n,k1*n + powerlaw(n,theta_pTRUE)*(N-n),lty=2)
lines(n,k1*n + true_k2_pTRUE(n)*(N-n),lty=3,lwd=3)
legend("topright",
  c(expression(mu(n)),
    expression(mu(n) %+-% 3*sqrt(psi(n))),
    "Par. est. cost",
    "True",
    "d",
    "Next pt."),
  lty=c(1,1,2,3,NA,NA),lwd=c(1,1,1,3,NA,NA),pch=c(NA,NA,NA,NA,16,124),pt.cex=c(NA,NA,NA,NA,1,1),
  col=c("blue","red","black","black","purple","black"),bg="white",border=NA)

abline(v=nset_pTRUE[np+1])



## Second panel
plot(0,xlim=range(n),ylim=yrange,type="n",
  xlab="Training/holdout set size",
  ylab="Total cost (= num. cases)")
lines(n,p_mu_pFALSE,col="blue")
lines(n,p_mu_pFALSE - 3*sqrt(pmax(0,p_var_pFALSE)),col="red")
lines(n,p_mu_pFALSE + 3*sqrt(pmax(0,p_var_pFALSE)),col="red")
points(nset_pFALSE[1:np],k1*nset_pFALSE[1:np] + d_pFALSE[1:np]*(N-nset_pFALSE[1:np]),pch=16,cex=1,col="purple")
lines(n,k1*n + powerlaw(n,theta_pFALSE)*(N-n),lty=2)
lines(n,k1*n + true_k2_pFALSE(n)*(N-n),lty=3,lwd=3)
legend("topright",
  c(expression(mu(n)),
    expression(mu(n) %+-% 3*sqrt(psi(n))),
    "Par. est. cost",
    "True",
    "d",
    "Next pt."),
  lty=c(1,1,2,3,NA,NA),lwd=c(1,1,1,3,NA,NA),pch=c(NA,NA,NA,NA,16,124),pt.cex=c(NA,NA,NA,NA,1,1),
  col=c("blue","red","black","black","purple","black"),bg="white",border=NA)

abline(v=nset_pFALSE[np+1])







################################################################################
## Next point using Bayesian emulation                                        ##
################################################################################

## Choose an initial five training sizes at which to evaluate k2
set.seed(32424) # start from same seed as before
nstart=5
nset0=round(runif(nstart,1000,N/2))
var_w0=runif(nstart,vwmin,vwmax)
d0_pTRUE=rnorm(nstart,mean=true_k2_pTRUE(nset0),sd=sqrt(var_w0))
d0_pFALSE=rnorm(nstart,mean=true_k2_pFALSE(nset0),sd=sqrt(var_w0))


# These are our sets of training sizes and k2 estimates, which will be built up.
nset_pTRUE=nset0
d_pTRUE=d0_pTRUE
var_w_pTRUE=var_w0

nset_pFALSE=nset0
d_pFALSE=d0_pFALSE
var_w_pFALSE=var_w0


# Go up to this many points
max_points=200

while(length(nset_pTRUE)<= max_points) {
  set.seed(46352 + length(nset_pTRUE))

  # Estimate parameters for parametric part of semi-parametric method
  theta_pTRUE=powersolve(nset_pTRUE,d_pTRUE,y_var=var_w_pTRUE,lower=theta_lower,upper=theta_upper,init=theta_init)$par
  theta_pFALSE=powersolve(nset_pTRUE,d_pFALSE,y_var=var_w_pTRUE,lower=theta_lower,upper=theta_upper,init=theta_init)$par

  # Mean and variance of emulator for cost function, parametric assumptions satisfied
  p_mu_pTRUE=mu_fn(n,nset=nset_pTRUE,d=d_pTRUE,var_w = var_w_pTRUE,theta=theta_pTRUE,N=N,k1=k1)
  p_var_pTRUE=psi_fn(n,nset=nset_pTRUE,N=N,var_w = var_w_pTRUE)

  # Mean and variance of emulator for cost function, parametric assumptions not satisfied
  p_mu_pFALSE=mu_fn(n,nset=nset_pFALSE,d=d_pFALSE,var_w = var_w_pFALSE,theta=theta_pFALSE,N=N,k1=k1)
  p_var_pFALSE=psi_fn(n,nset=nset_pFALSE,N=N,var_w = var_w_pFALSE)



  # Two panels
  par(mfrow=c(1,2))
  yrange=c(0,100000)


  ## First panel
  plot(0,xlim=range(n),ylim=yrange,type="n",
    xlab="Training/holdout set size",
    ylab="Total cost (= num. cases)")
  lines(n,p_mu_pTRUE,col="blue")
  lines(n,p_mu_pTRUE - 3*sqrt(pmax(0,p_var_pTRUE)),col="red")
  lines(n,p_mu_pTRUE + 3*sqrt(pmax(0,p_var_pTRUE)),col="red")
  points(nset_pTRUE,k1*nset_pTRUE + d_pTRUE*(N-nset_pTRUE),pch=16,col="purple")
  lines(n,k1*n + powerlaw(n,theta_pTRUE)*(N-n),lty=2)
  lines(n,k1*n + true_k2_pTRUE(n)*(N-n),lty=3,lwd=3)
  legend("topright",
    c(expression(mu(n)),
      expression(mu(n) %+-% 3*sqrt(psi(n))),
      "prior(n)",
      "True",
      "d",
      "Next pt."),
    lty=c(1,1,2,3,NA,NA),lwd=c(1,1,1,3,NA,NA),pch=c(NA,NA,NA,NA,16,124),pt.cex=c(NA,NA,NA,NA,1,1),
    col=c("blue","red","black","black","purple","black"),bg="white",border=NA)

  # Add vertical line at next suggested point
  exp_imp_em_pTRUE = exp_imp_fn(n,nset=nset_pTRUE,d=d_pTRUE,var_w = var_w_pTRUE, N=N,k1=k1)
  nextn_pTRUE = n[which.max(exp_imp_em_pTRUE)]
  abline(v=nextn_pTRUE)


  ## Second panel
  plot(0,xlim=range(n),ylim=yrange,type="n",
    xlab="Training/holdout set size",
    ylab="Total cost (= num. cases)")
  lines(n,p_mu_pFALSE,col="blue")
  lines(n,p_mu_pFALSE - 3*sqrt(pmax(0,p_var_pFALSE)),col="red")
  lines(n,p_mu_pFALSE + 3*sqrt(pmax(0,p_var_pFALSE)),col="red")
  points(nset_pFALSE,k1*nset_pFALSE + d_pFALSE*(N-nset_pFALSE),pch=16,col="purple")
  lines(n,k1*n + powerlaw(n,theta_pFALSE)*(N-n),lty=2)
  lines(n,k1*n + true_k2_pFALSE(n)*(N-n),lty=3,lwd=3)
  legend("topright",
    c(expression(mu(n)),
      expression(mu(n) %+-% 3*sqrt(psi(n))),
      "prior(n)",
      "True",
      "d",
      "Next pt."),
    lty=c(1,1,2,3,NA,NA),lwd=c(1,1,1,3,NA,NA),pch=c(NA,NA,NA,NA,16,124),pt.cex=c(NA,NA,NA,NA,1,1),
    col=c("blue","red","black","black","purple","black"),bg="white",border=NA)

  # Add vertical line at next suggested point
  exp_imp_em_pFALSE = exp_imp_fn(n,nset=nset_pFALSE,d=d_pFALSE,var_w = var_w_pFALSE, N=N,k1=k1)
  nextn_pFALSE = n[which.max(exp_imp_em_pFALSE)]
  abline(v=nextn_pFALSE)


  # New estimates of k2
  var_w_new_pTRUE=runif(1,vwmin,vwmax)
  d_new_pTRUE=rnorm(1,mean=true_k2_pTRUE(nextn_pTRUE),sd=sqrt(var_w_new_pTRUE))

  var_w_new_pFALSE=runif(1,vwmin,vwmax)
  d_new_pFALSE=rnorm(1,mean=true_k2_pFALSE(nextn_pFALSE),sd=sqrt(var_w_new_pFALSE))


  # Update data
  nset_pTRUE=c(nset_pTRUE,nextn_pTRUE)
  d_pTRUE=c(d_pTRUE,d_new_pTRUE)
  var_w_pTRUE=c(var_w_pTRUE,var_w_new_pTRUE)

  nset_pFALSE=c(nset_pFALSE,nextn_pFALSE)
  d_pFALSE=c(d_pFALSE,d_new_pFALSE)
  var_w_pFALSE=c(var_w_pFALSE,var_w_new_pFALSE)

  print(length(nset_pFALSE))


}

data_nextpoint_em=list(
  nset_pTRUE=nset_pTRUE,nset_pFALSE=nset_pFALSE,
  d_pTRUE=d_pTRUE,d_pFALSE=d_pFALSE,
  var_w_pTRUE=var_w_pTRUE,var_w_pFALSE=var_w_pFALSE)

save(data_nextpoint_em,file="data/data_nextpoint_em.RData")



################################################################################
## Analyse convergence rate of each algorithm                                 ##
################################################################################

# Function to resample values of d and regenerate OHS given nset and var_w
ntri=function(nset,var_w,k2,nx=100,method="MLE") {
  out=rep(0,nx)
  for (i in 1:nx) {
    d1=rnorm(length(nset),mean=k2(nset),sd=sqrt(var_w))
    theta1=powersolve(nset,d1,y_var=var_w,lower=theta_lower,upper=theta_upper,init=theta_true)$par
    if (method=="MLE") {
      out[i]=optimal_holdout_size(N,k1,theta1)$size
    } else {
      nn=seq(1000,N,length=1000)
      p_mu=mu_fn(nn,nset=nset,d=d1,var_w = var_w, N=N,k1=k1,theta=theta1)
      out[i]=nn[which.min(p_mu)]
    }
  }
  return(out)
}


# Load datasets of 'next point'
load("data/data_nextpoint_em.RData")
load("data/data_nextpoint_par.RData")

# Maximum number of training set sizes we will consider
n_iter=200

# Generate random 'next points'
set.seed(36279)
data_nextpoint_rand=list(
  nset_pTRUE=round(runif(n_iter,1000,N)),
  nset_pFALSE=round(runif(n_iter,1000,N)),
  var_w_pTRUE=runif(n_iter,vwmin,vwmax),
  var_w_pFALSE=runif(n_iter,vwmin,vwmax)
)

# Initialise matrices of records
# ohs_array[n,i,j,k,l] is
#  using the first n training set sizes
#  the ith resample
#  using the jth version of k2 (j=1: pTRUE, j=2: pFALSE)
#  using the kth algorithm (k=1: parametric, k=2: semiparametric/emulation)
#  using the lth method of selecting next points (l=1: random, l=2: systematic)
nr=200 # recalculate OHS this many times
ohs_array=array(NA,dim=c(n_iter,nr,2,2,2))

for (i in 5:n_iter) {
  set.seed(363726 + i)
  # Resamplings for parametric algorithm, random next point
  ohs_array[i,,1,1,1]=ntri(
    nset=data_nextpoint_rand$nset_pTRUE[1:i],
    var_w=data_nextpoint_rand$var_w_pTRUE[1:i],
    k2=true_k2_pTRUE,nx=nr,method="MLE")
  ohs_array[i,,2,1,1]=ntri(
    nset=data_nextpoint_rand$nset_pFALSE[1:i],
    var_w=data_nextpoint_rand$var_w_pFALSE[1:i],
    k2=true_k2_pFALSE,nx=nr,method="MLE")

  # Resamplings for semiparametric/emulation algorithm, random next point
  ohs_array[i,,1,2,1]=ntri(
    nset=data_nextpoint_rand$nset_pTRUE[1:i],
    var_w=data_nextpoint_rand$var_w_pTRUE[1:i],
    k2=true_k2_pTRUE,nx=nr,method="EM")
  ohs_array[i,,2,2,1]=ntri(
    nset=data_nextpoint_rand$nset_pFALSE[1:i],
    var_w=data_nextpoint_rand$var_w_pFALSE[1:i],
    k2=true_k2_pFALSE,nx=nr,method="EM")

  # Resamplings for parametric algorithm, nonrandom (systematic) next point
  ohs_array[i,,1,1,2]=ntri(
    nset=data_nextpoint_par$nset_pTRUE[1:i],
    var_w=data_nextpoint_par$var_w_pTRUE[1:i],
    k2=true_k2_pTRUE,nx=nr,method="MLE")
  ohs_array[i,,2,1,2]=ntri(
    nset=data_nextpoint_par$nset_pFALSE[1:i],
    var_w=data_nextpoint_par$var_w_pFALSE[1:i],
    k2=true_k2_pFALSE,nx=nr,method="MLE")

  # Resamplings for semiparametric/emulation algorithm, nonrandom (systematic) next point
  ohs_array[i,,1,2,2]=ntri(
    nset=data_nextpoint_em$nset_pTRUE[1:i],
    var_w=data_nextpoint_em$var_w_pTRUE[1:i],
    k2=true_k2_pTRUE,nx=nr,method="EM")
  ohs_array[i,,2,2,2]=ntri(
    nset=data_nextpoint_em$nset_pFALSE[1:i],
    var_w=data_nextpoint_em$var_w_pFALSE[1:i],
    k2=true_k2_pFALSE,nx=nr,method="EM")

  print(i)

  save(ohs_array,file="data/ohs_array.RData")
}







# Load data
data(ohs_array)
# ohs_array[n,i,j,k,l] is
#  using the first n training set sizes
#  the ith resample
#  using the jth version of k2 (j=1: pTRUE, j=2: pFALSE)
#  using the kth algorithm (k=1: parametric, k=2: semiparametric/emulation)
#  using the lth method of selecting next points (l=1: random, l=2: systematic)

# Settings
alpha=0.5; # Plot alpha/2,1-alpha/2 quantiles
dd=3 # horizontal line spacing
n_iter=dim(ohs_array)[1] # X axis range
ymax=80000 # Y axis range

# Plot drawing function
plot_ci_convergence=function(title,key,M1,M2,ohs_true) {

  # Set up plot parameters
  par(mar=c(1,4,4,0.1))
  layout(mat=rbind(matrix(1,4,4),matrix(2,2,4)))

  # Number of estimates
  n_iterx=dim(M1)[1]

  # Initialise
  plot(0,xlim=c(5,n_iterx),ylim=c(0,ymax),type="n",
    ylab="OHS and error",xaxt="n",main=title)
  abline(h=ohs_true,col="blue",lty=2)

  # Plot medians
  points(1:n_iterx,rowMedians(M1,na.rm=T),pch=16,cex=0.5,col="black")
  points(1:n_iterx,rowMedians(M2,na.rm=T),pch=16,cex=0.5,col="red")

  # CIs
  ci1=rbind(apply(M1,1,function(x) pmax(0,quantile(x,alpha/2,na.rm=T))),
    apply(M1,1,function(x) pmin(ymax,quantile(x,1-alpha/2,na.rm=T))))
  ci2=rbind(apply(M2,1,function(x) pmax(0,quantile(x,alpha/2,na.rm=T))),
    apply(M2,1,function(x) pmin(ymax,quantile(x,1-alpha/2,na.rm=T))))

  # Plot CIs
  segments(
    1:n_iterx,ci1[1,],
    1:n_iterx,ci1[2,],
    col="black"
  )
  segments(
    1:n_iterx + 1/dd,ci2[1,],
    1:n_iterx + 1/dd,ci2[2,],
    col="red"
  )

  # Add legend
  legend("topright",
    legend=key,
    col=c("black","red"),lty=1)


  # Bottom panel setup
  par(mar=c(4,4,0.1,0.1))
  plot(0,xlim=c(5,n_iterx),
       ylim=c(0,quantile(c(ci1[2,]-ci1[1,],ci2[2,]-ci2[1,]),0.95,na.rm=T)),
    type="n",ylab="IQR",xlab=expression(paste("Number of estimates of k"[2],"(n)")))

  # Draw lines
  lines(1:n_iterx,ci1[2,]-ci1[1,],col="black")
  lines(1:n_iterx,ci2[2,]-ci2[1,],col="red")

}



# Extract matrices from aray
M111=ohs_array[1:n_iter,,1,1,1] # pTRUE, param algorithm, random nextpoint
M112=ohs_array[1:n_iter,,1,1,2] # pTRUE, param algorithm, systematic nextpoint

M211=ohs_array[1:n_iter,,2,1,1] # pFALSE, param algorithm, random nextpoint
M212=ohs_array[1:n_iter,,2,1,2] # pFALSE, param algorithm, systematic nextpoint

M121=ohs_array[1:n_iter,,1,2,1] # pTRUE, emul algorithm, random nextpoint
M122=ohs_array[1:n_iter,,1,2,2] # pTRUE, emul algorithm, systematic nextpoint

M221=ohs_array[1:n_iter,,2,2,1] # pFALSE, emul algorithm, random nextpoint
M222=ohs_array[1:n_iter,,2,2,2] # pFALSE, emul algorithm, systematic nextpoint


# True OHS
nc=1000:N
nmax=50 # Plot up to this number of estimates.
true_ohs_pTRUE=nc[which.min(k1*nc + true_k2_pTRUE(nc)*(N-nc))]
true_ohs_pFALSE=nc[which.min(k1*nc + true_k2_pFALSE(nc)*(N-nc))]

if (save_pdf) pdf(paste0(pdf_path,"conv_comp_11.pdf"),width=4,height=4)
plot_ci_convergence("Params. satis, param. alg.",
  c("Rand. next n","Syst. next n"),M111[1:nmax,],M112[1:nmax,],true_ohs_pTRUE)
if (save_pdf) dev.off()

if (save_pdf) pdf(paste0(pdf_path,"conv_comp_12.pdf"),width=4,height=4)
plot_ci_convergence("Params. not satis, param. alg.",
  c("Rand. next n","Syst. next n"),M211[1:nmax,],M212[1:nmax,],true_ohs_pFALSE)
if (save_pdf) dev.off()


if (save_pdf) pdf(paste0(pdf_path,"conv_comp_21.pdf"),width=4,height=4)
plot_ci_convergence("Params. satis, emul. alg.",
  c("Rand. next n","Syst. next n"),M121[1:nmax,],M122[1:nmax,],true_ohs_pTRUE)
if (save_pdf) dev.off()

if (save_pdf) pdf(paste0(pdf_path,"conv_comp_22.pdf"),width=4,height=4)
plot_ci_convergence("Params. not satis, emul. alg.",
  c("Rand. next n","Syst. next n"),M221[1:nmax,],M222[1:nmax,],true_ohs_pFALSE)
if (save_pdf) dev.off()

