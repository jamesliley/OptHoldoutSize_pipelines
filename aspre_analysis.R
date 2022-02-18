################################################################################
## R script to estimate parameters of learning curve for ASPRE model          ##
################################################################################
##
## James Liley
## 21 October 2021
##
## Note: we do not have access to the script used to fit the model. THIS IS A
##  SIMULATION ONLY, INTENDED ONLY TO DEMONSTRATE HOW THE ALGORITHM COULD BE
##  APPLIED. DO NOT INTERPRET ANY RESULTS AS CLINICALLY USABLE.
##
## This simulation works by taking the final ASPRE model (Akolekar et al, 2013)
##  with >55,000 samples and considering it to have optimal performance.

## This script should be ran in the directory in which it is saved, or in 
##  some directory with subdirectories 'data', 'figures'.
## Not all figures are necessarily used in the manuscript.

######################################################################
## Scripts, switches and libraries                                  ##
######################################################################

# Set random seed
seed=463825
set.seed(seed)

# Libraries
library(mvtnorm)
library(matrixStats)
library(mle.tools)
library(OptHoldoutSize)

# Save plots to file, or not
save_plot=TRUE

# Force redo: set to TRUE to regenerate all datasets from scratch
force_redo=FALSE

# Include legends on plots for main paper
inc_legend=FALSE

# PDF dimensions in inches
pdf_dim=3.5

######################################################################
## Parameters                                                       ##
######################################################################


#### ASPRE-related settings

# Total individuals in trial; all data
n_aspre_total=58794

# Population untreated PRE prevalence
pi_PRE = 1426/58974

# Maximum score sensitivity amongst highest 10%: assumed to be that of ASPRE
sens_max = (138+194)/2707  # = 0.122645 , from abstract of Rolnik 2017 Ultrasound in O&G

# Intervene with aspirin on this proportion of individuals
pi_intervention=0.1

# Aspirin reduces PRE risk by approximately this much
alpha=0.37
SE_alpha=0.09

# Candidate values for n
nval=round(seq(500,30000,length=100))

# Variance and covariance parameters
var_u=1000
k_width=5000

######################################################################
## Mockup of real data                                              ##
######################################################################

# Parameters of true ASPRE dataset
data(params_aspre)

# Simulate random dataset
X=sim_random_aspre(n_aspre_total,params=params_aspre)
X1=add_aspre_interactions(X)

# Risk will be monotonic to ASPRE risk, but we will transform to match
#  population prevalence of PE and sensitivity of ASPRE score.
risk0=aspre(X1)

# Find a linear transformation ax+b of lrisk such that population prevalence
#  and expected sensitivity match. Suppose P(Y_i=1)=score_i
# Expected sensitivity = E_{Y|scores}(sens)
#                      = (1/(pi_intervention*n_aspre_total))*E{sum_{i:score(i)>thresh} [Y_i]}
#                      = (1/5879)*sum_{i:score(i)>thresh} [(score(i)])
lrisk0=logistic(risk0)
f_ab=function(ab) {
  a=ab[1]; b=ab[2]
  risk_ab=a*lrisk0 + b
  pop_prev=mean(logit(risk_ab))
  q_pi=quantile(risk_ab,0.9)
  sens=(1/(pi_intervention*n_aspre_total))*sum(logit(risk_ab)*(risk_ab>q_pi))
  return((pop_prev-pi_PRE)^2 + (sens - sens_max)^2)
}
abmin=optim(c(1,0),f_ab)$par
lrisk=abmin[1]*lrisk0 +abmin[2]
risk=logit(lrisk)

# PRE is a 0/1 variable indicating whether that simulated patient had PRE.
PRE=rbinom(n_aspre_total,1,prob=risk) # ASPRE=ground truth



######################################################################
## Approximate ASPRE model with logistic model                      ##
######################################################################

## ASPRE model performance is close enough to well-approximated by a logistic model
if (FALSE) {
  set.seed(seed)
  train=sample(n_aspre_total,40000); test=setdiff(1:n_aspre_total,train)
  aspre_simple=glm(PRE~.,data=cbind(X,PRE)[train,],family=binomial(link="logit"))
  ytest=predict(aspre_simple,X[test,],type="response")
  pre_test=PRE[test]
  sens_max_logistic=sens10(pre_test,ytest) # close enough to ASPRE sensitivity
}



######################################################################
## Values and errors for N and k1                                   ##
######################################################################

# Parameter calculation for N
N=400000; SE_N=1500

# Parameter calculation for k1
NICE_sensitivity=0.2
pi_1=NICE_sensitivity*(239/8875)/pi_intervention
pi_0=(1-NICE_sensitivity)*(239/8875)/(1-pi_intervention)
SE_pi_1=sqrt(pi_1*(1-pi_1)/(8875*0.1))
SE_pi_0=sqrt(pi_0*(1-pi_0)/(8875*0.9))
k1=pi_0*(1-pi_intervention) + pi_1*pi_intervention*alpha

# Standard error for k1
pi_1_s=rnorm(1000,mean=pi_1,sd=SE_pi_1)
pi_0_s=rnorm(1000,mean=pi_0,sd=SE_pi_0)
alpha_s=rnorm(1000,mean=alpha,sd=SE_alpha)
SE_k1=sd(pi_0_s*(1-pi_intervention) + pi_1_s*pi_intervention*alpha_s)




######################################################################
## Parametric approach: choose a set of values of n                 ##
######################################################################

if (!file.exists("data/aspre_parametric.RData") | force_redo) {

  set.seed(487276)

  # Start with estimates of k2 at 10 values of n
  nn_par=round(runif(20,20,150)^2)
  k2_par=0*nn_par;
  for (i in 1:length(nn_par)) {
    k2_par[i]=aspre_k2(nn_par[i],X,PRE)
  }

  # Starting value for theta
  theta=powersolve_general(nn_par,k2_par)$par
  theta_se=powersolve_se(nn_par,k2_par,init=theta)

  # Rough estimate for variance of k2
  dvar0=var(k2_par-powerlaw(nn_par,theta))
  s2_par=rep(dvar0,length(k2_par))

  ## Successively add new points
  for (i in 1:100) {
    nxn=next_n(nval,nn_par,k2_par,var_k2 = s2_par,N=N,k1=k1,nmed=10)
    if (any(is.finite(nxn))) n_new=nval[which.min(nxn)] else n_new=sample(nval,1)
    k2_new=aspre_k2(n_new,X,PRE)
    nn_par=c(nn_par,n_new)
    k2_par=c(k2_par,k2_new)
    s2_par=c(s2_par,dvar0)
    print(i)
  }

  # Resample k2(n), to avoid double-dipping effect
  for (i in 1:length(nn_par)) {
    k2_par[i]=aspre_k2(nn_par[i],X,PRE)
  }

  # Transform to total cost
  cc_par=k1*nn_par + k2_par*(N-nn_par)

  # Save
  aspre_parametric=list(nn_par=nn_par,k2_par=k2_par,s2_par=s2_par,cc_par=cc_par)
  save(aspre_parametric,file="data/aspre_parametric.RData")

} else load("data/aspre_parametric.RData")

for (i in 1:length(aspre_parametric)) assign(names(aspre_parametric)[i],aspre_parametric[[i]])

######################################################################
## Estimate power law parameters and variance                       ##
######################################################################

theta=powersolve_general(nn_par,k2_par)$par
theta_se=powersolve_se(nn_par,k2_par,init=theta)


######################################################################
## Estimate optimum holdout size and variance                       ##
######################################################################


# Optimal holdout set size and cost
optim_aspre=optimal_holdout_size(N,k1,theta)
OHS_ASPRE=optim_aspre$size
Min_cost_ASPRE=optim_aspre$cost

# Errors
cov_par=matrix(0,5,5);
cov_par[1,1]=SE_N^2; cov_par[2,2]=SE_k1^2
cov_par[3:5,3:5]=theta_se
CI_OHS_ASPRE=ci_ohs(N,k1,theta,sigma=cov_par,mode = "asymptotic",grad_nstar=grad_nstar_powerlaw,alpha = 0.1)



######################################################################
## Trace of OHS at different N using parametric algorithm           ##
######################################################################

if (!file.exists("data/aspre_trace_par.RData")) {
  
  ## Store OHS at various number of samples n
  ohs_trace_par=rep(NA,length(nn_par))
  mincost_trace_par=rep(NA,length(nn_par))
  ci_trace_par=list()
  
  ## Compute, starting with 5 points
  cov_par_sub=matrix(0,5,5);
  cov_par_sub[1,1]=SE_N^2; cov_par_sub[2,2]=SE_k1^2
  
  for (i in 5:length(nn_par)) {
    nnsub=nn_par[1:i]; k2sub=k2_par[1:i]; s2sub=s2_par[1:i]
    thetasub=powersolve_general(nnsub,k2sub)$par
    theta_sesub=powersolve_se(nnsub,k2sub,init=theta)
    
    # Optimal holdout set size and cost
    optim_aspre_sub=optimal_holdout_size(N,k1,thetasub)
    ohs_trace_par[i]=optim_aspre_sub$size
    mincost_trace_par[i]=optim_aspre_sub$cost
    
    # Errors
    cov_par_sub[3:5,3:5]=theta_sesub
    ci_trace_par[[i]]=ci_ohs(N,k1,thetasub,sigma=cov_par_sub,
                             mode = "asymptotic",
                             grad_nstar=grad_nstar_powerlaw,alpha = 0.1)
  }
  
  save(ohs_trace_par,mincost_trace_par,ci_trace_par,file="data/aspre_trace_par.RData")
  
} else load("data/aspre_trace_par.RData")



######################################################################
## Draw figure for cost function                                    ##
######################################################################


if (save_plot) pdf("figures/cost_function_estimate_param.pdf",width=pdf_dim,height=pdf_dim)

plot(0,xlim=range(nn_par),ylim=range(cc_par),type="n",
     xlab="Training set size",yaxt="n",
     ylab=expression(paste("Total. cost ", "(", "","",
                           phantom() %prop% phantom(), " sens.", ")", "")))
axis(2,las=2)
points(nn_par,cc_par,pch=16,cex=0.5)
lines(nval,k1*nval + powerlaw(nval,theta)*(N-nval))
e_min=min(CI_OHS_ASPRE); e_max=max(CI_OHS_ASPRE); c_min=min(cc_par); c_max=max(cc_par);
polygon(c(e_min,e_min,e_max,e_max),c(c_min,c_max,c_max,c_min),
        col=rgb(1,0,0,alpha=0.2),border=NA)
points(OHS_ASPRE,Min_cost_ASPRE,pch=16,col="red")

if (inc_legend) legend("topright",
       c("Cost function",
         "Est cost (d)",
         "OHS",
         "OHS err."),
       lty=c(1,NA,NA,NA),lwd=c(1,NA,NA,NA),pch=c(NA,16,16,16),pt.cex=c(NA,0.5,1,2),
       col=c("black","black","red",rgb(1,0,0,alpha=0.2)),bg="white",bty="n")

if (save_plot) dev.off()





######################################################################
## Emulation approach: choose a set of values of n                  ##
######################################################################

if (!file.exists("data/aspre_emulation.RData") | force_redo) {

  # Variance and covariance parameters
  var_u=1000
  k_width=5000

  # Begin as for parametric approach
  set.seed(487276)

  # Start with estimates of k2 at 10 values of n
  nn_emul=round(runif(20,20,150)^2)
  k2_emul=0*nn_emul;
  for (i in 1:length(nn_emul)) {
    k2_emul[i]=aspre_k2(nn_emul[i],X,PRE)
  }

  # Candidate values for n
  nval=round(seq(500,30000,length=100))

  # Starting value for theta
  theta=powersolve_general(nn_emul,k2_emul)$par

  # Rough estimate for variance of k2
  dvar0=var(k2_emul-powerlaw(nn_emul,theta))
  s2_emul=rep(dvar0,length(k2_emul))


  ## Successively add new points
  for (i in 1:100) {
    nxn = exp_imp_fn(nval,nset=nn_emul,d=k2_emul,var_k2=s2_emul, N=N,k1=k1,theta=theta,var_u=var_u,k_width=k_width)
    n_new = nval[which.max(nxn)]
    k2_new=aspre_k2(n_new,X,PRE)
    nn_emul=c(nn_emul,n_new)
    k2_emul=c(k2_emul,k2_new)
    s2_emul=c(s2_emul,dvar0)
    theta=powersolve_general(nn_emul,k2_emul)$par
    print(i)
  }

  # Transform estimated k2 to costs
  cc_emul=k1*nn_emul + k2_emul*(N-nn_emul)

  # Save
  aspre_emulation=list(nn_emul=nn_emul,k2_emul=k2_emul,s2_emul=s2_emul,cc_emul=cc_emul)
  save(aspre_emulation,file="data/aspre_emulation.RData")


} else load("data/aspre_emulation.RData")
for (i in 1:length(aspre_emulation)) assign(names(aspre_emulation)[i],aspre_emulation[[i]])


# Mean and variance of emulator for cost function, parametric assumptions satisfied
theta=powersolve_general(nn_emul,k2_emul)$par
p_mu=mu_fn(nval,nset=nn_emul,k2=k2_emul,var_k2 = s2_emul,theta=theta,
           N=N,k1=k1,var_u=var_u,k_width=k_width)
p_var=psi_fn(nval,nset=nn_emul,var_k2=s2_emul,N=N,var_u=var_u,
             k_width=k_width)


######################################################################
## Estimate optimum holdout size and measure of error               ##
######################################################################

OHS_ASPRE=nval[which.min(p_mu)]
MIN_COST_ASPRE=min(p_mu)
OHS_ERR=error_ohs_emulation(nn_emul,k2_emul,var_k2=s2_emul,N=N,k1=k1,alpha=0.1,
                            var_u=var_u,k_width=k_width,theta=theta)



######################################################################
## Trace of OHS at different N using emulation algorithm            ##
######################################################################

if (!file.exists("data/aspre_trace_emul.RData")) {
  
  ## Store OHS at various number of samples n
  ohs_trace_emul=rep(NA,length(nn_emul))
  mincost_trace_emul=rep(NA,length(nn_emul))
  error_trace_emul=list()
  
  ## Compute, starting with 5 points
  for (i in 5:length(nn_emul)) {
    nnsub=nn_emul[1:i]; k2sub=k2_emul[1:i]; s2sub=s2_emul[1:i]
    thetasub=powersolve_general(nnsub,k2sub)$par
    musub=mu_fn(nval,nset=nnsub,k2=k2sub,var_k2 = s2sub,theta=thetasub,
                N=N,k1=k1,var_u=var_u,k_width=k_width)
    psisub=psi_fn(nval,nset=nnsub,var_k2=s2sub,N=N,var_u=var_u,
                  k_width=k_width)
    
    ohs_trace_emul[i]=nval[which.min(musub)]
    mincost_trace_emul=min(musub)
    error_trace_emul[[i]]=error_ohs_emulation(nnsub,k2sub,var_k2=s2sub,N=N,k1=k1,alpha=0.1,
                                              var_u=var_u,k_width=k_width,theta=theta)
  }
  
  save(ohs_trace_emul,mincost_trace_emul,error_trace_emul,file="data/aspre_trace_emul.RData")
  
} else load("data/aspre_trace_emul.RData")



######################################################################
## Draw figure for cost fuction                                     ##
######################################################################

if (save_plot) pdf("figures/cost_function_estimate_emul.pdf",width=pdf_dim,height=pdf_dim)

plot(0,xlim=range(nn_emul),ylim=range(cc_emul),type="n",
     xlab="Training set size",yaxt="n",
     ylab=expression(paste("Total. cost ", "(", "","",
                           phantom() %prop% phantom(), " sens.", ")", "")))
axis(2,las=2)
points(nn_emul,cc_emul,pch=16,cex=0.5)
lines(nval,p_mu)
lines(nval,p_mu+3*sqrt(pmax(0,p_var)),col="blue")
lines(nval,p_mu-3*sqrt(pmax(0,p_var)),col="blue")
e_min=min(OHS_ERR); e_max=max(OHS_ERR); c_min=min(cc_emul); c_max=max(cc_emul);
polygon(c(e_min,e_min,e_max,e_max),c(c_min,c_max,c_max,c_min),
        col=rgb(0,0,1,alpha=0.2),border=NA)
points(OHS_ASPRE,MIN_COST_ASPRE,pch=16,col="red")

if (inc_legend) legend("topright",
       c(expression(mu(n)),
         expression(mu(n) %+-% 3*sqrt(psi(n))),
         "Est cost (d)",
         "OHS",
         "OHS err."),
       lty=c(1,1,NA,NA,NA),lwd=c(1,1,NA,NA,NA),pch=c(NA,NA,16,16,16),pt.cex=c(NA,NA,0.5,1,2),
       col=c("black","blue","black","red",rgb(0,0,1,alpha=0.2)),bg="white",bty="n")

if (save_plot) dev.off()




######################################################################
## Draw figure to track OHS at various N                            ##
######################################################################

if (save_plot) pdf("figures/aspre_track.pdf",width=pdf_dim,height=pdf_dim)

ymax=25000
par(mar=c(5.1, 4.8, 4.1, 2.1))
plot(0,type="n",xlim=c(0,length(nn_emul)),ylim=c(0,ymax),
     ylab="",xlab="|n|",yaxt="n")
axis(2,las=2); title(ylab = "Opt. holdout size", line = 3.5)
lines(1:length(nn_par),ohs_trace_par,type="l",col="red")
lines(1:length(nn_emul),ohs_trace_emul,type="l",col="blue",lty=2)
ip=c(); ie=c()
for (i in 1:length(nn_par)) {
  if (length(error_trace_emul[[i]]>0)) {
    ie=rbind(ie,range(error_trace_emul[[i]],na.rm=T))
  } else ie=rbind(ie,c(0,ymax))
  if (length(ci_trace_par[[i]]>0)) {
    ip=rbind(ip,range(ci_trace_par[[i]],na.rm=T))
  } else ip=rbind(ip,c(0,ymax))
}
ie[,1]=pmax(ie[,1],0); ip[,1]=pmax(ip[,1],0)
ie[,2]=pmin(ie[,2],ymax); ip[,2]=pmin(ip[,2],ymax)

# Emulation error
polygon(c(1:length(nn_par),length(nn_par):1),c(ie[,1],rev(ie[,2])),
        col=rgb(0,0,1,alpha=0.2),border=NA)
# Parametric error
polygon(c(1:length(nn_par),length(nn_par):1),c(ip[,1],rev(ip[,2])),
        col=rgb(1,0,0,alpha=0.2),border=NA)

if (inc_legend) legend("topright",c("Par. OHS","Em. OHS","Par. CI","Em. err."),
       lty=c(1,2,NA,NA),pch=c(NA,NA,16,16),pt.cex=c(NA,NA,2,2),
       col=c("red","blue",rgb(1,0,0,alpha=0.2),rgb(0,0,1,alpha=0.2)),
             bty="n")

if (save_plot) dev.off()
