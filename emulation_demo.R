################################################################################
## R script to demonstrate emulation approach to OHS estimation               ##
################################################################################
##
## Sam Emerson and James Liley
## 21 December 2021
##
## This script should be ran in the directory in which it is saved, or in 
##  some directory with subdirectories 'data', 'figures'.
## Not all figures are necessarily used in the manuscript.

######################################################################
## Scripts, switches and libraries                                  ##
######################################################################

# Set random seed
seed=958436
set.seed(seed)

# Libraries


# Save plots to file, or not
save_plot=FALSE

# Force redo: set to TRUE to regenerate all datasets from scratch
force_redo=FALSE


######################################################################
## Parameters                                                       ##
######################################################################

# Suppose we have population size and cost-per-sample without a risk score as follows
N=100000
k1=0.4

# Kernel width and variance for GP
k_width=5000
var_u=8000000


######################################################################
## Simulation                                                       ##
######################################################################


# We begin with k2() estimates at n-values
nset=c(10000,20000,30000)

# with cost-per-individual estimates 
# (note that since empirical k2(n) is non-monotonic, it cannot be perfectly 
#  approximated with a power-law function)
k2=c(0.35,0.26,0.28)

# and associated error on those estimates
var_k2=c(0.02^2,0.01^2,0.03^2)

# We estimate theta from these three points
theta=powersolve(nset,k2,y_var=var_k2)$par

# We will estimate the posterior at these values of n
n=seq(1000,50000,length=1000)

# Mean and variance
p_mu=mu_fn(n,nset=nset,k2=k2,var_k2 = var_k2, N=N,k1=k1,theta=theta,
           k_width=k_width,var_u=var_u)
p_var=psi_fn(n,nset=nset,N=N,var_k2 = var_k2,k_width=k_width,var_u=var_u)

# Expected improvement
exp_imp=exp_imp_fn(n,nset=nset,k2=k2,var_k2=var_k2,N=N,k1=k1,var_u=var_u,
                   k_width=k_width)

######################################################################
## Draw figures                                                     ##
######################################################################

if (!save_plot) par(mfrow=c(1,2))

if (save_plot) pdf("figures/initial_emulator.pdf",width=5,height=5)
  
# Plot empirical and approximated cost function
plot(0,xlim=range(n),ylim=c(20000,60000),type="n",
     xlab="Training/holdout set size",
     ylab="Total cost (= num. cases)")
lines(n,p_mu,col="blue")
lines(n,p_mu - 3*sqrt(p_var),col="red")
lines(n,p_mu + 3*sqrt(p_var),col="red")
points(nset,k1*nset + k2*(N-nset),pch=16,col="purple")
lines(n,k1*n + powerlaw(n,theta)*(N-n),lty=2)
segments(nset,k1*nset + (k2 - 3*sqrt(var_k2))*(N-nset),
         nset,k1*nset + (k2 + 3*sqrt(var_k2))*(N-nset))
legend("topright",
       c(expression(mu(n)),
         expression(mu(n) %+-% 3*sqrt(psi(n))),
         "m(n)",
         expression("d"^1),
         expression(paste("3SD(d"^1,")"))),
       lty=c(1,1,2,NA,NA),lwd=c(1,1,1,NA,NA),pch=c(NA,NA,NA,16,124),
       pt.cex=c(NA,NA,NA,1,1),bty="n",
       col=c("blue","red","black","purple","black"),bg="white")

if (save_plot) dev.off()



if (save_plot) pdf("figures/initial_expected_improvement.pdf",width=5,height=5)

# Plot expected improvement function
plot(0,xlim=range(n),ylim=range(exp_imp),type="n",
     xlab="Training/holdout set size",
     ylab="Expected improvement")
lines(n,exp_imp,col="black")
abline(v=n[which.max(exp_imp)],col="red")
points(nset,rep(0,length(nset)),pch=16)
abline(v=nset,col="gray",lty=2)

legend("topright",
       c("Exp. impr.",
         "Next n",
         "Current n"),
       lty=c(1,1,NA),lwd=c(1,1,NA),pch=c(NA,NA,16),bty="n",
       col=c("black","red","black"),bg="white")

if (save_plot) dev.off()