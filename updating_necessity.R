##****************************************************************************##
## R script to demonstrate necessity of managing model updating safely      ####
##****************************************************************************##
##
## James Liley
## 19 April 2022
##
## This work is loosely based on the ASPRE model for pre-eclampsia. THIS IS A
##  SIMULATION ONLY, INTENDED ONLY TO DEMONSTRATE HOW THE ALGORITHM COULD BE
##  APPLIED. DO NOT INTERPRET ANY RESULTS AS CLINICALLY USABLE.
##
## In this simulation, we demonstrate the necessity of developing a protocol for
##  safe updating of a predictive model. We use a setting similar to that for 
##  the ASPRE score for predicting pre-eclampsia.
##
## We presume 'drift' occurs in that the distribution of (Y|X) will change over 
##  time (Y: target; X: covariates). This necessitates model updating. We also
##  presume that interventions are made on the basis of the predicted score, 
##  which may be on the covariates X1 included in the predictive score or on 
##  a 'latent' covariate X2 representing other components to risk. 
##
## We compare three strategies:
##  1. Not updating a model at all
##  2. Updating a model 'naively'
##  3. Updating a model using a hold-out set
##

## This script should be ran in the directory in which it is saved, or in 
##  some directory with subdirectories 'data', 'figures'.
## Not all figures are necessarily used in the manuscript.


##*******************************************************************#
## Scripts, switches and libraries                                ####
##*******************************************************************#

# Set random seed
seed=857366
set.seed(seed)

# Libraries
library(mvtnorm)
library(matrixStats)
library(mle.tools)
library(OptHoldoutSize)
library(MASS)

# Save plots to file, or not
save_plot=TRUE

# Force redo: set to TRUE to regenerate all datasets from scratch
force_redo=FALSE

# Include legends on plots for main paper
inc_legend=FALSE

# PDF dimensions in inches
pdf_width=6
pdf_height=3.5

# Print progress
verbose=TRUE

# Force redo: if this script has already been ran, it will save progress,
#  speeding up recomputation. 
force_redo=TRUE

# Print session information
sink("data/updating_necessity_session_info.txt")
sessionInfo()
sink()



##*******************************************************************#
## Parameters                                                     ####
##*******************************************************************#

#### For simulated population

# Total individuals in population at any time
Nt=200000

# Number of 'timepoints' to simulate
n_timepoint=50

# Drift type: linear (coefficients change linearly) or GP (coefficients follow a Gaussian process)
drift_type="GP"

# Rate of drift: higher means faster. Gradient for linear, inverse of covariance for radial kernel for GP.
drift_rate=1/50

# Contribution of latent covariate to risk: higher is more, meaning less predictable outcome
latent_intensity=3

# Intercept term. A low intercept term will tend to mean more people are at low risk.
intercept=-5

# Frequency of updates (number of timepoints between updates)
update_frequency= 10

# Holdout set sizes
holdout_sizes=c(10000,20000)


### Simple intervention (give aspirin to 10% highest risk)

# Intervene (with aspirin) on this proportion of individuals
pi_intervention=0.2

# Aspirin reduces latent PRE risk by this much
latent_change=4


### Complex intervention

# Intensity of intervention: higher means a more drastic intervention
intervention_level=1




##*******************************************************************#
## Random matrix generator                                        ####
##*******************************************************************#

# ASPRE details for simulation
data(params_aspre)
set.seed(10+seed)
s0=cbind(sim_random_aspre(1000,params=params_aspre),Y=1)
s1=model.matrix(Y~.+0,data=s0)
mu_aspre=colMeans(s1)
sd_aspre=colSds(s1)

## Utility function to generate normalised ASPRE data matrix. Vectors mu_n
##  and sd_n indicate normalising values.
aspre_numeric=function(n, mu_n=mu_aspre, sd_n=sd_aspre) {
  X=cbind(sim_random_aspre(n,params=params_aspre),Y=1)
  Xn=model.matrix(Y~.+0,data=X)
  Xm=t((t(Xn)-mu_aspre)/sd_aspre)
  return(as.data.frame(Xm))
}


## Define timepoints at which an update will take place (first at update_frequency)
update_points=seq(1, n_timepoint-update_frequency,by=update_frequency) + update_frequency-1



##*******************************************************************#
## Coefficients governing true risk with drift                    ####
##*******************************************************************#

## Define coefficients governing true risk. Coefficients keep same sign over time.
coefs_sign=sample(c(-1,1),length(mu_aspre),rep=T)

# Initialise 
coefs=matrix(0,length(coefs_sign),n_timepoint)

# Simulate coefficients over time, including drift
if (drift_type=="linear") {
  
  # Define coefficients: 
  for (p in 1:length(coefs_sign)) 
    coefs[p,]=coefs_sign[p]*abs(rnorm(1) + drift_rate*(1:n_timepoint))
  
  # Latent coefficients: smoothly variable but between 0 and 2 times latent_intensity
  latent_coefs=abs(rnorm(1) + drift_rate*(1:n_timepoint))
  
  for (t in 1:n_timepoint)
    coefs[,t]=(coefs[,t]-mean(coefs[,t]))/(3*sd(coefs[,t]))
}
if (drift_type=="GP") {
  
  # Variance matrix for Gaussian process governing change in coefficients over time
  c_mat=outer(1:n_timepoint,1:n_timepoint,function(x,y) exp(-(x-y)^2 /(2*(1/drift_rate)^2)))
  
  # Define coefficients: 
  for (p in 1:length(coefs_sign)) 
    coefs[p,]=coefs_sign[p]*abs(mvrnorm(1,mu=rep(0,n_timepoint),Sigma=c_mat))/3
  
  # Latent coefficients: smoothly variable but between 0 and 2 times latent_intensity
  latent_coefs=abs(mvrnorm(1,mu=rep(0,n_timepoint),Sigma=c_mat))
}


##*******************************************************************#
## Intervention and true risk specification                       ####
##*******************************************************************#

## Intervention function
intervene=function(score_orig,X0,Xl,holdout=NULL) {
  
  # Non-holdout set
  intervention=setdiff(1:dim(X0)[1],holdout)
  
  # Intervene on intervention set. For proportion pi_intervention of samples at 
  #  highest assessed risk, intervene by reducing latent covariate value by 
  #  amount latent_change.
  
  # Samples at high assessed risk
  w=which(score_orig[intervention]> quantile(score_orig[intervention],1-pi_intervention))
  
  # Change latent covariates
  Xl[intervention[w]]=Xl[intervention[w]] - latent_change
  
  # Returned data frame: contains both X0 and Xl
  out=data.frame(X0,Xl)  
  return(out)  
}


## True risk function
true_risk=function(X1,coefs,latent_coef) {
  # Response and output
  resp=intercept + (as.matrix(X1) %*% c(coefs,latent_coef)) # Response
  out=logit(resp)
  return(out)
}



##*******************************************************************#
## Initialise simulation                                          ####
##*******************************************************************#

costs_file="data/updating_costs.RData"
if (!file.exists(costs_file) | force_redo) {
  
  # Initialise cost records
  cost_orig=rep(0,n_timepoint)
  cost_naive=rep(0,n_timepoint)
  cost_holdout=matrix(0,n_timepoint,length(holdout_sizes))
  cost_alt=rep(0,n_timepoint)
  cost_semioracle=rep(0,n_timepoint)
  cost_oracle=rep(0,n_timepoint)
  cost_null=rep(0,n_timepoint)
  
  
  # Original model; presumed fitted to independent study data
  X0_study=aspre_numeric(Nt)
  Xl_study=rnorm(Nt) # Latent covariate (not included in model)
  X1_study=data.frame(X0_study,Xl_study) # Make no interventions
  risk_study=true_risk(X1_study,coefs=coefs[,1],latent_coef=latent_coefs[1]) # True risks in study population
  events_study=rbinom(Nt,1,prob=risk_study) # True
  
  # Fit model to study data; if we do not update, this will be used throughout
  data_study=cbind(X0_study,events_study)
  model_orig=suppressWarnings(glm(events_study~.,data=data_study,family=binomial(link="logit")))
  predict_orig=function(x) suppressWarnings(predict(model_orig,newdata=x,type="response"))
  
  # Initialise predictive model for naive updating and holdout set
  predict_naive=predict_orig
  for (i in 1:length(holdout_sizes)) assign(paste0("predict_holdout",i),predict_orig)
  
  # Initialise alternative predictor
  predict_alt=function(x) true_risk(cbind(x,Xl=0),coefs=coefs[,1],latent_coef=latent_intensity)
  
  # Initialise semioracle predictor (just 'oracle' on legend and in paper)
  predict_semioracle=function(x) true_risk(cbind(x,Xl=0),coefs=coefs[,1],latent_coef=latent_intensity)
  
  # Initialise oracle predictor
  predict_oracle=function(x) true_risk(cbind(x,Xl=0),coefs=coefs[,1],latent_coef=latent_intensity)
  
  
  ##*******************************************************************#
  ## Run simulation                                                 ####
  ##*******************************************************************#
  
  for (e in 1:n_timepoint) {
    
    # Set random seed
    set.seed(seed + 352*e)
    
    # Define holdout set, if we are updating at this time point
    if (e %in% update_points) {
      for (i in 1:length(holdout_sizes)) assign(paste0("holdout_indices",i),sample(Nt,holdout_sizes[i]))
    } else {
      for (i in 1:length(holdout_sizes)) assign(paste0("holdout_indices",i),NULL)
    }
    
    ## Data for this timepoint. We presume we have a new population at each timepoint.
    X0=aspre_numeric(Nt)
    Xl=rnorm(Nt) # Latent covariate (not included in model)
    
    
    ## Compute risk scores for this timepoint
    score_orig=predict_orig(X0)        # No updating
    score_naive=predict_naive(X0)      # Naively updating
    for (i in 1:length(holdout_sizes)) {
      predict_holdout=get(paste0("predict_holdout",i))
      assign(paste0("score_holdout",i),predict_holdout(X0))  # Holdout set
    }
    score_alt=predict_alt(X0)    # Alternative strategy; include partial knowledge of intervention as covariate (oracle in paper)
    score_semioracle=predict_semioracle(X0)    # Semi-oracle (oracle in paper)
    score_oracle=predict_oracle(X0)    # True oracle (perfect knowledge of f_t)
    
    # Interventions: for holdout strategy, do not intervene on holdout set 
    X1_orig=intervene(score_orig,X0,Xl)
    X1_naive=intervene(score_naive,X0,Xl)
    for (i in 1:length(holdout_sizes)) {
      score_holdout=get(paste0("score_holdout",i))
      holdout_indices=get(paste0("holdout_indices",i))
      assign(paste0("X1_holdout",i),intervene(score_holdout,X0,Xl,holdout=holdout_indices))
    }
    X1_alt=intervene(score_alt,X0,Xl)
    X1_semioracle=intervene(score_semioracle,X0,Xl)
    X1_oracle=intervene(score_oracle,X0,Xl)
    
    # Record partial knowledge of intervention for alternative strategy
    alt_int=as.numeric(score_alt>quantile(score_alt,1-pi_intervention))
    
    
    # Risk under no risk-score based intervention (what we are aiming to estimate; f_t)
    risk0=true_risk(cbind(X0,Xl),coefs=coefs[,e],latent_coef=latent_coefs[e]) # Risk pre-intervention (all)
    
    # Risk under best possible risk-score based intervention (oracle)
    risk_oracle=true_risk(X1_oracle,coefs=coefs[,e],latent_coef=latent_coefs[e])
    
    # Risk when we use various risk scores
    risk_orig=true_risk(X1_orig,coefs=coefs[,e],latent_coef=latent_coefs[e])
    risk_naive=true_risk(X1_naive,coefs=coefs[,e],latent_coef=latent_coefs[e])
    for (i in 1:length(holdout_sizes)) {
      X1_holdout=get(paste0("X1_holdout",i))
      assign(paste0("risk_holdout",i),true_risk(X1_holdout,coefs=coefs[,e],latent_coef=latent_coefs[e]))
    }
    risk_alt=true_risk(X1_alt,coefs=coefs[,e],latent_coef=latent_coefs[e])
    risk_semioracle=true_risk(X1_semioracle,coefs=coefs[,e],latent_coef=latent_coefs[e])
    
    # True (observed) events under various risk score strategies
    events_orig=rbinom(Nt,1,prob=risk_orig) # No risk score
    events_naive=rbinom(Nt,1,prob=risk_naive)
    for (i in 1:length(holdout_sizes)) {
      risk_holdout=get(paste0("risk_holdout",i))
      assign(paste0("events_holdout",i),rbinom(Nt,1,prob=risk_holdout))
    }
    events_alt=rbinom(Nt,1,prob=risk_alt)
    events_semioracle=rbinom(Nt,1,prob=risk_semioracle)
    events_oracle=rbinom(Nt,1,prob=risk_oracle)
    
    # Add to costs
    cost_orig[e]=sum(risk_orig)
    cost_naive[e]=sum(risk_naive)
    for (i in 1:length(holdout_sizes)) 
      cost_holdout[e,i]=sum(get(paste0("risk_holdout",i)))
    cost_alt[e]=sum(risk_alt)
    cost_semioracle[e]=sum(risk_semioracle)
    cost_oracle[e]=sum(risk_oracle)
    cost_null[e]=sum(risk0) 
    
    # If we are to update this timepoint, then do so
    if (e %in% update_points) {
      
      # Original model: do not chance
      #predict_orig=predict_orig
      
      # Naive updating: fit to all data from previous time point
      data_naive=cbind(X0,events_naive)
      model_naive=suppressWarnings(glm(events_naive~.,data=data_naive,family=binomial(link="logit")))
      predict_naive=function(x) suppressWarnings(predict(model_naive,newdata=x,type="response"))
      
      # Holdout sets: fit to all data from holdout set at previous time point
      for (i in 1:length(holdout_sizes)) {
        events_holdout=get(paste0("events_holdout",i))
        holdout_indices=get(paste0("holdout_indices",i))
        data_holdout=cbind(X0,events_holdout)[holdout_indices,]
        assign(paste0("model_holdout",i),suppressWarnings(glm(events_holdout~.,data=data_holdout,family=binomial(link="logit"))))
        assign(paste0("predict_holdout",i),function(x) suppressWarnings(predict(get(paste0("model_holdout",i)),newdata=x,type="response")))
      }
      
      # Predict including partial knowledge of intervention as a covariate
      data_alt=cbind(X0,alt_int,events_alt)
      model_alt=suppressWarnings(glm(events_alt~.,data=data_alt,family=binomial(link="logit")))
      predict_alt=function(x) suppressWarnings(predict(model_alt,newdata=cbind(x,events_alt=rep(0.5,dim(x)[1])),type="response"))
      
      # Predict using a (non-existent) set of observations of mu_t, f_t(x), x~\mu_t (call this semi-oracle; corresponds to oracle in the paper; using term because we have an actual 'oracle' here)
      events_semioracle=rbinom(Nt,1,prob=risk0) # Note we are regenerating a non-existent set of patients, using true risk f_t (risk0)
      data_semioracle=cbind(X0,events_semioracle)
      model_semioracle=suppressWarnings(glm(events_semioracle~.,data=data_semioracle,family=binomial(link="logit")))
      predict_semioracle=function(x) suppressWarnings(predict(model_semioracle,newdata=x,type="response"))
      
    }
    
    # Oracle: best estimate of risk, without knowledge of latent variable. Update this whether 'updating' or not.
    predict_oracle=function(x) true_risk(cbind(x,Xl=0),coefs=coefs[,e],latent_coef=latent_intensity)
    
    if (verbose) print(e)
  }
  
  save(cost_orig,cost_naive,cost_holdout,cost_alt,cost_semioracle,cost_oracle,cost_null,
       file=costs_file)
  
} else load(costs_file)



##*******************************************************************#
## Draw plot                                                    ######
##*******************************************************************#

# Differences
df_orig=(cost_orig-cost_oracle)/Nt
df_naive=(cost_naive-cost_oracle)/Nt
df_holdout=(cost_holdout-cost_oracle)/Nt
df_null=(cost_null-cost_oracle)/Nt
df_alt=(cost_alt-cost_oracle)/Nt
df_semioracle=(cost_semioracle-cost_oracle)/Nt
yr=range(c(df_orig,df_naive,df_holdout))

if (save_plot) pdf("figures/updating_necessity.pdf",width=pdf_width,height=pdf_height)
par(mgp=c(4,1,0),mar=c(5,5,1,1))

# Setup
plot(0,type="n",xlim=c(0,n_timepoint),ylim=yr,xlab="Time point",
     ylab="Excess cost",yaxt="n")
pxx=pretty(yr); plab=c(pxx[1:(length(pxx)-1)],paste0("> ",max(pxx)))
axis(2,at=pxx,labels=plab,las=2)

# Updating points
sx=update_frequency*(1:floor(n_timepoint/update_frequency))
ax=1:(n_timepoint-1); axn=ax; axn[sx]=NA

# No updating
points(df_orig,cex=0.3)
segments(axn,df_orig[axn],axn+1,df_orig[axn+1])

# Naive updating
points(df_naive,cex=0.3,col="red")
segments(axn,df_naive[axn],axn+1,df_naive[axn+1],col="red")

# Alternative updating
points(df_alt,cex=0.3,col="orange")
segments(axn,df_alt[axn],axn+1,df_alt[axn+1],col="orange")

# Semioracle updating
points(df_semioracle,cex=0.3,col="purple")
segments(axn,df_semioracle[axn],axn+1,df_semioracle[axn+1],col="purple")


# Holdout sizes
cx=gray(0.5+(1:length(holdout_sizes))/(2*(1+length(holdout_sizes))))
ccex=seq(0.6,0.8,length=length(holdout_sizes))
for (i in 1:length(holdout_sizes)) {
  points(df_holdout[,i],cex=ccex[i],col=cx[i],lty=i)
  segments(axn,df_holdout[axn,i],axn+1,df_holdout[axn+1,],col=cx[i])
}  

# Draw update points=}|


segments(sx,rep(0,length(sx)),sx,rep(0.5*max(c(df_orig,df_holdout)),length(sx)),lty=2,col="blue")

# Legend
hslab=c("None","Naive")
for (i in 1:length(holdout_sizes)) 
  hslab=c(hslab,as.expression(bquote("H.S:"~.(floor(holdout_sizes[i])/10000)~"x10"^5)))
hslab=c(hslab,"Alternative","Oracle","Update")
legend("topleft",legend=hslab,
       pch=c(1,1,rep(1,length(holdout_sizes)),1,1,NA),
       pt.cex=c(0.3,0.5,ccex,0.3,0.3,NA),
       lty=c(1,1,rep(1,length(holdout_sizes)),1,1,2),
       col=c("black","red",cx,"orange","purple","blue"),bty="n",title=" ")

if (save_plot) dev.off()


##*******************************************************************#
## Draw plot of cumulative costs                                  ####
##*******************************************************************#

if (save_plot) pdf("figures/updating_necessity_cumulative.pdf",width=pdf_width,height=pdf_height)
par(mgp=c(4,1,0),mar=c(5,5,1,1))

# Setup
yr2=c(0,1.5*sum(df_orig))
plot(0,type="n",xlim=c(0,n_timepoint),ylim=yr2,xlab="Time point",
     ylab="Excess cost (cumulative)",yaxt="n")
pxx=pretty(yr2); plab=c(pxx[1:(length(pxx)-1)],paste0("> ",max(pxx)))
axis(2,at=pxx,labels=plab,las=2)

# Points and lines: the 'points' function does not seem to work
pointz=function(y,...) {
  x=1:length(y)
  ax=1:(length(x)-1)
  points(x,y,...)
  segments(x[ax],y[ax],x[1+ax],y[1+ax],...)
}

# Updating points
sx=update_frequency*(1:floor(n_timepoint/update_frequency))
ax=1:(n_timepoint-1); axn=ax; axn[sx]=NA

# Points and lines as above
pointz(cumsum(df_orig),cex=0.3) # No update
pointz(cumsum(df_naive),cex=0.3,col="red") # Naive update
pointz(cumsum(df_alt),cex=0.3,col="orange") # Alternative update
pointz(cumsum(df_semioracle),cex=0.3,col="purple") # Semioracle updating

# Holdout updating
cx=gray(0.5+(1:length(holdout_sizes))/(2*(1+length(holdout_sizes))))
ccex=seq(0.6,0.8,length=length(holdout_sizes))
for (i in 1:length(holdout_sizes)) {
  pointz(cumsum(df_holdout[,i]),cex=ccex[i],col=cx[i],lty=i)
}  

# Draw update points
segments(sx,rep(0,length(sx)),sx,rep(0.5*sum(df_orig),length(sx)),lty=2,col="blue")

# Legend
hslab=c("None","Naive")
for (i in 1:length(holdout_sizes)) 
  hslab=c(hslab,as.expression(bquote("H.S:"~.(floor(holdout_sizes[i])/10000)~"x10"^5)))
hslab=c(hslab,"Alternative","Oracle","Update")
legend("topleft",legend=hslab,
       pch=c(1,1,rep(1,length(holdout_sizes)),1,1,NA),
       pt.cex=c(0.3,0.5,ccex,0.3,0.3,NA),
       lty=c(1,1,rep(1,length(holdout_sizes)),1,1,2),
       col=c("black","red",cx,"orange","purple","blue"),bty="n",title=" ")

if (save_plot) dev.off()



##***************************************************************##
## Lipschitz and intervention parameters                     ######
##***************************************************************##

parameter_file="data/updating_parameters.RData"
if (!file.exists(parameter_file) | force_redo) {
  
  # Brief simulation to establish parameters
  nsim_param=200000
  xtest=as.matrix(aspre_numeric(nsim_param)) # random draws from mu_t
  ltest=rnorm(nsim_param) # Latent coefficients
  
  # Lipschitz constant for f_t
  scores=xtest %*% coefs # Scores for each individual over time (f_t values)
  dscores=(cbind(scores,scores[,dim(scores)[2]])-cbind(scores[,1],scores))[,1:dim(scores)[2]] # Differentials
  
  # Parameters governing relation of cost to accuracy
  ## Function takes a set of coefficients governing predictive score rho, and a set of 'true' coefficients
  ##  governing f_t, and returns xi_t(rho_t,f_t) = E((rho_t(x)-f_t(x))^2
  xi_t=function(coefs_pred,coefs_true) {
    rho_t=logit(xtest %*% coefs_pred)
    f_t=logit(xtest %*% coefs_true)
    return(mean((rho_t-f_t)^2))
  }
  cost_t=function(coefs_pred,coefs_true,latent_coef_true=mean(latent_coefs)) {
    X1_pred=intervene(logit(xtest %*% coefs_pred),xtest,ltest)
    X1_orac=intervene(logit(xtest %*% coefs_true),xtest,ltest)
    true_risk_pred=true_risk(X1_pred,coefs_true,latent_coef = latent_coef_true )
    true_risk_orac=true_risk(X1_orac,coefs_true,latent_coef = latent_coef_true)
    return(mean(true_risk_pred-true_risk_orac))
  }
  coef_true=rowMeans(coefs)
  ntrial=1000
  xi_ts=rep(0,ntrial); costs=xi_ts
  for (i in 1:ntrial) {
    coef_pred=coef_true+rnorm(length(coef_true),sd=0.25)
    xi_ts[i]=xi_t(coef_pred,coef_true)
    costs[i]=cost_t(coef_pred,coef_true)
    if ((i %% 100)==0) print(paste0(i," of ",ntrial))
  }
  # Plot (only if save_plot=TRUE)
  if (save_plot) {
    pdf("figures/cost_vs_xi.pdf",width=pdf_width,height=pdf_height)
    plot(xi_ts,costs,xlim=c(0,max(xi_ts)),ylim=c(0,max(costs)),
         xlab=expression(paste(xi[t]^2,"(f"[t],", ",rho[t],")")),
         ylab=expression("c"[t]),
         pch=16,cex=0.5)
    abline(0,max(costs/xi_ts),col="red")
    abline(0,min(costs/xi_ts),col="red")
    dev.off()
  }
  
  ## Lipschitz constants
  alpha_1 = 0 # No drift in mu_t, for simplicity
  alpha_2=max(abs(dscores)) # Drift in f_t (maximum δ f_t(x)/δt over x)
  
  delta=update_frequency
  s=1 
  
  # k1 (governs dependence of cost on E_x(f_t(x)-g_t(x)))
  k1=mean(df_null)
  
  # k2 upper and lower (governs dependence of cost on xi(f_t,rho_t))
  k2_u=max(costs/xi_ts)
  k2_l=min(costs/xi_ts)
  
  save(alpha_1,alpha_2,delta,s,k1,k2_u,k2_l,file=parameter_file)
  
} else load(parameter_file)

## Lipschitz parameters
print(paste0("α₁ = ",alpha_1))
print(paste0("α₂ = ",alpha_2))

## Updating parameters, set above
print(paste0("δ = ",delta))
print(paste0("s = ",s))

## Dependence of cost on E_x(f_t(x)-g_t(x))
print(paste0("k₁ = ",k1))

## Dependence of cost on xi(f_t,rho_t)
print(paste0("k₂ᵘ = ", k2_u))
print(paste0("k₂ˡ = ", k2_l))
