################################################################################
## R script to demonstrate necessity of managing model updating safely        ##
################################################################################
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


######################################################################
## Scripts, switches and libraries                                  ##
######################################################################

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
#force_redo=FALSE

# Include legends on plots for main paper
inc_legend=FALSE

# PDF dimensions in inches
pdf_width=6
pdf_height=3.5

# Print progress
verbose=TRUE


######################################################################
## Parameters                                                       ##
######################################################################

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

# Mode of intervention: simple or complex
intervention_mode="simple"

# Frequency of updates (number of timepoints between updates)
update_frequency= 10

# Holdout set sizes
holdout_sizes=c(10000,20000)


### Simple intervention (give aspirin to 10% highest risk)

# Intervene (with aspirin) on this proportion of individuals
pi_intervention=0.2

# Aspirin reduces latent PRE risk by this much
latent_change=2


### Complex intervention

# Intensity of intervention: higher means a more drastic intervention
intervention_level=1




######################################################################
## Random matrix generator                                          ##
######################################################################

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



######################################################################
## Coefficients governing true risk with drift                      ##
######################################################################

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


######################################################################
## Intervention and true risk specification                         ##
######################################################################

## Intervention function
intervene=function(score_orig,X0,Xl,mode=intervention_mode,holdout=NULL) {
  
  # Non-holdout set
  intervention=setdiff(1:dim(X0)[1],holdout)
  
  # Simple or random intervention
  if (intervention_mode=="simple") sc=1 else sc=0
  
  # Intervene on intervention set. For proportion pi_intervention of samples at 
  #  highest assessed risk, intervene by reducing latent covariate value by 
  #  amount latent_change.
  
  # Change covariates by this proportion
  prop_change=logit(-latent_change) 
  
  # Samples at high assessed risk
  w=which(score_orig[intervention]> quantile(score_orig[intervention],1-pi_intervention))
  
  # Change latent covariates
#  Xl[intervention[w]]=Xl[intervention[w]]*prop_change*runif(1,sc,1)
  Xl[intervention[w]]=Xl[intervention[w]] - 4
  
  # Change other covariate values
#  for (p in 1:length(coefs_sign)) {
#    if (coefs_sign[p]>0) X0[intervention[w],p]=X0[intervention[w],p]*prop_change*runif(1,sc,1)
#  }
  
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



######################################################################
## Initialise simulation                                            ##
######################################################################

# Initialise cost records
cost_orig=rep(0,n_timepoint)
cost_naive=rep(0,n_timepoint)
cost_holdout=matrix(0,n_timepoint,length(holdout_sizes))
cost_oracle=rep(0,n_timepoint)

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

# Initialise oracle predictor
predict_oracle=function(x) true_risk(cbind(x,Xl=0),coefs=coefs[,1],latent_coef=latent_intensity)



######################################################################
## Run simulation                                                   ##
######################################################################

### Tracer. This will indicate the overlap between the top-10% of 
###  samples according to each risk measure.
highrisk_track=c()


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
  score_oracle=predict_oracle(X0)    # Oracle
  
  # Interventions: for holdout strategy, do not intervene on holdout set 
  X1_orig=intervene(score_orig,X0,Xl,mode=intervention_mode)
  X1_naive=intervene(score_naive,X0,Xl,mode=intervention_mode)
  for (i in 1:length(holdout_sizes)) {
    score_holdout=get(paste0("score_holdout",i))
    holdout_indices=get(paste0("holdout_indices",i))
    assign(paste0("X1_holdout",i),intervene(score_holdout,X0,Xl,mode=intervention_mode,holdout=holdout_indices))
  }
  X1_oracle=intervene(score_oracle,X0,Xl,mode=intervention_mode)
  
  # True risks
  risk0=true_risk(cbind(X0,Xl),coefs=coefs[,e],latent_coef=latent_coefs[e]) # Risk pre-intervention (all)
  risk_oracle=true_risk(X1_oracle,coefs=coefs[,e],latent_coef=latent_coefs[e])
  
  # Approximated risk
  risk_orig=true_risk(X1_orig,coefs=coefs[,e],latent_coef=latent_coefs[e])
  risk_naive=true_risk(X1_naive,coefs=coefs[,e],latent_coef=latent_coefs[e])
  for (i in 1:length(holdout_sizes)) {
    X1_holdout=get(paste0("X1_holdout",i))
    assign(paste0("risk_holdout",i),true_risk(X1_holdout,coefs=coefs[,e],latent_coef=latent_coefs[e]))
  }
  
  w1=which(risk_oracle>quantile(risk_oracle,0.9))
  nwo=length(intersect(w1,which(risk_orig>quantile(risk_orig,0.9))))
  nwn=length(intersect(w1,which(risk_naive>quantile(risk_naive,0.9))))
  nwh=c(); 
  for (i in 1:length(holdout_sizes)) {
    risk_holdout=get(paste0("risk_holdout",i))
    nwh=c(nwh,length(intersect(w1,which(risk_holdout>quantile(risk_holdout,0.9)))))
  }
  highrisk_track=rbind(highrisk_track,c(nwo,nwn,nwh))
  if (e %in% c(1,4,15)) {
    #browser()
  }

  # True events
  events_orig=rbinom(Nt,1,prob=risk_orig)
  events_naive=rbinom(Nt,1,prob=risk_naive)
  for (i in 1:length(holdout_sizes)) {
    risk_holdout=get(paste0("risk_holdout",i))
    assign(paste0("events_holdout",i),rbinom(Nt,1,prob=risk_holdout))
  }
  events_oracle=rbinom(Nt,1,prob=risk_oracle)
  
  # Add to costs
  cost_orig[e]=sum(risk_orig)
  cost_naive[e]=sum(risk_naive)
  for (i in 1:length(holdout_sizes)) 
    cost_holdout[e,i]=sum(get(paste0("risk_holdout",i)))
  cost_oracle[e]=sum(risk_oracle)

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
  }
  
  # Oracle: best estimate of risk, without knowledge of latent variable. Update this whether 'updating' or not.
  predict_oracle=function(x) true_risk(cbind(x,Xl=0),coefs=coefs[,e],latent_coef=latent_intensity)
  
  if (verbose) print(e)
}






######################################################################
## Draw plot                                                        ##
######################################################################

# Differences
df_orig=cost_orig-cost_oracle
df_naive=cost_naive-cost_oracle
df_holdout=cost_holdout-cost_oracle
df_naive=pmin(df_naive,100*round(1.5*max(c(df_orig,df_holdout))/100))
yr=100*round(range(c(df_orig,df_naive,df_holdout))/100)


if (save_plot) pdf("figures/updating_necessity.pdf",width=pdf_width,height=pdf_height)

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
points(df_naive,cex=0.5,col="red")
segments(axn,df_naive[axn],axn+1,df_naive[axn+1],col="red")

# Holdout sizes
cx=gray(0.5+(1:length(holdout_sizes))/(2*(1+length(holdout_sizes))))
ccex=seq(0.6,0.8,length=length(holdout_sizes))
for (i in 1:length(holdout_sizes)) {
  points(df_holdout[,i],cex=ccex[i],col=cx[i],lty=i)
  segments(axn,df_holdout[axn,i],axn+1,df_holdout[axn+1,],col=cx[i])
}  

# Draw update points
segments(sx,rep(0,length(sx)),sx,rep(500,length(sx)),lty=2,col="blue")

# Legend
hslab=c("None","Naive")
for (i in 1:length(holdout_sizes)) 
  hslab=c(hslab,as.expression(bquote("H.S:"~.(floor(holdout_sizes[i])/10000)~"x10"^5)))
hslab=c(hslab,"Update")
legend(0,2300,legend=hslab,
       pch=c(1,1,rep(1,length(holdout_sizes)),NA),
       pt.cex=c(0.3,0.5,ccex,NA),
       lty=c(1,1,rep(1,length(holdout_sizes)),2),
       col=c("black","red",cx,"blue"),bty="n")

if (save_plot) dev.off()

