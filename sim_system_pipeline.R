##*****************************************************************************#
## R script to simulate the emergence of an optimal holdout set size on a   ####
##  population for which intervention is guided by a risk score             ####
##*****************************************************************************#
##
## Sami Haidar-Wehbe and James Liley
## 1 November 2021
##

## This script should be ran in the directory in which it is saved, or in 
##  some directory with subdirectories 'data', 'figures'.
## Not all figures are necessarily used in the manuscript.

##*****************************************************************************#
## Setup                                                                    ####
##*****************************************************************************#

# Package loading
library("dplyr")           # Syntax
library("ranger")          # Random forests
library("OptHoldoutSize")  # Functions for OHS estimation
library("showtext")        # Unicode in plots

## Set seed for reproducibility
seed=1234
set.seed(seed)

# Save plots to file, or not
save_plot=TRUE

# Force redo: set to TRUE to regenerate all datasets from scratch
force_redo=TRUE

# Print session information
sink("data/sim_system_pipeline_session_info.txt")
sessionInfo()
sink()



##*****************************************************************************#
## Set parameters                                                           ####
##*****************************************************************************#

# Initialisation of patient data
n_iter = 500           # Number of point estimates to be calculated
nobs = 5000            # Number of observations, i.e patients
npreds = 7             # Number of predictors

# We will loop over these parameters. Number of interaction terms governs
#  deviation from logistic regression assumptions.
families = c("log_reg","rand_forest")   # Model family
interactions = c(0,10)                  # Number of interaction terms

# We allow multiple levels of predictive accuracy at 'baseline' (e.g., no risk
#  score). We will characterise baseline behaviour as an oracle Bayes-optimal
#  predictor on some subset of variables. These flags and vars generate multiple
#  such predictive powers.
max_base_powers = 1 # When testing different Dr predictive powers, maximum power to be tested

# Checks
if(max_base_powers > npreds) stop("max_base_powers cannot be larger than npreds")


# Check holdout sizes within these two proportions of total
min_frac = 0.02
max_frac = 0.15
num_sizes = 50

# Fraction of patients assigned to the holdout set
frac_ho = seq(min_frac, max_frac, length.out = num_sizes)


# Variables for the threshold estimation
c_tn = 0 # Cost of true neg
c_tp = 0.5 # Cost of true pos
c_fp = 0.5 # Cost of false pos
c_fn = 1 # Cost of false neg
num_thresh = 10 # Number of probabilities to scan
prob_vals = seq(0, 1, length.out = num_thresh) # Vector with probability thresh to scan

# Set whether to estimate optimal classification thresholds, or just use 1/2
thresh_vals = seq()
k = 5 # Number of CV folds
set_thresh = FALSE # Flag to estimate classification threshold


# Cost estimation
costs_tot_resample = matrix(nrow = n_iter, ncol = num_sizes)
cost_type = "gen_cost" # gen_cost = generalised cost, deaths = deaths (i.e. fn)
cost_mat = rbind(c(c_tn, c_fp), c(c_fn, c_tp))


# Matrix with the number of deaths per epoch and holdout size
costs_inter_resample = costs_ho_resample = matrix(nrow = n_iter, ncol = num_sizes)



##*****************************************************************************#
## Simulate system dynamics                                                 ####
##*****************************************************************************#

# Generate coefficients for underlying model and for predictions in h.o. set
set.seed(seed + 1)

# Set ground truth coefficients, and the accuracy at baseline
coefs_general = rnorm(npreds,sd=1/sqrt(npreds))
coefs_base = gen_base_coefs(coefs_general, max_base_powers = max_base_powers)

old_progress=0
if (!file.exists("data/data_example_simulation.RData")|force_redo) {

  for (ninters in interactions) {

    # If ninters>0, append coefficients for interaction terms to coefficients
    if (ninters>0)
      coefs_generate=c(coefs_general,rnorm(ninters,sd=2/sqrt(npreds)))
    else
      coefs_generate=coefs_general

    for (family in families) {

      # Arrays to populate
      costs_tot_resample = array(0, dim = c(max_base_powers, n_iter, num_sizes))
      costs_sd = array(0, dim = c(max_base_powers, num_sizes))


      print(paste0("Starting with interactions = ",ninters," and family = ",family))
      
      # Sweep through h.o. set sizes of interest
      for (i in 1:num_sizes) {
        # Look over resamplings
        for (b in 0:n_iter) {

          progress = 100 * (((i - 1) * n_iter) + b) / (num_sizes * (n_iter + 1))
          if (abs(floor(old_progress) - floor(progress)) > 0) {
            cat(floor(progress), "%\n")
          }

          # Set random seed
          set.seed(seed + b + i*n_iter)

          # Generate dataset
          X = gen_preds(nobs, npreds, ninters)

          # Generate labels
          newdata = gen_resp(X, coefs = coefs_generate)
          Y = newdata$classes

          # This contains only non-interaction columns and is used to train models.
          Xtrain=X[,1:npreds]

          # Combined dataset
          pat_data = cbind(Xtrain, Y)
          pat_data$Y = factor(pat_data$Y)

          # For each holdout size, split data into intervention and holdout set
          mask = split_data(pat_data, frac_ho[i])
          data_interv = pat_data[!mask,]
          data_hold = pat_data[mask,]

          if (set_thresh) {
          #### Calculate optimal threshold
            indices = sample(1:nrow(data_hold)) # Shuffled indices of holdout set
            folds = cut(1:length(indices), breaks = k, labels = FALSE) # Mask to partition data in CV analysis
            cost_tot = numeric(num_thresh)

            for (f in 1:k) {
              # Train model for each fold of CV, then scan all probs to find loss
              val_indices = folds == f
              val_data = data_hold[val_indices,]
              partial_train_data = data_hold[!val_indices,]

              thresh_model = model_train(partial_train_data, model_family = family)

              thresh_pred = model_predict(val_data, thresh_model, return_type = "probs",
                                           model_family = family)


              for(p in 1:num_thresh) {
                num_tn = sum(as.numeric(val_data["Y"] == 0) & as.numeric(thresh_pred < prob_vals[p]))
                num_fn = sum(as.numeric(val_data["Y"] == 1) & as.numeric(thresh_pred < prob_vals[p]))
                num_fp = sum(as.numeric(val_data["Y"] == 0) & as.numeric(thresh_pred >= prob_vals[p]))
                num_tp = sum(as.numeric(val_data["Y"] == 1) & as.numeric(thresh_pred >= prob_vals[p]))

                cost_tot[p] = cost_tot[p] + c_tn * num_tn +
                  c_tp * num_tp +
                  c_fn * num_fn +
                  c_fp * num_fp
              }
            }

            # Rescale loss
            cost_tot = cost_tot / k
          }

          # Train model
          trained_model = model_train(data_hold, model_family = family)
          thresh = ifelse(set_thresh, prob_vals[which.min(cost_tot)], 0.5)

          # Predict
          class_pred = model_predict(data_interv, trained_model,
                                      return_type = "class",
                                      threshold = thresh, model_family = family)



          for (base_vars in 1:max_base_powers) { # sweep through different dr predictive powers
            base_pred = oracle_pred(data_hold, coefs_base[base_vars, ])

            # Those with disease, predicted not to die, will die
            # This "if clause" calculates the cost for each h.o. set size
            # as the number of deaths for each of them
            if (cost_type == "deaths") {
                costs_inter_resample[b, i] = sum(data_interv$Y == 1 & class_pred != 1)
                costs_ho_resample[b, i] = sum(data_hold$Y == 1 & base_pred != 1)
            }


            # Alternatively, use a generalised cost function with a specific cost for
            # each of fn, fp, tp and tn
            if (cost_type == "gen_cost"){
              # Generate confusion matrices
              confus_inter = table(factor(data_interv$Y, levels=0:1),
                                    factor(class_pred, levels=0:1))

              confus_hold = table(factor(data_hold$Y, levels=0:1),
                                   factor(base_pred, levels=0:1))

                costs_inter_resample[b, i] = sum(confus_inter * cost_mat)
                costs_ho_resample[b, i] = sum(confus_hold * cost_mat)
                costs_tot_resample[base_vars, b, i] = costs_ho_resample[b, i] + costs_inter_resample[b, i]
            }

          }

          old_progress = progress
        }

      }

      # Calculate Standard Deviations and Errors
      for (base_vars in 1:max_base_powers) {
        costs_sd[base_vars, ] =
          apply(costs_tot_resample[base_vars, , ], 2, sd)
      }
      costs_se = costs_sd / sqrt(n_iter)

      # Store data for this model family and number of interactions
      data_list=list(max_base_powers=max_base_powers,
                     frac_ho=frac_ho,
                     costs_tot_resample=costs_tot_resample,
                     costs_ho_resample=costs_ho_resample,
                     costs_inter_resample=costs_inter_resample,
                     base_vars=base_vars,
                     costs_sd=costs_sd,
                     costs_se=costs_se)
      assign(paste0("data_",family,"_inter",ninters),data_list)
    }
  }

  # Collate and save all data
  data_example_simulation=list()
  ind=1
  for (family in families) {
    for (ninters in interactions) {
      xname=paste0("data_",family,"_inter",ninters)
      data_example_simulation[[ind]]=get(xname)
      names(data_example_simulation)[ind]=xname
      ind=ind+1
    }
  }
  save(data_example_simulation,file="data/data_example_simulation.RData")

} else load("data/data_example_simulation.RData")





##*****************************************************************************#
## Draw plots                                                               ####
##*****************************************************************************#

if (!save_plot) par(mfrow=c(2,2))

# Loop through values of 'families' and 'interactions'
for (xf in 1:length(families)) {
  for (xi in 1:length(interactions)) {
    family=families[xf]
    ninters=interactions[xi]

    if (save_plot) {
      showtext_auto()
      pdf(paste0("figures/example_simulation_",family,"_int",ninters,".pdf"),
        width=4,height=4)
    }


    X=data_example_simulation[[paste0("data_",family,"_inter",ninters)]]
    for (i in 1:length(X)) assign(names(X)[i],X[[i]])

    # Plot title
    if (xi==1 & xf==1) ptitle=expression(paste("Lin. und., lin. ",rho))
    if (xi==2 & xf==1) ptitle=expression(paste("Non-lin. und., lin. ",rho))
    if (xi==1 & xf==2) ptitle=expression(paste("Lin. und., non-lin. ",rho))
    if (xi==2 & xf==2) ptitle=expression(paste("Non-lin. und., non-lin. ",rho))

    # Y axis range
    if (xi==1) yl=c(2100,2350) else yl=c(2450,2700)

    # Holdout set size in absolute terms
    n_ho=frac_ho*nobs

    # Colour intensities
    r_alpha=0.2
    b_alpha=0.6

    # Create figure
    par(mar=c(4,4,3,4))
    plot(0, 0, type = "n",
         ylab = expression(paste("\u2113(n)")),
         xlab = "Holdout set size (n)",
         main=ptitle,
         xlim=range(n_ho),
         ylim = yl
    )


  ##### Draw learning curve with axis on right

    # Compute
    k2=colMeans(costs_inter_resample)/(nobs-n_ho)

    # Scaling
    if (xi==2) scvec=c(0.48,20) else scvec=c(0.42,20)
    k2_sc=function(x) yl[1] + (yl[2]-yl[1])*(x-scvec[1])*scvec[2] # Scaling transform
    i_k2_sc=function(x) (x-yl[1])/(scvec[2]*(yl[2]-yl[1])) + scvec[1] # Inverse

    # Add standard deviation. Note this is scaled by (nobs-n_ho)
    sd_k2=colSds(costs_inter_resample)/(nobs-n_ho)
    polygon(c(n_ho, rev(n_ho)),
            k2_sc(c(k2 + sd_k2,rev(k2-sd_k2))),
            col = rgb(1,0,0,alpha=r_alpha),
            border = NA)

    # Add axis and label
    axis(4,at=k2_sc(pretty(i_k2_sc(yl))),labels=pretty(i_k2_sc(yl)))
    mtext(expression(paste("k"[2],"(n)")), side=4, line=3)


    # Estimate parameters.
    theta=c(1,1,1);
    for (i in 1:5) theta=powersolve(n_ho,k2,init=theta)$par
    k1=median(colMeans(costs_ho_resample)/(nobs*frac_ho))
    k2=theta[1]*(n_ho)^(-theta[2]) + theta[3]

    # Line for estimated k2 curve
    lines(n_ho,
          k2_sc(k2),
          col="red",lty=2)

    # Line for estimated cost function
    xcost=n_ho*k1 + k2*(nobs-n_ho)
    lines(n_ho,
          xcost,
          col="blue",lty=2)

    # Add minimum L and OHS
    ohs=n_ho[which.min(xcost)]
    minL=min(xcost)
    lines(c(ohs,ohs),c(0,minL),lty=2)
    abline(h=minL,lty=2)
    points(ohs,minL,pch=4)

    # Finally, add cost function (on top)
    # Plot Cost line
    lines(n_ho,
          colMeans(costs_tot_resample[base_vars, , ]),
          pch = 16,
          lwd = 1,
          col = "blue")

    # Standard Deviation Bands
    polygon(c(n_ho, rev(n_ho)),
            c(colMeans(costs_tot_resample[1,,]) - costs_sd[1, ],
              rev(colMeans(costs_tot_resample[1,,]) + costs_sd[1, ])),
            col = rgb(0,0,1,alpha=b_alpha),
            border = NA)


    # Add legend
    legend("topright", legend = c("\u2113(n)","Fitted",
                                  "OHS",
                                  expression(paste("k"[2],"(n)")), "Fitted",
                                  "Min. \u2113"),
           col = c("blue","blue","black",
                   "red","red","black"),
           lty=c(1,2,NA,
                 1,2,2),
           pch=c(NA,NA,4,
                 NA,NA,NA),
           bty="n",bg="white",
           ncol=2
    )

    if (save_plot) dev.off()
  }
}
