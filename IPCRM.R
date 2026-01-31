setwd('/rsrch8/home/biostatistics/syang10/IPDE-Project')
.libPaths("/rsrch8/home/biostatistics/syang10/R/x86_64-pc-linux-gnu-library/4.4")
##########################################################################
#                         IP-CRM Design simulations                      #
##########################################################################
# Algorithmic/structural implementation of Section 2.1 with
# patient-level dose/toxicity sequences stored for debugging.
#
# PI        --> True toxicity probabilities (length J)
# TARGET    --> Target toxicity rate phi
# p         --> Skeleton / initial guess (length J)
# COHORTSIZE--> cohort size c
# ncohort   --> number of cohorts (total N = COHORTSIZE * ncohort patients)
# K         --> max number of DIFFERENT doses a patient may receive
# c_stop    --> early stopping cutoff c (placeholder)
#
# NOTE: All CRM model fitting/posterior pieces are placeholders.
##########################################################################
library(coda)
library(rjags)
# ---- PLACEHOLDER: estimate MTD from all accumulated (y,d) ----
# Returns an integer in 1:J
posterior <- function(beta0,beta1,y,d,p,lambda = 1, sigma = 0.5, mu = 3){
  #likelihood
  if (beta1 < 0) return(0)
  lik = 1
  for(i in 1:length(y)){
    eta <- beta0 + beta1 * p[d[i]]
    pi <- plogis(eta)
    # pi = exp(p[d[i]] * beta1 + beta0)/(1+exp(p[d[i]] * beta1 + beta0))
    lik = lik * pi^y[i] * (1-pi)^(1-y[i])
  }
  #like * prior 
  #prior beta1 ~ exp(1)
  #prior beta0 ~ N(3,0.5^2) 
  lik = lik * lambda * exp(-lambda * beta1)
  lik = lik * exp(-0.5 * ((beta0 - mu)^2)/sigma)
  return(lik)
}
estimate_MTD <- function(y, d, p, TARGET, lambda, sigma, mu) {
  posttox <- rep(0, length(p))
  inner_den <- function(beta1) {
    sapply(beta1, function(b1) integrate(function(b0) posterior(b0,b1, y,d,p, lambda, sigma, mu), 
                                         -Inf, Inf)$value)
  }
  marginal <- integrate(inner_den, 0, Inf)$value
  for(i in 1:length(p)){
    inner_tox <- function(beta1){
      sapply(beta1, function(b1) integrate(function(b0) plogis(b0 + b1 * p[i])* posterior(b0, b1, y,d,p, lambda, sigma, mu), 
                                           -Inf, Inf)$value)
    }
    posttox[i] = (integrate(inner_tox, 0,Inf)$value)/marginal
  }
  # cat('tox',posttox, '\n')
  diff = abs(posttox-TARGET);
  dose.best = min(which(diff==min(diff)))
  return(list(MTD = dose.best, posttox = posttox))
}
overtox_prob <- function(y, d, p, TARGET, lambda, sigma, mu) {
  q <- qlogis(TARGET)   # logit(phi)
  p1 <- p[1]            # effective dose at level 1
  
  # denominator Z(beta1) = ∫ k(b0,b1) db0
  inner_den <- function(beta1) {
    sapply(beta1, function(b1) integrate(function(b0) posterior(b0,b1, y,d,p,lambda, sigma, mu), 
                                         -Inf, Inf)$value)
  }
  
  denom <- integrate(inner_den, lower = 0, upper = Inf)$value
  
  # numerator: integrate k(b0,b1) over region pi1>phi
  inner_num <- function(beta1) {
    sapply(beta1, function(b1) {
      lb <- q - b1 * p1  # boundary for b0 given b1
      integrate(function(b0)
        posterior(b0, b1, y, d, p, lambda, sigma, mu),
        lower = lb, upper = Inf
      )$value
    })
  }
  
  num <- integrate(inner_num, lower = 0, upper = Inf)$value
  
  return(num/denom)
}
select_dose <- function(post_prob, dose_trt,TARGET) {
  treated <- sort(unique(dose_trt))
  diff <- abs(post_prob[treated] - TARGET)
  cand    <- treated[diff[treated] == min(diff[treated])]
  min(cand)
}

# ------------------------------------------------------------
# Posterior MTD estimator (handles all 4 model types)
# ------------------------------------------------------------
estimate_MTD_JAGS <- function(y, d_level, p, Tcum, pid,
                              TARGET = 0.3,
                              cutoff = 0.96,
                              model_file = "logit.bug",
                              n.chains = 3,
                              n.adapt  = 1000,
                              n.burn   = 2000,
                              n.iter   = 5000,
                              thin     = 2,
                              # for CO / CO-RF posterior curve definition
                              T_ref_for_curve = 0) {
  
  stopifnot(length(y) == length(d_level),
            length(y) == length(Tcum),
            length(y) == length(pid))
  
  Nobs <- length(y)
  nPat <- if (Nobs == 0) 0L else max(pid)
  
  # d in the JAGS model is the numeric dose on model scale
  # Here: use standardized/working dose values p[d_level]
  d_numeric <- as.numeric(p[d_level])
  T_numeric <- as.numeric(Tcum)
  pid_int   <- as.integer(pid)
  
  data_jags <- list(
    Nobs = Nobs,
    nPat = nPat,
    y_bin = as.integer(y),
    d = d_numeric,
    T = T_numeric,
    pid = pid_int
  )
  
  jags <- rjags::jags.model(
    file   = model_file,
    data   = data_jags,
    n.chains = n.chains,
    n.adapt  = n.adapt,
    quiet = TRUE
  )
  
  update(jags, n.burn, progress.bar = "none")
  
  # Decide which parameters to monitor
  is_CO   <- model_file %in% c("logit_CO.bug", "logit_CORF.bug")
  is_RF   <- model_file %in% c("logit_RF.bug", "logit_CORF.bug")
  
  var.names <- c("beta0", "beta1")
  if (is_CO) var.names <- c(var.names, "beta2")
  if (is_RF) var.names <- c(var.names, "sigma") # u[] not needed for dose curve
  
  smp <- coda.samples(
    model          = jags,
    variable.names = var.names,
    n.iter         = n.iter,
    thin           = thin,
    progress.bar   = "none"
  )
  
  draws <- as.matrix(smp)
  b0 <- draws[, "beta0"]
  b1 <- draws[, "beta1"]
  b2 <- if (is_CO) draws[, "beta2"] else 0
  
  J <- length(p)
  
  # Dose-level posterior mean tox curve
  # - logit/RF: logit^{-1}(b0 + b1*p[j])
  # - CO/CO-RF: logit^{-1}(b0 + b1*p[j] + b2*T_ref_for_curve)
  posttox <- vapply(seq_len(J), function(j) {
    mean(plogis(b0 + b1 * p[j] + b2 * T_ref_for_curve))
  }, numeric(1))
  
  # MTD = closest to TARGET
  diff <- abs(posttox - TARGET)
  dose.best <- which(diff == min(diff))[1]
  
  # Early stop rule: Pr(p1 > TARGET | data) > cutoff
  # Define p1 consistently with curve definition (T_ref_for_curve).
  pi1_draws <- plogis(b0 + b1 * p[1] + b2 * T_ref_for_curve)
  prob_overtox <- mean(pi1_draws > TARGET)
  stop_flag <- as.integer(prob_overtox > cutoff)
  
  list(
    MTD = dose.best,
    posttox = posttox,
    prob_overtox = prob_overtox,
    stop = stop_flag
  )
}

# ------------------------------------------------------------
# A simple final dose selector (example)
# - You can replace with your own select_dose()
# ------------------------------------------------------------
select_dose_closest <- function(posttox, TARGET) {
  which.min(abs(posttox - TARGET))[1]
}


# ------------------------------------------------------------
# Main IPCRM simulator with:
# - observation-level y_all, d_all, T_all, pid_all
# - correct Tij: cumulative dose PRIOR to current cycle
# - CO standardization applied to p (only for CO/CO-RF)
# ------------------------------------------------------------
IPCRM <- function(
    PI,
    TARGET = 0.3,
    p,
    COHORTSIZE = 3,
    ncohort = 10,
    ntrial,
    K = 3,
    c_stop = 0.96,
    seed = 6,
    gamma = 0.2, 
    model_file = "logit.bug",
    T_ref_for_curve = 0,
    verbose = FALSE
) {
  stopifnot(length(PI) == length(p))
  stopifnot(COHORTSIZE >= 1, ncohort >= 1, ntrial >= 1, K >= 1)
  
  set.seed(seed)
  J <- length(p)
  
  is_CO <- model_file %in% c("logit_CO.bug", "logit_CORF.bug")
  
  # For CO / CO-RF: standardize dose values used in model
  # NOTE: PI is indexed by dose level and stays on level scale
  p_model <- p
  if (is_CO) {
    p_model <- (p_model - mean(p_model)) / (2 * sd(p_model))
  }
  
  dose.select <- numeric(J)
  nstop <- 0L
  
  # optional summaries (keep your originals if you like)
  nobs     <- rep(0L, J)
  ntox_obs <- rep(0L, J)
  
  patient_seq_by_trial <- vector("list", ntrial)
  
  for (trial in 1:ntrial) {
    
    # observation-level data for JAGS
    y_all   <- integer(0)
    d_all   <- integer(0)  # dose level index 1..J
    T_all   <- numeric(0)  # cumulative dose prior exposure (model scale)
    pid_all <- integer(0)  # patient id for each observation
    
    patient_seq <- list()
    patient_id  <- 0L
    stop_trial  <- FALSE
    
    # cohort starting dose rule trackers
    j_H <- 0L
    j_S_prev <- 1L
    j_S_curr <- 1L
    
    for (coh in 1:ncohort) {
      
      if (coh > 1 && length(y_all) > 0) {
        res <- estimate_MTD_JAGS(
          y = y_all,
          d_level = d_all,
          p = p_model,
          Tcum = T_all,
          pid = pid_all,
          TARGET = TARGET,
          cutoff = c_stop,
          model_file = model_file,
          T_ref_for_curve = T_ref_for_curve
        )
        j_MTD <- res$MTD
        
        if (j_MTD > j_H) {
          j_S_curr <- min(j_H + 1L, J)
        } else if (j_MTD < j_S_prev) {
          j_S_curr <- max(j_S_prev - 1L, 1L)
        } else {
          j_S_curr <- j_MTD
        }
      }
      
      j_S_prev <- j_S_curr
      j_H <- max(j_H, j_S_curr)
      
      # enroll COHORTSIZE new patients
      for (k in seq_len(COHORTSIZE)) {
        patient_id <- patient_id + 1L
        patient_seq[[patient_id]] <- list(
          id = patient_id, cohort = coh,
          doses = integer(0), tox = integer(0),
          stop = NA_character_
        )
      }
      pid_local <- (patient_id - COHORTSIZE + 1L):patient_id
      
      # within-cohort states
      curr_dose <- rep(j_S_curr, COHORTSIZE)  # dose level
      active    <- rep(TRUE,  COHORTSIZE)
      ncycle    <- rep(0L, COHORTSIZE)
      # Carryover state: previous cycle's *effective* toxicity probability c_{k-1}
      c_prev <- rep(NA_real_, COHORTSIZE)
      
      # cumulative dose PRIOR exposure on model scale (numeric)
      # This is what generates Tij
      Tprior <- rep(0, COHORTSIZE)
      
      repeat {
        idx <- which(active & (ncycle < K))
        if (length(idx) == 0) break
        
        d_cycle  <- curr_dose[idx]          # dose level indices
        dose_amt <- p_model[d_cycle]        # model scale numeric dose
        
        # ---- FIXED Tij ----
        # Tij for this observation is cumulative dose PRIOR to current cycle:
        T_cycle <- Tprior[idx]
        
        # --------------------------------------------
        # Carryover toxicity data generation (IPCRM)
        # c_k = min( PI[d_k] + gamma * c_{k-1}, 1 )
        # c_1 = PI[d_1]
        # --------------------------------------------
        pi_raw <- PI[d_cycle]
        
        p_eff <- ifelse(
          is.na(c_prev[idx]),
          pi_raw,
          pmin(pi_raw + gamma * c_prev[idx], 1)
        )
        
        y_cycle <- rbinom(length(idx), 1, p_eff)
        
        # update carryover state AFTER generating this cycle outcome
        c_prev[idx] <- p_eff
        
        # append observation-level data
        y_all   <- c(y_all, y_cycle)
        d_all   <- c(d_all, d_cycle)
        T_all   <- c(T_all, T_cycle)
        pid_all <- c(pid_all, pid_local[idx])
        
        # update patient sequence debug
        for (m in seq_along(idx)) {
          pg <- pid_local[idx[m]]
          patient_seq[[pg]]$doses <- c(patient_seq[[pg]]$doses, d_cycle[m])
          patient_seq[[pg]]$tox   <- c(patient_seq[[pg]]$tox,   y_cycle[m])
          if (y_cycle[m] == 1L && is.na(patient_seq[[pg]]$stop)) {
            patient_seq[[pg]]$stop <- "toxicity"
          }
        }
        
        # update counters AFTER recording
        ncycle[idx] <- ncycle[idx] + 1L
        Tprior[idx] <- Tprior[idx] + dose_amt
        
        # stop those with toxicity
        idx_tox <- idx[y_cycle == 1L]
        if (length(idx_tox) > 0) active[idx_tox] <- FALSE
        
        # update highest dose tried
        j_H <- max(j_H, max(d_cycle))
        
        # if no active patients left, stop
        if (!any(active)) break
        
        # update estimated MTD using all accumulated data
        res <- estimate_MTD_JAGS(
          y = y_all,
          d_level = d_all,
          p = p_model,
          Tcum = T_all,
          pid = pid_all,
          TARGET = TARGET,
          cutoff = c_stop,
          model_file = model_file,
          T_ref_for_curve = T_ref_for_curve
        )
        j_MTD <- res$MTD
        
        # decide next cycle for active patients
        idx2 <- which(active & (ncycle < K))
        if (length(idx2) == 0) break
        
        # escalation rule (your style): no "stay"; escalate only if allowed
        can_escalate <- (curr_dose[idx2] < J) & (curr_dose[idx2] < j_MTD)
        
        # stop if cannot escalate
        if (any(!can_escalate)) {
          idx_stop <- idx2[!can_escalate]
          for (ii in idx_stop) {
            pg <- pid_local[ii]
            if (is.na(patient_seq[[pg]]$stop)) patient_seq[[pg]]$stop <- "no_escalation"
          }
          active[idx_stop] <- FALSE
        }
        
        # escalate by 1 if allowed
        if (any(can_escalate)) {
          idx_up <- idx2[can_escalate]
          curr_dose[idx_up] <- pmin(curr_dose[idx_up] + 1L, J)
          j_H <- max(j_H, max(curr_dose[idx_up]))
        }
        
        # stop if reached max cycles
        idx_done <- which(active & (ncycle >= K))
        if (length(idx_done) > 0) {
          for (ii in idx_done) {
            pg <- pid_local[ii]
            if (is.na(patient_seq[[pg]]$stop)) patient_seq[[pg]]$stop <- "maxK"
          }
          active[idx_done] <- FALSE
        }
      } # repeat within cohort
      # end cohorts
    }
    # label remaining NA stops
    for (i in seq_along(patient_seq)) {
      if (is.na(patient_seq[[i]]$stop)) patient_seq[[i]]$stop <- "trial_end"
    }
    patient_seq_by_trial[[trial]] <- patient_seq
    
    if (stop_trial) {
      nstop <- nstop + 1L
    } else if (length(y_all) > 0) {
      res <- estimate_MTD_JAGS(
        y = y_all,
        d_level = d_all,
        p = p_model,
        Tcum = T_all,
        pid = pid_all,
        TARGET = TARGET,
        cutoff = c_stop,
        model_file = model_file,
        T_ref_for_curve = T_ref_for_curve
      )
      
      final_MTD <- select_dose_closest(res$posttox, TARGET)
      dose.select[final_MTD] <- dose.select[final_MTD] + 1
    }
    
    # cycle-level summaries
    for (j in seq_len(J)) {
      nobs[j]     <- nobs[j]     + sum(d_all == j)
      ntox_obs[j] <- ntox_obs[j] + sum(y_all[d_all == j] == 1L)
    }
  } # end trials
  
  list(
    sel_pct = dose.select / ntrial,
    stop_pct = nstop / ntrial,
    avg_cycles_at_dose = nobs / ntrial,
    avg_tox_events_at_dose = ntox_obs / ntrial,
    patient_seq_by_trial = patient_seq_by_trial,
    model_file = model_file,
    T_ref_for_curve = T_ref_for_curve
  )
}


#mu_beta - mean of prior 
backsol <- function(ske, mu_beta0, mu_beta1){
  dose = (log(ske/(1-ske)) - mu_beta0)/mu_beta1
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Please provide the job index as the first argument.")
job_i <- as.integer(args[1])
ske1 = c(0.02, 0.12, 0.3, 0.5, 0.65)
#backsolve for d_j
dose = backsol(ske1, mu_beta0 = 3, mu_beta1 = 1)
model_file = 'logit_CO.bug'
gamma = 0.2
# sce1 = c(0.02, 0.05, 0.08, 0.1, 0.3)
sce3 = c(0.05, 0.15, 0.3, 0.5, 0.8)
# sce1 = c(0.02,0.05,0.1,0.3,0.5)
cat('seed', job_i, '\n')
res <- IPCRM(
              sce3,
              TARGET = 0.3,
              dose,
              COHORTSIZE = 3,
              ncohort = 3,
              ntrial = 1,
              K = 3,
              c_stop = 0.96,
              seed = job_i,
              model_file = model_file,
              gamma = gamma
)
foldername = 'IPCRM/res3'
if (!dir.exists(paste0('results/',foldername))) {
  dir.create(paste0('results/',foldername), recursive = TRUE)}
saveRDS(res, paste0('results/',foldername,'/trial-',job_i))

sce1 = c(0.3, 0.5, 0.7, 0.8, 0.9)
res <- IPCRM(
  sce1,
  TARGET = 0.3,
  dose,
  COHORTSIZE = 3,
  ncohort = 3,
  ntrial = 1,
  K = 3,
  c_stop = 0.96,
  seed = job_i,
  model_file = model_file,
  gamma = gamma
)
foldername = 'IPCRM/res1'
if (!dir.exists(paste0('results/',foldername))) {
  dir.create(paste0('results/',foldername), recursive = TRUE)}
saveRDS(res, paste0('results/',foldername,'/trial-',job_i))
sce2 = c(0.1, 0.3, 0.55, 0.65, 0.75)
res <- IPCRM(
  sce2,
  TARGET = 0.3,
  dose,
  COHORTSIZE = 3,
  ncohort = 3,
  ntrial = 1,
  K = 3,
  c_stop = 0.96,
  seed = job_i,
  model_file = model_file,
  gamma = gamma
)
foldername = 'IPCRM/res2'
if (!dir.exists(paste0('results/',foldername))) {
  dir.create(paste0('results/',foldername), recursive = TRUE)}
saveRDS(res, paste0('results/',foldername,'/trial-',job_i))
sce4 = c(0.02, 0.05, 0.1, 0.3, 0.5)
res <- IPCRM(
  sce4,
  TARGET = 0.3,
  dose,
  COHORTSIZE = 3,
  ncohort = 3,
  ntrial = 1,
  K = 3,
  c_stop = 0.96,
  seed = job_i,
  model_file = model_file,
  gamma = gamma
)
foldername = 'IPCRM/res4'
if (!dir.exists(paste0('results/',foldername))) {
  dir.create(paste0('results/',foldername), recursive = TRUE)}
saveRDS(res, paste0('results/',foldername,'/trial-',job_i))
sce5 = c(0.02, 0.05, 0.08, 0.1, 0.3)
res <- IPCRM(
  sce5,
  TARGET = 0.3,
  dose,
  COHORTSIZE = 3,
  ncohort = 3,
  ntrial = 1,
  K = 3,
  c_stop = 0.96,
  seed = job_i,
  model_file = model_file,
  gamma = gamma
)
foldername = 'IPCRM/res5'
if (!dir.exists(paste0('results/',foldername))) {
  dir.create(paste0('results/',foldername), recursive = TRUE)}
saveRDS(res, paste0('results/',foldername,'/trial-',job_i))
sce6 = c(0.6, 0.7, 0.75, 0.8, 0.85)
res <- IPCRM(
  sce6,
  TARGET = 0.3,
  dose,
  COHORTSIZE = 3,
  ncohort = 3,
  ntrial = 1,
  K = 3,
  c_stop = 0.96,
  seed = job_i,
  model_file = model_file,
  gamma = gamma
)
foldername = 'IPCRM/res6'
if (!dir.exists(paste0('results/',foldername))) {
  dir.create(paste0('results/',foldername), recursive = TRUE)
}
saveRDS(res, paste0('results/',foldername,'/trial-',job_i))
