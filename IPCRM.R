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
    pi = exp(p[d[i]] * beta1 + beta0)/(1+exp(p[d[i]] * beta1 + beta0))
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
  diff = abs(posttox-TARGET);
  dose.best = min(which(diff==min(diff)))
  return(dose.best)
}

estimate_MTD_JAGS <- function(y, d, p,
                              TARGET = 0.3,
                              cutoff = 0.96,
                              model_file = "logit.bug",
                              n.adapt = 1000,
                              n.burn  = 1000,
                              n.iter  = 2000) {
  
  # d is integer dose level per observation; p is dtilde vector length J
  data_jags <- list(
    N     = length(y),
    y_bin = as.integer(y),
    d     = as.numeric(p[d])  # pass effective dose value for each observation
  )
  
  jags <- rjags::jags.model(
    file    = model_file,
    data    = data_jags,
    n.adapt = n.adapt,
    quiet   = TRUE
  )
  
  update(jags, n.burn, progress.bar = "none")
  
  samples <- coda::as.mcmc(do.call(rbind, coda.samples(
    model          = jags,
    variable.names = c("beta0", "beta1"),
    n.iter         = n.iter,
    progress.bar   = "none"
  )))
  
  draws <- as.matrix(samples)
  b0 <- draws[, "beta0"]
  b1 <- draws[, "beta1"]
  
  J <- length(p)
  
  # posterior mean tox probs: pi_hat[j] = E[ plogis(b0 + b1 * p[j]) | data ]
  posttox <- vapply(seq_len(J), function(j) {
    mean(plogis(b0 + b1 * p[j]))
  }, numeric(1))
  
  dose.best <- which.min(abs(posttox - TARGET))
  
  # early stopping probability: Pr(pi1 > TARGET | data)
  pi1_draws <- plogis(b0 + b1 * p[1])
  prob_overtox <- mean(pi1_draws > TARGET)
  
  stop <- as.integer(prob_overtox > cutoff)
  list(
    MTD = dose.best,
    posttox = posttox,
    prob_overtox = prob_overtox,
    stop = stop
  )
}


# ---- PLACEHOLDER: compute Pr(pi1 > TARGET | data) for early stopping ----
overtox_prob <- function(y, d, p, TARGET, lambda, sigma, mu) {
  q <- qlogis(TARGET)   # logit(phi)
  p1 <- p[1]            # effective dose at level 1
  
  # denominator Z(beta1) = ∫ k(b0,b1) db0
  browser()
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

IPCRM <- function(
    PI,
    TARGET = 0.3,
    p,
    COHORTSIZE = 3,
    ncohort = 10,
    ntrial,
    K = 3,
    c_stop = 0.96,
    seed = 6
){
  stopifnot(length(PI) == length(p))
  stopifnot(COHORTSIZE >= 1, ncohort >= 1, ntrial >= 1, K >= 1)
  set.seed(seed)
  
  J <- length(p)
  
  # summaries across trials
  dose.select <- numeric(J)   # final selected MTD counts
  ntox        <- rep(0, J)    # patient-level tox counts at each dose (count pt once per dose)
  ntrted      <- rep(0, J)    # patient-level treated counts at each dose (count pt once per dose)
  nstop       <- 0            # early stopped trials
  
  # for debugging: store patient sequences for EACH trial (can be large!)
  patient_seq_by_trial <- vector("list", ntrial)
  
  t.start <- Sys.time()
  
  for(trial in 1:ntrial) {
    
    # running observation-level data across entire trial
    y_all <- integer(0)  # toxicity outcomes per observation (cycle-dose)
    d_all <- integer(0)  # dose level per observation
    
    # patient-level debug list for this trial
    patient_seq <- list()
    patient_id  <- 0L
    
    # tracking for cohort starting dose rule
    j_H <- 0L        # highest tried dose level so far
    j_S_prev <- 1L   # previous cohort starting dose (initialize 1)
    stop_trial <- FALSE
    
    # ---- Step 1: first cohort starts at dose 1 ----
    j_S_curr <- 1L
    
    for(coh in 1:ncohort) {
      
      # ---- Step 3: update starting dose for new cohort (for coh > 1) ----
      if(coh > 1) {
        res <- estimate_MTD_JAGS(y_all, d_all, p, TARGET, c_stop)
        j_MTD <- res$MTD
        
        if(j_MTD > j_H) {
          j_S_curr <- min(j_H + 1L, J)
        } else if(j_MTD < j_S_prev) {
          j_S_curr <- max(j_S_prev - 1L, 1L)
        } else {
          j_S_curr <- j_MTD
        }
      }
      
      j_S_prev <- j_S_curr
      
      # ---- Enroll COHORTSIZE new patients; all start at j_S_curr ----
      for(k in seq_len(COHORTSIZE)) {
        patient_id <- patient_id + 1L
        patient_seq[[patient_id]] <- list(
          id     = patient_id,
          cohort = coh,
          doses  = integer(0),
          tox    = integer(0),
          stop   = NA_character_
        )
      }
      pid <- (patient_id - COHORTSIZE + 1L):patient_id
      
      # Patient state variables within cohort
      curr_dose <- rep(j_S_curr, COHORTSIZE)
      active    <- rep(TRUE,  COHORTSIZE)  # FALSE once patient stops
      # count of DIFFERENT doses taken so far (unique-dose count)
      ndoses_taken <- rep(1L, COHORTSIZE)
      
      # update highest tried dose
      j_H <- max(j_H, j_S_curr)
      
      # ---- Step 2: intra-patient escalation by cycles ----
      repeat {
        
        # stop if nobody is active
        if(!any(active)) break
        
        # stop if all active patients already reached K different doses
        if(all(!active | (ndoses_taken >= K))) {
          # mark those still active but reached maxK
          idx_maxk <- which(active & ndoses_taken >= K)
          if(length(idx_maxk) > 0) {
            for(ii in idx_maxk) {
              p_global <- pid[ii]
              if(is.na(patient_seq[[p_global]]$stop)) patient_seq[[p_global]]$stop <- "maxK"
            }
          }
          break
        }
        
        # administer current cycle doses to active patients
        idx <- which(active)
        d_cycle <- curr_dose[idx]
        y_cycle <- rbinom(length(idx), 1, PI[d_cycle])
        
        # append to accumulated observation-level data
        y_all <- c(y_all, y_cycle)
        d_all <- c(d_all, d_cycle)
        
        # record patient-level sequences + toxicity stop
        for(m in seq_along(idx)) {
          p_global <- pid[idx[m]]
          patient_seq[[p_global]]$doses <- c(patient_seq[[p_global]]$doses, d_cycle[m])
          patient_seq[[p_global]]$tox   <- c(patient_seq[[p_global]]$tox,   y_cycle[m])
          
          if(y_cycle[m] == 1L && is.na(patient_seq[[p_global]]$stop)) {
            patient_seq[[p_global]]$stop <- "toxicity"
          }
        }
        
        # toxicity => stop further treatment for that patient
        active[idx[y_cycle == 1L]] <- FALSE
        
        # update highest tried dose based on administered cycle
        j_H <- max(j_H, max(d_cycle))
        
        # update estimated MTD using all data
        res <- estimate_MTD_JAGS(y_all, d_all, p, TARGET, c_stop)
        j_MTD <- res$MTD
        
        # decide next cycle dose for those still active and not yet at K different doses
        idx2 <- which(active & (ndoses_taken < K))
        if(length(idx2) == 0) break
        
        next_dose <- pmin(curr_dose[idx2] + 1L, J)
        
        # escalation availability in YOUR stated IPCRM stopping rule:
        # - if cannot escalate, patient stops (NOT "stay and continue")
        can_escalate <- (next_dose <= j_MTD) & (curr_dose[idx2] < J)
        
        # patients who cannot escalate => stop due to "no IPE available"
        if(any(!can_escalate)) {
          idx_stop <- idx2[!can_escalate]
          for(ii in idx_stop) {
            p_global <- pid[ii]
            if(is.na(patient_seq[[p_global]]$stop)) patient_seq[[p_global]]$stop <- "no_escalation"
          }
          active[idx_stop] <- FALSE
        }
        
        # those who can escalate: move up by 1 and increment number of different doses taken
        if(any(can_escalate)) {
          idx_up <- idx2[can_escalate]
          curr_dose[idx_up] <- next_dose[can_escalate]
          ndoses_taken[idx_up] <- ndoses_taken[idx_up] + 1L
          j_H <- max(j_H, max(curr_dose[idx_up]))
        }
        
      } # end cycles within cohort
      
      # ---- Early termination rule (trial-level): Pr(p1 > TARGET | data) > c_stop ----
      overtox <- estimate_MTD_JAGS(y_all, d_all, p, TARGET, c_stop)
      if(isTRUE(overtox$stop == 1)) {
        stop_trial <- TRUE
        break
      }
      
    } # end cohorts
    
    # finalize stop labels for any remaining NA (e.g., still active when trial ended)
    for(i in seq_along(patient_seq)) {
      if(is.na(patient_seq[[i]]$stop)) patient_seq[[i]]$stop <- "trial_end"
    }
    
    patient_seq_by_trial[[trial]] <- patient_seq
    
    if(stop_trial) {
      nstop <- nstop + 1
    } else {
      res <- estimate_MTD_JAGS(y_all, d_all, p, TARGET, c_stop)
      final_MTD <- res$MTD
      dose.select[final_MTD] <- dose.select[final_MTD] + 1
    }
    
    # -------------------------------
    # PATIENT-LEVEL summary updating
    # -------------------------------
    for(pid_i in seq_along(patient_seq)) {
      doses_i <- patient_seq[[pid_i]]$doses
      tox_i   <- patient_seq[[pid_i]]$tox
      if(length(doses_i) == 0) next
      
      # patient treated at these dose levels (count once per dose)
      dset <- unique(doses_i)
      ntrted[dset] <- ntrted[dset] + 1
      
      # patient had >=1 toxicity while at dose j (count once per dose)
      tox_doses <- unique(doses_i[tox_i == 1L])
      if(length(tox_doses) > 0) ntox[tox_doses] <- ntox[tox_doses] + 1
    }
  } # end trials
  
  t.end <- Sys.time()
  
  cat("Selection probability (%): ", round(dose.select / ntrial * 100, 2), "\n")
  cat("Avg # pts treated at dose: ", round(ntrted / ntrial, 3), "\n")
  cat("Avg # pts tox at dose:     ", round(ntox / ntrial, 3), "\n")
  cat("Early stop rate (%):       ", round(nstop / ntrial * 100, 2), "\n")
  print(t.end - t.start)
  
  invisible(list(
    sel_pct   = dose.select / ntrial,
    avg_trted = ntrted / ntrial,
    avg_tox   = ntox / ntrial,
    stop_pct  = nstop / ntrial,
    patient_seq_by_trial = patient_seq_by_trial
  ))
}

#mu_beta - mean of prior 
backsol <- function(ske, mu_beta0, mu_beta1){
  dose = (log(ske/(1-ske)) - mu_beta0)/mu_beta1
}
ske1 = c(0.02, 0.12, 0.3, 0.5, 0.65)
#backsolve for d_j
dose = backsol(ske1, mu_beta0 = 3, mu_beta1 = 1)
sce1 = c(0.05, 0.15, 0.3, 0.5, 0.8)
res <- IPCRM(
              sce1,
              TARGET = 0.3,
              dose,
              COHORTSIZE = 3,
              ncohort = 3,
              ntrial = 2000,
              K = length(dose),
              c_stop = 0.96,
              seed = 6
)
