select_dose <- function(post_prob, dose_trt) {
  treated <- sort(unique(dose_trt))
  cand    <- treated[post_prob[treated] == max(post_prob[treated])]
  min(cand)
}
CRM <- function(
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
    
    # tracking for cohort starting dose rule
    j_H <- 0L        # highest tried dose level so far
    j_S_prev <- 1L   # previous cohort starting dose (initialize 1)
    stop_trial <- FALSE
    
    # ---- Step 1: first cohort starts at dose 1 ----
    j_curr <- 1L
    
    for(coh in 1:ncohort) {
      
      # ---- Step 3: update starting dose for new cohort (for coh > 1) ----
      if(coh > 1) {
        res <- estimate_MTD_JAGS(y_all, d_all, p, TARGET, c_stop)
        j_MTD <- res$MTD
        
        if(j_MTD > j_curr) {
          # j_S_curr <- min(j_H + 1L, J)
          j_curr <- min(j_curr + 1, J)
        } else if(j_MTD < j_curr) {
          j_curr <- max(j_curr - 1L, 1L)
        } else {
          j_curr <- j_MTD
        }
      }
      # ---- Enroll COHORTSIZE new patients; all start at j_S_curr ----
      y_coh <- rbinom(COHORTSIZE, 1, PI[j_curr])
      d_all <- c(d_all, rep(j_curr, COHORTSIZE))
      y_all <- c(y_all, y_coh)

      # ---- Early termination rule (trial-level): Pr(p1 > TARGET | data) > c_stop ----
      overtox <- estimate_MTD_JAGS(y_all, d_all, p, TARGET, c_stop)
      if(isTRUE(overtox$stop == 1)) {
        stop_trial <- TRUE
        break
      }
      
    } # end cohorts
    # cat('y_all', y_all, '\n')
    # cat('d_all', d_all, '\n')
    # finalize stop labels for any remaining NA (e.g., still active when trial ended)
    if(stop_trial) {
      nstop <- nstop + 1
    } else {
      res <- estimate_MTD_JAGS(y_all, d_all, p, TARGET, c_stop)
      # final_MTD <- select_dose(res$posttox, d_all)
      final_MTD <- res$MTD
      dose.select[final_MTD] <- dose.select[final_MTD] + 1
    }
    # -------------------------------
    # PATIENT-LEVEL summary updating
    # -------------------------------
    for(j in 1:J){
      ntrted[j] = ntrted[j] + sum(d_all == j)
    }
  } # end trials
  
  t.end <- Sys.time()
  
  cat("Selection probability (%): ", round(dose.select / ntrial * 100, 2), "\n")
  cat("Avg # pts treated at dose: ", round(ntrted / ntrial, 3), "\n")
  cat("Early stop rate (%):       ", round(nstop / ntrial * 100, 2), "\n")
  print(t.end - t.start)
  
  invisible(list(
    sel_pct   = dose.select / ntrial,
    avg_trted = ntrted / ntrial,
    stop_pct  = nstop / ntrial
  ))
}
backsol <- function(ske, mu_beta0, mu_beta1){
  dose = (log(ske/(1-ske)) - mu_beta0)/mu_beta1
}
ske1 = c(0.02, 0.12, 0.3, 0.5, 0.65)
#backsolve for d_j
dose = backsol(ske1, mu_beta0 = 3, mu_beta1 = 1)
sce1 = c(0.05, 0.15, 0.3, 0.5, 0.8)
res <- CRM(
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