draws_JAGS <- function(y, d, p, e,
                       model_file = "eff.bug",
                       n.chains = 3,
                       n.adapt  = 1000,
                       n.burn   = 2000,
                       n.iter   = 5000,
                       thin     = 2) {
  
  data_jags <- list(
    N     = length(y),
    y_bin = as.integer(y),
    e_bin = as.integer(e),
    dose_index = as.integer(d),
    D = length(p),
    d_j     = as.numeric(p[d])  # p is effective dose vector (dtilde)
  )
  
  jags <- rjags::jags.model(
    file    = model_file,
    data    = data_jags,
    n.chains = n.chains,
    n.adapt = n.adapt,
    quiet   = TRUE
  )
  
  update(jags, n.burn, progress.bar = "none")
  
  smp <- coda.samples(
    model          = jags,
    variable.names = c("beta0", "beta1", "gamma"),
    n.iter         = n.iter,
    thin           = thin,
    progress.bar   = "none"
  )
  draws <- as.matrix(smp)  # works for mcmc.list; stacks chains
  beta0 <- draws[, "beta0"]
  beta1 <- draws[, "beta1"]
  gamma <- draws[, 3:ncol(draws)]
  list(
    beta0 = beta0,
    beta1 = beta1,
    gamma = gamma
  )
}

get.MTD <- function(draws, T.TARGET,p){
  J = length(p)
  posttox <- vapply(seq_len(J), function(j) {
    mean(plogis(draws$beta0 + draws$beta1 * p[j]))
  }, numeric(1))
  diff <- abs(posttox - T.TARGET)
  dose.best <- which(diff == min(diff))[1]  # explicit tie rule
  return(list(MTD = dose.best, posttox = posttox))
}
overtox <- function(draws,T.TARGET,T.cutoff,p){
  prob_overtox = rep(0,length(p))
  for(j in 1:length(p)){
    pi1_draws <- plogis(draws$beta0 + draws$beta1 * p[j])
    prob_overtox[j] <- mean(pi1_draws > T.TARGET)
  }
  stop <- as.integer(prob_overtox > T.cutoff)
  return(stop)
}
futility <- function(draws, E.TARGET, E.cutoff, p){
  prob_futility = rep(0,length(p))
  for(j in 1:length(p)){
    pi1_draws <- pnorm(draws$gamma[,j])
    prob_futility[j] <- mean(pi1_draws < E.TARGET)
  }
  stop <- as.integer(prob_futility > E.cutoff)
  return(stop)
}

get.utility <- function(draws, uti, p){
  J = length(p)
  posttox <- vapply(seq_len(J), function(j) {
    mean(plogis(draws$beta0 + draws$beta1 * p[j]))
  }, numeric(1))
  posteff <- vapply(seq_len(J), function(j) {
    mean(pnorm(draws$gamma[,j]))
  }, numeric(1))
  postnotox <- 1 - posttox
  postnoeff <- 1 - posteff
  cat('posttox:', posttox, '\n')
  cat('posteff:', posteff, '\n')
  utility <- uti[1] * postnotox * posteff + uti[2] * postnotox * postnoeff + 
    uti[3] * posteff * posttox 
  return(utility)
}
# utility vector over doses -> average utility over all length-K trajectories -> pick start by softmax
select_traj_start_softmax <- function(draws, uti, p, K, admissible = NULL, tau = 1) {
  J <- length(p)
  
  # posterior mean utility per dose (your function)
  u <- get.utility(draws, uti, p)
  
  # default admissible: 1:J
  if (is.null(admissible)) admissible <- seq_len(J)
  
  # require consecutive admissible (your assumption)
  admissible <- sort(unique(admissible))
  if (!all(diff(admissible) == 1L)) stop("admissible must be consecutive indices.")
  
  topA <- max(admissible)
  if (K > length(admissible)) stop("K > length(admissible).")
  
  # all possible starting points (contiguous windows of length K)
  starts <- seq(from = min(admissible), to = topA - K + 1L, by = 1L)
  cat('starts:', starts, '\n')
  # average utility of each trajectory/window
  traj_avg <- vapply(starts, function(s) mean(u[s:(s + K - 1L)]), numeric(1))
  cat('traj_avg:', traj_avg, '\n' )
  # softmax with temperature tau (stable)
  z <- (traj_avg - max(traj_avg)) / tau
  prob <- exp(z)
  prob <- prob / sum(prob)
  
  # sample starting point
  start_sel <- sample(starts, size = 1, prob = prob)
  
  list(
    start = start_sel,
    starts = starts,
    traj_avg = traj_avg,
    prob = prob,
    utility = u
  )
}

select_OBD <- function(utility, dose_trt) {
  treated <- sort(unique(dose_trt))
  treated[which.max(utility[treated])]
}


IPCRM_TS <- function(
    T.PI,
    E.PI,
    T.TARGET = 0.3,
    E.TARGET = 0.2,
    p,
    COHORTSIZE = 3,
    ncohort = 10,
    s1.coh = 3,
    ntrial,
    K = 3,
    T.cutoff = 0.9,
    E.cutoff = 0.9,
    seed = 6
){
  stopifnot(length(T.PI) == length(p))
  stopifnot(COHORTSIZE >= 1, ncohort >= 1, ntrial >= 1, K >= 1)
  set.seed(seed)
  
  J <- length(p)
  levels <- 1:J
  # summaries across trials
  dose.select <- numeric(J)   # final selected MTD counts
  ntox        <- rep(0, J)    # patient-level tox counts at each dose (count pt once per dose)
  ntrted      <- rep(0, J)    # patient-level treated counts at each dose (count pt once per dose)
  nstop       <- 0            # early stopped trials
  nobs     <- rep(0L, J)  # total dose-cycles at each dose across trials
  ntox_obs <- rep(0L, J)  # total tox events (per cycle) at each dose across trials
  
  # for debugging: store patient sequences for EACH trial (can be large!)
  patient_seq_by_trial <- vector("list", ntrial)
  
  t.start <- Sys.time()
  
  for(trial in 1:ntrial) {
    
    # running observation-level data across entire trial
    y_all <- integer(0)  # toxicity outcomes per observation (cycle-dose)
    d_all <- integer(0)  # dose level per observation
    e_all <- integer(0)
    # patient-level debug list for this trial
    patient_seq <- list()
    patient_id  <- 0L
    
    # tracking for cohort starting dose rule
    j_H <- 0L        # highest tried dose level so far
    j_S_prev <- 1L   # previous cohort starting dose (initialize 1)
    stop_trial <- FALSE
    stageI <- TRUE
    # ---- Step 1: first cohort starts at dose 1 ----
    j_S_curr <- 1L
    
    for(coh in 1:ncohort) {
      cat('coh', coh ,'start', '\n')
      if(coh == s1.coh){
        stageI <- FALSE
      }
      if(stageI == TRUE){
        # ---- Step 3: update starting dose for new cohort (for coh > 1) ----
        if(coh > 1) {
          
          
          draws <- draws_JAGS(y_all, d_all, p, e_all)
          res <- get.MTD(draws, T.TARGET,p)
          j_MTD <- res$MTD
          cat("posttox=", paste(round(res$posttox,3), collapse=" "), " j_MTD=", j_MTD, "\n")
          if(j_MTD > j_H) {
            j_S_curr <- min(j_H + 1L, J)
          } else if(j_MTD < j_S_prev) {
            j_S_curr <- max(j_S_prev - 1L, 1L)
          } else {
            j_S_curr <- j_MTD
          }
        }
        cat('highest dose tried', j_H, '\n')
        cat('previous dose', j_S_prev, '\n')
        cat('curr dose', j_S_curr, '\n')
        j_S_prev <- j_S_curr
        
        # ---- Enroll COHORTSIZE new patients; all start at j_S_curr ----
        for(k in seq_len(COHORTSIZE)) {
          patient_id <- patient_id + 1L
          patient_seq[[patient_id]] <- list(
            id     = patient_id,
            cohort = coh,
            doses  = integer(0),
            tox    = integer(0),
            eff    = integer(0),
            stop   = NA_character_
          )
        }
        pid <- (patient_id - COHORTSIZE + 1L):patient_id
        
        # Patient state variables within cohort
        curr_dose <- rep(j_S_curr, COHORTSIZE)
        active    <- rep(TRUE,  COHORTSIZE)  # FALSE once patient stops
        # count of DIFFERENT doses taken so far (unique-dose count)
        # ndoses_taken <- rep(1L, COHORTSIZE)
        # NEW: count cycles received (per patient)
        ncycle <- rep(0L, COHORTSIZE)
        
        # update highest tried dose
        j_H <- max(j_H, j_S_curr)
        
        # ---- Step 2: intra-patient escalation by cycles ----
        repeat {
          cat('Intra-patient coh', coh, '\n')
          # stop if nobody is active
          if(!any(active)) break
          
          # stop if all active patients already reached K cycles
          if(all(!active | (ncycle >= K))) {
            idx_done <- which(active & ncycle >= K)
            if(length(idx_done) > 0) {
              for(ii in idx_done) {
                p_global <- pid[ii]
                if(is.na(patient_seq[[p_global]]$stop)) patient_seq[[p_global]]$stop <- "maxK"
              }
              active[idx_done] <- FALSE
            }
            break
          }
          
          # administer current cycle doses to active patients
          idx <- which(active & (ncycle < K))
          # idx <- which(active)
          if(length(idx) == 0) break
          # cat("d_cycle:", d_cycle, "\n")
          # cat("PI[d_cycle]:", PI[d_cycle], "\n")
          # cat("y_cycle:", y_cycle, "\n")
          d_cycle <- curr_dose[idx]
          y_cycle <- rbinom(length(idx), 1, T.PI[d_cycle])
          e_cycle <- rbinom(length(idx), 1, E.PI[d_cycle])
          
          cat("d_cycle:", d_cycle, "\n")
          cat("PI[d_cycle]:", T.PI[d_cycle], "\n")
          cat("y_cycle:", y_cycle, "\n")
          cat("e_cycle:", e_cycle, "\n")
          # update cycle counters
          ncycle[idx] <- ncycle[idx] + 1L
          
          # append to accumulated observation-level data
          y_all <- c(y_all, y_cycle)
          d_all <- c(d_all, d_cycle)
          e_all <- c(e_all, e_cycle)
          # record patient-level sequences + toxicity stop
          for(m in seq_along(idx)) {
            p_global <- pid[idx[m]]
            patient_seq[[p_global]]$doses <- c(patient_seq[[p_global]]$doses, d_cycle[m])
            patient_seq[[p_global]]$tox   <- c(patient_seq[[p_global]]$tox,   y_cycle[m])
            patient_seq[[p_global]]$eff   <- c(patient_seq[[p_global]]$eff,   e_cycle[m])
            
            if(y_cycle[m] == 1L && is.na(patient_seq[[p_global]]$stop)) {
              patient_seq[[p_global]]$stop <- "toxicity"
            }
          }
          
          # toxicity => stop further treatment
          idx_tox <- idx[y_cycle == 1L]
          if(length(idx_tox) > 0) active[idx_tox] <- FALSE
          
          # update highest tried dose
          j_H <- max(j_H, max(d_cycle))
          
          # update estimated MTD using all data
          
          draws <- draws_JAGS(y_all, d_all, p, e_all)
          res <- get.MTD(draws, T.TARGET,p)
          j_MTD <- res$MTD
          cat('y_all', y_all, '\n')
          cat('d_all', d_all, '\n')
          cat('e_all', e_all, '\n')
          cat("intra posttox=", paste(round(res$posttox,3), collapse=" "), " j_MTD=", j_MTD, "\n")
          # decide next cycle dose for those still active and still have cycles left
          idx2 <- which(active & (ncycle < K))
          if(length(idx2) == 0) break
          
          next_dose <- pmin(curr_dose[idx2] + 1L, J)
          
          
          # Escalation availability rule:
          # - must have a higher dose available (curr < J)
          # - and escalation must be allowed by the model (next_dose <= j_MTD)
          can_escalate <- (curr_dose[idx2] < J) & (curr_dose[idx2] <= j_MTD)
          
          # If cannot escalate => end treatment (NO "stay")
          if(any(!can_escalate)) {
            idx_stop <- idx2[!can_escalate]
            for(ii in idx_stop) {
              p_global <- pid[ii]
              if(is.na(patient_seq[[p_global]]$stop)) patient_seq[[p_global]]$stop <- "no_escalation"
            }
            active[idx_stop] <- FALSE
          }
          
          # If can escalate => move up by 1 for next cycle
          if(any(can_escalate)) {
            idx_up <- idx2[can_escalate]
            curr_dose[idx_up] <- curr_dose[idx_up] + 1L
            j_H <- max(j_H, max(curr_dose[idx_up]))
          }
          
          
          # after escalation attempt, immediately stop anyone who just reached K cycles
          idx_done2 <- which(active & (ncycle >= K))
          if(length(idx_done2) > 0) {
            for(ii in idx_done2) {
              p_global <- pid[ii]
              if(is.na(patient_seq[[p_global]]$stop)) patient_seq[[p_global]]$stop <- "maxK"
            }
            active[idx_done2] <- FALSE
          }
          
          
          cat('can escalation', can_escalate, '\n')
          cat('active', active, '\n')
          
        }
        
        
        # ---- Early termination rule (trial-level): Pr(p1 > TARGET | data) > c_stop ----
        tox <- overtox(draws,T.TARGET,T.cutoff,p)
        if(tox[1] == 1){
          stop_trial <- TRUE
          break
        }
        
      }
      else{
        #determine admissible dose
        draws <- draws_JAGS(y_all, d_all, p, e_all)
        res <- get.MTD(draws, T.TARGET,p)
        tox <- overtox(draws,T.TARGET,T.cutoff,p)
        fut <- futility(draws, E.TARGET, E.cutoff, p)
        MTD <- res$MTD
        good <- which(tox == 0 & fut == 0)
        if (length(good) == 0) { stop_trial <- TRUE; break }
        start_l <- min(good)
        admissible <- start_l:MTD
        cat('admissible:', admissible, '\n')
        if(length(admissible) > K){
          cat('> K', '\n')
          utility <- get.utility(draws, uti, p)
          cat('utility:', utility, '\n')
          cand <- select_traj_start_softmax(draws, uti, p, K, admissible = admissible, tau = 1)
          start <- cand$start  # from my function
          admissible <- start:(start + K - 1L) 
          cat('admissible',admissible,'\n')
        }
        j_curr <- admissible[1]
        topA <- max(admissible)
        # ---- Enroll COHORTSIZE new patients; all start at j_S_curr ----
        for(k in seq_len(COHORTSIZE)) {
          patient_id <- patient_id + 1L
          patient_seq[[patient_id]] <- list(
            id     = patient_id,
            cohort = coh,
            doses  = integer(0),
            tox    = integer(0),
            eff    = integer(0),
            stop   = NA_character_
          )
        }
        pid <- (patient_id - COHORTSIZE + 1L):patient_id
        # Patient state variables within cohort
        curr_dose <- rep(j_curr, COHORTSIZE)
        active    <- rep(TRUE,  COHORTSIZE)  # FALSE once patient stops
        repeat {
          # stop if nobody is active
          if (!any(active)) break
          
          # administer current-cycle doses to active patients
          idx <- which(active)
          if (length(idx) == 0) break
          
          d_cycle <- curr_dose[idx]
          y_cycle <- rbinom(length(idx), 1, T.PI[d_cycle])
          e_cycle <- rbinom(length(idx), 1, E.PI[d_cycle])
          
          # append accumulated observation-level data
          y_all <- c(y_all, y_cycle)
          d_all <- c(d_all, d_cycle)
          e_all <- c(e_all, e_cycle)
          cat('y_all', y_all, '\n')
          cat('d_all', d_all, '\n')
          cat('e_all', e_all, '\n')
          # record patient-level sequences + toxicity stop
          for (m in seq_along(idx)) {
            p_global <- pid[idx[m]]
            patient_seq[[p_global]]$doses <- c(patient_seq[[p_global]]$doses, d_cycle[m])
            patient_seq[[p_global]]$tox   <- c(patient_seq[[p_global]]$tox,   y_cycle[m])
            patient_seq[[p_global]]$eff   <- c(patient_seq[[p_global]]$eff,   e_cycle[m])
            
            if (y_cycle[m] == 1L && is.na(patient_seq[[p_global]]$stop)) {
              patient_seq[[p_global]]$stop <- "toxicity"
            }
          }
          
          # toxicity => stop further treatment
          idx_tox <- idx[y_cycle == 1L]
          if (length(idx_tox) > 0) active[idx_tox] <- FALSE
          
          # among those who remained non-toxic this cycle, if they were treated at topA, they finished trajectory
          idx_finish <- idx[(y_cycle == 0L) & (d_cycle == topA)]
          if (length(idx_finish) > 0) {
            for (pp in idx_finish) {
              p_global <- pid[pp]
              if (is.na(patient_seq[[p_global]]$stop)) patient_seq[[p_global]]$stop <- "topAdmissible"
            }
            active[idx_finish] <- FALSE
          }
          
          # escalate remaining active by 1 level, capped at topA
          idx2 <- which(active)
          if (length(idx2) == 0) break
          
          curr_dose[idx2] <- pmin(curr_dose[idx2] + 1L, topA)
        }
      }
    }
    
    # end cohorts
    
    # finalize stop labels for any remaining NA (e.g., still active when trial ended)
    for(i in seq_along(patient_seq)) {
      if(is.na(patient_seq[[i]]$stop)) patient_seq[[i]]$stop <- "trial_end"
    }
    
    patient_seq_by_trial[[trial]] <- patient_seq
    
    if(stop_trial) {
      nstop <- nstop + 1
    } else {
      #select the OBD among tried doses
      draws <- draws_JAGS(y_all, d_all, p, e_all)
      utility <- get.utility(draws, uti, p)
      final_OBD = select_OBD(utility, d_all)
      # final_MTD <- res$MTD
      dose.select[final_OBD] <- dose.select[final_OBD] + 1
      cat("utility=", paste(round(utility,3), collapse=" "), " final_OBD=", final_OBD, "\n")
      
      
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
    for(j in seq_len(J)) {
      nobs[j]     <- nobs[j]     + sum(d_all == j)
      ntox_obs[j] <- ntox_obs[j] + sum(y_all[d_all == j] == 1L)
    }
    
  } # end trials
  
  t.end <- Sys.time()
  
  cat("Selection probability (%): ", round(dose.select / ntrial * 100, 2), "\n")
  cat("Avg # pts treated at dose: ", round(nobs / ntrial, 3), "\n")
  cat("Avg # pts tox at dose:     ", round(ntox / ntrial, 3), "\n")
  cat("Early stop rate (%):       ", round(nstop / ntrial * 100, 2), "\n")
  print(t.end - t.start)
  
  invisible(list(
    sel_pct   = dose.select / ntrial,
    avg_trted = nobs / ntrial,
    avg_cycle = ntrted / ntrial,
    avg_tox   = ntox / ntrial,
    stop_pct  = nstop / ntrial,
    patient_seq_by_trial = patient_seq_by_trial
  ))
}