# Adapted from Decision Analysis R code/Functions_cSTM_3state_time.R.
# Retains calculate_ce_out() and generate_psa_params() interfaces.
# Model calculation blocks are taken from the supplied/adapted DARTH model.
# The authors thank the DARTH workgroup for providing the coding framework in R.
# Original templates are unchanged. Source references: Tunnel_parameter_extraction.md.
#
# >>> CHANGED (ageing) --------------------------------------------------------
# This version matches the main model with ageing
# (ECD_DCD_tunnel_cSTM_ageing__V2.Rmd, sections 03.1-08):
#   - background mortality from the Statistics Canada life table [15]
#     (ages 75-109; stored in get_base_params() as v_q_pop)
#   - death rate in each state = population rate at current age + extra rate
#       r_dial_exc: extra rate on dialysis (replaces r_dial_dead)
#       r_G_exc:    extra rate with a working graft, month 4+ (replaces r_G_dead)
#   - early death (graft months 1-3) = hr_early x dialysis rate at that age
#   - death in the first year after graft loss never below dialysis
#   - Simpson's 1/3 within-cycle correction and template discounting, as in
#     the main model (previously: no correction, effective annual discounting)
# As before, r_dial_exc and r_G_exc are calibrated at the base case (55% and
# 59% alive at 5 years; check with calibrate_base_params()) and then varied as
# rates in the sensitivity analyses; PSA draws are not re-calibrated.
# Uncertainty rules: SE from the published 95% CI where one exists, (upper - lower) / 3.92
# (get_ci_se()); assumed 10% of the mean where none is published (r_GL_y1, r_GL_y2p,
# c_dial); utility SD 0.04 not used as an SE.
#------------------------------------------------------------------------------

get_base_params <- function() list(
  # >>> CHANGED (ageing): extra death rates on top of the general population.
  # Values come from the base-case calibration (same as the main model,
  # Steps 1A and 6); calibrate_base_params() reproduces them from the inputs.
  r_dial_exc = 0.08923827135004225,      # extra rate on dialysis: 55% alive at 5 years [1]
  hr_early   = 3,                        # early hazard ratio vs dialysis [1]
  r_GL_y1    = -log(0.933),              # graft loss, year 1 [3]
  r_GL_y2p   = -log(0.866/0.933)/2,      # graft loss, month 13+ [13]
  r_G_exc    = 0.041944239502720235,     # extra rate with a working graft: 59% alive at 5 years [1]
  p_GL_dead_annual = 0.502,              # death in the first year after graft loss [4]
  c_dial = 8928.97, c_G1to3 = 14133.18,
  c_G4to12 = 4747.60, c_G13p = 3272.01, c_tx_once = 61041.80,
  u_dial = 0.71, u_graft = 0.82,
  early_utility_fraction = 0,
  d_c = 0.015, d_e = 0.015, cycle_length = 1/12, n_cycles = 420,
  c_D = 0, u_D = 0,
  wcc_method = "Simpson1/3",             # >>> CHANGED (ageing): as in the main model (was "none")
  # >>> NEW (ageing): Statistics Canada Table 13-10-0114-01, 2022/2024 [15],
  # both sexes, qx (death probability between age x and x+1), ages 75-109.
  # Same values as Life_Table_StatsCan_75years.csv used by the main model.
  age_start = 75,
  v_q_pop = c(0.02405, 0.02663, 0.02951, 0.03275, 0.03638, 0.04045, 0.04503,
              0.05017, 0.05596, 0.06249, 0.06985, 0.07816, 0.08754, 0.09816,
              0.11017, 0.12379, 0.13885, 0.15504, 0.17234, 0.19072, 0.21188,
              0.23216, 0.25330, 0.27515, 0.29752, 0.32023, 0.34306, 0.36580,
              0.38824, 0.41017, 0.43143, 0.45185, 0.47129, 0.48967, 0.50690)
)

#------------------------------------------------------------------------------#
#### >>> NEW (ageing): run the model and return the cohort traces          ####
#------------------------------------------------------------------------------#
# Same steps as the main model, sections 03.1 (probabilities), 04.3 (arrays)
# and 05 (model run). Used by calculate_ce_out() and calibrate_base_params().
run_model_traces <- function(l_params_all, verbose = FALSE) {
  l_params_all <- modifyList(get_base_params(), as.list(l_params_all))
  l_params_all$n_cycles <- as.integer(round(l_params_all$n_cycles))
  list2env(l_params_all, envir = environment())
  if (cycle_length != 1/12) stop("Tunnel structure requires monthly cycles.")
  if (n_cycles < 1L || n_cycles > length(v_q_pop) * 12) stop("Invalid horizon.")
  
  ## Step 0: background mortality, one rate per month (each age repeated 12 times)
  v_r_pop   <- -log(1 - v_q_pop)
  v_r_pop_m <- rep(v_r_pop, each = 12)[1:n_cycles]
  v_p_pop_m <- rate_to_prob(r = v_r_pop_m, time = cycle_length)
  
  ## Step 1: dialysis = population + extra rate; after 12+ months of graft loss = dialysis
  v_r_dial_dead  <- v_r_pop_m + r_dial_exc
  v_p_dial_dead  <- rate_to_prob(r = v_r_dial_dead, time = cycle_length)
  v_p_GL13p_dead <- v_p_dial_dead
  
  ## Step 2: graft months 1-3 = hr_early x dialysis rate at that age
  v_p_G1to3_dead <- rate_to_prob(r = hr_early * v_r_dial_dead, time = cycle_length)
  
  ## Steps 3-4: graft loss among people alive with a graft (no age effect)
  p_GL_y1  <- rate_to_prob(r = r_GL_y1,  time = cycle_length)
  p_GL_y2p <- rate_to_prob(r = r_GL_y2p, time = cycle_length)
  
  ## Step 5: first year after graft loss, never below dialysis at the same age
  p_GL1to12_dead   <- 1 - (1 - p_GL_dead_annual)^cycle_length
  v_p_GL1to12_dead <- pmax(p_GL1to12_dead, v_p_dial_dead)
  
  ## Step 6 (applied): graft month 4+ = population + extra rate with a graft
  v_p_G4p_dead <- rate_to_prob(r = v_r_pop_m + r_G_exc, time = cycle_length)
  
  ## Transition probabilities: 13 tunnel positions x n_cycles months
  v_p_G1to3_GL  <- (1 - v_p_G1to3_dead) * p_GL_y1
  v_p_G4to12_GL <- (1 - v_p_G4p_dead)   * p_GL_y1
  v_p_G13p_GL   <- (1 - v_p_G4p_dead)   * p_GL_y2p
  m_p_G_dead  <- rbind(matrix(v_p_G1to3_dead, nrow = 3,  ncol = n_cycles, byrow = TRUE),
                       matrix(v_p_G4p_dead,   nrow = 10, ncol = n_cycles, byrow = TRUE))
  m_p_G_GL    <- rbind(matrix(v_p_G1to3_GL,  nrow = 3, ncol = n_cycles, byrow = TRUE),
                       matrix(v_p_G4to12_GL, nrow = 9, ncol = n_cycles, byrow = TRUE),
                       v_p_G13p_GL)
  m_p_GL_dead <- rbind(matrix(v_p_GL1to12_dead, nrow = 12, ncol = n_cycles, byrow = TRUE),
                       v_p_GL13p_dead)
  
  ## Tunnels and state names
  n_tunnel_size  <- nrow(m_p_G_dead)                                  # 13
  v_Graft_tunnel <- c(paste("G_",  seq(1, n_tunnel_size - 1), "Mo", sep = ""), paste("G_",  n_tunnel_size, "Mo+", sep = ""))
  v_GL_tunnel    <- c(paste("GL_", seq(1, n_tunnel_size - 1), "Mo", sep = ""), paste("GL_", n_tunnel_size, "Mo+", sep = ""))
  v_names_states_tunnels_tx <- c(v_Graft_tunnel, v_GL_tunnel, "Dead")
  n_states_tunnels_tx       <- length(v_names_states_tunnels_tx)
  v_names_states_dial       <- c("Dialysis", "Dead")
  v_names_cycles            <- paste("cycle", 0:n_cycles)
  
  ## Transition probability arrays (different in each cycle because of ageing)
  a_P_dial <- array(0, dim = c(2, 2, n_cycles),
                    dimnames = list(v_names_states_dial, v_names_states_dial,
                                    v_names_cycles[-length(v_names_cycles)]))
  a_P_dial["Dialysis", "Dialysis", ] <- 1 - v_p_dial_dead
  a_P_dial["Dialysis", "Dead", ]     <-     v_p_dial_dead
  a_P_dial["Dead", "Dead", ]         <- 1
  
  a_P_tunnels_tx <- array(0, dim = c(n_states_tunnels_tx, n_states_tunnels_tx, n_cycles),
                          dimnames = list(v_names_states_tunnels_tx, v_names_states_tunnels_tx,
                                          v_names_cycles[-length(v_names_cycles)]))
  # Part 1: graft tunnels G_1Mo-G_12Mo
  for (i in 1:(n_tunnel_size - 1)) {
    a_P_tunnels_tx[v_Graft_tunnel[i], v_Graft_tunnel[i + 1], ] <- 1 - m_p_G_dead[i, ] - m_p_G_GL[i, ]
    a_P_tunnels_tx[v_Graft_tunnel[i], "GL_1Mo", ]              <-     m_p_G_GL[i, ]
    a_P_tunnels_tx[v_Graft_tunnel[i], "Dead", ]                <-     m_p_G_dead[i, ]
  }
  # Part 2: graft holding state G_13Mo+
  a_P_tunnels_tx[v_Graft_tunnel[n_tunnel_size], v_Graft_tunnel[n_tunnel_size], ] <- 1 - m_p_G_dead[n_tunnel_size, ] - m_p_G_GL[n_tunnel_size, ]
  a_P_tunnels_tx[v_Graft_tunnel[n_tunnel_size], "GL_1Mo", ]                      <-     m_p_G_GL[n_tunnel_size, ]
  a_P_tunnels_tx[v_Graft_tunnel[n_tunnel_size], "Dead", ]                        <-     m_p_G_dead[n_tunnel_size, ]
  # Part 3: graft-loss tunnels GL_1Mo-GL_12Mo
  for (i in 1:(n_tunnel_size - 1)) {
    a_P_tunnels_tx[v_GL_tunnel[i], v_GL_tunnel[i + 1], ] <- 1 - m_p_GL_dead[i, ]
    a_P_tunnels_tx[v_GL_tunnel[i], "Dead", ]             <-     m_p_GL_dead[i, ]
  }
  # Part 4: graft-loss holding state GL_13Mo+
  a_P_tunnels_tx[v_GL_tunnel[n_tunnel_size], v_GL_tunnel[n_tunnel_size], ] <- 1 - m_p_GL_dead[n_tunnel_size, ]
  a_P_tunnels_tx[v_GL_tunnel[n_tunnel_size], "Dead", ]                     <-     m_p_GL_dead[n_tunnel_size, ]
  a_P_tunnels_tx["Dead", "Dead", ] <- 1
  
  # Checks (stop with an error if any draw gives invalid probabilities)
  check_transition_probability(a_P_dial,       verbose = verbose, err_stop = TRUE)
  check_transition_probability(a_P_tunnels_tx, verbose = verbose, err_stop = TRUE)
  check_sum_of_transition_array(a_P_dial,       n_states = 2,                   n_cycles = n_cycles, verbose = verbose, err_stop = TRUE)
  check_sum_of_transition_array(a_P_tunnels_tx, n_states = n_states_tunnels_tx, n_cycles = n_cycles, verbose = verbose, err_stop = TRUE)
  
  ## Cohort traces (section 05)
  m_M_dial <- matrix(0, nrow = n_cycles + 1, ncol = 2,
                     dimnames = list(0:n_cycles, v_names_states_dial))
  m_M_tunnels_tx <- matrix(0, nrow = n_cycles + 1, ncol = n_states_tunnels_tx,
                           dimnames = list(0:n_cycles, v_names_states_tunnels_tx))
  m_M_dial[1, "Dialysis"]    <- 1
  m_M_tunnels_tx[1, "G_1Mo"] <- 1
  for (t in 1:n_cycles) {
    m_M_dial      [t + 1, ] <- m_M_dial      [t, ] %*% a_P_dial      [, , t]
    m_M_tunnels_tx[t + 1, ] <- m_M_tunnels_tx[t, ] %*% a_P_tunnels_tx[, , t]
  }
  list(m_M_dial = m_M_dial, m_M_tunnels_tx = m_M_tunnels_tx, v_p_pop_m = v_p_pop_m,
       n_cycles = n_cycles, v_Graft_tunnel = v_Graft_tunnel, v_GL_tunnel = v_GL_tunnel)
}

#------------------------------------------------------------------------------#
#### Costs, QALYs and NMB for both strategies                                ####
#------------------------------------------------------------------------------#
# Same interface as before. >>> CHANGED (ageing): the model block now calls
# run_model_traces(). Rewards depend only on the state occupied, so
# trace %*% rewards gives the same totals as the transition-dynamics arrays
# in section 08 of the main model (and runs faster in OWSA/PSA).
calculate_ce_out <- function(l_params_all, n_wtp=50000, verbose=FALSE) {
  l_params_all <- modifyList(get_base_params(), as.list(l_params_all))
  # >>> NEW (ageing): Simpson's 1/3 rule needs an even number of cycles, so in
  # the OWSA on n_cycles (60-420, 100 values) the horizon is rounded to an
  # even number of months (e.g. 63.6 -> 64). 420 is even, so the base case is unchanged.
  if (l_params_all$wcc_method != "none")
    l_params_all$n_cycles <- 2 * round(l_params_all$n_cycles / 2)
  l_out <- run_model_traces(l_params_all, verbose = verbose)
  list2env(l_params_all, envir = environment())
  n_cycles    <- l_out$n_cycles
  v_names_str <- c("Continued dialysis", "ECD/DCD transplant")
  
  v_wcc <- if (wcc_method == "none") rep(1, n_cycles + 1) else
    gen_wcc(n_cycles = n_cycles, method = wcc_method)
  # >>> CHANGED (ageing): discount weights as in the main model / template
  v_dwc <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))
  v_dwe <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
  
  # 0 = dialysis utility during recovery; 1 = full graft utility from month 1.
  u_G1to3 <- u_dial + early_utility_fraction * (u_graft - u_dial)
  u_GL    <- u_dial
  c_GL    <- c_dial
  
  # State rewards (section 07): utilities per month, costs already per month
  v_u_dial <- c(Dialysis = u_dial, Dead = u_D) * cycle_length
  v_c_dial <- c(Dialysis = c_dial, Dead = c_D)
  v_u_tx   <- c(rep(u_G1to3, 3), rep(u_graft, 10), rep(u_GL, 13), u_D) * cycle_length
  v_c_tx   <- c(rep(c_G1to3, 3), rep(c_G4to12, 9), c_G13p, rep(c_GL, 13), c_D)
  
  # Discounted totals with the within-cycle correction (section 08)
  v_tot_qaly <- c(sum((l_out$m_M_dial       %*% v_u_dial) * v_dwe * v_wcc),
                  sum((l_out$m_M_tunnels_tx %*% v_u_tx)   * v_dwe * v_wcc))
  v_tot_cost <- c(sum((l_out$m_M_dial       %*% v_c_dial) * v_dwc * v_wcc),
                  sum((l_out$m_M_tunnels_tx %*% v_c_tx)   * v_dwc * v_wcc) + c_tx_once)
  
  v_nmb <- v_tot_qaly*n_wtp-v_tot_cost
  df_ce <- data.frame(Strategy=v_names_str, Cost=as.numeric(v_tot_cost),
                      Effect=as.numeric(v_tot_qaly), NMB=as.numeric(v_nmb))
  if(any(!is.finite(as.matrix(df_ce[,-1])))) stop("Non-finite model outcome.")
  return(df_ce)
}

#------------------------------------------------------------------------------#
#### >>> NEW (ageing): reproduce the base-case calibration                   ####
#------------------------------------------------------------------------------#
# Recalculates r_dial_exc (closed form) and r_G_exc (uniroot) from the life
# table and the 5-year targets, using the same model as calculate_ce_out().
# Use it to check the values stored in get_base_params().
calibrate_base_params <- function(l_params_all = get_base_params(),
                                  surv5_dial = 0.55, surv5_tx = 0.59) {
  l <- modifyList(get_base_params(), as.list(l_params_all))
  v_r_pop_m  <- rep(-log(1 - l$v_q_pop), each = 12)
  r_dial_exc <- (-log(surv5_dial) - sum(v_r_pop_m[1:60]) * l$cycle_length) / 5
  f_surv5_tx <- function(r) {
    out <- run_model_traces(modifyList(l, list(r_dial_exc = r_dial_exc, r_G_exc = r, n_cycles = 60)))
    1 - out$m_M_tunnels_tx[61, "Dead"]
  }
  r_G_exc <- uniroot(function(r) f_surv5_tx(r) - surv5_tx, interval = c(0, 1), tol = 1e-10)$root
  c(r_dial_exc = r_dial_exc, r_G_exc = r_G_exc)
}

#------------------------------------------------------------------------------#
####             Generate a PSA input parameter dataset                     ####
#------------------------------------------------------------------------------#
# Same seed, distribution families and data-frame interface as DARTH template.
# Assumed raw-scale SE is explicit; lognormal arguments are moment matched.
# Reported utility SD 0.04 is NOT automatically interpreted as an SE.
# Shared dialysis/graft-failure costs and utilities are derived in the model.
generate_psa_params <- function(n_sim=1000, seed=71818,
                                relative_se=0.10,
                                params_basecase=get_base_params()) {
  set.seed(seed)
  spec <- get_sa_spec(params_basecase, relative_se)
  df_psa <- as.data.frame(setNames(lapply(seq_len(nrow(spec)),function(i) {
    x <- spec[i,]
    switch(x$distribution,
           lognormal=rlnorm(n_sim,meanlog=x$meanlog,sdlog=x$sdlog),
           gamma=rgamma(n_sim,shape=x$shape,scale=x$scale),
           beta=rbeta(n_sim,shape1=x$alpha,shape2=x$beta))
  }),spec$pars))
  return(df_psa)
}

# Addition: expose assumptions and gamma/beta/lognormal arguments for review.
# (With ageing, the rates selected are r_dial_exc, r_GL_y1, r_GL_y2p and
#  r_G_exc; the life table v_q_pop is not varied.)
#------------------------------------------------------------------------------#
####  Published 95% CIs and the standard errors (SE) used for the PSA        ####
#------------------------------------------------------------------------------#
# Returns one row per PSA parameter: base value (mean), lower and upper limit,
# SE and its basis. Lower/upper are the SAME limits as the OWSA table
# (df_owsa_lit in the Rmd). Where a published 95% CI exists: SE = (upper - lower) / 3.92.
# Where it does not (r_GL_y1, r_GL_y2p, c_dial: ranges across studies):
# SE = relative_se x mean, and lower/upper are NA.
get_ci_se <- function(params_basecase = get_base_params(), relative_se = 0.10) {
  pars <- names(params_basecase)
  pars <- pars[grepl("^(r_|hr_|p_GL_dead_annual$|c_|u_)", pars)]
  pars <- pars[!pars %in% c("c_dead", "u_dead", "c_D", "u_D")]
  m <- setNames(unlist(params_basecase[pars], use.names = FALSE), pars)
  
  # Relative width of the 95% CI of the 2-year deceased-donor transplant cost [6]
  v_tx_ci_rel <- c(114287, 127956) / 121121
  
  # The two calibrated death rates: convert the CI of 5-year survival into a rate
  cal_dial  <- calibrate_base_params(params_basecase, surv5_dial = 0.54)   # 54% alive at 5 y
  cal_tx_hi <- calibrate_base_params(params_basecase, surv5_tx = 0.66)     # 66% alive at 5 y (lower rate)
  cal_tx_lo <- calibrate_base_params(params_basecase, surv5_tx = 0.51)     # 51% alive at 5 y (higher rate)
  
  ci <- data.frame(
    pars  = c("r_dial_exc", "r_G_exc", "hr_early", "p_GL_dead_annual",
              "u_dial", "u_graft", "c_G1to3", "c_G4to12", "c_G13p", "c_tx_once"),
    lower = c(m[["r_dial_exc"]],         unname(cal_tx_hi["r_G_exc"]), 1.50, 0.474,
              0.82 - 0.15, 0.74,
              m[["c_G1to3"]] * v_tx_ci_rel[1], m[["c_G4to12"]] * v_tx_ci_rel[1],
              m[["c_G13p"]]  * v_tx_ci_rel[1], m[["c_tx_once"]] * v_tx_ci_rel[1]),
    upper = c(unname(cal_dial["r_dial_exc"]), unname(cal_tx_lo["r_G_exc"]), 3.49, 0.533,
              0.82 - 0.08, 0.90,
              m[["c_G1to3"]] * v_tx_ci_rel[2], m[["c_G4to12"]] * v_tx_ci_rel[2],
              m[["c_G13p"]]  * v_tx_ci_rel[2], m[["c_tx_once"]] * v_tx_ci_rel[2]),
    stringsAsFactors = FALSE)
  
  j   <- match(pars, ci$pars)
  has <- !is.na(j)
  lower <- upper <- rep(NA_real_, length(pars))
  lower[has] <- ci$lower[j[has]]
  upper[has] <- ci$upper[j[has]]
  
  se    <- ifelse(has, (upper - lower) / 3.92, relative_se * m)
  basis <- ifelse(has, "SE from published 95% CI: (upper - lower) / 3.92",
                  paste0("Assumed SE = ", 100 * relative_se, "% of mean"))
  
  data.frame(pars = pars, mean = unname(m), lower = lower, upper = upper,
             se = unname(se), basis = basis, stringsAsFactors = FALSE)
}

get_sa_spec <- function(params_basecase = get_base_params(), relative_se = 0.10) {
  ci    <- get_ci_se(params_basecase, relative_se)
  pars  <- ci$pars
  m     <- ci$mean
  se    <- ci$se
  basis <- ci$basis
  
  # >>> CHANGED: rates (r_) and the hazard ratio (hr_) use gamma instead of
  # lognormal, so every PSA input is beta (utilities, probabilities) or gamma
  # (costs, rates, hazard ratio), as in the SA input table. The lognormal
  # branch below is kept so it can be switched back.
  family <- ifelse(grepl("^c_|^r_|^hr_", pars), "gamma",
                   ifelse(grepl("^u_|^p_", pars), "beta", "lognormal"))
  x <- data.frame(pars = pars, mean = m, se = se, distribution = family, basis = basis,
                  stringsAsFactors = FALSE)
  x$shape <- x$scale <- x$alpha <- x$beta <- x$meanlog <- x$sdlog <- NA_real_
  for (i in seq_len(nrow(x))) {
    if (family[i] == "gamma") {x$shape[i] <- (m[i]/se[i])^2; x$scale[i] <- se[i]^2/m[i]}
    if (family[i] == "beta") {
      k <- m[i]*(1 - m[i])/se[i]^2 - 1
      if (!is.finite(k) || k <= 0) stop("Invalid beta uncertainty: ", pars[i])
      x$alpha[i] <- m[i]*k; x$beta[i] <- (1 - m[i])*k
    }
    if (family[i] == "lognormal") {
      x$sdlog[i]   <- sqrt(log1p((se[i]/m[i])^2))
      x$meanlog[i] <- log(m[i]) - x$sdlog[i]^2/2
    }
  }
  x
}

get_owsa_ranges <- function(params_basecase=get_base_params()) {
  x <- get_sa_spec(params_basecase)
  result <- data.frame(pars=x$pars,min=pmax(1e-8,x$mean-1.96*x$se),
                       max=x$mean+1.96*x$se)
  b <- x$distribution=="beta"
  result$max[b] <- pmin(1-1e-8,result$max[b])
  j <- match("p_GL_dead_annual",result$pars)
  if(!is.na(j)) {result$min[j]<-.474;result$max[j]<-.533}
  # Settings are explored deterministically, not sampled in PSA.
  rbind(result,data.frame(pars=c("d_c","d_e","n_cycles"),
                          min=c(0,0,5/params_basecase$cycle_length),
                          max=c(.05,.05,params_basecase$n_cycles)))
}