##################################################
# Simulation: recoverability of directed (top-down vs bottom-up) connectivity among 
# Primary Auditory Cortex (PAC) / Superior Temporal Gyrus (STG) / Inferior Frontal 
# Gyrus (IFG) after 64-channel -> 10-PC dimension reduction, for the AVH network 
# analysis.
# Matthieu Vignes
# August 2026
##################################################
# Method:
# 1. Simulate 3 latent Region of Interest (ROI) sources (PAC, STG, IFG) via a 
#    lag-delayed VAR process. Two generating hypotheses:
#       - bottom_up: PAC -> STG -> IFG dominant (feedforward)
#       - top_down : IFG -> STG -> PAC dominant (feedback)
# 2. Project analytically through a 64-channel mixing matrix + sensor noise, 
#    reduce to 10 PCA components, and reconstruct 3 ROI estimates via a linear 
#    MMSE (Wiener) filter. Linear-Gaussian chain -> this reduces to a fixed 3x3 
#    "gain" matrix G (true ROI -> reconstructed ROI) plus a 3x3 reconstruction
#    noise covariance.
# 3. From the reconstructed ROI estimates, compute a lagged cross-correlation 
#    asymmetry index (fast proxy for directed influence / Granger-style asymmetry).
# 4. Test recoverability of an AVH-locked top-down shift (conditions 4/6, 
#    AVH-present vs AVH-absent epochs) across N and AVH event rate, and estimate 
#    the false-positive rate in control conditions (1/2/silence, no true shift 
#    expected).

set.seed(2026)
library(ggplot2)

# parameters
n_channels <- 64
n_components <- 10
T_samples <- 200 # -200 to 600 ms @ 250 Hz
lag_delay <- 3 # ~12 ms transmission lag between ROIs
own_rho <- 0.5 # per-node local persistence (lag-1)
cf <- 0.26 # dominant-direction coupling
cb <- 0.16 # reciprocal coupling (moderate asymmetry, not a toy-easy contrast)
innovation_sd_base <- 1
channel_noise_sd_base <- 80 # single-trial sensor noise (channel units); calibrated so power
                            # is informative (not floor/ceiling) across the N x AVH-rate grid
n_design_trials <- 40 # trial count assumed when the spatial filter was calibrated

# Node order: 1 = PAC, 2 = STG, 3 = IFG
cross_matrix <- function(hypothesis = c("bottom_up", "top_down")){
  hypothesis <- match.arg(hypothesis)
  M <- matrix(0, 3, 3) # M[target, source]
  if(hypothesis == "bottom_up"){
    M[2, 1] <- cf; M[3, 2] <- cf # PAC->STG, STG->IFG (feedforward, strong)
    M[1, 2] <- cb; M[2, 3] <- cb # STG->PAC, IFG->STG (feedback, weak)
  }else{
    M[1, 2] <- cf; M[2, 3] <- cf # STG->PAC, IFG->STG (feedback, strong = top-down)
    M[2, 1] <- cb; M[3, 2] <- cb # PAC->STG, STG->IFG (feedforward, weak residual)
  }
  M
}

## Fixed reduction/reconstruction
L <- matrix(rnorm(n_channels * 3), n_channels, 3) # random leadfield-like mixing, fixed
noise_var_design <- (channel_noise_sd_base^2) / n_design_trials
Cov_S_design <- diag(3)
Cov_Y_design <- L %*% Cov_S_design %*% t(L) + noise_var_design * diag(n_channels)

eig <- eigen(Cov_Y_design, symmetric = TRUE)
V <- eig$vectors[, 1:n_components, drop = FALSE] # 64 x 10 PCA basis
var_explained <- sum(eig$values[1:n_components]) / sum(eig$values)

Cov_SC <- t(L) %*% V # 3 x 10
Cov_C  <- t(V) %*% Cov_Y_design %*% V # 10 x 10
W <- Cov_SC %*% solve(Cov_C) # 3 x 10 Wiener weights
Recon <- W %*% t(V) # 3 x 64 combined operator

G <- Recon %*% L # 3 x 3 true -> reconstructed gain / leakage matrix
Sigma_eps_per_trial <- Recon %*% t(Recon) * channel_noise_sd_base^2 # 3x3 recon-noise cov, single-trial scale

rownames(G) <- colnames(G) <- c("PAC", "STG", "IFG")
cat("=== Variance explained by top", n_components, "of", n_channels, "channels ===\n")
cat(round(100 * var_explained, 1), "%\n\n")
cat("=== Reconstruction gain / cross-talk matrix G (rows = recovered, cols = true source) ===\n")
print(round(G, 3))
cat("\n(Diagonal near 1 = faithful recovery; off-diagonal = cross-talk/leakage between ROIs\n",
    "introduced by reducing 64 channels to", n_components, "components.)\n\n")

chol_Sigma <- function(sigma_mat) tryCatch(chol(sigma_mat), error = function(e) chol(sigma_mat + diag(1e-6, 3)))

# Batch generation
# Returns: the lagged-correlation directionality index D:
# - positive = feedforward/bottom-up dominant,
# - negative = feedback/top-down dominant
row_cor <- function(X, Y){
  # row-wise Pearson correlation between two (n_reps x T) matrices
  n <- ncol(X)
  mx <- rowMeans(X); my <- rowMeans(Y)
  sx <- sqrt(rowSums((X - mx)^2))
  sy <- sqrt(rowSums((Y - my)^2))
  rowSums((X - mx) * (Y - my)) / (sx * sy)
}

simulate_D_batch <- function(n_reps, hypothesis, n_trials_avg){
  Cross <- cross_matrix(hypothesis)
  eff_innov_sd  <- innovation_sd_base # coupling signal is fixed
  Sigma_eps <- Sigma_eps_per_trial / n_trials_avg  # sensor/reconstruction noise averages down
  U <- chol_Sigma(Sigma_eps) # for MVN noise draws: eps = t(U) %*% z
  # history buffer: list of (lag_delay+1) matrices, each 3 x n_reps
  hist <- vector("list", lag_delay + 1)
  for(k in seq_along(hist)) hist[[k]] <- matrix(rnorm(3 * n_reps, 0, eff_innov_sd), 3, n_reps)

  Shat_PAC <- matrix(0, n_reps, T_samples)
  Shat_STG <- matrix(0, n_reps, T_samples)
  Shat_IFG <- matrix(0, n_reps, T_samples)

  state_now <- hist[[lag_delay + 1]]
  for(t in seq_len(T_samples)){
    state_lag1 <- hist[[lag_delay + 1]]
    state_lagK <- hist[[1]]
    innov <- matrix(rnorm(3 * n_reps, 0, eff_innov_sd), 3, n_reps)
    new_state <- own_rho * state_lag1 + Cross %*% state_lagK + innov

    z <- matrix(rnorm(3 * n_reps), 3, n_reps)
    eps <- t(U) %*% z
    Shat_t <- G %*% new_state + eps

    Shat_PAC[, t] <- Shat_t[1, ]
    Shat_STG[, t] <- Shat_t[2, ]
    Shat_IFG[, t] <- Shat_t[3, ]

    hist <- c(hist[-1], list(new_state))
  }

  idx_early <- 1:(T_samples - lag_delay)
  idx_late  <- (lag_delay + 1):T_samples

  ff1 <- row_cor(Shat_PAC[, idx_early, drop = FALSE], Shat_STG[, idx_late, drop = FALSE]) -
         row_cor(Shat_STG[, idx_early, drop = FALSE], Shat_PAC[, idx_late, drop = FALSE])
  ff2 <- row_cor(Shat_STG[, idx_early, drop = FALSE], Shat_IFG[, idx_late, drop = FALSE]) -
         row_cor(Shat_IFG[, idx_early, drop = FALSE], Shat_STG[, idx_late, drop = FALSE])

  ff1 + ff2 # D: positive = feedforward/bottom-up dominant
}

# Main grid: AVH-locked shift in conditions 4 / 6
N_grid <- c(20, 30, 40, 50)
avh_rate_grid <- c(0.10, 0.20, 0.35)
n_trials_condition <- 60
n_sims <- 1000

main_results <- list(); row_i <- 1
for(N in N_grid){
  for(rate in avh_rate_grid){
    n_avh <- max(round(n_trials_condition * rate), 3)
    n_non <- n_trials_condition - n_avh

    D_present <- simulate_D_batch(N * n_sims, "top_down",  n_avh) # AVH epochs: hypothesised top-down
    D_absent  <- simulate_D_batch(N * n_sims, "bottom_up", n_non) # non-AVH epochs: baseline bottom-up

    delta <- (-D_present) - (-D_absent) # flip sign so positive = "more top-down during AVH", the predicted direction
    delta_mat <- matrix(delta, nrow = N, ncol = n_sims)

    m  <- colMeans(delta_mat)
    sdv <- apply(delta_mat, 2, sd)
    tstat <- m / (sdv / sqrt(N))
    pval <- 2 * pt(-abs(tstat), df = N - 1)
    power <- mean(pval < .05 & m > 0)

    main_results[[row_i]] <- data.frame(N = N, avh_rate = rate, n_avh_trials = n_avh,
                 n_nonavh_trials = n_non, power = power, se = sqrt(power * (1 - power) / n_sims))
    row_i <- row_i + 1
  }
}
main_results <- do.call(rbind, main_results)

cat("=== Power to detect AVH-locked top-down shift (conditions 4/6) ===\n")
print(main_results, row.names = FALSE)

# Specificity check: control conditions (1/2/silence)
control_results <- list(); row_i <- 1
for(N in N_grid){
  D_half1 <- simulate_D_batch(N * n_sims, "bottom_up", 30)
  D_half2 <- simulate_D_batch(N * n_sims, "bottom_up", 30)
  delta <- (-D_half1) - (-D_half2)
  delta_mat <- matrix(delta, nrow = N, ncol = n_sims)
  m  <- colMeans(delta_mat)
  sdv <- apply(delta_mat, 2, sd)
  tstat <- m / (sdv / sqrt(N))
  pval <- 2 * pt(-abs(tstat), df = N - 1)
  fpr <- mean(pval < .05)
  control_results[[row_i]] <- data.frame(N = N, false_positive_rate = fpr,
                                          se = sqrt(fpr * (1 - fpr) / n_sims))
  row_i <- row_i + 1
}
control_results <- do.call(rbind, control_results)

cat("\n=== False-positive rate in control conditions (no true shift; nominal alpha = .05) ===\n")
print(control_results, row.names = FALSE)

write.csv(main_results, "roi_connectivity_power.csv", row.names = FALSE)
write.csv(control_results, "roi_connectivity_falsepositive.csv", row.names = FALSE)
write.csv(as.data.frame(G), "roi_leakage_matrix.csv", row.names = TRUE)

# Plot
main_results$avh_pct <- factor(paste0(main_results$avh_rate * 100, "% AVH epochs"),
                                levels = paste0(sort(unique(main_results$avh_rate)) * 100, "% AVH epochs"))

p <- ggplot(main_results, aes(x = N, y = power, colour = avh_pct)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = N_grid) +
  scale_colour_manual(values = c("orange", "purple", "green")) +
  labs(title = "Recoverability of an AVH-locked top-down connectivity shift",
       subtitle = paste0("PAC/STG/IFG directionality after 64-channel \u2192 ", n_components,
                          "-component reduction; 1000 sims/point; dashed = 80% benchmark"),
       x = "Sample size (N participants)", y = "Power", colour = "AVH event rate\n(of 60 trials)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 12.5, face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey"))

ggsave("roi_connectivity_power_plot.png", p, width = 7.5, height = 6, dpi = 300)
