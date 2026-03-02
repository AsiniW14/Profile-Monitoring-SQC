#libraries
library(ggplot2)

#global Parameters
set.seed(123) 
n <- 10        # number of observations per profile
x_values <- 1:n # fixed observation values for all profiles

# true in-control parameters
true_beta0_in_control <- 5
true_beta1_in_control <- 2
true_sigma_error <- 0.5 # standard deviation of the error term

#phase I data generation and establishment 
m1_phase1 <- 30 

# creating  a data frame to store phase I data
phase1_data <- data.frame(
  profile_id = rep(1:m1_phase1, each = n),
  x = rep(x_values, m1_phase1)
)
# generating the y values for phase I profiles
phase1_data$y <- true_beta0_in_control + true_beta1_in_control * phase1_data$x +
  rnorm(m1_phase1 * n, mean = 0, sd = true_sigma_error)

# phase I parameter estimation and control Limit calculation
# fit linear model to each phase I profile
phase1_models <- lapply(split(phase1_data, phase1_data$profile_id),
                        function(df) lm(y ~ x, data = df))

# extract parameter estimates (beta0_hat, beta1_hat) for each phase I profile
beta_hats_phase1 <- t(sapply(phase1_models, coef)) # transpose to get m1 x p matrix

# calculatingg the overall in-control parameter vector (mean of Phase I estimates)
# this is our estimate for beta_0 in the T^2 formula
beta_0_hat_overall <- colMeans(beta_hats_phase1)

#calculating the estimated in-control covariance matrix of the parameter estimates
# this is the estimate for Sigma_beta_0 in the T^2 formula
Sigma_beta_0_hat <- cov(beta_hats_phase1)

#calculate Hotelling's T^2 for for retrospective analysis
T2_phase1_values <- numeric(m1_phase1)
for (i in 1:m1_phase1) {
  # (beta_hat_i - beta_0_hat_overall)
  beta_diff <- matrix(beta_hats_phase1[i, ] - beta_0_hat_overall, ncol = 1)
  
  #T^2 = (beta_diff)^T %*% solve(Sigma_beta_0_hat) %*% beta_diff
  T2_phase1_values[i] <- t(beta_diff) %*% solve(Sigma_beta_0_hat) %*% beta_diff
}

# determine Hotelling's T^2 control limit for Phase II 
# p is the number of parameters being monitored
p <- 2 
alpha <- 0.0027 # standard alpha

UCL_T2 <- qchisq(1 - alpha, df = p)

# phase I plot Hotelling's T^2 for the retrospective analysis
plot_df_phase1 <- data.frame(Profile = 1:m1_phase1, T2 = T2_phase1_values)
ggplot(plot_df_phase1, aes(x = Profile, y = T2)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  geom_hline(yintercept = UCL_T2, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = max(plot_df_phase1$Profile) * 0.9, y = UCL_T2 * 1.05, 
           label = paste("UCL =", round(UCL_T2, 2)), color = "red", size = 4) +
  labs(title = "Phase I: Hotelling's T² Chart for Linear Profiles",
       y = "Hotelling's T² Statistic", x = "Profile Number") +
  theme_minimal() +
  ylim(0, max(T2_phase1_values, UCL_T2) * 1.1)

# phase II data generation with a shift 
m2_phase2 <- 20 # number of profiles in phase II
shift_start_profile <- 11 # shift starts at this profile number
shifted_beta1 <- 2.5    # new shifted value

# create a data frame for phase II
phase2_data <- data.frame(
  profile_id = rep(1:m2_phase2, each = n),
  x = rep(x_values, m2_phase2)
)

# generate y values for phase II profiles, introducing a shift
for (i in 1:m2_phase2) {
  current_beta0 <- true_beta0_in_control
  current_beta1 <- ifelse(i >= shift_start_profile, shifted_beta1, true_beta1_in_control)
  
  start_row <- (i - 1) * n + 1
  end_row <- i * n
  
  phase2_data$y[start_row:end_row] <- current_beta0 + current_beta1 * 
    phase2_data$x[start_row:end_row] + rnorm(n, mean = 0, sd = true_sigma_error)
}

# phase II: monitoring and plotting
# fit linear model to each Phase II profile
phase2_models <- lapply(split(phase2_data, phase2_data$profile_id),
                        function(df) lm(y ~ x, data = df))

# extract parameter estimates for each ehase II profile
beta_hats_phase2 <- t(sapply(phase2_models, coef))


# calculate Hotelling's T^2 for phase II profiles
T2_phase2_values <- numeric(m2_phase2)
for (i in 1:m2_phase2) {
  beta_diff <- matrix(beta_hats_phase2[i, ] - beta_0_hat_overall, ncol = 1)
  T2_phase2_values[i] <- t(beta_diff) %*% solve(Sigma_beta_0_hat) %*% beta_diff}

# phase II: plot Hotelling's T^2 
plot_df_phase2 <- data.frame(Profile = 1:m2_phase2, T2 = T2_phase2_values)
ggplot(plot_df_phase2, aes(x = Profile, y = T2)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  geom_hline(yintercept = UCL_T2, linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(xintercept = shift_start_profile - 0.5, linetype = "dotted", color = "darkgreen", linewidth = 1) +
  annotate("text", x = max(plot_df_phase2$Profile) * 0.9, y = UCL_T2 * 1.05, 
           label = paste("UCL =", round(UCL_T2, 2)), color = "red", size = 4) +
  annotate("text", x = shift_start_profile - 1, y = max(T2_phase2_values, UCL_T2) * 0.95, 
           label = "Shift Starts", color = "darkgreen", hjust = 1, size = 3) +
  labs(title = "Phase II: Hotelling's T² Chart for Linear Profiles (Shift in Slope)",
       y = "Hotelling's T² Statistic", x = "Profile Number") +
  theme_minimal() +
  ylim(0, max(T2_phase2_values, UCL_T2) * 1.1)

cat("Phase II Monitoring Results")
out_of_control_profiles <- which(T2_phase2_values > UCL_T2)
if (length(out_of_control_profiles) > 0) {
  cat("Out-of-control signals detected at profiles:", out_of_control_profiles)
} else { cat("No out-of-control signals detected in Phase II.")}


#summary of key estimates
cat("Summary of Phase I Estimates")
cat("Estimated In-Control Profile Parameters (beta_0, beta_1):")
print(beta_0_hat_overall)
cat("Estimated In-Control Covariance Matrix of Parameters:")
print(Sigma_beta_0_hat)
cat(paste("Calculated Hotelling's T^2 UCL for Phase II:", round(UCL_T2, 2)))

cat("Phase II Monitoring Results")
out_of_control_profiles <- which(T2_phase2_values > UCL_T2)
if (length(out_of_control_profiles) > 0) {
  cat("Out-of-control signals detected at profiles:", out_of_control_profiles)
} else {
  cat("No out-of-control signals detected in Phase II.")
}

