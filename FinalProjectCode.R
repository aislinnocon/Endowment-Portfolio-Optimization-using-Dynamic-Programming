############################################################
## CFRM 507 — FINAL PROJECT
## Aislinn O'Connell
############################################################

# generates multivariate random variables
library(MASS) 

# allows for sets to be reproducible
set.seed(123)

############################################################
#  Part 1 - Set up Problem
############################################################

T_horizon <- 10          # Years
W0 <- 40               # Initial wealth ($ millions)

# Spending and donations per year (millions)
spending  <- c(1.7, 1.8, 1.9, 2.0, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7)
donations <- c(2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 3.5, 3.4, 3.4, 3.4)

# Annual Return Distribution  
mu  <- c(0.057, 0.054, 0.052, 0.050, 0.033, 0.063, 0.028) # mean
sig <- c(0.176, 0.187, 0.243, 0.192, 0.037, 0.066, 0.056) # standard deviation

# Correlation Matrix
corr <- matrix(c(
  1.00,0.74,0.67,0.74,0.13,0.47,0.02,
  0.74,1.00,0.70,0.78,0.09,0.46,0.00,
  0.67,0.70,1.00,0.66,0.07,0.45,-0.03,
  0.74,0.78,0.66,1.00,0.10,0.37,-0.03,
  0.13,0.09,0.07,0.10,1.00,0.10,0.10,
  0.47,0.46,0.45,0.37,0.10,1.00,0.55,
  0.02,0.00,-0.03,-0.03,0.10,0.55,1.00
), byrow = TRUE, nrow = 7)

Sigma <- diag(sig) %*% corr %*% diag(sig)

# Fund Mixes Matrix
# columns correspond assets, & rows correspond to mixes
mix_weights <- rbind(
  c(0.00,0.00,0.00,0.00,0.70,0.08,0.22),
  c(0.00,0.00,0.00,0.00,0.57,0.43,0.00),
  c(0.00,0.00,0.00,0.00,0.34,0.66,0.00),
  c(0.00,0.00,0.00,0.00,0.12,0.88,0.00),
  c(0.03,0.00,0.16,0.00,0.00,0.81,0.00),
  c(0.06,0.00,0.42,0.00,0.00,0.52,0.00),
  c(0.09,0.00,0.69,0.00,0.00,0.22,0.00)
)

n_funds <- nrow(mix_weights)

#  Client Utility Function
u_terminal <- function(W) {
  u <- numeric(length(W)) #placeholder vector

  u[W < 95] <- W[W < 95] - 0.2 * (95 - W[W < 95])^2
  
  idx <- (W >= 95 & W < 100)
  u[idx] <- 1.01 * (W[idx] - exp(-0.001 * (100 - W[idx])))
  
  idx <- (W >= 100 & W < 110)
  u[idx] <- W[idx]
  
  idx <- (W >= 110)
  u[idx] <- 5.8 * (sqrt(W[idx]) - 1) / 0.5
  
  return(u)
}

############################################################
# Part 2 - Set up for Bellman's Backwards Recursion
############################################################

# Set up max & min wealth values
W_min  <- 0
W_max  <- 200
W_step <- 1

# create a grid of wealth values incremented by step
W_grid <- seq(W_min, W_max, by = W_step)
nW <- length(W_grid)

# Number of monte carlo scenarios
n_scenarios <- 5000

# Simulates jointly lognormal returns
log1p_R_assets <- mvrnorm(n_scenarios, mu = mu, Sigma = Sigma)
R_assets <- exp(log1p_R_assets) - 1

# Convert log returns to normal returns
R_funds <- R_assets %*% t(mix_weights)

# Stores most optimal solutions, and fund chosen at year t
V <- matrix(NA_real_, nrow = T_horizon + 1, ncol = nW)
policy <- matrix(NA_integer_, nrow = T_horizon, ncol = nW)

V[T_horizon + 1, ] <- u_terminal(W_grid)

############################################################
# Part 3 - Bellman's for Backwards Recursion
############################################################

# loop backwards
for (t in T_horizon:1) {
  St <- spending[t]
  Dt <- donations[t]
  V_next <- V[t + 1, ]
  
  #loop over each possible wealth value
  for (iw in 1:nW) {
    Wt <- W_grid[iw]
    
    # find how much can be invested at t
    investable <- Wt - St + Dt
    if (investable < 0) investable <- 0
    
    EV <- numeric(n_funds)
    
    #loop over mixed fund choices
    for (f in 1:n_funds) {
      #monte carlo small tree
      W_next_s <- investable * (1 + R_funds[, f])
      
      W_next_s <- pmin(pmax(W_next_s, W_min), W_max)
      
      #interpolate value of next year wealth
      V_interp <- approx(W_grid, V_next, xout = W_next_s, rule = 2)$y
      
      # find expected continuation value
      EV[f] <- mean(V_interp)
    }
    
    #pick the most optimal fund
    policy[t, iw] <- which.max(EV)
    V[t, iw] <- max(EV)
  }
}


############################################################
# Part 4 - Defines forward simulation
############################################################

# Applies optimal decisions and applies them to simulate outcomes
simulate_optimal_paths <- function(n_paths,
                                   W0,
                                   W_grid,
                                   policy,
                                   spending,
                                   donations,
                                   mu,
                                   Sigma,
                                   mix_weights) {
  
  T_horizon <- length(spending)
  n_funds <- nrow(mix_weights)
  
  # Store results of simulation
  W_path <- matrix(NA_real_,   nrow = T_horizon + 1, ncol = n_paths)
  fund_path <- matrix(NA_integer_, nrow = T_horizon, ncol = n_paths)
  
  W_path[1, ] <- W0
  
  # Loop over years t
  for (t in 1:T_horizon) {
    St <- spending[t]
    Dt <- donations[t]
    
    # Loop over simulations
    for (p in 1:n_paths) {
      # Current wealth
      Wt <- W_path[t, p]
      
      # Find optimal fund
      iw <- which.min(abs(W_grid - Wt))
      f  <- policy[t, iw]
      fund_path[t, p] <- f
      
      investable <- Wt - St + Dt
      if (investable < 0) investable <- 0
      
      # Simulate the asset returns
      logR <- mvrnorm(1, mu = mu, Sigma = Sigma)
      R_assets <- exp(logR) - 1 
      R_f <- sum(R_assets * mix_weights[f, ]) # compute the return of fund
      
      # Find next year wealth value
      W_path[t + 1, p] <- investable * (1 + R_f) 
    }
  }
  
  # return simulation results
  list(
    W_path = W_path,
    fund_path = fund_path,
    terminal_W = W_path[T_horizon + 1, ]
  )
}

############################################################
# Part 5 - Solves forward simulation
############################################################

# Solves Part 4
res <- simulate_optimal_paths(
  n_paths = 5000,
  W0 = W0,
  W_grid = W_grid,
  policy = policy,
  spending = spending,
  donations = donations,
  mu = mu,
  Sigma = Sigma,
  mix_weights = mix_weights
)

############################################################
# Part 6 - Collect Data
############################################################
# Finds mean of wealth
mean_terminal <- mean(res$terminal_W)
prob_100 <- mean(res$terminal_W >= 100) # probability it reaches W = 100
prob_110 <- mean(res$terminal_W >= 110)

cat("Mean terminal wealth:", mean_terminal, "\n")
cat("Prob(W_T >= 100):", prob_100, "\n")
cat("Prob(W_T >= 110):", prob_110, "\n")

# Wealth matrix
W_path    <- res$W_path
fund_path <- res$fund_path
terminal_W <- res$terminal_W

# Example: show wealth for first few paths
head(t(W_path))
head(t(fund_path))


# Finds wealth at each value t and the fund
avg_wealth <- rowMeans(W_path)

# Most common fund selected at each year (mode)
mode_fund <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
avg_fund_choice <- apply(fund_path, 1, mode_fund)

cat("\n==== Average Values Across All Paths ====\n")
for (t in 1:T_horizon) {
  cat(sprintf("Year %d: Avg Wealth = %.2f, Most-Chosen Fund = %d\n",
              t, avg_wealth[t], avg_fund_choice[t]))
}
cat("\nTerminal Avg Wealth:", avg_wealth[T_horizon + 1], "\n")

# chooses last wealth value
median_index <- which.min(abs(terminal_W - median(terminal_W)))

#creates table to print out wealth and fund
rep_wealth <- W_path[, median_index]        # length 11 (0..10)
rep_fund   <- fund_path[, median_index]     # length 10 (years 1..10)
rep_table <- data.frame(
  Year = 0:T_horizon,
  Wealth = rep_wealth,
  Fund_Chosen = c(NA, rep_fund)  # NA at year 0
)

cat("\n==== Representative Path (Median Terminal Wealth) ====\n")
print(rep_table, row.names = FALSE)
cat("\nTerminal Wealth (Representative):", rep_wealth[T_horizon + 1], "\n")

cat("\n==== Summary Statistics ====\n")
cat("Mean terminal wealth:", mean_terminal, "\n")
cat("Median terminal wealth:", median(terminal_W), "\n")
cat("Std dev terminal wealth:", sd(terminal_W), "\n")
cat("Prob terminal >= 100:", prob_100, "\n")
cat("Prob terminal >= 110:", prob_110, "\n")
