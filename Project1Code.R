# ldi_run.R
# Pension LDI LP — uses RealCashflow.txt now that it's uploaded.
# Outputs: bond_purchases.csv, portfolio_summary.csv, allocations_by_type.csv, credit_cap_sensitivity.csv

library(lpSolve)
# Set directory - will comment out for turning in
setwd("C:/Users/aisli/OneDrive/Documents/CFRM/Fall_2025/CFRM 507/Project1")

#######################################################
# CFRM 507 - Project 1
# by Aislinn O'Connell
#######################################################

# import necessary Libraries
suppressPackageStartupMessages({
  library(ROI)
  library(ROI.plugin.glpk)
})

# --- Load Data ---
benefit_lines <- scan("BenefitPayments.txt", what="", sep="\n", quiet=TRUE)
Q <- as.integer(trimws(benefit_lines[1]))
QP <- as.numeric(sapply(benefit_lines[2:(1+Q)], trimws))
QP <- QP / 1000  # convert benefits from dollars to 'per $1000 par' scale

bond_lines <- readLines("BondData.txt")
nb <- as.integer(trimws(strsplit(bond_lines[1], "\\s+")[[1]][1]))
bond_rows <- strsplit(bond_lines[2:(1+nb)], "\\s+")
bond_mat <- t(sapply(bond_rows, as.numeric))
colnames(bond_mat) <- c("type","rating","price","ub")
BOND <- as.data.frame(bond_mat)
types   <- BOND$type
ratings <- BOND$rating
prices  <- BOND$price
ubs     <- BOND$ub

nom_cf  <- as.matrix(read.table("NomCashflow.txt",  header=FALSE))
real_cf <- as.matrix(read.table("RealCashflow.txt", header=FALSE))
stopifnot(nrow(nom_cf)==Q, ncol(nom_cf)==nb)

# --- Parameters ---
r_cash_ann <- 0.02
r_q <- (1 + r_cash_ann)^(1/4) - 1
infl_ann <- 0.03
d_q <- (1 + infl_ann)^(- (0:(Q-1))/4)
hold_caps <- c(`1`=0.30, `2`=0.80, `3`=0.80, `4`=0.10)
credit_cap <- 0.4
nom_frac_cap <- 0.2

# --- Decision Variables ---
nv <- nb + (Q + 1)
idx_X <- function(b) b
idx_S <- function(t) nb + 1 + t

# Objective: Minimize Z = V + s(0) + s(80)
cvec <- numeric(nv)
cvec[1:nb] <- prices
cvec[idx_S(0)] <- 1
cvec[idx_S(Q)] <- 1

# --- Equality constraints: Cash balance each quarter ---
Aeq <- matrix(0, nrow=Q, ncol=nv)
beq <- numeric(Q)
for (t in 1:Q) {
  Aeq[t, idx_S(t)]   <- 1
  Aeq[t, idx_S(t-1)] <- -(1 + r_q)
  inflow <- real_cf[t, ] + d_q[t] * nom_cf[t, ]
  Aeq[t, 1:nb] <- -inflow
  beq[t] <- -QP[t]
}

# --- Inequality constraints (A_ub x <= b_ub) ---
Aub <- list(); bub <- c()

# Holding limits by type
for (g in names(hold_caps)) {
  cap <- hold_caps[[g]]
  row <- numeric(nv)
  row[1:nb] <- ((types == as.numeric(g)) - cap) * prices
  Aub[[length(Aub)+1]] <- row; bub <- c(bub, 0)
}

# Credit quality limit
row <- numeric(nv)
row[1:nb] <- (ratings - credit_cap) * prices
Aub[[length(Aub)+1]] <- row; bub <- c(bub, 0)

# Nominal cash flow fraction limit
nom_coeff <- colSums(nom_cf)
tot_coeff <- colSums(nom_cf + real_cf)
row <- numeric(nv)
row[1:nb] <- nom_coeff - nom_frac_cap * tot_coeff
Aub[[length(Aub)+1]] <- row; bub <- c(bub, 0)

# Upper bounds per bond
for (b in 1:nb) {
  row <- numeric(nv)
  row[b] <- prices[b]
  Aub[[length(Aub)+1]] <- row; bub <- c(bub, ubs[b])
}

Aub <- do.call(rbind, Aub)

# --- Variable bounds ---
lower <- rep(0, nv)
upper <- rep(Inf, nv)

# --- Solve LP using GLPK ---
Amat <- rbind(Aub, Aeq)
dir  <- c(rep("<=", nrow(Aub)), rep("==", nrow(Aeq)))
rhs  <- c(bub, beq)

op <- OP(
  objective = L_objective(cvec),
  constraints = L_constraint(L = Amat, dir = dir, rhs = rhs),
  bounds = V_bound(li = 1:nv, ui = 1:nv, lb = lower, ub = upper),
  maximum = FALSE
)

cat("Solving LP with GLPK... this may take a moment...\n")
sol <- ROI_solve(op, solver = "glpk")

if (sol$status$code != 0) stop(paste("GLPK solver failed:", sol$status$msg))

xopt <- solution(sol)

# --- Extract results ---
X <- xopt[1:nb]
S <- xopt[(nb+1):(nb+1+Q)]
V <- sum(prices * X)
Z <- V + S[1] + S[Q+1]

# --- Summaries ---
alloc_by_type <- tapply(prices*X, types, sum)
alloc_by_type <- alloc_by_type[as.character(1:4)]; alloc_by_type[is.na(alloc_by_type)] <- 0
credit_wavg <- if (V > 0) sum(ratings * prices * X)/V else 0
nom_total   <- sum(nom_coeff * X)
tot_total   <- sum(tot_coeff * X)
nom_share   <- if (tot_total > 0) 100 * nom_total / tot_total else 0

summary_df <- data.frame(
  Objective_Z = Z,
  V_PurchaseCost = V,
  S0_InitialCash = S[1],
  S80_LeftoverCash = S[Q+1],
  CreditWeightedAvg = credit_wavg,
  NominalSharePercent = nom_share
)

print(summary_df)
sum(prices * X)   # should be realistic, likely 5–50 million
sum(X > 0)        # how many bonds used


# --- Outputs ---
write.csv(summary_df, "solution_summary.csv", row.names = FALSE)
write.csv(data.frame(
  Type = c("CIPS","TIPS","Nominal","SMA"),
  Dollars = as.numeric(alloc_by_type),
  Cap = c(0.30,0.80,0.80,0.10)
), "allocations_by_type_solution.csv", row.names = FALSE)

write.csv(data.frame(
  Bond_Index = 1:nb,
  Type = types,
  Rating = ratings,
  Price = prices,
  Upper_Bound = ubs,
  Xb_Units = X,
  Purchase_Dollars = prices*X
), "bond_purchases_solution.csv", row.names = FALSE)

cat("\n=== SOLVED SUCCESSFULLY ===\n")
print(summary_df)
cat("\nResults saved:\n- solution_summary.csv\n- allocations_by_type_solution.csv\n- bond_purchases_solution.csv\n")