# Load/install necessary packages:
if (!suppressPackageStartupMessages(require(pacman))) install.packages("pacman")
pacman::p_load(data.table, ggplot2, plm)
pacman::p_load_gh("walshc/linearFeedbackModel")

library(linearFeedbackModel)
detach("package:linearFeedbackModel", unload = TRUE)
install.packages(".", repos = NULL, type = "source")
library(linearFeedbackModel)

set.seed(2435262)

gamma       <- 0.5
beta        <- 0.5
rho         <- 0.5
tau         <- 0.1
var_eta     <- 0.5
var_epsilon <- 0.5

N <- 10000
start <- -50
nT <- 10
df <- data.table(expand.grid(i = 1:N, t = start:nT))

df[, eta := rnorm(length(unique(i)), mean = 0, sd = sqrt(var_eta)), by = i]
df[, xi := rnorm(length(unique(i)), mean = 0,
                 sd = sqrt(var_epsilon / (1 - rho^2))), by = i]
df[t == start, x := (tau / (1 - rho)) * eta + xi]
df[t == start, y := rpois(.N, lambda = exp(x * beta + eta))]

for (j in (start+1):nT) {
  df[t == j, lag.x := df[t == j - 1, x]]
  df[t == j, lag.y := df[t == j - 1, y]]
  df[t == j, x := rho * lag.x + tau * eta + rnorm(.N)]
  df[t == j, y := rpois(.N, lambda = gamma * lag.y + exp(x * beta + eta))]
}
df <- df[t > 0]

print(lfm(y ~ lag(y, k = 1) + x | lag(y, k = 2:4) + lag(x, k = 1:4),
          data = df, effect = "individual", model = "twosteps",
          index = c("i", "t")))
