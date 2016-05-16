lfm <- function(formula, data, effect = "individual", model = "onestep") {
  if(!require(plm)) stop("Package 'plm' required.")
  if(!require(MASS)) stop("Package 'MASS' required.")
  if(!require(Rcpp)) stop("Package 'Rcpp' required.")
  require(plm)
  require(MASS)
  require(Rcpp)

  # Store the function call:
  cl <- match.call()

  # If the dependent variable is dep.var, the first RHS variable must be
  # lag(dep.var, k = 1) so check for lag and dep.var
  dep.var <- as.character(formula[[2]])
  first.indep.var <-   # First independent variable:
    strsplit(as.character(formula[[3]][2]), "+", fixed = TRUE)[[1]][1]

  if (!(grepl("lag", first.indep.var) & grepl(dep.var, first.indep.var))) {
    stop("First independent variable is not the first lag of the dependent
          variable. Use lag(depvar, k = 1) as the first variable after ~")
  }

  # Check if data is a pdata.frame:
  if ((inherits(data, "pdata.frame")) == 0)
    stop("'data' is not a pdata.frame. 
         Use pdata.frame() from package plm to convert.")

  # Check if data is balanced:
  if (pdim(data)[[3]] != 1)
    stop("'data' is not a balanced panel.
         Unbalanced panel data are not currently supported")

  # Check if effect option was used correctly:
  if (effect != "individual" & effect != "twoways")
    stop("'effect' must be either 'individual' or 'twoways'.
         Default is 'individual'.")

  # Check if model option was used correctly:
  if (model != "onestep" & model != "twosteps")
    stop("'model' must be either 'onestep' or 'twosteps'.
         Default is 'onestep'.")


  # Create the instrument matrix (W). The pgmm function does this (same
  # instruments as Arellano Bond). Use those coefficients as starting values
  # too.
  arellano.bond <- pgmm(formula, effect = effect, model = model, data = data)
  Z <- arellano.bond$W
  if (model == "onestep") {
    start <- arellano.bond$coefficients  # Use the AB as starting values
  } else {
    start <- arellano.bond$coefficients[[1]]
  }
  rm(arellano.bond)

  # Create a data set from the given formula (excluding the instruments):
  mdf <- model.frame(pFormula(formula), data = data, rhs = 1)

  N  <- length(unique(attr(data, "index")[[1]]))  # Number of individual obs.
  nT <- length(unique(attr(data, "index")[[2]]))  # Number of time periods.

  # Add ID and time identifiers to the data frame (in the 1st and 2nd column):
  mdf <- data.frame(id   = rep(1:N, each = nrow(mdf)/N),
                    time = rep(2:nT, N), mdf)
  colnames(mdf)[3:4] <- c("dep.var", "lag.dep.var")

  # Add time fixed effects if effect = "twoways":
  if (effect == "twoways") {
    for (t in 3:nT) {
      mdf[[paste0("time.", t)]] <- ifelse(mdf$time == t, 1, 0)
    }
  }

  # Number of regressors (exclude id, time, dependent variable and its lag):
  K <- ncol(mdf) - 4

  # Get the weighting matrix (use sum(Z'_i*Z_i) in the first step):
  W1 <- firstWeightMatrix(Z = do.call(rbind, Z))
  # Use identity matrix instead:
  # W1 <- diag(1, nrow = ncol(Z[[1]]), ncol = ncol(Z[[1]]))

  # Invert the weighting matrix (taking the general inverse if necessary):
  minevW1 <- min(eigen(W1)$values)
  if (minevW1 < .Machine$double.eps) {
    W1.inv <- ginv(W1)
    warning("The first-step matrix is singular, a general inverse is used")
  } else {
    W1.inv <- solve(W1)
  }

  GMMfirstStep <- function(theta) {
    GMM(theta = as.double(theta),
        idx   = as.matrix(mdf[, 1:2]),
        nT    = as.integer(max(mdf$t)),
        data  = as.matrix(mdf[, 3:ncol(mdf)]),
        Z     = as.matrix(do.call(rbind, Z)),
        W_inv = as.matrix(W1.inv))
  }

  # Obtain a first step estimate of theta from the initial weight matrix:
  first <- optim(start, GMMfirstStep)

  quasiDifference <- function(theta) {
    out <- qMu(theta = as.double(theta),
               idx   = as.matrix(mdf[, 1:2]),
               nT    = as.integer(max(mdf$t)),
               data  = as.matrix(mdf[, 3:ncol(mdf)]))
    names(out) <- c("q", "mu")
    return(out)
  }

  if (model == "twosteps") {
    # Find the efficient weight matrix (sum(Z'_i*q_i*q'_i,Z_i)) using
    # the estimates from the first step:
    W2 <- secondWeightMatrix(theta = as.double(first$par),
                             idx   = as.matrix(mdf[, 1:2]),
                             nT    = as.integer(max(mdf$t)),
                             data  = as.matrix(mdf[, 3:ncol(mdf)]),
                             Z     = do.call(rbind, Z))
    # Invert the second-step weight matrix:
    minevW2 <- min(eigen(W2)$values)
    if (minevW2 < .Machine$double.eps) {
      W2.inv <- ginv(W2)
      warning("The second-step matrix is singular, a general inverse is used")
    } else {
      W2.inv <- solve(W2)
    }

    GMMsecondStep <- function(theta) {
      GMM(theta = as.double(first$par),
          idx   = as.matrix(mdf[, 1:2]),
          nT    = as.integer(max(mdf$t)),
          data  = as.matrix(mdf[, 3:ncol(mdf)]),
          Z     = as.matrix(do.call(rbind, Z)),
          W_inv = as.matrix(W2.inv))
    }

    second <- optim(start, GMMsecondStep)
    result <- list(call = cl, coefficients = second$par, model = mdf,
                   first = first$par, W1 = W1, W2 = W2, Z = Z, call = cl)
  } else {
    result <- list(call = cl, coefficients = first$par, model = mdf,
                   W1 = W1, Z = Z)
  }

  # Get an estimate of the fixed effects:
  # Find X*Beta if there are any X's:
  xBeta <- rep(0, nrow(mdf))
  if (K > 0) {
    for (k in 1:K) {
      xBeta <- xBeta + result$coefficients[[k + 1]] * mdf[[k + 4]]
    }
  }
  mu <- aggregate((mdf$dep.var - result$coefficients[[1]] * mdf$lag.dep.var)/
                   exp(xBeta) ~ mdf$id, FUN = mean, na.rm = TRUE)
  result$fixed.effects <- data.frame(id = 1:N, fixed.effects = mu[[2]])

  # Get the fitted values:
  result$fitted.values <-
    data.frame(id = mdf$id, time = mdf$time,
               fitted.values = result$coefficients[[1]] * mdf$lag.dep.var +
                               exp(xBeta)*rep(mu[[2]], each = nT - 1))

  # Get the residuals:
  result$residuals <-
    data.frame(id = mdf$id, time = mdf$time,
               residuals =  mdf$dep.var - result$fitted.values$fitted.values)

  # Find the variance-covariance matrix:
  mdf$mu <- quasiDifference(result$coefficients)$mu
  q.prime <- matrix(0, N * (nT - 2), K + 1)

  # First column is derivative wrt the coefficient on the lagged dependent var:
  q.prime[, 1] <- na.omit(- mdf$lag.dep.var/mdf$mu -
                          - lag(mdf$lag.dep.var, k = 1)/lag(mdf$mu, k = 1))

  # Remaining columns are the derivatives w.r.t. the coefficients on covariates:
  for (k in 1:K) {
    q.prime[, k + 1] <-
      na.omit(- mdf[[k + 4]] * (mdf$dep.var -
              result$coefficients[[1]] * mdf$lag.dep.var) / mdf$mu +
              lag(mdf[[k + 4]], k = 1) * (lag(mdf$dep.var, k = 1) -
              result$coefficients[[1]] * lag(mdf$lag.dep.var, k = 1)) /
                                         lag(mdf$mu, k = 1))
  }

  # Turn the matrix into a list of matrices, split by individuals:
  id <- rep(1:N, each = nT - 2)
  q.prime.list <-
    lapply(split(q.prime, id),
           function(x) matrix(x, nT - 2, length(result$coefficients)))

  # Multiply by the instrument matrix:
  tmp <- as.list(data.frame(do.call(rbind, list(lapply(q.prime.list, t), Z))))
  q.prime.Z <- lapply(tmp, function(x) x[[1]] %*% x[[2]])

  # Average over individuals:
  D <- Reduce("+", q.prime.Z) / length(q.prime.Z)

  if (model == "onestep") {
    W <- W1
  } else {
    W <- W2
  }

  result$D <- D
  result$vcov <- (1/N) * solve(D %*% solve(W) %*% t(D))
  class(result) <- "lfm"
  return(result)
}

# This hides the large matrices from printing when the function is called:
print.lfm <- function(x) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}

# This shows the table when summary(lfm) is called (not finished yet):
summary.lfm <- function(x) {
  require(lmtest)
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(coeftest(x, vcov = x$vcov))
}
