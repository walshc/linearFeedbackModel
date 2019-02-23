lfm <- function(formula, data, effect = "individual", model = "onestep",
                weight.matrix = "identity", index = NULL) {

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
  if (is.null(index)) {
    if ((inherits(data, "pdata.frame")) == 0)
      stop("'data' is not a pdata.frame.
           Use 'index' option or convert to pdata.frame")
  } else {
    data <- pdata.frame(data, index = index)
  }

  # Check if data is balanced and if not balance it:
  if (!plm::pdim(data)$balanced) {
    un.id <- sort(unique(index(data, "id")))
    un.time <- sort(unique(index(data, "time")))
    rownames(data) <- paste(index(data, "id"), index(data, "time"), sep = ".")
    allRows <- as.character(t(outer(un.id, un.time, paste, sep = ".")))
    data <- data[allRows, ]
    rownames(data) <- allRows
    index <- data.frame(id = rep(un.id, each = length(un.time)),
                        time = rep(un.time, length(un.id)),
                        row.names = rownames(data))
    class(index) <- c("pindex", "data.frame")
    attr(data, "index") <- index
  }

  # Check if effect option was used correctly:
  if (effect != "individual" & effect != "twoways")
    stop("'effect' must be either 'individual' or 'twoways'.
         Default is 'individual'.")

  # Check if model option was used correctly:
  if (model != "onestep" & model != "twosteps")
    stop("'model' must be either 'onestep' or 'twosteps'.
         Default is 'onestep'.")

  if (weight.matrix != "identity" & weight.matrix != "instruments")
    stop("'weight.matrix' must be either 'identity' or 'instruments'.
         Default is 'identity'.")

  x <- formula
  if (!inherits(x, "Formula")) x <- Formula(formula)
  # gmm instruments : named list with the lags, names being the variables
  gmm.form <- formula(x, rhs = 2, lhs = 0)

  # Function from the plm package:
  getvar <- function(x) {
    x <- as.list(x)
    result <- lapply(x, function(y) {
      deb <- as.numeric(gregexpr("lag\\(", y)[[1]])
      if (deb == -1) {
          lags <- 0
          avar <- y
      } else {
        inspar <- substr(y, deb + 4, nchar(y) - 1)
        coma <- as.numeric(gregexpr(",", inspar)[[1]][1])
        if (coma == -1) {
            endvar <- nchar(inspar)
            lags <- 1
        } else {
          endvar <- coma - 1
          lags <- substr(inspar, coma + 1, nchar(inspar))
          lags <- eval(parse(text = lags))
        }
        avar <- substr(inspar, 1, endvar)
      }
      list(avar, lags)
    })
    nres <- sapply(result, function(x) x[[1]])
    result <- lapply(result, function(x) x[[2]])
    names(result) <- nres
    result
  }
  # Function from the plm package:
  dynterms <- function(x) {
    trms.lab <- attr(terms(x), "term.labels")
    result <- getvar(trms.lab)
    nv <- names(result)
    dn <- names(table(nv))[table(nv) > 1]
    un <- names(table(nv))[table(nv) == 1]
    resu <- result[un]
    for (i in dn) {
      v <- sort(unique(unlist(result[nv == i])))
      names(v) <- NULL
      resu[[i]] <- v
    }
    resu
  }
  gmm.lags <- dynterms(gmm.form)

  main.form <- formula(x, rhs = 1, lhs = 1)
  main.lags <- dynterms(main.form)

  # How many time series are lost ? May be the maximum number of lags
  # of any covariates + 1 because of first - differencing or the
  # largest minimum lag for any gmm or normal instruments
  gmm.minlag <- max(sapply(gmm.lags, min))
  inst.maxlag <- 0
  main.maxlag <- max(sapply(main.lags, max))
  TL1 <- max(main.maxlag + 1, gmm.minlag)
  TL2 <- max(main.maxlag, gmm.minlag - 1)

  gmm.form <- as.formula(paste("~", paste(names(gmm.lags), collapse = "+")))
  Form <- as.Formula(main.form, gmm.form)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("plm")
  mf$model <- NA
  mf$formula <- Form
  mf$na.action <- "na.pass"
  mf$subset <- NULL
  data <- eval(mf, parent.frame())

  attr(data, "formula") <- formula(main.form)

  extractData <- function(x, as.matrix = TRUE) {
    form <- attr(data, "formula")
    trms <- terms(form)
    has.response <- attr(trms, "response") == 1
    has.intercept <- attr(trms, "intercept") == 1
    if (has.intercept == 1) {
      form <- Formula(update(formula(form), ~. - 1))
    }
    index <- attr(data, "index")
    X <- model.matrix(form, data)
    if (has.response) {
      X <- cbind(data[[1]], X)
      colnames(X)[1] <- deparse(trms[[2]])
    }
    data <- split(as.data.frame(X), index[[1]])
    time <- split(index[[2]], index[[1]])
    data <- mapply(function(x, y) {
      rownames(x) <- y
      if (as.matrix) {
        x <- as.matrix(x)
      }
      x
    }, data, time, SIMPLIFY = FALSE)
    data
  }
  yX <- do.call(rbind, extractData(data))

  attr(data, "formula") <- gmm.form
  W <- extractData(data, as.matrix = FALSE)
  makegmm <- function(x, g, TL1) {
    nT <- length(x)
    rg <- range(g)
    z <- as.list((TL1 + 1):nT)
    x <- lapply(z, function(y) x[max(1, y - rg[2]):(y - rg[1])])
    lx <- sapply(x, length)
    n <- length(x)
    lxc <- cumsum(lx)
    before <- c(0, lxc[-n])
    after <- lxc[n] - sapply(x, length) - before
    result <- t(mapply(function(x, y, z) c(rep(0, y), x,
      rep(0, z)), x, before, after, SIMPLIFY = TRUE))
    result
  }
  Z <- lapply(W, function(x) {
    u <- mapply(makegmm, x, gmm.lags, TL1, SIMPLIFY = FALSE)
    matrix(unlist(u), nrow = nrow(u[[1]]))
  })

  N  <- length(unique(attr(data, "index")[[1]]))  # Number of individual obs.
  nT <- length(unique(attr(data, "index")[[2]]))  # Number of time periods.

  # Model data frame:
  row.names(yX) <- NULL
  mdf <- pdata.frame(cbind(index(data), yX), index = names(index(data)))
  names(mdf)[1:2] <- c("i", "t")

  # Convert indices from factor to integer:
  mdf$i <- as.integer(mdf$i)
  mdf$t <- as.integer(mdf$t)

  # Drop missings (from taking lags):
  mdf <- na.omit(mdf)

  # Add time fixed effects if effect = "twoways":
  if (effect == "twoways") {
    for (t in 3:nT) {
      mdf[[paste0("t.", t)]] <- ifelse(mdf$t == t, 1, 0)
    }
  }

  # Number of regressors (exclude id, time, dependent variable and its lag):
  K <- ncol(mdf) - 4

  if (weight.matrix == "identity") {
    # Identity matrix:
    W1 <- diag(1, nrow = ncol(Z[[1]]), ncol = ncol(Z[[1]]))
  } else if (weight.matrix == "instruments") {
    # Use sum(Z'_i*Z_i):
    W1 <- firstWeightMatrix(Z = do.call(rbind, Z))
  }

  # Invert the weighting matrix (taking the general inverse if necessary):
  W1.inv <- MASS::ginv(W1)
  if (min(eigen(W1)$values) < .Machine$double.eps) {
    warning("The first-step matrix is singular, a general inverse is used")
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
  if (K == 0) {
    first <- stats::optimize(f = GMMfirstStep, interval = c(-1e3, 1e3))
    first$par <- first$minimum
  } else {
    first <- stats::optim(rep(0, K + 1), GMMfirstStep)
  }
  names(first$par) <- names(mdf)[4:ncol(mdf)]

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
    W2.inv <- MASS::ginv(W2)
    if (min(eigen(W2)$values) < .Machine$double.eps) {
      warning("The second-step matrix is singular, a general inverse is used")
    }

    GMMsecondStep <- function(theta) {
      GMM(theta = as.double(theta),
          idx   = as.matrix(mdf[, 1:2]),
          nT    = as.integer(max(mdf$t)),
          data  = as.matrix(mdf[, 3:ncol(mdf)]),
          Z     = as.matrix(do.call(rbind, Z)),
          W_inv = as.matrix(W2.inv))
    }

    if (K == 0) {
      second <- stats::optimize(f = GMMsecondStep, interval = c(-1e3, 1e3))
      second$par <- second$minimum
    } else {
      second <- stats::optim(rep(0, K + 1), GMMsecondStep)
    }

    names(second$par) <- names(mdf)[4:ncol(mdf)]
    result <- list(call = cl, coefficients = second$par, model = mdf,
                   first = first$par, W1 = W1, W2 = W2, Z = Z)
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
  mu <- aggregate((mdf[[3]] - result$coefficients[[1]] * mdf[[4]])/
                   exp(xBeta) ~ mdf$i, FUN = mean, na.rm = TRUE)
  result$fixed.effects <- data.frame(id = 1:N, fixed.effects = mu[[2]])

  # Get the fitted values:
  result$fitted.values <-
    data.frame(id = mdf$i, time = mdf$t,
               fitted.values = result$coefficients[[1]] * mdf[[4]] +
                               exp(xBeta) * rep(mu[[2]], each = nT - 1))

  # Get the residuals:
  result$residuals <-
    data.frame(id = mdf$i, time = mdf$t,
               residuals =  mdf[[3]] - result$fitted.values$fitted.values)

  # Find the variance-covariance matrix:
  mdf$mu <- quasiDifference(result$coefficients)$mu
  q.prime <- matrix(0, N * (nT - 2), K + 1)

  # First column is derivative wrt the coefficient on the lagged dependent var:
  q.prime[, 1] <- stats::na.omit(- mdf[[4]] / mdf$mu -
                                 - lag(mdf[[4]], k = 1) /
                                   lag(mdf$mu, k = 1))

  # Remaining columns are the derivatives w.r.t. the coefficients on covariates:
  if (K > 0) {
    for (k in 1:K) {
      q.prime[, k + 1] <-
        stats::na.omit(- mdf[[k + 4]] * (mdf[[3]] -
                       result$coefficients[[1]] * mdf[[4]]) / mdf$mu +
                       lag(mdf[[k + 4]], k = 1) * (lag(mdf[[3]], k = 1) -
                       result$coefficients[[1]] * lag(mdf[[4]], k = 1)) /
                                                  lag(mdf$mu, k = 1))
    }
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
  result$vcov <- (1 / N) * MASS::ginv(D %*% MASS::ginv(W) %*% t(D))
  class(result) <- "lfm"
  return(result)
}
