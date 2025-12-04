ocpt.meanvar.initialise <- function(data,
                                    penalty         = "Manual",
                                    pen.value       = length(data),
                                    Q               = 5,
                                    test.stat       = "Normal",
                                    class           = TRUE,
                                    param.estimates = TRUE,
                                    shape           = 1,
                                    minseglen       = 2,
                                    alpha           = 1,
                                    verbose         = FALSE) {
  checkData(data)
  
  ## ECP branch: unchanged
  if (test.stat == "ECP") {
    ecpans <- e.cp3o_delta.online.initialise(
      Z       = data,
      K       = Q,
      delta   = minseglen + 1,
      alpha   = alpha,
      verbose = verbose
    )
    return(ecpans)
  }
  
  ## Enforce minimum segment length
  if (minseglen < 2) {
    minseglen <- 2
    warning("Minimum segment length for a change in mean and variance is 2, ",
            "automatically changed to be 2.")
  }
  
  ## CROPS is not available for online version
  if (penalty == "CROPS") {
    stop("Use cpt.meanvar from changepoint as CROPS is not available with changepoint.online")
  }
  
  ## ------------------------------------------------------------
  ## Choose cost function name (Normal / Exp / Gamma / Poisson)
  ## ------------------------------------------------------------
  if (penalty == "MBIC") {
    if (test.stat == "Normal") {
      cost_func <- "meanvar.norm.mbic"
    } else if (test.stat == "Exponential") {
      cost_func <- "meanvar.exp.mbic"
    } else if (test.stat == "Gamma") {
      cost_func <- "meanvar.gamma.mbic"
    } else if (test.stat == "Poisson") {
      cost_func <- "meanvar.poisson.mbic"
    } else {
      stop("Not a valid test statistic. Must be Normal, Exponential, Gamma or Poisson")
    }
  } else {
    if (test.stat == "Normal") {
      cost_func <- "meanvar.norm"
    } else if (test.stat == "Exponential") {
      cost_func <- "meanvar.exp"
    } else if (test.stat == "Gamma") {
      cost_func <- "meanvar.gamma"
    } else if (test.stat == "Poisson") {
      cost_func <- "meanvar.poisson"
    } else {
      stop("Not a valid test statistic. Must be Normal, Exponential, Gamma or Poisson")
    }
  }
  
  ## ------------------------------------------------------------
  ## Build sufficient statistics matrix
  ## ------------------------------------------------------------
  y  <- coredata(data)
  mu <- mean(y)
  sumstat <- cbind(
    c(0, cumsum(y)),
    c(0, cumsum(y^2)),
    cumsum(c(0, (y - mu)^2))
  )
  
  ## ------------------------------------------------------------
  ## Map (penalty, pen.value) -> numeric penalty (includes ARL)
  ## ------------------------------------------------------------
  pen.value <- online.decision(
    penalty   = penalty,
    pen.value = pen.value,
    n         = length(data),
    diffparam = 1,
    asymcheck = cost_func,
    method    = "AMOC"
  )
  
  method <- "PELT"
  ans <- PELT.online.initialise(
    sumstat   = sumstat,
    pen       = pen.value,
    cost_func = cost_func,
    shape     = shape,
    minseglen = minseglen
  )
  
  ## ------------------------------------------------------------
  ## Wrap as ocpt object
  ## IMPORTANT: pass the *matrix* sumstat, not ans$sumstat
  ## ------------------------------------------------------------
  online.class_input(
    sumstat         = sumstat,
    cpttype         = "mean and variance",
    method          = method,
    test.stat       = test.stat,
    penalty         = penalty,
    pen.value       = ans$penalty,
    minseglen       = minseglen,
    param.estimates = param.estimates,
    out             = sort(ans$cptsout[ans$cptsout > 0]),
    shape           = ans$shape,
    Q               = Q,
    lastchangelike  = ans$lastchangelike,
    lastchangecpts  = ans$lastchangecpts,
    checklist       = ans$checklist,
    nchecklist      = ans$nchecklist,
    ndone           = ans$ndone,
    nupdate         = ans$nupdate,
    cost_func       = ans$cost_func
  )
}

ocpt.meanvar.initialize <- function(data,
                                    penalty         = "Manual",
                                    pen.value       = length(data),
                                    Q               = 5,
                                    test.stat       = "Normal",
                                    class           = TRUE,
                                    param.estimates = TRUE,
                                    shape           = 1,
                                    minseglen       = 2,
                                    alpha           = 1,
                                    verbose         = FALSE) {
  ocpt.meanvar.initialise(
    data            = data,
    penalty         = penalty,
    pen.value       = pen.value,
    Q               = Q,
    test.stat       = test.stat,
    class           = class,
    param.estimates = param.estimates,
    shape           = shape,
    minseglen       = minseglen,
    alpha           = alpha,
    verbose         = verbose
  )
}

ocpt.meanvar.update <- function(previousanswer, newdata) {
  checkData(newdata)
  
  ## ECP branch: unchanged
  if (class(previousanswer) == "ecp.ocpt") {
    ecpans <- e.cp3o_delta.online.update(
      previousanswer,
      newdata,
      K = 2 * (previousanswer@number)
    )
    return(ecpans)
  }
  
  param.estimates <- is.list(previousanswer@param.est)
  
  nextans <- PELT.online.update(
    previousanswer = previousanswer,
    newdata        = newdata
  )
  
  online.class_input(
    sumstat         = nextans$sumstat,
    cpttype         = "mean and variance",
    method          = previousanswer@method,
    test.stat       = previousanswer@test.stat,
    penalty         = previousanswer@pen.type,
    pen.value       = nextans$penalty,
    minseglen       = nextans$minseglen,
    param.estimates = param.estimates,
    out             = sort(nextans$cptsout[nextans$cptsout > 0]),
    shape           = previousanswer@shape,
    lastchangelike  = nextans$lastchangelike,
    lastchangecpts  = nextans$lastchangecpts,
    checklist       = nextans$checklist,
    nchecklist      = nextans$nchecklist,
    ndone           = nextans$ndone,
    nupdate         = nextans$nupdate,
    cost_func       = previousanswer@cost_func
  )
}