ocpt.var.initialise <- function(
    data,
    penalty         = "Manual",
    pen.value       = length(data),
    know.mean       = FALSE,
    mu              = NA,
    Q               = 5,
    test.stat       = "Normal",
    class           = TRUE,
    param.estimates = TRUE,
    shape           = 1,
    minseglen       = 1,
    alpha           = 1,
    verbose         = FALSE
) {
  
  checkData(data)

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
  
  ## Basic checks
  if (minseglen < 1) {
    minseglen <- 1
    warning("Minimum segment length for a change in variance is 1, automatically changed to 1.")
  }
  
  if (length(data) <= 2 * minseglen) {
    stop("Data length must be larger than 2*minseglen.")
  }
  
  if (penalty == "CROPS") {
    stop("Use cpt.var from changepoint as CROPS is not available with changepoint.online")
  }
  
  ## Cost function choice
  if (penalty == "MBIC") {
    if (test.stat == "Normal") {
      cost_func <- "var.norm.mbic"
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
      cost_func <- "var.norm"
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
  
  ## Known / unknown mean
  core <- coredata(data)
  if (!know.mean) {
    mu <- mean(core)
  }
  
  ## Sufficient statistics
  sumstat <- cbind(
    c(0, cumsum(core)),
    c(0, cumsum(core^2)),
    cumsum(c(0, (core - mu)^2))
  )
  
  ## Map (penalty, pen.value) -> numeric penalty via online.decision
  pen.num <- online.decision(
    penalty   = penalty,
    pen.value = pen.value,
    n         = length(data),
    diffparam = 1,
    asymcheck = cost_func,
    method    = "AMOC"   # matches original var code
  )
  
  method <- "PELT"
  
  ## Initialise online PELT with numeric penalty
  ans <- PELT.online.initialise(
    sumstat   = sumstat,
    pen       = pen.num,
    cost_func = cost_func,
    shape     = shape,
    minseglen = minseglen
  )
  
  ## Build & return ocpt object
  return(
    online.class_input(
      sumstat        = sumstat,
      cpttype        = "variance",
      method         = method,
      test.stat      = test.stat,
      penalty        = penalty,        
      pen.value      = ans$penalty,    
      minseglen      = minseglen,
      param.estimates = param.estimates,
      out            = sort(ans$cptsout[ans$cptsout > 0]),
      shape          = shape,
      Q              = Q,
      lastchangelike = ans$lastchangelike,
      lastchangecpts = ans$lastchangecpts,
      checklist      = ans$checklist,
      nchecklist     = ans$nchecklist,
      ndone          = ans$ndone,
      nupdate        = ans$nupdate,
      cost_func      = ans$cost_func
    )
  )
}

## Backwards-compatible alias
ocpt.var.initialize <- function(
    data,
    penalty         = "Manual",
    pen.value       = length(data),
    know.mean       = FALSE,
    mu              = NA,
    Q               = 5,
    test.stat       = "Normal",
    class           = TRUE,
    param.estimates = TRUE,
    shape           = 1,
    minseglen       = 1,
    alpha           = 1,
    verbose         = FALSE
) {
  ocpt.var.initialise(
    data            = data,
    penalty         = penalty,
    pen.value       = pen.value,
    know.mean       = know.mean,
    mu              = mu,
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

## Update method
ocpt.var.update <- function(previousanswer, newdata) {
  
  checkData(newdata)
  
  if (class(previousanswer) == "ecp.ocpt") {
    ecpans <- e.cp3o_delta.online.update(
      previousanswer,
      newdata,
      K = 2 * (previousanswer@number)
    )
    return(ecpans)
  }
  
  if (class(previousanswer@param.est) == "list") {
    param.estimates <- TRUE
  } else {
    param.estimates <- FALSE
  }
  
  nextans <- PELT.online.update(previousanswer = previousanswer,
                                newdata        = newdata)
  
  return(
    online.class_input(
      sumstat        = nextans$sumstat,
      cpttype        = "variance",
      method         = previousanswer@method,
      test.stat      = previousanswer@test.stat,
      penalty        = previousanswer@pen.type,
      pen.value      = nextans$penalty,
      minseglen      = nextans$minseglen,
      param.estimates = param.estimates,
      out            = sort(nextans$cptsout[nextans$cptsout > 0]),
      shape          = previousanswer@shape,
      lastchangelike = nextans$lastchangelike,
      lastchangecpts = nextans$lastchangecpts,
      checklist      = nextans$checklist,
      nchecklist     = nextans$nchecklist,
      ndone          = nextans$ndone,
      nupdate        = nextans$nupdate,
      cost_func      = previousanswer@cost_func
    )
  )
}