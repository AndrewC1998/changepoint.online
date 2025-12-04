online.decision <- function(penalty,
                            pen.value,
                            n,
                            diffparam,
                            asymcheck,
                            method) {
  
  if (penalty == "NoPenalty") penalty <- "None"
  
  ## ARL penalty
  if (penalty == "ARL") {
    if (!is.numeric(pen.value) || length(pen.value) != 1L || pen.value <= 1) {
      stop("For penalty='ARL', pen.value must be a single numeric ARL > 1.")
    }
    runlen <- pen.value
    alpha  <- 1 / runlen
  
    asym_name <- sub("\\.mbic$", "", asymcheck)
    
    pen_asym <- suppressWarnings(
      online.decision(
        penalty   = "Asymptotic",
        pen.value = alpha,
        n         = n,
        diffparam = diffparam,
        asymcheck = asym_name,
        method    = method
      )
    )
    
    if (!is.finite(pen_asym)) {
      warning("ARL to penalty calibration was numerically unstable, ",
              "falling back to log(ARL).")
      pen.return <- log(runlen)
    } else {
      pen.return <- pen_asym
    }
    
    return(pen.return)
  }

  ## Existing logic for standard penalties + 'Manual' + 'Asymptotic'
  if ((penalty == "SIC0") || (penalty == "BIC0")) {
    pen.return <- diffparam * log(n)
  } else if ((penalty == "SIC") || (penalty == "BIC")) {
    pen.return <- (diffparam + 1) * log(n)
  } else if (penalty == "MBIC") {
    pen.return <- (diffparam + 2) * log(n)
  } else if ((penalty == "SIC1") || (penalty == "BIC1")) {
    stop("SIC1 and BIC1 have been depreciated, use SIC or BIC for the same result.")
  } else if (penalty == "AIC0") {
    pen.return <- 2 * diffparam
  } else if (penalty == "AIC") {
    pen.return <- 2 * (diffparam + 1)
  } else if (penalty == "AIC1") {
    stop("AIC1 has been depreciated, use AIC for the same result.")
  } else if (penalty == "Hannan-Quinn0") {
    pen.return <- 2 * diffparam * log(log(n))
  } else if (penalty == "Hannan-Quinn") {
    pen.return <- 2 * (diffparam + 1) * log(log(n))
  } else if (penalty == "Hannan-Quinn1") {
    stop("Hannan-Quinn1 has been depreciated, use Hannan-Quinn for the same result.")
  } else if (penalty == "None") {
    pen.return <- 0
  } else if ((penalty != "Manual") && (penalty != "Asymptotic")) {
    stop("Unknown penalty type.")
  }
  
  ## Manual penalty 
  if (penalty == "Manual" && !is.numeric(pen.value)) {
    pen.value <- try(eval(parse(text = paste(pen.value))), silent = TRUE)
    if (inherits(pen.value, "try-error")) {
      stop("Your manual penalty cannot be evaluated")
    } else {
      pen.return <- pen.value
    }
  }
  if (penalty == "Manual" && is.numeric(pen.value)) {
    pen.return <- pen.value
  }
  
  ## Asymptotic penalties
  if (penalty == "Asymptotic") {
    if ((pen.value <= 0) || (pen.value > 1)) {
      stop("Asymptotic penalty values must be > 0 and <= 1")
    }
    if (method != "AMOC") {
      warning("Asymptotic penalty value is not accurate for multiple changes, ",
              "it should be treated the same as a manual penalty choice.")
    }
    if (asymcheck == "mean.norm") {
      alpha <- pen.value
      alogn <- (2 * log(log(n)))^(-1/2)
      blogn <- (alogn^(-1)) + 0.5 * alogn * log(log(log(n)))
      pen.return <- (-alogn *
                       log(log((1 - alpha +
                                  exp(-2 * (pi^(1/2)) *
                                        exp(blogn / alogn)))^
                                 (-1 / (2 * (pi^(1/2)))))) +
                       blogn)^2
    } else if (asymcheck == "var.norm") {
      alpha <- pen.value
      alogn <- sqrt(2 * log(log(n)))
      blogn <- 2 * log(log(n)) +
        (log(log(log(n)))) / 2 - log(gamma(1 / 2))
      pen.return <- (-(log(log((1 - alpha +
                                  exp(-2 * exp(blogn)))^(-1 / 2)))) / alogn +
                       blogn / alogn)^2
    } else if (asymcheck == "meanvar.norm") {
      alpha <- pen.value
      alogn <- sqrt(2 * log(log(n)))
      blogn <- 2 * log(log(n)) + log(log(log(n)))
      pen.return <- (-(log(log((1 - alpha +
                                  exp(-2 * exp(blogn)))^(-1 / 2)))) / alogn +
                       blogn / alogn)^2
    } else if (asymcheck == "reg.norm") {
      alpha <- pen.value
      top <- -(log(log((1 - alpha +
                          exp(-2 * exp(2 * (log(log(n))) + (diffparam / 2) *
                                         (log(log(log(n)))) -
                                         log(gamma(diffparam / 2)))))^(-1 / 2)))) +
        2 * (log(log(n))) + (diffparam / 2) * (log(log(log(n)))) -
        log(gamma(diffparam / 2))
      bottom <- (2 * log(log(n)))^(1 / 2)
      pen.return <- (top / bottom)^2
    } else if (asymcheck == "var.css") {
      if (pen.value == 0.01) pen.return <- 1.628
      else if (pen.value == 0.05) pen.return <- 1.358
      else if (pen.value == 0.1) pen.return <- 1.224
      else if (pen.value == 0.25) pen.return <- 1.019
      else if (pen.value == 0.5) pen.return <- 0.828
      else if (pen.value == 0.75) pen.return <- 0.677
      else if (pen.value == 0.9) pen.return <- 0.571
      else if (pen.value == 0.95) pen.return <- 0.520
      else stop("Only alpha values of 0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95 are valid for CSS")
    } else if (asymcheck == "mean.cusum") {
      stop("Asymptotic penalties have not been implemented yet for CUSUM")
    } else if (asymcheck == "meanvar.gamma") {
      stop("Asymptotic penalties for the Gamma test statistic are not defined, please choose an alternative penalty type")
    } else if (asymcheck == "meanvar.exp") {
      alpha <- pen.value
      an <- (2 * log(log(n)))^(1 / 2)
      bn <- 2 * log(log(n)) + 0.5 * log(log(log(n))) - 0.5 * log(pi)
      pen.return <- (-1 / an) * log(-0.5 * log(1 - alpha)) + bn
      if (alpha == 1) pen.return <- 1.42417
    } else if (asymcheck == "meanvar.poisson") {
      stop("Asymptotic penalties for the Poisson test statistic are not available yet, please choose an alternative penalty type")
    }
  }
  
  ## sanity check
  if (is.finite(pen.return) && pen.return < 0) {
    stop("pen.value cannot be negative, please change your penalty value")
  }
  
  return(pen.return)
}