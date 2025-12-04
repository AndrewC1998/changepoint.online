.compress_lastchange <- function(lastchangelike,
                                 lastchangecpts,
                                 checklist,
                                 nchecklist,
                                 ndone,
                                 nupdate) {
  # If we don't have a meaningful nupdate yet, just return empty matrices
  total_len <- length(lastchangelike)
  if (total_len == 0L) {
    return(list(
      lastchangelike = matrix(numeric(0), ncol = 2,
                              dimnames = list(NULL, c("like", "idx"))),
      lastchangecpts = matrix(integer(0), ncol = 2,
                              dimnames = list(NULL, c("cpt", "idx")))
    ))
  }
  
  # Time index runs 1:(ndone + nupdate + 1) in the C code
  total_T <- ndone + nupdate + 1L
  total_T <- min(total_T, total_len)  # safety if things are slightly off
  
  # Indices we care about:
  #  - 1 (initial)
  #  - checklist entries (where > 0)
  #  - final index total_T
  cl_valid <- integer(0)
  if (nchecklist > 0L && length(checklist) > 0L) {
    cl_valid <- as.integer(checklist[seq_len(nchecklist)])
    cl_valid <- cl_valid[cl_valid >= 1L & cl_valid <= total_T]
  }
  keep_idx <- sort(unique(c(1L, cl_valid, total_T)))
  
  # Guard against weirdness
  keep_idx <- keep_idx[keep_idx >= 1L & keep_idx <= total_len]
  
  # Build matrices
  ll_mat <- cbind(
    like = as.numeric(lastchangelike[keep_idx]),
    idx  = keep_idx
  )
  cpt_mat <- cbind(
    cpt = as.integer(lastchangecpts[keep_idx]),
    idx = keep_idx
  )
  
  colnames(ll_mat)   <- c("like", "idx")
  colnames(cpt_mat)  <- c("cpt",  "idx")
  
  list(
    lastchangelike = ll_mat,
    lastchangecpts = cpt_mat
  )
}

online.class_input <- function(sumstat, cpttype, method, test.stat,
                               penalty, pen.value, minseglen,
                               param.estimates, out = list(),
                               Q = NA, shape = NA,
                               lastchangelike = c(0),
                               lastchangecpts = c(0),
                               checklist = c(0),
                               nchecklist = 0,
                               ndone = 0,
                               nupdate = length(data),
                               cost_func) {
  
  ## ------------------------------------------------------------
  ## Ensure lastchangelike / lastchangecpts are 2-col matrices
  ##   col1 = value (likelihood or cpt)
  ##   col2 = original index (for later trimming/reconstruction)
  ## ------------------------------------------------------------
  if (is.null(dim(lastchangelike))) {
    idx_lik <- seq_along(lastchangelike)
    lastchangelike <- cbind(
      like = as.numeric(lastchangelike),
      idx  = as.integer(idx_lik)
    )
  }
  
  if (is.null(dim(lastchangecpts))) {
    idx_cp <- seq_along(lastchangecpts)
    lastchangecpts <- cbind(
      cpt = as.integer(lastchangecpts),
      idx = as.integer(idx_cp)
    )
  }
  
  ## ------------------------------------------------------------
  ## Choose class (ocpt.range vs ocpt)
  ## ------------------------------------------------------------
  if (method == "BinSeg" || method == "SegNeigh" || penalty == "CROPS") {
    ans <- new("ocpt.range")
  } else {
    ans <- new("ocpt")
  }
  
  ## core slots (using accessors where they exist)
  sumstat(ans)        <- sumstat
  cpttype(ans)        <- cpttype
  method(ans)         <- method
  test.stat(ans)      <- test.stat
  pen.type(ans)       <- penalty
  pen.value(ans)      <- pen.value
  minseglen(ans)      <- minseglen
  ndone(ans)     <- ndone
  nupdate(ans)   <- nupdate
  checklist(ans) <- checklist
  ans@nchecklist <- nchecklist
  lc <- .compress_lastchange(
    lastchangelike = lastchangelike,
    lastchangecpts = lastchangecpts,
    checklist      = checklist,
    nchecklist     = nchecklist,
    ndone          = ndone,
    nupdate        = nupdate
  )
  lastchangelike(ans) <- lc$lastchangelike
  lastchangecpts(ans) <- lc$lastchangecpts
  shape(ans)          <- shape
  cost_func(ans)      <- cost_func
  
  ## ------------------------------------------------------------
  ## Single optimal segmentation (non-CROPS)
  ## ------------------------------------------------------------
  if (penalty != "CROPS") {  # CROPS is the only one that doesn't give a single set of cpts
    cpts(ans) <- out
    
    if (isTRUE(param.estimates)) {
      if (test.stat(ans) == "Gamma") {
        ans <- param(ans, shape)
      } else {
        ans <- param(ans)
      }
    }
  }
  
  ## ------------------------------------------------------------
  ## ncpts.max depending on method
  ## ------------------------------------------------------------
  if (method == "PELT") {
    ncpts.max(ans) <- Inf
  } else if (method == "AMOC") {
    ncpts.max(ans) <- 1
  } else {
    ncpts.max(ans) <- Q
  }
  
  ## ------------------------------------------------------------
  ## Range methods (BinSeg / SegNeigh / CROPS)
  ## ------------------------------------------------------------
  if (method == "BinSeg") {
    l <- list()
    for (i in 1:(length(out$cps) / 2)) {
      l[[i]] <- out$cps[1, 1:i]
    }
    m <- t(sapply(l, "[", 1:max(sapply(l, length))))
    
    cpts.full(ans)     <- m
    pen.value.full(ans) <- out$cps[2, ]
    
  } else if (method == "SegNeigh") {
    
    cpts.full(ans)      <- out$cps[-1, ]
    pen.value.full(ans) <- -diff(out$like.Q)
    
  } else if (penalty == "CROPS") {
    
    m <- t(sapply(out[[2]], "[", 1:max(sapply(out[[2]], length))))
    
    cpts.full(ans)      <- m
    pen.value.full(ans) <- out[[1]][1, ]
    if (test.stat(ans) == "Gamma") param.est(ans)$shape <- shape
  }
  
  return(ans)
}

online.ecp.class_input <- function(number, estimates, GofM, delta, alpha, verbose, csum, dll, dlr, drr, left, right, datalength, functime, width, cpLoc){
    
    ans = new("ecp.ocpt")
    
    number(ans)=number; estimates(ans)=estimates; GofM(ans)=GofM; delta(ans)=delta; alpha(ans)=alpha; verbose(ans)=verbose; csum(ans)=csum; dll(ans)=dll; dlr(ans)=dlr; drr(ans)=drr; left(ans)=left; right(ans)=right; datalength(ans)=datalength; functime(ans)=functime; width(ans)=width; cpLoc(ans)=cpLoc;
    return(ans)
}

