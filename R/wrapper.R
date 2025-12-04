expand_lastchange <- function(lastlike, lastcpt, N) {
  # Expand compact representations of lastchangelike / lastchangecpts
  # into full vectors of length N.
  #
  # lastlike, lastcpt can be:
  #  - 2-column matrices with columns c("like","idx") or c("cpt","idx")
  #  - full numeric vectors (for backwards compatibility)
  
  lastchangelike_full <- numeric(N)
  lastchangecpts_full <- integer(N)
  
  ## lastchangelike
  if (is.matrix(lastlike)) {
    if (nrow(lastlike) > 0L) {
      idx <- as.integer(lastlike[, "idx"])
      lastchangelike_full[idx] <- lastlike[, "like"]
    }
  } else if (is.numeric(lastlike)) {
    lastchangelike_full[seq_len(min(length(lastlike), N))] <- lastlike
  }
  
  ## lastchangecpts
  if (is.matrix(lastcpt)) {
    if (nrow(lastcpt) > 0L) {
      idx <- as.integer(lastcpt[, "idx"])
      lastchangecpts_full[idx] <- as.integer(lastcpt[, "cpt"])
    }
  } else if (is.numeric(lastcpt)) {
    lastchangecpts_full[seq_len(min(length(lastcpt), N))] <- as.integer(lastcpt)
  }
  
  list(
    lastchangelike = lastchangelike_full,
    lastchangecpts = lastchangecpts_full
  )
}

compress_lastchange <- function(lastchangelike_full,
                                lastchangecpts_full,
                                checklist_full,
                                cpts_full,
                                ndone,
                                nupdate) {
  # Tree-aware compression of the DP state.
  #  - Keep:
  #      * new data indices (latest block),
  #      * all checklist entries,
  #      * all current changepoints,
  #      * and their ancestors in the parent tree.
  #  - Return compact 2-col matrices for lastchangelike / lastchangecpts,
  #    plus a trimmed checklist of length N, and nchecklist.
  
  N <- ndone + nupdate + 1L
  
  if (length(lastchangelike_full) != N || length(lastchangecpts_full) != N) {
    stop("compress_lastchange: length mismatch in DP state.")
  }
  
  keep <- rep(FALSE, N)
  
  ## New data region indices (completed indices: ndone+1 .. ndone+nupdate)
  new_start <- ndone + 1L
  new_end   <- ndone + nupdate
  if (new_end >= new_start) {
    keep[new_start:new_end] <- TRUE
  }
  
  ## Checklist entries (where checklist_full > 0)
  if (length(checklist_full) == N) {
    keep[which(checklist_full > 0L)] <- TRUE
  }
  
  ## Current changepoints (non-zero entries in cpts_full)
  current_cpts <- cpts_full[cpts_full > 0L]
  if (length(current_cpts)) {
    keep[current_cpts] <- TRUE
  }
  
  ## Close under the parent-pointer tree: follow lastchangecpts_full[j]
  ## backwards until 0, keeping everything on each path.
  for (t in which(keep)) {
    j <- lastchangecpts_full[t]
    while (j > 0L && !keep[j]) {
      keep[j] <- TRUE
      j <- lastchangecpts_full[j]
    }
  }
  
  keep_idx <- which(keep)
  
  lastlike_mat_trim <- cbind(
    like = lastchangelike_full[keep_idx],
    idx  = keep_idx
  )
  colnames(lastlike_mat_trim) <- c("like", "idx")
  
  lastcpt_mat_trim <- cbind(
    cpt = lastchangecpts_full[keep_idx],
    idx = keep_idx
  )
  colnames(lastcpt_mat_trim) <- c("cpt", "idx")
  
  checklist_trim <- integer(N)
  if (length(checklist_full) == N) {
    checklist_trim[keep_idx] <- checklist_full[keep_idx]
  }
  nchecklist <- sum(checklist_trim > 0L)
  
  list(
    lastchangelike = lastlike_mat_trim,
    lastchangecpts = lastcpt_mat_trim,
    checklist      = checklist_trim,
    nchecklist     = nchecklist
  )
}

# PELT.online.initialise

PELT.online.initialise <- function(sumstat,
                                   pen       = 0,
                                   cost_func = "mean.norm",
                                   shape     = 1,
                                   minseglen = 1) {
  # Initialisation for online PELT.
  # sumstat: n x 3 matrix of sufficient statistics
  #          (with first row = 0s for AMOC-style cumulative representation)
  
  if (is.null(dim(sumstat))) {
    sumstat <- as.matrix(sumstat, ncol = 1L)
  }
  
  ndone   <- 0L
  nupdate <- nrow(sumstat) - 1L  # data points (excluding initial 0-row)
  
  lastchangelike <- rep(0, nupdate + ndone + 1L)
  storage.mode(lastchangelike) <- "double"
  
  lastchangecpts <- rep(0L, nupdate + ndone + 1L)
  storage.mode(lastchangecpts) <- "integer"
  
  checklist <- rep(0L, nupdate + ndone + 1L)
  storage.mode(checklist) <- "integer"
  
  nchecklist <- 0L
  
  cptsout <- rep(0L, ndone + nupdate + 1L)
  storage.mode(cptsout) <- "integer"
  
  answer <- list()
  answer[[7]] <- 1L
  on.exit(.C("FreePELT", answer[[7]]))
  
  error <- 0L
  
  answer <- .C(
    "PELT_online",
    cost_func      = as.character(cost_func),
    sumstat        = as.double(sumstat),
    ndone          = as.integer(ndone),
    nupdate        = as.integer(nupdate),
    penalty        = as.double(pen),
    cptsout        = cptsout,
    error          = as.integer(error),
    shape          = as.double(shape),
    minseglen      = as.integer(minseglen),
    lastchangelike = lastchangelike,
    lastchangecpts = lastchangecpts,
    checklist      = checklist,
    nchecklist     = as.integer(nchecklist)
  )
  
  names(answer) <- c(
    "cost_func", "sumstat", "ndone", "nupdate",
    "penalty", "cptsout", "error", "shape", "minseglen",
    "lastchangelike", "lastchangecpts", "checklist", "nchecklist"
  )
  
  if (answer$error > 0L) {
    stop("C code error:", answer$error, call. = FALSE)
  }
  
  ## Compress DP state into 2-col matrices + trimmed checklist
  comp <- compress_lastchange(
    lastchangelike_full = answer$lastchangelike,
    lastchangecpts_full = answer$lastchangecpts,
    checklist_full      = answer$checklist,
    cpts_full           = answer$cptsout,
    ndone               = answer$ndone,
    nupdate             = answer$nupdate
  )
  
  answer$lastchangelike <- comp$lastchangelike
  answer$lastchangecpts <- comp$lastchangecpts
  answer$checklist      <- comp$checklist
  answer$nchecklist     <- as.integer(comp$nchecklist)
  
  return(answer)
}

# PELT.online.update

PELT.online.update <- function(previousanswer, newdata) {
  # Online update step for PELT.
  # previousanswer: "ocpt" S4 object (or compatible list from initialiser)
  # newdata: numeric vector of new observations
  
  checkData(newdata)
  
  ## Completed length so far
  ndone_prev   <- previousanswer@ndone + previousanswer@nupdate
  nupdate      <- length(newdata)
  ndone        <- ndone_prev           # for clarity; C will see this as "completed"
  
  ## -----------------------------------
  ## Update sumstat with the new data
  ## -----------------------------------
  y_new <- coredata(newdata)
  mu    <- mean(y_new)
  
  end.sumstat <- previousanswer@sumstat[nrow(previousanswer@sumstat), ]
  
  sumstat_new <- rbind(
    previousanswer@sumstat,
    matrix(
      c(
        end.sumstat[1] + cumsum(y_new),
        end.sumstat[2] + cumsum(y_new^2),
        end.sumstat[3] + cumsum((y_new - mu)^2)
      ),
      ncol = 3L
    )
  )
  storage.mode(sumstat_new) <- "double"
  
  ## -----------------------------------------------
  ## Expand lastchangelike / lastchangecpts to full
  ## vectors of length ndone_prev + 1
  ## -----------------------------------------------
  N_prev   <- ndone_prev + 1L
  expanded <- expand_lastchange(
    lastlike = previousanswer@lastchangelike,
    lastcpt  = previousanswer@lastchangecpts,
    N        = N_prev
  )
  
  full_ll_prev  <- expanded$lastchangelike
  full_cpt_prev <- expanded$lastchangecpts
  
  ## Extend for the new data
  lastchangelike <- c(full_ll_prev, numeric(nupdate))
  storage.mode(lastchangelike) <- "double"
  
  lastchangecpts <- c(full_cpt_prev, integer(nupdate))
  storage.mode(lastchangecpts) <- "integer"
  
  ## Checklist: keep entire previous vector, then append zeros for new data
  checklist <- c(previousanswer@checklist, integer(nupdate))
  storage.mode(checklist) <- "integer"
  nchecklist <- as.integer(sum(checklist > 0L))
  
  ## Output changepoint vector for C
  cptsout <- rep(0L, ndone + nupdate + 1L)
  storage.mode(cptsout) <- "integer"
  
  error  <- 0L
  answer <- list()
  answer[[7]] <- 1L
  on.exit(.C("FreePELT", answer[[7]]))
  
  answer <- .C(
    "PELT_online",
    cost_func      = as.character(previousanswer@cost_func),
    sumstat        = sumstat_new,
    ndone          = as.integer(ndone),
    nupdate        = as.integer(nupdate),
    penalty        = as.double(previousanswer@pen.value),
    cptsout        = cptsout,
    error          = as.integer(error),
    shape          = previousanswer@shape,
    minseglen      = as.integer(previousanswer@minseglen),
    lastchangelike = lastchangelike,
    lastchangecpts = lastchangecpts,
    checklist      = checklist,
    nchecklist     = nchecklist
  )
  
  if (answer$error > 0L) {
    stop("C code error:", answer$error, call. = FALSE)
  }
  
  names(answer) <- c(
    "cost_func", "sumstat", "ndone", "nupdate",
    "penalty", "cptsout", "error", "shape", "minseglen",
    "lastchangelike", "lastchangecpts", "checklist", "nchecklist"
  )
  
  ## Compress DP state again
  comp <- compress_lastchange(
    lastchangelike_full = answer$lastchangelike,
    lastchangecpts_full = answer$lastchangecpts,
    checklist_full      = answer$checklist,
    cpts_full           = answer$cptsout,
    ndone               = answer$ndone,
    nupdate             = answer$nupdate
  )
  
  answer$lastchangelike <- comp$lastchangelike
  answer$lastchangecpts <- comp$lastchangecpts
  answer$checklist      <- comp$checklist
  answer$nchecklist     <- as.integer(comp$nchecklist)
  
  return(answer)
}