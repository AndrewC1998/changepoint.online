## ============================================================
## S4 classes and methods for online changepoint objects
## ============================================================

## Core classes

setClass("ocpt",
         slots = list(
           sumstat        = "array",
           cpttype        = "character",
           method         = "character",
           test.stat      = "character",
           pen.type       = "character",
           pen.value      = "numeric",
           minseglen      = "numeric",
           cpts           = "numeric",
           ncpts.max      = "numeric",
           param.est      = "list",
           date           = "character",
           version        = "character",
           lastchangelike = "matrix",   
           lastchangecpts = "matrix",   
           checklist      = "numeric",
           nchecklist     = "numeric",
           ndone          = "numeric",
           nupdate        = "numeric",
           cost_func      = "character",
           shape          = "numeric"
         ),
         prototype = prototype(
           cpttype        = "Not Set",
           date           = date(),
           version        = as(packageVersion("changepoint.online"), "character")
         )
)

setClass("ocpt.reg",
         slots = list(
           sumstat        = "matrix",
           cpttype        = "character",
           method         = "character",
           test.stat      = "character",
           pen.type       = "character",
           pen.value      = "numeric",
           minseglen      = "numeric",
           cpts           = "numeric",
           ncpts.max      = "numeric",
           param.est      = "list",
           date           = "character",
           version        = "character",
           lastchangelike = "matrix",
           lastchangecpts = "matrix",
           checklist      = "numeric",
           nchecklist     = "numeric",
           ndone          = "numeric",
           nupdate        = "numeric",
           cost_func      = "character",
           shape          = "numeric"
         ),
         prototype = prototype(
           cpttype        = "Not Set",
           date           = date(),
           version        = as(packageVersion("changepoint.online"), "character")
         )
)

setClass(
  "ocpt.range",
  slots    = list(cpts.full = "matrix", pen.value.full = "numeric"),
  prototype = prototype(),
  contains  = "ocpt"
)

## ecp.ocpt class (extends ocpt)

setClass("ecp.ocpt",
         slots = list(
           number     = "numeric",
           estimates  = "numeric",
           GofM       = "numeric",
           delta      = "numeric",
           alpha      = "numeric",
           verbose    = "logical",
           csum       = "numeric",
           dll        = "numeric",
           dlr        = "numeric",
           drr        = "numeric",
           left       = "matrix",
           right      = "matrix",
           datalength = "numeric",
           functime   = "numeric",
           width      = "numeric",
           cpLoc      = "list"
         ),
         contains = "ocpt"
)

## Generics

## sumstat
setGeneric("sumstat", function(object, ...) standardGeneric("sumstat"))

setMethod("sumstat", "ocpt",     function(object, ...) object@sumstat)
setMethod("sumstat", "ocpt.reg", function(object, ...) object@sumstat)

if (!isGeneric("cpttype")) {
  setGeneric("cpttype", function(object) standardGeneric("cpttype"))
}
setMethod("cpttype", "ocpt",     function(object) object@cpttype)
setMethod("cpttype", "ocpt.reg", function(object) object@cpttype)

if (!isGeneric("method")) {
  setGeneric("method", function(object) standardGeneric("method"))
}
setMethod("method", "ocpt",     function(object) object@method)
setMethod("method", "ocpt.reg", function(object) object@method)

## distribution
if (!isGeneric("distribution")) {
  setGeneric("distribution", function(object) standardGeneric("distribution"))
}
setMethod("distribution", "ocpt",     function(object) object@test.stat)
setMethod("distribution", "ocpt.reg", function(object) object@test.stat)

if (!isGeneric("test.stat")) {
  setGeneric("test.stat", function(object) standardGeneric("test.stat"))
}
setMethod("test.stat", "ocpt",     function(object) object@test.stat)
setMethod("test.stat", "ocpt.reg", function(object) object@test.stat)

if (!isGeneric("pen.type")) {
  setGeneric("pen.type", function(object) standardGeneric("pen.type"))
}
setMethod("pen.type", "ocpt",     function(object) object@pen.type)
setMethod("pen.type", "ocpt.reg", function(object) object@pen.type)

if (!isGeneric("pen.value")) {
  setGeneric("pen.value", function(object) standardGeneric("pen.value"))
}
setMethod("pen.value", "ocpt",     function(object) object@pen.value)
setMethod("pen.value", "ocpt.reg", function(object) object@pen.value)

if (!isGeneric("pen.value.full")) {
  setGeneric("pen.value.full", function(object) standardGeneric("pen.value.full"))
}
setMethod("pen.value.full", "ocpt.range", function(object) object@pen.value.full)

if (!isGeneric("minseglen")) {
  setGeneric("minseglen", function(object) standardGeneric("minseglen"))
}
setMethod("minseglen", "ocpt", function(object) object@minseglen)

if (!isGeneric("cpts")) {
  setGeneric("cpts", function(object) standardGeneric("cpts"))
}
setMethod("cpts", "ocpt",     function(object) object@cpts[-length(object@cpts)])
setMethod("cpts", "ocpt.reg", function(object) object@cpts[-length(object@cpts)])

if (!isGeneric("cpts.full")) {
  setGeneric("cpts.full", function(object) standardGeneric("cpts.full"))
}
setMethod("cpts.full", "ocpt.range", function(object) object@cpts.full)

if (!isGeneric("ncpts.max")) {
  setGeneric("ncpts.max", function(object) standardGeneric("ncpts.max"))
}
setMethod("ncpts.max", "ocpt",     function(object) object@ncpts.max)
setMethod("ncpts.max", "ocpt.reg", function(object) object@ncpts.max)

if (!isGeneric("param.est")) {
  setGeneric("param.est", function(object) standardGeneric("param.est"))
}
setMethod("param.est", "ocpt",     function(object) object@param.est)
setMethod("param.est", "ocpt.reg", function(object) object@param.est)

setMethod("coef", "ocpt",     function(object, ...) object@param.est)
setMethod("coef", "ocpt.reg", function(object, ...) object@param.est)

## ncpts
if (!isGeneric("ncpts")) {
  setGeneric("ncpts", function(object) standardGeneric("ncpts"))
}
setMethod("ncpts", "ocpt",     function(object) length(cpts(object)))
setMethod("ncpts", "ocpt.reg", function(object) length(cpts(object)))

## seg.len
if (!isGeneric("seg.len")) {
  setGeneric("seg.len", function(object) standardGeneric("seg.len"))
}
setMethod("seg.len", "ocpt",
          function(object) object@cpts - c(0, object@cpts[-length(object@cpts)]))
setMethod("seg.len", "ocpt.reg",
          function(object) object@cpts - c(0, object@cpts[-length(object@cpts)]))

## nseg
if (!isGeneric("nseg")) {
  setGeneric("nseg", function(object) standardGeneric("nseg"))
}
setMethod("nseg", "ocpt",     function(object) ncpts(object) + 1L)
setMethod("nseg", "ocpt.reg", function(object) ncpts(object) + 1L)

## lastchangelike

if (!isGeneric("lastchangelike")) {
  setGeneric("lastchangelike", function(object) standardGeneric("lastchangelike"))
}

setMethod(
  "lastchangelike",
  signature(object = "ocpt"),
  function(object) {
    ll <- object@lastchangelike
    if (is.matrix(ll) && ncol(ll) >= 1L) {
      if (!is.null(colnames(ll)) && "like" %in% colnames(ll)) {
        as.numeric(ll[, "like"])
      } else {
        as.numeric(ll[, 1L])
      }
    } else {
      as.numeric(ll)
    }
  }
)

setMethod(
  "lastchangelike",
  signature(object = "ocpt.reg"),
  function(object) {
    ll <- object@lastchangelike
    if (is.matrix(ll) && ncol(ll) >= 1L) {
      if (!is.null(colnames(ll)) && "like" %in% colnames(ll)) {
        as.numeric(ll[, "like"])
      } else {
        as.numeric(ll[, 1L])
      }
    } else {
      as.numeric(ll)
    }
  }
)

if (!isGeneric("lastchangelike<-")) {
  setGeneric("lastchangelike<-",
             function(object, value) standardGeneric("lastchangelike<-"))
}

setReplaceMethod(
  "lastchangelike",
  signature(object = "ocpt", value = "numeric"),
  function(object, value) {
    if (length(value) == 0L) {
      mat <- matrix(numeric(0), ncol = 2L)
      colnames(mat) <- c("like", "idx")
    } else {
      mat <- cbind(
        like = as.numeric(value),
        idx  = seq_along(value)
      )
    }
    object@lastchangelike <- mat
    validObject(object)
    object
  }
)

setReplaceMethod(
  "lastchangelike",
  signature(object = "ocpt.reg", value = "numeric"),
  function(object, value) {
    if (length(value) == 0L) {
      mat <- matrix(numeric(0), ncol = 2L)
      colnames(mat) <- c("like", "idx")
    } else {
      mat <- cbind(
        like = as.numeric(value),
        idx  = seq_along(value)
      )
    }
    object@lastchangelike <- mat
    validObject(object)
    object
  }
)

setReplaceMethod(
  "lastchangelike",
  signature(object = "ocpt", value = "matrix"),
  function(object, value) {
    mat <- value
    if (ncol(mat) == 1L) {
      mat <- cbind(
        like = as.numeric(mat[, 1L]),
        idx  = seq_len(nrow(mat))
      )
    } else {
      mat <- mat[, 1:2, drop = FALSE]
      colnames(mat) <- c("like", "idx")
      mat[, "like"] <- as.numeric(mat[, "like"])
      mat[, "idx"]  <- as.integer(mat[, "idx"])
    }
    object@lastchangelike <- mat
    validObject(object)
    object
  }
)

setReplaceMethod(
  "lastchangelike",
  signature(object = "ocpt.reg", value = "matrix"),
  function(object, value) {
    mat <- value
    if (ncol(mat) == 1L) {
      mat <- cbind(
        like = as.numeric(mat[, 1L]),
        idx  = seq_len(nrow(mat))
      )
    } else {
      mat <- mat[, 1:2, drop = FALSE]
      colnames(mat) <- c("like", "idx")
      mat[, "like"] <- as.numeric(mat[, "like"])
      mat[, "idx"]  <- as.integer(mat[, "idx"])
    }
    object@lastchangelike <- mat
    validObject(object)
    object
  }
)

## lastchangecpts

if (!isGeneric("lastchangecpts")) {
  setGeneric("lastchangecpts", function(object) standardGeneric("lastchangecpts"))
}

setMethod(
  "lastchangecpts",
  signature(object = "ocpt"),
  function(object) {
    lc <- object@lastchangecpts
    if (is.matrix(lc) && ncol(lc) >= 1L) {
      if (!is.null(colnames(lc)) && "cpt" %in% colnames(lc)) {
        as.numeric(lc[, "cpt"])
      } else {
        as.numeric(lc[, 1L])
      }
    } else {
      as.numeric(lc)
    }
  }
)

setMethod(
  "lastchangecpts",
  signature(object = "ocpt.reg"),
  function(object) {
    lc <- object@lastchangecpts
    if (is.matrix(lc) && ncol(lc) >= 1L) {
      if (!is.null(colnames(lc)) && "cpt" %in% colnames(lc)) {
        as.numeric(lc[, "cpt"])
      } else {
        as.numeric(lc[, 1L])
      }
    } else {
      as.numeric(lc)
    }
  }
)

if (!isGeneric("lastchangecpts<-")) {
  setGeneric("lastchangecpts<-",
             function(object, value) standardGeneric("lastchangecpts<-"))
}

setReplaceMethod(
  "lastchangecpts",
  signature(object = "ocpt", value = "numeric"),
  function(object, value) {
    if (length(value) == 0L) {
      mat <- matrix(integer(0), ncol = 2L)
      colnames(mat) <- c("cpt", "idx")
    } else {
      mat <- cbind(
        cpt = as.integer(value),
        idx = seq_along(value)
      )
    }
    object@lastchangecpts <- mat
    validObject(object)
    object
  }
)

setReplaceMethod(
  "lastchangecpts",
  signature(object = "ocpt.reg", value = "numeric"),
  function(object, value) {
    if (length(value) == 0L) {
      mat <- matrix(integer(0), ncol = 2L)
      colnames(mat) <- c("cpt", "idx")
    } else {
      mat <- cbind(
        cpt = as.integer(value),
        idx = seq_along(value)
      )
    }
    object@lastchangecpts <- mat
    validObject(object)
    object
  }
)

setReplaceMethod(
  "lastchangecpts",
  signature(object = "ocpt", value = "matrix"),
  function(object, value) {
    mat <- value
    if (ncol(mat) == 1L) {
      mat <- cbind(
        cpt = as.integer(mat[, 1L]),
        idx = seq_len(nrow(mat))
      )
    } else {
      mat <- mat[, 1:2, drop = FALSE]
      colnames(mat) <- c("cpt", "idx")
      mat[, "cpt"] <- as.integer(mat[, "cpt"])
      mat[, "idx"] <- as.integer(mat[, "idx"])
    }
    object@lastchangecpts <- mat
    validObject(object)
    object
  }
)

setReplaceMethod(
  "lastchangecpts",
  signature(object = "ocpt.reg", value = "matrix"),
  function(object, value) {
    mat <- value
    if (ncol(mat) == 1L) {
      mat <- cbind(
        cpt = as.integer(mat[, 1L]),
        idx = seq_len(nrow(mat))
      )
    } else {
      mat <- mat[, 1:2, drop = FALSE]
      colnames(mat) <- c("cpt", "idx")
      mat[, "cpt"] <- as.integer(mat[, "cpt"])
      mat[, "idx"] <- as.integer(mat[, "idx"])
    }
    object@lastchangecpts <- mat
    validObject(object)
    object
  }
)

## checklist
if (!isGeneric("checklist")) {
  setGeneric("checklist", function(object) standardGeneric("checklist"))
}
setMethod("checklist", "ocpt",     function(object) object@checklist)
setMethod("checklist", "ocpt.reg", function(object) object@checklist)

## ndone
if (!isGeneric("ndone")) {
  setGeneric("ndone", function(object) standardGeneric("ndone"))
}
setMethod("ndone", "ocpt",     function(object) object@ndone)
setMethod("ndone", "ocpt.reg", function(object) object@ndone)

## nupdate
if (!isGeneric("nupdate")) {
  setGeneric("nupdate", function(object) standardGeneric("nupdate"))
}
setMethod("nupdate", "ocpt",     function(object) object@nupdate)
setMethod("nupdate", "ocpt.reg", function(object) object@nupdate)

## cost_func (fixed fun <- cost_func)
if (!isGeneric("cost_func")) {
  setGeneric("cost_func", function(object) standardGeneric("cost_func"))
}
setMethod("cost_func", "ocpt",     function(object) object@cost_func)
setMethod("cost_func", "ocpt.reg", function(object) object@cost_func)

## shape
if (!isGeneric("shape")) {
  setGeneric("shape", function(object) standardGeneric("shape"))
}
setMethod("shape", "ocpt",     function(object) object@shape)
setMethod("shape", "ocpt.reg", function(object) object@shape)

## Replacement methods

setGeneric("sumstat<-", function(object, value) standardGeneric("sumstat<-"))
setReplaceMethod("sumstat", "ocpt",     function(object, value) { object@sumstat <- value; object })
setReplaceMethod("sumstat", "ocpt.reg", function(object, value) { object@sumstat <- value; object })

setGeneric("cpttype<-", function(object, value) standardGeneric("cpttype<-"))
setReplaceMethod("cpttype", "ocpt",     function(object, value) { object@cpttype <- value; object })
setReplaceMethod("cpttype", "ocpt.reg", function(object, value) { object@cpttype <- value; object })

setGeneric("method<-", function(object, value) standardGeneric("method<-"))
setReplaceMethod("method", "ocpt",     function(object, value) { object@method <- value; object })
setReplaceMethod("method", "ocpt.reg", function(object, value) { object@method <- value; object })

setGeneric("distribution<-", function(object, value) standardGeneric("distribution<-"))
setReplaceMethod("distribution", "ocpt",     function(object, value) { object@test.stat <- value; object })
setReplaceMethod("distribution", "ocpt.reg", function(object, value) { object@test.stat <- value; object })

setGeneric("test.stat<-", function(object, value) standardGeneric("test.stat<-"))
setReplaceMethod("test.stat", "ocpt",     function(object, value) { object@test.stat <- value; object })
setReplaceMethod("test.stat", "ocpt.reg", function(object, value) { object@test.stat <- value; object })

setGeneric("pen.type<-", function(object, value) standardGeneric("pen.type<-"))
setReplaceMethod("pen.type", "ocpt",     function(object, value) { object@pen.type <- value; object })
setReplaceMethod("pen.type", "ocpt.reg", function(object, value) { object@pen.type <- value; object })

setGeneric("pen.value<-", function(object, value) standardGeneric("pen.value<-"))
setReplaceMethod("pen.value", "ocpt",     function(object, value) { object@pen.value <- value; object })
setReplaceMethod("pen.value", "ocpt.reg", function(object, value) { object@pen.value <- value; object })

setGeneric("minseglen<-", function(object, value) standardGeneric("minseglen<-"))
setReplaceMethod("minseglen", "ocpt",       function(object, value) { object@minseglen <- value; object })
setReplaceMethod("minseglen", "ocpt.reg",   function(object, value) { object@minseglen <- value; object })
setReplaceMethod("minseglen", "ocpt.range", function(object, value) { object@minseglen <- value; object })

setGeneric("cpts<-", function(object, value) standardGeneric("cpts<-"))
setReplaceMethod("cpts", "ocpt",
                 function(object, value) {
                   if ((cpttype(object) == "meanar") || (cpttype(object) == "trendar")) {
                     n <- nrow(object@sumstat) - 2L
                   } else {
                     n <- nrow(object@sumstat) - 1L
                   }
                   if (value[length(value)] == n) {
                     object@cpts <- value
                   } else {
                     object@cpts <- c(value, n)
                   }
                   object
                 }
)
setReplaceMethod("cpts", "ocpt.reg",
                 function(object, value) {
                   n <- nrow(object@sumstat)
                   if (value[length(value)] == n) {
                     object@cpts <- value
                   } else {
                     object@cpts <- c(value, n)
                   }
                   object
                 }
)

setGeneric("ncpts.max<-", function(object, value) standardGeneric("ncpts.max<-"))
setReplaceMethod("ncpts.max", "ocpt",     function(object, value) { object@ncpts.max <- value; object })
setReplaceMethod("ncpts.max", "ocpt.reg", function(object, value) { object@ncpts.max <- value; object })

setGeneric("param.est<-", function(object, value) standardGeneric("param.est<-"))
setReplaceMethod("param.est", "ocpt",     function(object, value) { object@param.est <- value; object })
setReplaceMethod("param.est", "ocpt.reg", function(object, value) { object@param.est <- value; object })

setGeneric("cpts.full<-", function(object, value) standardGeneric("cpts.full<-"))
setReplaceMethod("cpts.full", "ocpt.range", function(object, value) { object@cpts.full <- value; object })

setGeneric("pen.value.full<-", function(object, value) standardGeneric("pen.value.full<-"))
setReplaceMethod("pen.value.full", "ocpt.range", function(object, value) { object@pen.value.full <- value; object })

setGeneric("checklist<-", function(object, value) standardGeneric("checklist<-"))
setReplaceMethod("checklist", "ocpt",     function(object, value) { object@checklist <- value; object })
setReplaceMethod("checklist", "ocpt.reg", function(object, value) { object@checklist <- value; object })

setGeneric("ndone<-", function(object, value) standardGeneric("ndone<-"))
setReplaceMethod("ndone", "ocpt",     function(object, value) { object@ndone <- value; object })
setReplaceMethod("ndone", "ocpt.reg", function(object, value) { object@ndone <- value; object })

setGeneric("nupdate<-", function(object, value) standardGeneric("nupdate<-"))
setReplaceMethod("nupdate", "ocpt",     function(object, value) { object@nupdate <- value; object })
setReplaceMethod("nupdate", "ocpt.reg", function(object, value) { object@nupdate <- value; object })

setGeneric("cost_func<-", function(object, value) standardGeneric("cost_func<-"))
setReplaceMethod("cost_func", "ocpt",     function(object, value) { object@cost_func <- value; object })
setReplaceMethod("cost_func", "ocpt.reg", function(object, value) { object@cost_func <- value; object })

setGeneric("shape<-", function(object, value) standardGeneric("shape<-"))
setReplaceMethod("shape", "ocpt",     function(object, value) { object@shape <- value; object })
setReplaceMethod("shape", "ocpt.reg", function(object, value) { object@shape <- value; object })

## param() generics and methods

setGeneric("param", function(object, ...) standardGeneric("param"))

setMethod("param", "ocpt", function(object, shape, ...) {
  
  param.mean <- function(object) {
    ss   <- sumstat(object)
    cpts <- c(0, object@cpts)
    nseg <- length(cpts) - 1L
    
    if (is.matrix(ss)) {
      S  <- ss[, 1L]
      mu <- numeric(nseg)
      for (j in seq_len(nseg)) {
        a   <- cpts[j]
        b   <- cpts[j + 1L]
        n   <- b - a
        Sy  <- S[b + 1L] - S[a + 1L]
        mu[j] <- Sy / n
      }
      mu
    } else {
      data <- ss
      mu   <- numeric(nseg)
      for (j in seq_len(nseg)) {
        mu[j] <- mean(data[(cpts[j] + 1L):cpts[j + 1L]])
      }
      mu
    }
  }
  
  param.var <- function(object) {
    ss   <- sumstat(object)
    cpts <- c(0, object@cpts)
    nseg <- length(cpts) - 1L
    
    if (is.matrix(ss)) {
      if (cpttype(object) == "variance") {
        Qmu <- ss[, 3L]
        v   <- numeric(nseg)
        for (j in seq_len(nseg)) {
          a   <- cpts[j]
          b   <- cpts[j + 1L]
          n   <- b - a
          SSE <- Qmu[b + 1L] - Qmu[a + 1L]
          v[j] <- SSE / n
        }
        v
      } else {
        S  <- ss[, 1L]
        Q  <- ss[, 2L]
        v  <- numeric(nseg)
        for (j in seq_len(nseg)) {
          a   <- cpts[j]
          b   <- cpts[j + 1L]
          n   <- b - a
          Sy  <- S[b + 1L] - S[a + 1L]
          Sq  <- Q[b + 1L] - Q[a + 1L]
          v[j] <- (Sq - Sy^2 / n) / n
        }
        v
      }
    } else {
      data   <- ss
      seglen <- seg.len(object)
      v      <- numeric(nseg)
      for (j in seq_len(nseg)) {
        v[j] <- var(data[(cpts[j] + 1L):cpts[j + 1L]])
      }
      v * (seglen - 1) / seglen
    }
  }
  
  param.scale <- function(object, shape) {
    cpts <- c(0, object@cpts)
    data <- sumstat(object)
    y    <- c(0, cumsum(data))
    tmpscale <- numeric(nseg(object))
    for (j in seq_len(nseg(object))) {
      tmpscale[j] <- (y[cpts[j + 1] + 1L] - y[cpts[j] + 1L]) /
        ((cpts[j + 1] - cpts[j]) * shape)
    }
    tmpscale
  }
  
  param.trend <- function(object) {
    cpts   <- c(0, object@cpts)
    seglen <- seg.len(object)
    data   <- sumstat(object)
    n      <- length(data)
    ss     <- cbind(cumsum(c(0, data)), cumsum(c(0, data * (1:n))))
    cptss  <- matrix(ss[object@cpts + 1L, ] - ss[c(0, cpts(object)) + 1L, ], ncol = 2)
    cptss[, 2L] <- cptss[, 2L] - cptss[, 1L] * c(0, cpts(object))
    
    thetaS <- (2 * cptss[, 1L] * (2 * seglen + 1) - 6 * cptss[, 2L]) /
      (2 * seglen * (2 * seglen + 1) - 3 * seglen * (seglen + 1))
    thetaT <- (6 * cptss[, 2L]) / ((seglen + 1) * (2 * seglen + 1)) +
      thetaS * (1 - (3 * seglen) / ((2 * seglen) + 1))
    cbind(thetaS, thetaT)
  }
  
  param.meanar <- function(object) {
    seglen <- seg.len(object)
    data   <- sumstat(object)
    n      <- length(data) - 1L
    ss     <- cbind(
      cumsum(c(0, data[-1L])),
      cumsum(c(0, data[-(n + 1L)])),
      cumsum(c(0, data[-1L] * data[-(n + 1L)])),
      cumsum(c(0, data[-1L]^2)),
      cumsum(c(0, data[-(n + 1L)]^2))
    )
    cptss <- matrix(ss[object@cpts + 1L, ] - ss[c(0, cpts(object)) + 1L, ], ncol = 5)
    beta2 <- (2 * seglen * cptss[, 3L] - cptss[, 1L] * cptss[, 2L]) /
      (2 * seglen * cptss[, 5L] * (1 - cptss[, 2L]^2))
    beta1 <- (2 * cptss[, 1L] - beta2 * cptss[, 2L]) / (2 * seglen)
    cbind(beta1, beta2)
  }
  
  param.trendar <- function(object) {
    seglen <- seg.len(object)
    data   <- sumstat(object)
    n      <- length(data) - 1L
    ss     <- cbind(
      cumsum(c(0, data[-1L])),
      cumsum(c(0, data[-(n + 1L)])),
      cumsum(c(0, data[-1L] * data[-(n + 1L)])),
      cumsum(c(0, data[-1L] * (1:n))),
      cumsum(c(0, data[-(n + 1L)] * 0:(n - 1L))),
      cumsum(c(0, data[-1L]^2)),
      cumsum(c(0, data[-(n + 1L)]^2))
    )
    cptss <- matrix(ss[object@cpts + 1L, ] - ss[c(0, cpts(object)) + 1L, ], ncol = 7)
    cptss[, 4L] <- cptss[, 4L] - cptss[, 1L] * c(0, cpts(object))
    cptss[, 5L] <- cptss[, 5L] - cptss[, 2L] * c(0, cpts(object))
    
    betatop <- seglen * (seglen - 1) * (
      seglen * (seglen - 1) * cptss[, 3L] +
        2 * (2 * seglen + 1) * cptss[, 1L] * (cptss[, 5L] - seglen * cptss[, 2L]) +
        6 * cptss[, 4L] * (cptss[, 2L] - cptss[, 5L])
    )
    betabottom <- seglen * (seglen - 1) * cptss[, 7L] +
      2 * (2 * seglen + 1) * cptss[, 2L] * (seglen * cptss[, 2L] - cptss[, 5L]) +
      6 * cptss[, 5L] * (cptss[, 5L] - cptss[, 2L])
    beta    <- betatop / betabottom
    thetajp <- (6 * (seglen + 2) * (cptss[, 4L] - beta * cptss[, 5L])) /
      ((seglen + 1) * (2 * seglen + 1)) -
      2 * (cptss[, 1L] - beta * cptss[, 2L])
    thetaj  <- (2 * (2 * seglen + 1) * (cptss[, 1L] - beta * cptss[, 2L]) -
                  6 * (cptss[, 4L] - beta * cptss[, 5L])) / (seglen - 1)
    cbind(beta, thetajp, thetaj)
  }
  
  ## Dispatch by cpttype/test.stat
  
  if (cpttype(object) == "mean") {
    param.est(object) <- list(mean = param.mean(object))
  } else if (cpttype(object) == "variance") {
    param.est(object) <- list(variance = param.var(object))
  } else if (cpttype(object) == "mean and variance") {
    if (test.stat(object) == "Normal") {
      param.est(object) <- list(
        mean     = param.mean(object),
        variance = param.var(object)
      )
    } else if (test.stat(object) == "Gamma") {
      param.est(object) <- list(scale = param.scale(object, shape = shape), shape = shape)
    } else if (test.stat(object) == "Exponential") {
      param.est(object) <- list(rate = 1 / param.mean(object))
    } else if (test.stat(object) == "Poisson") {
      param.est(object) <- list(lambda = param.mean(object))
    } else {
      stop("Unknown test statistic for a change in mean and variance")
    }
  } else if (cpttype(object) == "trend") {
    if (test.stat(object) == "Normal") {
      tmp <- param.trend(object)
      param.est(object) <- list(thetaS = tmp[, 1L], thetaT = tmp[, 2L])
    } else {
      stop("Unknown test statistic for a change in trend")
    }
  } else if (cpttype(object) == "trendar") {
    if (test.stat(object) == "Normal") {
      tmp <- param.trendar(object)
      param.est(object) <- list(beta = tmp[, 1L], thetajpo = tmp[, 2L], thetaj = tmp[, 3L])
    } else {
      stop("Unknown test statistic for a change in trend+ar")
    }
  } else if (cpttype(object) == "meanar") {
    if (test.stat(object) == "Normal") {
      tmp <- param.meanar(object)
      param.est(object) <- list(beta1 = tmp[, 1L], beta2 = tmp[, 2L])
    } else {
      stop("Unknown test statistic for a change in mean+ar")
    }
  } else {
    stop("Unknown changepoint type, must be 'mean', 'variance', 'mean and variance', 'trend', 'meanar' or 'trendar'.")
  }
  
  object
})

setMethod("param", "ocpt.range", function(object, ncpts = NA, shape, ...) {
  
  if (is.na(ncpts)) {
    cpts <- c(0, object@cpts)
  } else {
    ncpts.full <- apply(cpts.full(object), 1L, function(x) sum(x > 0, na.rm = TRUE))
    row <- try(which(ncpts.full == ncpts), silent = TRUE)
    if (inherits(row, "try-error") || length(row) == 0L) {
      stop("Your input object doesn't have a segmentation with the requested number of changepoints.")
    }
    cpts <- c(0, cpts.full(object)[row, 1:ncpts], as.integer(ndone(object) + nupdate(object)))
  }
  
  param.mean <- function(object, cpts) {
    nseg <- length(cpts) - 1L
    data <- sumstat(object)
    tmp  <- numeric(nseg)
    for (j in seq_len(nseg)) {
      tmp[j] <- mean(data[(cpts[j] + 1L):cpts[j + 1L]])
    }
    tmp
  }
  
  param.var <- function(object, cpts) {
    nseg   <- length(cpts) - 1L
    data   <- sumstat(object)
    seglen <- seg.len(object)
    tmp    <- numeric(nseg)
    for (j in seq_len(nseg)) {
      tmp[j] <- var(data[(cpts[j] + 1L):cpts[j + 1L]])
    }
    tmp * (seglen - 1) / seglen
  }
  
  param.scale <- function(object, cpts, shape) {
    nseg <- length(cpts) - 1L
    data <- sumstat(object)
    y    <- c(0, cumsum(data))
    tmp  <- numeric(nseg)
    for (j in seq_len(nseg)) {
      tmp[j] <- (y[cpts[j + 1] + 1L] - y[cpts[j] + 1L]) /
        ((cpts[j + 1] - cpts[j]) * shape)
    }
    tmp
  }
  
  
  if (cpttype(object) == "mean") {
    pe <- list(mean = param.mean(object, cpts))
  } else if (cpttype(object) == "variance") {
    pe <- list(variance = param.var(object, cpts))
  } else if (cpttype(object) == "mean and variance") {
    if (test.stat(object) == "Normal") {
      pe <- list(
        mean     = param.mean(object, cpts),
        variance = param.var(object, cpts)
      )
    } else if (test.stat(object) == "Gamma") {
      pe <- list(scale = param.scale(object, cpts, shape = shape), shape = shape)
    } else if (test.stat(object) == "Exponential") {
      pe <- list(rate = 1 / param.mean(object, cpts))
    } else if (test.stat(object) == "Poisson") {
      pe <- list(lambda = param.mean(object, cpts))
    } else {
      stop("Unknown test statistic for a change in mean and variance")
    }
  } else if (cpttype(object) %in% c("trend", "trendar", "meanar")) {
    stop("param() for '", cpttype(object), "' not fully implemented for ocpt.range in this trimmed version.")
  } else {
    stop("Unknown changepoint type, must be 'mean', 'variance' or 'mean and variance'")
  }
  
  if (is.na(ncpts)) {
    param.est(object) <- pe
    object
  } else {
    out <- new("ocpt.range")
    param.est(out) <- pe
    out
  }
})

setMethod("param", "ocpt.reg", function(object, shape, ...) {
  param.norm <- function(object) {
    cpts <- c(0, object@cpts)
    data <- sumstat(object)
    p    <- ncol(data) - 1L
    nb   <- nseg(object)
    tmpbeta  <- matrix(NA_real_, ncol = p, nrow = nb)
    tmpsigma <- numeric(nb)
    for (j in seq_len(nb)) {
      rhs <- paste0("data[", cpts[j] + 1L, ":", cpts[j + 1L], ",2]")
      if (p > 1L) {
        for (i in 2:p) {
          rhs <- paste0(rhs, "+data[", cpts[j] + 1L, ":", cpts[j + 1L], ",", i + 1L, "]")
        }
      }
      form   <- paste0("data[", cpts[j] + 1L, ":", cpts[j + 1L], ",1] ~ -1 + ", rhs)
      tmpfit <- eval(parse(text = paste0("lm(", form, ")")))
      tmpbeta[j, ]  <- tmpfit$coefficients
      tmpsigma[j]   <- var(tmpfit$residuals)
    }
    list(beta = tmpbeta, sig2 = tmpsigma)
  }
  
  if (test.stat(object) == "Normal") {
    param.est(object) <- param.norm(object)
  } else {
    stop("Unknown test statistic, must be 'Normal'")
  }
  object
})

## summary

setMethod("summary", "ocpt", function(object) {
  cat("Created Using changepoint.online version", object@version, "\n")
  cat("Changepoint type      : Change in", cpttype(object), "\n")
  cat("Method of analysis    :", method(object), "\n")
  cat("Test Statistic        :", test.stat(object), "\n")
  cat("Type of penalty       :", pen.type(object), "with value", pen.value(object), "\n")
  cat("Minimum Segment Length:", minseglen(object), "\n")
  cat("Maximum no. of cpts   :", ncpts.max(object), "\n")
  if (length(cpts(object)) <= 20L) {
    cat("Changepoint Locations :", cpts(object), "\n")
  } else {
    cat("Number of changepoints:", ncpts(object), "\n")
  }
  cat("ndone                 :", ndone(object), "\n")
  cat("nupdate               :", nupdate(object), "\n")
})

setMethod("summary", "ocpt.range", function(object) {
  cat("Created Using changepoint.online version", object@version, "\n")
  cat("Changepoint type      : Change in", cpttype(object), "\n")
  cat("Method of analysis    :", method(object), "\n")
  cat("Test Statistic        :", test.stat(object), "\n")
  cat("Type of penalty       :", pen.type(object), "with value", pen.value(object), "\n")
  cat("Minimum Segment Length:", minseglen(object), "\n")
  cat("Maximum no. of cpts   :", ncpts.max(object), "\n")
  if (length(cpts(object)) <= 20L) {
    cat("Changepoint Locations :", cpts(object), "\n")
  } else {
    cat("Number of changepoints:", ncpts(object), "\n")
  }
  cf <- cpts.full(object)
  if ((nrow(cf) <= 5L) && (ncol(cf) <= 20L)) {
    cat("Range of segmentations:\n")
    print(cf)
    cat("\nFor penalty values:", pen.value.full(object), "\n")
  } else {
    cat("Number of segmentations recorded:", nrow(cf),
        "with between",
        sum(cf[nrow(cf), ] > 0, na.rm = TRUE), "and",
        sum(cf[1, ] > 0, na.rm = TRUE), "changepoints.\nPenalty value ranges from:",
        min(pen.value.full(object)), "to", max(pen.value.full(object)), "\n")
  }
  cat("ndone                 :", ndone(object), "\n")
  cat("nupdate               :", nupdate(object), "\n")
})

setMethod("summary", "ocpt.reg", function(object) {
  cat("Created Using changepoint.online version", object@version, "\n")
  cat("Changepoint type      : Change in", cpttype(object), "\n")
  cat("Method of analysis    :", method(object), "\n")
  cat("Test Statistic        :", test.stat(object), "\n")
  cat("Type of penalty       :", pen.type(object), "with value", pen.value(object), "\n")
  cat("Maximum no. of cpts   :", ncpts.max(object), "\n")
  if (length(cpts(object)) <= 20L) {
    cat("Changepoint Locations :", cpts(object), "\n")
  } else {
    cat("Number of changepoints:", ncpts(object), "\n")
  }
  cat("ndone                 :", ndone(object), "\n")
  cat("nupdate               :", nupdate(object), "\n")
})

setMethod("summary", "ecp.ocpt", function(object) {
  cat("Number of Changepoints :", number(object), "\n")
  cat("Estimate Locations     :", estimates(object), "\n")
  cat("Goodness of Fit Model  :", GofM(object), "\n")
  cat("Delta                  :", delta(object), "\n")
  cat("Alpha                  :", alpha(object), "\n")
  cat("Verbose                :", verbose(object), "\n")
  cat("Number of Data points  :", datalength(object), "\n")
  cat("Calculation Time       :", functime(object), "\n")
})

setMethod("show", "ocpt", function(object) {
  cat("Class 'ocpt' : Changepoint Object\n")
  cat("       ~~   : S4 class containing", length(attributes(object)) - 1L,
      "slots with names\n")
  cat("             ", names(attributes(object))[1:(length(attributes(object)) - 1L)], "\n\n")
  cat("Created on  :", object@date, "\n\n")
  cat("summary(.)  :\n----------\n")
  summary(object)
})

setMethod("show", "ocpt.reg", function(object) {
  cat("Class 'ocpt.reg' : Changepoint Regression Object\n")
  cat("       ~~   : S4 class containing", length(attributes(object)) - 1L,
      "slots with names\n")
  cat("             ", names(attributes(object))[1:(length(attributes(object)) - 1L)], "\n\n")
  cat("Created on  :", object@date, "\n\n")
  cat("summary(.)  :\n----------\n")
  summary(object)
})

setMethod("show", "ecp.ocpt", function(object) {
  cat("Class 'ecp.ocpt' : Changepoint Object\n")
  cat("summary(.)  :\n----------\n")
  summary(object)
})

## Plot methods

setMethod(
  "plot", "ocpt",
  function(x, data, start = 1, cpt.col = "red", cpt.width = 1,
           cpt.style = 1, window = NA, ...) {
    
    if (length(param.est(x)) == 0L) {
      cat("Calculating parameter estimates...")
      x <- param(x)
      cat("done.\n")
    }
    
    if (is.na(window)) {
      plot(data, ...)
    } else {
      if (window > (ndone(x) + nupdate(x))) {
        window <- ndone(x) + nupdate(x)
      }
      start <- ndone(x) + nupdate(x) - window
      plot(data[start:(ndone(x) + nupdate(x))], ...)
    }
    
    if (cpttype(x) == "variance") {
      abline(v = cpts(x), col = cpt.col, lwd = cpt.width, lty = cpt.style)
      
    } else if (cpttype(x) == "mean" || cpttype(x) == "mean and variance") {
      
      cpts.full <- c(0, x@cpts)
      if (test.stat(x) %in% c("Normal", "CUSUM")) {
        means <- param.est(x)$mean
      } else if (test.stat(x) == "Gamma") {
        means <- param.est(x)$scale * param.est(x)$shape
      } else if (test.stat(x) == "Exponential") {
        means <- 1 / param.est(x)$rate
      } else if (test.stat(x) == "Poisson") {
        means <- param.est(x)$lambda
      } else {
        stop("Invalid Changepoint test statistic")
      }
      
      for (i in seq_len(nseg(x))) {
        segments(
          cpts.full[i] + 1L, means[i],
          cpts.full[i + 1L], means[i],
          col = cpt.col, lwd = cpt.width, lty = cpt.style
        )
      }
      
    } else if (cpttype(x) == "trend") {
      
      cpts.full <- c(0, x@cpts)
      seglen    <- x@cpts - c(0, cpts(x))
      intercept <- rep(param.est(x)$thetaS, seglen)
      slope     <- rep(param.est(x)$thetaT - param.est(x)$thetaS, seglen) /
        rep(seglen, seglen)
      cptn <- rep(c(0, cpts(x)), seglen)
      n    <- nrow(sumstat(x)) - 1L
      means <- intercept + slope * ((seq_len(n)) - cptn)
      
      for (i in seq_len(nseg(x))) {
        segments(
          cpts.full[i] + 1L, means[cpts.full[i] + 1L],
          cpts.full[i + 1L], means[cpts.full[i + 1L]],
          col = cpt.col, lwd = cpt.width, lty = cpt.style
        )
      }
      
    } else {
      stop("Invalid Changepoint Type for plotting.\n",
           "Can only plot mean, variance, mean and variance, or trend.")
    }
  }
)

setMethod(
  "plot", "ocpt.range",
  function(x, data, ncpts = NA, diagnostic = FALSE,
           cpt.col = "red", cpt.width = 1, cpt.style = 1, ...) {
    
    if (diagnostic) {
      ncpts.full <- apply(cpts.full(x), 1L, function(z) sum(z > 0, na.rm = TRUE))
      return(
        plot(
          ncpts.full,
          pen.value.full(x),
          type = "l",
          xlab = "Number of Changepoints",
          ylab = "Difference in Test Statistic",
          ...
        )
      )
    }
    
    plot(data, ...)
    
    if (is.na(ncpts)) {
      if (pen.type(x) == "CROPS") {
        stop("CROPS does not supply an optimal set of changepoints.\n",
             "Set ncpts to the desired segmentation to plot or use diagnostic=TRUE")
      }
      cpts.to.plot <- cpts(x)
      pe           <- x
    } else {
      ncpts.full <- apply(cpts.full(x), 1L, function(z) sum(z > 0, na.rm = TRUE))
      row        <- which(ncpts.full == ncpts)
      if (length(row) == 0L) {
        stop(
          "Your input object doesn't have a segmentation with the requested number of changepoints.\n",
          "Possible ncpts are: ", paste(ncpts.full, collapse = ",")
        )
      }
      cpts.to.plot <- cpts.full(x)[row, 1:ncpts]
      if (test.stat(x) == "Gamma") {
        pe <- param(x, ncpts, shape = param.est(x)$shape)
      } else {
        pe <- param(x, ncpts)
      }
    }
    
    if (cpttype(x) == "variance") {
      abline(v = cpts.to.plot, col = cpt.col, lwd = cpt.width, lty = cpt.style)
      
    } else if (cpttype(x) == "mean" || cpttype(x) == "mean and variance") {
      
      if (test.stat(x) %in% c("Normal", "CUSUM")) {
        means <- param.est(pe)$mean
      } else if (test.stat(x) == "Gamma") {
        means <- param.est(pe)$scale * param.est(pe)$shape
      } else if (test.stat(x) == "Exponential") {
        means <- 1 / param.est(pe)$rate
      } else if (test.stat(x) == "Poisson") {
        means <- param.est(pe)$lambda
      } else {
        stop("Invalid Changepoint test statistic")
      }
      
      nseg.plot   <- ncpts + 1L
      cpts.ext    <- c(0, cpts.to.plot, nrow(sumstat(x)) - 1L)
      for (i in seq_len(nseg.plot)) {
        segments(
          cpts.ext[i] + 1L, means[i],
          cpts.ext[i + 1L], means[i],
          col = cpt.col, lwd = cpt.width, lty = cpt.style
        )
      }
      
    } else {
      stop("Invalid Changepoint Type for plotting.\n",
           "Can only plot mean, variance, mean and variance.")
    }
  }
)

setMethod(
  "plot", "ocpt.reg",
  function(x, data, cpt.col = "red", cpt.width = 1, cpt.style = 1, ...) {
    
    if (length(param.est(x)) == 0L) {
      cat("Calculating parameter estimates...")
      x <- param(x)
      cat("done.\n")
    }
    
    plot(data[, 1L], type = "l", ...)
    
    if (test.stat(x) == "Normal") {
      cpts.full <- c(0, x@cpts)
      betas     <- param.est(x)$beta
      for (i in seq_len(nseg(x))) {
        idx    <- (cpts.full[i] + 1L):cpts.full[i + 1L]
        fitted <- as.numeric(betas[i, ] %*% t(data[idx, -1L, drop = FALSE]))
        lines(idx, fitted, col = cpt.col, lwd = cpt.width, lty = cpt.style)
      }
    } else {
      stop("Invalid Changepoint test statistic for ocpt.reg plot (must be Normal).")
    }
  }
)

## logLik methods (negative log-likelihood style)

setMethod(
  "logLik", "ocpt",
  function(object) {
    
    if (length(param.est(object)) == 0L) {
      cat("Calculating parameter estimates...")
      object <- param(object)
      cat("done.\n")
    }
    
    add_penalty <- function(base_like) {
      if (pen.type(object) == "MBIC") {
        base_like + (nseg(object) - 2L) * pen.value(object) + sum(log(seg.len(object)))
      } else {
        base_like + (nseg(object) - 1L) * pen.value(object)
      }
    }
    
    if (test.stat(object) == "Normal") {
      
      if (cpttype(object) == "mean") {
        
        means <- rep(param.est(object)$mean, seg.len(object))
        y     <- as.numeric(sumstat(object))
        rss   <- sum((y - means)^2)
        n     <- length(y)
        base  <- n * (log(2 * pi) + log(rss / n) + 1)
        pen   <- add_penalty(base)
        
      } else if (cpttype(object) == "variance") {
        
        y      <- as.numeric(sumstat(object))
        rss    <- c(0, cumsum((y - param.est(object)$mean)^2))
        cpts.f <- c(0, object@cpts)
        n      <- length(y)
        seglen <- seg.len(object)
        sigmas <- (rss[cpts.f[-1L] + 1L] - rss[cpts.f[-length(cpts.f)] + 1L]) / seglen
        base   <- n * log(2 * pi) + sum(seglen * log(sigmas)) + n
        pen    <- add_penalty(base)
        
      } else if (cpttype(object) == "mean and variance") {
        
        y      <- as.numeric(sumstat(object))
        means  <- rep(param.est(object)$mean, seg.len(object))
        rss    <- sum((y - means)^2)
        n      <- length(y)
        seglen <- seg.len(object)
        sigmas <- param.est(object)$variance
        base   <- n * log(2 * pi) + sum(seglen * log(sigmas)) + n
        pen    <- add_penalty(base)
        
      } else if (cpttype(object) == "trend") {
        
        cpts.f   <- c(0, object@cpts)
        seglen   <- seg.len(object)
        n        <- nrow(sumstat(object)) - 1L
        intercept <- rep(param.est(object)$thetaS, seglen)
        slope     <- rep(param.est(object)$thetaT - param.est(object)$thetaS, seglen) /
          rep(seglen, seglen)
        cptn  <- rep(c(0, cpts(object)), seglen)
        mu    <- intercept + slope * ((seq_len(n)) - cptn)
        y     <- as.numeric(sumstat(object))
        rss   <- sum((y - mu)^2)
        base  <- n * (log(2 * pi) + log(rss / n) + 1)
        pen   <- add_penalty(base)
        
      } else if (cpttype(object) == "trendar") {
        
        seglen   <- seg.len(object)
        n        <- nrow(sumstat(object)) - 1L
        intercept <- rep(param.est(object)$thetaj, seglen)
        slope     <- rep(param.est(object)$thetajpo - param.est(object)$thetaj, seglen) /
          rep(seglen, seglen)
        ar       <- rep(param.est(object)$beta, seglen)
        cptn     <- rep(c(0, cpts(object)), seglen)
        y        <- as.numeric(sumstat(object))
        mu       <- numeric(n)
        mu[1]    <- 0
        for (i in 2:n) {
          mu[i] <- intercept[i] + slope[i] * (i - cptn[i]) + ar[i] * mu[i - 1L]
        }
        mu    <- mu[-1L]
        rss   <- sum((y[-1L] - mu)^2)
        base  <- n * (log(2 * pi) + log(rss / n) + 1)
        pen   <- add_penalty(base)
        
      } else if (cpttype(object) == "meanar") {
        
        seglen   <- seg.len(object)
        n        <- nrow(sumstat(object)) - 1L
        intercept <- rep(param.est(object)$beta1, seglen)
        ar       <- rep(param.est(object)$beta2, seglen)
        cptn     <- rep(c(0, cpts(object)), seglen)
        y        <- as.numeric(sumstat(object))
        mu       <- numeric(n)
        mu[1]    <- 0
        for (i in 2:n) {
          mu[i] <- intercept[i] + ar[i] * mu[i - 1L]
        }
        mu    <- mu[-1L]
        rss   <- sum((y[-1L] - mu)^2)
        base  <- n * (log(2 * pi) + log(rss / n) + 1)
        pen   <- add_penalty(base)
        
      } else {
        stop("Unknown changepoint type for test.stat='Normal'.",
             " Must be 'mean', 'variance', 'mean and variance', 'trend', 'meanar' or 'trendar'.")
      }
      
    } else if (test.stat(object) == "Gamma") {
      
      if (cpttype(object) != "mean and variance") {
        stop("Unknown changepoint type for test.stat='Gamma', must be 'mean and variance'")
      }
      
      warning("Not changed to be -2*logLik; using original mean/var contrast.")
      mll.meanvarg <- function(x, n, shape) {
        n * shape * log(n * shape) - n * shape * log(x)
      }
      y     <- sumstat(object)[, 1L]
      shape <- param.est(object)$shape
      cpts.f <- c(0, object@cpts)
      nseg.o <- nseg(object)
      base   <- 0
      for (j in seq_len(nseg.o)) {
        base <- base + mll.meanvarg(
          y[cpts.f[j + 1L] + 1L] - y[cpts.f[j] + 1L],
          cpts.f[j + 1L] - cpts.f[j],
          shape
        )
      }
      pen <- add_penalty(base)
      
    } else if (test.stat(object) == "Exponential") {
      
      if (cpttype(object) != "mean and variance") {
        stop("Unknown changepoint type for test.stat='Exponential', must be 'mean and variance'")
      }
      
      warning("Not changed to be -2*logLik; using original mean contrast.")
      mll.meanvare <- function(x, n) {
        n * log(n) - n * log(x)
      }
      y      <- sumstat(object)[, 1L]
      cpts.f <- c(0, object@cpts)
      base   <- 0
      for (j in seq_len(nseg(object))) {
        base <- base + mll.meanvare(
          y[cpts.f[j + 1L] + 1L] - y[cpts.f[j] + 1L],
          cpts.f[j + 1L] - cpts.f[j]
        )
      }
      pen <- add_penalty(base)
      
    } else if (test.stat(object) == "Poisson") {
      
      if (cpttype(object) != "mean and variance") {
        stop("Unknown changepoint type for test.stat='Poisson', must be 'mean and variance'")
      }
      
      warning("Not changed to be -2*logLik; using original mean contrast.")
      mll.meanvarp <- function(x, n) {
        x * log(x) - x * log(n)
      }
      y      <- sumstat(object)[, 1L]
      cpts.f <- c(0, object@cpts)
      base   <- 0
      for (j in seq_len(nseg(object))) {
        base <- base + mll.meanvarp(
          y[cpts.f[j + 1L] + 1L] - y[cpts.f[j] + 1L],
          cpts.f[j + 1L] - cpts.f[j]
        )
      }
      pen <- add_penalty(base)
      
    } else {
      stop("logLik is only valid for distributional assumptions, not CUSUM or CSS")
    }
    
    like <- c(`-like` = base, `-likepen` = pen)
    like
  }
)

setMethod(
  "logLik", "ocpt.range",
  function(object, ncpts = NA) {
    
    warning("Not changed to be -2*logLik; using original contrast scale.")
    
    if (is.na(ncpts)) {
      if (pen.type(object) == "CROPS") {
        stop("CROPS does not supply an optimal set of changepoints.\n",
             "Set ncpts or use diagnostic=TRUE.")
      }
      cpts.f   <- c(0, object@cpts)
      pen.val  <- pen.value(object)
    } else {
      ncpts.full <- apply(cpts.full(object), 1L, function(z) sum(z > 0, na.rm = TRUE))
      row        <- which(ncpts.full == ncpts)
      if (length(row) == 0L) {
        stop(
          "Your input object doesn't have a segmentation with the requested number of changepoints.\n",
          "Possible ncpts are: ", paste(ncpts.full, collapse = ",")
        )
      }
      cpts.f  <- c(0, cpts.full(object)[row, 1:ncpts], as.integer(ndone(object) + nupdate(object)))
      pen.val <- pen.value.full(object)[row]
    }
    
    nseg.o <- length(cpts.f) - 1L
    
    add_penalty <- function(base_like) {
      if (pen.type(object) == "MBIC") {
        base_like + (nseg.o - 2L) * pen.val + sum(log(seg.len(object)))
      } else {
        base_like + (nseg.o - 1L) * pen.val
      }
    }
    
    if (test.stat(object) == "Normal") {
      
      if (cpttype(object) == "mean") {
        
        mll.mean <- function(x2, x, n) {
          x2 - (x^2) / n
        }
        y2   <- sumstat(object)[, 2L]
        y    <- sumstat(object)[, 1L]
        base <- 0
        for (j in seq_len(nseg.o)) {
          base <- base + mll.mean(
            y2[cpts.f[j + 1L] + 1L] - y2[cpts.f[j] + 1L],
            y[cpts.f[j + 1L] + 1L] - y[cpts.f[j] + 1L],
            cpts.f[j + 1L] - cpts.f[j]
          )
        }
        pen <- add_penalty(base)
        
      } else if (cpttype(object) == "variance") {
        
        mll.var <- function(x, n) {
          neg      <- x <= 0
          x[neg]   <- 1e-11
          n * (log(2 * pi) + log(x / n) + 1)
        }
        y2   <- sumstat(object)[, 3L]
        base <- 0
        for (j in seq_len(nseg.o)) {
          base <- base + mll.var(
            y2[cpts.f[j + 1L] + 1L] - y2[cpts.f[j] + 1L],
            cpts.f[j + 1L] - cpts.f[j]
          )
        }
        pen <- add_penalty(base)
        
      } else if (cpttype(object) == "mean and variance") {
        
        mll.meanvar <- function(x2, x, n) {
          sigmasq <- (1 / n) * (x2 - (x^2) / n)
          neg     <- sigmasq <= 0
          sigmasq[neg] <- 1e-11
          n * (log(2 * pi) + log(sigmasq) + 1)
        }
        y2   <- sumstat(object)[, 2L]
        y    <- sumstat(object)[, 1L]
        base <- 0
        for (j in seq_len(nseg.o)) {
          base <- base + mll.meanvar(
            y2[cpts.f[j + 1L] + 1L] - y2[cpts.f[j] + 1L],
            y[cpts.f[j + 1L] + 1L] - y[cpts.f[j] + 1L],
            cpts.f[j + 1L] - cpts.f[j]
          )
        }
        pen <- add_penalty(base)
        
      } else {
        stop("Unknown changepoint type, must be 'mean', 'variance' or 'mean and variance'")
      }
      
    } else if (test.stat(object) == "Gamma") {
      
      if (cpttype(object) != "mean and variance") {
        stop("Unknown changepoint type for test.stat='Gamma', must be 'mean and variance'")
      }
      mll.meanvarg <- function(x, n, shape) {
        n * shape * log(n * shape) - n * shape * log(x)
      }
      y     <- sumstat(object)[, 1L]
      shape <- param.est(object)$shape
      base  <- 0
      for (j in seq_len(nseg.o)) {
        base <- base + mll.meanvarg(
          y[cpts.f[j + 1L] + 1L] - y[cpts.f[j] + 1L],
          cpts.f[j + 1L] - cpts.f[j],
          shape
        )
      }
      pen <- add_penalty(base)
      
    } else if (test.stat(object) == "Exponential") {
      
      if (cpttype(object) != "mean and variance") {
        stop("Unknown changepoint type for test.stat='Exponential', must be 'mean and variance'")
      }
      mll.meanvare <- function(x, n) {
        n * log(n) - n * log(x)
      }
      y    <- sumstat(object)[, 1L]
      base <- 0
      for (j in seq_len(nseg.o)) {
        base <- base + mll.meanvare(
          y[cpts.f[j + 1L] + 1L] - y[cpts.f[j] + 1L],
          cpts.f[j + 1L] - cpts.f[j]
        )
      }
      pen <- add_penalty(base)
      
    } else if (test.stat(object) == "Poisson") {
      
      if (cpttype(object) != "mean and variance") {
        stop("Unknown changepoint type for test.stat='Poisson', must be 'mean and variance'")
      }
      mll.meanvarp <- function(x, n) {
        x * log(x) - x * log(n)
      }
      y    <- sumstat(object)[, 1L]
      base <- 0
      for (j in seq_len(nseg.o)) {
        base <- base + mll.meanvarp(
          y[cpts.f[j + 1L] + 1L] - y[cpts.f[j] + 1L],
          cpts.f[j + 1L] - cpts.f[j]
        )
      }
      pen <- add_penalty(base)
      
    } else {
      stop("logLik is only valid for distributional assumptions, not CUSUM or CSS")
    }
    
    like <- c(`-like` = base, `-likepen` = pen)
    like
  }
)

setMethod(
  "logLik", "ocpt.reg",
  function(object) {
    
    if (length(param.est(object)) == 0L) {
      cat("Calculating parameter estimates...")
      object <- param(object)
      cat("done.\n")
    }
    
    if (test.stat(object) != "Normal") {
      stop("logLik is only valid for Normal distributional assumption for ocpt.reg.")
    }
    
    cpts.f  <- c(0, object@cpts)
    seglen  <- seg.len(object)
    data    <- sumstat(object)
    beta    <- param.est(object)$beta
    sigmas  <- param.est(object)$sig2
    rss.vec <- numeric(length(seglen))
    
    for (i in seq_along(seglen)) {
      idx <- (cpts.f[i] + 1L):cpts.f[i + 1L]
      mu  <- as.numeric(data[idx, -1L, drop = FALSE] %*% beta[i, ])
      rss.vec[i] <- sum((data[idx, 1L] - mu)^2)
    }
    
    base <- sum(seglen * log(2 * pi * sigmas)) + sum(rss.vec / sigmas)
    
    if (pen.type(object) == "MBIC") {
      pen <- base + (nseg(object) - 2L) * pen.value(object) + sum(log(seg.len(object)))
    } else {
      pen <- base + (nseg(object) - 1L) * pen.value(object)
    }
    
    like <- c(`-like` = base, `-likepen` = pen)
    like
  }
)

if (!isGeneric("likelihood")) {
  setGeneric("likelihood", function(object) standardGeneric("likelihood"))
}

setMethod("likelihood", "ocpt", function(object) {
  logLik(object)
})

setGeneric("acf", function(object, ...) standardGeneric("acf"))

setMethod(
  "acf", "ocpt",
  function(object, lag.max = NULL, ...) {
    cpts.f <- c(0, object@cpts)
    nseg.o <- nseg(object)
    data   <- as.numeric(sumstat(object))
    for (i in seq_len(nseg.o)) {
      stats::acf(
        data[(cpts.f[i] + 1L):cpts.f[i + 1L]],
        main = paste("Series part:", (cpts.f[i] + 1L), ":", cpts.f[i + 1L]),
        lag.max = lag.max,
        ...
      )
    }
  }
)

setMethod(
  "acf", "ocpt.reg",
  function(object, lag.max = NULL, ...) {
    cpts.f <- c(0, object@cpts)
    nseg.o <- nseg(object)
    data   <- sumstat(object)[, 1L]
    for (i in seq_len(nseg.o)) {
      stats::acf(
        data[(cpts.f[i] + 1L):cpts.f[i + 1L]],
        main = paste("Series part:", (cpts.f[i] + 1L), "-", cpts.f[i + 1L]),
        lag.max = lag.max,
        ...
      )
    }
  }
)

## ecp.ocpt accessors and replacement methods

if (!isGeneric("number")) {
  setGeneric("number", function(object) standardGeneric("number"))
}
setMethod("number", "ecp.ocpt", function(object) object@number)

if (!isGeneric("estimates")) {
  setGeneric("estimates", function(object) standardGeneric("estimates"))
}
setMethod("estimates", "ecp.ocpt", function(object) object@estimates)

if (!isGeneric("GofM")) {
  setGeneric("GofM", function(object) standardGeneric("GofM"))
}
setMethod("GofM", "ecp.ocpt", function(object) object@GofM)

if (!isGeneric("delta")) {
  setGeneric("delta", function(object) standardGeneric("delta"))
}
setMethod("delta", "ecp.ocpt", function(object) object@delta)

if (!isGeneric("alpha")) {
  setGeneric("alpha", function(object) standardGeneric("alpha"))
}
setMethod("alpha", "ecp.ocpt", function(object) object@alpha)

if (!isGeneric("verbose")) {
  setGeneric("verbose", function(object) standardGeneric("verbose"))
}
setMethod("verbose", "ecp.ocpt", function(object) object@verbose)

if (!isGeneric("csum")) {
  setGeneric("csum", function(object) standardGeneric("csum"))
}
setMethod("csum", "ecp.ocpt", function(object) object@csum)

if (!isGeneric("dll")) {
  setGeneric("dll", function(object) standardGeneric("dll"))
}
setMethod("dll", "ecp.ocpt", function(object) object@dll)

if (!isGeneric("dlr")) {
  setGeneric("dlr", function(object) standardGeneric("dlr"))
}
setMethod("dlr", "ecp.ocpt", function(object) object@dlr)

if (!isGeneric("drr")) {
  setGeneric("drr", function(object) standardGeneric("drr"))
}
setMethod("drr", "ecp.ocpt", function(object) object@drr)

if (!isGeneric("left")) {
  setGeneric("left", function(object) standardGeneric("left"))
}
setMethod("left", "ecp.ocpt", function(object) object@left)

if (!isGeneric("right")) {
  setGeneric("right", function(object) standardGeneric("right"))
}
setMethod("right", "ecp.ocpt", function(object) object@right)

if (!isGeneric("datalength")) {
  setGeneric("datalength", function(object) standardGeneric("datalength"))
}
setMethod("datalength", "ecp.ocpt", function(object) object@datalength)

if (!isGeneric("functime")) {
  setGeneric("functime", function(object) standardGeneric("functime"))
}
setMethod("functime", "ecp.ocpt", function(object) object@functime)

if (!isGeneric("width")) {
  setGeneric("width", function(object) standardGeneric("width"))
}
setMethod("width", "ecp.ocpt", function(object) object@width)

if (!isGeneric("cpLoc")) {
  setGeneric("cpLoc", function(object) standardGeneric("cpLoc"))
}
setMethod("cpLoc", "ecp.ocpt", function(object) object@cpLoc)

## Replacement methods for ecp.ocpt slots -----------------------

setGeneric("number<-", function(object, value) standardGeneric("number<-"))
setReplaceMethod("number", "ecp.ocpt", function(object, value) {
  object@number <- value
  object
})

setGeneric("estimates<-", function(object, value) standardGeneric("estimates<-"))
setReplaceMethod("estimates", "ecp.ocpt", function(object, value) {
  object@estimates <- value
  object
})

setGeneric("GofM<-", function(object, value) standardGeneric("GofM<-"))
setReplaceMethod("GofM", "ecp.ocpt", function(object, value) {
  object@GofM <- value
  object
})

setGeneric("delta<-", function(object, value) standardGeneric("delta<-"))
setReplaceMethod("delta", "ecp.ocpt", function(object, value) {
  object@delta <- value
  object
})

setGeneric("alpha<-", function(object, value) standardGeneric("alpha<-"))
setReplaceMethod("alpha", "ecp.ocpt", function(object, value) {
  object@alpha <- value
  object
})

setGeneric("verbose<-", function(object, value) standardGeneric("verbose<-"))
setReplaceMethod("verbose", "ecp.ocpt", function(object, value) {
  object@verbose <- value
  object
})

setGeneric("csum<-", function(object, value) standardGeneric("csum<-"))
setReplaceMethod("csum", "ecp.ocpt", function(object, value) {
  object@csum <- value
  object
})

setGeneric("dll<-", function(object, value) standardGeneric("dll<-"))
setReplaceMethod("dll", "ecp.ocpt", function(object, value) {
  object@dll <- value
  object
})

setGeneric("dlr<-", function(object, value) standardGeneric("dlr<-"))
setReplaceMethod("dlr", "ecp.ocpt", function(object, value) {
  object@dlr <- value
  object
})

setGeneric("drr<-", function(object, value) standardGeneric("drr<-"))
setReplaceMethod("drr", "ecp.ocpt", function(object, value) {
  object@drr <- value
  object
})

setGeneric("left<-", function(object, value) standardGeneric("left<-"))
setReplaceMethod("left", "ecp.ocpt", function(object, value) {
  object@left <- value
  object
})

setGeneric("right<-", function(object, value) standardGeneric("right<-"))
setReplaceMethod("right", "ecp.ocpt", function(object, value) {
  object@right <- value
  object
})

setGeneric("datalength<-", function(object, value) standardGeneric("datalength<-"))
setReplaceMethod("datalength", "ecp.ocpt", function(object, value) {
  object@datalength <- value
  object
})

setGeneric("functime<-", function(object, value) standardGeneric("functime<-"))
setReplaceMethod("functime", "ecp.ocpt", function(object, value) {
  object@functime <- value
  object
})

setGeneric("width<-", function(object, value) standardGeneric("width<-"))
setReplaceMethod("width", "ecp.ocpt", function(object, value) {
  object@width <- value
  object
})

setGeneric("cpLoc<-", function(object, value) standardGeneric("cpLoc<-"))
setReplaceMethod("cpLoc", "ecp.ocpt", function(object, value) {
  object@cpLoc <- value
  object
})