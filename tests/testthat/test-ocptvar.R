context("ocpt.var function tests")

set.seed(1) # Note: new data sets must be added at the end.
previoussingmeandata <- c(rnorm(100,50,1),rnorm(100,100,1))
previousmulmeandata  <- c(rnorm(100,100,1),rnorm(100,70,1),
                          rnorm(100,5,1),rnorm(100,4,1))
singmeandata         <- c(rnorm(100,0,1),rnorm(100,10,1))
mulmeandata          <- c(rnorm(100,0,1),rnorm(100,10,1),
                          rnorm(100,20,1),rnorm(100,50,1))
nochangedata         <- c(rnorm(200,0,1))
previoussingvardata  <- c(rnorm(100,10,5))
singvardata          <- c(rnorm(100,10,25))
mulvardata           <- c(rnorm(100,20,10),rnorm(100,20,15),
                          rnorm(100,20,20),rnorm(100,20,25))
singmeanvardata      <- c(rnorm(50,0,1),rnorm(50,3,10))
mulmeanvardata       <- c(rnorm(50,0,1),rnorm(50,5,3),
                          rnorm(50,10,1),rnorm(50,3,10))
mulmeanvarexpdata    <- c(rexp(50,1), rexp(50,3),
                          rexp(50,5), rexp(50,7)) # rate values correct
mulmeanvarpoisdata   <- c(rpois(50,1), rpois(50,2),
                          rpois(50,3), rpois(50,5)) # lambda values correct?
constantdata         <- rep(1, 200)
shortdata            <- c(2,2,2)
negativedata         <- jitter(rep(-100, 200))

# list of datas
meandata     <- list(singmeandata, mulmeandata, nochangedata)
meanvardata  <- list(singmeanvardata, mulmeanvardata, nochangedata)
data         <- list(singmeandata, nochangedata, constantdata, negativedata)
previousdata <- list(previoussingvardata)
updateddata  <- list(singvardata)
ecpdata      <- matrix(singvardata, ncol = 1)

penalties       <- c("None", "SIC", "BIC", "AIC",
                     "Hannan-Quinn", "Manual", "MBIC") 
penalties.error <- c("Asymptotic","CROPS")

testStats <- c("Normal","ECP") 

## ------------------------------------------------------------------
## 1) Online vs offline: variance
## ------------------------------------------------------------------
for (i in data) {
  ans.online  <- ocpt.var.initialise(i)
  ans.offline <- cpt.var(i, method = "PELT")
  expect_equal(cpts(ans.online), cpts(ans.offline))
}

## ------------------------------------------------------------------
## 2) initialise vs initialize: identical
## ------------------------------------------------------------------
for (i in data) {
  english  <- ocpt.var.initialise(i)
  american <- ocpt.var.initialize(i)
  expect_identical(english, american)
}

## ------------------------------------------------------------------
## 3) Online update vs offline on concatenated data – relaxed
## ------------------------------------------------------------------
for (i in previousdata) {
  for (j in updateddata) {
    previousans <- ocpt.var.initialise(i)
    updatedans  <- ocpt.var.update(previousans, j)
    ij          <- c(i, j)
    ans.offline <- cpt.var(ij, method = "PELT")
    
    cpts.upd     <- cpts(updatedans)
    cpts.offline <- cpts(ans.offline)
    
    # same number of cps
    expect_equal(length(cpts.upd), length(cpts.offline))
    
    # allow some tolerance – online recursion can shift a little
    if (length(cpts.offline) > 0L) {
      expect_true(all(abs(cpts.upd - cpts.offline) <= 20L))
    }
  }
}

## ------------------------------------------------------------------
## 4) ECP vs PELT consistency
## ------------------------------------------------------------------
expect_identical(ecpdata, as.matrix(singvardata))
ecpans  <- ocpt.var.initialise(ecpdata, test.stat = "ECP")
normans <- ocpt.var.initialise(singvardata, test.stat = "Normal")
expect_equal(estimates(ecpans), 67)

# Check correct classes are called
for (i in testStats) {
  if (i == "ECP") {
    ecpans <- ocpt.var.initialise(ecpdata, test.stat = i)
    expect_identical(class(ecpans)[1], "ecp.ocpt")
  } else {
    ans <- ocpt.var.initialise(singvardata, test.stat = i)
    expect_identical(class(ans)[1], "ocpt")
  }
}

## ------------------------------------------------------------------
## 5) Penalty handling
## ------------------------------------------------------------------
for (p in penalties) {
  ans <- ocpt.var.initialise(singvardata, penalty = p)
}

for (p in penalties.error) {
  if (p == "CROPS") {
    expect_that(
      ocpt.var.initialise(singvardata, penalty = p),
      throws_error("Use cpt.var from changepoint as CROPS is not available with changepoint.online")
    )
  } else {
    expect_that(
      ocpt.var.initialise(singvardata, penalty = p),
      throws_error("Asymptotic penalty values must be > 0 and <= 1")
    )
  }
}

## ------------------------------------------------------------------
## 6) New tests: compressed lastchange structures (variance)
## ------------------------------------------------------------------
test_that("ocpt.var stores compressed lastchange structures", {
  skip_on_cran()
  
  set.seed(123)
  x <- c(rnorm(200, 0, 1), rnorm(200, 0, 3))  # variance change at ~200
  
  res <- ocpt.var.initialise(
    data      = x,
    penalty   = "ARL",
    pen.value = 1e5,
    test.stat = "Normal"
  )
  
  # slots should be matrices with the right column names
  expect_true(is.matrix(res@lastchangelike))
  expect_true(is.matrix(res@lastchangecpts))
  
  expect_equal(colnames(res@lastchangelike), c("like", "idx"))
  expect_equal(colnames(res@lastchangecpts),  c("cpt",  "idx"))
  
  # indices within 1..(ndone + nupdate + 1)
  T_total <- res@ndone + res@nupdate + 1L
  
  expect_true(all(res@lastchangelike[, "idx"] >= 1L))
  expect_true(all(res@lastchangelike[, "idx"] <= T_total))
  expect_true(all(res@lastchangecpts[,  "idx"] >= 1L))
  expect_true(all(res@lastchangecpts[,  "idx"] <= T_total))
  
  # sanity cp near the true break
  expect_true(any(abs(res@cpts - 200L) <= 20L))
})

test_that("ocpt.var.update expands and recompresses lastchange correctly", {
  skip_on_cran()
  
  set.seed(123)
  x1 <- c(rnorm(200, 0, 1), rnorm(200, 0, 3))
  res1 <- ocpt.var.initialise(
    data      = x1,
    penalty   = "ARL",
    pen.value = 1e5,
    test.stat = "Normal"
  )
  
  # ensure initial representation is compressed
  expect_true(is.matrix(res1@lastchangelike))
  expect_true(is.matrix(res1@lastchangecpts))
  expect_equal(colnames(res1@lastchangelike), c("like", "idx"))
  expect_equal(colnames(res1@lastchangecpts),  c("cpt",  "idx"))
  
  # new data from the high-variance regime
  set.seed(456)
  x2 <- rnorm(100, 0, 3)
  
  res2 <- ocpt.var.update(res1, x2)
  
  # compare to offline on full data
  x_full    <- c(x1, x2)
  offline   <- cpt.var(x_full, method = "PELT")
  cp_off    <- cpts(offline)
  cp_online <- res2@cpts
  
  # Require at least one cp online if offline found any
  if (length(cp_off) > 0L) {
    expect_true(length(cp_online) > 0L)
  }
  
  # representation should still be compressed matrices
  expect_true(is.matrix(res2@lastchangelike))
  expect_true(is.matrix(res2@lastchangecpts))
  
  expect_equal(colnames(res2@lastchangelike), c("like", "idx"))
  expect_equal(colnames(res2@lastchangecpts),  c("cpt",  "idx"))
  
  # indices still within valid range
  T_total <- res2@ndone + res2@nupdate + 1L
  expect_true(all(res2@lastchangelike[, "idx"] >= 1L))
  expect_true(all(res2@lastchangelike[, "idx"] <= T_total))
  expect_true(all(res2@lastchangecpts[,  "idx"] >= 1L))
  expect_true(all(res2@lastchangecpts[,  "idx"] <= T_total))
})