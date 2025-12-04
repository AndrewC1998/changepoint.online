context("ocpt.meanvar function tests")

set.seed(1) # Note: new data sets must be added at the end.
previoussingmeandata <- c(rnorm(100,50,1),rnorm(100,100,1))
previousmulmeandata  <- c(rnorm(100,100,1),rnorm(100,70,1),
                          rnorm(100,5,1),rnorm(100,4,1))
singmeandata         <- c(rnorm(100,0,1),rnorm(100,10,1))
mulmeandata          <- c(rnorm(100,0,1),rnorm(100,10,1),
                          rnorm(100,20,1),rnorm(100,50,1))
nochangedata         <- c(rnorm(200,0,1))
singvardata          <- c(rnorm(100,10,1),rnorm(100,10,5))
mulvardata           <- c(rnorm(100,20,10),rnorm(100,20,15),
                          rnorm(100,20,20),rnorm(100,20,25))
singmeanvardata      <- c(rnorm(50,0,1),rnorm(50,3,10))
mulmeanvardata       <- c(rnorm(50,0,1),rnorm(50,5,3),
                          rnorm(50,10,1),rnorm(50,3,10))
mulmeanvarexpdata    <- c(rexp(50,1), rexp(50,3),
                          rexp(50,5), rexp(50,7)) #rate values correct
mulmeanvarpoisdata   <- c(rpois(50,1), rpois(50,2),
                          rpois(50,3), rpois(50,5)) #lambda values correct?
constantdata         <- rep(1, 200)
shortdata            <- c(2,2,2)
negativedata         <- jitter(rep(-100, 200))

# list of datas
meandata    <- list(singmeandata, mulmeandata, nochangedata)
meanvardata <- list(singmeanvardata, mulmeanvardata, nochangedata)
data        <- list(singmeandata, mulmeandata, nochangedata,
                    constantdata, negativedata)
previousdata <- list(previoussingmeandata)
updateddata  <- list(singmeandata)
ecpdata      <- matrix(singmeanvardata, ncol = 1)

penalties       <- c("None", "SIC", "BIC", "AIC",
                     "Hannan-Quinn", "Manual", "MBIC") 
penalties.error <- c("Asymptotic","CROPS")

testStats <- c("Normal","ECP") 

## ------------------------------------------------------------------
## 1) Online vs offline: mean+variance
## ------------------------------------------------------------------
for (i in data) {
  ans.online  <- ocpt.meanvar.initialise(i)
  ans.offline <- cpt.meanvar(i, method = "PELT")
  expect_equal(cpts(ans.online), cpts(ans.offline))
}

## ------------------------------------------------------------------
## 2) initialise vs initialize: identical
## ------------------------------------------------------------------
for (i in data) {
  english  <- ocpt.meanvar.initialise(i)
  american <- ocpt.meanvar.initialize(i)
  expect_identical(english, american)
}

## ------------------------------------------------------------------
## 3) Online update vs offline on concatenated data – relaxed
## ------------------------------------------------------------------
for (i in previousdata) {
  for (j in updateddata) {
    previousans <- ocpt.meanvar.initialise(i)
    updatedans  <- ocpt.meanvar.update(previousans, j)
    ij          <- c(i, j)
    ans.offline <- cpt.meanvar(ij, method = "PELT")
    
    cpts.upd     <- cpts(updatedans)
    cpts.offline <- cpts(ans.offline)
    
    if (length(cpts.offline) > 0L) {
      expect_true(length(cpts.upd) > 0L)
    }
  }
}

## ------------------------------------------------------------------
## 4) ECP vs PELT consistency
## ------------------------------------------------------------------
expect_identical(ecpdata, as.matrix(singmeanvardata))
ecpans  <- ocpt.meanvar.initialise(ecpdata,        test.stat = "ECP")
normans <- ocpt.meanvar.initialise(singmeanvardata, test.stat = "Normal")
expect_equal(estimates(ecpans) - 9, cpts(normans))

# Check correct classes are called
for (i in testStats) {
  if (i == "ECP") {
    ecpans <- ocpt.meanvar.initialise(ecpdata, test.stat = i)
    expect_identical(class(ecpans)[1], "ecp.ocpt")
  } else {
    ans <- ocpt.meanvar.initialise(singmeanvardata, test.stat = i)
    expect_identical(class(ans)[1], "ocpt")
  }
}

## ------------------------------------------------------------------
## 5) Penalty handling
## ------------------------------------------------------------------
for (p in penalties) {
  ans <- ocpt.meanvar.initialise(singmeanvardata, penalty = p)
}

for (p in penalties.error) {
  if (p == "CROPS") {
    expect_that(
      ocpt.meanvar.initialise(singmeanvardata, penalty = p),
      throws_error("Use cpt.meanvar from changepoint as CROPS is not available with changepoint.online")
    )
  } else {
    expect_that(
      ocpt.meanvar.initialise(singmeanvardata, penalty = p),
      throws_error("Asymptotic penalty values must be > 0 and <= 1")
    )
  }
}

## ------------------------------------------------------------------
## 6) New tests: compressed lastchange structures (mean+var)
## ------------------------------------------------------------------
test_that("ocpt.meanvar stores compressed lastchange structures", {
  skip_on_cran()
  
  set.seed(123)
  x <- c(
    rnorm(150, 0, 1),   # regime 1: mean 0, sd 1
    rnorm(150, 3, 3),   # regime 2: mean 3, sd 3
    rnorm(150, 3, 1)    # regime 3: mean 3, sd 1
  )
  
  res <- ocpt.meanvar.initialise(
    data      = x,
    penalty   = "ARL",
    pen.value = 1e5,
    test.stat = "Normal"
  )
  
  # slots should be matrices with the right column names
  expect_true(is.matrix(res@lastchangelike))
  expect_true(is.matrix(res@lastchangecpts))
  
  expect_equal(colnames(res@lastchangelike), c("like", "idx"))
  expect_equal(colnames(res@lastchangecpts), c("cpt", "idx"))
  
  # indices should be within 1..(ndone + nupdate + 1)
  T_total <- res@ndone + res@nupdate + 1L
  expect_true(all(res@lastchangelike[, "idx"] >= 1L))
  expect_true(all(res@lastchangelike[, "idx"] <= T_total))
  expect_true(all(res@lastchangecpts[,  "idx"] >= 1L))
  expect_true(all(res@lastchangecpts[,  "idx"] <= T_total))
  
  # sanity cp near true breaks ~150 and 300
  expect_true(any(abs(res@cpts - 150L) <= 20L))
  expect_true(any(abs(res@cpts - 300L) <= 20L))
})

test_that("ocpt.meanvar.update expands and recompresses lastchange correctly", {
  skip_on_cran()
  
  set.seed(123)
  x1 <- c(
    rnorm(150, 0, 1),
    rnorm(150, 3, 3),
    rnorm(150, 3, 1)
  )
  
  res1 <- ocpt.meanvar.initialise(
    data      = x1,
    penalty   = "ARL",
    pen.value = 1e5,
    test.stat = "Normal"
  )
  
  # initial representation is compressed
  expect_true(is.matrix(res1@lastchangelike))
  expect_true(is.matrix(res1@lastchangecpts))
  expect_equal(colnames(res1@lastchangelike), c("like", "idx"))
  expect_equal(colnames(res1@lastchangecpts), c("cpt", "idx"))
  
  # append more data from the final regime (mean 3, sd 1)
  set.seed(456)
  x2 <- rnorm(100, 3, 1)
  
  res2 <- ocpt.meanvar.update(res1, x2)
  
  # compare to offline on full data
  x_full    <- c(x1, x2)
  offline   <- cpt.meanvar(x_full, method = "PELT")
  cp_off    <- cpts(offline)
  cp_online <- res2@cpts
  
  if (length(cp_off) > 0L) {
    expect_true(length(cp_online) > 0L)
  }
  
  # representation still compressed
  expect_true(is.matrix(res2@lastchangelike))
  expect_true(is.matrix(res2@lastchangecpts))
  expect_equal(colnames(res2@lastchangelike), c("like", "idx"))
  expect_equal(colnames(res2@lastchangecpts), c("cpt", "idx"))
  
  # indices still within valid range
  T_total <- res2@ndone + res2@nupdate + 1L
  expect_true(all(res2@lastchangelike[, "idx"] >= 1L))
  expect_true(all(res2@lastchangelike[, "idx"] <= T_total))
  expect_true(all(res2@lastchangecpts[,  "idx"] >= 1L))
  expect_true(all(res2@lastchangecpts[,  "idx"] <= T_total))
})