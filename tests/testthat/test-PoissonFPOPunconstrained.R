library(testthat)
library(PeakSegOptimal)

select_segmentor_model <- function(seg, penalty) {
  costs <- seg@likelihood
  K_values <- seq_len(nrow(costs))
  penalized <- costs[,1] + penalty * (K_values - 1)
  best_K <- which.min(penalized)
  breaks <- Segmentor3IsBack::getBreaks(seg)[best_K, seq_len(best_K)]
  means  <- Segmentor3IsBack::getParameters(seg)[best_K, seq_len(best_K)]
  list(K = best_K, breaks = breaks, means = means,
       cost = costs[best_K, 1])
}

fpop_to_breaks <- function(result, n) {
  K <- result$n.segments
  if(K == 1) return(n)
  c(result$seg.end[2:K], n)
}

expand_means <- function(result, n) {
  fitted <- numeric(n)
  seg.start <- 1L
  for(k in seq_len(result$n.segments)) {
    if(k < result$n.segments) {
      seg.stop <- result$seg.end[k + 1L]
    } else {
      seg.stop <- n
    }
    fitted[seg.start:seg.stop] <- result$seg.mean[k]
    seg.start <- seg.stop + 1L
  }
  fitted
}

if(requireNamespace("Segmentor3IsBack", quietly = TRUE)) {

  test_that("3-segment Poisson data matches Segmentor", {
    set.seed(42)
    data_vec <- as.integer(c(rpois(30, 5), rpois(30, 20), rpois(30, 5)))
    n <- length(data_vec)
    seg <- Segmentor3IsBack::Segmentor(data_vec, model = 1, Kmax = 30)
    for (pen in c(5, 10, 50, 100)) {
      fpop <- PoissonFPOPunconstrained(data_vec, penalty = pen)
      segm <- select_segmentor_model(seg, pen)
      fpop_breaks <- fpop_to_breaks(fpop, n)
      expect_identical(as.integer(fpop$n.segments), as.integer(segm$K),
                       info = paste("n.segments at penalty =", pen))
      expect_identical(as.integer(fpop_breaks), as.integer(segm$breaks),
                       info = paste("breaks at penalty =", pen))
    }
  })

  test_that("single changepoint data matches Segmentor", {
    set.seed(123)
    data_vec <- as.integer(c(rpois(50, 3), rpois(50, 15)))
    n <- length(data_vec)
    seg <- Segmentor3IsBack::Segmentor(data_vec, model = 1, Kmax = 20)
    for (pen in c(10, 50, 200)) {
      fpop <- PoissonFPOPunconstrained(data_vec, penalty = pen)
      segm <- select_segmentor_model(seg, pen)
      fpop_breaks <- fpop_to_breaks(fpop, n)
      expect_identical(as.integer(fpop$n.segments), as.integer(segm$K),
                       info = paste("n.segments at penalty =", pen))
      expect_identical(as.integer(fpop_breaks), as.integer(segm$breaks),
                       info = paste("breaks at penalty =", pen))
    }
  })

  test_that("constant-rate data gives 1 segment at high penalty", {
    set.seed(7)
    data_vec <- as.integer(rpois(80, 10))
    n <- length(data_vec)
    seg <- Segmentor3IsBack::Segmentor(data_vec, model = 1, Kmax = 10)
    fpop <- PoissonFPOPunconstrained(data_vec, penalty = 500)
    segm <- select_segmentor_model(seg, 500)
    expect_identical(as.integer(fpop$n.segments), 1L)
    expect_identical(as.integer(fpop$n.segments), as.integer(segm$K))
  })

  test_that("many-changepoint data matches Segmentor at low penalty", {
    set.seed(99)
    data_vec <- as.integer(c(
      rpois(20, 2), rpois(20, 15), rpois(20, 4),
      rpois(20, 25), rpois(20, 3)
    ))
    n <- length(data_vec)
    seg <- Segmentor3IsBack::Segmentor(data_vec, model = 1, Kmax = 30)
    for (pen in c(1, 5, 20, 100)) {
      fpop <- PoissonFPOPunconstrained(data_vec, penalty = pen)
      segm <- select_segmentor_model(seg, pen)
      fpop_breaks <- fpop_to_breaks(fpop, n)
      expect_identical(as.integer(fpop$n.segments), as.integer(segm$K),
                       info = paste("n.segments at penalty =", pen))
      expect_identical(as.integer(fpop_breaks), as.integer(segm$breaks),
                       info = paste("breaks at penalty =", pen))
    }
  })

} else {
  message("Segmentor3IsBack not installed, skipping comparison tests")
}

test_that("penalty=0 gives every point its own segment", {
  data_vec <- as.integer(c(3, 10, 3, 10, 3))
  fpop <- PoissonFPOPunconstrained(data_vec, penalty = 0.0)
  expect_identical(as.integer(fpop$n.segments), length(data_vec))
})

test_that("n=2 works for both 1-segment and 2-segment outcomes", {
  data_vec <- as.integer(c(1, 100))
  fpop2 <- PoissonFPOPunconstrained(data_vec, penalty = 0.0)
  expect_identical(as.integer(fpop2$n.segments), 2L)
  fpop1 <- PoissonFPOPunconstrained(data_vec, penalty = 1e10)
  expect_identical(as.integer(fpop1$n.segments), 1L)
})

test_that("rejects negative data", {
  expect_error(PoissonFPOPunconstrained(as.integer(c(-1, 2, 3)), penalty = 5))
})
