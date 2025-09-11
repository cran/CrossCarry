test_that("Kronecker correlation matrices are valid", {
  skip_on_cran()
  data("Arterial", package = "CrossCarry")
  cc <- createCarry(Arterial, "Treatment", "Period", "Subject", TRUE)

  fit <- CrossGEEKron(response="Pressure", treatment="Treatment",
                      period="Period", id="Subject",
                      time = "Time", carry = cc$carryover, data = cc$data,
                      correlation = "AR-M", Mv = 1)

  Rw <- fit$Within
  Rb <- fit$Between

  for (M in list(Rw, Rb)) {
    expect_true(is.matrix(M))
    expect_equal(diag(M), rep(1, nrow(M)))
    expect_lt(max(abs(M - t(M))), 1e-8)        # symmetric
    expect_true(all(M >= -1 & M <= 1))
  }
})
