test_that("computeqic returns proper vector", {
  data("Arterial", package = "CrossCarry")
  cc <- createCarry(Arterial, "Treatment", "Period", "Subject", TRUE)
  fit <- CrossGEE("Pressure", "Treatment", "Period", "Subject",
                  carry = cc$carryover, data = cc$data,
                  correlation = "exchangeable", Mv = 1)
  q <- fit$QIC
  expect_type(q, "list")
  expect_true(all(!is.na(q)))
  expect_true(length(q) >= 1)
  expect_true(!is.null(names(q)))
})
