test_that("Edge cases are handled gracefully", {
  data("Arterial", package = "CrossCarry")
  cc <- createCarry(Arterial, "Treatment", "Period", "Subject", TRUE)

  # One subject only
  d1 <- subset(Arterial, Subject == 1)
  cc1 <- createCarry(d1, "Treatment", "Period", "Subject", TRUE)
  expect_error(
    CrossGEE("Pressure", "Treatment", "Period", "Subject",
             carry = cc1$carryover, data = cc1$data,
             correlation = "exchangeable", Mv = 1)
  )

  # Missing values in response
  d2 <- Arterial
  d2$Pressure[1:3] <- NA
  cc2 <- createCarry(d2, "Treatment", "Period", "Subject", TRUE)
  expect_no_error(
    CrossGEE("Pressure", "Treatment", "Period", "Subject",
             carry = cc2$carryover, data = cc2$data,
             correlation = "AR-M", Mv = 1)
  )
})
